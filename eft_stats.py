"""Get raster of EFT types."""
from datetime import datetime
import argparse
import glob
import logging
import multiprocessing
import os

from osgeo import gdal
from ecoshard import geoprocessing
from ecoshard import taskgraph
import numpy

gdal.SetCacheMax(2**26)
logging.basicConfig(
    level=logging.DEBUG,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(funcName)s:%(lineno)d] %(message)s'))
logging.getLogger('ecoshard.taskgraph').setLevel(logging.INFO)
LOGGER = logging.getLogger(__name__)

NODATA = -1


def _add_thresholds(base_raster_path_list, target_threshold_count_raster_path):
    """Add all the base counting >0 as 1 to target threshold count."""
    def _add_threshold_op(*base_array_list):
        result = numpy.full(
            base_array_list[0].shape, NODATA, dtype=numpy.int32)
        for base_array in base_array_list:
            result[base_array > 0] += 1
        return result
    geoprocessing.raster_calculator(
        [(path, 1) for path in base_raster_path_list],
        _add_threshold_op,
        target_threshold_count_raster_path, gdal.GDT_Int32, NODATA)


def _basename(path):
    """Return basename without extension."""
    return os.path.basename(os.path.splitext(path)[0])


def _mask_by_value(base_raster_path, mask_value, target_raster_path):
    """Mask base by value to target."""
    nodata = geoprocessing.get_raster_info(base_raster_path)['nodata'][0]
    target_nodata = 2
    def _mask_op(array):
        """Mask array by mask value but keep nodata."""
        result = numpy.full(array.shape, target_nodata, dtype=numpy.int8)
        if nodata is not None:
            valid_mask = (array != nodata)
            print('nnot non')
        else:
            valid_mask = numpy.isfinite(array)
        result[valid_mask] = (array[valid_mask] == mask_value).astype(numpy.int8)
        return result

    geoprocessing.raster_calculator(
        [(base_raster_path, 1)], _mask_op,
        target_raster_path, gdal.GDT_Byte, target_nodata)


def _make_radius_kernel(n_pixels, target_kernel_path):
    """Create a flat radius kernel at target kernel path."""
    nodata = -1
    kernel = numpy.arange(0, n_pixels*2)
    circle_mask = (
        (kernel[numpy.newaxis, :]-n_pixels)**2 +
        (kernel[:, numpy.newaxis]-n_pixels)**2) < n_pixels**2
    geoprocessing.numpy_array_to_raster(
        circle_mask, nodata, (1, -1), (1, 1), None,
        target_kernel_path)


def main():
    """Entry point."""
    parser = argparse.ArgumentParser(
        description='EFT counter')
    parser.add_argument(
        'raster_pattern', type=str, nargs='+',
        help='path to raster or wildcards to process')
    parser.add_argument(
        '--search_radius', type=float, nargs='+', required=True,
        help="List of search radii to search in projected units")
    parser.add_argument(
        '--force', action='store_true', help=(
            'override any errors, if you do not know what this is you do '
            'not need it'))
    args = parser.parse_args()

    workspace_dir = '_eft_stats_workspace'
    os.makedirs(workspace_dir, exist_ok=True)
    task_graph = taskgraph.TaskGraph(
        workspace_dir, multiprocessing.cpu_count(), 10.0,
        parallel_mode='thread')

    # TODO: was I going to make stats?
    stats_path = (
        f'eft_stats_{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}.csv')
    raster_stats_list = []
    eft_raster_path_list = [
        raster_path
        for raster_pattern in args.raster_pattern
        for raster_path in glob.glob(raster_pattern)]
    kernel_lookup_by_n_pixels = {}
    for eft_raster_path in eft_raster_path_list:
        # get unique values
        unique_task = task_graph.add_task(
            func=geoprocessing.get_unique_values,
            args=((eft_raster_path, 1),),
            store_result=True,
            task_name=f'unique values for {eft_raster_path}')

        local_working_dir = os.path.join(
            workspace_dir, _basename(eft_raster_path))
        os.makedirs(local_working_dir, exist_ok=True)

        # build kernel for each search radius
        kernel_path_list = []
        for search_radius in args.search_radius:
            raster_info = geoprocessing.get_raster_info(eft_raster_path)
            min_raster_size = raster_info['pixel_size'][0]*min(
                *raster_info['raster_size'])
            if not args.force and (search_radius > .05*min_raster_size):
                raise ValueError(
                    f'the search radius {search_radius} is quite large '
                    f'compared to the raster size {min_raster_size}. I think '
                    f'this means that you put in meters when the raster is '
                    f'in WGS84 projection or some other issue. Either put '
                    f'radius in same units as raster or pass the --force '
                    f'flag if you think this was a mistake.')
            n_pixels = round(search_radius/raster_info['pixel_size'][0])
            kernel_path = os.path.join(
                local_working_dir, f'kernel_{n_pixels}.tif')
            if n_pixels not in kernel_lookup_by_n_pixels:
                kernel_task = task_graph.add_task(
                    func=_make_radius_kernel,
                    args=(n_pixels, kernel_path),
                    target_path_list=[kernel_path],
                    task_name=f'make kernel {n_pixels} {kernel_path}')
                kernel_lookup_by_n_pixels[n_pixels] = (
                    kernel_path, kernel_task)
            kernel_path = kernel_lookup_by_n_pixels[n_pixels][0]
            kernel_path_list.append((search_radius, kernel_path, kernel_task))

        # save relevant information so all kernels and unique tasks can be
        # scheduled before masking and convolution starts
        raster_stats_list.append(
            (eft_raster_path, unique_task, kernel_path_list,
             local_working_dir))

    for eft_raster_path, unique_task, kernel_path_list, local_working_dir in \
            raster_stats_list:
        unique_set = unique_task.get()
        for unique_value in unique_set:
            # create a mask of that value/raster
            mask_raster_path = os.path.join(
                local_working_dir, f'{unique_value}.tif')
            mask_task = task_graph.add_task(
                transient_run=True,
                func=_mask_by_value,
                args=(eft_raster_path, unique_value, mask_raster_path),
                target_path_list=[mask_raster_path],
                task_name=f'mask {mask_raster_path}')

            # do a convolution on that mask for the distance desired
            convolution_raster_list = []
            convolution_task_list = []
            for search_radius, kernel_path, kernel_task in kernel_path_list:
                convolution_raster_path = os.path.join(
                    local_working_dir, f'''{_basename(mask_raster_path)[0]}_{
                    unique_value}_{search_radius}.tif''')
                convolution_task = task_graph.add_task(
                    func=geoprocessing.convolve_2d,
                    args=((mask_raster_path, 1), (kernel_path, 1),
                          convolution_raster_path),
                    kwargs={
                        'working_dir': local_working_dir,
                        'target_datatype': gdal.GDT_Int32,
                        'target_nodata': NODATA,
                        },
                    dependent_task_list=[mask_task, kernel_task],
                    target_path_list=[convolution_raster_path],
                    task_name=f'convolve {convolution_raster_path}'
                    )
                convolution_task_list.append(convolution_task)
                convolution_raster_list.append(convolution_raster_path)

        # add all thresholds together, anything > 1 counts as 1
        eft_count_raster_path = os.path.join(
            workspace_dir, f'eft_count_{_basename(eft_raster_path)}.tif')
        task_graph.add_task(
            func=_add_thresholds,
            args=(convolution_raster_list, eft_count_raster_path),
            target_path_list=[eft_count_raster_path],
            dependent_task_list=convolution_task_list,
            task_name=f'adding thresholds for {eft_count_raster_path}')

    task_graph.join()
    task_graph.close()
    task_graph = None


if __name__ == '__main__':
    main()
