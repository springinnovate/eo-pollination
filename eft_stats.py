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

gdal.SetCacheMax(2**26)
logging.basicConfig(
    level=logging.DEBUG,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(funcName)s:%(lineno)d] %(message)s'))
logging.getLogger('ecoshard.taskgraph').setLevel(logging.INFO)
LOGGER = logging.getLogger(__name__)


def _mask_by_value(base_raster_path, mask_value, target_raster_path):
    """Mask base by value to target."""
    geoprocessing.raster_calculator(
        [(base_raster_path, 1)], lambda a: (a == mask_value),
        target_raster_path, gdal.GDT_Byte, None)


def _make_radius_kernel(n_pixels, target_kernel_path):
    """Create a flat radius kernel at target kernel path."""
    pass


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

    stats_path = (
        f'eft_stats_{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}.csv')
    raster_stats_list = []
    raster_path_list = [
        raster_path
        for raster_pattern in args.raster_pattern
        for raster_path in glob.glob(raster_pattern)]
    for raster_path in raster_path_list:
        local_working_dir = os.path.join(
            workspace_dir, os.path.basename(
                os.path.splitext(raster_path)[0]))

        kernel_path_list = []
        for search_radius in args.search_radius:
            raster_info = geoprocessing.get_raster_info(raster_path)
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
            kernel_path = os.path.join(
                local_working_dir, f'kernel_{search_radius}.tif')
            n_pixels = round(search_radius/raster_info['pixel_size'])
            kernel_task = task_graph.add_task(
                func=_make_radius_kernel,
                args=(n_pixels, kernel_path),
                target_path_list=[kernel_path],
                task_name=f'make kernel {n_pixels} {kernel_path}')
            kernel_path_list.append((search_radius, kernel_path))

        unique_task = task_graph.add_task(
            func=geoprocessing.get_unique_values,
            args=((raster_path, 1),),
            store_result=True,
            task_name=f'unique values for {raster_path}')
        raster_stats_list.append((raster_path, unique_task, kernel_path_list))

    for raster_path, unique_task, kernel_path_list in raster_stats_list:
        unique_set = unique_task.get()
        local_working_dir = os.path.join(
            workspace_dir, os.path.basename(os.path.splitext(raster_path)[0]))
        for unique_value in unique_set:
            # create a mask of that value/raster
            mask_raster_path = os.path.join(
                local_working_dir, f'{unique_value}.tif')
            mask_task = task_graph.add_task(
                func=_mask_by_value,
                args=(raster_path, unique_value, mask_raster_path),
                target_path_list=[mask_raster_path],
                task_name=f'mask {mask_raster_path}')
            # do a convolution on that mask for the distance desired
            # threshold the convolution so >0 = 1
        # add all thresholds together

    task_graph.join()
    task_graph.close()
    task_graph = None


if __name__ == '__main__':
    main()
