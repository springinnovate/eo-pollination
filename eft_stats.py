"""Get raster of EFT types."""
import argparse
import datetime
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
logging.getLogger('ecoshard.taskgraph').setLevel(logging.WARN)
LOGGER = logging.getLogger(__name__)


def _mask_by_value(base_raster_path, mask_value, target_raster_path):
    """Mask base by value to target."""
    def _mask(array, value):
        result = (array == value)
        return result

    geoprocessing.raster_calculator(
        [(base_raster_path, 1)], _mask, target_raster_path, gdal.GDT_Byte,
        None)


def main():
    """Entry point."""
    parser = argparse.ArgumentParser(
        description='EFT counter')
    parser.add_argument(
        'raster_pattern', type=str, nargs='+',
        help='path to raster or wildcards to process')
    args = parser.parse_args()

    workspace_dir = '_eft_stats_workspace'
    os.makedirs(workspace_dir, exist_ok=True)
    task_graph = taskgraph.TaskGraph(
        workspace_dir, multiprocessing.cpu_count(), 10.0)

    stats_path = (
        f'eft_stats_{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}.csv')
    raster_stats_list = []
    for raster_pattern in args.raster_pattern:
        for raster_path in glob.glob(raster_pattern):
            task = task_graph.add_task(
                func=geoprocessing.get_unique_values,
                args=((raster_path, 1),),
                store_result=True,
                task_name=f'unique values for {raster_path}')
            raster_stats_list.append((raster_path, task))

    for raster_path, unique_task in raster_stats_list:
        unique_set = unique_task.get()
        local_working_dir = os.path.join(
            workspace_dir, os.path.basename(os.path.splitext(raster_path)))
        for unique_value in unique_set:
            # create a mask of that value/raster
            mask_raster_path = os.path.join(
                local_working_dir, f'{unique_value}.tif')
            task_graph.add_task(
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
