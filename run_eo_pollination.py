"""
EO Pollination test script.

docker run --rm -it -v "%CD%":/usr/local/workspace therealspring/inspring:latest eo_test.py
"""
import logging
import sys
import multiprocessing
import inspring
from osgeo import gdal

gdal.SetCacheMax(2**27)

logging.basicConfig(
    level=logging.DEBUG,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)

LOGGER = logging.getLogger(__name__)
logging.getLogger('taskgraph').setLevel(logging.INFO)

args = {
    'eft_raster_path': 'eo_pollination_data/EFT1-45_CR_md5_f98029e71fbfc0f9e823dfc27c3b323d.tif',
    'guild_table_path': 'eo_pollination_data/guildtable_protobee_take2.csv',
    'farm_vector_path': 'eo_pollination_data/coffee_iCAFE_2012_md5_8f3644b1353587ef9e30fb61c0386f72.gpkg',
    'farm_floral_eft_field': 'eft',
    'results_suffix': 'eo_modification_take2',
    'workspace_dir': 'eo_pollination_data/eo_workspace',
    'n_workers': multiprocessing.cpu_count(),
    'pixel_scale': 1.0
}


if __name__ == '__main__':
    inspring.eo_pollination.execute(args)