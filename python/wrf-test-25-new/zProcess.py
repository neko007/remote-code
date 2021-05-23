import geopandas
import salem
import xarray as xr
import os
import warnings
warnings.filterwarnings("ignore")

def load_NamCo_shp():
    f_in = '/home/zzhzhao/code/shpfiles/boundary/hyd1_4m/hyd1_4p.shp'
    shp = geopandas.read_file(f_in, encoding='gbk')
    shp = shp.loc[shp['NAME'].isin(["纳木错"])]
    return shp

def mask_lake(data_dir, shp, testname, domain):
    # geo_path = f'/home/zzhzhao/Model/tests/{testname}/WPS/geo_em.d{domain:0>2d}.nc'
    geo_path = f'/home/zzhzhao/Model/tests/test-25/WPS/geo_em.d{domain:0>2d}.nc'
    lu = salem.open_wrf_dataset(os.path.join(data_dir, geo_path))['LU_INDEX'].isel(time=0)
    lu_lake = lu.salem.roi(shape=shp)
    mask = xr.where(lu_lake.notnull(), True, False)
    return mask