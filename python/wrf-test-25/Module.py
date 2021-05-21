import wrf as w
import numpy as np
import xarray as xr
from netCDF4 import Dataset
import pandas as pd 
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs
import cmaps
import os 
import salem 
import warnings
warnings.filterwarnings("ignore")

def load_wrfdata(data_dir, domain):
    wrf_files = [f for f in os.listdir(data_dir) if f[11]==f'{domain}']
    wrflist = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files] # 

    tsk = w.getvar(wrflist, 'TSK', timeidx=w.ALL_TIMES, method='cat')

    # ### 强制附上salem的grid属性，才能施以salem独有的方法
    # tsk.attrs['pyproj_srs'] = salem.open_wrf_dataset(os.path.join(data_dir, wrf_files[0])).attrs['pyproj_srs']

    # del tsk.attrs['projection']
    # del tsk.attrs['coordinates']

    lats, lons = w.latlon_coords(tsk)
    time = tsk.Time.to_index() 

    return tsk, lats, lons, time 

def load_NamCo_shp():
    import geopandas
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

def load_modis(file_path, date_start, date_end, day_or_night):
    modis = salem.open_xr_dataset(file_path)
    time = modis.indexes['time'].to_datetimeindex()
    modis.coords['time'] = time
    modis = modis.sel(time=slice(date_start,date_end))
    modis_lake = modis[f'LST_{day_or_night}_1km'].salem.roi(shape=load_NamCo_shp()) 

    ### 质量控制
    qc = modis[f'QC_{day_or_night}'].sel(time=pd.date_range(date_start,date_end))
    qc_lake = qc.salem.roi(shape=load_NamCo_shp())
    total_lakegrid = qc_lake.isel(time=0).count().values # 湖泊格点总数
    modis_lake_qc = modis_lake.where(modis_lake.count(dim=['lon','lat']) >= total_lakegrid/10.) # 有效湖泊格点数>=湖泊总格点数的1/2
    # modis_lake_qc = modis_lake_qc.where(qc_lake<127) # 63:error<1K; 127:error<2K

    ### 湖面空间平均
    # modis_lake_mean = modis_lake.mean(dim=['lon', 'lat'])
    modis_lake_mean = modis_lake_qc.mean(dim=['lon', 'lat'])
    return modis_lake_mean

def load_cmfd(date_start, date_end, lat_range, lon_range):
    file_path = "/home/zzhzhao/data/CMFD_Prec_TP/CMFD_Prec_TP_1979-2018.nc"

    cmfd = xr.open_dataset(file_path)['prec'].sel(time=slice(date_start,date_end)) * 3
    
    cmfd = cmfd.sel(lat=slice(lat_range[0],lat_range[1]), lon=slice(lon_range[0],lon_range[1])) 

    return cmfd