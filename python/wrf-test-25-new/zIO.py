import wrf as w
import numpy as np
import xarray as xr
import salem
from netCDF4 import Dataset
import pandas as pd 
import os 
import warnings
warnings.filterwarnings("ignore")
from zProcess import load_NamCo_shp

def load_wrflist(data_dir, domain):
    wrf_files = [f for f in os.listdir(data_dir) if f[11]==f'{domain}']
    wrflist = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files] 
    return wrflist

def load_t(data_dir, domain, t_name='TSK'): # t_name: 'TSK' or 'T2' or 'T_LAKE3D'
    wrf_files = [f for f in os.listdir(data_dir) if f[11]==f'{domain}']
    wrflist = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files] # 

    t = w.getvar(wrflist, t_name, timeidx=w.ALL_TIMES, method='cat')

    lats, lons = w.latlon_coords(t)
    time = t.Time.to_index() 
    return t, lats, lons, time 

def load_prec(data_dir, domain):
    wrf_files = [f for f in os.listdir(data_dir) if f[11]==f'{domain}']
    wrflist = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files]

    rainc = w.getvar(wrflist, 'RAINC', timeidx=w.ALL_TIMES, method='cat')
    rainnc = w.getvar(wrflist, 'RAINNC', timeidx=w.ALL_TIMES, method='cat')
    total_rain = rainc + rainnc

    prec = total_rain.diff('Time', 1)
    lats, lons = w.latlon_coords(prec)
    time = total_rain.Time.to_index() 
    return prec, lats, lons, time

def load_wind(data_dir, domain):
    wrf_files = [f for f in os.listdir(data_dir) if f[11]==f'{domain}']
    wrf_list = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files] 


    var_list = ['U10', 'V10']
    u10, v10 = [w.getvar(wrf_list, var, timeidx=w.ALL_TIMES, method='cat') for var in var_list]

    lats, lons = w.latlon_coords(u10)
    time = u10.Time.to_index() 
    return u10, v10, lats, lons, time

def load_cmfd(date_start, date_end, lat_range, lon_range):
    file_path = "/home/zzhzhao/data/CMFD_Prec_TP/CMFD_Prec_TP_1979-2018.nc"
    cmfd = xr.open_dataset(file_path)['prec'].sel(time=slice(date_start,date_end)) * 3
    cmfd = cmfd.sel(lat=slice(lat_range[0],lat_range[1]), lon=slice(lon_range[0],lon_range[1])) 
    return cmfd

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
    modis_lake_qc = modis_lake.where(modis_lake.count(dim=['lon','lat']) >= total_lakegrid/3.) # 有效湖泊格点数>=湖泊总格点数的1/2
    modis_lake_qc = modis_lake_qc.where(qc_lake<127) # 63:error<1K; 127:error<2K

    ### 湖面空间平均
    # modis_lake_mean = modis_lake.mean(dim=['lon', 'lat'])
    modis_lake_mean = modis_lake_qc.mean(dim=['lon', 'lat'])
    return modis_lake_mean