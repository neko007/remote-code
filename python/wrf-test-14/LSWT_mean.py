#%%
import wrf as w
import numpy as np
import xarray as xr
from netCDF4 import Dataset
import pandas as pd 
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs
import cmaps
import os 
import sys
sys.path.append('/home/zzhzhao/code/python/wrf-test-14')
from zMap import set_grid, add_NamCo
import salem 
import warnings
warnings.filterwarnings("ignore")

def load_wrfdata(data_dir):
    wrf_files = [f for f in os.listdir(data_dir) if f[11]=='2']
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

def mask_lake(data_dir, shp):
    geo_path = '/home/zzhzhao/Model/tests/test-14/WPS/geo_em.d02.nc'
    lu = salem.open_wrf_dataset(os.path.join(data_dir, geo_path))['LU_INDEX'].isel(time=0)
    lu_lake = lu.salem.roi(shape=shp)
    mask = xr.where(lu_lake.notnull(), True, False)
    return mask

def load_modis(file_path):
    modis = salem.open_xr_dataset(file_path)
    time = modis.indexes['time'].to_datetimeindex()
    modis.coords['time'] = time
    modis = modis.sel(time=pd.date_range('2017-06-01','2017-06-30'))
    modis_lake = modis['LST_Day_1km'].salem.roi(shape=load_NamCo_shp()) 

    ### 质量控制
    qc = modis['QC_Day'].sel(time=pd.date_range('2017-06-01','2017-06-30'))
    qc_lake = qc.salem.roi(shape=load_NamCo_shp())
    total_lakegrid = qc_lake.isel(time=0).count().values # 湖泊格点总数
    modis_lake_qc = modis_lake.where(modis_lake.count(dim=['lon','lat']) >= total_lakegrid/2.) # 有效湖泊格点数>=湖泊总格点数的1/2
    modis_lake_qc = modis_lake_qc.where(qc_lake<63) # 63:error<1K; 127:error<2K

    ### 湖面空间平均
    # modis_lake_mean = modis_lake.mean(dim=['lon', 'lat'])
    modis_lake_mean = modis_lake_qc.mean(dim=['lon', 'lat'])
    return modis_lake_mean

#%%
if __name__ == '__main__':
    data_dir1 = '/home/zzhzhao/Model/wrfout/test-14-oriLD'
    data_dir2 = '/home/zzhzhao/Model/wrfout/test-15'
    # data_dir1 = '/home/zzhzhao/Model/wrfout/test-14-oriLD'
    # data_dir2 = '/home/zzhzhao/Model/wrfout/test-14-nolake-oriLD'
    tsk1, lats, lons, time = load_wrfdata(data_dir1)
    tsk2, lats, lons, time = load_wrfdata(data_dir2) 

    mask = mask_lake(data_dir1, load_NamCo_shp())
    tsk1_lake = tsk1.where(mask) # 切出NamCo范围
    tsk1_lake_mean = tsk1_lake.mean(dim=['west_east','south_north'])
    tsk2_lake = tsk2.where(mask) # 切出NamCo范围
    tsk2_lake_mean = tsk2_lake.mean(dim=['west_east','south_north'])

    ### MODIS
    file_path = "/home/Public_Data/MODIS/MOD11A1/MOD11A1_NamCo_2017.nc"
    modis_lake_mean = load_modis(file_path)

    ### 取3 UTC和6 UTC的平均
    hour_list = [3, 6]
    tsk1_lake_daily = tsk1_lake_mean.sel(Time=tsk1_lake_mean.Time.dt.hour.isin(hour_list)).resample(Time='D').mean()
    tsk2_lake_daily = tsk2_lake_mean.sel(Time=tsk2_lake_mean.Time.dt.hour.isin(hour_list)).resample(Time='D').mean()

    ### 取3 UTC
    # hour_list = [6]
    # tsk1_lake_daily = tsk1_lake_mean.sel(Time=tsk1_lake_mean.Time.dt.hour.isin(hour_list))
    # tsk2_lake_daily = tsk2_lake_mean.sel(Time=tsk2_lake_mean.Time.dt.hour.isin(hour_list))

    fig, ax = plt.subplots(dpi=100, figsize=(7,4))
    tsk1_lake_daily.plot.line('r^', label='Ctrl', ax=ax)
    tsk2_lake_daily.plot.line('bo', label='noLake', ax=ax)
    # tsk1_lake_mean.plot.line('r^', label='Ctrl', ax=ax)
    # tsk2_lake_mean.plot.line('bo', label='noLake', ax=ax)
    modis_lake_mean.plot.line('g+', label='Modis', ax=ax)
    ax.legend(frameon=False, ncol=3)
    ax.set_ylabel('LSWT / K')
    # ax.grid(alpha=0.5, ls='--')
    import matplotlib.dates as mdate    
    ax.xaxis.set_major_formatter(mdate.DateFormatter('%m-%d'))
    
    # fig.savefig('/home/zzhzhao/code/python/wrf-test-14/fig/lswt_oriLD.jpg', dpi=300, bbox_inches='tight', pad_inches=0.1)
    