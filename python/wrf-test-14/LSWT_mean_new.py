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
    data_dir = '/home/zzhzhao/Model/wrfout'
    testname_list = [
        'modis',
        'test-14',
        'test-14-oriLD',
        'test-15',
        'test-15-oriLD',
        'test-17',
        'test-18',
        ]
    N_test = len(testname_list)

    tsk_list = dict()
    tsk_NamCo_daily_list = dict()
    for testname in testname_list:
        if testname == 'modis':
            file_path = "/home/Public_Data/MODIS/MOD11A1/MOD11A1_NamCo_2017.nc"
            tsk_list[testname] = load_modis(file_path)
            tsk_NamCo_daily_list[testname] = tsk_list[testname]
        else:
            data_path = os.path.join(data_dir, testname)
            tsk, lats, lons, time = load_wrfdata(data_path)
            tsk = xr.where(tsk>0, tsk, np.nan)
            tsk_list[testname] = tsk

            mask = mask_lake(data_path, load_NamCo_shp())
            tsk_NamCo = tsk.where(mask) # 切出NamCo范围
            tsk_NamCo_mean = tsk_NamCo.mean(dim=['west_east','south_north'])
            # tsk_NamCo_mean = tsk_NamCo.isel(west_east=74, south_north=39)
            

            ### 取3 UTC和6 UTC的平均
            hour_list = [3, 6]
            tsk_NamCo_daily = tsk_NamCo_mean.sel(Time=tsk_NamCo_mean.Time.dt.hour.isin(hour_list)).resample(Time='D').mean()
            tsk_NamCo_daily_list[testname] = tsk_NamCo_daily


#%%
    labels = [
        'Modis',
        'Wuyang_90m', 
        'Wuyang_0.5m',
        'Default_50m', 
        'Default_0.5m',
        'Wuyang_50m', 
        'Wuyang_20m'
        ]
    markers = list('P^.sxD+')
    fig, ax = plt.subplots(dpi=100)
    for i, testname in enumerate(testname_list):
        var = tsk_NamCo_daily_list[testname]

        var.plot.line(lw=0, marker=markers[i], mfc='none', label=labels[i], ax=ax)
        ax.legend(loc=2, bbox_to_anchor=(1.0,1.0), borderaxespad=0, frameon=False)
        import matplotlib.dates as mdate  
        ax.xaxis.set_major_formatter(mdate.DateFormatter('%m-%d'))

    fig.savefig('fig/lwst_alltest.jpg', dpi=300, bbox_inches='tight', pad_inches=0.1)

    