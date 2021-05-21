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

def load_wrfdata1(data_dir, domain):
    wrf_files = [f for f in os.listdir(data_dir) if f[11]==f'{domain}']
    wrflist = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files] # 

    tsk = w.getvar(wrflist, 'TSK', timeidx=w.ALL_TIMES, method='cat')

    lats, lons = w.latlon_coords(tsk)
    time = tsk.Time.to_index() 

    return tsk, lats, lons, time 

def load_wrfdata2(data_dir, domain):
    wrf_files = [f for f in os.listdir(data_dir) if f[11]==f'{domain}']
    wrflist = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files] # 

    tsk = w.getvar(wrflist, 'T2', timeidx=w.ALL_TIMES, method='cat')

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
    geo_path = f'/home/zzhzhao/Model/tests/test-25/WPS/geo_em.d{domain:0>2d}.nc'
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


if __name__ == '__main__':
    data_dir = '/home/zzhzhao/Model/wrfout'
    testname_list = [
        'test-25-pre',
        'test-25',
        'test-25-WY',
        'test-25-WY2',
        'test-25-WY3',
        'test-25-WY4',
        ]
    N_test = len(testname_list)

    t_diff_list = dict()
    tsk_list = dict()
    t2_list = dict()

    for testname in testname_list:
        data_path = os.path.join(data_dir, testname)
        domain = 1 
        tsk, lats, lons, time = load_wrfdata1(data_path, domain)
        t2, lats, lons, time = load_wrfdata2(data_path, domain)
        tsk = xr.where(tsk>0, tsk, np.nan)
        t2 = xr.where(t2>0, t2, np.nan)
        # t_diff = t2 - tsk       

        mask = mask_lake(data_path, load_NamCo_shp(), testname, domain)
        # t_diff_NamCo = t_diff.where(mask) # 切出NamCo范围
        # t_diff_NamCo_mean = t_diff_NamCo.mean(dim=['west_east','south_north']).resample(Time='D').mean()

        # t_diff_list[testname] = t_diff_NamCo_mean

        t2_NamCo = t2.where(mask) # 切出NamCo范围
        t2_NamCo_mean = t2_NamCo.mean(dim=['west_east','south_north'])#.resample(Time='D').mean()
        tsk_NamCo = tsk.where(mask) # 切出NamCo范围
        tsk_NamCo_mean = tsk_NamCo.mean(dim=['west_east','south_north'])#.resample(Time='D').mean()
        t2_list[testname] = t2_NamCo_mean
        tsk_list[testname] = tsk_NamCo_mean

        ### 取3 UTC和6 UTC的平均
        # hour_list = [3, 6]
        # tsk_NamCo_daily = tsk_NamCo_mean.sel(Time=tsk_NamCo_mean.Time.dt.hour.isin(hour_list)).resample(Time='D').mean()
        # tsk_NamCo_daily_list[testname] = tsk_NamCo_daily


#%%
    labels = [
        'CTL_285K',
        'CTL',
        'WY',
        'WY2',
        'WY3',
        'WY4',
        ]

    ylen = np.ceil(np.sqrt(N_test)).astype(int); xlen = np.ceil(N_test/ylen).astype(int)
    default_len = 5

    fig = plt.figure(figsize=(xlen*default_len, ylen*default_len), dpi=100)
    fig.subplots_adjust(hspace=0.3, wspace=0.2)
    for i, testname in enumerate(testname_list):
        ax = fig.add_subplot(ylen, xlen, i+1)

        var1 = t2_list[testname]
        var2 = tsk_list[testname]
        var1.plot.line('r', mfc='none', label='T2', ax=ax)
        var2.plot.line('b', mfc='none', label='TSK', ax=ax)
        ax.set_title(labels[i])

        import matplotlib.dates as mdate  
        ax.xaxis.set_major_formatter(mdate.DateFormatter('%m-%d'))


    