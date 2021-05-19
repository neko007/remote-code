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
from zMap import set_grid, add_NamCo
import warnings
warnings.filterwarnings("ignore")

def load_wrfdata(data_dir):
    wrf_files = [f for f in os.listdir(data_dir) if f[11]=='2']
    wrflist = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files]

    rainc = w.getvar(wrflist, 'RAINC', timeidx=w.ALL_TIMES, method='cat')
    rainnc = w.getvar(wrflist, 'RAINNC', timeidx=w.ALL_TIMES, method='cat')
    total_rain = rainc + rainnc

    prec = total_rain.diff('Time', 1)
    lats, lons = w.latlon_coords(prec)
    time = total_rain.Time.to_index() 

    return prec, lats, lons, time, wrflist

if __name__ == '__main__':
    data_dir = '/home/zzhzhao/Model/wrfout'
    testname_list = [
        'cmfd',
        # 'obs',
        'test-14',
        # 'test-14-oriLD',
        'test-19',
        # 'test-15',
        # 'test-15-oriLD',
        # 'test-17',
        # 'test-18',
        'test-20',
        'test-23',
        'test-24',
        ]
    
    N_test = len(testname_list)
    # NamCo_latlon = (30.75, 90.8)
    lat_range = (30.5, 31)
    lon_range = (90.2, 91)
    date_start = '2017-06-01'
    date_end = '2017-06-30'

    prec_mean_list = dict()
    for testname in testname_list:
        if testname == 'cmfd':
            file_path = "/home/zzhzhao/code/python/wrf-test-10/data/prec_CMFD_201706.nc"
            cmfd = xr.open_dataset(file_path)['prec']
            cmfd_mean = cmfd.sel(lat=slice(lat_range[0],lat_range[1]), lon=slice(lon_range[0],lon_range[1])).mean(dim=['lat', 'lon']).resample(time='1D').sum() * 3
            prec_mean_list[testname] = cmfd_mean
        elif testname == 'obs':
            df = pd.read_excel('data/纳木错站2017-2018.xlsx', index_col=0)
            obs_sta = df.loc[pd.date_range(date_start, date_end)]['降水量']
            prec_mean_list[testname] = obs_sta
        else:
            data_path = os.path.join(data_dir, testname)
            prec, lats, lons, time, wrflist = load_wrfdata(data_path)
            prec = xr.where(prec>0, prec, np.nan)
            ### 纳木错站点插值
            left_bottom = w.ll_to_xy(wrflist, lat_range[0], lon_range[0]).values
            right_top = w.ll_to_xy(wrflist, lat_range[1], lon_range[1]).values
            prec_mean = prec.sel(south_north=slice(left_bottom[1], right_top[1]), west_east=slice(left_bottom[0], right_top[0])).mean(dim=['south_north', 'west_east']).resample(Time='1D').sum()
            prec_mean_list[testname] = prec_mean


    labels = [
        'CMFD',
        # 'StaObs',
        'Wuyang_90m', 
        # 'Wuyang_0.5m',
        'Default_90m',
        # 'Default_50m', 
        # 'Default_0.5m',
        # 'Wuyang_50m', 
        # 'Wuyang_20m',
        'Wuyang_90m_Update',
        'Wuyang_90m_279.5K',
        'Default_90m_277K',
        ]
    markers = list('PX^.sxD+*p')
    fig, ax = plt.subplots(dpi=200)
    for i, testname in enumerate(testname_list):
        var = prec_mean_list[testname]

        # var.plot.line(lw=0, marker=markers[i], mfc='none', label=labels[i], ax=ax)
        var.plot.line(lw=1.2, mfc='none', label=labels[i], ax=ax)
        ax.legend(loc=2, bbox_to_anchor=(1.0,1.0), borderaxespad=0, frameon=False)
        import matplotlib.dates as mdate  
        ax.xaxis.set_major_formatter(mdate.DateFormatter('%m-%d'))

    # fig.savefig('fig/prec_sta_alltest.jpg', dpi=300, bbox_inches='tight', pad_inches=0.1)