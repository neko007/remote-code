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
from Module import load_cmfd
import warnings
warnings.filterwarnings("ignore")

def load_wrfdata(data_dir, domain):
    wrf_files = [f for f in os.listdir(data_dir) if f[11]==f'{domain}']
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
        # 'test-25-pre',
        # 'test-25',
        # 'test-25-WY',
        'test-25-WY4',
        ]
    
    N_test = len(testname_list)
    # NamCo_latlon = (30.75, 90.8)
    lat_range = (30.5, 31)
    lon_range = (90, 91)
    date_start = '2013-08-23'
    date_end = '2013-09-01'

    prec_mean_list = dict()
    for testname in testname_list:
        if testname == 'cmfd':
            cmfd = load_cmfd(date_start, date_end, lat_range, lon_range)
            prec_mean_list[testname] = cmfd.mean(dim=['lat', 'lon']).resample(time='1D').sum()
        elif testname == 'obs':
            df = pd.read_excel('data/纳木错站2017-2018.xlsx', index_col=0)
            obs_sta = df.loc[pd.date_range(date_start, date_end)]['降水量']
            prec_mean_list[testname] = obs_sta
        else:
            data_path = os.path.join(data_dir, testname)
            domain = 1
            prec, lats, lons, time, wrflist = load_wrfdata(data_path, domain)
            prec = xr.where(prec>0, prec, np.nan)
            ### 纳木错站点插值
            left_bottom = w.ll_to_xy(wrflist, lat_range[0], lon_range[0]).values
            right_top = w.ll_to_xy(wrflist, lat_range[1], lon_range[1]).values
            prec_mean = prec.sel(south_north=slice(left_bottom[1], right_top[1]), west_east=slice(left_bottom[0], right_top[0])).mean(dim=['south_north', 'west_east']).resample(Time='1D').sum()
            prec_mean_list[testname] = prec_mean


    labels = [
        'CMFD',
        # 'CTL_285K',
        # 'CTL',
        # 'WY',
        'WY4',
        ]
    markers = list('PX^.sxD+*p')
    fig, ax = plt.subplots(dpi=200)
    for i, testname in enumerate(testname_list):
        var = prec_mean_list[testname]

        # var.plot.line(lw=0, marker=markers[i], mfc='none', label=labels[i], ax=ax)
        var.plot.step(lw=1.2, mfc='none', label=labels[i], ax=ax)
        ax.legend(loc=2, bbox_to_anchor=(1.0,1.0), borderaxespad=0, frameon=False)
        import matplotlib.dates as mdate  
        ax.xaxis.set_major_formatter(mdate.DateFormatter('%m-%d'))

    # fig.savefig('fig/prec_sta_alltest.jpg', dpi=300, bbox_inches='tight', pad_inches=0.1)