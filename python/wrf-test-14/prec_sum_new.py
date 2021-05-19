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
import warnings
warnings.filterwarnings("ignore")

def load_wrfdata(data_dir):
    wrf_files = [f for f in os.listdir(data_dir) if f[11]=='2']
    wrflist = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files]

    rainc = w.getvar(wrflist, 'RAINC', timeidx=w.ALL_TIMES, method='cat')
    rainnc = w.getvar(wrflist, 'RAINNC', timeidx=w.ALL_TIMES, method='cat')
    total_rain = rainc + rainnc

    prec = total_rain.diff('Time', 1)
    # prec = total_rain.isel(Time=-1)
    lats, lons = w.latlon_coords(prec)
    time = total_rain.Time.to_index() 

    return prec, lats, lons, time 

#%%
if __name__ == '__main__':
    data_dir = '/home/zzhzhao/Model/wrfout'
    testname_list = [
        'cmfd',
        'test-14',
        # 'test-14-oriLD',
        'test-19',
        'test-15',
        # 'test-15-oriLD',
        'test-17',
        'test-18',
        'test-20',
        'test-24-pre',
        'test-24',
        ]
    
    N_test = len(testname_list)
    lat_range = (28, 34)
    lon_range = (86, 94)

    prec_list = dict()
    prec_sum_list = dict()
    for testname in testname_list:
        if testname == 'cmfd':
            file_path = "/home/zzhzhao/code/python/wrf-test-10/data/prec_CMFD_201706.nc"
            cmfd = xr.open_dataset(file_path)['prec'].sel(lat=slice(lat_range[0],lat_range[1]), lon=slice(lon_range[0],lon_range[1])) * 3
            prec_list[testname] = cmfd
            prec_sum_list[testname] = cmfd.sum(dim='time')
        else:
            data_path = os.path.join(data_dir, testname)
            prec, lats, lons, time = load_wrfdata(data_path)
            prec = xr.where(prec>0, prec, np.nan)
            prec_list[testname] = prec
            prec_sum_list[testname] = prec.sum(dim='Time')

#%%
    ### 累计降水分布
    labels = [
        'CMFD',
        'Wuyang_90m', 
        # 'Wuyang_0.5m',
        'Default_90m',
        'Default_50m', 
        # 'Default_0.5m',
        'Wuyang_50m', 
        'Wuyang_20m',
        'Wuyang_90m_Update',
        'Default_90m_279.5K',
        'Default_90m_277K',
        ]
    proj = ccrs.PlateCarree()
    cmap = cmaps.WhiteBlueGreenYellowRed
    ylen = np.ceil(np.sqrt(N_test)).astype(int); xlen = np.ceil(N_test/ylen).astype(int)
    default_len = 5

    fig = plt.figure(figsize=(xlen*default_len, ylen*default_len), dpi=100)
    fig.subplots_adjust(hspace=0.3, wspace=0.2)
    axes = []
    for i, testname in enumerate(testname_list):
        ax = fig.add_subplot(ylen, xlen, i+1, projection=proj)
        var = prec_sum_list[testname]
        if testname == 'cmfd':
            lat, lon = prec_sum_list[testname].lat, prec_sum_list[testname].lon 
        else:
            lat, lon = lats, lons
    
        set_grid(ax, lat=[30, 31.5], lon=[90, 91.5], span=.5)
        # set_grid(ax, lat=[29, 32], lon=[89, 92], span=.5)
        add_NamCo(ax)
        c = ax.pcolor(lon, lat, var, cmap=cmap, vmin=0, vmax=250, transform=proj)
        ax.set_title(labels[i])
     
        # axes.append(ax)
        
    # cb = fig.colorbar(c, ax=np.array(axes).all(), orientation='horizontal', pad=0.05, shrink=0.9, aspect=35)
    # cb.set_label('Precipitation / mm', fontsize=14)

