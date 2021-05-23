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
import salem 
import warnings
warnings.filterwarnings("ignore")

def load_wrfdata(data_dir, domain):
    wrf_files = [f for f in os.listdir(data_dir) if f[11]==f'{domain}']
    wrflist = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files] # 

    T_lake3d = w.getvar(wrflist, 'T_LAKE3D', timeidx=w.ALL_TIMES, method='cat')

    lats, lons = w.latlon_coords(T_lake3d)
    time = T_lake3d.Time.to_index() 

    return T_lake3d, lats, lons, time, wrflist

def load_z(testname):
    version = 'WRFV3' if 'WY' in testname else 'WRF'
    if testname == 'test-25-pre': testname = 'test-25'
    z_lake3d = w.getvar(Dataset(f'/home/zzhzhao/Model/tests/{testname}/{version}/run/wrfout_d01_2013-08-23_00:00:00'), 'Z_LAKE3D', timeidx=0)
    return z_lake3d

def load_sta():
    ds = xr.open_dataset('/home/Public_Data/NamCo_Station_Data/NamCo_2013_LSWT.nc')
    return ds.TLake

if __name__ == '__main__':
    data_dir = '/home/zzhzhao/Model/wrfout'
    testname_list = [
        'sta',
        'test-25',
        'test-25-pre',
        # 'test-25-WY',
        # 'test-25-WY2',
        'test-25-WY3',
        'test-25-WY4',
        'test-25-3',
        ]
    N_test = len(testname_list)

    lat = 30.75
    lon = 90.75

    tl_list = dict()
    zl_list = dict()

    for testname in testname_list:
        if testname == 'sta':
            TLake = load_sta().sel(time=slice(20130823,20130901))
            tl_list[testname] = TLake
        else:
            data_path = os.path.join(data_dir, testname)
            domain = 1
            tl, lats, lons, time, wrflist = load_wrfdata(data_path, domain)
            zl = load_z(testname)

            tl = xr.where(tl>0, tl, np.nan)

            x, y = w.ll_to_xy(wrflist, lat, lon)

            tl_list[testname] = tl.sel(west_east=x, south_north=y)
            zl_list[testname] = zl.sel(west_east=x, south_north=y)

        
        
#%%
    crange = np.arange(276, 286+1, 1)
    cmap = 'rainbow'
    titles = [
        'Obs',
        'CTL',
        'CTL_285K',
        # 'WY',
        # 'WY2',
        'WY3',
        'WY4',
        'CTL3_285K',
        ]

    ylen = np.ceil(np.sqrt(N_test)).astype(int); xlen = np.ceil(N_test/ylen).astype(int)
    default_len = 5

    fig = plt.figure(figsize=(xlen*default_len, ylen*default_len), dpi=100)
    fig.subplots_adjust(hspace=0.3, wspace=0.2)
    for i, testname in enumerate(testname_list):
        ax = fig.add_subplot(ylen, xlen, i+1)
        var = tl_list[testname].T

        if testname == 'sta':
            t = pd.to_datetime(var.time)
            p = ax.contourf(t, var.lake_level, var+273.15, levels=crange, cmap=cmap)
        else:
        # var.plot.contourf(levels=crange, cmap=cmap, add_labels=False, ax=ax)
            var2 = var.resample(Time='D').mean()
            t = var2.Time
            p = ax.contourf(t, zl_list[testname].values, var2, levels=crange, cmap=cmap)
        plt.colorbar(p)

        ax.set_ylim([0, 50])
        ax.invert_yaxis()
        ax.set_title(titles[i])

        
        import matplotlib.dates as mdate  
        ax.xaxis.set_major_formatter(mdate.DateFormatter('%m-%d'))
    
    # fig.savefig('fig/Tlake3d_alltest.jpg', dpi=300, bbox_inches='tight', pad_inches=0.1)