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
    wrf_list = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files] 


    var_list = ['U10', 'V10']
    u10, v10 = [w.getvar(wrf_list, var, timeidx=w.ALL_TIMES, method='cat') for var in var_list]

    lats, lons = w.latlon_coords(u10)

    return u10, v10, lats, lons

def add_terrain(ax, geo_path):
    ds = xr.open_dataset(geo_path)
    hgt = ds.HGT_M.squeeze()
    lats = ds.XLAT_M.squeeze()
    lons = ds.XLONG_M.squeeze()

    crange = np.arange(3000, 7000+100, 100)
    # p = ax.pcolor(lons, lats, hgt, vmin=3000, vmax=7000, transform=ccrs.PlateCarree(), cmap='terrain')
    p = ax.contourf(lons, lats, hgt, levels=crange, extend='both', transform=ccrs.PlateCarree(), cmap='terrain')

#%%  
if __name__ == '__main__':
    # data_dir1 = '/home/zzhzhao/Model/wrfout/test-14'
    # data_dir2 = '/home/zzhzhao/Model/wrfout/test-14-nolake'
    data_dir1 = '/home/zzhzhao/Model/wrfout/test-14-oriLD'
    data_dir2 = '/home/zzhzhao/Model/wrfout/test-14-nolake-oriLD'
    geo_path = '/home/zzhzhao/Model/tests/test-14/WPS/geo_em.d02.nc'
    u101, v101, lats, lons = load_wrfdata(data_dir1)
    u102, v102, _, _ = load_wrfdata(data_dir2)

    u101_mean, v101_mean, u102_mean, v102_mean = [var.mean(dim='Time') for var in [u101, v101, u102, v102]]

    ### ERA5 （CMFD只有全风速）
    file_path = '/home/zzhzhao/code/python/wrf-test-10/data/wind_era5_201706_monthly.nc'
    u10_era5 = xr.open_dataset(file_path)['u10'].squeeze()
    v10_era5 = xr.open_dataset(file_path)['v10'].squeeze()
    lat, lon = u10_era5.latitude, u10_era5.longitude

#%%

    proj = ccrs.PlateCarree()
    labels = ['WRF', 'WRF-nolake', 'ERA5', 'Difference']

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10,10), dpi=100, subplot_kw={'projection':proj})
    for i in range(2):
        set_grid(axes[0][i], lat=[30, 31.5], lon=[90, 91.5], span=.5)
        add_NamCo(axes[0][i])
        add_terrain(axes[0][i], geo_path)
        axes[0][i].set_title(labels[i], fontsize=14, weight='bold')

        q = axes[0][i].quiver(lons.values, lats.values, eval(f'u10{i+1}_mean.values'), eval(f'v10{i+1}_mean.values'), color='b', transform=ccrs.PlateCarree(), scale=35, regrid_shape=15, width=0.0035, headwidth=5)
        axes[0][i].quiverkey(q, X=0.85, Y=1.05, U=2, label='2 m/s', labelpos='E', fontproperties={'size':12})

    set_grid(axes[1][0], lat=[30, 31.5], lon=[90, 91.5], span=.5)
    add_NamCo(axes[1][0])
    add_terrain(axes[1][0], geo_path)
    axes[1][0].set_title(labels[2], fontsize=14, weight='bold')

    q = axes[1][0].quiver(lon.values, lat.values, u10_era5.values, v10_era5.values, color='b', transform=ccrs.PlateCarree(), scale=35, regrid_shape=15, width=0.0035, headwidth=5)
    axes[1][0].quiverkey(q, X=0.85, Y=1.05, U=2, label='2 m/s', labelpos='E', fontproperties={'size':12})

    axes[1][1].set_visible(False)
    
    fig.savefig('/home/zzhzhao/code/python/wrf-test-14/fig/wind_oriLD.jpg', dpi=300)
        