#%%
import wrf as w
import xarray as xr
from netCDF4 import Dataset
import pandas as pd 
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs
import cmaps
import os 
import sys
sys.path.append('/home/zzhzhao/code/python/wrf-test-10')
from zMap import set_grid, add_NamCo
import warnings
warnings.filterwarnings("ignore")

def load_wrfdata(data_dir):
    wrf_files = [f for f in os.listdir(data_dir) if f[9]=='2']
    wrf_list = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files] 

    var_list = ['U10', 'V10']
    u10, v10 = [w.getvar(wrf_list, var, timeidx=w.ALL_TIMES, method='cat') for var in var_list]

    lats, lons = w.latlon_coords(u10)

    return u10, v10, lats, lons

if __name__ == '__main__':
    data_dir1 = '/home/zzhzhao/Model/wrfout/test-10'
    data_dir2 = '/home/zzhzhao/Model/wrfout/test-10-removelake'
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
        axes[0][i].set_title(labels[i], fontsize=14, weight='bold')

        q = axes[0][i].quiver(lons.values, lats.values, eval(f'u10{i+1}_mean.values'), eval(f'v10{i+1}_mean.values'), color='b', transform=ccrs.PlateCarree(), scale=35, regrid_shape=15, width=0.0035, headwidth=5)
        axes[0][i].quiverkey(q, X=0.85, Y=1.05, U=2, label='2 m/s', labelpos='E', fontproperties={'size':12})

    set_grid(axes[1][0], lat=[30, 31.5], lon=[90, 91.5], span=.5)
    add_NamCo(axes[1][0])
    axes[1][0].set_title(labels[2], fontsize=14, weight='bold')

    q = axes[1][0].quiver(lon.values, lat.values, u10_era5.values, v10_era5.values, color='b', transform=ccrs.PlateCarree(), scale=35, regrid_shape=15, width=0.0035, headwidth=5)
    axes[1][0].quiverkey(q, X=0.85, Y=1.05, U=2, label='2 m/s', labelpos='E', fontproperties={'size':12})

    axes[1][1].set_visible(False)
    
    fig.savefig('/home/zzhzhao/code/python/wrf-test-10/fig/wind.jpg', dpi=300)
        