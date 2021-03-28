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
    wrflist = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files]

    rainc = w.getvar(wrflist, 'RAINC', timeidx=w.ALL_TIMES, method='cat')
    rainnc = w.getvar(wrflist, 'RAINNC', timeidx=w.ALL_TIMES, method='cat')
    total_rain = rainc + rainnc

    prec = total_rain.diff('Time', 1)#.sel(Time=pd.date_range('2017-06-01 3:00:00', '2017-06-8 00:00:00', freq='3H'))
    # prec = total_rain.isel(Time=-1)
    lats, lons = w.latlon_coords(prec)
    time = total_rain.Time.to_index() 

    return prec, lats, lons, time 



#%%
if __name__ == '__main__':
    data_dir1 = '/home/zzhzhao/Model/wrfout/test-10'
    data_dir2 = '/home/zzhzhao/Model/wrfout/test-10-removelake'
    prec1, lats, lons, time = load_wrfdata(data_dir1)
    prec2, lats, lons, time = load_wrfdata(data_dir2) 

    lat_range = (28, 34)
    lon_range = (86, 94)

    ### CMFD资料
    file_path = '/home/zzhzhao/code/python/wrf-test-10/data/prec_CMFD_201706.nc'
    cmfd = xr.open_dataset(file_path)['prec'].sel(lat=slice(lat_range[0],lat_range[1]), lon=slice(lon_range[0],lon_range[1])) * 3
    lat, lon =cmfd.lat, cmfd.lon

    ### 累计降水
    # prec_sum = prec.sel(Time=second_period).sum(dim='Time')
    cmfd_sum = cmfd.sum(dim='time')
    prec1_sum = prec1.sum(dim='Time')
    prec2_sum = prec2.sum(dim='Time')


#%%
    ### 累计降水分布
    proj = ccrs.PlateCarree()
    # crange = np.arange(0, 200+10, 10)
    labels = ['WRF', 'WRF-nolake', 'CMFD', 'Difference']
    fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(9,11), subplot_kw={'projection':proj})
    fig.subplots_adjust(hspace=0.2, wspace=0.15)
    for i in range(2):
        for j in range(2):
            set_grid(axes[i, j], lat=[30, 31.5], lon=[90, 91.5], span=.5)
            add_NamCo(axes[i, j])
    for j, prec_sum in enumerate([prec1_sum, prec2_sum]):
        c = axes[0][j].pcolor(lons, lats, prec_sum, vmin=0, vmax=250, cmap=cmaps.WhiteBlueGreenYellowRed, transform=proj)
        axes[0][j].set_title(labels[j], fontsize=14, weight='bold')
    axes[1][0].pcolor(lon, lat, cmfd_sum, vmin=0, vmax=250, cmap=cmaps.WhiteBlueGreenYellowRed, transform=proj)
    axes[1][0].set_title(labels[2], fontsize=14, weight='bold')
    
    c2 = axes[1][1].pcolor(lons, lats, prec1_sum-prec2_sum, vmin=-60, vmax=60, cmap='RdBu', transform=proj)
    axes[1][1].set_title(labels[3], fontsize=14, weight='bold')
    cb = fig.colorbar(c, ax=axes, orientation='horizontal', pad=0.05, shrink=0.9, aspect=35)
    cb.set_label('Precipitation / mm', fontsize=14)
    
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    axins = inset_axes(axes[1][1],
                  width="5%", # width = 10% of parent_bbox width
                  height="100%", # height : 50%
                  loc=6,
                  bbox_to_anchor=(1.05, 0., 1, 1),
                  bbox_transform=axes[1][1].transAxes,
                  borderpad=0,
              )
    
    cb2 = fig.colorbar(c2, cax=axins)#, orientation='vertical', shrink=0.6, aspect=25)
    # cb2.set_label('Precipitation / mm', fontsize=14)
    # axes[0][1].set_visible(False)

    fig.savefig('/home/zzhzhao/code/python/wrf-test-10/fig/prec.jpg', dpi=300)

