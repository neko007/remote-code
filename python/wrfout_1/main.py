'''
对test-5的ERA5 reanalysis和降尺度wrfout输出进行对比
'''
#%%
import numpy as np
import xarray as xr 
import pandas as pd
from netCDF4 import Dataset 
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs 
from cartopy.feature import NaturalEarthFeature
import cmaps
import wrf as w
import os 
import xesmf as xe
### 隐藏警告
import warnings
warnings.filterwarnings("ignore")

#%%

if __name__ == '__main__':
    # 再分析资料
    re_dir = '/home/zzhzhao/Model/tests/test-5/data'
    wrf_dir = '/home/zzhzhao/Model/tests/test-5/WRF/run'

    re_file = 'ERA5-20200801-20200802-pl.grib'
    ds_re = xr.open_dataset(os.path.join(re_dir, re_file), engine='cfgrib')
    z_500_re = ds_re['z'].sel(isobaricInhPa=500) / 9.80665
    lat1d, lon1d = z_500_re.latitude, z_500_re.longitude

    # wrfout数据
    start_date = '20200801 00:00:00'
    end_date = '20200802 18:00:00'
    date_range = pd.date_range(start_date, end_date, freq='6H')
    
    wrf_files = [os.path.join(wrf_dir, f'wrfout_d01_{date.year}-{date.month:0>2d}-{date.day:0>2d}_{date.hour:0>2d}:00:00') for date in date_range]
    ds_wrf = [Dataset(wrf_file) for wrf_file in wrf_files]
    z = w.getvar(ds_wrf, 'z', timeidx=w.ALL_TIMES)
    lat2d, lon2d = w.latlon_coords(z)

    p = w.getvar(ds_wrf, "pressure")
    z_500 = w.interplevel(z, p, '500').sel(Time=z_500_re.time)
#%%
    # 插值 era5 插到 wrf 网格上
    src_grid = {'lat':lat1d, 'lon':lon1d}
    dst_grid = {'lat':lat2d, 'lon':lon2d}
    regridder = xe.Regridder(src_grid, dst_grid, 'bilinear')
    z_500_re_interp = regridder(z_500_re)
    regridder.clean_weight_file()

#%%
    # 绘图比较
    diff = z_500 - z_500_re_interp
    var_z = [z_500_re_interp, z_500, diff]

    # proj = ccrs.PlateCarree()
    # extent = [90,130,25,50]
    proj = w.get_cartopy(z_500)
    crange = [np.arange(5550,6000,25)] * 2 + [np.arange(-130,130+20,20)]
    cmap = [cmaps.rainbow] * 2 + ['RdBu_r']
    titles = ['ERA5 Reanalysis', 'WRF output', 'WRF - ERA5']


    fig, axes = plt.subplots(nrows=8, ncols=3, figsize=(12,24), subplot_kw={'projection':proj})
    fig.subplots_adjust(hspace=0.12,) 

    states = NaturalEarthFeature(category="cultural", scale="50m",
                             facecolor="none",
                             name="admin_0_boundary_lines_land")

    ps = []
    for i in range(8):
        for j in range(3):
            axes[i][j].coastlines(resolution='50m', lw=.8)
            # axes[0].set_extent(extent, crs=proj)
            axes[i][j].add_feature(states, linewidth=.5, edgecolor="black")
            p = axes[i][j].contourf(lon2d, lat2d, var_z[j].isel(time=i), levels=crange[j], transform=ccrs.PlateCarree(), cmap=cmap[j])
            axes[i][j].gridlines(color='k', linestyle='dotted')
            axes[0][j].set_title(titles[j], loc='left', fontsize=16)

            axes[i][j].set_xticks([], crs=proj)
            axes[i][j].set_yticks([], crs=proj)
            axes[i][0].set_ylabel(z_500_re.time.to_index()[i].strftime('%m-%d %H:%M'), fontsize=14)
            ps.append(p)


    cb1 = plt.colorbar(ps[0], ax=axes[:, :2], orientation='horizontal', shrink=0.8, aspect=30, pad=0.02)
    cb2 = plt.colorbar(ps[-1], ax=axes[:, 2], orientation='horizontal', shrink=0.8, aspect=15, pad=0.02)
    # # cb.set_label('z')


    fig.savefig('comp.jpg', dpi=300, bbox_inches='tight', pad_inches=0.1)

# %%
