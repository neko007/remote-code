#%%
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature import NaturalEarthFeature
from netCDF4 import Dataset
import wrf as w
import cmaps
import warnings
warnings.filterwarnings("ignore")

# plt.style.use(['science','ieee'])

def add_artist(ax, proj, lat=[20,55], lon=[70,135]):
    from cartopy.io.shapereader import Reader
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import matplotlib.ticker as mticker
    '''
    添加地图要素 → 海岸线，限定区域，网格，省界
    '''
    lw = 1
    ax.coastlines(resolution='50m', lw=lw)
    ax.set_extent([70, 130, 20, 55], crs=proj)
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), lw=lw)

    ### 限定区域
    extent = [np.min(lon), np.max(lon), np.min(lat), np.max(lat)]
    ax.set_extent(extent, crs=proj)

    ### 设置网格
    lat_span = 10
    lon_span = 15

    ax.set_xticks([], crs=ccrs.PlateCarree())
    ax.set_yticks([], crs=ccrs.PlateCarree())

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color=None, alpha=0, linestyle='--')
    gl.xlabels_top = False  # 关闭顶端的经纬度标签
    gl.ylabels_right = False  # 关闭右侧的经纬度标签
    gl.xformatter = LONGITUDE_FORMATTER  # x轴设为经度的格式
    gl.yformatter = LATITUDE_FORMATTER  # y轴设为纬度的格式
    gl.xlocator = mticker.FixedLocator(np.arange(extent[0], extent[1]+lon_span, lon_span))
    gl.ylocator = mticker.FixedLocator(np.arange(extent[2], extent[3]+lat_span, lat_span))
    gl.xlabel_style = {'size':12}
    gl.ylabel_style = {'size':12}
    return ax

def add_province(ax):
    '''
    添加省界
    '''
    import geopandas
    f_in = '../../boundary/province/Province_9.shp'
    shp = geopandas.read_file(f_in)
    ax.add_geometries(shp.geometry, crs=ccrs.PlateCarree(), edgecolor='k', alpha=0.5, facecolor='none', lw=0.8)
    return ax

def add_citylabel(ax, proj):
    city_names = ['YL', 'Urumqi', 'YC', 'TY', 'BJ']
    city_latlons = [(43.9404, 81.2815), (43.8303, 87.5801), (38.4975, 106.2328), (37.7124, 112.469), (40.0031, 116.407)]
    for i in range(len(city_names)):
        lat, lon = city_latlons[i]
        ax.plot(lon, lat, c='b', marker='.', transform=proj)
        ax.text(lon, lat+1, city_names[i], c='b', fontsize=12, weight=600, horizontalalignment='center')
    return ax

#%%

if __name__ == '__main__':
    # wrfout
    data_dir = '/home/zzhzhao/Model/tests/test-chem-6/WRF/run'
    date1 = '20190508 00:00:00'
    date2 = '20190519 18:00:00'
    date_range = pd.date_range(date1, date2, freq='1D')

    wrf_files = [f'{data_dir}/wrfout_d01_{date.year}-{date.month:0>2d}-{date.day:0>2d}_{date.hour:0>2d}:00:00' for date in date_range]

    wrflist = [Dataset(wrf_file) for wrf_file in wrf_files]

    u10_wrf, v10_wrf, pm10_wrf = [w.getvar(wrflist, var_name, timeidx=w.ALL_TIMES, method='cat') for var_name in ['U10', 'V10', 'PM10']]
    pm10_wrf = pm10_wrf.isel(bottom_top=0)
    
    lats, lons = w.latlon_coords(u10_wrf)

    date_range = pd.date_range(start='2019-05-10 00:00:00', end='2019-05-18 00:00:00', freq='6H')
#%%
    proj = ccrs.PlateCarree()
    crange = np.arange(0, 1000+25, 25)
    cmap = cmaps.WhiteBlueGreenYellowRed
    
    for t in date_range:
        fig, ax = plt.subplots(figsize=(8,6), subplot_kw={'projection':proj})

        ax = add_artist(ax, proj)
        # ax = add_province(ax)
        ax = add_citylabel(ax, proj)

        p = ax.pcolor(lons, lats, pm10_wrf.sel(Time=t), vmax=1000, vmin=0, transform=ccrs.PlateCarree(), cmap=cmap)

        cb = fig.colorbar(p, ax=ax, orientation='horizontal', shrink=1, pad=0.07, aspect=25)
        cb.set_label('$\mathrm{\mu g/m^3}$', fontsize=14)
        cb.ax.tick_params(labelsize=12)

        q = ax.quiver(lons.values, lats.values, u10_wrf.sel(Time=t).values, v10_wrf.sel(Time=t).values, color='k', transform=ccrs.PlateCarree(), scale=250, regrid_shape=20)
        ax.quiverkey(q, X=0.9, Y=1.04, U=12, label='12 m/s', labelpos='E', fontproperties={'size':12})

        title_time = t.strftime('%m-%d %H')
        ax.set_title(f"Wrtout 10m Wind & Surface PM10 \n{title_time}UTC", fontsize=14)
        # break

        fig.savefig(f'fig/surface_{title_time}.jpg', dpi=200)
