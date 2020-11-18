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
import xesmf as xe
import cmaps
import warnings
warnings.filterwarnings("ignore")

plt.style.use(['science','ieee'])

def nearest_position(latlon, lats, lons):

    lat, lon = latlon
    difflat = lat - lats
    difflon = lon - lons
    rad = np.multiply(difflat,difflat)+np.multiply(difflon , difflon)#difflat * difflat + difflon * difflon
    aa=np.where(rad==np.min(rad))
    ind=np.squeeze(np.array(aa))
    return tuple(ind)

def add_artist(ax, proj, lat=[22, 40], lon=[112, 124]):
    from cartopy.io.shapereader import Reader
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import matplotlib.ticker as mticker
    '''
    添加地图要素 → 海岸线，限定区域，网格，省界
    '''

    ax.coastlines(resolution='50m', lw=0.8)
    ax.set_extent([70, 130, 20, 55], crs=proj)
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), lw=0.8)

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
    gl.xlabel_style = {'size':8}
    gl.ylabel_style = {'size':8}
    return ax

#%%

if __name__ == '__main__':

    data_dir = '/home/zzhzhao/Model/tests/test-chem-6/WRF/run'
    date1 = '20190508 00:00:00'
    date2 = '20190519 18:00:00'
    date_range = pd.date_range(date1, date2, freq='1D')

    wrf_files = [f'{data_dir}/wrfout_d01_{date.year}-{date.month:0>2d}-{date.day:0>2d}_{date.hour:0>2d}:00:00' for date in date_range]

    wrflist = [Dataset(wrf_file) for wrf_file in wrf_files]

    rainc = w.getvar(wrflist, 'RAINC', timeidx=w.ALL_TIMES, method='cat')

    lats, lons = w.latlon_coords(rainc)
    time = rainc.Time.to_index()
    time_modified = time + pd.Timedelta('8 hours')


    city_names = ['YL', 'Urumqi', 'YC', 'TY', 'BJ']
    city_latlons = [(43.9404, 81.2815), (43.8303, 87.5801), (38.4975, 106.2328), (37.7124, 112.469), (40.0031, 116.407)]
    n_city = len(city_names)

    city_position = [nearest_position(city_latlons[i], lats, lons) for i in range(n_city)]
    rainc_city = [rainc.sel(south_north=city_position[i][0], west_east=city_position[i][1]) for i in range(n_city)]

#%%
    # 城市PM10时间序列
    fig, axes = plt.subplots(nrows=n_city, ncols=1, figsize=(4,1.5*n_city), sharex=True)
    fig.subplots_adjust(hspace=0.08)
    for i in range(n_city):
        axes[i].plot(time_modified, rainc_city[i], lw=0.6, c='r', marker='o', mfc=None, markersize=0.8, label='WRF-Chem')
        axes[i].set_title(city_names[i], loc='left', y=0.75, x=0.02)
        axes[i].set_ylabel('RAINC / mm')
    # axes[0].legend(fontsize=8, loc='upper right')
    
    import matplotlib.dates as mdate    
    axes[-1].xaxis.set_major_formatter(mdate.DateFormatter('%m-%d'))

    # fig.savefig('fig/pm10_ts.jpg')
