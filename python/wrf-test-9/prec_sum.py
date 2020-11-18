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
import os 
import warnings
warnings.filterwarnings("ignore")

def nearest_position(latlon, lats, lons):

    lat, lon = latlon
    difflat = lat - lats
    difflon = lon - lons
    rad = np.multiply(difflat,difflat)+np.multiply(difflon , difflon)#difflat * difflat + difflon * difflon
    aa=np.where(rad==np.min(rad))
    ind=np.squeeze(np.array(aa))
    return tuple(ind)

def add_artist(ax, proj, lat=[29,33], lon=[87,93]):
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import matplotlib.ticker as mticker
    '''
    添加地图要素 → 海岸线，限定区域，网格，省界
    '''
    linewidth = 1
    ax.coastlines(resolution='10m', lw=linewidth)
    ax.add_feature(cfeature.BORDERS.with_scale('10m'), lw=linewidth)

    ### 限定区域
    extent = [np.min(lon), np.max(lon), np.min(lat), np.max(lat)]
    ax.set_extent(extent, crs=proj)

    ### 设置网格
    lat_span = 1
    lon_span = 2

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

def add_NamCo(ax):
    import geopandas
    f_in = '/home/zzhzhao/code/shpfiles/boundary/hyd1_4m/hyd1_4p.shp'
    shp = geopandas.read_file(f_in, encoding='gbk')
    shp2 = pd.concat([shp[shp['NAME']=="纳木错"], shp[shp['NAME']=="色林错"]], axis=0)
    ax.add_geometries(shp2.geometry, crs=ccrs.PlateCarree(), edgecolor='k',
                      alpha=1, facecolor='none', lw=1)
    return ax
#%%
if __name__ == '__main__':
    data_dir = '/home/zzhzhao/Model/wrfout/test-9.1'
    wrf_files = [f for f in os.listdir(data_dir) if f[9]=='2']

    wrflist = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files]

    rainc = w.getvar(wrflist, 'RAINC', timeidx=w.ALL_TIMES, method='cat')
    rainnc = w.getvar(wrflist, 'RAINNC', timeidx=w.ALL_TIMES, method='cat')
    total_rain = rainc + rainnc

    prec = total_rain.diff('Time', 1).sel(Time=pd.date_range('2017-06-04 3:00:00', '2017-06-8 00:00:00', freq='3H'))
    lats, lons = w.latlon_coords(prec)
    time = prec.Time.to_index()

    ### 卫星资料
    trmm_file = 'data/obs.nc'
    trmm = xr.open_dataset(trmm_file)['PRECIPITATION'].sel(TIME=time) * 3
    lat, lon = trmm.LAT, trmm.LON

    ### 累计降水
    # second_period = pd.date_range('2017-06-02 00:00:00', '2017-08-02 00:00:00', freq='6H')
    # prec_sum = prec.sel(Time=second_period).sum(dim='Time')
    # trmm_sum = trmm.sel(TIME=second_period).sum(dim='TIME')
    prec_sum = prec.sum(dim='Time')
    trmm_sum = trmm.sum(dim='TIME')

#%%

    proj = ccrs.PlateCarree()
    crange = np.arange(0, 30+2.5, 2.5)
    labels = ['TRMM', 'WRF']
    fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(10,6), subplot_kw={'projection':proj})
    # fig.subplots_adjust(hspace=0.03)
    for i in range(2):
        axes[i] = add_artist(axes[i], proj)
        axes[i] = add_NamCo(axes[i])
    c = axes[0].contourf(lon, lat, trmm_sum, vmin=2.5, levels=crange, extend='max', cmap=cmaps.WhiteBlueGreenYellowRed, transform=proj, add_labels=False)
    c = axes[1].contourf(lons, lats, prec_sum, vmin=2.5, levels=crange, extend='max', cmap=cmaps.WhiteBlueGreenYellowRed, transform=proj, add_labels=False)

    cb = fig.colorbar(c, ax=axes, orientation='horizontal', pad=0.1, shrink=0.8, aspect=35)
    cb.set_label('Precipitation / mm', fontsize=12)
# %%
