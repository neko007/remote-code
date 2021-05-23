import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature import NaturalEarthFeature
import geopandas
import salem

### 用于替换老旧的add_artist()
def set_grid(ax, lat=[29,32.5], lon=[87,93], span=1.):
    '''
    魔改matplotlib添加经纬度坐标

    Parameters
    ----------
    lat, lon: 经纬度范围
    '''
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

    ### 限定区域
    extent = [np.min(lon), np.max(lon), np.min(lat), np.max(lat)]
    ax.set_extent(extent, crs=ccrs.PlateCarree())

    lat_span = span # 经纬度间隔
    lon_span = span
    ax.set_xticks(np.arange(lon[0], lon[1]+lon_span, lon_span))
    ax.set_yticks(np.arange(lat[0], lat[1]+lat_span, lat_span))
    ax.tick_params(top=True,bottom=True,left=True,right=True)

    ### 取消坐标单位度°的小圆点
    lon_formatter = LongitudeFormatter(degree_symbol='')
    lat_formatter = LatitudeFormatter(degree_symbol='')
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(axis='both',labelsize=12, direction='out')

    ### 国界线
    linewidth = 0.8
    # ax.add_feature(cfeature.BORDERS.with_scale('10m'), lw=linewidth) # 藏南被吃了，慎用
    f_in = '/home/zzhzhao/code/shpfiles/boundary/bou1_4m/bou1_4l.shp'
    shp = geopandas.read_file(f_in, encoding='gbk')
    ax.add_geometries(shp.geometry, crs=ccrs.PlateCarree(), edgecolor='k', alpha=0.4, facecolor='none', lw=linewidth)
    
    ### 海岸线
    ax.coastlines(resolution='10m', lw=linewidth) 

def add_NamCo(ax, linewidth=1):
    '''
    添加纳木错湖泊
    '''
    f_in = '/home/zzhzhao/code/shpfiles/boundary/hyd1_4m/hyd1_4p.shp'
    shp = geopandas.read_file(f_in, encoding='gbk')
    shp2 = pd.concat([shp[shp['NAME']=="纳木错"], shp[shp['NAME']=="色林错"]], axis=0)
    ax.add_geometries(shp2.geometry, crs=ccrs.PlateCarree(), edgecolor='k',
                      alpha=1, facecolor='none', lw=linewidth)
