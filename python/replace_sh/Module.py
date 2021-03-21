### 不同脚本都会用到的包
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
import geopandas
import os 
import warnings
import salem
warnings.filterwarnings("ignore")

def add_artist(ax, proj, lat=[29,32.5], lon=[87,93]):
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

def add_NamCo(ax):
    f_in = '/home/zzhzhao/code/shpfiles/boundary/hyd1_4m/hyd1_4p.shp'
    shp = geopandas.read_file(f_in, encoding='gbk')
    shp2 = pd.concat([shp[shp['NAME']=="纳木错"], shp[shp['NAME']=="色林错"]], axis=0)
    ax.add_geometries(shp2.geometry, crs=ccrs.PlateCarree(), edgecolor='k',
                      alpha=1, facecolor='none', lw=1)

### 用于替换老旧的add_artist()
def set_grid(ax, lat=[29,32.5], lon=[87,93]):
    '''
    魔改matplotlib添加经纬度坐标

    Parameters
    ----------
    lat, lon: 经纬度范围
    '''
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter,AutoMinorLocator

    # 限定区域
    extent = [np.min(lon), np.max(lon), np.min(lat), np.max(lat)]
    ax.set_extent(extent, crs=ccrs.PlateCarree())

    lat_span = .5 # 经纬度间隔
    lon_span = .5
    ax.set_xticks(np.arange(lon[0], lon[1]+lon_span, lon_span))
    ax.set_yticks(np.arange(lat[0], lat[1]+lat_span, lat_span))
    ax.tick_params(top=True,bottom=True,left=True,right=True)

    # 取消坐标单位度°的小圆点
    lon_formatter = LongitudeFormatter(degree_symbol='')
    lat_formatter = LatitudeFormatter(degree_symbol='')
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(axis='both',labelsize=12, direction='out')

    # 国界线
    linewidth = 0.8
    ax.add_feature(cfeature.BORDERS.with_scale('10m'), lw=linewidth)
    ax.coastlines(resolution='10m', lw=linewidth)
