# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 17:04:52 2020

@author: zzz
"""

import numpy as np
import pandas as pd
import xarray as xr
import geopandas
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# from matplotlib import rcParams
# rcParams['font.family'] = 'Arial'


def add_artist(ax, proj, lat=[27,35], lon=[85,95]):
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import matplotlib.ticker as mticker
    '''
    添加地图要素 → 海岸线，限定区域，网格，省界
    '''
    linewidth = 1
    ax.coastlines(resolution='10m', lw=linewidth)
    ax.set_extent([70, 130, 20, 55], crs=proj)
    ax.add_feature(cfeature.BORDERS.with_scale('10m'), lw=linewidth)

    ### 限定区域
    extent = [np.min(lon), np.max(lon), np.min(lat), np.max(lat)]
    ax.set_extent(extent, crs=proj)

    ### 设置网格
    lat_span = 2
    lon_span = 2

    ax.set_xticks([], crs=ccrs.PlateCarree())
    ax.set_yticks([], crs=ccrs.PlateCarree())

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color=None, alpha=0, linestyle='--')
    gl.xlabels_bottom = False  # 关闭顶端的经纬度标签
    gl.ylabels_right = False  # 关闭右侧的经纬度标签
    gl.xformatter = LONGITUDE_FORMATTER  # x轴设为经度的格式
    gl.yformatter = LATITUDE_FORMATTER  # y轴设为纬度的格式
    gl.xlocator = mticker.FixedLocator(np.arange(extent[0], extent[1]+lon_span, lon_span))
    gl.ylocator = mticker.FixedLocator(np.arange(extent[2], extent[3]+lat_span, lat_span))
    gl.xlabel_style = {'size':14}
    gl.ylabel_style = {'size':14}
    return ax

def add_NamCo(ax):
    f_in = '/home/zzhzhao/code/shpfiles/boundary/hyd1_4m/hyd1_4p.shp'
    shp = geopandas.read_file(f_in, encoding='gbk')
    shp2 = pd.concat([shp[shp['NAME']=="纳木错"], shp[shp['NAME']=="色林错"]], axis=0)
    ax.add_geometries(shp2.geometry, crs=ccrs.PlateCarree(), edgecolor='k',
                      alpha=1, facecolor='deepskyblue', lw=0)
    return ax

def add_terrain(fig, ax):
    ds = xr.open_dataset('geo_output.nc')
    hgt = ds.HGT_M.squeeze()
    lats = ds.XLAT_M.squeeze()
    lons = ds.XLONG_M.squeeze()

    crange = np.arange(3000, 7000+100, 100)
    # p = ax.pcolor(lons, lats, hgt, vmin=3000, vmax=7000, transform=ccrs.PlateCarree(), cmap='terrain')
    p = ax.contourf(lons, lats, hgt, levels=crange, extend='both', transform=ccrs.PlateCarree(), cmap='terrain')

    position = fig.add_axes([0.16, 0.06, 0.57, 0.035])
    cb = plt.colorbar(p, cax=position, orientation='horizontal')
    # cb.set_label('height / m', ha='left', va='baseline', fontsize=12, position=(1.1, -1), rotation=0)
    fig.text(0.75, 0.065, 'height / m', fontsize=14, weight='bold')

    return ax, p
    
if __name__ == '__main__':

    proj = ccrs.PlateCarree()
    fig, ax = plt.subplots(figsize=(8,6), subplot_kw={'projection':proj})
    ax = add_artist(ax, proj)
    ax = add_NamCo(ax)

    ax, p = add_terrain(fig, ax)

    ax.text(90.7, 31, 'Nam Co', fontsize=14, fontweight=600)
    ax.text(89, 32, 'Sering Co', fontsize=14, fontweight=600)

    ax.plot(90.9, 30.75, c='k', transform=proj, marker='o', markersize=5)

    # fig.savefig('fig/terrain2.png', dpi=300)
    # fig.savefig('fig/terrain2.eps', dpi=300, format='eps')






