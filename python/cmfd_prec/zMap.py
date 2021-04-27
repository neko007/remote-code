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

    # ### 国界线
    # linewidth = 1
    # # ax.add_feature(cfeature.BORDERS.with_scale('10m'), lw=linewidth) # 藏南被吃了，慎用
    # f_in = '/home/zzhzhao/code/shpfiles/boundary/bou1_4m/bou1_4l.shp'
    # shp = geopandas.read_file(f_in, encoding='gbk')
    # ax.add_geometries(shp.geometry, crs=ccrs.PlateCarree(), edgecolor='k', alpha=0.4, facecolor='none', lw=linewidth)
    
    # ### 海岸线
    # ax.coastlines(resolution='10m', lw=linewidth) 

    ### 青藏高原
    add_TP(ax, linewidth=1)


def add_TP(ax, linewidth=1):
    file_path = "/home/zzhzhao/code/shpfiles/DBATP/DBATP_Line.shp"
    shp = geopandas.read_file(file_path, encoding='gbk')
    ax.add_geometries(shp.geometry, crs=ccrs.PlateCarree(), edgecolor='k',
                      alpha=1, facecolor='none', lw=linewidth)
    

def _add_subaxes(ax, rect):
    '''
    添加图中图
    '''
    inset_ax = inset_axes(ax, width='18%', height='30%', loc=4)
    inset_ax.set_xticks([ ])
    inset_ax.set_yticks([ ])
    inset_ax.set_xlim(rect[:2])
    inset_ax.set_ylim(rect[2:])
    return inset_ax

# def add_southsea(ax):
#     '''
#     添加南海子图
#     '''
#     f_in = 'boundary/border/国界线.shp'
#     shp = geopandas.read_file(f_in)
#     shp = shp.to_crs(epsg=4326)

#     southsea_range = [105, 125, 0, 25]
#     inset_ax = _add_subaxes(ax, rect=southsea_range)
#     shp.plot(ax=inset_ax, color='k', lw=0.8)

def add_NamCo(ax, linewidth=1):
    '''
    添加纳木错湖泊
    '''
    f_in = '/home/zzhzhao/code/shpfiles/boundary/hyd1_4m/hyd1_4p.shp'
    shp = geopandas.read_file(f_in, encoding='gbk')
    shp2 = pd.concat([shp[shp['NAME']=="纳木错"], shp[shp['NAME']=="色林错"]], axis=0)
    ax.add_geometries(shp2.geometry, crs=ccrs.PlateCarree(), edgecolor='k',
                      alpha=1, facecolor='none', lw=linewidth)

# def add_rivers(ax):
#     '''
#     添加长江、黄河 
#     '''
#     f_in = './boundary/hyd1_4m/hyd1_4l.shp'
#     shp = geopandas.read_file(f_in, encoding='gbk')
#     shp2 = pd.concat([shp[shp['NAME']=='黄河'], shp[shp['NAME']=='长江'],
#                       shp[shp['NAME']=='金沙江'], shp[shp['NAME']=="沱沱河(玛曲)"], shp[shp['NAME']=='通天河']], axis=0)
#     shp2.plot(ax=ax, color='royalblue', lw=0.7)
#     # ax.add_geometries(shp.geometry, crs=ccrs.PlateCarree(), alpha=0.5, edgecolor='b', lw=0.8)
#     return ax

# def add_province(ax):
#     '''
#     添加省界
#     '''
#     f_in = './boundary/province/Province_9.shp'
#     shp = geopandas.read_file(f_in)
#     ax.add_geometries(shp.geometry, crs=ccrs.PlateCarree(), edgecolor='k', alpha=0.4, facecolor='none', lw=0.6)
#     return ax

### 藏南有点问题
def country_mask(ax, c, res='50m', nations=['China', 'Taiwan'], boundary=True):
    '''
    调用cartopy内建的shpfile，对指定国家作白化处理

    Parameters
    ----------
    ax : TYPE
        axes.
    c : TYPE
        contourf的返回值.
    res : {'10m', '50m', '110m'}
        地图分辨率. The default is '50m'.
    nations : list, optional
        需要白化的国家列表. The default is ['China', 'Taiwan'].
    boundary : bool, optional
        是否需要绘制国界线. The default is True.

    Returns
    -------
    ax : TYPE
        axes.

    '''
    shpfilename = shpreader.natural_earth(resolution=res,
                                          category='cultural',
                                          name='admin_0_countries')

    world = geopandas.read_file(shpfilename)
    ### 提取指定国家
    df_nation = world[world['NAME'].map(lambda n: n in nations)]

    multipoly = df_nation.geometry.values
    paths = geos_to_path(list(multipoly))

    ###
    ver = np.concatenate([path.vertices for path in paths])
    cod = np.concatenate([path.codes for path in paths])
    path = Path(ver, cod)

    plate_carre_data_transform = ccrs.PlateCarree()._as_mpl_transform(ax)

    for col in c.collections:
        col.set_clip_path(path, plate_carre_data_transform)

    ### 国界线
    if boundary == True:
        ax.add_geometries(df_nation.geometry, crs=ccrs.PlateCarree(), facecolor='none', edgecolor='k', lw=0.8)

    return ax

### 老版本的set_grid，不再使用
def add_artist(ax, proj, lat=[29,32.5], lon=[87,93]):
    '''
    添加地图要素 → 海岸线，限定区域，网格，省界
    '''
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import matplotlib.ticker as mticker

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

def load_TPshp():
    file_path = "/home/zzhzhao/code/shpfiles/Tibet/Tibet.shp"
    shp = geopandas.read_file(file_path, encoding='gbk')
    shp.crs = "EPSG:4326"
    return shp

def add_TP(ax, linewidth=1):
    shp = load_TPshp()
    ax.add_geometries(shp.geometry, crs=ccrs.PlateCarree(), edgecolor='k', facecolor='k', lw=linewidth)

if __name__ == '__main__':
    fig = plt.figure(figsize=(8,8))
    ax = plt.axes(projection=ccrs.PlateCarree())

    # ax = add_southsea(ax)
    ax.set_extent([70, 140, 10, 55])
    ax.gridlines()
    # ax = add_rivers(ax)

    x = np.arange(70, 140, 0.1)

    y = np.arange(15, 55, 0.1)
    lon, lat = np.meshgrid(x, y)
    z = np.sin(lon) + np.cos(lat)

    # c = ax.contourf(lon, lat, z, transform=ccrs.PlateCarree())
    # ax = country_mask(ax, c)

    # ax = add_rivers(ax)
    # ax = add_province(ax)
    file_path = "/home/zzhzhao/code/shpfiles/Tibet/Tibet.shp"
    shp = geopandas.read_file(file_path, encoding='gbk')
    # shp.crs = "EPSG:4326"
    shp = shp.to_crs(epsg=4326)
    ax.add_geometries(shp.geometry, crs=ccrs.PlateCarree(), edgecolor='k', facecolor='none', lw=1)
