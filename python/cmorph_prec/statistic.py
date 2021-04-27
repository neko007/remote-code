#%%
import numpy as np 
import xarray as xr
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs 
import cartopy.feature as cfeature
import geopandas
import salem
import cmaps
import warnings
warnings.filterwarnings("ignore")

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

    add_TP(ax, linewidth=1)

def load_TPshp():
    file_path = "/home/zzhzhao/code/shpfiles/Tibet/Tibet.shp"
    shp = geopandas.read_file(file_path, encoding='gbk')
    shp = shp.to_crs(epsg=4326)
    return shp

def add_TP(ax, linewidth=1):
    shp = load_TPshp()
    ax.add_geometries(shp.geometry, crs=ccrs.PlateCarree(), edgecolor='k',
                      alpha=1, facecolor='none', lw=linewidth)

def add_lakes(ax):
    '''
    添加青藏高原湖泊
    '''
    f_in = '/home/zzhzhao/code/shpfiles/TP_lake/青藏高原主要湖泊.shp'
    shp = geopandas.read_file(f_in, encoding='gbk')
    shp1 = shp.to_crs(epsg=4326)
    # shp1.plot(ax=ax, lw=0.3, facecolor='None', edgecolor='k')
    shp1.plot(ax=ax, facecolor='blue', edgecolor='None')

### 区域滑动距平
def _region_rolling_mean(var, rolling_window):
    return var.rolling(lat=rolling_window, lon=rolling_window, center=True).mean()

def region_rolling_anom(var, rolling_window=3):
    return var - _region_rolling_mean(var, rolling_window)

def draw(var, figname, cmap=cmaps.WhiteBlueGreenYellowRed):
    '''
    绘制降水气候态或者频率分布
    '''
    proj = ccrs.PlateCarree()
    # lat_min, lat_max = 25, 43
    # lon_min, lon_max = 72, 105

    fig, ax = plt.subplots(figsize=(8,4), dpi=300, subplot_kw={'projection':proj})
    # set_grid(ax, lat=[lat_min, lat_max], lon=[lon_min, lon_max], span=3)
    ax.set_extent([73, 105, 25, 40], crs=ccrs.PlateCarree())
    add_TP(ax)
    

    p = var.plot.contourf(ax=ax, levels=20, cmap=cmap, robust=True, transform=proj, add_labels=False, add_colorbar=False)
    fig.colorbar(p, ax=ax, orientation='horizontal', aspect=35, shrink=0.75, pad=0.01)

    add_lakes(ax)
    ### 去除地图边框
    ax.background_patch.set_visible(False)
    ax.outline_patch.set_visible(False)
    
    fig.savefig(f'fig/{figname}.jpg', dpi=300, bbox_inches='tight', pad_inches=0.1)

if __name__ == '__main__':
    data_dir = '/home/zzhzhao/data/CMORPH_Prec_TP'
    with xr.open_mfdataset(f'{data_dir}/*.nc')['prec'] as prec:
        prec = prec.where(prec < 200) # 小时降水200以上舍去
        ### 筛选出季节
        # prec_daily = prec.resample(time='D').mean() * 24 # 转换为日平均
        # prec_seasons = dict(prec_daily.groupby('time.season'))
        prec_seasons = dict(prec.groupby('time.season'))
#%% 
    seasons = ['MAM', 'JJA', 'SON', 'DJF']
    # seasons = ['JJA']
    for season in seasons:
        print(f'** {season} **')
        ### 提取青藏高原，各个季节平均
        print('>> 开始计算降水气候态 <<')
        prec_season_mean = prec_seasons[season].mean(dim='time')
        prec_season_mean_TP = prec_season_mean.salem.roi(shape=load_TPshp())

        var1 = prec_season_mean_TP * 24 * 92 # 季度降水
        print('>> 开始绘制降水气候态 <<')
        draw(var1, f'{season}_mean')

        ### 计算降水时次 / 总时次
        print('>> 开始计算降水频率 <<')
        israin = 1e-1
        prec_season_israin = xr.where(prec_seasons[season] >= israin, prec_seasons[season], np.nan).count(dim='time')
        prec_season_count = prec_seasons[season].count(dim='time')
        prec_season_freq = prec_season_israin / prec_season_count

        print('>> 开始绘制降水频率 <<')
        var2 = prec_season_freq.salem.roi(shape=load_TPshp()) * 100
        draw(var2, f'{season}_freq')

#%%
        # 九点距平
        print('>> 计算和绘制气候态的9点滑动区域距平 <<')
        var1_9anom = region_rolling_anom(var1, rolling_window=3)
        draw(var1_9anom, f'{season}_mean_9point', cmap='RdBu_r')
        print('>> 计算和绘制频率的9点滑动区域距平 <<')
        var2_9anom = region_rolling_anom(var2, rolling_window=3)
        draw(var2_9anom, f'{season}_freq_9point', cmap='RdBu_r')





