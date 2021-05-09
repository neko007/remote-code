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
from statistic import load_TPshp, add_TP, draw
warnings.filterwarnings("ignore")

def cut_NamCo(var):
    lat_min, lat_max = 30, 31.5 
    lon_min, lon_max = 90, 91.5 
    var_NamCo = var.sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))
    return var_NamCo

def add_lakes(ax):
    '''
    添加青藏高原湖泊
    '''
    f_in = '/home/zzhzhao/code/shpfiles/TP_lake/青藏高原主要湖泊.shp'
    shp = geopandas.read_file(f_in, encoding='gbk')
    shp1 = shp.to_crs(epsg=4326)
    # shp1.plot(ax=ax, lw=0.3, facecolor='None', edgecolor='k')
    shp1.plot(ax=ax, facecolor='deepskyblue', edgecolor='None')

def convert_color(var):
    colors = xr.where((var<=3) | (var>21), 'b', var)
    colors = xr.where((var<=9) & (var>3), 'g', colors)
    colors = xr.where((var<=15) & (var>9), 'm', colors)
    colors = xr.where((var<=21) & (var>15), 'r', colors)
    return colors

if __name__ == '__main__':
    data_dir = '/home/zzhzhao/data/CMORPH_Prec_TP'
    with xr.open_mfdataset(f'{data_dir}/*.nc')['prec'] as prec:
        prec = prec.where(prec < 100) # 小时降水100以上舍去
        # prec = cut_NamCo(prec) # 纳木错
        prec_seasons = dict(prec.groupby('time.season'))

#%% 小时降水强度的计算和绘制
    seasons = ['MAM', 'JJA', 'SON', 'DJF']
    proj = ccrs.PlateCarree()
    cmap = cmaps.WhiteBlueGreenYellowRed
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10,8), dpi=300, subplot_kw={'projection':proj})
    fig.subplots_adjust(hspace=0.1)
    rows = [0, 0, 1, 1]
    cols = [0, 1, 0, 1]
    # seasons = ['JJA']
    for i, season in enumerate(seasons):
        print(f'** {season} **')
        ### 提取青藏高原，各个季节平均
        print('>> 开始计算夏季总降水 <<')
        prec_season_mean = prec_seasons[season].sum(dim='time')

        ### 计算降水总时次
        print('>> 开始计算降水总时次 <<')
        israin = 1e-1
        prec_season_israin = xr.where(prec_seasons[season] >= israin, prec_seasons[season], np.nan).count(dim='time')

        ### 计算小时降水强度
        print('>> 计算并绘制小时降水强度 <<')
        prec_hour_intensity = prec_season_mean / prec_season_israin
        var = prec_hour_intensity.salem.roi(shape=load_TPshp())
        
        ax = axes[rows[i]][cols[i]]
        
        ax.set_extent([73, 105, 25, 40], crs=ccrs.PlateCarree())
        add_TP(ax)
        
        crange = np.arange(0, 1.8+0.1, 0.1)
        p = var.plot.contourf(ax=ax, levels=crange, cmap=cmap, transform=proj, add_labels=False, add_colorbar=False)

        add_lakes(ax)
        ### 去除地图边框
        ax.background_patch.set_visible(False)
        ax.outline_patch.set_visible(False)
        ax.set_title(season, fontsize=12)

    fig.colorbar(p, ax=axes, orientation='horizontal', aspect=35, shrink=0.75, pad=0.01)
    fig.savefig(f'fig2/hour_intensity.jpg', dpi=300, bbox_inches='tight', pad_inches=0.1)

#%% 日变化的相位 （还是要把计算和画图模块分开，尤其是面对计算量大的场景）
    seasons = ['MAM', 'JJA', 'SON', 'DJF']

    colors_dict = dict()
    u_dict, v_dict = dict(), dict()
    for season in seasons:
        print(f'** {season} **')
        prec_hourly_cli = prec_seasons[season].groupby('time.hour').mean()
        prec_hourly_cli.coords['hour'] = prec_hourly_cli.coords['hour'].roll(hour=-6, roll_coords=True) # 转换成当地时
        prec_phase = prec_hourly_cli.idxmax(dim='hour').load() # 此处不宜用argmax，因为其返回的是索引

        angle = - (2*np.pi/24 * prec_phase - np.pi/2)
        u, v = np.cos(angle), np.sin(angle)
        colors = convert_color(prec_phase)

        colors_dict.update({season:colors})
        u_dict.update({season:u})
        v_dict.update({season:v})
  
    lon, lat = u_dict[season].coords['lon'], u_dict[season].coords['lat']
    LON, LAT = np.meshgrid(lon.values, lat.values)
#%%
    print('>> Visualization <<')
    proj = ccrs.PlateCarree()
    c = list('bgmr')
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8,6), dpi=300, subplot_kw={'projection':proj})
    fig.subplots_adjust(hspace=0.15, wspace=0.1)
    rows = [0, 0, 1, 1]
    cols = [0, 1, 0, 1]

    for i, season in enumerate(seasons):        
        ax = axes[rows[i]][cols[i]]
        colors, u, v = [var[season] for var in [colors_dict, u_dict, v_dict]]

        ax.set_extent([73, 105, 25, 40], crs=proj) # 青藏高原
        # ax.set_extent([90, 91.5, 30, 31.5], crs=proj) # 纳木错

        add_TP(ax)
        add_lakes(ax)
        ### 去除地图边框
        ax.background_patch.set_visible(False)
        ax.outline_patch.set_visible(False)
        ax.set_title(season, fontsize=12) # 标题
        
        for k in range(4):
            span = 5
            mask = xr.where(colors==c[k], True, False)[::span, ::span]

            u_mask, v_mask = [xr.where(mask, var, np.nan).salem.roi(shape=load_TPshp()).values for var in [u[::span, ::span], v[::span, ::span]]] # 筛选,裁切

            # q = ax.quiver(LON[::span, ::span], LAT[::span, ::span], u_mask, v_mask, color=c[k], transform=proj, scale=30, width=0.005, pivot='mid', headwidth=5, headlength=4)
            q = ax.quiver(LON[::span, ::span], LAT[::span, ::span], u_mask, v_mask, color=c[k], transform=proj, scale=100, width=0.0015, headwidth=4, headlength=4, pivot='mid')

        # ax.quiverkey(q, X=0.85, Y=1.05, U=10, label='10 m/s', labelpos='E', fontproperties={'size':12})
    fig.savefig(f'fig2/phase.jpg', dpi=500, bbox_inches='tight', pad_inches=0.1)
    # fig.savefig(f'fig2/phase_NamCo2.jpg', dpi=500, bbox_inches='tight', pad_inches=0.1)
