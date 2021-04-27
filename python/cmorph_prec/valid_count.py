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
from statistic import add_TP, add_lakes, load_TPshp
warnings.filterwarnings("ignore")

if __name__ == '__main__':
    data_dir = '/home/zzhzhao/data/CMORPH_Prec_TP'
    with xr.open_mfdataset(f'{data_dir}/*.nc')['prec'] as prec:
        prec = prec.where(prec < 200) # 小时降水200以上舍去
        ### 筛选出季节
        prec_seasons = dict(prec.groupby('time.season'))

        seasons = ['MAM', 'JJA', 'SON', 'DJF']
        proj = ccrs.PlateCarree()
        cmap = cmaps.WhiteBlueGreenYellowRed
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10,8), dpi=300, subplot_kw={'projection':proj})
        rows = [0, 0, 1, 1]
        cols = [0, 1, 0, 1]
        for i, season in enumerate(seasons):
            print(f'** {season} **')
            valid_freq = prec_seasons[season].count(dim='time') / prec_seasons[season].shape[0]

            print('>> 开始绘图 <<')
            var = valid_freq.salem.roi(shape=load_TPshp()) * 100
            ax = axes[rows[i]][cols[i]]
            
            ax.set_extent([73, 105, 25, 40], crs=ccrs.PlateCarree())
            add_TP(ax)
            
            crange = np.arange(80, 100+1, 1)
            p = var.plot.contourf(ax=ax, levels=crange, cmap=cmap, transform=proj, add_labels=False, add_colorbar=False)

            add_lakes(ax)
            ### 去除地图边框
            ax.background_patch.set_visible(False)
            ax.outline_patch.set_visible(False)
            ax.set_title(season, fontsize=12)

        fig.colorbar(p, ax=axes, orientation='horizontal', aspect=35, shrink=0.75, pad=0.01)
        fig.savefig(f'fig2/valid.jpg', dpi=300, bbox_inches='tight', pad_inches=0.1)
