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

def area_mean(da, latlon):
    lat, lon = latlon
    d = 0.5
    extent = [lon-d, lon+d, lat-d, lat+d]
    da_area = da.where(lons>extent[0]).where(lons<extent[1]).where(lats>extent[2]).where(lats<extent[3])
    da_area_sum = da_area.mean(dim='south_north').sum(dim='west_east')

    return da_area_sum

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

    ax.coastlines(resolution='50m', lw=0.7)
    ax.set_extent([70, 130, 20, 55], crs=proj)
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), lw=0.7)

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

def add_province(ax):
    '''
    添加省界
    '''
    import geopandas
    f_in = 'boundary/province/Province_9.shp'
    shp = geopandas.read_file(f_in)
    ax.add_geometries(shp.geometry, crs=ccrs.PlateCarree(), edgecolor='k', alpha=0.4, facecolor='none', lw=0.6)
    return ax

def add_citylabel(ax, proj):
    city_names = ['YL', 'Urumqi', 'YC', 'TY', 'BJ']
    city_latlons = [(43.9404, 81.2815), (43.8303, 87.5801), (38.4975, 106.2328), (37.7124, 112.469), (40.0031, 116.407)]
    for i in range(len(city_names)):
        lat, lon = city_latlons[i]
        ax.plot(lon, lat, c='b', marker='.', transform=proj)
        ax.text(lon, lat+1, city_names[i], c='b', fontsize=8, weight=600, horizontalalignment='center')
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

    gh_wrf = w.getvar(wrflist, 'z', timeidx=w.ALL_TIMES, method='cat')
    t_wrf = w.getvar(wrflist, 'temp', timeidx=w.ALL_TIMES, method='cat')
    p = w.getvar(wrflist, "pressure")
    gh_wrf = w.interplevel(gh_wrf, p, 500)
    t_wrf = w.interplevel(t_wrf, p, 500)

    lats, lons = w.latlon_coords(gh_wrf)
    time = gh_wrf.Time.to_index()
    time_modified = time + pd.Timedelta('8 hours')

    ## 再分析资料
    fnl_dir = '/home/zzhzhao/Model/tests/test-chem-6/WPS/'
    date_range = pd.date_range(start='2019-05-10 00:00:00', end='2019-05-18 00:00:00', freq='12H')
    fnl_files = [f'{fnl_dir}met_em.d01.2019-05-{date.day:0>2d}_{date.hour:0>2d}:00:00.nc' for date in date_range]
    fnllist= [Dataset(fnl_file) for fnl_file in fnl_files]

    gh_fnl = w.getvar(fnllist, 'z', timeidx=w.ALL_TIMES, method='cat')
    t_fnl = w.getvar(fnllist, 'TT', timeidx=w.ALL_TIMES, method='cat')

    p_fnl = w.getvar(fnllist, "pressure").isel(west_east=0, south_north=0)
    gh_fnl = gh_fnl.where(p_fnl==500, drop=True).squeeze()
    t_fnl = t_fnl.where(p_fnl==500, drop=True).squeeze()




#%%
    # 500hPa高度场
    gh_list = [gh_fnl, gh_wrf]
    t_list = [t_fnl, t_wrf]

    proj = ccrs.PlateCarree()
    proj2 = w.get_cartopy(gh_wrf)
    crange1 = np.arange(5200, 5920+40, 40)
    crange2 = np.arange(232, 276+4, 4)

    # date = date_range[0]
    for date in date_range:
        fig, axes = plt.subplots(nrows=2, figsize=(5,7), subplot_kw={'projection':proj})
        fig.subplots_adjust(hspace=0.03)
        for i in range(2):
            axes[i] = add_artist(axes[i], proj, lat=[20,55], lon=[70,130])
            # axes[i] = add_province(axes[i])
            axes[i] = add_citylabel(axes[i], proj)

            c = axes[i].contour(lons, lats, gh_list[i].sel(Time=date), levels=crange1, colors='purple', linewidths=1.2, transform=proj, add_labels=False)
            cf = axes[i].contourf(lons, lats, t_list[i].sel(Time=date), levels=crange2, cmap='coolwarm', transform=proj, add_labels=False, add_colorbar=False)
            
            axes[i].clabel(c, inline=1, fontsize=6, fmt='%d')
            title_date = (date+pd.Timedelta('8H')).strftime('%m-%d %H:00:00')
        axes[0].set_title(title_date, fontsize=12, y=0.98)

        cb = fig.colorbar(cf, ax=axes, orientation='horizontal', pad=0.04, shrink=0.95, aspect=30)
        cb.set_label('Temperature (K)', fontsize=10)
        # break
        fig.savefig(f'fig/500hPa_{title_date}.jpg', dpi=300)

    # from concurrent.futures import ThreadPoolExecutor
    # with ThreadPoolExecutor(max_workers=2) as executor:
    #     executor.map(add_province, axes.ravel())
