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

def add_province(ax):
    '''
    添加省界
    '''
    import geopandas
    f_in = 'boundary/province/Province_9.shp'
    shp = geopandas.read_file(f_in)
    ax.add_geometries(shp.geometry, crs=ccrs.PlateCarree(), edgecolor='k', alpha=0.4, facecolor='none', lw=0.6)
    return ax

#%%

if __name__ == '__main__':

    data_dir = '/home/zzhzhao/Model/tests/test-chem-6/WRF/run'
    date1 = '20190508 00:00:00'
    date2 = '20190519 18:00:00'
    date_range = pd.date_range(date1, date2, freq='1D')

    wrf_files = [f'{data_dir}/wrfout_d01_{date.year}-{date.month:0>2d}-{date.day:0>2d}_{date.hour:0>2d}:00:00' for date in date_range]

    wrflist = [Dataset(wrf_file) for wrf_file in wrf_files]


    vars_dust = [f'DUST_{i+1}' for i in range(5)]
    dust_list = [w.getvar(wrflist, vars_dust[i], timeidx=w.ALL_TIMES, method='cat').isel(bottom_top=0) for i in range(5)]
    dust = xr.concat(dust_list, dim='cat_dust').rename('dust').sum(dim='cat_dust')

    pm10 = w.getvar(wrflist, 'PM10', timeidx=w.ALL_TIMES, method='cat').isel(bottom_top=0)

    lats, lons = w.latlon_coords(dust)
    time = dust.Time.to_index()
    time_modified = time + pd.Timedelta('8 hours')
    


    city_names = ['YL', 'Urumqi', 'YC', 'TY', 'BJ']
    city_latlons = [(43.9404, 81.2815), (43.8303, 87.5801), (38.4975, 106.2328), (37.7124, 112.469), (40.0031, 116.407)]
    n_city = len(city_names)

    city_position = [nearest_position(city_latlons[i], lats, lons) for i in range(n_city)]
    pm10_city = [pm10.sel(south_north=city_position[i][0], west_east=city_position[i][1]) for i in range(n_city)]

#%%
    ## 观测数据
    f_in = 'PM10data.xlsx'
    sheet_order = ['BJ', 'Urumqi', 'TY', 'YC', 'YL']
    data = pd.read_excel(f_in, sheet_name=[0,1,2,3,4], names=['year', 'month', 'day', 'hour', 'pm10'])
    for i in range(n_city):
        df = data[i]
        periods = pd.PeriodIndex(year=df['year'], month=df['month'],day=df['day'], hour=df['hour'], freq="H").to_timestamp()
        data[i].set_index(periods, inplace=True)
        data[i].loc[pd.Timestamp('2019-05-10 17:00:00')] = [2019, 5, 10, 17, 0] # 缺失值
    data_re = {sheet_order[i]:data[i].loc[time_modified] for i in range(n_city)}

#%%
    ## 再分析资料
    fnl_dir = '/home/zzhzhao/Model/tests/test-chem-6/data/'
    date_range = pd.date_range(start='2019-05-11 12:00:00', end='2019-05-13 00:00:00', freq='12H')
    da_list = []
    for date in date_range:
        ds_fnl = xr.open_mfdataset(f'{fnl_dir}fnl_2019{date.month:0>2d}{date.day:0>2d}_{date.hour:0>2d}_00.grib2',
            engine='cfgrib', 
            # parallel=True,
            # concat_dim='time',
            backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {
                "shortName": "gh",
                "typeOfLevel": "isobaricInhPa"
            },
        },)
        da_list.append(ds_fnl.gh.sel(isobaricInhPa=500))
    gh = xr.concat(da_list, dim='time')
    

    # 模式资料
    gh_wrf = w.getvar(wrflist, 'z', timeidx=w.ALL_TIMES, method='cat')
    p = w.getvar(wrflist, "pressure")
    gh_wrf = w.interplevel(gh_wrf, p, 500).sel(Time=gh.time)

#%%
    # 500hPa高度场
    proj = ccrs.PlateCarree()
    proj2 = w.get_cartopy(gh_wrf)
    crange = np.arange(5200, 5920+40, 40)

    fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(13,5), subplot_kw={'projection':proj})
    fig.subplots_adjust(hspace=0.03)

    for i in range(2):
        for j in range(4):
            axes[i][j] = add_artist(axes[i][j], proj, lat=[20,55], lon=[70,130])
            
            if i == 0:
                c = axes[i][j].contourf(gh.longitude, gh.latitude, gh.isel(time=j), levels=crange, cmap='coolwarm', transform=proj, add_labels=False, add_colorbar=False)
                axes[i][j].contour(gh.longitude, gh.latitude, gh.isel(time=j), levels=crange, colors='purple', linewidths=0.6, transform=proj, add_labels=False, add_colorbar=False)
            else:
                c = axes[i][j].contourf(lons, lats, gh_wrf.isel(time=j), levels=crange, cmap='coolwarm', transform=proj, add_labels=False, add_colorbar=False)
                axes[i][j].contour(lons, lats, gh_wrf.isel(time=j), levels=crange, colors='purple', linewidths=0.6, transform=proj, add_labels=False,)
            t = (gh.time+pd.Timedelta('8H')).to_index().strftime('%m-%d %H:00:00')[j]
            axes[0][j].set_title(t)

    axes[0][0].set_ylabel('FNL', labelpad=30, fontsize=10)
    axes[1][0].set_ylabel('WRF-Chem', labelpad=30, fontsize=10)

    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=8) as executor:
        executor.map(add_province, axes.ravel())

    cb = fig.colorbar(c, ax=axes, orientation='horizontal', pad=0.06, shrink=0.9, aspect=45)
    cb.set_label('Geopotential Height (gpm)', fontsize=12)

    fig.savefig('fig/gh.jpg')

#%%
    # 城市PM10时间序列
    fig, axes = plt.subplots(nrows=n_city, ncols=1, figsize=(4,1.5*n_city), sharex=True)
    fig.subplots_adjust(hspace=0.08)
    for i in range(n_city):
        axes[i].plot(time_modified, pm10_city[i], lw=0.6, c='r', marker='o', mfc=None, markersize=0.8, label='WRF-Chem')
        axes[i].plot(time_modified, data_re[city_names[i]]['pm10'], lw=0.6, c='b', label='Obs', marker='o', mfc=None, markersize=0.8,)
        axes[i].set_title(city_names[i], loc='left', y=0.75, x=0.02)
        axes[i].set_ylabel('PM10 ($\mathrm{\mu g / m^3}$)')
    axes[0].legend(fontsize=8, loc='upper right')
    
    import matplotlib.dates as mdate    
    axes[-1].xaxis.set_major_formatter(mdate.DateFormatter('%m-%d'))

    fig.savefig('fig/pm10_ts.jpg')
