#%%
import sys
sys.path.append("/home/zzhzhao/code/python/wrfchem_2")
from module import *

def load_wrf(date_dir, date1, date2):
    date_range = pd.date_range(date1, date2, freq='1D')
    wrf_files = [f'{data_dir}/wrfout_d01_{date.year}-{date.month:0>2d}-{date.day:0>2d}_{date.hour:0>2d}:00:00' for date in date_range]
    wrflist = [Dataset(wrf_file) for wrf_file in wrf_files]
    return wrflist

def add_citylabel(ax, proj):
    city_names = ['TJ']
    city_latlons = [(117.3, 39.3)]
    color = 'k'
    for i in range(len(city_names)):
        lon, lat = city_latlons[i]
        ax.plot(lon, lat, c=color, marker='.', transform=proj, markersize=10)
        ax.text(lon, lat+0.5, city_names[i], c=color, fontsize=16, weight=600, horizontalalignment='center')
    return ax

def add_artist(ax, proj, lat=[25, 45], lon=[110, 125]):
    from cartopy.io.shapereader import Reader
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import matplotlib.ticker as mticker
    '''
    添加地图要素 → 海岸线，限定区域，网格，省界
    '''

    ax.coastlines(resolution='10m', lw=1)
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), lw=1)

    ### 限定区域
    extent = [np.min(lon), np.max(lon), np.min(lat), np.max(lat)]
    ax.set_extent(extent, crs=proj)

    ### 设置网格
    lat_span = 5
    lon_span = 5

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

if __name__ == '__main__':
    # wrfout
    data_dir = '/home/zzhzhao/code/python/wrfchem_2/data2'
    date1 = '20170625 00:00:00'
    date2 = '20170705 00:00:00'
    wrflist = load_wrf(data_dir, date1, date2)
    
    u10, v10, t2m, rh2, o3, no, no2 = [w.getvar(wrflist, var_name, timeidx=w.ALL_TIMES, method='cat') for var_name in ['U10', 'V10', 'T2', 'rh2', 'o3', 'no', 'no2']]
    o3 = o3.isel(bottom_top=0) * 1e3
    nox = (no + no2).isel(bottom_top=0) * 1e3
    t2m = t2m - 273.15
    lats, lons = w.latlon_coords(u10)

    # era5
    data_dir = '/home/zzhzhao/code/python/wrfchem_2/era5_data'
    u10_era, v10_era, t2m_era = [xr.open_dataset(data_dir+'/20170630.nc')[var_name].sel(longitude=slice(110, 125), latitude=slice(45, 25)) for var_name in ['u10', 'v10', 't2m']]
    lat, lon = u10_era.latitude, u10_era.longitude
    t2m_era = t2m_era.copy() - 273.15
#%% Wind & O3
    proj = ccrs.PlateCarree()
    # crange = np.arange(0, 1000+25, 25)
    cmap = cmaps.WhiteBlueGreenYellowRed

    for i in [3, 6, 9, 12]:   
        fig, ax = plt.subplots(figsize=(6,6), subplot_kw={'projection':proj})
        ax = add_artist(ax, proj)
        ax = add_citylabel(ax, proj)
        t = pd.Timestamp(f'20170630 {i:0>2d}:00:00')
        p = ax.pcolor(lons, lats, o3.sel(Time=t), vmin=0, vmax=120, transform=ccrs.PlateCarree(), cmap=cmap)

        cb = fig.colorbar(p, ax=ax, orientation='vertical', shrink=1, pad=0.05, aspect=25)
        cb.set_label('O3 mixing ratio / ppbv', fontsize=12)

        q = ax.quiver(lons.values, lats.values, u10.sel(Time=t).values, v10.sel(Time=t).values, color='k', 
                        transform=ccrs.PlateCarree(), scale=100, regrid_shape=15, width=0.004, headwidth=5)
        ax.quiverkey(q, X=0.9, Y=1.02, U=10, label='10 m/s', labelpos='E', fontproperties={'size':12})

        title_time = (t+pd.Timedelta('8H')).strftime('%m-%d %H')
        ax.set_title(f"WRF 10m Wind & Surface O3 \n{title_time}CST", fontsize=14)
        fig.savefig(f"/home/zzhzhao/code/python/wrfchem_2/fig/surface_o3&wind_{(t+pd.Timedelta('8H')).strftime('%H')}.jpg", dpi=300)
        
#%% WRF Wind & T2m
    proj = ccrs.PlateCarree()
    # crange = np.arange(0, 1000+25, 25)
    cmap = cmaps.MPL_OrRd

    for i in [3, 6, 9, 12]:   
        fig, ax = plt.subplots(figsize=(6,6), subplot_kw={'projection':proj})
        ax = add_artist(ax, proj)
        ax = add_citylabel(ax, proj)
        t = pd.Timestamp(f'20170630 {i:0>2d}:00:00')
        p = ax.pcolor(lons, lats, t2m.sel(Time=t), vmin=20, vmax=40, transform=ccrs.PlateCarree(), cmap=cmap)

        cb = fig.colorbar(p, ax=ax, orientation='vertical', shrink=1, pad=0.05, aspect=25)
        cb.set_label('T2m / $\mathrm{^\circ C}$', fontsize=12)

        q = ax.quiver(lons.values, lats.values, u10.sel(Time=t).values, v10.sel(Time=t).values, color='k', 
                        transform=ccrs.PlateCarree(), scale=100, regrid_shape=15, width=0.004, headwidth=5)
        ax.quiverkey(q, X=0.9, Y=1.02, U=10, label='10 m/s', labelpos='E', fontproperties={'size':12})

        title_time = (t+pd.Timedelta('8H')).strftime('%m-%d %H')
        ax.set_title(f"WRF 10m Wind & T2m \n{title_time}CST", fontsize=14)
        fig.savefig(f"/home/zzhzhao/code/python/wrfchem_2/fig/surface_t2m&wind_{(t+pd.Timedelta('8H')).strftime('%H')}.jpg", dpi=300)
#%% WRF Wind & RH2
    proj = ccrs.PlateCarree()
    cmap = cmaps.cmocean_deep

    for i in [3, 6, 9, 12]:   
        fig, ax = plt.subplots(figsize=(6,6), subplot_kw={'projection':proj})
        ax = add_artist(ax, proj)
        ax = add_citylabel(ax, proj)
        t = pd.Timestamp(f'20170630 {i:0>2d}:00:00')
        p = ax.pcolor(lons, lats, rh2.sel(Time=t), vmin=0, vmax=100, transform=ccrs.PlateCarree(), cmap=cmap)

        cb = fig.colorbar(p, ax=ax, orientation='vertical', shrink=1, pad=0.05, aspect=25)
        cb.set_label('2m Relative Humidity / %', fontsize=12)

        q = ax.quiver(lons.values, lats.values, u10.sel(Time=t).values, v10.sel(Time=t).values, color='k', transform=ccrs.PlateCarree(), scale=200, regrid_shape=20)
        ax.quiverkey(q, X=0.9, Y=1.02, U=10, label='10 m/s', labelpos='E', fontproperties={'size':12})

        title_time = (t+pd.Timedelta('8H')).strftime('%m-%d %H')
        ax.set_title(f"WRF 10m Wind & 2m RH \n{title_time}CST", fontsize=14)
        fig.savefig(f"/home/zzhzhao/code/python/wrfchem_2/fig/surface_RH2&wind_{(t+pd.Timedelta('8H')).strftime('%H')}.jpg", dpi=300)

#%% WRF Wind & NOx
    proj = ccrs.PlateCarree()
    cmap = cmaps.MPL_YlGn

    for i in [3, 6, 9, 12]:   
        fig, ax = plt.subplots(figsize=(6,6), subplot_kw={'projection':proj})
        ax = add_artist(ax, proj)
        ax = add_citylabel(ax, proj)
        t = pd.Timestamp(f'20170630 {i:0>2d}:00:00')
        p = ax.pcolor(lons, lats, nox.sel(Time=t), vmin=0, vmax=100, transform=ccrs.PlateCarree(), cmap=cmap)

        cb = fig.colorbar(p, ax=ax, orientation='vertical', shrink=1, pad=0.05, aspect=25)
        cb.set_label('NOx mix ratio / ppbv', fontsize=12)

        q = ax.quiver(lons.values, lats.values, u10.sel(Time=t).values, v10.sel(Time=t).values, color='k', 
                        transform=ccrs.PlateCarree(), scale=100, regrid_shape=15, width=0.004, headwidth=5)
        ax.quiverkey(q, X=0.9, Y=1.02, U=10, label='10 m/s', labelpos='E', fontproperties={'size':12})

        title_time = (t+pd.Timedelta('8H')).strftime('%m-%d %H')
        ax.set_title(f"WRF 10m Wind & NOx \n{title_time}CST", fontsize=14)
        fig.savefig(f"/home/zzhzhao/code/python/wrfchem_2/fig/surface_NOx&wind_{(t+pd.Timedelta('8H')).strftime('%H')}.jpg", dpi=300)

#%% ERA5 Wind & T2m
    proj = ccrs.PlateCarree()
    # crange = np.arange(0, 1000+25, 25)
    cmap = cmaps.MPL_OrRd

    for i in [3, 6, 9, 12]:   
        fig, ax = plt.subplots(figsize=(6,6), subplot_kw={'projection':proj})
        ax = add_artist(ax, proj)
        ax = add_citylabel(ax, proj)
        t = pd.Timestamp(f'20170630 {i:0>2d}:00:00')
        p = ax.pcolor(lon, lat, t2m_era.sel(time=t), vmin=20, vmax=40, transform=ccrs.PlateCarree(), cmap=cmap)

        cb = fig.colorbar(p, ax=ax, orientation='vertical', shrink=1, pad=0.05, aspect=25)
        cb.set_label('T2m / $\mathrm{^\circ C}$', fontsize=12)

        q = ax.quiver(lon.values, lat.values, u10_era.sel(time=t).values, v10_era.sel(time=t).values, color='k', 
                        transform=ccrs.PlateCarree(), scale=100, regrid_shape=15, width=0.004, headwidth=5)
        ax.quiverkey(q, X=0.9, Y=1.02, U=10, label='10 m/s', labelpos='E', fontproperties={'size':12})

        title_time = (t+pd.Timedelta('8H')).strftime('%m-%d %H')
        ax.set_title(f"ERA5 10m Wind & T2m \n{title_time}CST", fontsize=14)
        fig.savefig(f"/home/zzhzhao/code/python/wrfchem_2/fig/surface_t2m&wind_{(t+pd.Timedelta('8H')).strftime('%H')}_era.jpg", dpi=300)
        # break