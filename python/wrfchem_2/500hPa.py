#%%
from module import *

def load_wrf(date_dir, date1, date2):
    date_range = pd.date_range(date1, date2, freq='1D')
    wrf_files = [f'{data_dir}/wrfout_d01_{date.year}-{date.month:0>2d}-{date.day:0>2d}_{date.hour:0>2d}:00:00' for date in date_range]
    wrflist = [Dataset(wrf_file) for wrf_file in wrf_files]
    return wrflist

def add_citylabel(ax, proj):
    city_names = ['TJ']
    city_latlons = [(117.3616, 39.3434)]
    color = 'k'
    for i in range(len(city_names)):
        lon, lat = city_latlons[i]
        ax.plot(lon, lat, c=color, marker='.', transform=proj, markersize=10)
        ax.text(lon, lat+0.5, city_names[i], c=color, fontsize=16, weight=600, horizontalalignment='center')
    return ax

def add_artist(ax, proj, lat=[25, 50], lon=[105, 130]):
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
    data_dir = '/home/zzhzhao/code/python/wrfchem_2/data'
    date1 = '20170625 00:00:00'
    date2 = '20170705 00:00:00'
    wrflist = load_wrf(data_dir, date1, date2)

    gh = w.getvar(wrflist, 'z', timeidx=w.ALL_TIMES, method='cat')
    p = w.getvar(wrflist, 'pressure')
    gh = w.interplevel(gh, p, 500)

    lats, lons = w.latlon_coords(gh)
    time = gh.Time.to_index()
    time_bj = time + pd.Timedelta('8 hours')

    # gh_mean = gh.mean(dim='Time')
    gh_mean = gh.sel(Time='20170630 06:00:00')
#%%
    proj = ccrs.PlateCarree()
    crange = np.arange(5600, 5920+40, 40)
    fig, ax = plt.subplots(figsize=(6,6), subplot_kw={'projection':proj})
    ax = add_artist(ax, proj)
    from scipy.ndimage.filters import gaussian_filter # 平滑等高线
    sigma = 0
    data = gaussian_filter(gh_mean, sigma)
    # y = gaussian_filter(lats, sigma)
    # x = gaussian_filter(lons, sigma)
    # c = ax.contour(x, y, data, levels=crange, transform=proj, colors='b', ls=0.8)
    c = ax.contourf(lons, lats, gh_mean, levels=crange, transform=proj, ls=0.8)
    plt.colorbar(c)
    ax.clabel(c, inline=1, fontsize=12, fmt='%d')
    # ax = add_citylabel(ax, proj)

#%%
    fnl_dir = './data2/'
    date_range = pd.date_range(start='2017-06-30 06:00:00', end='2017-06-30 06:00:00')
    fnl_files = [f'{fnl_dir}met_em.d01.2017-06-{date.day:0>2d}_{date.hour:0>2d}:00:00.nc' for date in date_range]
    fnllist= [Dataset(fnl_file) for fnl_file in fnl_files]

    gh_fnl = w.getvar(fnllist, 'z', timeidx=w.ALL_TIMES, method='cat')
    # t_fnl = w.getvar(fnllist, 'TT', timeidx=w.ALL_TIMES, method='cat')

    p_fnl = w.getvar(fnllist, "pressure").isel(west_east=0, south_north=0)
    gh_fnl = gh_fnl.where(p_fnl==500, drop=True).squeeze()
