#%%
from Module import *

def nearest_position(latlon, lats, lons):

    lat, lon = latlon
    difflat = lat - lats
    difflon = lon - lons
    rad = np.multiply(difflat,difflat)+np.multiply(difflon , difflon)#difflat * difflat + difflon * difflon
    aa=np.where(rad==np.min(rad))
    ind=np.squeeze(np.array(aa))
    return tuple(ind)

def add_artist(ax, proj, lat=[29,33], lon=[87,93]):
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
    return ax

def add_NamCo(ax):
    import geopandas
    f_in = '/home/zzhzhao/code/shpfiles/boundary/hyd1_4m/hyd1_4p.shp'
    shp = geopandas.read_file(f_in, encoding='gbk')
    shp2 = shp.loc[shp['NAME'].isin(["纳木错", "色林错"])]
    ax.add_geometries(shp2.geometry, crs=ccrs.PlateCarree(), edgecolor='k',
                      alpha=1, facecolor='none', lw=1)
    return ax

def load_wrfdata(data_dir):
    wrf_files = [f for f in os.listdir(data_dir) if f[9]=='2']
    wrflist = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files]

    rainc = w.getvar(wrflist, 'RAINC', timeidx=w.ALL_TIMES, method='cat')
    rainnc = w.getvar(wrflist, 'RAINNC', timeidx=w.ALL_TIMES, method='cat')
    total_rain = rainc + rainnc

    prec = total_rain.diff('Time', 1)#.sel(Time=pd.date_range('2017-06-01 3:00:00', '2017-06-8 00:00:00', freq='3H'))
    # prec = total_rain.isel(Time=-1)
    lats, lons = w.latlon_coords(prec)
    time = total_rain.Time.to_index() 

    return prec, lats, lons, time 

#%%
if __name__ == '__main__':
    data_dir1 = '/home/zzhzhao/Model/wrfout/test-9.4'
    data_dir2 = '/home/zzhzhao/Model/wrfout/test-9.4-nolake'
    prec1, lats, lons, time = load_wrfdata(data_dir1)
    prec2, lats, lons, time = load_wrfdata(data_dir2) 

    lat_range = (28, 34)
    lon_range = (86, 94)
    ### TRMM资料
    period = pd.date_range('2017-06-01 00:00:00', '2017-06-30 21:00:00', freq='3H')
    trmm_file = '/home/Public_Data/TRMM/TRMM_3hr_3B42_v7/2017.nc'
    trmm = xr.open_dataset(trmm_file)['PRECIPITATION']\
        .sel(TIME=period).sel(LAT=slice(lat_range[0],lat_range[1]), LON=slice(lon_range[0],lon_range[1])) * 3
    lat, lon = trmm.LAT, trmm.LON

    ### 累计降水
    # prec_sum = prec.sel(Time=second_period).sum(dim='Time')
    trmm_sum = trmm.sel(TIME=period).sum(dim='TIME')
    prec1_sum = prec1.sum(dim='Time')
    prec2_sum = prec2.sum(dim='Time')


#%%
    ### 累计降水分布
    proj = ccrs.PlateCarree()
    crange = np.arange(0, 200+10, 10)
    labels = ['TRMM', 'WRF', 'WRF-LakeTurnOff']
    fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(10,10), subplot_kw={'projection':proj})
    fig.subplots_adjust(hspace=0.01)
    for i in range(2):
        for j in range(2):
            axes[i, j] = add_artist(axes[i, j], proj)
            axes[i, j] = add_NamCo(axes[i, j])
    c = axes[0][0].pcolor(lon, lat, trmm_sum, vmin=0, vmax=200, cmap=cmaps.WhiteBlueGreenYellowRed, transform=proj)
    axes[0][0].set_title(labels[0], fontsize=14, weight='bold')
    for j, prec_sum in enumerate([prec1_sum, prec2_sum]):
        c = axes[1][j].pcolor(lons, lats, prec_sum, vmin=0, vmax=200, cmap=cmaps.WhiteBlueGreenYellowRed, transform=proj)
        axes[i][j].set_title(labels[j+1], fontsize=14, weight='bold')
    cb = fig.colorbar(c, ax=axes, orientation='horizontal', pad=0.05, shrink=0.9, aspect=35)
    cb.set_label('Precipitation / mm', fontsize=14)
    axes[0][1].set_visible(False)

    # fig.savefig('fig-test-9.4/prec.jpg', dpi=300)

