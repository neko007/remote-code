#%%
from Module import *
from prec_sum import add_artist, add_NamCo

def load_wrfdata_T2m(data_dir):
    wrf_files = [f for f in os.listdir(data_dir) if f[9]=='2']
    wrflist = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files]

    t2m = w.getvar(wrflist, 'T2', timeidx=w.ALL_TIMES, method='cat') - 273.15
    lats, lons = w.latlon_coords(t2m)
    time = t2m.Time.to_index() 

    return t2m, lats, lons, time

def load_wrfdata_T500(data_dir):
    wrf_files = [f for f in os.listdir(data_dir) if f[9]=='2']
    wrflist = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files]

    t = w.getvar(wrflist, 'temp', timeidx=w.ALL_TIMES, method='cat', units='degC')
    p = w.getvar(wrflist, 'pressure')
    t500 = w.interplevel(t, p, 500)
    lats, lons = w.latlon_coords(t500)
    time = t500.Time.to_index() 

    return t500, lats, lons, time
    
if __name__ == '__main__':
    data_dir1 = '/home/zzhzhao/Model/wrfout/test-9.4-initLSWT'
    data_dir2 = '/home/zzhzhao/Model/wrfout/test-9.4-initLSWT-laketurnoff'
    t2m1, lats, lons, time = load_wrfdata_T2m(data_dir1)
    t2m2, lats, lons, time = load_wrfdata_T2m(data_dir2)

    # ERA5资料
    era5_file = 'data/2017.06.nc'
    t2m_era5 = xr.open_dataset(era5_file)['t2m'] - 273.15
    lat, lon = t2m_era5.latitude, t2m_era5.longitude

    ### 平均表面温度
    # tsk = tsk.sel(Time=pd.date_range('2017-6-01 18:00:00', periods=30, freq='1D'))
    t2m1_mean = t2m1.mean(dim='Time')
    t2m2_mean = t2m2.mean(dim='Time')
    t2m_era5_mean = t2m_era5.mean(dim='time')
    xr.Dataset({'t2m1':t2m1_mean, 't2m2':t2m2_mean, 'lat':lats, 'lon':lons}).to_netcdf('data/t2m_mean.nc')
#%%
    ### 表面温度分布
    proj = ccrs.PlateCarree()
    labels = ['ERA5', 'WRF', 'WRF-LakeTurnOff']
    fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(10,10), subplot_kw={'projection':proj})
    fig.subplots_adjust(hspace=0.01)
    for i in range(2):
        for j in range(2):
            axes[i, j] = add_artist(axes[i, j], proj)
            axes[i, j] = add_NamCo(axes[i, j])
    c = axes[0][0].pcolor(lon, lat, t2m_era5_mean, vmin=-10, vmax=20, cmap=cmaps.WhiteBlueGreenYellowRed, transform=proj)
    # c = axes[0][0].pcolor(era5_skt.longitude, era5_skt.latitude, era5_skt, vmin=-10, vmax=20, cmap=cmaps.WhiteBlueGreenYellowRed, transform=proj)
    axes[0][0].set_title(labels[0], fontsize=14, weight='bold')
    for j, tsk_mean in enumerate([t2m1_mean, t2m2_mean]):
        c = axes[1][j].pcolor(lons, lats, tsk_mean, vmin=-10, vmax=20, cmap=cmaps.WhiteBlueGreenYellowRed, transform=proj)
        axes[i][j].set_title(labels[j+1], fontsize=14, weight='bold')
    cb = fig.colorbar(c, ax=axes, orientation='horizontal', pad=0.05, shrink=0.9, aspect=35)
    cb.set_label('2m temperature / $\mathrm{^\circ C}$', fontsize=14)
    axes[0][1].set_visible(False)

    fig.savefig('fig-test-9.4/t2m2.jpg', dpi=300)

#%%
    t5001, lats, lons, time = load_wrfdata_T500(data_dir1)
    t5002, lats, lons, time = load_wrfdata_T500(data_dir2)
    t5001_mean = t5001.mean(dim='Time')
    t5002_mean = t5002.mean(dim='Time')

    # ERA5资料
    era5_file2 = 'data/t_500hpa_201706.nc'
    t500_era5 = xr.open_dataset(era5_file2)['t'] - 273.15
    t500_era5_mean = t500_era5.mean(dim='time')

#%%
    ### 500hPa温度分布
    proj = ccrs.PlateCarree()
    labels = ['ERA5', 'WRF', 'WRF-LakeTurnOff']
    fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(10,10), subplot_kw={'projection':proj})
    fig.subplots_adjust(hspace=0.01)
    for i in range(2):
        for j in range(2):
            axes[i, j] = add_artist(axes[i, j], proj)
            axes[i, j] = add_NamCo(axes[i, j])
    c = axes[0][0].pcolor(lon, lat, t500_era5_mean, vmin=-4, vmax=1, cmap='Spectral_r', transform=proj)
    # c = axes[0][0].pcolor(era5_skt.longitude, era5_skt.latitude, era5_skt, vmin=-10, vmax=20, cmap=cmaps.WhiteBlueGreenYellowRed, transform=proj)
    axes[0][0].set_title(labels[0], fontsize=14, weight='bold')
    for j, tsk_mean in enumerate([t5001_mean, t5002_mean]):
        c = axes[1][j].pcolor(lons, lats, tsk_mean, vmin=-4, vmax=1, cmap='Spectral_r', transform=proj)
        axes[i][j].set_title(labels[j+1], fontsize=14, weight='bold')
    cb = fig.colorbar(c, ax=axes, orientation='horizontal', pad=0.05, shrink=0.9, aspect=35)
    cb.set_label('temperature @500hPa / $\mathrm{^\circ C}$', fontsize=14)
    axes[0][1].set_visible(False)

    fig.savefig('fig-test-9.4/t5002.jpg', dpi=300)

#%%
    ### 表面-500hPa温度分布
    proj = ccrs.PlateCarree()
    labels = ['ERA5', 'WRF', 'WRF-LakeTurnOff']
    fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(10,10), subplot_kw={'projection':proj})
    fig.subplots_adjust(hspace=0.01)
    for i in range(2):
        for j in range(2):
            axes[i, j] = add_artist(axes[i, j], proj)
            axes[i, j] = add_NamCo(axes[i, j])
    c = axes[0][0].pcolor(lon, lat, t2m_era5_mean-t500_era5_mean, vmin=-10, vmax=20, cmap=cmaps.WhiteBlueGreenYellowRed, transform=proj)
    # c = axes[0][0].pcolor(era5_skt.longitude, era5_skt.latitude, era5_skt, vmin=-10, vmax=20, cmap=cmaps.WhiteBlueGreenYellowRed, transform=proj)
    axes[0][0].set_title(labels[0], fontsize=14, weight='bold')
    for j, tsk_mean in enumerate([t2m1_mean-t5001_mean, t2m2_mean-t5002_mean]):
        c = axes[1][j].pcolor(lons, lats, tsk_mean, vmin=-10, vmax=20, cmap=cmaps.WhiteBlueGreenYellowRed, transform=proj)
        axes[i][j].set_title(labels[j+1], fontsize=14, weight='bold')
    cb = fig.colorbar(c, ax=axes, orientation='horizontal', pad=0.05, shrink=0.9, aspect=35)
    cb.set_label('2m-500hPa temperature / $\mathrm{^\circ C}$', fontsize=14)
    axes[0][1].set_visible(False)

    fig.savefig('fig-test-9.4/diff2.jpg', dpi=300)