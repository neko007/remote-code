#%%
import wrf as w
import xarray as xr
from netCDF4 import Dataset
import pandas as pd 
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs
import cmaps
import os 
from zMap import set_grid, add_NamCo
import warnings
warnings.filterwarnings("ignore")

def load_wrfdata(data_dir):
    wrf_files = [f for f in os.listdir(data_dir) if f[9]=='2']
    wrflist = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files]

    rainc = w.getvar(wrflist, 'RAINC', timeidx=w.ALL_TIMES, method='cat')
    rainnc = w.getvar(wrflist, 'RAINNC', timeidx=w.ALL_TIMES, method='cat')
    total_rain = rainc + rainnc

    prec = total_rain.diff('Time', 1)
    # prec = total_rain.isel(Time=-1)
    lats, lons = w.latlon_coords(prec)
    time = total_rain.Time.to_index() 

    return prec, lats, lons, time, wrflist

if __name__ == '__main__':
    test_number = 10
    data_dir1 = f'/home/zzhzhao/Model/wrfout/test-{test_number}'
    data_dir2 = f'/home/zzhzhao/Model/wrfout/test-{test_number}-removelake'
    prec1, lats, lons, time, wrflist = load_wrfdata(data_dir1)
    prec2, lats, lons, time, _ = load_wrfdata(data_dir2) 

    lat_range = (28, 34)
    lon_range = (86, 94)
    period = pd.date_range('2017-06-01 00:00:00', '2017-06-30 21:00:00', freq='3H')

    # ### TRMM资料
    # trmm_file = '/home/Public_Data/TRMM/TRMM_3hr_3B42_v7/2017.nc'
    # trmm = xr.open_dataset(trmm_file)['PRECIPITATION']\
    #     .sel(TIME=period).sel(LAT=slice(lat_range[0],lat_range[1]), LON=slice(lon_range[0],lon_range[1])) * 3
    # lat, lon = trmm.LAT, trmm.LON

    ### CMFD资料
    file_path = '/home/zzhzhao/code/python/wrf-test-10/data/prec_CMFD_201706.nc'
    cmfd = xr.open_dataset(file_path)['prec'].sel(lat=slice(lat_range[0],lat_range[1]), lon=slice(lon_range[0],lon_range[1])) * 3
    lat, lon =cmfd.lat, cmfd.lon

    ### 站点观测
    df = pd.read_excel('data/纳木错站2017-2018.xlsx', index_col=0)
    obs_NamCo = df.loc[pd.date_range('2017-06-01','2017-06-30')]['降水量']

#%%
    # NamCo station
    NamCo_latlon = (30.75, 90.9)
    NamCo_position = w.ll_to_xy(wrflist, NamCo_latlon[0], NamCo_latlon[1])
    prec1_NamCo = prec1.sel(south_north=NamCo_position[0], west_east=NamCo_position[1]).resample(Time='1D').sum()
    prec2_NamCo = prec2.sel(south_north=NamCo_position[0], west_east=NamCo_position[1]).resample(Time='1D').sum()
    # trmm_NamCo = trmm.sel(LAT=NamCo_latlon[0], LON=NamCo_latlon[1], method='ffill').resample(TIME='1D').sum()
    cmfd_NamCo = cmfd.sel(lat=NamCo_latlon[0], lon=NamCo_latlon[1], method='ffill').resample(time='1D').sum()

#%%
    fig, ax = plt.subplots(figsize=(9,4))
    ax.plot(obs_NamCo.index, obs_NamCo, lw=1.4, c='k', marker='o', mfc=None, markersize=3.5, label='OBS')
    ax.plot(prec1_NamCo.Time, prec1_NamCo, lw=1.4, c='r', marker='o', mfc=None, markersize=3.5, label='Ctrl')
    ax.plot(prec2_NamCo.Time, prec2_NamCo, lw=1.4, c='g', marker='o', mfc=None, markersize=3.5, label='nolake')
    # ax.plot(trmm_NamCo.TIME, trmm_NamCo, lw=1.4, c='b', marker='o', mfc=None, markersize=3.5, label='TRMM')
    ax.plot(cmfd_NamCo.time, cmfd_NamCo, lw=1.4, c='b', marker='o', mfc=None, markersize=3.5, label='CMFD')
    ax.set_title('NamCo Station', loc='left', y=0.9, x=0.02, fontsize=14, weight='bold')
    ax.set_ylabel('Precipitation $\mathrm{mmd^{-1}}$', fontsize=14)
    ax.legend(fontsize=13, loc='upper right', ncol=2, frameon=False)
    ax.set_ylim([-1, 20])

    import matplotlib.dates as mdate    
    ax.xaxis.set_major_formatter(mdate.DateFormatter('%m-%d'))

    fig.savefig('fig/TS.jpg', dpi=300)
    # fig.savefig('TS.eps',dpi=300,format='eps')
