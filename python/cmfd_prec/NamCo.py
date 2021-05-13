'''
纳木错地区区域平均降水时间序列
'''
#%%
import xarray as xr 
import numpy as np
import matplotlib.pyplot as plt

def threshold(prec):
    prec_seasons = dict(prec.groupby('time.season'))
    seasons = ['MAM', 'JJA', 'SON', 'DJF']
    threshold_seasons = dict()
    for season in seasons:
        var = prec_seasons[season].values.ravel()
        threshold_seasons[season] = np.percentile(var, 85)
    return threshold_seasons

def pick_2017(prec_NamCo_daily):
    prec_NamCo_2017 = prec_NamCo_daily.sel(time=slice('2017-5', '2017-10'))
    prec_NamCo_2017.plot.step()
#%%
if __name__ == '__main__':
    data_path = '/home/zzhzhao/data/CMFD_Prec_TP/CMFD_Prec_TP_1979-2018.nc'
    lat_min, lat_max = 29, 32.5 
    lon_min, lon_max = 87, 93 
    with xr.open_dataset(data_path)['prec'] as prec:
        ### 纳木错区域平均
        prec_NamCo = prec.sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))
        prec_NamCo_mean = prec_NamCo.mean(dim=['lat', 'lon'])
        prec_NamCo_daily = prec_NamCo_mean.resample(time='D').mean() * 24
#%%
    threshold_seasons = threshold(prec_NamCo_daily)
    ### 极端降水检测，有为1，无为0
    prec_extreme = prec_NamCo_daily.copy()
    prec_extreme[:] = 0

    for season in ['MAM', 'JJA', 'SON', 'DJF']:
        extreme_index = xr.where((prec_NamCo_daily.time.dt.season==season) & (prec_NamCo_daily>=threshold_seasons[season]), True, False)
        prec_extreme = xr.where(extreme_index, 1, prec_extreme)

#%%
    fig, ax = plt.subplots(figsize=(8,4), dpi=100)
    ax.bar(np.arange(1,367), height=prec_extreme.groupby('time.dayofyear').sum())
    ax.set_xlabel('Dayofyear', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)