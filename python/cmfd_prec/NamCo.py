'''
纳木错地区区域平均降水时间序列
'''
#%%
import xarray as xr 

if __name__ == '__main__':
    data_path = '/home/zzhzhao/data/CMFD_Prec_TP/CMFD_Prec_TP_1979-2018.nc'
    lat_min, lat_max = 29, 32.5 
    lon_min, lon_max = 87, 93 
    with xr.open_dataset(data_path)['prec'] as prec:
        prec_NamCo = prec.sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))
        prec_NamCo_mean = prec_NamCo.mean(dim=['lat', 'lon'])
#%%
        prec_NamCo_daily = prec_NamCo_mean.resample(time='D').mean() * 24