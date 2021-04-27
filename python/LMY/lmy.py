#%%
import xarray as xr 
import pandas as pd
#%%
if __name__ == '__main__':
    data_path1 = '/home/zzhzhao/data/BCC-CSM2-MR/*.nc'
    with xr.open_mfdataset(data_path1, concat_dim='time') as ds:
        time_index = ds.indexes['time'].to_datetimeindex()
        ds.coords['time'] = time_index
        ds_monthly = ds.resample(time='M', keep_attrs=True).mean()
        ds_monthly.to_netcdf('sithick_monthly_1970-2014.nc', engine='scipy', format='NETCDF3_CLASSIC')
    
#%%
    data_path2 = '/home/zzhzhao/data/sithick_grid_1979_2019_2.nc'
    with xr.open_dataset(data_path2) as ds:
        ds.coords['day'] = pd.date_range('2001-01-01', periods=365)
        ds_monthly = ds.resample(day='M', keep_attrs=True).mean()
        del ds_monthly.coords['day']
        ds_monthly.to_netcdf('sithick_monthly_1979-2019.nc', engine='scipy', format='NETCDF3_CLASSIC')