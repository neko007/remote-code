#%%
import xarray as xr 

if __name__ == '__main__':
    data_dir = '/home/zzhzhao/data/CMFD_Prec'

    with xr.open_mfdataset(f'{data_dir}/*.nc', concat_dim='time', parallel=True) as ds:
        ### 青藏高原范围
        lat_min, lat_max = 27, 40
        lon_min, lon_max = 73, 105
        ### 提取该范围
        ds_TP = ds.sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))
    ### 输出成新文件
    output_dir = '/home/zzhzhao/data/CMFD_Prec_TP'
    print('>> 正在输出netcdf <<')
    ds_TP.to_netcdf(f'{output_dir}/CMFD_Prec_TP_1979-2018.nc', mode='w', format='NETCDF4')