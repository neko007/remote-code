import xarray as xr 
import geopandas

if __name__ == '__main__':
    lake_init = 277.

    data_dir = '/home/zzhzhao/Model/tests/test-9.4/WRF/run/'
    geo_file = '/home/zzhzhao/Model/tests/test-9.4/WPS/geo_em.d02.nc'
    wrfinput_file = 'wrfinput_d02'
    with xr.open_dataset(data_dir+wrfinput_file) as ds:
        
        lu = xr.open_dataset(geo_file)['LU_INDEX']
        lake_mask = xr.where(lu==21, True, False).squeeze()           

        ### 替换
        lswt = ds['TSK']
        lswt_new = lswt.copy()
        lswt_new.loc[dict(Time=0)] = xr.where(lake_mask, lake_init, lswt.loc[dict(Time=0)])

        del ds['TSK']
        ds = ds.assign(TSK=lswt_new)
        ds.to_netcdf(f'/home/zzhzhao/code/python/replace_sh/{wrfinput_file}', engine='scipy') 

        print('** OK **')
