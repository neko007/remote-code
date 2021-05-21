import xarray as xr 
import geopandas

if __name__ == '__main__':
    lake_init = 277.
    lake_init_flag = 2
    Temp = 3 # Kelvin

    data_dir = '/home/zzhzhao/Model/tests/test-9.4/WRF/run/'
    geo_file = '/home/zzhzhao/Model/tests/test-9.4/WPS/geo_em.d02.nc'
    wrfinput_file = 'wrfinput_d02'
    with xr.open_dataset(data_dir+wrfinput_file) as ds:
        
        lu = xr.open_dataset(geo_file)['LU_INDEX']
        lake_mask = xr.where(lu==21, True, False).squeeze()           

        lswt = ds['TSK']
        lswt_new = lswt.copy()

        if lake_init_flag == 1:
            ### 替换
            lswt_new.loc[dict(Time=0)] = xr.where(lake_mask, lake_init, lswt.loc[dict(Time=0)])
        elif lake_init_flag == 2:
            t2 = ds['T2']
            lswt_new.loc[dict(Time=0)] = xr.where(lake_mask, t2.loc[dict(Time=0)]-Temp, lswt.loc[dict(Time=0)])


        del ds['TSK']
        ds = ds.assign(TSK=lswt_new)
        ds.to_netcdf(f'/home/zzhzhao/code/python/replace_sh/{wrfinput_file}', engine='scipy') 

        print('** OK **')
