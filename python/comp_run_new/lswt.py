import xarray as xr 
import geopandas
from modify_nml2 import max_dom, wps_dir, wrf_dir, lswt_init 
import os 

if __name__ == '__main__':

    for i in range(1, max_dom):
        geo_file      = os.path.join(wps_dir, f'geo_em.d0{i+1}.nc')
        wrfinput_file = os.path.join(wrf_dir, 'run', f'wrfinput_d0{i+1}')
        with xr.open_dataset(wrfinput_file) as ds:
            
            lu = xr.open_dataset(geo_file)['LU_INDEX']
            lake_mask = xr.where(lu==21, True, False).squeeze()           

            ### 替换
            lswt = ds['TSK']
            lswt_new = lswt.copy()
            lswt_new.loc[dict(Time=0)] = xr.where(lake_mask, lswt_init, lswt.loc[dict(Time=0)])

            del ds['TSK']
            ds = ds.assign(TSK=lswt_new)

            if os.path.exists(wrfinput_file):
                print(f'  **** delete {wrfinput_file} ****')
                ds.to_netcdf(wrfinput_file, engine='scipy') # scipy引擎才能输出netcdf3，默认是输出hdf5包装下的netcdf4
                print(f'  **** rewrite {wrfinput_file} ****')