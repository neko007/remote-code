import xarray as xr 
import numpy as np
import geopandas
from modify_nml2 import max_dom, wps_dir, wrf_dir, lswt_init, lswt_init_flag
import os 

def nearest_position(latlon, lats, lons):
    lat, lon = latlon
    difflat = lat - lats
    difflon = lon - lons
    rad = np.multiply(difflat,difflat)+np.multiply(difflon , difflon)#difflat * difflat + difflon * difflon
    aa=np.where(rad==np.min(rad))
    ind=np.squeeze(np.array(aa))
    return list(ind)

if __name__ == '__main__':
    Temp = 3
    # 纳木错范围
    lat_min, lat_max = 30.5, 31. 
    lon_min, lon_max = 90.2, 91.

    for i in range(max_dom):
        geo_file      = os.path.join(wps_dir, f'geo_em.d0{i+1}.nc')
        wrfinput_file = os.path.join(wrf_dir, 'run', f'wrfinput_d0{i+1}')
        with xr.open_dataset(wrfinput_file) as ds:
            
            lats, lons = ds['XLAT'].squeeze(), ds['XLONG'].squeeze()
            latlon_min = nearest_position([lat_min, lon_min], lats.values, lons.values)
            latlon_max = nearest_position([lat_max, lon_max], lats.values, lons.values)

            lu = xr.open_dataset(geo_file)['LU_INDEX']
            lake_mask = xr.where(lu==21, True, False).squeeze()           
            NamCo_mask = lu.where(np.logical_and(latlon_min[0]<lu.south_north, lu.south_north<latlon_max[0]) & np.logical_and(latlon_min[1]<lu.west_east, lu.west_east<latlon_max[1])).squeeze() 
            NamCo_mask = ~ NamCo_mask.isnull()

            lake_mask2 = lake_mask * NamCo_mask       

            ### 替换
            lswt = ds['TSK']
            lswt_new = lswt.copy()
            if lswt_init_flag == 1:
                ### 替换
                print(f'  **** lswt_init_flag: {lswt_init_flag} ****')
                print(f'  **** NamCo init LSWT: {lswt_init} ****')
                lswt_new.loc[dict(Time=0)] = xr.where(lake_mask2, lswt_init, lswt.loc[dict(Time=0)])
            elif lswt_init_flag == 2:
                print(f'  **** lswt_init_flag: {lswt_init_flag} ****')
                t2 = ds['T2']
                lswt_new.loc[dict(Time=0)] = xr.where(lake_mask, t2.loc[dict(Time=0)]+Temp, lswt.loc[dict(Time=0)])

            del ds['TSK']
            ds = ds.assign(TSK=lswt_new)

            if os.path.exists(wrfinput_file):
                print(f'  **** delete {wrfinput_file} ****')
                ds.to_netcdf(wrfinput_file, engine='scipy') # scipy引擎才能输出netcdf3，默认是输出hdf5包装下的netcdf4
                print(f'  **** rewrite {wrfinput_file} ****')