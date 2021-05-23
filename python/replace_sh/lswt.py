import xarray as xr 
import geopandas
import numpy as np

def nearest_position(latlon, lats, lons):
    lat, lon = latlon
    difflat = lat - lats
    difflon = lon - lons
    rad = np.multiply(difflat,difflat)+np.multiply(difflon , difflon)#difflat * difflat + difflon * difflon
    aa=np.where(rad==np.min(rad))
    ind=np.squeeze(np.array(aa))
    return list(ind)

if __name__ == '__main__':
    lake_init = 277.
    lake_init_flag = 1
    Temp = 3 # Kelvin
    
    # 纳木错范围
    lat_min, lat_max = 30.5, 31. 
    lon_min, lon_max = 90.2, 91.
    
    data_dir = '/home/zzhzhao/Model/tests/test-9.4/WRF/run/'
    geo_file = '/home/zzhzhao/Model/tests/test-9.4/WPS/geo_em.d02.nc'
    wrfinput_file = 'wrfinput_d02'
    with xr.open_dataset(data_dir+wrfinput_file) as ds:
        lats, lons = ds['XLAT'].squeeze(), ds['XLONG'].squeeze()
        latlon_min = nearest_position([lat_min, lon_min], lats.values, lons.values)
        latlon_max = nearest_position([lat_max, lon_max], lats.values, lons.values)

        lu = xr.open_dataset(geo_file)['LU_INDEX']
        lake_mask = xr.where(lu==21, True, False).squeeze()           
        NamCo_mask = lu.where(np.logical_and(latlon_min[0]<lu.south_north, lu.south_north<latlon_max[0]) & np.logical_and(latlon_min[1]<lu.west_east, lu.west_east<latlon_max[1])).squeeze() 
        NamCo_mask = ~ NamCo_mask.isnull()

        lake_mask2 = lake_mask * NamCo_mask 

        lswt = ds['TSK']
        lswt_new = lswt.copy()

        if lake_init_flag == 1:
            ### 替换
            lswt_new.loc[dict(Time=0)] = xr.where(lake_mask2, lake_init, lswt.loc[dict(Time=0)])
        elif lake_init_flag == 2:
            t2 = ds['T2']
            lswt_new.loc[dict(Time=0)] = xr.where(lake_mask, t2.loc[dict(Time=0)]-Temp, lswt.loc[dict(Time=0)])


        del ds['TSK']
        ds = ds.assign(TSK=lswt_new)
        ds.to_netcdf(f'/home/zzhzhao/code/python/replace_sh/{wrfinput_file}', engine='scipy') 

        print('** OK **')
