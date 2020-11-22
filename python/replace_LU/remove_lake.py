#%%
from Module import *
import salem # salem用于后处理没毛病，用于中间处理容易产生投影、CRS的问题

def load_NamCo_shp():
    f_in = '/home/zzhzhao/code/shpfiles/boundary/hyd1_4m/hyd1_4p.shp'
    shp = geopandas.read_file(f_in, encoding='gbk')
    shp = shp.loc[shp['NAME'].isin(["纳木错"])]
    return shp

if __name__ == '__main__':
    data_dir = '/home/zzhzhao/Model/tests/test-9.4/WPS/'
    geo_file = 'geo_em.d01.nc'
    
    with xr.open_dataset(data_dir+geo_file) as ds1:
        with salem.open_xr_dataset(data_dir+geo_file) as ds2: # 只能以salem这种方式来读取，不然会损失grid信息
            shp = load_NamCo_shp()
            lu = ds2['LU_INDEX']
            lu_lake = lu.salem.roi(shape=shp).isel(Time=0)

            mask = xr.where(lu_lake==21, True, False)
            del mask.coords['west_east'], mask.coords['south_north']

        # 变量的metadata和coords要跟原始一致
        lu = ds1['LU_INDEX']
        lu_new = lu.copy()
        lu_new.loc[dict(Time=0)] = xr.where(mask, 10, lu.loc[dict(Time=0)])

        luf = ds1['LANDUSEF']
        luf_new = luf.copy()
        luf_new.loc[dict(land_cat=20, Time=0)] = xr.where(mask, 0., luf_new.sel(land_cat=20, Time=0))
        luf_new.loc[dict(land_cat=9, Time=0)] = luf_new.sel(land_cat=9, Time=0) + xr.where(mask, luf.sel(land_cat=20, Time=0), 0.)

        del ds1['LU_INDEX'], ds1['LANDUSEF']
        ds1 = ds1.assign(LU_INDEX=lu_new, LANDUSEF=luf_new) 

        ds1.to_netcdf(geo_file, engine='scipy') # scipy引擎才能输出netcdf3，默认是输出hdf5包装下的netcdf4