#%%
import xarray as xr
import geopandas
import salem # salem用于后处理没毛病，用于中间处理容易产生投影、CRS的问题
from modify_nml2 import max_dom, wps_dir
import os 

def load_NamCo_shp():
    f_in = '/home/zzhzhao/code/shpfiles/boundary/hyd1_4m/hyd1_4p.shp'
    shp = geopandas.read_file(f_in, encoding='gbk')
    shp = shp.loc[shp['NAME'].isin(["纳木错"])]
    return shp

if __name__ == '__main__':
    data_dir = wps_dir
    for i in range(max_dom):
        geo_file = f'geo_em.d0{i+1}.nc'
        
        with xr.open_dataset(f'{data_dir}/{geo_file}') as ds1:
            with salem.open_xr_dataset(f'{data_dir}/{geo_file}') as ds2: # 只能以salem这种方式来读取，不然会损失grid信息
            ### 以两种方式读取的原因是修改salem（内层）读出来的数据，然后填充进纯净的xarray读出的数据（外层）中
                shp = load_NamCo_shp()
                lu = ds2['LU_INDEX']
                lu_lake = lu.salem.roi(shape=shp).isel(Time=0)

                # mask = xr.where(lu_lake==21, True, False)
                mask = lu_lake.notnull()
                del mask.coords['west_east'], mask.coords['south_north']

            ### 变量的metadata和coords要跟原始一致
            ### 修改LU_INDEX、LANDUSEF、LANDMASK
            lu = ds1['LU_INDEX']
            lu_new = lu.copy()
            lu_new.loc[dict(Time=0)] = xr.where(mask, 10, lu.loc[dict(Time=0)])

            luf = ds1['LANDUSEF']
            luf_new = luf.copy()
            luf_new.loc[dict(land_cat=20, Time=0)] = xr.where(mask, 0., luf_new.sel(land_cat=20, Time=0))
            luf_new.loc[dict(land_cat=9, Time=0)] = luf_new.sel(land_cat=9, Time=0) + xr.where(mask, luf.sel(land_cat=20, Time=0), 0.)

            landmask = ds1['LANDMASK']
            landmask_new = landmask.copy()
            landmask_new.loc[dict(Time=0)] = xr.where(mask, 1, landmask.loc[dict(Time=0)])

            del ds1['LU_INDEX'], ds1['LANDUSEF'], ds1['LANDMASK']
            ds1 = ds1.assign(LU_INDEX=lu_new, LANDUSEF=luf_new, LANDMASK=landmask_new) 
            
            if os.path.exists(f'{data_dir}/{geo_file}'):
                print(f'  **** delete {geo_file} ****')
            ds1.to_netcdf(f'{data_dir}/{geo_file}', engine='scipy') # scipy引擎才能输出netcdf3，默认是输出hdf5包装下的netcdf4
            print(f'  **** rewrite {geo_file} ****')
