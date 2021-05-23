#%%
import xarray as xr
import geopandas
import salem
from modify_nml2 import max_dom, wps_dir
import os 

def load_mountain_shp():
    f_in = '/home/zzhzhao/code/shpfiles/Nq_shp/Nq.shp'
    shp = geopandas.read_file(f_in, encoding='gbk')
    shp = shp.to_crs(epsg=4326)
    return shp

if __name__ == '__main__':
    data_dir = wps_dir

    for i in range(max_dom):
        geo_file = f'geo_em.d0{i+1}.nc'
        with xr.open_dataset(os.path.join(data_dir, geo_file)) as ds1:
            with salem.open_xr_dataset(os.path.join(data_dir, geo_file)) as ds2: # 只能以salem这种方式来读取，不然会损失grid信息
            ### 以两种方式读取的原因是修改salem（内层）读出来的数据，然后填充进纯净的xarray读出的数据（外层）中
                shp = load_mountain_shp()
                hgt = ds2['HGT_M']
                hgt_mountain = hgt.salem.roi(shape=shp).isel(Time=0)

                # mask = xr.where(lu_lake==21, True, False)
                mask = hgt_mountain.notnull()
                del mask.coords['west_east'], mask.coords['south_north']

            ### 变量的metadata和coords要跟原始一致
            ### 修改HGT_M
            hgt = ds1['HGT_M']
            hgt_new = hgt.copy()
            scale = 0.1
            base_hight = 4718
            new_values = (hgt.loc[dict(Time=0)] - base_hight) * scale + base_hight
            hgt_new.loc[dict(Time=0)] = xr.where(mask, new_values, hgt.loc[dict(Time=0)])

            del ds1['HGT_M']
            ds1 = ds1.assign(HGT_M=hgt_new)

            if os.path.exists(f'{data_dir}/{geo_file}'):
                print(f'  **** delete {geo_file} ****')
            ds1.to_netcdf(f'{data_dir}/{geo_file}', engine='scipy') # scipy引擎才能输出netcdf3，默认是输出hdf5包装下的netcdf4
            print(f'  **** rewrite {geo_file} ****')