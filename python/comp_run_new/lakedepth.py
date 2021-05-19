import xarray as xr 
import xesmf
import salem 
import geopandas
import os 
from modify_nml2 import max_dom, wps_dir
import warnings
warnings.filterwarnings("ignore")

def load_NamCo_shp():
    f_in = '/home/zzhzhao/code/shpfiles/boundary/hyd1_4m/hyd1_4p.shp'
    shp = geopandas.read_file(f_in, encoding='gbk')
    shp = shp.loc[shp['NAME'].isin(["纳木错"])]
    return shp

def regrid(lakedepth, bathy):
    '''
    将obs数据插值到WRF网格上
    '''
    target = lakedepth.isel(Time=0)
    regridder = xesmf.Regridder(bathy, target, 'bilinear')
    bathy_regrid = regridder(bathy)
    regridder.clean_weight_file() # 清理中间文件

    return bathy_regrid

if __name__ == '__main__':
    file_path = '/home/Public_Data/NamCo_Station_Data/NamCo_Bathymetry.nc'
    bathy = salem.open_xr_dataset(file_path)['Bathy']

    data_dir = wps_dir
    for i in range(max_dom):
        geo_file = f'geo_em.d0{i+1}.nc'
        with xr.open_dataset(os.path.join(data_dir, geo_file)) as ds1:
            with salem.open_xr_dataset(os.path.join(data_dir, geo_file)) as ds2:
                lakedepth = ds2['LAKE_DEPTH']
                lakedepth_NamCo = lakedepth.salem.roi(shape=load_NamCo_shp()).isel(Time=0)
                mask = lakedepth_NamCo.notnull()
                del mask.coords['west_east'], mask.coords['south_north']
                
                ### 不加坐标不给插值啊
                lat, lon = ds2['XLAT_M'], ds2['XLONG_M']
                lakedepth.coords['lat'] = lat 
                lakedepth.coords['lon'] = lon

                ### 将观测插值到WRF网格
                bathy_regrid = regrid(lakedepth, bathy)
                ### 填充缺测
                bathy_regrid = bathy_regrid.fillna(0)

                ### 替换
                lakedepth = ds1['LAKE_DEPTH']
                lakedepth_new = lakedepth.copy()
                lakedepth_new.loc[dict(Time=0)] = xr.where(mask, bathy_regrid, lakedepth.loc[dict(Time=0)])

                del ds1['LAKE_DEPTH']
                ds1 = ds1.assign(LAKE_DEPTH=lakedepth_new)

                if os.path.exists(os.path.join(data_dir, geo_file)):
                    print(f'  **** delete {geo_file} ****')
                ds1.to_netcdf(f'{data_dir}/{geo_file}', engine='scipy') # scipy引擎才能输出netcdf3，默认是输出hdf5包装下的netcdf4
                print(f'  **** rewrite {geo_file} ****')
