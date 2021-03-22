#%%
import xarray as xr
import geopandas
import salem 
import xesmf

def load_NamCo_shp():
    f_in = '/home/zzhzhao/code/shpfiles/boundary/hyd1_4m/hyd1_4p.shp'
    shp = geopandas.read_file(f_in, encoding='gbk')
    shp = shp.loc[shp['NAME'].isin(["纳木错"])]
    return shp

def regrid(tsk, modis_lst):
    '''
    将modis数据插值到WRF网格(tsk)上
    '''
    target = tsk.isel(Time=0).rename({'XLONG': 'lon', 'XLAT': 'lat'})
    regridder = xesmf.Regridder(modis_lst.isel(time=0), target, 'nearest_s2d')
    modis_regrid = regridder(modis_lst)
    regridder.clean_weight_file() # 清理中间文件
    modis_regrid = modis_regrid.rename({'lon': 'XLONG', 'lat': 'XLAT'})

    return modis_regrid


if __name__ == '__main__':
    data_dir = '/home/zzhzhao/Model/tests/test-9.4/WRF/run/'
    wrfinput_file = 'wrfinput_d01'

    with xr.open_dataset(data_dir+wrfinput_file) as ds1:
        with salem.open_xr_dataset(data_dir+wrfinput_file) as ds2:
            shp = load_NamCo_shp()
            tsk = ds2['TSK']
            tsk_lake = tsk.salem.roi(shape=shp).isel(Time=0)

            mask = tsk_lake.notnull()
            del mask.coords['west_east'], mask.coords['south_north']
    
    ### 读取MOD11A1
    modis = xr.open_dataset('/home/zzhzhao/code/python/replace_sh/MOD11A2.006_1km_aid0001.nc')
    modis_lst = modis['LST_Day_1km']

    ### 将MODIS插值到WRF网格
    modis_regrid = regrid(tsk, modis_lst)
    
    ### 替换
    tsk = ds1['TSK']
    tsk_new = tsk.copy()
    tsk_new.loc[dict(Time=0)] = xr.where(mask, modis_regrid.isel(time=0), tsk.loc[dict(Time=0)])

    del ds1['TSK']
    ds1 = ds1.assign(TSK=tsk_new)
    ds1.to_netcdf(f'/home/zzhzhao/code/python/replace_sh/{wrfinput_file}', engine='scipy') 

# %%
