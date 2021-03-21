#%%
import xarray as xr
import geopandas
import salem 

def load_NamCo_shp():
    f_in = '/home/zzhzhao/code/shpfiles/boundary/hyd1_4m/hyd1_4p.shp'
    shp = geopandas.read_file(f_in, encoding='gbk')
    shp = shp.loc[shp['NAME'].isin(["纳木错"])]
    return shp


if __name__ == '__main__':
    data_dir = '/home/zzhzhao/Model/tests/test-9.4/WRF/run/'
    wrfinput_file = 'wrfinput_d02'

    with xr.open_dataset(data_dir+wrfinput_file) as ds1:
        with salem.open_xr_dataset(data_dir+wrfinput_file) as ds2:
            shp = load_NamCo_shp()
            tsk = ds2['TSK']
            tsk_lake = tsk.salem.roi(shape=shp).isel(Time=0)

            mask = tsk_lake.notnull()
            # del mask.coords['west_east'], mask.coords['south_north']