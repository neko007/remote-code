#%%
from Module import *
import salem

if __name__ == '__main__':
    geo_file = '/home/zzhzhao/Model/tests/test-9.4/WPS/geo_em.d02.nc'
    geo = salem.open_xr_dataset(geo_file) # 只能以这种方式来读取，不然会损失grid信息
    lu = geo.LU_INDEX.isel(Time=0)
    lats, lons = w.latlon_coords(w.getvar(geo._file_obj.ds, 'LU_INDEX'))

    f_in = '/home/zzhzhao/code/shpfiles/boundary/hyd1_4m/hyd1_4p.shp'
    shp = geopandas.read_file(f_in, encoding='gbk')
    shp = shp.loc[shp['NAME'].isin(["纳木错"])]
    lu = lu.salem.roi(shape=shp)
    mask = xr.where(lu==21, True, False)
    mask.name = 'mask'
    mask.coords['lats'] = lats 
    mask.coords['lons'] = lons
    mask.attrs['long_name'] = 'NamCo Lake Mask (NamCo: True, Others: False)'

    # ds = xr.Dataset({'mask': lu})
    mask.to_netcdf('NamCo_mask.nc', 'w', 'netcdf4')
