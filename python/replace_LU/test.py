#%%
import salem
from salem.utils import get_demo_file
import geopandas

# ds = salem.open_xr_dataset('/home/zzhzhao/Model/wrfout/test-9.4/wrfout_d02_2017-07-01_00:00:00.nc')
# ds = salem.open_xr_dataset(get_demo_file('wrfout_d01.nc'))
# t2 = ds.T2.isel(Time=2)
# t2_sub = t2.salem.subset(corners=((77., 25.), (100., 35.)), crs=salem.wgs84)

# f_in = '/home/zzhzhao/code/shpfiles/boundary/hyd1_4m/hyd1_4p.shp'
# shp = geopandas.read_file(f_in, encoding='gbk')
# shp = shp.loc[shp['NAME'].isin(["纳木错", "色林错"])]
# t2_sub = t2_sub.salem.roi(shape=shp)  # add 2 grid points

geo_file = '/home/zzhzhao/Model/tests/test-9.4/WPS/geo_em.d02.nc'
ds = xr.open_dataset(geo_file)

