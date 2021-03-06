#%%
import wrf as w
import numpy as np
import xarray as xr
from netCDF4 import Dataset
import pandas as pd 
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs
import cmaps
import os 
import sys
sys.path.append('/home/zzhzhao/code/python/wrf-test-14')
from zMap import set_grid, add_NamCo
import salem 
import warnings
warnings.filterwarnings("ignore")

def load_wrfdata(data_dir, domain):
    wrf_files = [f for f in os.listdir(data_dir) if f[11]==f'{domain}']
    wrflist = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files] # 

    T_lake3d = w.getvar(wrflist, 'T_LAKE3D', timeidx=w.ALL_TIMES, method='cat')

    lats, lons = w.latlon_coords(T_lake3d)
    time = T_lake3d.Time.to_index() 

    return T_lake3d, lats, lons, time 

def load_NamCo_shp():
    import geopandas
    f_in = '/home/zzhzhao/code/shpfiles/boundary/hyd1_4m/hyd1_4p.shp'
    shp = geopandas.read_file(f_in, encoding='gbk')
    shp = shp.loc[shp['NAME'].isin(["纳木错"])]
    return shp

def mask_lake(data_dir, shp, testname, domain):
    geo_path = f'/home/zzhzhao/Model/tests/{testname}/WPS/geo_em.d{domain:0>2d}.nc'
    lu = salem.open_wrf_dataset(os.path.join(data_dir, geo_path))['LU_INDEX'].isel(time=0)
    lu_lake = lu.salem.roi(shape=shp)
    mask = xr.where(lu_lake.notnull(), True, False)
    return mask

if __name__ == '__main__':
    data_dir = '/home/zzhzhao/Model/wrfout'
    testname_list = [
        'test-14',
        # 'test-14-oriLD',
        'test-19',
        'test-15',
        # 'test-15-oriLD',
        'test-17',
        # 'test-18',
        'test-20',
        'test-21',
        'test-22',
        'test-23',
        'test-24',
        'test-24-ERA5',
        ]
    N_test = len(testname_list)

    tl_list = dict()
    tl_NamCo_mean_list = dict()
    for testname in testname_list:
        data_path = os.path.join(data_dir, testname)
        domain = 1 if 'ERA5' in testname else 2
        tl, lats, lons, time = load_wrfdata(data_path, domain)
        tl = xr.where(tl>0, tl, np.nan)
        tl_list[testname] = tl

        mask = mask_lake(data_path, load_NamCo_shp(), testname, domain)
        tl_NamCo = tl.where(mask) # 切出NamCo范围
        tl_NamCo_mean = tl_NamCo.mean(dim=['west_east','south_north'])
        # tl_NamCo_mean = tl_NamCo.isel(west_east=74, south_north=39)
        tl_NamCo_mean_list[testname] = tl_NamCo_mean
        
#%%
    crange = np.arange(276, 282+.5, .5)
    cmap = 'rainbow'
    titles = [
        'Wuyang_90m', 
        # 'Wuyang_0.5m',
        'Default_90m', 
        'Default_50m', 
        # 'Default_0.5m',
        'Wuyang_50m', 
        # 'Wuyang_20m',
        'Wuyang_90m_Update',
        'Wuyang_90m_Update2',
        'Wuyang_90m_277K',
        'Wuyang_90m_279.5K',
        'Default_90m_277K',
        'Default_90m_277K_ERA5',
        ]
    ylen = np.ceil(np.sqrt(N_test)).astype(int); xlen = np.ceil(N_test/ylen).astype(int)
    default_len = 5

    fig = plt.figure(figsize=(xlen*default_len, ylen*default_len), dpi=100)
    fig.subplots_adjust(hspace=0.3, wspace=0.2)
    for i, testname in enumerate(testname_list):
        ax = fig.add_subplot(ylen, xlen, i+1)
        var = tl_NamCo_mean_list[testname].T

        var.plot.contourf(levels=crange, cmap=cmap, add_labels=False, ax=ax)
        ax.invert_yaxis()
        ax.set_title(titles[i])
        import matplotlib.dates as mdate  
        ax.xaxis.set_major_formatter(mdate.DateFormatter('%m-%d'))
    
    fig.savefig('fig/Tlake3d_alltest.jpg', dpi=300, bbox_inches='tight', pad_inches=0.1)