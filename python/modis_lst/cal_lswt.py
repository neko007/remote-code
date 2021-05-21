'''
计算某日00UTC的纳木错平均温度
1. 如果有数据格点占到整个湖泊的1/3以上，计算Day&Night平均值
2. 如果数据大面积缺失,向后延一天
'''
#%%
import numpy as np
import xarray as xr
import salem 
import pandas as pd
import geopandas 
import warnings
warnings.filterwarnings("ignore")

def load_NamCo_shp():
    import geopandas
    f_in = '/home/zzhzhao/code/shpfiles/boundary/hyd1_4m/hyd1_4p.shp'
    shp = geopandas.read_file(f_in, encoding='gbk')
    shp = shp.loc[shp['NAME'].isin(["纳木错"])]
    return shp

def load_modis(data_path, date, day_or_night='Day'):
    modis = salem.open_xr_dataset(data_path)
    time = modis.indexes['time'].to_datetimeindex()
    modis.coords['time'] = time
    modis = modis.sel(time=date)
    modis_lake = modis[f'LST_{day_or_night}_1km'].salem.roi(shape=load_NamCo_shp()) 

    ### 质量控制
    qc = modis[f'QC_{day_or_night}']
    qc_lake = qc.salem.roi(shape=load_NamCo_shp())
    total_lakegrid = qc_lake.count().values # 湖泊格点总数
    valid_lakegrid = modis_lake.count().values
    if valid_lakegrid / total_lakegrid > 1/2:
        print(f' * {day_or_night}数据有效')
        flag = True
    else:
        print(f' * {day_or_night}数据无效')
        flag = False

    return modis_lake, flag

if __name__ == '__main__':
    date = pd.Timestamp('2013-08-23')
    data_path = f'/home/Public_Data/MODIS/MOD11A1/MOD11A1_NamCo_{date.year}.nc'
    
    while(True):
        print(f"> 今天是{date.strftime('%Y年%m月%d日')}")
        modis_day, flag_day = load_modis(data_path, date, day_or_night='Day')
        modis_night, flag_night = load_modis(data_path, date-pd.Timedelta('1D'), day_or_night='Night')
        if flag_day == True or flag_night == True:
            modis_mean = (modis_day - modis_night) * (10.5-6)/(6+1.5) + modis_night 
            break
        else:
            date += pd.Timedelta('1D')
    print(f'$ 平均湖泊温度 {modis_mean.mean().values:.1f}K $')