'''
读取缺失时间csv -> 复制每年第一个文件填补缺失文件（包括大小不匹配的文件） -> 利用ctl读取grd -> 缺测值处理（-999）-> 青藏高原区域裁切 -> 缺测值处理（缺失文件所在时间赋值nan）-> 导出
'''
#%%
from xgrads import CtlDescriptor, open_CtlDataset
import xarray as xr
import numpy as np
import os
import pandas as pd
import subprocess as sp

if __name__ == '__main__':
    year = 2008 # 平年8760，闰年8784
    ctl_path = 'cmorph.ctl'
    data_dir = f'/home/zzhzhao/data/CMORPH_sta/{year}'
    lossdate_dir = '/home/zzhzhao/code/python/cmorph_prec/loss_date'
    head = 'SEVP'
    copy_file = f'{head}_CLI_CHN_MERGE_CMP_PRE_HOUR_GRID_0.10-{year}010100.grd'

    # ctl = CtlDescriptor(file=ctl_path)

    ### 循环填补空缺文件以及大小不为2464000的文件
    losstime_list = pd.read_csv(os.path.join(lossdate_dir, f'{year}.csv'), index_col=None)['date'].values

    for time in pd.DatetimeIndex(losstime_list):
        source_file = os.path.join(data_dir, copy_file)
        target_file = os.path.join(data_dir, f"{head}_CLI_CHN_MERGE_CMP_PRE_HOUR_GRID_0.10-{time.strftime('%Y%m%d%H')}.grd")

        print(f"** 正在复制 {time.strftime('%Y%m%d%H')} **")
        sp.run(f'cp {source_file} {target_file}', shell=True)
    print('>> 复制完毕 <<')

#%%
    with open_CtlDataset(ctl_path) as ds:
        print('>> 正在读取 grd <<')
        rain = ds.crain.squeeze()
        rain = rain.where(rain>-999.).load()
        ### 青藏高原范围
        lat_min, lat_max = 27, 40
        lon_min, lon_max = 73, 105
        ### 提取该范围
        rain_TP = rain.sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))

        ## 将缺失文件填补成nan
        print(f'>> 正在填补缺测值 <<')
        time_index = pd.DatetimeIndex(losstime_list)
        rain_TP.loc[dict(time=time_index)] = np.nan

        ### 输出成新文件
        output_dir = '/home/zzhzhao/data/CMORPH_Prec_TP'
        print(f'>> 正在导出 {year} netcdf <<')
        ds_TP = xr.Dataset({'prec': rain_TP})
        ds_TP.to_netcdf(f'{output_dir}/CMORPH_Prec_TP_{year}.nc', mode='w', format='NETCDF4')
        print(f'>> 导出 {year} netcdf 完成 <<')