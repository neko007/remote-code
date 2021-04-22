'''
有些年份grd文件不全，先用每年第一个grd的复制填上，等xgrads读进来再修改对应时间为nan
'''
#%%
import pandas as pd
import numpy as np 
import os 

if __name__ == '__main__':
    # year = 2016
    for year in range(2008, 2009):
        file_dir = f'/home/zzhzhao/data/CMORPH_sta/{year}'
        time_range = pd.date_range(f'{year}-01-01 00:00:00', f'{year}-12-31 23:00:00 ', freq='h')
        losstime_list = []
        for time in time_range:
            # grd_file = f"SURF_CLI_CHN_MERGE_CMP_PRE_HOUR_GRID_0.10-{time.strftime('%Y%m%d%H')}.grd"
            grd_file = f"SEVP_CLI_CHN_MERGE_CMP_PRE_HOUR_GRID_0.10-{time.strftime('%Y%m%d%H')}.grd"
            grd_path = os.path.join(file_dir, grd_file)
            if os.path.exists(grd_path) == 0 or (os.path.getsize(grd_path) != 2464000):
                losstime_list.append(time)

        pd.DataFrame(losstime_list, index=None, columns=['date']).to_csv(f'loss_date/{year}.csv')

#%%
    ### 2015年检索出了点问题   
    pd.concat([
        pd.read_csv('loss_date/2015_1.csv',index_col=0), 
        pd.read_csv('loss_date/2015.csv',index_col=0)
        ], ignore_index=True).to_csv('loss_date/2015_2.csv')
    