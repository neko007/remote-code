#%%
import sys
sys.path.append("/home/zzhzhao/code/python/wrfchem_2")
from module import *

def load_wrf(date_dir, date1, date2):
    date_range = pd.date_range(date1, date2, freq='1D')
    wrf_files = [f'{data_dir}/wrfout_d01_{date.year}-{date.month:0>2d}-{date.day:0>2d}_{date.hour:0>2d}:00:00' for date in date_range]
    wrflist = [Dataset(wrf_file) for wrf_file in wrf_files]
    return wrflist

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

if __name__ == '__main__':
    fig_dir = '/home/zzhzhao/code/python/wrfchem_2/fig/'
    data_dir = '/home/zzhzhao/code/python/wrfchem_2/data2'
    date1 = '20170625 00:00:00'
    date2 = '20170705 00:00:00'
    wrflist = load_wrf(data_dir, date1, date2)

    var_names = ['o3', 'no', 'no2']
    var_dict = {var_name: w.getvar(wrflist, var_name, timeidx=w.ALL_TIMES, method='cat').isel(bottom_top=0) * 1e3 for var_name in var_names}

    var2_names = ['U10', 'V10', 'T2']
    var2_dict = {var_name: w.getvar(wrflist, var_name, timeidx=w.ALL_TIMES, method='cat') for var_name in var2_names}
    windspeed = np.sqrt(var2_dict['U10']**2 + var2_dict['V10']**2)

    lats, lons = w.latlon_coords(var_dict[var_names[0]])
    time = var_dict[var_names[0]].Time.to_index()
    time_bj = time + pd.Timedelta('8H')

    TJ_ll = (117.3, 39.3)
    TJ_xy = w.ll_to_xy(wrflist, TJ_ll[1], TJ_ll[0])
    var_TJ_dict = {var_name: var_dict[var_name].sel(west_east=TJ_xy[0], south_north=TJ_xy[1]) for var_name in var_names}
    windspeed_TJ = windspeed.sel(west_east=TJ_xy[0], south_north=TJ_xy[1])
    T2_TJ = var2_dict['T2'].sel(west_east=TJ_xy[0], south_north=TJ_xy[1]) - 273.15
   
    ### 观测
    f_obs = '/home/zzhzhao/code/python/wrfchem_2/data2/data_obs.xlsx'
    data = pd.read_excel(f_obs, sheet_name=['O3', 'NO2'])
    def get_date(x):
        return pd.Timestamp(f"{x['date']} {x['hour']}:00:00")
    for var in ['O3', 'NO2']: 
        date_index = data[var].apply(get_date, axis=1)
        data[var].set_index(date_index, inplace=True) 
    sta_index = '1021A'
    cf = 1.97
    obs = {var: data[var].loc[time_bj][sta_index] / cf for var in ['O3', 'NO2']}

#%%
# NO, NO2
    for var in ['o3', 'no2']:
        fig, ax = plt.subplots(figsize=(6,3), dpi=300)
        
        ax.plot(time_bj, var_TJ_dict[var], c='red', lw=0.9, marker='o', markersize=3, label='WRF')
        ax.plot(time_bj, obs[var.upper()], c='blue', lw=0.9, marker='o', markersize=3, label='OBS')
        
        ax.set_ylabel('mixing ratio / ppbv', fontsize=10)
        ax.legend(loc=1, ncol=2, fontsize=10, frameon=False) 
        ax.set_title('Tianjin', fontsize=10, loc='left', weight='bold')
        ax.set_title(var.upper(), fontsize=10, loc='right')
        import matplotlib.dates as mdate    
        ax.xaxis.set_major_formatter(mdate.DateFormatter('%m-%d'))
        fig.savefig(f'{fig_dir}ts_{var.upper()}.jpg', dpi=300)
#%% 14:00 时间序列
    date1_day = '20170625 06:00:00'
    date2_day = '20170705 06:00:00'
    date_range = pd.date_range(date1_day, date2_day, freq='1D')
    var_TJ_day_dict = {var_name: var_TJ_dict[var_name].sel(Time=date_range) for var_name in var_names}
    var_TJ_day_dict['NOx'] = var_TJ_day_dict['no'] + var_TJ_day_dict['no2']

    from scipy.stats import pearsonr
    corr_no2, p = pearsonr(var_TJ_day_dict['o3'].values, var_TJ_day_dict['no2'].values)
    corr_t2, p = pearsonr(var_TJ_day_dict['o3'].values, T2_TJ.sel(Time=date_range).values)
    corr_ws, p = pearsonr(var_TJ_day_dict['o3'].values, windspeed_TJ.sel(Time=date_range).values)

    ### NO2, Wind
    
    fig, ax = plt.subplots(figsize=(6,3), dpi=300)
    # ax.plot(date_range+pd.Timedelta('8H'), var_TJ_day_dict['o3'], c='k', marker='o', markersize=2, label='O3')
    ax.bar(date_range+pd.Timedelta('8H'), var_TJ_day_dict['o3'], color='royalblue', label='O3')
    ax.bar(date_range+pd.Timedelta('8H'), var_TJ_day_dict['no2'], color='goldenrod', label='NO2')

   
    ax.set_ylabel('mixing ratio / ppbv', fontsize=10)
    ax.legend(loc=0, ncol=1, fontsize=10, frameon=False) 
    ax.set_ylim([0, 110])
    ax.set_title('14 CST', fontsize=10, loc='right')
    ax.set_title('Tianjin', fontsize=10, loc='left', weight='bold')
    

    ax1 = ax.twinx()
    l1 = ax1.plot(date_range+pd.Timedelta('8H'), windspeed_TJ.sel(Time=date_range), c='k', marker='o', markersize=4, label='Windspeed')
    # ax1.legend(loc='upper left', ncol=1, fontsize=10, frameon=False) 
    # ax1.tick_params(axis='y', labelcolor='g')
    ax1.set_ylabel('Windspeed / $\mathrm{ms^{-1}}$', c='k', fontsize=10)
    ax1.set_ylim([0, 7])
    ax1.grid(False)

    ax2 = ax.twinx()
    ax2.spines["right"].set_position(("axes", 1.12))
    l2 = ax2.plot(date_range+pd.Timedelta('8H'), T2_TJ.sel(Time=date_range), c='red', ls='--', marker='o', markersize=4, label='T2m')
    ax2.tick_params(axis='y', labelcolor='red')
    ax2.set_ylabel('T2m / $\mathrm{^\circ C}$', color='red', fontsize=10)
    ax2.set_ylim([30, 38])
    # make_patch_spines_invisible(ax2)
    # ax2.spines["right"].set_visible(True)
    ls = l1 + l2
    labs = [l.get_label() for l in ls]
    ax2.legend(ls, labs, loc='upper left', fontsize=10, frameon=False) 

    import matplotlib.dates as mdate  
    ax.xaxis.set_major_formatter(mdate.DateFormatter('%m-%d'))

    # 相关系数
    fig.text(0.12, -0.04, f'Correlation : NO2: {corr_no2:.2f}    Windspeed: {corr_ws:.2f}    T2m: {corr_t2:.2f}', c='blue', fontsize=12)
    fig.savefig(f'{fig_dir}ts.jpg', dpi=300, bbox_inches='tight', pad_inches=0.15)

#%% CO, NO2
    plt.style.use(['science', 'ieee'])
    fig, axes = plt.subplots(figsize=(6,3))
    color = ['k', 'r', 'b']
    ls = ['-', '--', '-.']
    l_list = []
    for i in range(3):
        if var_names[i] == 'co': 
            ax = axes.twinx()
            ax.tick_params(axis='y', labelcolor=color[i])
            ax.set_ylabel('CO mixing ratio / ppmv', c='r')
        else:
            ax = axes
            ax.tick_params(axis='y', labelcolor=color[i])
            ax.set_ylabel('O3 & NO2 mixing ratio / ppmv', c='b')
        l = ax.plot(time_bj, var_NJ_dict[var_names[i]], c=color[i], ls=ls[i], marker='o', markersize=2, label=var_names[i].upper())
        l_list.append(l)
    ls = l_list[0] + l_list[1] + l_list[2]
    labs = [l.get_label() for l in ls]
    ax.legend(ls, labs, loc=1, ncol=2, fontsize=10) 
    # fig.savefig('fig/ts.jpg', dpi=300)