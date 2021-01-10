#%%
import sys
sys.path.append("/home/zzhzhao/code/python/wrfchem_2")
from module import *

def load_wrf(date_dir, date1, date2):
    date_range = pd.date_range(date1, date2, freq='1D')
    wrf_files = [f'{data_dir}/wrfout_d01_{date.year}-{date.month:0>2d}-{date.day:0>2d}_{date.hour:0>2d}:00:00' for date in date_range]
    wrflist = [Dataset(wrf_file) for wrf_file in wrf_files]
    return wrflist

if __name__ == '__main__':
    fig_dir = '/home/zzhzhao/code/python/wrfchem_2/fig'
    data_dir = '/home/zzhzhao/code/python/wrfchem_2/data2'
    date1 = '20170629 00:00:00'
    date2 = '20170630 00:00:00'
    wrflist = load_wrf(data_dir, date1, date2)

    date_range = pd.date_range('20170629 18:00:00', '20170630 15:00:00', freq='3H')
    o3 = w.getvar(wrflist, 'o3', timeidx=w.ALL_TIMES, method='cat').sel(Time=date_range, bottom_top=slice(0,12)) * 1e3
    no2 = w.getvar(wrflist, 'no2', timeidx=w.ALL_TIMES, method='cat').sel(Time=date_range, bottom_top=slice(0,12)) * 1e3
    z = w.getvar(wrflist, 'height', timeidx=w.ALL_TIMES, method='cat').isel(Time=0, bottom_top=slice(0,12))
    W = w.getvar(wrflist, 'wa', timeidx=w.ALL_TIMES, method='cat').sel(Time=date_range, bottom_top=slice(0,12))
    
    TJ_ll = (117.3, 39.3)
    TJ_xy = w.ll_to_xy(wrflist, TJ_ll[1], TJ_ll[0])
    o3_TJ = o3.sel(west_east=TJ_xy[0], south_north=TJ_xy[1]) 
    no2_TJ = no2.sel(west_east=TJ_xy[0], south_north=TJ_xy[1]) 
    z_TJ = z.sel(west_east=TJ_xy[0], south_north=TJ_xy[1])
    w_TJ = W.sel(west_east=TJ_xy[0], south_north=TJ_xy[1])

#%% O3
    cmap = 'Blues'
    fig, ax = plt.subplots(dpi=300, figsize=(5,3))
    c = ax.contourf(date_range+pd.Timedelta('8H'), z_TJ, o3_TJ.T, cmap=cmap)
    plt.colorbar(c, label='O3 mix ratio / ppbv')

    ax.set_ylim([0, 3000])
    ax.set_xlim([date_range[0]+pd.Timedelta('7H'), date_range[-1]+pd.Timedelta('9H')])
    ax.set_xticks(date_range+pd.Timedelta('8H'))

    x, y = np.meshgrid(date_range+pd.Timedelta('8H'),z_TJ )
    q = ax.quiver(x, y, 0, w_TJ.T.values, color='k', scale=1.5, width=0.004, headwidth=3.5)#, regrid_shape=20)
    ax.quiverkey(q, X=0.85, Y=1.04, U=0.1, label='10 m/s', labelpos='E', fontproperties={'size':10})

    ax.set_title('2017-06-30 Tianjin', loc='left', fontsize=10, weight='bold')
    ax.set_xlabel('CST', fontsize=10)
    ax.set_ylabel('Height / m', fontsize=10)
    import matplotlib.dates as mdate  
    ax.xaxis.set_major_formatter(mdate.DateFormatter('%H'))

    fig.savefig(f'{fig_dir}/o3_t_h.jpg', dpi=300, bbox_inches='tight', pad_inched=0.1)

#%% NO2
    cmap = cmaps.MPL_YlGn
    fig, ax = plt.subplots(dpi=300, figsize=(5,3))
    c = ax.contourf(date_range+pd.Timedelta('8H'), z_TJ, no2_TJ.T, cmap=cmap)
    plt.colorbar(c, label='NO2 mix ratio / ppbv')

    ax.set_ylim([0, 3000])
    ax.set_xlim([date_range[0]+pd.Timedelta('7H'), date_range[-1]+pd.Timedelta('9H')])
    ax.set_xticks(date_range+pd.Timedelta('8H'))

    x, y = np.meshgrid(date_range+pd.Timedelta('8H'),z_TJ )
    q = ax.quiver(x, y, 0, w_TJ.T.values, color='k', scale=1.5, width=0.004, headwidth=3.5)#, regrid_shape=20)
    ax.quiverkey(q, X=0.85, Y=1.04, U=0.1, label='10 m/s', labelpos='E', fontproperties={'size':10})

    ax.set_title('2017-06-30 Tianjin', loc='left', fontsize=10, weight='bold')
    ax.set_xlabel('CST', fontsize=10)
    ax.set_ylabel('Height / m', fontsize=10)
    import matplotlib.dates as mdate  
    ax.xaxis.set_major_formatter(mdate.DateFormatter('%H'))

    fig.savefig(f'{fig_dir}/no2_t_h.jpg', dpi=300, bbox_inches='tight', pad_inched=0.1)