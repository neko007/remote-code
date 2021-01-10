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
    data_dir = '/home/zzhzhao/code/python/wrfchem_2/data2'
    date1 = '20170625 00:00:00'
    date2 = '20170705 00:00:00'
    wrflist = load_wrf(data_dir, date1, date2)
    
    o3 = w.getvar(wrflist, 'o3', timeidx=w.ALL_TIMES, method='cat').sel(bottom_top=slice(0, 12)) * 1e3
    
    h = w.getvar(wrflist, 'height', timeidx=-1).sel(bottom_top=slice(0, 12))
    vv = w.getvar(wrflist, 'V', timeidx=w.ALL_TIMES, method='cat').sel(bottom_top=slice(0, 12))
    ww = w.getvar(wrflist, 'W', timeidx=w.ALL_TIMES, method='cat').sel(bottom_top_stag=slice(0, 12))

    lats, lons = w.latlon_coords(o3)
    TJ_ll = (117.3616, 39.3434)
    TJ_xy = w.ll_to_xy(wrflist, TJ_ll[1], TJ_ll[0])

    start_point = w.CoordPair(lat=TJ_ll[1]-5, lon=TJ_ll[0])
    end_point = w.CoordPair(lat=TJ_ll[1]+5, lon=TJ_ll[0])

    # var_TJ_dict = {var_name: var_dict[var_name].sel(west_east=TJ_xy[0], bottom_top=slice(0, 12)) for var_name in ['V', 'o3', 'height']}
    # var_TJ_dict['W'] = var_dict['W'].sel(west_east=TJ_xy[0], bottom_top_stag=slice(0, 12)) 

    o3_cross = w.vertcross(o3, h, wrfin=wrflist, start_point=start_point,
                       end_point=end_point, latlon=True, meta=True)
    v_cross = w.vertcross(vv, h, wrfin=wrflist, start_point=start_point,
                       end_point=end_point, latlon=True, meta=True)
    w_cross = w.vertcross(ww, h, wrfin=wrflist, start_point=start_point,
                       end_point=end_point, latlon=True,  meta=True)

    


#%%
    crange = np.arange(0, 135+15, 15)
    for i in [3, 6, 9, 12]:   
        fig, ax = plt.subplots(figsize=(6,4), dpi=300)
        t = pd.Timestamp(f'20170630 {i:0>2d}:00:00')
        x = [pair.lat for pair in w.to_np(o3_cross.coords['xy_loc'])]
        y = w.to_np(o3_cross.coords['vertical'])
        var = o3_cross.sel(Time=t)

        p = ax.contourf(x, y, var, levels=crange, cmap='Blues')
        cb = fig.colorbar(p, ax=ax, orientation='vertical', shrink=1, pad=0.03, aspect=25)
        cb.set_label('O3 mixing ratio / ppbv', fontsize=10)

        ax.set_ylim([0, 3000])
        
        ax.plot(TJ_ll[1], ax.get_ylim()[0]+20, c='r', markersize=8, marker='^')
        ax.text(TJ_ll[1], ax.get_ylim()[0]-180, 'TJ', c='r', ha='center', weight='bold', fontsize=12)
        ax.set_ylabel('Height / m', fontsize=12)
        ax.set_xlabel('Latitide', fontsize=12)

        # Lat = np.tile(lat, (12,1))
        var1 = v_cross.sel(Time=t)[::3,::3]
        var2 = w_cross.sel(Time=t)[::3,::3] * 25 # w太小需要乘以一个系数
        q = ax.quiver(x[::3], y[::3], var1, var2, color='k', scale=150)
        # ax.quiverkey(q, X=0.9, Y=1.02, U=5, label='5 m/s', labelpos='E', fontproperties={'size':10})

        title_time = (t+pd.Timedelta('8H')).strftime('%m-%d %H')
        ax.set_title(f"WRF O3 profile \n{title_time}CST", fontsize=12)
        # fig.savefig(f"/home/zzhzhao/code/python/wrfchem_2/fig/profile_o3_{(t+pd.Timedelta('8H')).strftime('%H')}.jpg", dpi=300, bbox_inches='tight', pad_inched=0.15)
        # break

