#%%
import sys
sys.path.append('/home/zzhzhao/code/python/replace_LU')
from Module import *

def load_wrfdata_windspeed(data_dir):
    wrf_files = [f for f in os.listdir(data_dir) if f[9]=='2']
    wrflist = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files]

    uu, ww, rh = [w.getvar(wrflist, var, timeidx=w.ALL_TIMES, method='cat') for var in ['ua', 'wa', 'rh']]
    h, ter = [w.getvar(wrflist, var, timeidx=-1) for var in ['height', 'ter']]

    lats, lons = w.latlon_coords(uu)
    time = uu.Time.to_index() 

    return uu, ww, rh, h, ter, lats, lons, time, wrflist

def xy(var_cross):
        x = [pair.lon for pair in w.to_np(var_cross.coords['xy_loc'])]
        y = w.to_np(var_cross.coords['vertical'])
        return x, y

def fill_gap(var_cross):
    #* 消除contourf与地形填色之间的空隙     
    var_cross_filled = np.ma.copy(w.to_np(var_cross))
    for i in range(var_cross_filled.shape[-1]):
        column_vals = var_cross_filled[:, i]
        first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
        var_cross_filled[0:first_idx, i] = var_cross_filled[first_idx, i]
    
    return var_cross_filled

#%%
if __name__ == '__main__':
    data_dir1 = '/home/zzhzhao/Model/wrfout/test-9.4-initLSWT'
    data_dir2 = '/home/zzhzhao/Model/wrfout/test-9.4-initLSWT-laketurnoff'
    u1, w1, rh1, h, ter, lats, lons, time, wrflist = load_wrfdata_windspeed(data_dir1) 

    NamCo_ll = (30.75, 90.9)
    # NamCo_xy = w.ll_to_xy(wrflist, NamCo_ll[0], NamCo_ll[1])

    start_point = w.CoordPair(lat=NamCo_ll[0], lon=NamCo_ll[1]-1)
    end_point = w.CoordPair(lat=NamCo_ll[0], lon=NamCo_ll[1]+1)
    
    u1_cross = w.vertcross(u1, h, wrfin=wrflist, start_point=start_point,
                       end_point=end_point, latlon=True, meta=True)
    w1_cross = w.vertcross(w1, h, wrfin=wrflist, start_point=start_point,
                       end_point=end_point, latlon=True,  meta=True)
    rh1_cross = w.vertcross(rh1, h, wrfin=wrflist, start_point=start_point,
                       end_point=end_point, latlon=True,  meta=True)
    ter_line = w.interpline(ter, wrfin=wrflist, start_point=start_point,
                      end_point=end_point)

#%%
    fig, ax = plt.subplots(figsize=(8,4), dpi=300)
    t = pd.Timestamp(f'20170610 12:00:00')
    x, y = xy(u1_cross)     

    # add NamCo station
    label_color = 'white'
    ax.plot(NamCo_ll[1], 4760, c=label_color, markersize=5, marker='^')
    ax.text(NamCo_ll[1], 4580, 'NamCo', c=label_color, ha='center',weight='bold', fontsize=8)

    # RH
    var = fill_gap(rh1_cross.sel(Time=t))
    p = ax.contourf(x, y, var, levels=np.arange(0, 100+5, 5), cmap='Blues')
    cb = fig.colorbar(p, ax=ax, orientation='vertical', shrink=1, pad=0.03, aspect=25)
    cb.set_label('Relative humidity / %')

    # UW
    span = 1
    var1 = u1_cross.sel(Time=t)[::span, ::span]
    var2 = w1_cross.sel(Time=t)[::span,::span] * 25
    q = ax.quiver(x[::span], y[::span], var1, var2, color='k', scale=300)

    # add terrain
    ht_fill = ax.fill_between(x, 0, w.to_np(ter_line),
                                facecolor="saddlebrown")

    # adjust 
    ax.set_ylim([4500, 7500])
    ax.set_ylabel('Height / m')
    ax.set_xlabel('Longitude')
