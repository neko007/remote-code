#%%
import sys
import wrf as w
import numpy as np
import xarray as xr
from netCDF4 import Dataset
import pandas as pd 
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs
import cmaps
import os 
from zMap import set_grid, add_NamCo
import warnings
warnings.filterwarnings("ignore")

def load_wrfdata_windspeed(data_dir):
    wrf_files = [f for f in os.listdir(data_dir) if f[11]=='2']
    wrflist = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files]

    uu, vv, ww, qv = [w.getvar(wrflist, var, timeidx=w.ALL_TIMES, method='cat') for var in ['ua', 'va', 'wa', 'QVAPOR']]
    h, ter = [w.getvar(wrflist, var, timeidx=-1) for var in ['height', 'ter']]

    lats, lons = w.latlon_coords(uu)
    time = uu.Time.to_index() 

    return uu, vv, ww, qv, h, ter, lats, lons, time, wrflist

def xy(var_cross):
        x_labels = [pair.latlon_str() for pair in w.to_np(var_cross.coords['xy_loc'])]
        y = w.to_np(var_cross.coords['vertical'])
        return x_labels, y

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
    test_number = 15
    data_dir1 = f'/home/zzhzhao/Model/wrfout/test-{test_number}'
    data_dir2 = f'/home/zzhzhao/Model/wrfout/test-{test_number}-nolake'
    u1, v1, w1, qv1, h, ter, lats, lons, time, wrflist = load_wrfdata_windspeed(data_dir1) 

    NamCo_ll = (30.75, 90.9)
    # NamCo_xy = w.ll_to_xy(wrflist, NamCo_ll[0], NamCo_ll[1])

    start_point = w.CoordPair(lat=NamCo_ll[0]+0.5, lon=NamCo_ll[1]-0.5)
    end_point = w.CoordPair(lat=NamCo_ll[0]-0.5, lon=NamCo_ll[1]+0.5)
    
    u1_cross = w.vertcross(u1, h, wrfin=wrflist, start_point=start_point,
                       end_point=end_point, latlon=True, meta=True)
    v1_cross = w.vertcross(v1, h, wrfin=wrflist, start_point=start_point,
                       end_point=end_point, latlon=True, meta=True)
    w1_cross = w.vertcross(w1, h, wrfin=wrflist, start_point=start_point,
                       end_point=end_point, latlon=True,  meta=True)
    qv1_cross = w.vertcross(qv1, h, wrfin=wrflist, start_point=start_point,
                       end_point=end_point, latlon=True,  meta=True)
    ter_line = w.interpline(ter, wrfin=wrflist, start_point=start_point,
                      end_point=end_point)

#%%
    fig, ax = plt.subplots(figsize=(15,3), dpi=100)
    t = pd.Timestamp(f'20170610 12:00:00')
    x_labels, y = xy(u1_cross)   
    xs = np.arange(0, u1_cross.shape[-1], 1)

    # add NamCo station
    # label_color = 'white'
    # ax.plot(NamCo_ll[1], 4760, c=label_color, markersize=5, marker='^')
    # ax.text(NamCo_ll[1], 4580, 'NamCo', c=label_color, ha='center',weight='bold', fontsize=8)

    # QV
    var = fill_gap(qv1_cross.sel(Time=t)) * 1000
    p = ax.contourf(xs, y, var, levels=np.arange(3, 8+0.25, 0.25), cmap='Blues')
    cb = fig.colorbar(p, ax=ax, orientation='vertical', shrink=1, pad=0.01, aspect=15)
    cb.set_label('Mixing ratio / $g\cdot kg^{-1}$')

    # UW
    # span = 1
    # var1 = u1_cross.sel(Time=t)[::span, ::span]
    # var2 = w1_cross.sel(Time=t)[::span,::span] * 20
    # q = ax.quiver(x[::span], y[::span], var1, var2, color='k', scale=300)


    # x-label
    num_ticks = 10
    thin = int((len(xs) / num_ticks) + .5)
    ax.set_xticks(xs[::thin])
    ax.set_xticklabels(x_labels[::thin], rotation=45, fontsize=8)

    # add terrain
    ht_fill = ax.fill_between(xs, 0, w.to_np(ter_line),
                                facecolor="saddlebrown")

    # adjust 
    ax.set_ylim([4200, 6000])
    ax.set_ylabel('Height / m')
    # ax.set_xlabel('Longitude')
