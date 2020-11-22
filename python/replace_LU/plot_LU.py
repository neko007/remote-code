#%%
from Module import *
import salem

lu_name = ['Evergreen Needleleaf Forest',
    'Evergreen Broadleaf Forest',
    'Deciduous Needleleaf Forest',
    'Deciduous Broadleaf Forest',
    'Mixed Forests',
    'Closed Shrublands',
    'Open Shrublands',
    'Woody Savannas',
    'Savannas',
    'Grasslands',
    'Permanent Wetlands',
    'Croplands',
    'Urban and Built-Up',
    'Cropland/Natural Vegetation Mosaic',
    'Snow and Ice',
    'Barren or Sparsely Vegetated',
    'Water',
    'Wooded Tundra',
    'Mixed Tundra',
    'Barren Tundra',
    'Lake',
    ]
lu_name.reverse()

if __name__ == '__main__':
    geo_file = '/home/zzhzhao/Model/tests/test-9.4/WPS/geo_em.d02.nc'
    geo = salem.open_xr_dataset(geo_file)
    lu = geo.LU_INDEX.isel(Time=0)
    lats, lons = w.latlon_coords(w.getvar(geo._file_obj.ds, 'LU_INDEX'))
    
    lu = lu.salem.roi(shape=shp)
    lu = np.abs(lu - 22)

    proj = ccrs.PlateCarree()
    fig, ax = plt.subplots(figsize=(10,8), subplot_kw={'projection':proj})
    ax = add_artist(ax, proj)
    ax = add_NamCo(ax)

    p = ax.pcolor(lons, lats, lu, vmin=0.5, vmax=21.5, cmap=cmaps.t2m_29lev[:21], transform=proj)
    cb = fig.colorbar(p, orientation='vertical', ticks=np.arange(1,22), shrink=0.62, pad=0.03, aspect=25)
    cb.ax.set_yticklabels(lu_name, rotation=0)
