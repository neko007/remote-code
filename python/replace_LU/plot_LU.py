#%%
import sys
sys.path.append('/home/zzhzhao/code/python/replace_LU')
from Module import *
import salem

lu_name_color = np.array([
    ['Evergreen Needleleaf Forest', 'darkgreen'],
    ['Evergreen Broadleaf Forest', 'forestgreen'],
    ['Deciduous Needleleaf Forest', 'yellowgreen'],
    ['Deciduous Broadleaf Forest', 'salmon'],
    ['Mixed Forests', 'darkolivegreen'],
    ['Closed Shrublands', 'yellow'],
    ['Open Shrublands', 'lawngreen'],
    ['Woody Savannas', 'olivedrab'],
    ['Savannas', 'coral'],
    ['Grasslands', 'limegreen'],
    ['Permanent Wetlands', 'turquoise'],
    ['Croplands', 'gold'],
    ['Urban and Built-Up', 'red'],
    ['Cropland/Natural Vegetation Mosaic', 'orange'],
    ['Snow and Ice', 'white'],
    ['Barren or Sparsely Vegetated', 'seashell'],
    ['Water', 'dodgerblue'],
    ['Wooded Tundra', 'plum'],
    ['Mixed Tundra', 'violet'],
    ['Barren Tundra', 'purple'],
    ['Lake', 'cyan'],
    ])
lu_name = [item[0] for item in lu_name_color]
lu_cmap = [item[1] for item in lu_name_color]

if __name__ == '__main__':
    geo_file = '/home/zzhzhao/Model/tests/test-9.4/WPS/geo_em.d02.nc'
    fig_dir = '/home/zzhzhao/code/python/replace_LU/fig'
    # geo_file = 'geo_em.d02.nc'
    geo = salem.open_xr_dataset(geo_file)
    lu = geo.LU_INDEX.isel(Time=0)
    lats, lons = w.latlon_coords(w.getvar(geo._file_obj.ds, 'LU_INDEX'))

    
    from matplotlib.colors import LinearSegmentedColormap
    cmap_modis = LinearSegmentedColormap.from_list('modis_21_landuse', lu_cmap)
    proj = ccrs.PlateCarree()

    fig, ax = plt.subplots(figsize=(13,9), subplot_kw={'projection':proj})
    ax = add_artist(ax, proj)
    ax = add_NamCo(ax)
    
    ax.pcolor(lons, lats, lu, cmap=cmap_modis, transform=proj)

    import matplotlib.patches as mpatches
    plt.legend([mpatches.Patch(color=i) for i in lu_cmap], lu_name, frameon=False, bbox_to_anchor=(0, -0.05), loc='upper left', ncol=4)

    fig.savefig(f'{fig_dir}/landuse_new.jpg', dpi=300, bbox_inches='tight', pad_inches=0.1)
