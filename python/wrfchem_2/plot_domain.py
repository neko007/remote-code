#%%
from salem import geogrid_simulator
import matplotlib.pyplot as plt

if __name__ == '__main__':  
    geo_path = '/home/zzhzhao/code/python/wrfchem_2/data2/namelist.wps'
    g, maps = geogrid_simulator(geo_path)
    maps[0].set_rgb(natural_earth='lr')  # add a background image
    maps[0].visualize(title='Domains 1')
    plt.savefig('/home/zzhzhao/code/python/wrfchem_2/fig/domain.jpg', dpi=300)
# %%
