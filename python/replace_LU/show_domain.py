# %%
from salem import geogrid_simulator
import pylab as plt

fpath = '/home/zzhzhao/Model/tests/test-9.4/WPS/namelist.wps'
g, maps = geogrid_simulator(fpath)

maps[0].set_rgb(natural_earth='lr') 
maps[0].visualize(title='Domains 1 & 2')

plt.savefig('fig/domain.jpg', dpi=300)