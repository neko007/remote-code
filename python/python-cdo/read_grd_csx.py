#%%
from xgrads import CtlDescriptor, open_CtlDataset
import xarray as xr
import numpy as np

if __name__ == '__main__':
    year = 2009
    ctl_path = f'{year}.ctl'
    ctl = CtlDescriptor(file=ctl_path)
    ds = open_CtlDataset(ctl_path)
    rain = ds.r.squeeze()
    rain = xr.where(rain<0, np.nan, rain)
#%%
    import cartopy.crs as ccrs 
    import matplotlib.pyplot as plt 
    import pandas as pd 

    time_range = pd.date_range(f'{year}-06-01',f'{year}-07-01',freq='h')

    fig, ax = plt.subplots(dpi=100, subplot_kw={'projection':ccrs.PlateCarree()})
    rain.sel(time=time_range).sum(dim='time').plot.contourf(robust=True, ax=ax)

    # rain.sum(dim='time').plot.contourf(robust=True, transfrom=ccrs.PlateCarree(), ax=ax)

    ax.coastlines()

