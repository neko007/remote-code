#%%
from xgrads import CtlDescriptor, open_CtlDataset

if __name__ == '__main__':
    year = 2000
    ctl_path = 'CHN_PRCP_HOUR_MERG_DISPLAY_0.1deg.lnx.ctl'
    ctl = CtlDescriptor(file=ctl_path)
    ds = open_CtlDataset(ctl_path)
    rain = ds.crain.squeeze()
    rain = rain.where(rain>-999.)
#%%
    import cartopy.crs as ccrs 
    import matplotlib.pyplot as plt 
    import pandas as pd 

    time_range = pd.date_range(f'{year}-06-01',f'{year}-07-01',freq='h')

    fig, ax = plt.subplots(dpi=100, subplot_kw={'projection':ccrs.PlateCarree()})
    rain.sel(time=time_range).sum(dim='time').plot.contourf(robust=True, transfrom=ccrs.PlateCarree(), ax=ax)
    ax.coastlines()

