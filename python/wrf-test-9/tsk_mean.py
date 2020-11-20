#%%
from Module import *
from prec_sum import add_artist, add_NamCo
import gdal

def load_wrfdata_TSK(data_dir):
    wrf_files = [f for f in os.listdir(data_dir) if f[9]=='2']
    wrflist = [Dataset(os.path.join(data_dir, wrf_file)) for wrf_file in wrf_files]

    tsk = w.getvar(wrflist, 'TSK', timeidx=w.ALL_TIMES, method='cat') - 273.15
    lats, lons = w.latlon_coords(tsk)
    time = tsk.Time.to_index() 

    return tsk, lats, lons, time

class DDataset:
    def __init__(self, in_file):
        self.in_file = in_file  # Tiff或者ENVI文件

        dataset = gdal.Open(self.in_file)
        self.XSize = dataset.RasterXSize  # 网格的X轴像素数量
        self.YSize = dataset.RasterYSize  # 网格的Y轴像素数量
        self.GeoTransform = dataset.GetGeoTransform()  # 投影转换信息
        self.ProjectionInfo = dataset.GetProjection()  # 投影信息
    
    def get_data(self, band):
        dataset = gdal.Open(self.in_file)
        band = dataset.GetRasterBand(band)
        data = band.ReadAsArray()
        return data
    
    def get_lon_lat(self):
        gtf = self.GeoTransform
        x_range = range(0, self.XSize)
        y_range = range(0, self.YSize)
        x, y = np.meshgrid(x_range, y_range)
        lon = gtf[0] + x * gtf[1] + y * gtf[2]
        lat = gtf[3] + x * gtf[4] + y * gtf[5]
        return lon, lat

def get_MODIS(modis_file):
    dataset = DDataset(modis_file)
    band = 1
    modis = dataset.get_data(band)  # 获取第一个通道的数据
    lon, lat = dataset.get_lon_lat()  # 获取经纬度信息
    modis = xr.DataArray(modis, coords=[lat[:,0], lon[0,:]], dims=['lat', 'lon'])
    modis = modis.sel(lat=slice(34, 28), lon=slice(86,94))
    modis = xr.where(modis>1e4, np.nan, modis)
    return modis
    
if __name__ == '__main__':
    data_dir1 = '/home/zzhzhao/Model/wrfout/test-9.4'
    data_dir2 = '/home/zzhzhao/Model/wrfout/test-9.4-nolake'
    tsk1, lats, lons, time = load_wrfdata_TSK(data_dir1)
    tsk2, lats, lons, time = load_wrfdata_TSK(data_dir2)

    ### MODIS资料
    modis_file_day = 'data/MOD_LSTD_M_2017-06-01_rgb_3600x1800.FLOAT.TIFF'
    modis_file_night = 'data/MOD_LSTN_M_2017-06-01_rgb_3600x1800.FLOAT.TIFF'
    modis_day = get_MODIS(modis_file_day)
    modis_night = get_MODIS(modis_file_night)
    modis = (modis_day + modis_night) / 2

    ### 平均表面温度
    # tsk = tsk.sel(Time=pd.date_range('2017-6-01 18:00:00', periods=30, freq='1D'))
    tsk1_mean = tsk1.mean(dim='Time')
    tsk2_mean = tsk2.mean(dim='Time')
#%%
    ### 表面温度分布
    proj = ccrs.PlateCarree()
    labels = ['MODIS', 'WRF', 'WRF-LakeTurnOff']
    fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(10,10), subplot_kw={'projection':proj})
    fig.subplots_adjust(hspace=0.01)
    for i in range(2):
        for j in range(2):
            axes[i, j] = add_artist(axes[i, j], proj)
            axes[i, j] = add_NamCo(axes[i, j])
    c = axes[0][0].pcolor(modis.lon, modis.lat, modis, vmin=-10, vmax=20, cmap=cmaps.WhiteBlueGreenYellowRed, transform=proj)
    # c = axes[0][0].pcolor(era5_skt.longitude, era5_skt.latitude, era5_skt, vmin=-10, vmax=20, cmap=cmaps.WhiteBlueGreenYellowRed, transform=proj)
    axes[0][0].set_title(labels[0], fontsize=14, weight='bold')
    for j, tsk_mean in enumerate([tsk1_mean, tsk2_mean]):
        c = axes[1][j].pcolor(lons, lats, tsk_mean, vmin=-10, vmax=20, cmap=cmaps.WhiteBlueGreenYellowRed, transform=proj)
        axes[i][j].set_title(labels[j+1], fontsize=14, weight='bold')
    cb = fig.colorbar(c, ax=axes, orientation='horizontal', pad=0.05, shrink=0.9, aspect=35)
    cb.set_label('Skin temperature / $\mathrm{^\circ C}$', fontsize=14)
    axes[0][1].set_visible(False)

    # fig.savefig('fig-test-9.4/tsk.jpg', dpi=300)