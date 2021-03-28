import os
import f90nml 

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 需要手动修改的部分
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# 运行设置
core_num = 16

# 基本路径
model_dir = '/home/zzhzhao/Model'
test_number = ''

# 驱动场资料名称
data_file = '驱动场'
sst_file = '海温场'

# 驱动场资料时间间隔
interval_hours = 6
interval_seconds = interval_hours * 3600 

# 驱动场资料层数
num_metgrid_levels = 32

# 模拟的起止时间
y_start, m_start, d_start, h_start = 2020, 9, 1, 00
y_end, m_end, d_end, h_end = 2020, 9, 3, 00
run_days = 2
run_hours = 0

# restart设置
restart = False
restart_interval = 1440 # minute

# 文件输出时间间隔(min)、文件打包个数
history_interval = [60, 60, 60]
frames_per_outfile = [6, 6, 100]

# 嵌套区域相关配置
max_dom = 1
parent_ids = [1, 1]
parent_grid_ratios = [1, 3]

# 每个区域的起止经纬度， 按照d01、d02...的顺序排列
lat_mins = [
    20, 
    30
    ]
lat_maxs = [
    55, 
    45
    ]
lon_mins = [
    80, 
    105
    ]
lon_maxs = [
    145, 
    125
    ]

# 父区域分辨率
dx = 30000
dy = 30000

# 引入海温场
sst_flag = 0

# lake
sf_lake_physics = 1

# chem
chem = 0

# **********
# 参数化方案
# **********

mp_physics = 3
cu_physics = 1
ra_lw_physics = 1
ra_sw_physics = 1
bl_pbl_physics = 1
sf_sfclay_physics = 1
sf_surface_physics = 2

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 不经常修改的部分
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# 路径
wps_nml_name = 'namelist.wps'
wrf_nml_name = 'namelist.input'
tests_dir = os.path.join(model_dir, 'tests')
wps_dir = os.path.join(tests_dir, test_number, 'WPS')
wrf_dir = os.path.join(tests_dir, test_number, 'WRF')
data_dir = '../data' # 驱动场资料路径
geog_data_path = os.path.join(model_dir, 'Build_WRF/WPS_GEOG') # 静态地形资料的路径

# 投影
map_proj = 'lambert'
truelat1 = 30. 
truelat2 = 60. 

# 地形资料分辨率
geog_data_res = 'modis_lake+default'

# ungrib
prefix = 'FILE'

# 驱动场资料类型
Vtable_type = 'Vtable.GFS'

# 积分步长
time_step = dx / 1000 * 5

# wrfout输出路径
wrfout_path = f'/home/zzhzhao/Model/wrfout/{test_number}'
wrfout_file = 'wrfout_d<domain>_<date>.nc'

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 自动计算部分
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# 计算起止时间
start_date = '%d-%02d-%02d_%02d:00:00' % (y_start, m_start, d_start, h_start)
end_date = '%d-%02d-%02d_%02d:00:00' % (y_end, m_end, d_end, h_end)

# 计算网格配置
dxs = [] 
dys = []
e_wes = []
e_sns = []
i_parent_starts = []
j_parent_starts = []
for i in range(max_dom):
    scale_cof = parent_grid_ratios[i] * parent_grid_ratios[parent_ids[i]-1]  
    e_we = int((lon_maxs[i] - lon_mins[i]) * 110 * scale_cof / (dx * 1e-3))
    e_sn = int((lat_maxs[i] - lat_mins[i]) * 110 * scale_cof / (dx * 1e-3))
    if i+1 == 1:
        i_parent_start = 1 
        j_parent_start = 1
    else:
        k = parent_ids[i] - 1
        scale_cof_parent = parent_grid_ratios[k] * parent_grid_ratios[parent_ids[k]-1]
        i_parent_start = int((lon_mins[i] - lon_mins[k]) * 110 * scale_cof_parent / (dx * 1e-3))
        j_parent_start = int((lat_mins[i] - lat_mins[k]) * 110 * scale_cof_parent / (dx * 1e-3))
        while(e_we%3 != 1):
            e_we += 1
        while(e_sn%3 != 1):
            e_sn += 1
    dxs.append(int(dx/scale_cof))
    dys.append(int(dy/scale_cof))
    e_wes.append(e_we)
    e_sns.append(e_sn)
    i_parent_starts.append(i_parent_start)
    j_parent_starts.append(j_parent_start)

ref_lat = (lat_maxs[0] + lat_mins[0]) / 2.
ref_lon = (lon_maxs[0] + lon_mins[0]) / 2. 
stand_lon = ref_lon

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

def modify_wps_nml(sst_update):
    '''
    读取并修改namelist.wps
    '''
    nml_wps = f90nml.read(os.path.join(wps_dir, wps_nml_name))
    nml_wps['share']['max_dom'] = max_dom
    nml_wps['share']['start_date'] = [start_date for i in range(max_dom)]
    nml_wps['share']['end_date'] = [end_date for i in range(max_dom)] 
    nml_wps['share']['interval_seconds'] = interval_seconds

    nml_wps['geogrid']['parent_id'] = parent_ids
    nml_wps['geogrid']['parent_grid_ratio'] = parent_grid_ratios
    nml_wps['geogrid']['i_parent_start'] = i_parent_starts
    nml_wps['geogrid']['j_parent_start'] = j_parent_starts
    nml_wps['geogrid']['e_we'] = e_wes 
    nml_wps['geogrid']['e_sn'] = e_sns 
    nml_wps['geogrid']['dx'] = dx
    nml_wps['geogrid']['dy'] = dy 
    nml_wps['geogrid']['map_proj'] = map_proj
    nml_wps['geogrid']['ref_lat'] = ref_lat 
    nml_wps['geogrid']['ref_lon'] = ref_lon
    nml_wps['geogrid']['truelat1'] = truelat1
    nml_wps['geogrid']['truelat2'] = truelat2
    nml_wps['geogrid']['stand_lon'] = stand_lon 
    nml_wps['geogrid']['geog_data_path'] = geog_data_path
    nml_wps['geogrid']['geog_data_res'] = [geog_data_res] * max_dom

    nml_wps['ungrib']['prefix'] = prefix
    if sst_update == 1:
        nml_wps['ungrib']['prefix'] = 'SST'
        nml_wps['metgrid']['fg_name'] = ['FILE', 'SST']
    
    return nml_wps
    
def modify_wrf_nml(sst_update):
    '''
    读取并修改namelist.input 
    '''
    nml_wrf = f90nml.read(os.path.join(wrf_dir, 'run', wrf_nml_name))

    # 时间控制
    nml_wrf['time_control']['run_days'] = run_days
    nml_wrf['time_control']['run_hours'] = run_hours
    nml_wrf['time_control']['start_year'] = [y_start for i in range(max_dom)]
    nml_wrf['time_control']['start_month'] = [m_start for i in range(max_dom)]
    nml_wrf['time_control']['start_day'] = [d_start for i in range(max_dom)]
    nml_wrf['time_control']['start_hour'] = [h_start for i in range(max_dom)]
    nml_wrf['time_control']['end_year'] = [y_end for i in range(max_dom)]
    nml_wrf['time_control']['end_month'] = [m_end for i in range(max_dom)]
    nml_wrf['time_control']['end_day'] = [d_end for i in range(max_dom)]
    nml_wrf['time_control']['end_hour'] = [h_end for i in range(max_dom)]
    nml_wrf['time_control']['interval_seconds'] = interval_seconds
    nml_wrf['time_control']['history_interval'] = history_interval
    nml_wrf['time_control']['frames_per_outfile'] = frames_per_outfile
    nml_wrf['time_control']['restart'] = restart
    nml_wrf['time_control']['restart_interval'] = restart_interval
    nml_wrf['time_control'].update({'history_outname':os.path.join(wrfout_path, wrfout_file)})

    # 区域设置
    nml_wrf['domains']['max_dom'] = max_dom
    nml_wrf['domains']['time_step'] = time_step
    nml_wrf['domains']['e_we'] = e_wes
    nml_wrf['domains']['e_sn'] = e_sns
    nml_wrf['domains']['num_metgrid_levels'] = num_metgrid_levels
    nml_wrf['domains']['dx'] = dxs 
    nml_wrf['domains']['dy'] = dys
    nml_wrf['domains']['parent_id'] = [0, 1, 2] # warning
    nml_wrf['domains']['i_parent_start'] = i_parent_starts
    nml_wrf['domains']['j_parent_start'] = j_parent_starts
    nml_wrf['domains']['parent_grid_ratio'] = parent_grid_ratios
    nml_wrf['domains']['parent_time_step_ratio'] = parent_grid_ratios
    nml_wrf['domains'].update({'sfcp_to_sfcp':True})

    # 参数化方案
    nml_wrf['physics']['mp_physics'] = [mp_physics] * max_dom
    nml_wrf['physics']['cu_physics'] = [cu_physics] * max_dom
    nml_wrf['physics']['ra_lw_physics'] = [ra_lw_physics] * max_dom
    nml_wrf['physics']['ra_sw_physics'] = [ra_sw_physics] * max_dom
    nml_wrf['physics']['bl_pbl_physics'] = [bl_pbl_physics] * max_dom
    nml_wrf['physics']['sf_sfclay_physics'] = [sf_sfclay_physics] * max_dom
    nml_wrf['physics']['sf_surface_physics'] = [sf_surface_physics] * max_dom

    if sst_update == 1:
        nml_wrf['time_control'].update({'auxinput4_inname':'wrflowinp_d<domain>', 'auxinput4_interval':360, 'io_form_auxinput4':2})
        nml_wrf['physics'].update({'sst_update':sst_update})
    
    # lake
    if sf_lake_physics == 1:
        nml_wrf['physics'].update({'sf_lake_physics':[1] * max_dom})
        nml_wrf['physics'].update({'lakedepth_default':[50] * max_dom})
        nml_wrf['physics'].update({'lake_min_elev':[5] * max_dom})
        nml_wrf['physics'].update({'use_lakedepth':[0] * max_dom})

    # chem
    if chem == 1:
        # 将namelist.input.chem的&chem导入到namelist.input
        nml_wrf_chem = f90nml.read(os.path.join(wrf_dir, 'test', 'em_real', 'namelist.input.chem'))
        nml_wrf['chem'] = nml_wrf_chem['chem']

        nml_wrf['chem']['chem_opt'] = [300] * max_dom
        nml_wrf['chem']['kemit'] = [1] * max_dom 
        nml_wrf['chem']['chemdt'] = [20] * max_dom 
        nml_wrf['chem']['emiss_inpt_opt'] = [0] * max_dom
        nml_wrf['chem']['emiss_opt'] = [0] * max_dom
        nml_wrf['chem']['io_style_emissions'] = [0] * max_dom
        nml_wrf['chem']['phot_opt'] = [0] * max_dom 
        nml_wrf['chem']['bio_emiss_opt'] = [0] * max_dom
        nml_wrf['chem'].update({'ne_area':[200] * max_dom})
        nml_wrf['chem'].update({'depo_fact':[0.25] * max_dom})
        nml_wrf['chem']['dust_opt'] = [1] * max_dom
        nml_wrf['chem'].update({'aer_op_opt':[1] * max_dom})
        nml_wrf['chem'].update({'opt_pars_out':[1] * max_dom})
        nml_wrf['domains']['e_vert'] = [33] * max_dom # namelist.input.chem的e_vert默认是20


    return nml_wrf