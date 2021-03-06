import os
import f90nml 

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 需要经常修改的部分
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# 核数设置
core_num = 10

# 测试路径
test_number = '（记得修改）'

# 驱动场资料名称
data_file = '驱动场（记得修改）'
sst_file  = '海温场' # 有则修改

# 驱动场资料类型
Vtable_type = 'Vtable.GFS'

# 驱动场资料时间间隔
interval_hours     = 6
interval_seconds   = interval_hours * 3600 
# 驱动场资料层数
num_metgrid_levels =  32

# 模拟的起止时间
y_start, m_start, d_start, h_start = 2017, 6, 1, 00
y_end, m_end, d_end, h_end         = 2017, 7, 1, 00
run_days                           = 30
run_hours                          = 0

# restart设置
restart = False
restart_interval = 1440 # minute

# 文件输出时间间隔(min)、文件打包个数
history_interval   = 1440
frames_per_outfile = 1

# 辅助文件输出时间间隔(min)、文件打包个数
auxhist2_interval   = 180
frames_per_auxhist2 = 8

# 嵌套区域相关配置
max_dom            = 2
parent_ids         = [1, 1]
parent_grid_ratios = [1, 3]

# 每个区域的起止经纬度， 按照d01、d02...的顺序排列
lat_mins = [
    19, 
    29
    ]
lat_maxs = [
    43, 
    33
    ]
lon_mins = [
    77, 
    87
    ]
lon_maxs = [
    103, 
    93
    ]

# 父区域分辨率
dx = 15000
dy = 15000

### lake option
sf_lake_physics   = 1 # WRF-Lake开关
alternative_lswt  = 1 # 用WRF默认程式替代湖温（效果不佳）
alternative_lake  = 1 # 用相邻下垫面替换湖泊
use_lakedepth     = 1 # 是否使用静态地形数据中的湖深
md_lakedepth      = 1 # 替换纳木错湖泊深度（一般与alternative_lake相反）
lakedepth_default = 50. # 默认湖泊深度，use_lakedepth=0时生效

### zzz强制更改wrfinput_d0x湖温
lswt_init_flag       = False # 0: 关闭; 1:固定值替代湖温; 2:气温加减某值替代湖温
lswt_init            = 277.

### lake option based on Wuyang (part)
### WRF_version  = 'WRFV3' 时方可生效
tlake_init_flag      = True
tlake_init_value     = 276.05
eta_flag             = True # True: Wu; False: Gu
eta_scale_yw         = 0.575 # default: 0.6
eta_yw               = 0.1
diffusivity_flag     = True # True: Wu; False: Gu
diffusivity_index_yw = 4 # default: 1
tdmax                = 274.2
mixing_factor        = 100 # default: 40
mixing_factor_ked    = 40
roughness_flag       = True # True: Wu; False: Gu
coszen_flag          = True # True: Wu; False: Gu
albedo_flag          = True # True: Wu; False: Gu

### mountain
md_mountainHeight = 0 # 修改念青唐古拉山高度

# 引入海温场
sst_flag = 0

# chem
chem = 0

# **********
# 参数化方案
# **********

mp_physics         = 3
cu_physics         = 1
ra_lw_physics      = 1
ra_sw_physics      = 1
bl_pbl_physics     = 1
sf_sfclay_physics  = 1
sf_surface_physics = 2

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 不经常修改的部分
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# 路径
WRF_version    = 'WRFV3'
wps_nml_name   = 'namelist.wps'
wrf_nml_name   = 'namelist.input'
model_dir      = '/home/zzhzhao/Model'
tests_dir      = os.path.join(model_dir, 'tests')
root_dir       = os.path.join(tests_dir, test_number)
comp_run_dir   = os.path.join(root_dir, 'comp_run_new')
wps_dir        = os.path.join(root_dir, 'WPS')
wrf_dir        = os.path.join(root_dir, WRF_version) # 针对WRF3.9.1进行修改
data_dir       = '../data' # 驱动场资料路径
geog_data_path = os.path.join(model_dir, 'Build_WRF/WPS_GEOG') # 静态地形资料的路径

# 投影
map_proj = 'lambert'
truelat1 = 30. 
truelat2 = 60. 

# 地形资料分辨率
geog_data_res = 'modis_lake+default'

# ungrib
prefix = 'FILE'

# 积分步长
time_step = dx / 1000 * 5

# wrfout输出路径 （暂时不用）
wrfout_path = f'/home/zzhzhao/Model/wrfout/{test_number}'
wrfout_file = 'wrfout_d<domain>_<date>.nc'

# 辅助文件（自定义输出）
auxhist2_file    = 'auxhist2_d<domain>_<date>.nc'
io_form_auxhist2 = 2
iofile_name      = 'iofile'

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 自动计算部分
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# 计算起止时间
start_date = '%d-%02d-%02d_%02d:00:00' % (y_start, m_start, d_start, h_start)
end_date   = '%d-%02d-%02d_%02d:00:00' % (y_end, m_end, d_end, h_end)

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

def modify_wps_nml():
    '''
    读取并修改namelist.wps
    '''
    nml_wps                              = f90nml.read(os.path.join(wps_dir, wps_nml_name))
    nml_wps['share']['max_dom']          = max_dom
    nml_wps['share']['start_date']       = [start_date for i in range(max_dom)]
    nml_wps['share']['end_date']         = [end_date for i in range(max_dom)] 
    nml_wps['share']['interval_seconds'] = interval_seconds

    nml_wps['geogrid']['parent_id']         = parent_ids
    nml_wps['geogrid']['parent_grid_ratio'] = parent_grid_ratios
    nml_wps['geogrid']['i_parent_start']    = i_parent_starts
    nml_wps['geogrid']['j_parent_start']    = j_parent_starts
    nml_wps['geogrid']['e_we']              = e_wes 
    nml_wps['geogrid']['e_sn']              = e_sns 
    nml_wps['geogrid']['dx']                = dx
    nml_wps['geogrid']['dy']                = dy 
    nml_wps['geogrid']['map_proj']          = map_proj
    nml_wps['geogrid']['ref_lat']           = ref_lat 
    nml_wps['geogrid']['ref_lon']           = ref_lon
    nml_wps['geogrid']['truelat1']          = truelat1
    nml_wps['geogrid']['truelat2']          = truelat2
    nml_wps['geogrid']['stand_lon']         = stand_lon 
    nml_wps['geogrid']['geog_data_path']    = geog_data_path
    nml_wps['geogrid']['geog_data_res']     = [geog_data_res] * max_dom

    nml_wps['ungrib']['prefix'] = prefix
    if sst_flag == 1:
        nml_wps['ungrib']['prefix']   = 'SST'
        nml_wps['metgrid']['fg_name'] = ['FILE', 'SST']
    if alternative_lswt == 1:
        nml_wps['metgrid']['constants_name'] = 'TAVGSFC'
    
    return nml_wps
    
def modify_wrf_nml():
    '''
    读取并修改namelist.input 
    '''
    nml_wrf = f90nml.read(os.path.join(wrf_dir, 'run', wrf_nml_name))

    ### 时间控制
    nml_wrf['time_control']['run_days']           = run_days
    nml_wrf['time_control']['run_hours']          = run_hours
    nml_wrf['time_control']['start_year']         = [y_start for i in range(max_dom)]
    nml_wrf['time_control']['start_month']        = [m_start for i in range(max_dom)]
    nml_wrf['time_control']['start_day']          = [d_start for i in range(max_dom)]
    nml_wrf['time_control']['start_hour']         = [h_start for i in range(max_dom)]
    nml_wrf['time_control']['end_year']           = [y_end for i in range(max_dom)]
    nml_wrf['time_control']['end_month']          = [m_end for i in range(max_dom)]
    nml_wrf['time_control']['end_day']            = [d_end for i in range(max_dom)]
    nml_wrf['time_control']['end_hour']           = [h_end for i in range(max_dom)]
    nml_wrf['time_control']['interval_seconds']   = interval_seconds
    nml_wrf['time_control']['history_interval']   = [history_interval] * max_dom
    nml_wrf['time_control']['frames_per_outfile'] = [frames_per_outfile] * max_dom
    nml_wrf['time_control']['restart']            = restart
    nml_wrf['time_control']['restart_interval']   = restart_interval
    # nml_wrf['time_control'].update({'history_outname':os.path.join(wrfout_path, wrfout_file)})
    # 自定义输出
    nml_wrf['time_control'].update({'iofields_filename':[f'{iofile_name}_d0{i+1}.txt' for i in range(max_dom)]})
    nml_wrf['time_control'].update({'ignore_iofields_warning':True})
    nml_wrf['time_control'].update({'auxhist2_outname':os.path.join(wrfout_path, auxhist2_file)})
    nml_wrf['time_control'].update({'auxhist2_interval':[auxhist2_interval] * max_dom})
    nml_wrf['time_control'].update({'frames_per_auxhist2':[frames_per_auxhist2] * max_dom})
    nml_wrf['time_control'].update({'io_form_auxhist2':io_form_auxhist2})

    ### 区域设置
    nml_wrf['domains']['max_dom']            = max_dom
    nml_wrf['domains']['time_step']          = time_step
    nml_wrf['domains']['e_we']               = e_wes
    nml_wrf['domains']['e_sn']               = e_sns
    nml_wrf['domains']['num_metgrid_levels'] = num_metgrid_levels
    nml_wrf['domains']['dx']                 = dxs 
    nml_wrf['domains']['dy']                 = dys
    nml_wrf['domains']['parent_id']          = [0, 1, 2] # warning
    nml_wrf['domains']['i_parent_start']     = i_parent_starts
    nml_wrf['domains']['j_parent_start']     = j_parent_starts
    nml_wrf['domains']['parent_grid_ratio']  = parent_grid_ratios
    nml_wrf['domains']['parent_time_step_ratio'] = parent_grid_ratios
    nml_wrf['domains'].update({'sfcp_to_sfcp':True})

    ### 参数化方案
    nml_wrf['physics']['mp_physics']         = [mp_physics] * max_dom
    nml_wrf['physics']['cu_physics']         = [cu_physics] * max_dom
    nml_wrf['physics']['ra_lw_physics']      = [ra_lw_physics] * max_dom
    nml_wrf['physics']['ra_sw_physics']      = [ra_sw_physics] * max_dom
    nml_wrf['physics']['bl_pbl_physics']     = [bl_pbl_physics] * max_dom
    nml_wrf['physics']['sf_sfclay_physics']  = [sf_sfclay_physics] * max_dom
    nml_wrf['physics']['sf_surface_physics'] = [sf_surface_physics] * max_dom

    if sst_flag == 1:
        nml_wrf['time_control'].update({'auxinput4_inname':'wrflowinp_d<domain>', 'auxinput4_interval':360, 'io_form_auxinput4':2})
        nml_wrf['physics'].update({'sst_update':sst_update})
    
    ### lake
    if sf_lake_physics == 1:
        nml_wrf['physics'].update({'sf_lake_physics':[1] * max_dom})
        nml_wrf['physics'].update({'lakedepth_default':[lakedepth_default] * max_dom})
        nml_wrf['physics'].update({'lake_min_elev':[5] * max_dom})
        nml_wrf['physics'].update({'use_lakedepth':[use_lakedepth] * max_dom})
        if WRF_version == 'WRFV3': # 仅在WRF3.9.1中生效
            nml_wrf['physics'].update({'tlake_init_flag':[tlake_init_flag] * max_dom})
            nml_wrf['physics'].update({'tlake_init_value':[tlake_init_value] * max_dom})
            nml_wrf['physics'].update({'eta_flag':[eta_flag] * max_dom})
            nml_wrf['physics'].update({'eta_scale_yw':[eta_scale_yw] * max_dom})
            nml_wrf['physics'].update({'eta_yw':[eta_yw] * max_dom})
            nml_wrf['physics'].update({'diffusivity_flag':[diffusivity_flag] * max_dom})
            nml_wrf['physics'].update({'diffusivity_index_yw':[diffusivity_index_yw] * max_dom})
            nml_wrf['physics'].update({'tdmax':[tdmax] * max_dom})
            nml_wrf['physics'].update({'mixing_factor':[mixing_factor] * max_dom})
            nml_wrf['physics'].update({'mixing_factor_ked':[mixing_factor_ked] * max_dom})
            nml_wrf['physics'].update({'roughness_flag':[roughness_flag] * max_dom})
            nml_wrf['physics'].update({'coszen_flag':[coszen_flag] * max_dom})
            nml_wrf['physics'].update({'albedo_flag':[albedo_flag] * max_dom})

    ### 运行WRF3.9.1版本需要添加的地方，否则会在./wrf.exe时报错
    if WRF_version == 'WRFV3':
        nml_wrf['dynamics'].update({'max_rot_angle_gwd':100})

    ### chem
    if chem == 1:
        # 将namelist.input.chem的&chem导入到namelist.input
        nml_wrf_chem    = f90nml.read(os.path.join(wrf_dir, 'test', 'em_real', 'namelist.input.chem'))
        nml_wrf['chem'] = nml_wrf_chem['chem']

        nml_wrf['chem']['chem_opt']           = [300] * max_dom
        nml_wrf['chem']['kemit']              = [1] * max_dom 
        nml_wrf['chem']['chemdt']             = [20] * max_dom 
        nml_wrf['chem']['emiss_inpt_opt']     = [0] * max_dom
        nml_wrf['chem']['emiss_opt']          = [0] * max_dom
        nml_wrf['chem']['io_style_emissions'] = [0] * max_dom
        nml_wrf['chem']['phot_opt']           = [0] * max_dom 
        nml_wrf['chem']['bio_emiss_opt']      = [0] * max_dom
        nml_wrf['chem'].update({'ne_area'     : [200] * max_dom})
        nml_wrf['chem'].update({'depo_fact'   : [0.25] * max_dom})
        nml_wrf['chem']['dust_opt']           = [1] * max_dom
        nml_wrf['chem'].update({'aer_op_opt'  : [1] * max_dom})
        nml_wrf['chem'].update({'opt_pars_out': [1] * max_dom})
        nml_wrf['domains']['e_vert']          = [33] * max_dom # namelist.input.chem的e_vert默认是20


    return nml_wrf