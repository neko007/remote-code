import os
import f90nml 

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 需要手动修改的部分
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def revise_nml():
    # 基本路径
    model_dir = '/home/zzhzhao/Model'
    test_number = 'test-2'
    tests_dir = os.path.join(model_dir, 'tests')
    wps_dir = os.path.join(tests_dir, test_number, 'WPS')
    wrf_dir = os.path.join(tests_dir, test_number, 'WRF')
    wps_nml_name = 'namelist.wps'
    wrf_nml_name = 'namelist.input'

    # 再分析资料路径
    data_dir = '../data'
    data_file = 'GFS_*'

    # 模拟的起止时间
    y_start, m_start, d_start, h_start = 2020, 9, 1, 00
    y_end, m_end, d_end, h_end = 2020, 9, 3, 00
    run_days = 2
    run_hours = 0

    # 再分析资料时间间隔
    interval_hours = 3
    interval_seconds = interval_hours * 3600 

    # restart设置
    restart = False
    restart_interval = 1440 # minute

    # 文件输出时间间隔(min)、文件打包个数
    history_interval = [120, 120, 60]
    frames_per_outfile = [3, 3, 100]

    # 静态地形资料的路径
    geog_data_path = os.path.join(model_dir, 'Build_WRF/WPS_GEOG')
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

    # 积分步长
    time_step = dx / 1000 * 4

    # 投影
    map_proj = 'lambert'
    truelat1 = 30. 
    truelat2 = 60. 
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

    # 读取并修改namelist.wps
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

    # 读取并修改namelist.input 
    nml_wrf = f90nml.read(os.path.join(wrf_dir, 'run', wrf_nml_name))
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

    nml_wrf['domains']['max_dom'] = max_dom
    nml_wrf['domains']['time_step'] = time_step
    nml_wrf['domains']['e_we'] = e_wes
    nml_wrf['domains']['e_sn'] = e_sns
    nml_wrf['domains']['num_metgrid_levels'] = 34 # warning
    nml_wrf['domains']['dx'] = dxs 
    nml_wrf['domains']['dy'] = dys
    nml_wrf['domains']['parent_id'] = [0, 1, 2] # warning
    nml_wrf['domains']['i_parent_start'] = i_parent_starts
    nml_wrf['domains']['j_parent_start'] = j_parent_starts
    nml_wrf['domains']['parent_grid_ratio'] = parent_grid_ratios
    nml_wrf['domains']['parent_time_step_ratio'] = parent_grid_ratios
    nml_wrf['domains'].update({'sfcp_to_sfcp':True})

    return nml_wps, nml_wrf, wps_dir, wps_nml_name, wrf_dir, wrf_nml_name, data_dir, data_file, test_number
