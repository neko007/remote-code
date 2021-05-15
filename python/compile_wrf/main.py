import subprocess as sp 
import os 
import glob

if __name__ == '__main__':
    Model_dir = '/home/zzhzhao/Model'
    test_dir = os.path.join(Model_dir, 'test')
    source_name = 'original-WRF3.9.1'
    target_name = 'original-WRF3.9.1-YW_Lake-comp3'
    WRF_version = 'WRFV3'
    Modified_wrf_files_path = '/home/zzhzhao/code/python/compile_wrf/Modified_WRF_Files'

    wrf_compile_option = 15
    wps_compile_option = 17

    WRF_path = os.path.join(test_dir, target_name, WRF_version)
    WPS_path = os.path.join(test_dir, target_name, 'WPS')

    ### 复制相应文件
    print('>>>> Copy source file <<<<')
    os.chdir(test_dir)
    sp.run(f'cp -r {source_name} {target_name}', shell=True)

    ### 复制吴阳程序
    sp.run(f'cp {Modified_wrf_files_path}/YangW_Registry.EM_COMMON {WRF_path}/Registry/Registry.EM_COMMON', shell=True)
    sp.run(f'cp {Modified_wrf_files_path}/YangW_module_sf_lake.F {WRF_path}/phys/module_sf_lake.F', shell=True)
    sp.run(f'cp {Modified_wrf_files_path}/YangW_registry.lake {WRF_path}/Registry/registry.lake', shell=True)
    sp.run(f'cp {Modified_wrf_files_path}/YangW_registry.dimspec {WRF_path}/Registry/registry.dimspec', shell=True)
    sp.run(f'cp {Modified_wrf_files_path}/YangW_module_surface_driver.F {WRF_path}/phys/module_surface_driver.F', shell=True)
    sp.run(f'cp {Modified_wrf_files_path}/YangW_module_first_rk_step_part1.F {WRF_path}/dyn_em/module_first_rk_step_part1.F', shell=True)
    sp.run(f'cp {Modified_wrf_files_path}/YangW_module_physics_init.F {WRF_path}/phys/module_physics_init.F', shell=True)
    sp.run(f'cp {Modified_wrf_files_path}/YangW_start_em.F {WRF_path}/dyn_em/start_em.F', shell=True)

    ### 编译WRF
    os.chdir(WRF_path)
    print('>>>> Configure WRF <<<<')
    sp.run(f'echo {wrf_compile_option} | ./configure > log.configure', shell=True)
    print('>>>> Compile WRF <<<<')
    sp.run('./compile em_real >& log.compile', shell=True) 

    wrf_exefile = ['wrf.exe', 'real.exe', 'ndown.exe', 'tc.exe'].sort()
    if glob.glob(os.path.join(WRF_path, 'main, *.exe')).sort() == wrf_exefile:
        print('>>>> WRF Compile Success <<<<')
    else: 
        print('xxxx Error: WRF Compile Fail xxxx')
        os._exit(0)

    ### 编译WPS
    os.chdir(WPS_path)
    print('>>>> Configure WPS <<<<')
    sp.run(f'echo {wps_compile_option} | ./configure > log.configure', shell=True)
    print('>>>> Compile WPS <<<<')
    sp.run('./compile >& log.compile', shell=True) 

    wps_exefile = ['geogrid.exe', 'ungrib.exe', 'metgrid.exe'].sort()
    if glob.glob(os.path.join(WPS_path, '*.exe')).sort() == wps_exefile:
        print('>>>> WPS Compile Success <<<<')
    else: 
        print('xxxx Error: WPS Compile Fail xxxx')
        os._exit(0)

    print('**** Compile Success! ****')