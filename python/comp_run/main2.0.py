'''
需要准备好WPS和WRF均编译完成的原始文件
'''

import os
import subprocess as sp
import numpy as np
from modify_nml2 import *

# 打印基本信息
def describe(test_number):
    print(f'**** {test_number} start ****')

# 编译运行wps
def run_wps(sst_flag):
    # 配置编译
    # print('>>>> configure wps <<<<')
    # sp.call('echo 17 | ./configure > log.configure', shell=True)
    # print('>>>> compile wps <<<<')
    # sp.run('./compile >& log.compile', shell=True)
    
    # 修改namelist.wps
    nml_wps = modify_wps_nml(0)
    # 复写 namelist.wps
    rewrite_namelist(wps_dir, wps_nml_name, nml_wps, 1)
    # geogrid 地形插值
    print('>>>> geogrid.exe <<<<')
    sp.run('./geogrid.exe >& geogrid.log', shell=True)
    # 链接数据
    print('>>>> link data <<<<')
    sp.run('ln -sf ungrib/Variable_Tables/Vtable.GFS Vtable', shell=True)
    sp.run('./link_grib.csh ' + os.path.join(data_dir, data_file), shell=True)
    # ungrid 解码
    print('>>>> ungrib.exe <<<<')
    sp.run('./ungrib.exe >& ungrib.log', shell=True)
    if sst_flag == 1:
        # 链接海温场
        print('>>>> prepare SST data <<<<')
        sp.run('ln -sf ungrib/Variable_Tables/Vtable.SST Vtable', shell=True)
        sp.run('./link_grib.csh ' + os.path.join(data_dir, sst_file), shell=True)
        # 修改namelist.wps
        nml_wps = modify_wps_nml(1)
        # 复写 namelist.wps
        rewrite_namelist(wps_dir, wps_nml_name, nml_wps, 1)
        # ungrid 解码
        print('>>>> ungrib.exe <<<<')
        sp.run('./ungrib.exe >& ungrib.log', shell=True)
    # metgrid 气象场插值
    print('>>>> metgrid.exe <<<<')
    sp.run('./metgrid.exe > metgrid.log', shell=True)


# 运行WRF
def run_wrf(sst_flag):
    # 修改namelist.input
    nml_wrf = modify_wrf_nml(sst_flag)
    # 复写namelist.input
    rewrite_namelist(os.path.join(wrf_dir, 'run'), wrf_nml_name, nml_wrf, 2)
    # 链接气象场
    print('>>>> link met_em* <<<<')
    sp.run('ln -sf ../../WPS/met_em* .', shell=True)
    # 生成初始场
    print('>>>> real.exe <<<<')
    sp.run('mpirun -np 1 ./real.exe', shell=True)
    # 运行wrf
    print('>>>> wrf.exe <<<<')
    sp.run(f'mpirun -np {core_num} ./wrf.exe', shell=True)

# 复写namelist
def rewrite_namelist(target_dir, nml_name, nml, flag):
    if(flag == 1):
        print('>>>> enter WPS <<<<')
    else:
        print('>>>> enter WRF/run <<<<')
    os.chdir(target_dir)
    os.remove(nml_name)
    print('>>>> rewrite namelist.wps <<<<')
    nml.write(nml_name)

if __name__ == '__main__':
    # 打印基本信息
    describe(test_number)

    # 运行WPS
    run_wps(sst_flag)

    # 运行WRF
    run_wrf(sst_flag)
