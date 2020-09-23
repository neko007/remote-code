import os
import subprocess as sp
import numpy as np
from revise_nml2 import *

# 打印基本信息
def describe(test_number):
    print(f'**** {test_number} start ****')

# 编译运行wps
def cc_wps(sst_update):
    print('>>>> configure wps <<<<')
    sp.call('echo 17 | ./configure > log.configure', shell=True)
    print('>>>> compile wps <<<<')
    sp.run('./compile >& log.compile', shell=True)
    print('>>>> geogrid.exe <<<<')
    sp.run('./geogrid.exe >& log.geogrid', shell=True)
    # 链接数据
    print('>>>> link data <<<<')
    sp.run('ln -sf ungrib/Variable_Tables/Vtable.GFS Vtable', shell=True)
    sp.run('./link_grib.csh ' + os.path.join(data_dir, data_file), shell=True)
     # ungrid 解码
    print('>>>> ungrib.exe <<<<')
    sp.run('./ungrib.exe >& log.ungrid', shell=True)
    if sst_update == 1:
        print('>>>> prepare SST data <<<<')
        sp.run('ln -sf ungrib/Variable_Tables/Vtable.SST Vtable', shell=True)
        sp.run('./link_grib.csh ' + os.path.join(data_dir, sst_file), shell=True)
        # 修改namelist.wps
        nml_wps, nml_wrf = revise_nml(sst_update=1)
        rewrite_namelist(wps_dir, wps_nml_name, nml_wps, 1)
        print('>>>> ungrib.exe <<<<')
        sp.run('./ungrib.exe >& log.ungrid', shell=True)
    # 将气象场插值到地形上
    print('>>>> metgrid.exe <<<<')
    sp.run('./metgrid.exe >& log.metgrid', shell=True)


# 运行WRF
def cc_wrf():
    print('>>>> link met_em* <<<<')
    sp.run('ln -sf ../../WPS/met_em* .', shell=True)
    print('>>>> real.exe <<<<')
    sp.run('mpirun -np 1 ./real.exe', shell=True)
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
    # 修改namelist
    nml_wps, nml_wrf, sst_update = revise_nml(sst_update=0)

    describe(test_number)

    # 重新写namelist.wps
    rewrite_namelist(wps_dir, wps_nml_name, nml_wps, 1)

    # 配置编译运行WPS
    cc_wps(sst_update)

    # 重新写namlist.input
    rewrite_namelist(os.path.join(wrf_dir, 'run'), wrf_nml_name, nml_wrf, 2)

    # 运行WRF
    cc_wrf()
