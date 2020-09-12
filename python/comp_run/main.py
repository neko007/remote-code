import os
import subprocess as sp
import numpy as np
from revise_nml import revise_nml

# 打印基本信息
def describe(test_number):
    print(f'**** {test_number} start ****')

# 编译运行wps
def cc_wps(data_dir, data_file):
    print('>>>> configure wps <<<<')
    sp.call('echo 17 | ./configure > log.configure', shell=True)
    print('>>>> compile wps <<<<')
    sp.run('./compile >& log.compile', shell=True)
    print('>>>> geogrid.exe <<<<')
    sp.run('./geogrid.exe >& log.geogrid', shell=True)
    print('>>>> link data <<<<')
    sp.run('./link_grib.csh ' + os.path.join(data_dir, data_file), shell=True)
    sp.run('ln -sf ungrib/Variable_Tables/Vtable.GFS Vtable', shell=True)
    print('>>>> ungrib.exe <<<<')
    sp.run('./ungrib.exe >& log.ungrid', shell=True)
    print('>>>> metgrid.exe <<<<')
    sp.run('./metgrid.exe >& log.metgrid', shell=True)

# 运行WRF
def cc_wrf():
    print('>>>> link met_em* <<<<')
    sp.run('ln -sf ../../WPS/met_em* .', shell=True)
    print('>>>> real.exe <<<<')
    sp.run('mpirun -np 1 ./real.exe', shell=True)
    print('>>>> wrf.exe <<<<')
    sp.run('mpirun -np 20 ./wrf.exe', shell=True)

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
    # # 修改namelist
    nml_wps, nml_wrf, wps_dir, wps_nml_name, wrf_dir, wrf_nml_name, data_dir, data_file, test_number \
         = revise_nml()

    describe(test_number)

    # 重新写namelist.wps
    rewrite_namelist(wps_dir, wps_nml_name, nml_wps, 1)

    # 配置编译运行WPS
    cc_wps(data_dir, data_file)

    # 重新写namlist.input
    rewrite_namelist(os.path.join(wrf_dir, 'run'), wrf_nml_name, nml_wrf, 2)

    # 运行WRF
    cc_wrf()

