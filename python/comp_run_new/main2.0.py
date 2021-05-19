'''
需要准备好WPS和WRF均编译完成的原始文件
'''

import os
import subprocess as sp
import numpy as np
from modify_nml2 import *

# 打印基本信息
def describe():
    print(f'**** {test_number} start ****')

# 编译运行wps
def run_wps():
    # 修改namelist.wps
    nml_wps = modify_wps_nml()
    # 复写 namelist.wps
    rewrite_namelist(wps_dir, wps_nml_name, nml_wps, 1)
    # geogrid 地形插值
    print('>>>> geogrid.exe <<<<')
    sp.run('./geogrid.exe >& geogrid.log', shell=True)
    # 替换湖泊为相邻下垫面
    remove_lake(alternative_lake)
    # 替换湖泊深度 
    modify_lakedepth(md_lakedepth, alternative_lake)
    # 去除念青唐古拉山
    modify_mountainHeight(md_mountainHeight)
    # 链接数据
    print('>>>> link data <<<<')
    sp.run(f'ln -sf ungrib/Variable_Tables/{Vtable_type} Vtable', shell=True)
    sp.run('./link_grib.csh ' + os.path.join(data_dir, data_file), shell=True)
    # ungrid 解码
    print('>>>> ungrib.exe <<<<')
    sp.run('./ungrib.exe >& ungrib.log', shell=True)
    # 更新海温场
    sst_update(sst_flag)
    # 初始湖温替换
    alter_lswt(alternative_lswt)
    # metgrid 气象场插值
    print('>>>> metgrid.exe <<<<')
    sp.run('./metgrid.exe > metgrid.log', shell=True)

# 运行WRF
def run_wrf():
    # 修改namelist.input
    nml_wrf = modify_wrf_nml()
    # 复制iofield文件到WRF/run
    copy_iofiles()
    # 复写namelist.input
    rewrite_namelist(os.path.join(wrf_dir, 'run'), wrf_nml_name, nml_wrf, 2)
    # 链接气象场    
    print('>>>> link met_em* <<<<')
    sp.run('ln -sf ../../WPS/met_em* .', shell=True)
    # 生成初始场
    print('>>>> real.exe <<<<')
    sp.run('mpirun -np 1 ./real.exe', shell=True)
    # 修改初始湖温
    modify_lswt(lswt_init_flag)
    # 新建wrfout输出文件夹
    if os.path.exists(wrfout_path):
        print('*** Warning: path exists ***')
    else:
        print('>>>> mkdir wrfout_path <<<<')
        sp.run(f'mkdir {wrfout_path}', shell=True)
    # 运行wrf
    print('>>>> wrf.exe <<<<')
    sp.run(f'nohup mpirun -np {core_num} ./wrf.exe 2>&1 &', shell=True)

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

def sst_update(sst_flag):
    '''
    用额外的海温场更新初始/边界条件
    '''
    if sst_flag == 1:
        # 链接海温场
        print('>>>> prepare SST data <<<<')
        sp.run('ln -sf ungrib/Variable_Tables/Vtable.SST Vtable', shell=True)
        sp.run('./link_grib.csh ' + os.path.join(data_dir, sst_file), shell=True)
        # 修改namelist.wps
        nml_wps = modify_wps_nml()
        # 复写 namelist.wps
        rewrite_namelist(wps_dir, wps_nml_name, nml_wps, 1)
        # ungrid 解码
        print('>>>> ungrib.exe <<<<')
        sp.run('./ungrib.exe >& ungrib.log', shell=True)
    else:
        return 0

def alter_lswt(alternative_lswt):
    '''
    用WRF默认程式替换初始场湖温（WRFUsersGuide 3-28）
    '''
    if alternative_lswt == 1:
        print('>>>> cal avg_tsfc <<<<')
        sp.run('./util/avg_tsfc.exe > avg_tsfc.out', shell=True)
    else:
        return 0

def remove_lake(alternative_lake):
    '''
    用自己编写的脚本将湖泊（纳木错）替换为相邻土地利用类型
    '''
    if alternative_lake == 1:
        print('>>>> remove lake <<<<')
        # sp.run('conda activate geocat', shell=True)
        sp.run(f'python {comp_run_dir}/remove_lake.py', shell=True)
    else:
        return 0

def modify_lakedepth(md_lakedepth, alternative_lake):
    '''
    用观测湖盆深度替换默认的NamCo湖深
    '''
    if alternative_lake == 1:
        print('*** Warning: no need for modifing lakedepth ***')
    elif md_lakedepth == 1:
        print('>>>> modify lakedepth <<<<')
        # sp.run('conda activate geocat', shell=True)
        sp.run(f'python {comp_run_dir}/lakedepth.py', shell=True)
    else:
        return 0

def modify_mountainHeight(md_mountainHeight):
    '''
    去除念青唐古拉山
    '''
    if md_mountainHeight == 1:
        print('>>>> modify mountainHeight <<<<')
        # sp.run('conda activate geocat', shell=True)
        sp.run(f'python {comp_run_dir}/remove_mountain.py', shell=True)
    else:
        return 0

def modify_lswt(lswt_init_flag):
    '''
    去除念青唐古拉山
    '''
    if lswt_init_flag == 1:
        print('>>>> modify LSWT <<<<')
        sp.run(f'python {comp_run_dir}/lswt.py', shell=True)
    else:
        return 0

def copy_iofiles():
    '''
    拷贝自定义输入输出文件到WRF/run目录下
    '''
    source_path = os.path.join(root_dir, 'comp_run_new', f'{iofile_name}_d0x.txt')
    for i in range(max_dom):
        target_path = os.path.join(wrf_dir, 'run', f'{iofile_name}_d0{i+1}.txt')
        sp.run(f'cp {source_path} {target_path}', shell=True)
        print(f'>>>> copy {iofile_name}_d0{i+1}.txt complete <<<<')

if __name__ == '__main__':
    # 打印基本信息
    describe()

    # 运行WPS
    run_wps()

    # 运行WRF
    run_wrf()
