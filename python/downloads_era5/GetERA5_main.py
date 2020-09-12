import subprocess as sp 
import os 
from params import *

if __name__ == '__main__':
    if os.path.exists(Target_dir) == False: os.mkdir(target_dir)
    
    print(f'*** Downloading {DATE1}-{DATE2} ERA5 pressure levels data ***')
    sp.run('python GetERA5_pl.py', shell=True)
    print(f'*** Downloading {DATE1}-{DATE2} ERA5 surface data ***')
    sp.run('python GetERA5_sl.py', shell=True)
    

