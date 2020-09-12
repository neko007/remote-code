import subprocess as sp
import pandas as pd 

def curl(target_dir, year, month, day, hour, count):
    sp.run(f'curl -s --disable-epsv --connect-timeout 30 -m 60 -u anonymous:USER_ID@INSTITUTION -o {target_dir}/GFS_{count} ftp://ftpprd.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.{year}{month:0>2d}{day:0>2d}/00/gfs.t00z.pgrb2.0p50.f0{hour:0>2d}', shell=True)

if __name__ == '__main__':
    target_dir = '../data'
    start_date = '2020-09-01 00:00:00'
    end_date = '2020-09-03 00:00:00'
    date_range = pd.date_range(start_date, end_date, freq='3H')
    for i, date in enumerate(date_range):
        print(f'>>> Download {date} Data <<<')
        year, month, day, hour = [date.year, date.month, date.day, date.hour]
        curl(target_dir, year, month, day, hour, i)


