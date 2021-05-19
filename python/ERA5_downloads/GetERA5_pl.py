import cdsapi
from params import *

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type':'reanalysis',
        'format':'grib',
        'pressure_level':[
            '1','2','3',
            '5','7','10',
            '20','30','50',
            '70','100','125',
            '150','175','200',
            '225','250','300',
            '350','400','450',
            '500','550','600',
            '650','700','750',
            '775','800','825',
            '850','875','900',
            '925','950','975',
            '1000'
        ],
        'date':f'{DATE1}/{DATE2}',
        'area':f'{Nort}/{West}/{Sout}/{East}',
        'time':f'00/to/23/by/{time_interval}',
        'variable':[
            'geopotential','relative_humidity','specific_humidity',
            'temperature','u_component_of_wind','v_component_of_wind'
        ]
    },
    f'{Target_dir}/ERA5-{DATE1}-{DATE2}-pl.grib')
