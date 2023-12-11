import cdsapi
from datetime import datetime, timedelta
import os
import sys
from pathlib import Path
import yaml
import sys

c = cdsapi.Client()

years = [2009]


months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
#months = ['09']


days_30 = ['09', '04', '06', '11']
days_28 = ['02']


for year in years: 
    os.chdir('/disk1/Data/CAMS/aerosol_levels/{}'.format(year))
    for m in months: 

        dest_path = Path(m)
        if not dest_path.exists():
            dest_path.mkdir()
            os.chdir(dest_path)
            print('Making new folder:', dest_path)
        else:
            print('This folder already exists')
            continue

        if m in days_30:
            days = 30

        elif m in days_28:
            days = 28
        else:
            days = 31

        starting = datetime(year,int(m),1)    
        starting = starting.strftime('%Y-%m-%d')
        end = datetime(year,int(m),days) 
        end = end.strftime('%Y-%m-%d')
        
        filename = starting + '-' + end + '.nc'
        #Account for end of month...
        #days = (end - starting).days + 1       
        #for one_day in range(0,days):      
#             this_date = (starting + timedelta(days = one_day)).strftime('%Y-%m-%d')
#             filename = "%s.nc"%this_date  
#             day = one_day + 1
#             d = str(day).zfill(2)
            
        c.retrieve(
            'cams-global-reanalysis-eac4',
            {
                'date': f'{starting}/{end}',
                'format': 'netcdf',
                'variable': [
                    'dust_aerosol_0.03-0.55um_mixing_ratio', 'dust_aerosol_0.9-20um_mixing_ratio', 'hydrophilic_black_carbon_aerosol_mixing_ratio',
                    'hydrophilic_organic_matter_aerosol_mixing_ratio', 'hydrophobic_black_carbon_aerosol_mixing_ratio', 'hydrophobic_organic_matter_aerosol_mixing_ratio',
                    'sea_salt_aerosol_0.03-0.5um_mixing_ratio', 'sea_salt_aerosol_0.5-5um_mixing_ratio', 'sea_salt_aerosol_5-20um_mixing_ratio',
                    'sulphate_aerosol_mixing_ratio', 'sulphur_dioxide'
                ],
                'pressure_level': [
                    '100', '200', '300',
                    '400', '500', '600',
                    '700', '800', '900',
                    '1000',
                ],
                'time': [
                    '00:00', '03:00', '06:00',
                    '09:00', '12:00', '15:00',
                    '18:00', '21:00',
                ],
                'area': [
                    90, -180, -90,
                    180,
                ],
            },
            filename)
        os.chdir('..')
