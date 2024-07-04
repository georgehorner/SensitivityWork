import netCDF4 as nc
import xarray as xr
import numpy as np
import logging
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import os
import math
import scipy.stats as stats
#import imageio
from scipy.ndimage.interpolation import map_coordinates
from pprint import pprint
from glob import glob
from scipy import interpolate
import scipy.signal as sp
from scipy.interpolate import RegularGridInterpolator

import logging

# Create a logging instance
#logger = logging.getLogger('TSC_long')
#logger.setLevel(logging.INFO) # you can set this to be DEBUG, INFO, ERROR
# Assign a file-handler to that instance
#fh = logging.FileHandler("file_dir.txt")
#fh.setLevel(logging.INFO) # again, you can set this differently
# Format your logs (optional)
#fh.setFormatter(formatter) # This will set the format to the file handler
##formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# Add the handler to your logging instance
#logger.addHandler(fh)
#print = logger.info # this will make 'print' call the logger.info function




#####################################################
print('running')
def round_three(n):
    '''round nearest integer to multiple of three'''
    rounded = np.int(3*np.round(n/3.0))
    if n<2:
        return rounded
    else:
        return int(rounded / 3)

########################################################
def latlongrid(initlat,finlat,initlon,finlon,step):
    '''Create empty latlon grid'''
    latlen = finlat - initlat
    lonlen = finlon - initlon 
    
    tlats = np.arange(initlat,finlat+1,step)
    tlons = np.arange(initlon,finlon+1,step)
    
    #initialise empty lat and lon array values for tracers
    ilats = np.zeros((len(tlats),len(tlons)))
    ilons = np.zeros((len(tlats),len(tlons)))

    #populate initial tracer positions (-29.5,29.5 lat / 0.5 360.5 lon) from lat/lon arrays
    for i in range(len(tlons)):
        ilats[:,i] = tlats
    for i in range(len(tlats)):
        ilons[i,:] = tlons
    return np.array([ilats,ilons])
########################################################################
def nanmean(data, *args, **kwargs):
    '''Returns of mean of a numpy array ignoring missing data'''
    return np.ma.filled(np.ma.array(
        data, mask=np.bitwise_not(np.isfinite(data))).mean(*args, **kwargs),
                        np.nan)

########################################################################



def find_closest(A, target):
    '''Find closest advected pixel to gridbox'''
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx

def round_to_index(lati,loni,value):
    '''rounds float to index of array'''
    latlon = np.empty((1801,3600))
    latlon[:] = np.nan
    
    index_lat_f = 10*x_round(lati)
    index_lon_f = 10*x_round(loni)
        
    #index_lat_f = index_lat[~np.isnan(index_lat)]
    #index_lon_f = index_lon[~np.isnan(index_lon)]
    
    
    index_lat_f = index_lat_f + 900
    index_lat_f = np.where(index_lat_f < 1800,index_lat_f,np.nan)
    index_lon_f = np.where((index_lon_f>0) & (index_lon_f<3600),index_lon_f,index_lon_f % 3600)

    latlon[index_lat_f.astype(int),index_lon_f.astype(int)] = value
    return latlon

def resizetsc(inputtsc):
    '''resizes array from 1 to 0.1 degree'''
    s1 = reduceres(inputtsc)
    s2 = regrid(s1,641,3600)
    return s2

def x_round(x):
    '''rounds float to int'''
    return np.round(x*10)/10

def reduceres(inputdata):
    '''reduces array resolution from 0.1 to 1 degree'''
    test1 = reduce_res_axis(inputdata,5,1)
    test2 = reduce_res_axis(test1[0:640],5,0)
    return test2

def reduce_res_axis(data, res, axis=0, func=nanmean):
    '''Reduce the resolution of an array along a specified axis (axis) by
    a factor (res), using a function (func).'''
    if not isinstance(res, int):
        print('res converted to int')
        res = int(res)
    dshape = data.shape
    if axis < 0:
        axis = len(dshape) + axis
    dshape = np.concatenate((dshape[:axis], np.array(
        [dshape[axis]//res, res]), dshape[(axis+1):]))
    return func(data.reshape(dshape.astype(int)), axis=axis+1)

def resize(array):
    '''resizes DCC cores from 1deg to 0.1deg'''
    #270,1020
    new_dims = []
    for original_length, new_length in zip(array.shape, (600,3600)):
        new_dims.append(np.linspace(0, original_length-1, new_length))

    coords = np.meshgrid(*new_dims, indexing='ij')
    dcc_scale = map_coordinates(array, coords)
    return dcc_scale

def regrid(data, out_x, out_y):
    m = max(data.shape[0], data.shape[1])
    y = np.linspace(0, 1.0/m, data.shape[0])
    x = np.linspace(0, 1.0/m, data.shape[1])
    interpolating_function = RegularGridInterpolator((y, x), data)

    yv, xv = np.meshgrid(np.linspace(0, 1.0/m, out_y), np.linspace(0, 1.0/m, out_x))

    return interpolating_function((xv, yv))

def reduceres1deg(inputdata):
    test1 = reduce_res_axis(inputdata,10,1)
    test2 = reduce_res_axis(test1[0:640],10,0)
    return test2


def reduceres1deg_mode(inputdata):
    test1 = reduce_res_axis_MODE(inputdata,10,1)
    test2 = reduce_res_axis_MODE(test1[0:640],10,0)
    return test2

############################################################

def tracer_grid(inittime, time_u,latgrid,longrid):
    '''Main tracer function. Takes grid of lat/lon points and advects them in windfield'''
    
    time = np.arange(time_u)
    
    #initialise empty lat and lon array values for tracers
    ilats = np.zeros((len(time),len(latgrid[:,0]),len(longrid[0,:])))
    ilons = np.zeros((len(time),len(latgrid[:,0]),len(longrid[0,:])))

    #populate initial tracer positions (-29.5,29.5 lat / 0.5 360.5 lon) from lat/lon arrays
    
    ilats[0] = latgrid
    ilons[0] = longrid
                     
        #np.fromfunction <- arrays from a function
        
        #loop through each time step
    for i in range(len(time)-1):
    
        #conversion from km to 1deg of lat/lon
        lat1deg = 111
        lon1deg = np.zeros((len(time),60,360))
        #lon conversion depends on latitude
        lon1deg = 111*np.cos(np.radians(ilats))
    
        #print('latitudes are: ' +str(ilats[i]))
        #print('longitudes are: ' +str(ilons[i]))
            
        #if tracer longitude has been advected beyond 360, reset to -360 (bugs out otherwise)
        #modulus np.mod function
            
        #find the closest value between advected lat/lon to the original lat/lon arrays and 
        #return the index values of the lat/lon array
        inlats = find_closest(lat,-ilats[i])
        inlons = find_closest(lon,ilons[i])
    
        #print('index of lat values: ' +str(inlats))
        #print('index of lon values: ' +str(inlons))
        
        uwindmean = (u300[inittime+i]+ u200[inittime+i]) / 2
        vwindmean = (v300[inittime+i]+ v200[inittime+i]) / 2
        
        #npuwind = np.array(uwind[inittime+i])
        #npvwind = np.array(vwind[inittime+i])
        
        npuwind = np.array(uwindmean)
        npvwind = np.array(vwindmean)
            
        #define advection using these new defined inlat/inlot values
        #d_u = 3.6 * uwind[inittime+i,inlats,inlons]
        #d_v = 3.6 * vwind[inittime+i,inlats,inlons]
        
        d_u = 3.6 * npuwind[inlats,inlons]
        d_v = 3.6 * npvwind[inlats,inlons]
        
        dlat =  (d_v/lat1deg)
        dlon = (d_u/lon1deg[i])
        
        #advect lat/lons to to new timestep
        
        #print('ilats MIN before advection: ' +str(np.nanmin(ilats[i])))
        #print('ilats MAX before advection: ' +str(np.nanmax(ilats[i])))
    
        ilats[i+1] = ilats[i] + dlat
        ilons[i+1] = ilons[i] + dlon
        
        #print('ilats MIN after advection: ' +str(np.nanmin(ilats[i+1])))
        #print('ilats MAX after advection: ' +str(np.nanmax(ilats[i+1])))

        ilons = np.where((ilons>0) & (ilons<360),ilons,ilons % 360)
        ilats = np.where(ilats>-32,ilats,-32)
        ilats = np.where(ilats<32,ilats,32)
        #print('ilats MIN after advection and rounding: ' +str(np.nanmin(ilats[i+1])))
        #print('ilats MAX after advection and rounding: ' +str(np.nanmax(ilats[i+1])))
        #if i % 50 == 0:
            #print(str(i) +' time stamps done!')
            
    #return the array of each tracer position at each time step
    #print('done!')
    return np.array([ilats,ilons])

def round_to_index_REDUCE(lati,loni,value):
    '''index rounding function that deals with two pixels advecting into same gridbox'''
    latlon = np.empty((641,3600))
    latlon[:] = np.nan
    
    index_lat_f = 10*x_round(lati)
    index_lon_f = 10*x_round(loni)
        
    #index_lat_f = index_lat[~np.isnan(index_lat)]
    #index_lon_f = index_lon[~np.isnan(index_lon)]
    
    index_lat_f = index_lat_f + 320
    index_lat_f = np.where(index_lat_f > 640,640,index_lat_f)
    index_lat_f = np.where(index_lat_f > 0,index_lat_f,0)

    index_lon_f = np.where((index_lon_f>0) & (index_lon_f<3600),index_lon_f,index_lon_f % 3600)
    
    latlon[index_lat_f.astype(int),index_lon_f.astype(int)] = value
    latlon[638:640,:] = np.nan
    latlon[0:2,:] = np.nan
    latlonlowest = np.where(latlon[index_lat_f.astype(int),index_lon_f.astype(int)] > value,value,latlon)
    
    return latlonlowest

def convolve(input):
    '''interpolation to fill in missing data from divergent pixels'''
    interp1 = sp.convolve2d(np.where(~np.isfinite(input),0,input),np.ones((3,3)),boundary='fill')
    interp2 = sp.convolve2d(np.isfinite(input),np.ones((3,3)),boundary='fill')
    final = (interp1/interp2)[1:-1,1:-1]
    return final

########### CODE STARTS HERE ##################

#load in dataset used for land/ocean flag
fn2 = nc.Dataset('/net/seldon/disk2/Data/ECMWF/ERA5/SST_surf_1grid_timed/2020/SST_surf_001.nc')
sst = fn2['sst'][:]

import time
print('starting TSC calcs')
start_time = time.time()

#initialise lat lon array
itrace_lat = latlongrid(-32,31.1,0.0,359.0,0.1)[0]
itrace_lon = latlongrid(-32,31.1,0.0,359.0,0.1)[1]

lats = np.arange(-32,32.1,0.1)
lons = np.arange(0,360,0.1)
onedeglats = np.arange(-31.5,32.5,1)
onedeglons = np.arange(0.5,360.5,1)

#create empty arrays used to hold flag values/tsc/advected values etc
oldconv = np.zeros((641,3600))
advectconv = np.zeros((641,3600))
#tsc = np.ones((641,3600))
#load new tsc
tsc = np.load('temp_tsc_2008_365.npy')
convorigin = np.zeros((641,3600))
#convaod = np.zeros((641,3600))
convpc = np.zeros((641,3600))

convcir = np.zeros((641,3600))
convlifetime = np.zeros((641,3600))

newconv = np.zeros((641,3600))
newconv[:] = np.nan

oceanconv = np.zeros((641,3600))
oceanconv[:] = np.nan
landconv = np.zeros((641,3600))
landconv[:] = np.nan

h_aodconv = np.zeros((641,3600))
h_aodconv[:] = np.nan
l_aodconv = np.zeros((641,3600))
l_aodconv[:] = np.nan

pclevels = np.zeros((4,641,3600))
pclevels[:] = np.nan

convlocorigin = np.zeros((2,641,3600))
convlocorigin[:] = np.nan


###land/ocean mask
maskl = np.where(sst[0]<0,0,1)
mask = np.kron(maskl,np.ones((10,10)))[600:1200]

#initialise start loop
i=0
j=0
k=0
hour=0
dcchour=0
isccphour=0

year = 2009


#Import ISCCP data
#isccp = xr.open_mfdataset('/net/seldon/disk2/Users/gah20/ISCCP/access/isccp-basic/hgg/200[7-9]*/200[7-9]*.nc',combine='nested',concat_dim='time')
isccp = xr.open_mfdataset('/disk1/Data/ISCCP/access/isccp-basic/hgg/'+str(year)+'*/'+str(year)+'*.nc',combine='nested',concat_dim='time')
taupc = isccp['n_pctaudist']
levtau = isccp['levtau']
levpc = isccp['levpc']
tauraw = isccp['tau']
pcraw = isccp['pc']
isccptime = isccp['time']

## Import the DCC data
fn_cf = '/disk1/Users/gah20/DCC/'+str(year)+'/*.nc'
ds_cf = xr.open_mfdataset(fn_cf,combine='nested',concat_dim='time')
ctp_dcc = ds_cf['DCC'][:,60:120]
lentime = ds_cf['Time']


fu2 = nc.Dataset('/net/hardin/disk1/Data/ECMWF/ERA5/U-wind-component_200hPa_0.25grid_timed/'+str(year)+'/U-wind-component_200hPa_'+str(j+1).zfill(3)+'.nc')
fv2 = nc.Dataset('/net/hardin/disk1/Data/ECMWF/ERA5/V-wind-component_200hPa_0.25grid_timed/'+str(year)+'/V-wind-component_200hPa_'+str(j+1).zfill(3)+'.nc')
u200 = fu2['u']
v200 = fv2['v']
lat = fu2['latitude'][:]*(-1)
lon = fu2['longitude'][:]

fu3 = nc.Dataset('/net/hardin/disk1/Data/ECMWF/ERA5/U-wind-component_300hPa_0.25grid_timed/'+str(year)+'/U-wind-component_300hPa_'+str(j+1).zfill(3)+'.nc')
fv3 = nc.Dataset('/net/hardin/disk1/Data/ECMWF/ERA5/V-wind-component_300hPa_0.25grid_timed/'+str(year)+'/V-wind-component_300hPa_'+str(j+1).zfill(3)+'.nc')
u300 = fu3['u']
v300 = fv3['v']
timerec = fv3['time']

#################################################################
#ft = nc.Dataset('TSC/'+str(year)+'/TEMPDATASET.nc','w',format='NETCDF4')
#Create netcdf file for TSC output
ft = nc.Dataset('TSC_2origin/'+str(year)+'/TSC_000.nc','w',format='NETCDF4')
ft.createDimension('lat', 64)
ft.createDimension('lon', 360)
ft.createDimension('time', None)
ft.createDimension('latlon',2)


latitude = ft.createVariable('Latitude', 'f4', 'lat')
longitude = ft.createVariable('Longitude', 'f4', 'lon')
ntsc = ft.createVariable('TSC', 'f4', ('time', 'lat', 'lon',))
nconvorigin = ft.createVariable('ConvOrigin', 'f4', ('time', 'lat', 'lon',))
nconvcir = ft.createVariable('ConvCir','f4',('time','lat','lon'))
nconvpc = ft.createVariable('ConvPc','f4',('time','lat','lon'))
nconvlocorigin = ft.createVariable('LocOrigin','f4',('time','latlon','lat','lon'))


timen = ft.createVariable('Time', 'i4', 'time')
lifetimen = ft.createVariable('Lifetime', 'f4', ('time','lat','lon',))
longitude[:] = onedeglons
latitude[:] = onedeglats
##################################################################

#9504
#17520



for i in range(0,1000000):
    if (hour>0) & (i%24 == 0):
        k=0
        j=j+1
        ft.close()
        print('File ' +str(j-1)+ ' completed!')
        
        if (year % 4 == 0) or (year == 2000): #deal with leap year
            if (j>0) & (j%366 == 0):
                print('Year ' +str(year)+ ' completed!') 
                year+=1
                j=0

                fn_cf = '/disk1/Users/gah20/DCC/'+str(year)+'/*.nc'
                ds_cf = xr.open_mfdataset(fn_cf,combine='nested',concat_dim='time')
                ctp_dcc = ds_cf['DCC'][:,60:120]
                lentime = ds_cf['Time']

                isccp = xr.open_mfdataset('/disk1/Data/ISCCP/access/isccp-basic/hgg/'+str(year)+'*/'+str(year)+'*.nc',combine='nested',concat_dim='time')
                taupc = isccp['n_pctaudist']
                levtau = isccp['levtau']
                levpc = isccp['levpc']
                tauraw = isccp['tau']
                pcraw = isccp['pc']
                isccptime = isccp['time']

                hour = 0
                dcchour=0
                isccphour=0

        else:
            if (j>0) & (j%365 == 0):
                year+=1
                j=0

                fn_cf = '/disk1/Users/gah20/DCC/'+str(year)+'/*.nc'
                ds_cf = xr.open_mfdataset(fn_cf,combine='nested',concat_dim='time')
                ctp_dcc = ds_cf['DCC'][:,60:120]
                lentime = ds_cf['Time']

                isccp = xr.open_mfdataset('/disk1/Data/ISCCP/access/isccp-basic/hgg/'+str(year)+'*/'+str(year)+'*.nc',combine='nested',concat_dim='time')
                taupc = isccp['n_pctaudist']
                levtau = isccp['levtau']
                levpc = isccp['levpc']
                tauraw = isccp['tau']
                pcraw = isccp['pc']
                isccptime = isccp['time']
                hour=0
                isccphour = 0 

                hour = 0
                dcchour = 0 
                isccphour = 0 
                #print('Year ' +str(year)+ ' completed!')

                
        #import WINDdata for each timestep

        
        fu2 = nc.Dataset('/net/hardin/disk1/Data/ECMWF/ERA5/U-wind-component_200hPa_0.25grid_timed/'+str(year)+'/U-wind-component_200hPa_'+str(j+1).zfill(3)+'.nc')
        fv2 = nc.Dataset('/net/hardin/disk1/Data/ECMWF/ERA5/V-wind-component_200hPa_0.25grid_timed/'+str(year)+'/V-wind-component_200hPa_'+str(j+1).zfill(3)+'.nc')
        u200 = fu2['u']
        v200 = fv2['v']
        
        fu3 = nc.Dataset('/net/hardin/disk1/Data/ECMWF/ERA5/U-wind-component_300hPa_0.25grid_timed/'+str(year)+'/U-wind-component_300hPa_'+str(j+1).zfill(3)+'.nc')
        fv3 = nc.Dataset('/net/hardin/disk1/Data/ECMWF/ERA5/V-wind-component_300hPa_0.25grid_timed/'+str(year)+'/V-wind-component_300hPa_'+str(j+1).zfill(3)+'.nc')
        u300 = fu3['u']
        v300 = fv3['v']
        
        
        ft = nc.Dataset('TSC_2origin/'+str(year)+'/TSC_'+(str(j).zfill(3))+'.nc','w',format='NETCDF4')

        ft.createDimension('lat', 64)
        ft.createDimension('lon', 360)
        ft.createDimension('time', None)
        ft.createDimension('latlon',2)

        latitude = ft.createVariable('Latitude', 'f4', 'lat')
        longitude = ft.createVariable('Longitude', 'f4', 'lon')
        ntsc = ft.createVariable('TSC', 'f4', ('time', 'lat', 'lon',))
        nconvorigin = ft.createVariable('ConvOrigin','f4',('time','lat','lon',))
        nconvpc = ft.createVariable('ConvPc','f4',('time','lat','lon'))
        nconvcir = ft.createVariable('ConvCir','f4',('time','lat','lon'))
        nconvlocorigin = ft.createVariable('LocOrigin','f4',('time','latlon','lat','lon'))


        timen = ft.createVariable('Time', 'i4', 'time')
        lifetimen = ft.createVariable('Lifetime','f4',('time','lat','lon',))
        longitude[:] = onedeglons
        latitude[:] = onedeglats
    
    if (round_three(isccphour) >= len(pcraw)):
        print('hour is: ' +str(round_three(hour)))
        print('len is: ' +str(len(pcraw)))

        isccp = xr.open_mfdataset('/disk1/Data/ISCCP/access/isccp-basic/hgg/'+str(year+1)+'*/'+str(year+1)+'*.nc',combine='nested',concat_dim='time')
        taupc = isccp['n_pctaudist']
        levtau = isccp['levtau']
        levpc = isccp['levpc']
        tauraw = isccp['tau']
        pcraw = isccp['pc']
        isccptime = isccp['time']
        hour=0
        isccphour = 0 

    if round_three(dcchour) >= len(ctp_dcc):
        fn_cf = '/disk1/Users/gah20/DCC/'+str(year+1)+'/*.nc'
        ds_cf = xr.open_mfdataset(fn_cf,combine='nested',concat_dim='time')
        ctp_dcc = ds_cf['DCC'][:,60:120]
        lentime = ds_cf['Time']
        hour=0
        dcchour = 0

    #create DCC flag
    dcc_norm = (ctp_dcc[round_three(dcchour)]/ctp_dcc[round_three(dcchour)])


    dcc_scale =np.kron(dcc_norm,np.ones((10,10)))
    newconv[20:620] = dcc_scale/dcc_scale
    #aod4_scale = np.kron(aodgr[6+i],np.ones((10,10)))

    
    oceanconv[20:620] = np.where(mask==1,newconv[20:620],np.nan)
    landconv[20:620] = np.where(mask==0,newconv[20:620],np.nan)

    
    #h_aodconv[20:620] = np.where(aod4_scale>0.3,newconv[20:620],np.nan)
    #l_aodconv[20:620] = np.where(aod4_scale<0.3,newconv[20:620],np.nan)
    #print('hour is: ' +str(hour))

    #deals with rounding to next dataset - temporarily sets hour back to zero for next year

    pc = np.kron(pcraw[round_three(isccphour),60:120],np.ones((10,10)))


    pclevels[0,20:620] = np.where((pc<200)&(pc>175),newconv[20:620],np.nan)
    pclevels[1,20:620] = np.where((pc<175)&(pc>150),newconv[20:620],np.nan)
    pclevels[2,20:620] = np.where((pc<150)&(pc>125),newconv[20:620],np.nan)
    pclevels[3,20:620] = np.where((pc<125)&(pc>100),newconv[20:620],np.nan)

    #create insitu v real flag
    cf = np.nansum(np.nansum(taupc[round_three(isccphour),:,:3],axis=0),axis=0)
    cf = np.where(cf<0,np.nan,cf)
    cirrus = np.where(cf[60:120]>10,1,0)
    cirrus = np.where(np.isnan(cf[60:120]),np.nan,cirrus)
    cirrus_scale = np.kron(cirrus,np.ones((10,10)))
    
    #advect
    trace = tracer_grid(k,2,itrace_lat, itrace_lon)


    tsc[~np.isnan(newconv)] = 0
    convorigin[~np.isnan(oceanconv)] = -100
    convorigin[~np.isnan(landconv)] = 100
    #[~np.isnan(h_aodconv)] = 100
    #convaod[~np.isnan(l_aodconv)] = np.nan
    convpc[~np.isnan(pclevels[0])] = 1500
    convpc[~np.isnan(pclevels[1])] = 500
    convpc[~np.isnan(pclevels[2])] = -500
    convpc[~np.isnan(pclevels[3])] = -1500
    convlocorigin[:,~np.isnan(newconv)] = np.where(~np.isnan(newconv))

    convcir[~np.isnan(newconv)] = 100
    
    advectconv = round_to_index_REDUCE(trace[0,1],trace[1,1],tsc)
    advectconvorigin = round_to_index_REDUCE(trace[0,1],trace[1,1],convorigin)
    #advectconvaod = round_to_index_REDUCE(trace[0,1],trace[1,1],convaod)
    advectconvpc = round_to_index_REDUCE(trace[0,1],trace[1,1],convpc)
    advectconvcir = round_to_index_REDUCE(trace[0,1],trace[1,1],convcir)
    advectlocoriginlat = round_to_index_REDUCE(trace[0,1],trace[1,1],convlocorigin[0])
    advectlocoriginlon = round_to_index_REDUCE(trace[0,1],trace[1,1],convlocorigin[1])


    #if currently detrained but now CF<0.1, switch to detrained
    convlifetime[20:620] = np.where((cirrus_scale==0)&(advectconvcir[20:620]>0),advectconv[20:620],np.nan)
    advectconvcir[20:620] = np.where((cirrus_scale==0)&(advectconvcir[20:620]>0),-100,advectconvcir[20:620])

    #print(convlifetime)
    ## at this point, note down the TSC of the locations where this happens 
    ### SAVE EACH TSC TO A NETCDF FILE HERE;
    
    timen[k] = timerec[k]
    tsc1deg = reduceres1deg(tsc)
    lifetimen[k] = reduceres1deg(convlifetime)
    ntsc[k] = tsc1deg
    
    convorigin1deg = reduceres1deg(convorigin)
    convpc1deg = reduceres1deg(convpc)
    convcir1deg = reduceres1deg(convcir)
    convlocorigin1deglat = reduceres1deg(convlocorigin[0])
    convlocorigin1deglon = reduceres1deg(convlocorigin[1])

    nconvorigin[k] = convorigin1deg
    nconvpc[k] = convpc1deg
    nconvcir[k] = convcir1deg
    nconvlocorigin[k,0] = convlocorigin1deglat
    nconvlocorigin[k,1] = convlocorigin1deglon

    k=k+1 ##reset array position for nc file to 0 when new file is created

    tsc = advectconv + 1
    convorigin = advectconvorigin
    convpc = advectconvpc
    convcir = advectconvcir
    convlocorigin[0] = advectlocoriginlat
    convlocorigin[1] = advectlocoriginlon

    tsc_filler = convolve(tsc)
    conv_filler = convolve(convorigin)
    convloc_fillerlat = convolve(convlocorigin[0])
    convloc_fillerlon = convolve(convlocorigin[1])

    convpc_filler = convolve(convpc)
    convcir_filler = convolve(convcir)
    
    tsc = np.where(~np.isfinite(tsc),tsc_filler,tsc)
    convorigin = np.where(~np.isfinite(convorigin),conv_filler,convorigin)
    #convaod = np.where(~np.isfinite(convaod),convaod_filler,convaod)
    convpc = np.where(~np.isfinite(convpc),convpc_filler,convpc)
    convlocorigin[0] = np.where(~np.isfinite(convlocorigin[0]),convloc_fillerlat,convlocorigin[0])
    convlocorigin[1] = np.where(~np.isfinite(convlocorigin[1]),convloc_fillerlon,convlocorigin[1])

    convcir = np.where(~np.isfinite(convcir),convcir_filler,convcir)

    
    if (k==23):
        if os.path.exists('temp_tsc_'+str(year)+'_'+str(j-1).zfill(3)+'.npy'):
            #print('PATH EXISTS')
            os.remove('temp_tsc_'+str(year)+'_'+str(j-1).zfill(3)+'.npy')
        np.save('temp_tsc_'+str(year)+'_'+str(j).zfill(3)+'.npy',tsc)
    #print('hour is: ' +str(hour))

    hour = hour+1
    dcchour = dcchour + 1
    isccphour = isccphour + 1

ft.close()
print('Time elapsed: ' + str(time.time() - start_time)/3600)
