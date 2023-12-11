import netCDF4 as nc
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import os
import math
#import imageio
from scipy.ndimage.interpolation import map_coordinates
from pprint import pprint
from glob import glob
from scipy import interpolate
import scipy.signal as sp
from scipy.interpolate import RegularGridInterpolator

print('running')
def round_three(n):
    '''round nearest integer to multiple of three'''
    rounded = np.int(3*np.round(n/3.0))
    if n<2:
        return rounded
    else:
        return int(rounded / 3)


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

def nanmean(data, *args, **kwargs):
    '''Returns of mean of a numpy array ignoring missing data'''
    return np.ma.filled(np.ma.array(
        data, mask=np.bitwise_not(np.isfinite(data))).mean(*args, **kwargs),
                        np.nan)

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
        
        ilats[i+1] = ilats[i] + dlat
        ilons[i+1] = ilons[i] + dlon
        
        ilons = np.where((ilons>0) & (ilons<360),ilons,ilons % 360)
            
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
    index_lat_f = np.where(index_lat_f < 640,index_lat_f,640)
    index_lon_f = np.where((index_lon_f>0) & (index_lon_f<3600),index_lon_f,index_lon_f % 3600)
    
    latlon[index_lat_f.astype(int),index_lon_f.astype(int)] = value
    latlon[640,:] = np.nan
    latlon[0,:] = np.nan
    latlonlowest = np.where(latlon[index_lat_f.astype(int),index_lon_f.astype(int)] > value,value,latlon)
    
    return latlonlowest

def convolve(input):
    '''interpolation to fill in missing data from divergent pixels'''
    interp1 = sp.convolve2d(np.where(~np.isfinite(input),0,input),np.ones((3,3)),boundary='wrap')
    interp2 = sp.convolve2d(np.isfinite(input),np.ones((3,3)),boundary='wrap')
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
tsc = np.zeros((641,3600))
convorigin = np.zeros((641,3600))
convaod = np.zeros((641,3600))
convpc = np.zeros((641,3600))
convcir = np.zeros((641,3600))
convlifetime = np.zeros((641,3600))
#convsa = np.zeros((641,3600))
#convaf = np.zeros((641,3600))
#convmc = np.zeros((641,3600))

#tsc = last day of 2007 data
#tsc = nc.Dataset('/net/seldon/disk2/Users/gah20/TSC_Young2/2007/TSC_1094.nc')['TSC'][23]
    #tsc[:] = np.nan
newconv = np.zeros((641,3600))
newconv[:] = np.nan

oceanconv = np.zeros((641,3600))
oceanconv[:] = np.nan
landconv = np.zeros((641,3600))
landconv[:] = np.nan


#conv_sa = np.zeros((641,3600))
#conv_nsa = np.zeros((641,3600))
#conv_sa[:] = np.nan
#conv_nsa[:] = np.nan
#conv_af = np.zeros((641,3600))
#conv_af[:] = np.nan
#conv_naf = np.zeros((641,3600))
#conv_naf[:] = np.nan
#conv_mc = np.zeros((641,3600))
#conv_mc[:] = np.nan
#conv_nmc = np.zeros((641,3600))
#conv_nmc[:] = np.nan

#conv_af_int = np.zeros((641,3600))
#conv_af_int[:] = np.nan
#conv_mc_int = np.zeros((641,3600))
#conv_mc_int[:] = np.nan

h_aodconv = np.zeros((641,3600))
h_aodconv[:] = np.nan
l_aodconv = np.zeros((641,3600))
l_aodconv[:] = np.nan

#regionalflag = np.zeros((641,3600))
#regionalflag[:] = np.nan
#regionalflag[100:400,-850:-450] = 1
#regionalflag[100:450,100:400] = 2
#regionalflag[200:500,700:1550] = 3

highpc = np.zeros((641,3600))
highpc[:] = np.nan
lowpc = np.zeros((641,3600))
lowpc[:] = np.nan

###land/ocean mask
maskl = np.where(sst[0]<0,0,1)
mask = np.kron(maskl,np.ones((10,10)))[600:1200]

## Import the DCC data
fn_cf = '/net/seldon/disk2/Users/gah20/DCC/ISCCP/Total7_10.nc'
ds_cf = nc.Dataset(fn_cf)
ctp_dcc = ds_cf['DCC'][:,60:120]
lentime = ds_cf['Time']

#initialise start loop
i=0
j=0
k=0
year = 2007


## Import AOD data
fn_aod = xr.open_mfdataset('/net/seldon/disk2/Users/gah20/MODIS/subset/*/aod_int_8hr*.nc',combine='nested',concat_dim='time')['AOD_Terra']
aodgr = np.roll(np.flip(fn_aod[:],axis=1),180,axis=2)
aodmean = np.nanmean(aodgr,axis=0)
aodstd = np.nanstd(aodgr,axis=0)
aod2sigp = np.where(aodgr>(aodmean + 2*aodstd),aodgr,np.nan)
aod2sign = np.where(aodgr<(aodmean - 2*aodstd),aodgr,np.nan)

#Import ISCCP data
isccp = xr.open_mfdataset('/net/seldon/disk2/Users/gah20/ISCCP/access/isccp-basic/hgg/200[7-9]*/200[7-9]*.nc',combine='nested',concat_dim='time')
taupc = isccp['n_pctaudist']
levtau = isccp['levtau']
levpc = isccp['levpc']
tauraw = isccp['tau']
pcraw = isccp['pc']
isccptime = isccp['time']

## Import the windfield lat/londata
fn_wuLAT = '/net/seldon/disk1/Users/gah20/data/ECMWF/era5_u_component_of_wind_'+str(year)+'_hourly.nc'
ds_wL = nc.Dataset(fn_wuLAT)
lat = ds_wL['latitude'][:]*(-1)
lon = ds_wL['longitude'][:]

fu2 = nc.Dataset('/net/seldon/disk2/Data/ECMWF/ERA5/U-wind-component_200hPa_0.25grid_timed/'+str(year)+'/U-wind-component_200hPa_'+str(j+1).zfill(3)+'.nc')
fv2 = nc.Dataset('/net/seldon/disk2/Data/ECMWF/ERA5/V-wind-component_200hPa_0.25grid_timed/'+str(year)+'/V-wind-component_200hPa_'+str(j+1).zfill(3)+'.nc')
u200 = fu2['u']
v200 = fv2['v']

fu3 = nc.Dataset('/net/seldon/disk2/Data/ECMWF/ERA5/U-wind-component_300hPa_0.25grid_timed/'+str(year)+'/U-wind-component_300hPa_'+str(j+1).zfill(3)+'.nc')
fv3 = nc.Dataset('/net/seldon/disk2/Data/ECMWF/ERA5/V-wind-component_300hPa_0.25grid_timed/'+str(year)+'/V-wind-component_300hPa_'+str(j+1).zfill(3)+'.nc')
u300 = fu3['u']
v300 = fv3['v']
timerec = fv3['time']

#################################################################
#Create netcdf file for TSC output
ft = nc.Dataset('TSC/'+str(year)+'/TSC_000.nc','w',format='NETCDF4')

ft.createDimension('lat', 64)
ft.createDimension('lon', 360)
ft.createDimension('time', None)

latitude = ft.createVariable('Latitude', 'f4', 'lat')
longitude = ft.createVariable('Longitude', 'f4', 'lon')
ntsc = ft.createVariable('TSC', 'f4', ('time', 'lat', 'lon',))
nconvorigin = ft.createVariable('ConvOrigin', 'f4', ('time', 'lat', 'lon',))
nconvaod = ft.createVariable('ConvAOD','f4',('time','lat','lon',))
nconvpc = ft.createVariable('ConvPc','f4',('time','lat','lon'))
nconvcir = ft.createVariable('ConvCir','f4',('time','lat','lon'))
#nconvsa = ft.createVariable('Conv_SA','f4',('time','lat','lon'))
#nconvaf = ft.createVariable('Conv_AF','f4',('time','lat','lon'))
#nconvmc = ft.createVariable('Conv_MC','f4',('time','lat','lon'))

timen = ft.createVariable('Time', 'i4', 'time')
lifetimen = ft.createVariable('Lifetime', 'f4', ('time','lat','lon',))
longitude[:] = onedeglons
latitude[:] = onedeglats
##################################################################

#9504
#17520

for i in range(0,35064):
    if (i>0) & (i%24 == 0):
        k=0
        j=j+1
        ft.close()
        print('File ' +str(j-1)+ ' completed!')
        
        if year == 2008: #deal with leap year
            if (j>0) & (j%366 == 0):
                print('Year ' +str(year)+ ' completed!') 
                year+=1
                j=0
        if year != 2008:
            if (j>0) & (j%365 == 0):
                year+=1
                j=0
                print('Year ' +str(year)+ ' completed!')
                
        #import WINDdata for each timestep
        #load datasets
        fn_cf = '/net/seldon/disk2/Users/gah20/DCC/ISCCP/Total7_10.nc'
        ds_cf = nc.Dataset(fn_cf)
        ctp_dcc = ds_cf['DCC'][:,60:120]
        lentime = ds_cf['Time']
        
        fu2 = nc.Dataset('/net/seldon/disk2/Data/ECMWF/ERA5/U-wind-component_200hPa_0.25grid_timed/'+str(year)+'/U-wind-component_200hPa_'+str(j+1).zfill(3)+'.nc')
        fv2 = nc.Dataset('/net/seldon/disk2/Data/ECMWF/ERA5/V-wind-component_200hPa_0.25grid_timed/'+str(year)+'/V-wind-component_200hPa_'+str(j+1).zfill(3)+'.nc')
        u200 = fu2['u']
        v200 = fv2['v']
        
        fu3 = nc.Dataset('/net/seldon/disk2/Data/ECMWF/ERA5/U-wind-component_300hPa_0.25grid_timed/'+str(year)+'/U-wind-component_300hPa_'+str(j+1).zfill(3)+'.nc')
        fv3 = nc.Dataset('/net/seldon/disk2/Data/ECMWF/ERA5/V-wind-component_300hPa_0.25grid_timed/'+str(year)+'/V-wind-component_300hPa_'+str(j+1).zfill(3)+'.nc')
        u300 = fu3['u']
        v300 = fv3['v']
        
        #fu4 = nc.Dataset('/net/seldon/disk2/Data/ECMWF/ERA5/U-wind-component_400hPa_0.25grid_timed/'+str(year)+'/U-wind-component_400hPa_'+str(j+1).zfill(3)+'.nc')
        #fv4 = nc.Dataset('/net/seldon/disk2/Data/ECMWF/ERA5/V-wind-component_400hPa_0.25grid_timed/'+str(year)+'/V-wind-component_400hPa_'+str(j+1).zfill(3)+'.nc')
        #u400 = fu4['u']
        #v400 = fv4['v']
        #timerec = fv4['time']
        
        ft = nc.Dataset('TSC/'+str(year)+'/TSC_'+(str(j).zfill(3))+'.nc','w',format='NETCDF4')

        ft.createDimension('lat', 64)
        ft.createDimension('lon', 360)
        ft.createDimension('time', None)

        latitude = ft.createVariable('Latitude', 'f4', 'lat')
        longitude = ft.createVariable('Longitude', 'f4', 'lon')
        ntsc = ft.createVariable('TSC', 'f4', ('time', 'lat', 'lon',))
        nconvorigin = ft.createVariable('ConvOrigin','f4',('time','lat','lon',))
        nconvaod = ft.createVariable('ConvAOD','f4',('time','lat','lon',))
        nconvpc = ft.createVariable('ConvPc','f4',('time','lat','lon'))
        nconvcir = ft.createVariable('ConvCir','f4',('time','lat','lon'))
        #nconvsa = ft.createVariable('Conv_SA','f4',('time','lat','lon'))
        #nconvaf = ft.createVariable('Conv_AF','f4',('time','lat','lon'))
        #nconvmc = ft.createVariable('Conv_MC','f4',('time','lat','lon'))

        timen = ft.createVariable('Time', 'i4', 'time')
        lifetimen = ft.createVariable('Lifetime','f4',('time','lat','lon',))
        longitude[:] = onedeglons
        latitude[:] = onedeglats

    #create DCC flag
    dcc_norm = (ctp_dcc[round_three(i)]/ctp_dcc[round_three(i)])
    dcc_scale =np.kron(dcc_norm,np.ones((10,10)))
    newconv[20:620] = dcc_scale/dcc_scale
    aod4_scale = np.kron(aodgr[6+i],np.ones((10,10)))

    
    oceanconv[20:620] = np.where(mask==1,newconv[20:620],np.nan)
    landconv[20:620] = np.where(mask==0,newconv[20:620],np.nan)
    #conv_nsa[20:620] = np.where(regionalflag[20:620]==1,np.nan,newconv[20:620])
    #conv_sa[20:620] = np.where(regionalflag[20:620]==1,newconv[20:620],np.nan)
    #conv_af[20:620] = np.where(regionalflag[20:620]==2,newconv[20:620],np.nan)
    #conv_mc[20:620] = np.where(regionalflag[20:620]==3,newconv[20:620],np.nan)
    #conv_naf[20:620] = np.where(regionalflag[20:620]==2,np.nan,newconv[20:620])
    #conv_nmc[20:620] = np.where(regionalflag[20:620]==3,np.nan,newconv[20:620])

    
    h_aodconv[20:620] = np.where(aod4_scale>0.3,newconv[20:620],np.nan)
    l_aodconv[20:620] = np.where(aod4_scale<0.3,newconv[20:620],np.nan)

    pc = np.kron(pcraw[round_three(i),60:120],np.ones((10,10)))
    highpc[20:620] = np.where(pc>180,newconv[20:620],np.nan)
    lowpc[20:620] = np.where(pc<130,newconv[20:620],np.nan)

    #create insitu v real flag
    cf = np.nansum(np.nansum(taupc[round_three(i),:,:3],axis=0),axis=0)
    cf = np.where(cf<0,np.nan,cf)
    cirrus = np.where(cf[60:120]>10,1,0)
    cirrus = np.where(np.isnan(cf[60:120]),np.nan,cirrus)
    cirrus_scale = np.kron(cirrus,np.ones((10,10)))
    
    #advect
    trace = tracer_grid(k,2,itrace_lat, itrace_lon)


    tsc[~np.isnan(newconv)] = 0
    convorigin[~np.isnan(oceanconv)] = -100
    convorigin[~np.isnan(landconv)] = 100
    convaod[~np.isnan(h_aodconv)] = 100
    convaod[~np.isnan(l_aodconv)] = np.nan
    convpc[~np.isnan(highpc)] = 100
    convpc[~np.isnan(lowpc)] = -100
    convcir[~np.isnan(newconv)] = 100
    #convcir[~np.isnan(insitu)] = -100
    
    #convsa[~np.isnan(conv_sa)] = 100
    #convsa[~np.isnan(conv_nsa)] = -100
    #convaf[~np.isnan(conv_af)] = 100
    #convaf[~np.isnan(conv_naf)] = -100
    #convmc[~np.isnan(conv_mc)] = 100
    #convmc[~np.isnan(conv_nmc)] = -100



    #oldconv = round_to_index(trace[0,0],trace[1,0],tsc)
    advectconv = round_to_index_REDUCE(trace[0,1],trace[1,1],tsc)
    advectconvorigin = round_to_index_REDUCE(trace[0,1],trace[1,1],convorigin)
    advectconvaod = round_to_index_REDUCE(trace[0,1],trace[1,1],convaod)
    advectconvpc = round_to_index_REDUCE(trace[0,1],trace[1,1],convpc)
    advectconvcir = round_to_index_REDUCE(trace[0,1],trace[1,1],convcir)
    #advectconvaf = round_to_index_REDUCE(trace[0,1],trace[1,1],convaf)
    #advectconvsa = round_to_index_REDUCE(trace[0,1],trace[1,1],convsa)
    #advectconvmc = round_to_index_REDUCE(trace[0,1],trace[1,1],convmc)

    #if currently detrained but now CF<0.1, switch to detrained
    convlifetime[20:620] = np.where((cirrus_scale==0)&(advectconvcir[20:620]>0),advectconv[20:620],np.nan)
    advectconvcir[20:620] = np.where((cirrus_scale==0)&(advectconvcir[20:620]>0),-100,advectconvcir[20:620])
    #convlifetime[20:620] = np.where((cirrus_scale==0)&(advectconvcir[20:620]>0),100,0)
    #convlifetime[20:620] = np.where((cirrus_scale==0)&(advectconvcir[20:620]>0),5,0)
    print(convlifetime)
    ## at this point, note down the TSC of the locations where this happens 
    ### SAVE EACH TSC TO A NETCDF FILE HERE;
    
    timen[k] = timerec[k]
    tsc1deg = reduceres1deg(tsc)
    lifetimen[k] = reduceres1deg(convlifetime)
    ntsc[k] = tsc1deg
    
    convorigin1deg = reduceres1deg(convorigin)
    convaod1deg = reduceres1deg(convaod)
    convpc1deg = reduceres1deg(convpc)
    convcir1deg = reduceres1deg(convcir)
    #convsa1deg = reduceres1deg(convsa)
    #convaf1deg = reduceres1deg(convaf)
    #convmc1deg = reduceres1deg(convmc)

    nconvaod[k] = convaod1deg
    nconvorigin[k] = convorigin1deg
    nconvpc[k] = convpc1deg
    nconvcir[k] = convcir1deg
    #nconvsa[k] = convsa1deg
    #nconvaf[k] = convaf1deg
    #nconvmc[k] = convmc1deg
    
    k=k+1 ##reset array position for nc file to 0 when new file is created

    tsc = advectconv + 1
    convorigin = advectconvorigin
    convaod = advectconvaod
    convpc = advectconvpc
    convcir = advectconvcir
    #convsa = advectconvsa
    #convaf = advectconvaf
    #convmc = advectconvmc

    tsc_filler = convolve(tsc)
    conv_filler = convolve(convorigin)
    convaod_filler = convolve(convaod)
    convpc_filler = convolve(convpc)
    convcir_filler = convolve(convcir)

    #convaf_filler = convolve(convaf)
    #convsa_filler = convolve(convsa)
    #convmc_filler = convolve(convmc)
    
    tsc = np.where(~np.isfinite(tsc),tsc_filler,tsc)
    convorigin = np.where(~np.isfinite(convorigin),conv_filler,convorigin)
    convaod = np.where(~np.isfinite(convaod),convaod_filler,convaod)
    convpc = np.where(~np.isfinite(convpc),convpc_filler,convpc)
    convcir = np.where(~np.isfinite(convcir),convcir_filler,convcir)

    #convaf = np.where(~np.isfinite(convaf),convaf_filler,convaf)
    #convsa = np.where(~np.isfinite(convsa),convsa_filler,convsa)
    #convmc = np.where(~np.isfinite(convmc),convmc_filler,convmc)

    if (i % 10) ==0:
        print(str(i) +' done')
ft.close()
print('Time elapsed: ' + str(time.time() - start_time)/3600)
