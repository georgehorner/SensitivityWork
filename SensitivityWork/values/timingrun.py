import netCDF4 as nc
import xarray as xr
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import LinearSegmentedColormap
import cartopy.crs as ccrs
import os
import math
import csat2
#import imageio
from scipy.ndimage.interpolation import map_coordinates
from pprint import pprint
import glob
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator
from copy import copy
from copy import deepcopy

from scipy import stats
import multiprocessing as mp
import matplotlib.patches as patches

_wbgyr_cdict = {'red': ((0.0,  255./255, 255./255),
                        (0.125, 173./255, 173./255),
                        (0.25,  95./255,  95./255),
                        (0.375, 73./255,  73./255),
                        (0.5,  165./255, 164./255),
                        (0.625, 248./255, 248./255),
                        (0.75, 236./255, 236./255),
                        (0.875, 200./255, 200./255),
                        (1.0,  146./255, 146./255)),
                'green': ((0.0,  255./255, 255./255),
                          (0.125, 224./255, 224./255),
                          (0.25, 163./255, 163./255),
                          (0.375, 166./255, 166./255),
                          (0.5,  207./255, 207./255),
                          (0.625, 184./255, 184./255),
                          (0.75,  86./255,  86./255),
                          (0.875, 29./255,  29./255),
                          (1.0,   21./255,  21./255)),
                'blue': ((0.0,  255./255, 255./255),
                         (0.125, 248./255, 248./255),
                         (0.25, 214./255, 214./255),
                         (0.375, 120./255, 120./255),
                         (0.5,   81./255,  81./255),
                         (0.625, 73./255,  73./255),
                         (0.75,  41./255,  41./255),
                         (0.875, 38./255,  38./255),
                         (1.0,   25./255,  25./255))}
# And reverse
_wbgyr_cdict_r = deepcopy(_wbgyr_cdict)
for i in _wbgyr_cdict_r.keys():
    _wbgyr_cdict_r[i] = [(1-j[0], j[1], j[2]) for j in _wbgyr_cdict_r[i]]
    _wbgyr_cdict_r[i].reverse()

cmap = LinearSegmentedColormap('WBGYR', _wbgyr_cdict)
cmap.set_bad('#D2D2D2')
plt.register_cmap(cmap=cmap)

cmap_r = LinearSegmentedColormap('WBGYR', _wbgyr_cdict_r)
cmap_r.set_bad('#D2D2D2')
plt.register_cmap(cmap=cmap_r)


def circular_mean(data):
    # Convert times to angles in radians
    angles = data * (2 * np.pi / 24)

    # Compute the mean of sine and cosine components
    mean_sin = np.nanmean(np.sin(angles), axis=0)
    mean_cos = np.nanmean(np.cos(angles), axis=0)

    # Compute the circular mean angle
    mean_angle = np.arctan2(mean_sin, mean_cos)

    # Convert mean angle back to time
    mean_time = mean_angle * (24 / (2 * np.pi))

    # Ensure times are within the range [0, 24]
    mean_time = np.mod(mean_time, 24)

    mean_time
    return mean_time


lst_conv = np.zeros((2,5,60,360))

for year in range(2000,2017):

    tsc = xr.open_mfdataset('/net/hardin/disk1/Users/gah20/TSC_JUNE/'+str(year)+'/*.nc',combine='nested',concat_dim='time')

    timearray = np.asarray(np.repeat(tsc['Time']%24,360)).reshape(len(tsc['Time']),360)
    for j in range(0,100,1):
        for i in range(0,360,1):
            timearray[j,i] = csat2.misc.time.utc_to_lst(timearray[j,i],i)
    timearray = np.repeat(timearray,60).reshape(len(tsc['Time']),360,60)

    time_lst = np.zeros((len(tsc['Time'][:]),60,360))
    for i in range(len(timearray)):
        time_lst[i] = timearray[i].T

    allconv = np.where(tsc['TSC']==0,1,0)

    strongconv = np.where(tsc['ConvPc']<-750,allconv,np.nan)
    strongishconv = np.where((tsc['ConvPc']<0)&(tsc['ConvPc']>-750),allconv,np.nan)
    weakishconv = np.where((tsc['ConvPc']>0)&(tsc['ConvPc']<750),allconv,np.nan)
    weakconv = np.where((tsc['ConvPc']>750),allconv,np.nan)

    all_lst = np.where(allconv[:,2:62]==1,time_lst,np.nan)
    strong_lst = np.where(strongconv[:,2:62]==1,time_lst,np.nan)
    strongish_lst = np.where(strongishconv[:,2:62]==1,time_lst,np.nan)
    weakish_lst = np.where(weakishconv[:,2:62]==1,time_lst,np.nan)
    weak_lst = np.where(weakconv[:,2:62]==1,time_lst,np.nan)


    lst_conv[0,0] = circular_mean(all_lst)
    lst_conv[0,1] = circular_mean(strong_lst)
    lst_conv[0,2] = circular_mean(strongish_lst)
    lst_conv[0,3] = circular_mean(weakish_lst)
    lst_conv[0,4] = circular_mean(weak_lst)

    lst_conv[1,0] = np.nansum(allconv[:,2:62],axis=0)
    lst_conv[1,1] = np.nansum(strongconv[:,2:62],axis=0)
    lst_conv[1,2] = np.nansum(strongishconv[:,2:62],axis=0)
    lst_conv[1,3] = np.nansum(weakishconv[:,2:62],axis=0)
    lst_conv[1,4] = np.nansum(weakconv[:,2:62],axis=0)
    np.save('/net/hardin/disk1/Users/gah20/SensitivityWork/values/LSTs/'+str(year)+'_lst_conv.npy',lst_conv)
    print(year)
