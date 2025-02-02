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

for year in range(2014,2017,1):
    tsc = xr.open_mfdataset('/net/hardin/disk1/Users/gah20/TSC_2origin/'+str(year)+'/*.nc',combine='nested',concat_dim='time')
    detrained = xr.open_mfdataset('/net/hardin/disk1/Users/gah20/TSC_2origin/'+str(year)+'/*.nc',combine='nested',concat_dim='time')['ConvCir']

    loclat = tsc['LocOrigin'][:,0,2:62]
    loclon = tsc['LocOrigin'][:,1,2:62]



    lats_d = np.where(detrained[:,2:62]>0,loclat,np.nan)
    lons_d = np.where(detrained[:,2:62]>0,loclon,np.nan)

    ilat_d = np.floor(lats_d/10)
    ilon_d= np.floor(lons_d/10)
    

    heatmap_d = np.zeros((60,360,60,360))

    for i in range(0,60,1):
        for j in range(0,360,1):
            heatmap_d[i,j] = np.nansum(np.where((ilat_d[:]==i)&(ilon_d[:]==j),1,0),axis=0)
            if j%10==0:
                print(str(j) + ' done')
        print(str(i) + ' done!')
    np.save('/net/hardin/disk1/Users/gah20/Heatmaps/heatmap_1deg_d'+str(year)+'.npy',heatmap_d)
    print(str(year) + ' done')

