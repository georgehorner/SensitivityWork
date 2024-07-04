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
#from mpl_toolkits.axes_grid1 import make
#_axes_locatable
from scipy import stats
import multiprocessing as mp
import matplotlib.patches as patches


import warnings
print('is it doing?')
timelen=499
yedges_net = np.arange(-650,650,1)
yedges_lw = np.arange(-500,0,1)
yedges_sw = np.arange(-100,900,1)
yedges_solar = np.arange(0,1400,5)
xedges = np.arange(0,timelen+1,1)

cereshour = 8760
k=0


for year in range(2000,2017,1):
    isccp = xr.open_mfdataset('/net/hardin/disk1/Data/ISCCP/access/isccp-basic/hgg/'+str(year)+'*/'+str(year)+'*.nc',combine='nested',concat_dim='time')
    taupc = isccp['n_pctaudist'][:]
    levtau = isccp['levtau'][:]
    levpc = isccp['levpc'][:]
    tauraw = isccp['tau'][:]
    pcraw = isccp['pc'][:]
    isccptime = isccp['time'][:]

    tsc = xr.open_mfdataset('/disk1/Users/gah20/TSC_JUNE/'+str(year)+'/TSC*.nc',combine='nested',concat_dim='time')['TSC'][:,2:62]
    landocean = xr.open_mfdataset('/disk1/Users/gah20/TSC_JUNE/'+str(year)+'/TSC*.nc',combine='nested',concat_dim='time')['ConvOrigin'][:,2:62]
    hlpc = xr.open_mfdataset('/disk1/Users/gah20/TSC_JUNE/'+str(year)+'/TSC*.nc',combine='nested',concat_dim='time')['ConvPc'][:,2:62]
    #aodtsc = xr.open_mfdataset('/disk1/Users/gah20/TSC_FULL/'+str(year)+'/TSC*.nc',combine='nested',concat_dim='time')['ConvAOD'][:,2:62]
    detinsit = xr.open_mfdataset('/disk1/Users/gah20/TSC_JUNE/'+str(year)+'/TSC*.nc',combine='nested',concat_dim='time')['ConvCir'][:,2:62]

    ds_w = xr.open_mfdataset('/disk1/Data/CERES/CERES_year'+str(year)+'.nc')
    swall = ds_w['toa_sw_all_1h'][:,60:120]
    swclr = ds_w['toa_sw_clr_1h'][:,60:120]
    lwall = ds_w['toa_lw_all_1h'][:,60:120]
    lwclr = ds_w['toa_lw_clr_1h'][:,60:120]

    ds_alb = xr.open_mfdataset('/disk1/Data/CERES/2024_dwld/CERES_year'+str(year)+'.nc')
    solar = ds_alb['toa_solar_all_1h'][:,60:120]
    alb_all = ds_alb['toa_alb_all_1h'][:,60:120]
    alb_clr = ds_alb['toa_alb_clr_1h'][:,60:120]
    timeceres = ds_alb['time'][:]

    offset = np.argmin(np.abs(isccp['time'][:].values - ds_w['time'][0].values))
    if len(isccptime[offset:]) == len(tsc[offset*3:])/3:

        swcldycalc = -(swall[0::3][0:len(isccp['time'][offset:])] - swclr[0::3][0:len(isccp['time'][offset:])])



        meanlw = np.nanmean(lwclr[:],axis=0)
        levels = 3
        lowcldamt = 1

        lowcloudraw = np.nansum(isccp['n_pcdist'][offset:,:levels,60:120],axis=1)
        lowcloud_day = np.where(lowcloudraw<0,np.nan,lowcloudraw)
        #lowcloud_day = np.where(lowcloud_day[:]>100,0,lowcloud_day[:])

        solar_t = solar[0::3][0:len(isccp['time'][offset:])]
        alb_all_t = alb_all[0::3][0:len(isccp['time'][offset:])]
        alb_clr_t = alb_clr[0::3][0:len(isccp['time'][offset:])]

        alb_all_t_sm = np.where(solar_t<700,np.nan,alb_all_t) 
        alb_clr_t_sm = np.where(solar_t<700,np.nan,alb_clr_t)
        solar_t_sm = np.where(solar_t<700,np.nan,solar_t)

        low_alb = np.where(lowcloud_day[:]<1,alb_all_t_sm - alb_clr_t_sm,np.nan)
        bkg_alb = np.where(lowcloud_day[:]<1,alb_all_t_sm,np.nan)


        meanbkg = np.zeros((len(alb_all_t),60,360))
        for i in range(0,12):
            meanbkg[i*238:238*(i+1)] = np.nanmean(bkg_alb[i*238:238*(i+1),:,:],axis=0)

        from csat2 import misc
        lstbkgalb = np.zeros((345,25,60,360))
        transpose_bkgalb = np.transpose(np.asarray(bkg_alb), (1, 2, 0))
        for j in range(0,345,1):
            try:
                for i in range(0,25,3):
                    lstbkgalb[j,i] = misc.time.toLocalSolarTime(i,np.arange(0,25,3),np.arange(0,360,1),np.asarray(transpose_bkgalb[:,:,j*8:(j*8)+9]),interpolation='nearest')
            except:
                lstbkgalb[j,i] = lstbkgalb[j-1,i]
            print(j)
        lstbkgalb = np.where(lstbkgalb==0,np.nan,lstbkgalb)

        alb_high_low = alb_all_t - (alb_clr_t + np.nanmean(low_alb,axis=0))
        alb_high_low = np.where(np.isnan(alb_high_low),0,alb_high_low)
        #alb_high_bkg = alb_all_t - meanbkg
        alb_high_bkg = alb_all_t - np.nanmean(np.nanmean(lstbkgalb,axis=0),axis=0)
        alb_high_bkg = np.where(np.isnan(alb_high_bkg),0,alb_high_bkg)

        lwclr_f = np.where(np.isnan(lwclr[0::3][0:len(isccp['time'][offset:])]),meanlw,lwclr[0::3][0:len(isccp['time'][offset:])])
        lwcldy = (lwall[0::3][0:len(isccp['time'][offset:])] - lwclr_f) 
        swcldy = solar_t*alb_high_low
        swcldy_high = solar_t*alb_high_bkg
        swcldy_high_zero = np.where(swcldy_high<0,0,swcldy_high)
        

        tsc_d = np.asarray(np.where(detinsit[offset*3::3]>=0,tsc[offset*3::3],np.nan))
        #tsc_d = np.asarray(np.where(landocean[offset*3::3]>0,tsc_d,np.nan))

        tsc_d_1 = np.asarray(np.where((hlpc[offset*3::3]<-750),tsc_d,np.nan))
        tsc_d_2 = np.asarray(np.where((hlpc[offset*3::3]<0)&(hlpc[offset*3::3]>-750),tsc_d,np.nan))
        tsc_d_3 = np.asarray(np.where((hlpc[offset*3::3]<750)&(hlpc[offset*3::3]>0),tsc_d,np.nan))
        tsc_d_4 = np.asarray(np.where((hlpc[offset*3::3]>750),tsc_d,np.nan))

        H_lw_d_1, xedges, yedges_lw = np.histogram2d(tsc_d_1.flatten(),np.asarray(lwcldy).flatten(), bins=(xedges,yedges_lw))
        H_lw_d_2, xedges, yedges_lw = np.histogram2d(tsc_d_2.flatten(),np.asarray(lwcldy).flatten(), bins=(xedges,yedges_lw))
        H_lw_d_3, xedges, yedges_lw = np.histogram2d(tsc_d_3.flatten(),np.asarray(lwcldy).flatten(), bins=(xedges,yedges_lw))
        H_lw_d_4, xedges, yedges_lw = np.histogram2d(tsc_d_4.flatten(),np.asarray(lwcldy).flatten(), bins=(xedges,yedges_lw))

        H_sw_d_1, xedges, yedges_sw = np.histogram2d(tsc_d_1.flatten(),np.asarray(swcldy_high).flatten(), bins=(xedges,yedges_sw))
        H_sw_d_2, xedges, yedges_sw = np.histogram2d(tsc_d_2.flatten(),np.asarray(swcldy_high).flatten(), bins=(xedges,yedges_sw))
        H_sw_d_3, xedges, yedges_sw = np.histogram2d(tsc_d_3.flatten(),np.asarray(swcldy_high).flatten(), bins=(xedges,yedges_sw))
        H_sw_d_4, xedges, yedges_sw = np.histogram2d(tsc_d_4.flatten(),np.asarray(swcldy_high).flatten(), bins=(xedges,yedges_sw))

        print('histograms done')


        H_norm_lw_1_det = np.zeros((timelen,499))
        H_norm_lw_2_det = np.zeros((timelen,499))
        H_norm_lw_3_det = np.zeros((timelen,499))
        H_norm_lw_4_det = np.zeros((timelen,499))

        H_norm_sw_1_det = np.zeros((timelen,999))
        H_norm_sw_2_det = np.zeros((timelen,999))
        H_norm_sw_3_det = np.zeros((timelen,999))
        H_norm_sw_4_det = np.zeros((timelen,999))



        for i in range(timelen):

            H_norm_sw_1_det[i] = H_sw_d_1[i] / np.sum(H_sw_d_1[i])
            H_norm_sw_2_det[i] = H_sw_d_2[i] / np.sum(H_sw_d_2[i])
            H_norm_sw_3_det[i] = H_sw_d_3[i] / np.sum(H_sw_d_3[i])
            H_norm_sw_4_det[i] = H_sw_d_4[i] / np.sum(H_sw_d_4[i])
            
            H_norm_lw_1_det[i] = H_lw_d_1[i] / np.sum(H_lw_d_1[i])
            H_norm_lw_2_det[i] = H_lw_d_2[i] / np.sum(H_lw_d_2[i])
            H_norm_lw_3_det[i] = H_lw_d_3[i] / np.sum(H_lw_d_3[i])
            H_norm_lw_4_det[i] = H_lw_d_4[i] / np.sum(H_lw_d_4[i])


        print('doing averages')


        LW_AVERAGE_1_det = np.zeros((timelen))
        LW_AVERAGE_2_det = np.zeros((timelen))
        LW_AVERAGE_3_det = np.zeros((timelen))
        LW_AVERAGE_4_det = np.zeros((timelen))

        SW_AVERAGE_1_det = np.zeros((timelen))
        SW_AVERAGE_2_det = np.zeros((timelen))
        SW_AVERAGE_3_det = np.zeros((timelen))
        SW_AVERAGE_4_det = np.zeros((timelen))

        for i in range(timelen):

        
            LW_AVERAGE_1_det[i] = np.average(np.arange(499.5,0.5,-1), weights=H_norm_lw_1_det[i])
            LW_AVERAGE_2_det[i] = np.average(np.arange(499.5,0.5,-1), weights=H_norm_lw_2_det[i])
            LW_AVERAGE_3_det[i] = np.average(np.arange(499.5,0.5,-1), weights=H_norm_lw_3_det[i])
            LW_AVERAGE_4_det[i] = np.average(np.arange(499.5,0.5,-1), weights=H_norm_lw_4_det[i])

            SW_AVERAGE_1_det[i] = np.average(np.arange(99.5,-899.5,-1), weights=H_norm_sw_1_det[i])
            SW_AVERAGE_2_det[i] = np.average(np.arange(99.5,-899.5,-1), weights=H_norm_sw_2_det[i])
            SW_AVERAGE_3_det[i] = np.average(np.arange(99.5,-899.5,-1), weights=H_norm_sw_3_det[i])
            SW_AVERAGE_4_det[i] = np.average(np.arange(99.5,-899.5,-1), weights=H_norm_sw_4_det[i])
            
        print('averages done')

        counts = np.histogram(np.asarray(tsc).flatten(),bins=np.arange(0,500,1))[0]
        counts_det = np.histogram(np.asarray(tsc_d).flatten(),bins=np.arange(0,500,1))[0]
        counts_det_1 = np.histogram(np.asarray(tsc_d_1).flatten(),bins=np.arange(0,500,1))[0]
        counts_det_2 = np.histogram(np.asarray(tsc_d_2).flatten(),bins=np.arange(0,500,1))[0]
        counts_det_3 = np.histogram(np.asarray(tsc_d_3).flatten(),bins=np.arange(0,500,1))[0]
        counts_det_4 = np.histogram(np.asarray(tsc_d_4).flatten(),bins=np.arange(0,500,1))[0]
        

        LWDETsum_1 = np.nansum(LW_AVERAGE_1_det[:]*(counts_det_1[:]))/np.nansum(counts_det_1[:])
        SWDETsum_1 = np.nansum(SW_AVERAGE_1_det[:]*(counts_det_1[:]))/np.nansum(counts_det_1[:])
        LWDETsum_2 = np.nansum(LW_AVERAGE_2_det[:]*(counts_det_2[:]))/np.nansum(counts_det_2[:])
        SWDETsum_2 = np.nansum(SW_AVERAGE_2_det[:]*(counts_det_2[:]))/np.nansum(counts_det_2[:])
        LWDETsum_3 = np.nansum(LW_AVERAGE_3_det[:]*(counts_det_3[:]))/np.nansum(counts_det_3[:])
        SWDETsum_3 = np.nansum(SW_AVERAGE_3_det[:]*(counts_det_3[:]))/np.nansum(counts_det_3[:])
        LWDETsum_4 = np.nansum(LW_AVERAGE_4_det[:]*(counts_det_4[:]))/np.nansum(counts_det_4[:])
        SWDETsum_4 = np.nansum(SW_AVERAGE_4_det[:]*(counts_det_4[:]))/np.nansum(counts_det_4[:])

        CRE_TSC = np.zeros((2,4,499))

        CRE_TSC[0,0] = LW_AVERAGE_1_det
        CRE_TSC[0,1] = LW_AVERAGE_2_det
        CRE_TSC[0,2] = LW_AVERAGE_3_det
        CRE_TSC[0,3] = LW_AVERAGE_4_det
        CRE_TSC[1,0] = SW_AVERAGE_1_det
        CRE_TSC[1,1] = SW_AVERAGE_2_det
        CRE_TSC[1,2] = SW_AVERAGE_3_det
        CRE_TSC[1,3] = SW_AVERAGE_4_det
        np.save('/disk1/Users/gah20/SensitivityWork/values/CRES/CRE_TSC_'+str(year),CRE_TSC)
        #np.save('/disk1/Users/gah20/SensitivityWork/values/CRES/CRE_TSC_ocean_2017',CRE_TSC)

        forcings = np.zeros((2,4))
        forcings[0,0] = LWDETsum_1
        forcings[0,1] = LWDETsum_2
        forcings[0,2] = LWDETsum_3
        forcings[0,3] = LWDETsum_4

        forcings[1,0] = SWDETsum_1
        forcings[1,1] = SWDETsum_2
        forcings[1,2] = SWDETsum_3
        forcings[1,3] = SWDETsum_4

        print(str(year)+ ' forcings are: '+str(forcings))
        np.save('/disk1/Users/gah20/SensitivityWork/values/Forcings/forcings_'+str(year),forcings)