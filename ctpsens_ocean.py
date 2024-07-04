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


for year in range(2013,2017,1):
    isccp = xr.open_mfdataset('/net/hardin/disk1/Data/ISCCP/access/isccp-basic/hgg/'+str(year)+'*/'+str(year)+'*.nc',combine='nested',concat_dim='time')
    taupc = isccp['n_pctaudist']
    levtau = isccp['levtau']
    levpc = isccp['levpc']
    tauraw = isccp['tau']
    pcraw = isccp['pc']
    isccptime = isccp['time']

    tsc = xr.open_mfdataset('/disk1/Users/gah20/TSC_2origin/'+str(year)+'/TSC*.nc',combine='nested',concat_dim='time')['TSC'][:,2:62]
    landocean = xr.open_mfdataset('/disk1/Users/gah20/TSC_2origin/'+str(year)+'/TSC*.nc',combine='nested',concat_dim='time')['ConvOrigin'][:,2:62]
    hlpc = xr.open_mfdataset('/disk1/Users/gah20/TSC_2origin/'+str(year)+'/TSC*.nc',combine='nested',concat_dim='time')['ConvPc'][:,2:62]
    #aodtsc = xr.open_mfdataset('/disk1/Users/gah20/TSC_FULL/'+str(year)+'/TSC*.nc',combine='nested',concat_dim='time')['ConvAOD'][:,2:62]
    detinsit = xr.open_mfdataset('/disk1/Users/gah20/TSC_2origin/'+str(year)+'/TSC*.nc',combine='nested',concat_dim='time')['ConvCir'][:,2:62]

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

        alb_high_low = alb_all_t - (alb_clr_t + np.nanmean(low_alb,axis=0))
        alb_high_low = np.where(np.isnan(alb_high_low),0,alb_high_low)
        alb_high_bkg = alb_all_t - meanbkg
        alb_high_bkg = np.where(np.isnan(alb_high_bkg),0,alb_high_bkg)

        lwclr_f = np.where(np.isnan(lwclr[0::3][0:len(isccp['time'][offset:])]),meanlw,lwclr[0::3][0:len(isccp['time'][offset:])])
        lwcldy = (lwall[0::3][0:len(isccp['time'][offset:])] - lwclr_f) 
        swcldy = solar_t*alb_high_low
        swcldy_high = solar_t*alb_high_bkg
        swcldy_high_zero = np.where(swcldy_high<0,0,swcldy_high)

        tsc = np.where(landocean<0,tsc,np.nan)
        tsc_t = (np.asarray(tsc[offset*3::3]))
        tsc_d = np.asarray(np.where(detinsit[offset*3::3]>=0,tsc[offset*3::3],np.nan))
        tsc_i = np.asarray(np.where(np.isnan(tsc_d),tsc[offset*3::3],np.nan))

        tsc_ctp_1 = np.asarray(np.where(hlpc[offset*3::3]<-1000,tsc[offset*3::3],np.nan))
        tsc_ctp_2 = np.asarray(np.where((hlpc[offset*3::3]<0)&(hlpc[offset*3::3]>-1000),tsc[offset*3::3],np.nan))
        tsc_ctp_3 = np.asarray(np.where((hlpc[offset*3::3]<1000)&(hlpc[offset*3::3]>0),tsc[offset*3::3],np.nan))
        tsc_ctp_4 = np.asarray(np.where((hlpc[offset*3::3]>1000),tsc[offset*3::3],np.nan))

        tsc_d_1 = np.asarray(np.where(hlpc[offset*3::3]<-1000,tsc_d,np.nan))
        tsc_d_2 = np.asarray(np.where((hlpc[offset*3::3]<0)&(hlpc[offset*3::3]>-1000),tsc_d,np.nan))
        tsc_d_3 = np.asarray(np.where((hlpc[offset*3::3]<1000)&(hlpc[offset*3::3]>0),tsc_d,np.nan))
        tsc_d_4 = np.asarray(np.where((hlpc[offset*3::3]>1000),tsc_d,np.nan))

        tsc_i_1 = np.asarray(np.where(hlpc[offset*3::3]<-1000,tsc_i,np.nan))
        tsc_i_2 = np.asarray(np.where((hlpc[offset*3::3]<0)&(hlpc[offset*3::3]>-1000),tsc_i,np.nan))
        tsc_i_3 = np.asarray(np.where((hlpc[offset*3::3]<1000)&(hlpc[offset*3::3]>0),tsc_i,np.nan))
        tsc_i_4 = np.asarray(np.where((hlpc[offset*3::3]>1000),tsc_i,np.nan))


        #tsc_o_d = np.asarray(np.where(hlpc[offset*3::3]<-50,tsc_d,np.nan))
        #tsc_o_i = np.asarray(np.where(hlpc[offset*3::3]<-50,tsc_i,np.nan))

        #tsc_l_d = np.asarray(np.where(hlpc[offset*3::3]>50,tsc_d,np.nan))
        #tsc_l_i = np.asarray(np.where(hlpc[offset*3::3]>50,tsc_i,np.nan))

        H_lw, xedges, yedges_lw = np.histogram2d(tsc_t.flatten(),np.asarray(lwcldy).flatten(), bins=(xedges,yedges_lw)) 
        H_sw, xedges, yedges_sw = np.histogram2d(tsc_t.flatten(),np.asarray(swcldy_high).flatten(), bins=(xedges,yedges_sw))

        H_sw_zero, xedges, yedges_sw = np.histogram2d(tsc_t.flatten(),np.asarray(swcldy_high_zero).flatten(), bins=(xedges,yedges_sw))
        H_sw_det_zero, xedges, yedges_sw = np.histogram2d(tsc_d.flatten(),np.asarray(swcldy_high_zero).flatten(), bins=(xedges,yedges_sw))
        H_sw_ins_zero, xedges, yedges_sw = np.histogram2d(tsc_i.flatten(),np.asarray(swcldy_high_zero).flatten(), bins=(xedges,yedges_sw))

        H_sw_det, xedges, yedges_sw = np.histogram2d(tsc_d.flatten(),np.asarray(swcldy_high).flatten(), bins=(xedges,yedges_sw))
        H_sw_ins, xedges, yedges_sw = np.histogram2d(tsc_i.flatten(),np.asarray(swcldy_high).flatten(), bins=(xedges,yedges_sw))

        H_lw_det, xedges, yedges_lw = np.histogram2d(tsc_d.flatten(),np.asarray(lwcldy).flatten(), bins=(xedges,yedges_lw))
        H_lw_ins, xedges, yedges_lw = np.histogram2d(tsc_i.flatten(),np.asarray(lwcldy).flatten(), bins=(xedges,yedges_lw))

        H_sw_ctp_1, xedges, yedges_sw = np.histogram2d(tsc_ctp_1.flatten(),np.asarray(swcldy_high).flatten(), bins=(xedges,yedges_sw))
        H_sw_ctp_2, xedges, yedges_sw = np.histogram2d(tsc_ctp_2.flatten(),np.asarray(swcldy_high).flatten(), bins=(xedges,yedges_sw))
        H_sw_ctp_3, xedges, yedges_sw = np.histogram2d(tsc_ctp_3.flatten(),np.asarray(swcldy_high).flatten(), bins=(xedges,yedges_sw))
        H_sw_ctp_4, xedges, yedges_sw = np.histogram2d(tsc_ctp_4.flatten(),np.asarray(swcldy_high).flatten(), bins=(xedges,yedges_sw))

        H_lw_ctp_1, xedges, yedges_lw = np.histogram2d(tsc_ctp_1.flatten(),np.asarray(lwcldy).flatten(), bins=(xedges,yedges_lw))
        H_lw_ctp_2, xedges, yedges_lw = np.histogram2d(tsc_ctp_2.flatten(),np.asarray(lwcldy).flatten(), bins=(xedges,yedges_lw))
        H_lw_ctp_3, xedges, yedges_lw = np.histogram2d(tsc_ctp_3.flatten(),np.asarray(lwcldy).flatten(), bins=(xedges,yedges_lw))
        H_lw_ctp_4, xedges, yedges_lw = np.histogram2d(tsc_ctp_4.flatten(),np.asarray(lwcldy).flatten(), bins=(xedges,yedges_lw))

        H_lw_d_1, xedges, yedges_lw = np.histogram2d(tsc_d_1.flatten(),np.asarray(lwcldy).flatten(), bins=(xedges,yedges_lw))
        H_lw_d_2, xedges, yedges_lw = np.histogram2d(tsc_d_2.flatten(),np.asarray(lwcldy).flatten(), bins=(xedges,yedges_lw))
        H_lw_d_3, xedges, yedges_lw = np.histogram2d(tsc_d_3.flatten(),np.asarray(lwcldy).flatten(), bins=(xedges,yedges_lw))
        H_lw_d_4, xedges, yedges_lw = np.histogram2d(tsc_d_4.flatten(),np.asarray(lwcldy).flatten(), bins=(xedges,yedges_lw))

        H_lw_i_1, xedges, yedges_lw = np.histogram2d(tsc_i_1.flatten(),np.asarray(lwcldy).flatten(), bins=(xedges,yedges_lw))
        H_lw_i_2, xedges, yedges_lw = np.histogram2d(tsc_i_2.flatten(),np.asarray(lwcldy).flatten(), bins=(xedges,yedges_lw))
        H_lw_i_3, xedges, yedges_lw = np.histogram2d(tsc_i_3.flatten(),np.asarray(lwcldy).flatten(), bins=(xedges,yedges_lw))
        H_lw_i_4, xedges, yedges_lw = np.histogram2d(tsc_i_4.flatten(),np.asarray(lwcldy).flatten(), bins=(xedges,yedges_lw))

        H_sw_d_1, xedges, yedges_sw = np.histogram2d(tsc_d_1.flatten(),np.asarray(swcldy_high).flatten(), bins=(xedges,yedges_sw))
        H_sw_d_2, xedges, yedges_sw = np.histogram2d(tsc_d_2.flatten(),np.asarray(swcldy_high).flatten(), bins=(xedges,yedges_sw))
        H_sw_d_3, xedges, yedges_sw = np.histogram2d(tsc_d_3.flatten(),np.asarray(swcldy_high).flatten(), bins=(xedges,yedges_sw))
        H_sw_d_4, xedges, yedges_sw = np.histogram2d(tsc_d_4.flatten(),np.asarray(swcldy_high).flatten(), bins=(xedges,yedges_sw))

        H_sw_i_1, xedges, yedges_sw = np.histogram2d(tsc_i_1.flatten(),np.asarray(swcldy_high).flatten(), bins=(xedges,yedges_sw))
        H_sw_i_2, xedges, yedges_sw = np.histogram2d(tsc_i_2.flatten(),np.asarray(swcldy_high).flatten(), bins=(xedges,yedges_sw))
        H_sw_i_3, xedges, yedges_sw = np.histogram2d(tsc_i_3.flatten(),np.asarray(swcldy_high).flatten(), bins=(xedges,yedges_sw))
        H_sw_i_4, xedges, yedges_sw = np.histogram2d(tsc_i_4.flatten(),np.asarray(swcldy_high).flatten(), bins=(xedges,yedges_sw))

        print('histograms done')
        H_norm_lw_all = np.zeros((timelen,499))
        H_norm_sw_all = np.zeros((timelen,999))

        H_norm_lw_1 = np.zeros((timelen,499))
        H_norm_sw_1 = np.zeros((timelen,999))

        H_norm_lw_2 = np.zeros((timelen,499))
        H_norm_sw_2 = np.zeros((timelen,999))

        H_norm_lw_3 = np.zeros((timelen,499))
        H_norm_sw_3 = np.zeros((timelen,999))

        H_norm_lw_4 = np.zeros((timelen,499))
        H_norm_sw_4 = np.zeros((timelen,999))

        H_norm_sw_det = np.zeros((timelen,999))
        H_norm_sw_ins = np.zeros((timelen,999))
        H_norm_lw_det = np.zeros((timelen,499))
        H_norm_lw_ins = np.zeros((timelen,499))

        H_norm_lw_1_det = np.zeros((timelen,499))
        H_norm_lw_1_ins = np.zeros((timelen,499))
        H_norm_lw_2_det = np.zeros((timelen,499))
        H_norm_lw_2_ins = np.zeros((timelen,499))
        H_norm_lw_3_det = np.zeros((timelen,499))
        H_norm_lw_3_ins = np.zeros((timelen,499))
        H_norm_lw_4_det = np.zeros((timelen,499))
        H_norm_lw_4_ins = np.zeros((timelen,499))

        H_norm_sw_1_det = np.zeros((timelen,999))
        H_norm_sw_1_ins = np.zeros((timelen,999))
        H_norm_sw_2_det = np.zeros((timelen,999))
        H_norm_sw_2_ins = np.zeros((timelen,999))
        H_norm_sw_3_det = np.zeros((timelen,999))
        H_norm_sw_3_ins = np.zeros((timelen,999))
        H_norm_sw_4_det = np.zeros((timelen,999))
        H_norm_sw_4_ins = np.zeros((timelen,999))



        for i in range(timelen):

            H_norm_lw_all[i] = H_lw[i] / np.sum(H_lw[i])
            H_norm_sw_all[i] = H_sw[i] / np.sum(H_sw[i])

            H_norm_lw_1[i] = H_lw_ctp_1[i] / np.sum(H_lw_ctp_1[i])
            H_norm_sw_1[i] = H_sw_ctp_1[i] / np.sum(H_sw_ctp_1[i])

            H_norm_lw_2[i] = H_lw_ctp_2[i] / np.sum(H_lw_ctp_2[i])
            H_norm_sw_2[i] = H_sw_ctp_2[i] / np.sum(H_sw_ctp_2[i])

            H_norm_lw_3[i] = H_lw_ctp_3[i] / np.sum(H_lw_ctp_3[i])
            H_norm_sw_3[i] = H_sw_ctp_3[i] / np.sum(H_sw_ctp_3[i])

            H_norm_lw_4[i] = H_lw_ctp_4[i] / np.sum(H_lw_ctp_4[i])
            H_norm_sw_4[i] = H_sw_ctp_4[i] / np.sum(H_sw_ctp_4[i])

            H_norm_sw_det[i] = H_sw_det[i] / np.sum(H_sw_det[i])
            H_norm_sw_ins[i] = H_sw_ins[i] / np.sum(H_sw_ins[i])
            H_norm_lw_det[i] = H_lw_det[i] / np.sum(H_lw_det[i])
            H_norm_lw_ins[i] = H_lw_ins[i] / np.sum(H_lw_ins[i])


            H_norm_sw_1_det[i] = H_sw_d_1[i] / np.sum(H_sw_d_1[i])
            H_norm_sw_1_ins[i] = H_sw_i_1[i] / np.sum(H_sw_i_1[i])
            H_norm_sw_2_det[i] = H_sw_d_2[i] / np.sum(H_sw_d_2[i])
            H_norm_sw_2_ins[i] = H_sw_i_2[i] / np.sum(H_sw_i_2[i])
            H_norm_sw_3_det[i] = H_sw_d_3[i] / np.sum(H_sw_d_3[i])
            H_norm_sw_3_ins[i] = H_sw_i_3[i] / np.sum(H_sw_i_3[i])
            H_norm_sw_4_det[i] = H_sw_d_4[i] / np.sum(H_sw_d_4[i])
            H_norm_sw_4_ins[i] = H_sw_i_4[i] / np.sum(H_sw_i_4[i])
            
            H_norm_lw_1_det[i] = H_lw_d_1[i] / np.sum(H_lw_d_1[i])
            H_norm_lw_1_ins[i] = H_lw_i_1[i] / np.sum(H_lw_i_1[i])
            H_norm_lw_2_det[i] = H_lw_d_2[i] / np.sum(H_lw_d_2[i])
            H_norm_lw_2_ins[i] = H_lw_i_2[i] / np.sum(H_lw_i_2[i])
            H_norm_lw_3_det[i] = H_lw_d_3[i] / np.sum(H_lw_d_3[i])
            H_norm_lw_3_ins[i] = H_lw_i_3[i] / np.sum(H_lw_i_3[i])
            H_norm_lw_4_det[i] = H_lw_d_4[i] / np.sum(H_lw_d_4[i])
            H_norm_lw_4_ins[i] = H_lw_i_4[i] / np.sum(H_lw_i_4[i])


        print('doing averages')
        LW_AVERAGE_all = np.zeros((timelen))
        SW_AVERAGE_all = np.zeros((timelen))

        LW_AVERAGE_1 = np.zeros((timelen))
        SW_AVERAGE_1 = np.zeros((timelen))

        LW_AVERAGE_2 = np.zeros((timelen))
        SW_AVERAGE_2 = np.zeros((timelen))

        LW_AVERAGE_3 = np.zeros((timelen))
        SW_AVERAGE_3 = np.zeros((timelen))

        LW_AVERAGE_4 = np.zeros((timelen))
        SW_AVERAGE_4 = np.zeros((timelen))


        SW_AVERAGE_det = np.zeros((timelen))
        SW_AVERAGE_ins = np.zeros((timelen))

        LW_AVERAGE_det = np.zeros((timelen))
        LW_AVERAGE_ins = np.zeros((timelen))

        LW_AVERAGE_1_det = np.zeros((timelen))
        LW_AVERAGE_1_ins = np.zeros((timelen))
        LW_AVERAGE_2_det = np.zeros((timelen))
        LW_AVERAGE_2_ins = np.zeros((timelen))
        LW_AVERAGE_3_det = np.zeros((timelen))
        LW_AVERAGE_3_ins = np.zeros((timelen))
        LW_AVERAGE_4_det = np.zeros((timelen))
        LW_AVERAGE_4_ins = np.zeros((timelen))

        SW_AVERAGE_1_det = np.zeros((timelen))
        SW_AVERAGE_1_ins = np.zeros((timelen))
        SW_AVERAGE_2_det = np.zeros((timelen))
        SW_AVERAGE_2_ins = np.zeros((timelen))
        SW_AVERAGE_3_det = np.zeros((timelen))
        SW_AVERAGE_3_ins = np.zeros((timelen))
        SW_AVERAGE_4_det = np.zeros((timelen))
        SW_AVERAGE_4_ins = np.zeros((timelen))

        for i in range(timelen):

            LW_AVERAGE_all[i] = np.average(np.arange(499.5,0.5,-1), weights=H_norm_lw_all[i])
            SW_AVERAGE_all[i] = np.average(np.arange(99.5,-899.5,-1), weights=H_norm_sw_all[i])

            LW_AVERAGE_1[i] = np.average(np.arange(499.5,0.5,-1), weights=H_norm_lw_1[i])
            SW_AVERAGE_1[i] = np.average(np.arange(99.5,-899.5,-1), weights=H_norm_sw_1[i])
            LW_AVERAGE_2[i] = np.average(np.arange(499.5,0.5,-1), weights=H_norm_lw_2[i])
            SW_AVERAGE_2[i] = np.average(np.arange(99.5,-899.5,-1), weights=H_norm_sw_2[i])
            LW_AVERAGE_3[i] = np.average(np.arange(499.5,0.5,-1), weights=H_norm_lw_3[i])
            SW_AVERAGE_3[i] = np.average(np.arange(99.5,-899.5,-1), weights=H_norm_sw_3[i])
            LW_AVERAGE_4[i] = np.average(np.arange(499.5,0.5,-1), weights=H_norm_lw_4[i])
            SW_AVERAGE_4[i] = np.average(np.arange(99.5,-899.5,-1), weights=H_norm_sw_4[i])

            SW_AVERAGE_det[i] = np.average(np.arange(99.5,-899.5,-1), weights=H_norm_sw_det[i])
            SW_AVERAGE_ins[i] = np.average(np.arange(99.5,-899.5,-1), weights=H_norm_sw_ins[i])

            LW_AVERAGE_det[i] = np.average(np.arange(499.5,0.5,-1), weights=H_norm_lw_det[i])
            LW_AVERAGE_ins[i] = np.average(np.arange(499.5,0.5,-1), weights=H_norm_lw_ins[i])

            LW_AVERAGE_1_det[i] = np.average(np.arange(499.5,0.5,-1), weights=H_norm_lw_1_det[i])
            LW_AVERAGE_1_ins[i] = np.average(np.arange(499.5,0.5,-1), weights=H_norm_lw_1_ins[i])
            LW_AVERAGE_2_det[i] = np.average(np.arange(499.5,0.5,-1), weights=H_norm_lw_2_det[i])
            LW_AVERAGE_2_ins[i] = np.average(np.arange(499.5,0.5,-1), weights=H_norm_lw_2_ins[i])
            LW_AVERAGE_3_det[i] = np.average(np.arange(499.5,0.5,-1), weights=H_norm_lw_3_det[i])
            LW_AVERAGE_3_ins[i] = np.average(np.arange(499.5,0.5,-1), weights=H_norm_lw_3_ins[i])
            LW_AVERAGE_4_det[i] = np.average(np.arange(499.5,0.5,-1), weights=H_norm_lw_4_det[i])
            LW_AVERAGE_4_ins[i] = np.average(np.arange(499.5,0.5,-1), weights=H_norm_lw_4_ins[i])

            SW_AVERAGE_1_det[i] = np.average(np.arange(99.5,-899.5,-1), weights=H_norm_sw_1_det[i])
            SW_AVERAGE_1_ins[i] = np.average(np.arange(99.5,-899.5,-1), weights=H_norm_sw_1_ins[i])
            SW_AVERAGE_2_det[i] = np.average(np.arange(99.5,-899.5,-1), weights=H_norm_sw_2_det[i])
            SW_AVERAGE_2_ins[i] = np.average(np.arange(99.5,-899.5,-1), weights=H_norm_sw_2_ins[i]) 
            SW_AVERAGE_3_det[i] = np.average(np.arange(99.5,-899.5,-1), weights=H_norm_sw_3_det[i])
            SW_AVERAGE_3_ins[i] = np.average(np.arange(99.5,-899.5,-1), weights=H_norm_sw_3_ins[i])
            SW_AVERAGE_4_det[i] = np.average(np.arange(99.5,-899.5,-1), weights=H_norm_sw_4_det[i])
            SW_AVERAGE_4_ins[i] = np.average(np.arange(99.5,-899.5,-1), weights=H_norm_sw_4_ins[i]) 
            
        print('averages done')

        counts = np.histogram(np.asarray(tsc_t).flatten(),bins=np.arange(0,500,1))[0]
        counts_det = np.histogram(np.asarray(tsc_d).flatten(),bins=np.arange(0,500,1))[0]
        counts_ins = np.histogram(np.asarray(tsc_i).flatten(),bins=np.arange(0,500,1))[0]

        counts_1 = np.histogram(np.asarray(tsc_ctp_1).flatten(),bins=np.arange(0,500,1))[0]
        counts_2 = np.histogram(np.asarray(tsc_ctp_2).flatten(),bins=np.arange(0,500,1))[0]
        counts_3 = np.histogram(np.asarray(tsc_ctp_3).flatten(),bins=np.arange(0,500,1))[0]
        counts_4 = np.histogram(np.asarray(tsc_ctp_4).flatten(),bins=np.arange(0,500,1))[0]

        counts_det_1 = np.histogram(np.asarray(tsc_d_1).flatten(),bins=np.arange(0,500,1))[0]
        counts_ins_1 = np.histogram(np.asarray(tsc_i_1).flatten(),bins=np.arange(0,500,1))[0]
        counts_det_2 = np.histogram(np.asarray(tsc_d_2).flatten(),bins=np.arange(0,500,1))[0]
        counts_ins_2 = np.histogram(np.asarray(tsc_i_2).flatten(),bins=np.arange(0,500,1))[0]
        counts_det_3 = np.histogram(np.asarray(tsc_d_3).flatten(),bins=np.arange(0,500,1))[0]
        counts_ins_3 = np.histogram(np.asarray(tsc_i_3).flatten(),bins=np.arange(0,500,1))[0]
        counts_det_4 = np.histogram(np.asarray(tsc_d_4).flatten(),bins=np.arange(0,500,1))[0]
        counts_ins_4 = np.histogram(np.asarray(tsc_i_4).flatten(),bins=np.arange(0,500,1))[0]

        weightedlw = (LW_AVERAGE_det * (counts_det / counts)) + (LW_AVERAGE_ins * (counts_ins/counts))
        weightedsw = (SW_AVERAGE_det * (counts_det / counts)) + (SW_AVERAGE_ins * (counts_ins/counts))

        weightedlw_1 = (LW_AVERAGE_1_det * (counts_det_1 / counts_1)) + (LW_AVERAGE_1_ins * (counts_ins_1/counts_1))
        weightedsw_1 = (SW_AVERAGE_1_det * (counts_det_1 / counts_1)) + (SW_AVERAGE_1_ins * (counts_ins_1/counts_1))
        weightedlw_2 = (LW_AVERAGE_2_det * (counts_det_2 / counts_2)) + (LW_AVERAGE_2_ins * (counts_ins_2/counts_2))
        weightedsw_2 = (SW_AVERAGE_2_det * (counts_det_2 / counts_2)) + (SW_AVERAGE_2_ins * (counts_ins_2/counts_2))
        weightedlw_3 = (LW_AVERAGE_3_det * (counts_det_3 / counts_3)) + (LW_AVERAGE_3_ins * (counts_ins_3/counts_3))
        weightedsw_3 = (SW_AVERAGE_3_det * (counts_det_3 / counts_3)) + (SW_AVERAGE_3_ins * (counts_ins_3/counts_3))
        weightedlw_4 = (LW_AVERAGE_4_det * (counts_det_4 / counts_4)) + (LW_AVERAGE_4_ins * (counts_ins_4/counts_4))
        weightedsw_4 = (SW_AVERAGE_4_det * (counts_det_4 / counts_4)) + (SW_AVERAGE_4_ins * (counts_ins_4/counts_4))

        weightedlw = np.where(np.isnan(weightedlw),LW_AVERAGE_all,weightedlw)
        weightedsw = np.where(np.isnan(weightedsw),SW_AVERAGE_all,weightedsw)

        weightedlw_1 = np.where(np.isnan(weightedlw_1),LW_AVERAGE_1,weightedlw_1)
        weightedsw_1 = np.where(np.isnan(weightedsw_1),SW_AVERAGE_1,weightedsw_1)
        weightedlw_2 = np.where(np.isnan(weightedlw_2),LW_AVERAGE_2,weightedlw_2)
        weightedsw_2 = np.where(np.isnan(weightedsw_2),SW_AVERAGE_2,weightedsw_2)
        weightedlw_3 = np.where(np.isnan(weightedlw_3),LW_AVERAGE_3,weightedlw_3)
        weightedsw_3 = np.where(np.isnan(weightedsw_3),SW_AVERAGE_3,weightedsw_3)
        weightedlw_4 = np.where(np.isnan(weightedlw_4),LW_AVERAGE_4,weightedlw_4)
        weightedsw_4 = np.where(np.isnan(weightedsw_4),SW_AVERAGE_4,weightedsw_4)

        lwsum = np.nansum(weightedlw[:]*(counts[:]))/np.nansum(counts[:])
        swsum = np.nansum(weightedsw[:]*(counts[:]))/np.nansum(counts[:])

        LWDETsum_1 = np.nansum(LW_AVERAGE_1_det[:]*(counts_det_1[:]))/np.nansum(counts_det_1[:])
        SWDETsum_1 = np.nansum(SW_AVERAGE_1_det[:]*(counts_det_1[:]))/np.nansum(counts_det_1[:])
        LWDETsum_2 = np.nansum(LW_AVERAGE_2_det[:]*(counts_det_2[:]))/np.nansum(counts_det_2[:])
        SWDETsum_2 = np.nansum(SW_AVERAGE_2_det[:]*(counts_det_2[:]))/np.nansum(counts_det_2[:])
        LWDETsum_3 = np.nansum(LW_AVERAGE_3_det[:]*(counts_det_3[:]))/np.nansum(counts_det_3[:])
        SWDETsum_3 = np.nansum(SW_AVERAGE_3_det[:]*(counts_det_3[:]))/np.nansum(counts_det_3[:])
        LWDETsum_4 = np.nansum(LW_AVERAGE_4_det[:]*(counts_det_4[:]))/np.nansum(counts_det_4[:])
        SWDETsum_4 = np.nansum(SW_AVERAGE_4_det[:]*(counts_det_4[:]))/np.nansum(counts_det_4[:])

        forcings = np.zeros((4))
        forcings[0] = LWDETsum_1+SWDETsum_1
        forcings[1] = LWDETsum_2+SWDETsum_2
        forcings[2] = LWDETsum_3+SWDETsum_3
        forcings[3] = LWDETsum_4+SWDETsum_4
        print(str(year)+ ' forcings are: '+str(forcings))
        np.save('/disk1/Users/gah20/forcings_ocean_'+str(year),forcings)
