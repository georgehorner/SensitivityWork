import numpy as np
import netCDF4 as nc
import xarray as xr
import kdtree as kd
import pykdtree
import sklearn.neighbors as sk
import scipy.spatial as sp
import glob

def clusteredmap(time,regime,property):

    newmap = np.ma.masked_all((len(lats),len(lons)))

    for lon in range(len(lons)):
        for lat in range(len(lats)):
            pixel = (alb[time, lat, lon], ctp[time, lat, lon], cf[time, lat, lon])
            group = kd.kdtree_closest_point(kdtree,pixel)

            if (group == regime):
                newmap[lat,lon] = property[time,lat,lon]

    return newmap

def regimes(prop):
    output_array = np.ma.masked_all((8,len(time),len(lat),len(lon)),dtype=np.float32)
    for i in range(7):
        prod = (prop.flatten() * (np.ma.masked_not_equal(nnint, i) + 1)) / (i+1)
        output_array[i] = prod.reshape(cf.shape)
    return output_array

for year in range(2000,2017):
    for month in range(1,13):

        #import monthly isccp files
        isccp = xr.open_mfdataset('/disk1/Data/ISCCP/access/isccp-basic/hgg/'+str(year)+str(month).zfill(2)+'/*.nc')
        print('loaded in isccp file: ' + str(year) + str(month).zfill(2))
        #normalise and QC isccp files
        ctp = np.where(isccp['pc']<0,np.nan,isccp['pc'])
        cf = np.where(isccp['cldamt']>100,np.nan,isccp['cldamt'])
        tau = np.where(isccp['tau']<0,np.nan,isccp['tau'])
        albedo = tau ** 0.895 / ((tau**0.895) + 6.82)
        time = isccp['time']
        lat = isccp['lat']
        lon = isccp['lon']

        ctp_norm = ctp/1025 
        cf_norm = cf/100
        alb_norm = albedo/1

        #ISCCP-H WS CENTROIDS

        #r = [albedo,ctp,cf]
        r1 = (0.546,242.6/1100,99.5/100) # Deep Convective - RFO-9.56%
        r2 = (0.543,433.6/1100,99.2/100) 
        r3 = (0.147,316.3/1100,79.9/100)
        r4 = (0.229,395.6/1100,84.5/100)
        r5 = (0.524,606.9/1100,97.2/100)
        r6 = (0.293,645.1/1100,40.0/100)
        r7 = (0.336,840.1/1100,79.6/100)
        r8 = (0.432,725.5/1100,90.7/100)
        datapoints = [r1,r2,r3,r4,r5,r6,r7,r8]

        skdtree = sp.cKDTree(datapoints)
        print('skdtree done')

        data = np.transpose([alb_norm.flatten(),ctp_norm.flatten(),cf_norm.flatten()])
        nndist, nnint  = skdtree.query(data,k=1,n_jobs=10)
        cot_r = regimes(alb_norm)
        print('regimes created')

        dccraw = np.where(isccp['tc']<220,cot_r[0].filled(np.nan),np.nan)
        dcc_final = np.where(dccraw<0.5,np.nan,dccraw)

        ft = nc.Dataset('/disk1/Users/gah20/DCC/'+str(year)+'/'+str(month).zfill(2)+'.nc','w',format='NETCDF4')
        ft.createDimension('lat', 180)
        ft.createDimension('lon', 360)
        ft.createDimension('time', None)

        latitude = ft.createVariable('Latitude', 'f4', 'lat')
        longitude = ft.createVariable('Longitude', 'f4', 'lon')
        DCC = ft.createVariable('DCC', 'f4', ('time', 'lat', 'lon',))
        timen = ft.createVariable('Time', 'i4', 'time')

        longitude[:] = lon
        latitude[:] = lat
        DCC[:] = dcc_final
        timen[:] = time
        ft.close()

        print('file' + str(year) + str(month).zfill(2) + ' done')