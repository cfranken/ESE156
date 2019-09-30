#!/usr/bin/env python

import os
import sys
from optparse import OptionParser
import numpy as np
import h5py
import time
from netCDF4 import Dataset
from dateutil import rrule
from datetime import *
import time
import glob
from numba import jit
from collections import OrderedDict

dict_names = {
    'source': 'TROPOMI L1B',
    'references': '...',
    'Conventions': 'CF-1.6',
    'product_version': 'v1.0',
    'summmary': 'Fraunhofer-line based SIF retrievals',
    'keyword': 'satellite, OCO-2, Solar Induced Fluorescence, SIF',
    'keywords_vocabulary': 'NASA Global Change Master Directory (GCMD)',
    'cdm_data_type': 'grids',
    'comment': 'These data were produced at JPL/Caltech',
    'date_created': 'Created ' + time.ctime(time.time()),
    'creator_name': 'Caltech, Christian Frankenberg',
    'creator_email': 'cfranken@caltech.edu',
    'project': 'OCO-2 NASA/JPL',
    'geospatial_lat_min': '-90.0f; // float',
    'geospatial_lat_max': '90.0f; // float',
    'geospatial_lat_units': 'degrees_north',
    'geospatial_lon_min': '-180.0f; // float',
    'geospatial_lon_max': '180.0f; // float',
    'geospatial_lon_units': 'degrees_east',
    'geospatial_vertical_min': '0.0f; // float',
    'geospatial_vertical_max': '100000.0; // float',
    'standard_name_vocabulary': 'NetCDF Climate and Forecast (CF) Metadata Conventions Version 1.6',
    'platform': 'OCO-2',
    'sensor': 'OCO-2',
    'spatial_resolution': '2km x 1.3km at nadir (along-track x across-track)',
    '                   CoordSysBuilder': 'ucar.nc2.dataset.conv.CF1Convention',
    }

dict_l2 = OrderedDict({
   #'variable_name': 'standard_name_from_netcdf',
    'co_column': 'target_product/co_column',
    'co_column_precision': 'target_product/co_column_precision',
    'sza': 'instrument/solar_zenith_angle',
#    'vza': 'instrument/viewing_zenith_angle',
    'vaa': 'instrument/viewing_azimuth_angle',
#    'altitude': 'meteo/altitude_levels',
    'surface_height': 'meteo/surface_altitude',
'surface_height_sigma': 'meteo/surface_altitude_stdv',
#    'p': 'meteo/pressure_levels',
    'Tskin': 'meteo/skin_temperature',
#    'temperature': 'meteo/temperature',
    })
print(dict_l2)


def grid(options, args):
    # Get some indices (could change)
#    i_DC = list(dict_l2.keys()).index("daily_corr") 
    i_xco = list(dict_l2.keys()).index("co_column")
    i_co_sigma = list(dict_l2.keys()).index("co_column_precision") 
    # Define spatial grid:
    dlat = np.arange(options.latMin, options.latMax+1e-8, options.dLat)
    dlon = np.arange(options.lonMin, options.lonMax+1e-8, options.dLon)
    lat = np.arange(options.latMin+options.dLat/2., options.latMax+1e-8-options.dLat/2., options.dLat)
    lon = np.arange(options.lonMin+options.dLon/2., options.lonMax+1e-8-options.dLat/2., options.dLon)
    
    # how many time-slices for now
    start = datetime.strptime(options.start, "%Y-%m-%d").date()
    stop = datetime.strptime(options.stop, "%Y-%m-%d").date()
    print(start,stop)
    nT = 0
    # Check out times to be interated
    for dt in rrule.rrule(rrule.DAILY, interval=options.dTime, dtstart=start, until=stop):
      #print dt
      nT+=1

    print('Time dimension ' , nT)
    
    # Generate simple numpy arrays for averaging
    mat_data=np.zeros((1+len(dict_l2),len(lat),len(lon)))

    # Prestore dataset here (netCDF write can be slow!)
    #list_data = np.zeros((nT, len(lat),len(lon),1+len(dict_l2)))

    # create netCDF4 file:
    f = Dataset(options.outFile, 'w', format='NETCDF4')
    time = f.createDimension('time', None)
    lati = f.createDimension('lat', len(lat))
    loni = f.createDimension('lon', len(lon))
    times = f.createVariable('time','f8',('time',))
    latitudes = f.createVariable('lat','f8',('lat'))
    longitudes = f.createVariable('lon','f8',('lon'))
    latitudes[:] =  lat
    longitudes[:] = lon
    # Add units, long names, etc.
    times.units = 'days since 2014-1-1 0:0:0'
    latitudes.units="degrees_north"
    longitudes.units="degrees_east"
    latitudes.standard_name = "latitude"
    longitudes.standard_name="longitude"
    latitudes.axis = "Y"
    longitudes.axis="X"
    times.long_name = "time"
    latitudes.long_name="latitude"
    longitudes.long_name="longitude"

    # Create all NetCDF variables and store in list:
    list_nc = []

###NN: looping through the keys of the in dict_l2
### fills up the list with the variables 
    for i in dict_l2:
        list_nc.append(f.createVariable(i,"f4",("time","lat","lon",),zlib=True, least_significant_digit=4,fill_value=-999.9,chunksizes=(1, 100, 100)))
    
    # Create some extra ones here:     
    n = f.createVariable('N','f4',('time','lat','lon',), zlib=True,least_significant_digit=10,fill_value=-999.9, chunksizes=(1, 100, 100))
    XCO = f.createVariable('SIF_avg','f4',('time','lat','lon',), zlib=True,least_significant_digit=4,fill_value=-999.9, chunksizes=(1, 100, 100))
    SIF_norm = f.createVariable('SIF_avg_daily','f4',('time','lat','lon',), zlib=True,least_significant_digit=4,fill_value=-999.9, chunksizes=(1, 100, 100))
    
    # loop over files:
    counter = 0
    counter_time = 0

    # Start looping over files:
    counter = -1
    points = np.zeros((options.subGrid,options.subGrid,2))
    files = []
    for dt in (rrule.rrule(rrule.DAILY, interval=options.dTime,  dtstart=start, until=stop)):
#    for orbit_number in range(630, 670):
        # Set to 0 again
        mat_data[:]=0

#        files = glob.glob(options.folder + dt.strftime('%Y/%m/')+ 'tropomi_*_??'+dt.strftime('%m%d')+'_*.nc4')
        files = glob.glob('/export/data2/projects/TROPOMI_GASES/co/7_7/s5p_l2_co_0007_006*.nc')
#        print(filename) 
        #date_center = dt + timedelta(days=np.floor(options.dTime/2.))
        date_start = datetime(2014,1,1)
        delta = dt - date_start
        
        for file in files:
            print(file)
            fin = h5py.File(file,'r')

            lat_in = fin['instrument/latitude_center'][:]
            lon_in = fin['instrument/longitude_center'][:]
            lat_in_ = fin['instrument/latitude_corners'][:]
            lon_in_ = fin['instrument/longitude_corners'][:]
# NN May not need this for TROPOMI
            #To kick out soundings where a corner is missing (happens once in a while)
            test = np.min(lat_in_,1)
            test2 = np.min(lon_in_,1)

#            biome = fin['IGBP_index'][:]
#            mode = fin['measurement_mode'][:]
            eps =0.01
            ind2 = (lat_in>options.latMin+eps)&(lat_in<options.latMax-eps)&(lon_in<options.lonMax-eps)&(lon_in>options.lonMin+eps) 
# NN checking if there are valid points in the box 
            ind = np.array(np.where(ind2)[0])
            print(ind)
            print(file, counter , 'ind is equal to', len(ind),)
            if len(ind)>10:
                

                counter +=1
                times[counter] = delta.days
                # print len(ind)
                iLat = np.asarray(np.floor(((np.asarray(lat_in[ind])-options.latMin)/(options.latMax-options.latMin)*len(lat))),dtype=int)
                iLon = np.asarray(np.floor(((np.asarray(lon_in[ind])-options.lonMin)/(options.lonMax-options.lonMin)*len(lon))),dtype=int)
                mat_in = np.zeros((len(lat_in),len(dict_l2)+1))
                iT = counter
                #index_vector = np.asarray((iLon*len(lat)+iLat), dtype=int);
                co = 0
                # Read all data into matrix:
                for i in dict_l2:
                    print(dict_l2[i]) 
                    mat_in[:,co]=fin[dict_l2[i]][:].flatten()
 
                    co =co+1
                print('READ')
                
                mat_in[:,-1]=1
                if options.subGrid>1: #NN make use of verticies?
                    favg_all(mat_data, lat_in_[ind,:],lon_in_[ind,:],mat_in[ind,:],len(ind),mat_in.shape[1],options.subGrid, options.latMin, options.latMax, options.lonMin,options.lonMax, len(lat), len(lon), points ) 
                else:
                    favg(mat_data, iLat,iLon,ind,mat_in,len(ind),mat_in.shape[1]) 
                #favg(vec_SIF_avg,iLat,iLon,ind,0.5*sif_757nm_in+1.5/2.*sif_771nm_in,len(ind))
	#favg(vec_SIF_avg,iLat,iLon,ind,0.5*sif_757nm_in+1.5/2.*sif_771nm_in,len(ind))

                
                print('.. averaged')

            fin.close()
        # Demand a minimum of 5 points per grid cell
        # Norm Mat data
        if len(ind)>10:
            aveMat(mat_data,len(dict_l2))
            #for i in range(len(dict_l2)):    
            #    mat_data[i,:,:] = mat_data[i,:,:]/mat_data[-1,:,:]
            for i in range(len(dict_l2)):
                list_nc[i][counter,:,:]=-999.9
            # Find non-empty longitudes:
            ix = np.where(np.sum(mat_data[-1,:,:],axis=1)>2)[0]
            print(len(ix))
            #rint(len(np.where(np.sum(mat_data[-1,:,:],axis=0)>1)[0]))
            for x in range(mat_data.shape[2]):
                wo = np.where(mat_data[-1,:,x]>1)[0]
                #print(x,wo)
                if len(wo)>0:
                    for i in range(len(dict_l2)):
                        list_nc[i][counter,wo,x]=mat_data[i,wo,x]
                    n[counter,wo,x] = mat_data[-1,wo,x]/options.subGrid**2
#


    f.close()
      #  print ' done'

@jit()
def aveMat(mat_data,leng):
    for i in range(leng):
        for x in range(mat_data.shape[1]):
            wo = np.where(mat_data[-1,x,:]>1)[0]
            mat_data[i,x,wo] = mat_data[i,x,wo]/mat_data[-1,x,wo]


def getIndices(lat_in, lon_in, latMin, latMax, lonMin, lonMax, nLat, nLon):
    iLat = np.asarray(np.floor(((np.asarray(lat_in)-latMin)/(latMax-latMin)*nLat)),dtype=int)
    iLon = np.asarray(np.floor(((np.asarray(lon_in)-lonMin)/(lonMax-lonMin)*nLon)),dtype=int)
    return iLat, iLon

@jit(nopython=True, parallel=True)
def favg(arr,ix,iy,iz,inp,s,s2):
    for i in range(s):
        for z in range(s2):
            arr[z,ix[i],iy[i]]+=inp[iz[i],z]
    return arr


@jit(nopython=True, parallel=True)
def favg_all(arr,lat,lon,inp,s,s2,n, latMin, latMax, lonMin, lonMax, nLat, nLon, points):
    #w = 1/n**2
    for i in range(s):
        po = getPoints(points,lat[i,:],lon[i,:],n)
        lons = po[:,:,1].flatten()
        lats = po[:,:,0].flatten()
        ix = ((lats-latMin)/(latMax-latMin)*nLat)
        iy = ((lons-lonMin)/(lonMax-lonMin)*nLon)
        #print(ix,iy)
        #print(i,len(ix))
        for j in range(len(ix)):
            for z in range(s2):
                arr[z,int(ix[j]),int(iy[j])]+=inp[i,z]
    return arr

@jit(nopython=True)
def divLine(lat1,lon1,lat2,lon2,n):
    dLat = (lat2-lat1)/(2*n)
    dLon = (lon2-lon1)/(2*n)
    lats = np.linspace(lat1+dLat,lat2-dLat,n)
    lons = np.linspace(lon1+dLon,lon2-dLon,n)
    return lats,lons
## NN look here
# Divide each polygon into multiple points
@jit(nopython=True, parallel=True)
def getPoints(points, vert_lat, vert_lon, n):
# NN: Creates lines for the left and right edges
    # Get reference points for two lines at the extremes:
    lats_0, lons_0 = divLine(vert_lat[0],vert_lon[0],vert_lat[1],vert_lon[1],n)
    lats_1, lons_1 = divLine(vert_lat[3],vert_lon[3],vert_lat[2],vert_lon[2],n)
    for i in range(n):
        points[i,:,0], points[i,:,1] = divLine(lats_0[i], lons_0[i] ,lats_1[i], lons_1[i],n)
    return points

def standalone_main():
    parser = OptionParser(usage="usage: %prog l2_file2")
    parser.add_option( "-o","--outFile", dest="outFile",
                       default='co_tropomi_l3.nc',
                       help="output filename (default OCO2_SIF_map.nc)")
    parser.add_option( "--latMin", dest="latMin",
                       type=float,
                       default=-90,
                       help="min latitude region")
    parser.add_option( "--dLat", dest="dLat",
                       type=float,
                       default=1,
                       help="latitude resolution (1 degree default)")
    parser.add_option( "--dLon", dest="dLon",
                       type=float,
                       default=1,
                       help="longitude resolution (1 degree default)")
    parser.add_option( "--startTime", dest="start",
                       default='2017-11-08',
                       help="default 2014-09-06")
    parser.add_option( "--stopTime", dest="stop",
                       default='2017-12-08',
                       help="default 2015-01-01")
    parser.add_option( "--subGrid", dest="subGrid",
                       type=int,
                       default=1,
                       help="How many subgrid points to use for gridding (default=1 (just center pixel)")
    parser.add_option( "--dTime", dest="dTime",
                       default=1,
                       type=int,
                       help="default 1 day (determines time window size for each time step)")
    parser.add_option('--mode', dest='mode', type=int, default=-1,
                      help='mode (0=ND, 1=GL, 2=TG, etc, -1 for all)')

    parser.add_option('--biome', dest='biome', type=int, default=-1,
                      help='IGBP biome type (-1 default for ALL)')

    parser.add_option( "--latMax", dest="latMax",
                       type=float,
                       default=90,
                       help="max latitude region")
    parser.add_option( "--lonMin", dest="lonMin",
                       type=float,
                       default=-180,
                       help="min longitude region")
    parser.add_option( "--lonMax", dest="lonMax",
                       type=float,
                       default=180,
                       help="max longitude region")
    parser.add_option( "--folder", dest="folder",
                       default='/export/data2/projects/TROPOMI_GASES/co/7_7',
                       help="Default folder DIR root where data is stored in DIR/YYYY/MM/DD/LtSIF/oco2*.nc4")
   # /data/oco2/scf/product/B7???r/r0?/'
    # Parse command line arguments
    (options, args) = parser.parse_args()
    
    # start gridding
    grid(options, args)

if __name__ == "__main__":
    standalone_main()



