from __future__ import print_function
import numpy as np
import xarray as xr
# For the hdf file
from pyhdf.HDF import *
from pyhdf.VS import *
from pyhdf.SD import *

def get_lat_lon_lwp_iwp_from_hdf(fn):
    f = HDF(fn)
    vs = f.vstart()
    #print(vs.vdatainfo())

    Latitude = vs.attach('Latitude')
    Longitude = vs.attach('Longitude')
    #Profile_time = vs.attach('Profile_time')
    #UTC_start = vs.attach('UTC_start')
    #DEM_elevation = vs.attach('DEM_elevation')

    RO_liq_water_path = vs.attach('RO_liq_water_path')
    RO_ice_water_path = vs.attach('RO_ice_water_path')

    lats = np.squeeze(Latitude[:])
    lons = np.squeeze(Longitude[:])
    lons = np.where(lons<0, lons+360, lons)
    # profile_time= np.squeeze(Profile_time[:])
    # utc_start= np.squeeze(UTC_start[:])
    # dem= np.squeeze(DEM_elevation[:])

    lwp = np.squeeze(RO_liq_water_path[:])
    iwp = np.squeeze(RO_ice_water_path[:])

    lwp[lwp<-10] = np.nan
    iwp[iwp<-10] = np.nan

    Latitude.detach()
    Longitude.detach()
    RO_liq_water_path.detach()
    RO_ice_water_path.detach()
    vs.end()
    f.close()

    return lats, lons, lwp, iwp

def global_mean(dt):
    dt_ma = np.ma.array(dt, mask=np.isnan(dt))
    lats = dt.lat
    coslat = np.cos(np.deg2rad(lats))
    dt_gm = np.ma.average(np.ma.average(dt_ma, axis=1), axis=0, weights=coslat)
    return dt_gm

def yearly_sum(dt):
    dt_ma = np.ma.array(dt, mask=np.isnan(dt))
    dt_ym = np.ma.sum(dt_ma, axis=0)
    dt_ym = xr.DataArray(dt_ym, dims=('lat', 'lon'), coords=(dt.lat, dt.lon))
    return dt_ym
