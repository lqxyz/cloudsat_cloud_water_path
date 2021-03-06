from __future__ import print_function
import xarray as xr
import os
import numpy as np
import sys
import glob
import timeit
import warnings
warnings.simplefilter(action='ignore')
from calendar import monthrange
from datetime import date
from functions import get_lat_lon_lwp_iwp_from_hdf

def gridding_data():
    for i in range(nlats):
        for j in range(nlons):
            slat = bd_lats[i]
            nlat = bd_lats[i+1]
            llon = bd_lons[j]
            rlon = bd_lons[j+1]

            ind = (lats>=slat) & (lats<=nlat) & (lons>=llon) & (lons<=rlon)
            if np.sum(ind)>0:
                lwp_grid_dt[i,j] = np.nansum([ lwp_grid_dt[i,j], np.nanmean(lwp[ind]) ])
                iwp_grid_dt[i,j] = np.nansum([ iwp_grid_dt[i,j], np.nanmean(iwp[ind]) ])
                visit_dt[i,j] = visit_dt[i,j] + 1
            else:
                lwp_grid_dt[i,j] = np.nansum([lwp_grid_dt[i,j], np.nan])
                iwp_grid_dt[i,j] = np.nansum([iwp_grid_dt[i,j], np.nan])

if __name__ == '__main__':
    P = os.path.join
    saved_dir = './output_data/lwp_iwp_year_dt_flip_lon_t42_2012_2014'
    if not os.path.exists(saved_dir):
        os.makedirs(saved_dir)
    
    ## Read t42 lat lon boundary and lat lons from dataset
    t42_fn = 'atmos_monthly.nc'
    ds_t42 = xr.open_dataset(t42_fn, decode_times=False)
    bd_lats = ds_t42.latb.values
    bd_lons = ds_t42.lonb.values
    lats1 = ds_t42.lat.values
    lons1 = ds_t42.lon.values
    nlats = len(lats1)
    nlons = len(lons1)
    
    '''
    bd_lats = np.arange(-90, 91, 1)
    bd_lons = np.arange(-180, 181, 1)
    nlats = len(bd_lats)-1
    nlons = len(bd_lons)-1
    lats1 = (bd_lats[0:-1]+bd_lats[1:]) / 2.0
    lons1 = (bd_lons[0:-1]+bd_lons[1:]) / 2.0
    '''

    basedir = 'ftp.cloudsat.cira.colostate.edu/2B-CWC-RO.P1_R05/'

    # input year from command line
    year = int(sys.argv[1]) # e.g. 2012

    start = timeit.default_timer()

    lwp_grid_dt = np.zeros((nlats, nlons))
    iwp_grid_dt = np.zeros((nlats, nlons))
    visit_dt = np.zeros((nlats, nlons))

    for d in range(366):
        day = d + 1
        day_str = str(day).zfill(3)

        start_d = timeit.default_timer()

        fns = sorted(glob.glob(P(basedir, str(year), day_str, '*.hdf')))
        if fns != []:
            for kk, fn in enumerate(fns):

                start1 = timeit.default_timer()
                print(os.path.basename(fn))

                lats, lons, lwp, iwp = get_lat_lon_lwp_iwp_from_hdf(fn)
                print('mean lwp/iwp:', np.nanmean(lwp), np.nanmean(iwp))

                gridding_data()

                stop = timeit.default_timer()
                print(kk, ' Time has passed: ', stop - start1)

            stop = timeit.default_timer()
            print('Day ', day, ' Time has passed: ', stop - start_d)
            print('Day ', day, ' Total time has passed: ', stop - start)
            print('')
        else:
            print('Day ', day, ' is empty.')

    stop = timeit.default_timer()
    print('Total time has passed: ', stop - start)

    # ====================== save data sets ======================    
    visit_dt = xr.DataArray(visit_dt, dims=('lat', 'lon'), coords=(lats1, lons1))
    lwp_grid_dt = xr.DataArray(lwp_grid_dt, dims=('lat', 'lon'), coords=(lats1, lons1))
    iwp_grid_dt = xr.DataArray(iwp_grid_dt, dims=('lat', 'lon'), coords=(lats1, lons1))

    coords = {'lat': lats1, 'lon': lons1}

    visit_ds = xr.Dataset({'visit_num':visit_dt}, coords=coords)
    visit_ds.to_netcdf(P(saved_dir, 'visit_num_'+str(year)+'.nc'), format='NETCDF4_CLASSIC', mode='w')

    lwp_ds = xr.Dataset({'lwp_sum':lwp_grid_dt}, coords=coords)
    lwp_ds.to_netcdf(P(saved_dir, 'lwp_sum_'+str(year)+'.nc'), format='NETCDF4_CLASSIC', mode='w')

    iwp_ds = xr.Dataset({'iwp_sum':iwp_grid_dt}, coords=coords)
    iwp_ds.to_netcdf(P(saved_dir, 'iwp_sum_'+str(year)+'.nc'), format='NETCDF4_CLASSIC', mode='w')

    print('dataset saved...')
