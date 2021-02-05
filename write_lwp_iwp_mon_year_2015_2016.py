import xarray as xr
import os
import sys
import numpy as np
from functions import global_mean, yearly_sum

if __name__ == '__main__':
    P = os.path.join
    dt_dir = './output_data/lwp_iwp_year_dt_flip_lon_t42_2015_2016'
    saved_dir = './output_data'
    if not os.path.exists(saved_dir):
        os.makedirs(saved_dir)

    for year in range(2015, 2017): 
        print(year)

        # Read year data
        ds = xr.open_dataset(P(dt_dir, 'visit_num_'+str(year)+'.nc'), decode_times=False)
        visit_num = yearly_sum(ds.visit_num)

        nlats = len(ds.lat)
        nlons = len(ds.lon)

        lwp = np.zeros((1, nlats, nlons))
        iwp = np.zeros((1, nlats, nlons))
        cwp = np.zeros((1, nlats, nlons))

        ds = xr.open_dataset(P(dt_dir, 'lwp_sum_'+str(year)+'.nc'), decode_times=False)
        lwp_sum = yearly_sum(ds.lwp_sum)
        lwp[0,:,:] = np.where(visit_num > 0, lwp_sum / visit_num, np.nan) 

        lwp = xr.DataArray(lwp, dims=('year','lat', 'lon'), coords=(np.array([year]), ds.lat, ds.lon))
        #lwp_gm = global_mean(lwp[0,:,:])
        #lwp_gm = xr.DataArray(np.array([lwp_gm]), dims=('year'), coords={'year':np.array(year)})

        ds = xr.open_dataset(P(dt_dir, 'iwp_sum_'+str(year)+'.nc'), decode_times=False)
        iwp_sum = yearly_sum(ds.iwp_sum)
        iwp[0,:,:] = np.where(visit_num > 0, iwp_sum / visit_num, np.nan) 
        iwp = xr.DataArray(iwp, coords=(np.array([year]), ds.lat, ds.lon), dims=('year', 'lat', 'lon'))
        #iwp_gm = global_mean(iwp[0,:,:])
        #iwp_gm = xr.DataArray(np.array([iwp_gm]), dims=('year'), coords={'year':np.array(year)})
       
        cwp[0,:,:] = lwp[0,:,:] + iwp[0,:,:]
        cwp = xr.DataArray(cwp, dims=('year','lat', 'lon'), coords=(np.array([year]), ds.lat, ds.lon))
        #cwp_gm = lwp_gm.values + iwp_gm.values
        #cwp_gm = xr.DataArray(cwp_gm, dims=('year'), coords={'year':np.array(year)})

        ds_out = xr.Dataset({'lwp':lwp, 'iwp':iwp, 'cwp':cwp}, coords={'year':np.array(year), 'lat':ds.lat, 'lon':ds.lon})

        ds_out.to_netcdf(os.path.join(saved_dir, 'cld_water_path_data_t42_'+str(year)+'.nc'), format='NETCDF4_CLASSIC', mode='w')
       
        '''
        lwp_mon = np.zeros((12, nlats, nlons))
        iwp_mon = np.zeros((12, nlats, nlons))
        cwp_mon = np.zeros((12, nlats, nlons))
       
        lwp_gm_mon = np.zeros(12)
        iwp_gm_mon = np.zeros(12)
        cwp_gm_mon = np.zeros(12)

        for mm, mon in enumerate(range(12)):
            print(mon)
            mon_str = str(mon+1).zfill(2)

            ds = xr.open_dataset(P(dt_dir,'visit_num_'+str(year)+mon_str+'.nc'), decode_times=False)
            visit_num = ds.visit_num

            ds = xr.open_dataset(P(dt_dir,'lwp_sum_'+str(year)+mon_str+'.nc'), decode_times=False)
            lwp_sum = ds.lwp_sum

            lwp = np.where(visit_num > 0, lwp_sum / visit_num, np.nan) 
            lwp = xr.DataArray(lwp, coords=(ds.lat, ds.lon), dims=('lat', 'lon'))
            lwp_gm = global_mean(lwp)
            print('lwp global_mean=', lwp_gm)

            ds = xr.open_dataset(P(dt_dir,'iwp_sum_'+str(year)+mon_str+'.nc'), decode_times=False)
            iwp_sum = ds.iwp_sum
            iwp = np.where(visit_num > 0, iwp_sum / visit_num, np.nan) 
            iwp = xr.DataArray(iwp, coords=(ds.lat, ds.lon), dims=('lat', 'lon'))
            iwp_gm = global_mean(iwp)
            print('iwp global_mean=', iwp_gm)

            lwp_mon[mm,:,:] = lwp 
            iwp_mon[mm,:,:] = iwp 
            cwp_mon[mm,:,:] = lwp + iwp
            lwp_gm_mon[mm] = lwp_gm
            iwp_gm_mon[mm] = iwp_gm
            cwp_gm_mon[mm] = lwp_gm + iwp_gm

        mons = np.arange(1, 13)
        lwp_mon = xr.DataArray(lwp_mon, coords=(mons, ds.lat, ds.lon), dims=('month', 'lat', 'lon'))
        iwp_mon = xr.DataArray(iwp_mon, coords=(mons, ds.lat, ds.lon), dims=('month', 'lat', 'lon'))
        cwp_mon = xr.DataArray(cwp_mon, coords=(mons, ds.lat, ds.lon), dims=('month', 'lat', 'lon'))
        lwp_gm_mon = xr.DataArray(lwp_gm_mon, coords={'month':mons}, dims=('month'))
        iwp_gm_mon = xr.DataArray(iwp_gm_mon, coords={'month':mons}, dims=('month'))
        cwp_gm_mon = xr.DataArray(cwp_gm_mon, coords={'month':mons}, dims=('month'))
        #iwp_gm_mon = xr.DataArray(iwp_gm_mon, coords=(mons), dims=('month'))
        #cwp_gm_mon = xr.DataArray(cwp_gm_mon, coords=(mons), dims=('month'))

        coords={'month':mons, 'lat':ds.lat, 'lon':ds.lon}
        var_dict = {'lwp':lwp_mon, 'iwp':iwp_mon, 'cwp':cwp_mon,'lwp_gm':lwp_gm_mon, 'iwp_gm':iwp_gm_mon, 'cwp_gm':cwp_gm_mon}

        ds_out = xr.Dataset(var_dict, coords={'month':mons, 'lat':ds.lat, 'lon':ds.lon})
        ds_out.to_netcdf(os.path.join(saved_dir, 'cld_water_path_data_t42_monthly_'+str(year)+'.nc'), format='NETCDF4_CLASSIC', mode='w')
        '''

