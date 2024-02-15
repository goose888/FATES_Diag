#!/usr/bin/env python
"""
This package contains function for global FATES output analysis
Currently all functions are created for the use of postprocessing output from ELM-FATES 
with wood harvest activity.
"""

from numba import njit
import numpy as np
import copy
import netCDF4 as nc4

def retrieve_hrvrate(fpath, area_4x5, denflag=0, getmask=True, hr_lev=0.0, yr_beg=1850):
    """
    Return harvest rate and natural PFT fraction from upscaled 4x5 LUH2 harvest rate data
    # fpath - path to the file contains harvest rate data
    # area_4x5 - grid area at the same resolution as fpath
    # denflag - 0, kgC total amount; 1, kgC per m2 density value
    # getmask - if output mask only contains grids with harvest rate hihger than "hr_lev"
    # hr_lev - value used as threshold for mask creation
    # yr_beg - begin year of all cases, 1850 is the default value
    """
    nclu = nc4.Dataset(fpath)
    prif_hv = nclu['HARVEST_VH1'][:]
    prin_hv = nclu['HARVEST_VH2'][:]
    secmf_hv = nclu['HARVEST_SH1'][:]
    secyf_hv = nclu['HARVEST_SH2'][:]
    secn_hv = nclu['HARVEST_SH3'][:]
    pct_nat_pft = nclu['PCT_NAT_PFT'][:]
    nclu.close()
    # Obtain N year from an arbitrary variable
    nyr = np.shape(prif_hv)[0] + (1700-yr_beg)

    # Transfer into per m2 density value
    if(denflag == 1):
        prif_hv_den = prif_hv / (area_4x5*1e6)
        prin_hv_den = prin_hv / (area_4x5*1e6)
        secmf_hv_den = secmf_hv / (area_4x5*1e6)
        secyf_hv_den = secyf_hv / (area_4x5*1e6)
        secn_hv_den = secn_hv / (area_4x5*1e6)

        # Aggregate the number to certain year
        luh2_pri_hv = np.sum(prif_hv_den[:, :, :], 0) + np.sum(prin_hv_den[:, :, :], 0)
        luh2_sec_m_hv = np.sum(secmf_hv_den[:, :, :], 0)
        luh2_sec_y_hv = np.sum(secyf_hv_den[:, :, :], 0) + np.sum(secn_hv_den[:, :, :], 0)
        luh2_tot_hv = luh2_pri_hv + luh2_pri_hv + luh2_sec_y_hv
    else:   
        # Aggregate the harvest amount
        luh2_pri_hv = prif_hv + prin_hv
        luh2_sec_m_hv = secmf_hv
        luh2_sec_y_hv = secyf_hv + secn_hv
        luh2_tot_hv = prif_hv + prin_hv + secmf_hv + secyf_hv + secn_hv

    if(getmask):
        # mask value 1.0 represents nan
        pri_mask = copy.deepcopy(luh2_pri_hv)
        pri_mask[pri_mask>hr_lev] = -1.0
        pri_mask = pri_mask + 1.0
        luh2_pri_hv.mask = pri_mask
    
        secm_mask = copy.deepcopy(luh2_sec_m_hv)
        secm_mask[secm_mask>hr_lev] = -1.0
        secm_mask = secm_mask + 1.0
        luh2_sec_m_hv.mask = secm_mask
    
        secy_mask = copy.deepcopy(luh2_sec_y_hv)
        secy_mask[secy_mask>hr_lev] = -1.0
        secy_mask = secy_mask + 1.0
        luh2_sec_y_hv.mask = secy_mask
    
        tothrv_mask = copy.deepcopy(luh2_tot_hv)
        tothrv_mask[tothrv_mask>hr_lev] = -1.0
        tothrv_mask = tothrv_mask + 1.0
        luh2_tot_hv.mask = tothrv_mask

    # Crop and bareground fraction
    pct_crop = copy.deepcopy(pct_nat_pft[(yr_beg-1700):(yr_beg-1700+nyr),15,:,:])
    pct_bg = copy.deepcopy(pct_nat_pft[(yr_beg-1700):(yr_beg-1700+nyr),0,:,:])
    pct_nat = 100.0 - pct_crop - pct_bg
    
    return luh2_pri_hv, luh2_sec_m_hv, luh2_sec_y_hv, luh2_tot_hv, pct_nat

def retrieve_secfrac_luh2(fpath, verbose_output=False):
    """
    Only retreive secondary forest area time series
    Can return the secondary forest fraction and grid area time series if set verbose_output to True
    # fpath - path to LUH2 file contains secondat forest fraction data
    # verbose_output - output contains grid level map and grid area, otherwise only timeseries of global
    # total fraction
    """
    nclu=nc4.Dataset(fpath)
    secf_luh2 = nclu['secdf'][:]
    secn_luh2 = nclu['secdn'][:]
    lat_luh2 = nclu['lat'][:]
    lon_luh2 = nclu['lon'][:]
    nclu.close()
    
    # Global grid area for 0.25 deg for LUH2 datasets
    # Since LUH2 dataset has fixed resolution, nlat and nlon are now hard-coded here.
    nlon_loc = 1440
    nlat_loc = 720
    nyr = 166
    grid_area = np.ones((nlat_loc, nlon_loc))*-9999.
    earth_area = 5.096e14
    lat = np.arange(-89.875, 90.125, 0.25)
    res = 0.25;
    for i in np.arange(0,nlat_loc):
        for j in np.arange(0,nlon_loc):
            grid_area[i,j] = (earth_area/2)*abs(np.sin((lat[i] - res/2)*np.pi/180) -
                        np.sin((lat[i] + res/2)*np.pi/180))/(360/res)
    
    sec_tot = secf_luh2+secn_luh2
    sec_tot[secf_luh2<=0.0] = np.float('nan')
    seca_luh2 = copy.deepcopy(sec_tot[(1165-nyr):1165,:,:])
    for i in np.arange(0,nyr):
        seca_luh2[i,:,:] = sec_tot[(1165-nyr+i),:,:] * grid_area
    # m2 to km2
    seca_luh2_tot = np.nansum(np.nansum(seca_luh2, 1), 1)/1e6
    
    if(verbose_output):
        return seca_luh2_tot, sec_tot, grid_area
    else:
        return seca_luh2_tot

@njit
def apply_mask_ts(ts, mts, ts_dim, fdir, thold=0.0):
    """
    Apply a mask timeseries to the timeseries with the same dimension
    # ts - [array] timeseries (list_len * year * lat * lon)
    # mts - [array] timeseries of the mask, same dimension as ts
    # ts_dim - [array, 3*1] size of ts, time * lon * lat 
    # fdir - [scalar] the direction of filter (0 - less equal; 1 - larger equal)
    """
    time_len = ts_dim[0]
    lat_len = ts_dim[1]
    lon_len = ts_dim[2]
    for iyr in np.arange(0, time_len):
        for j in np.arange(0, lat_len):
            for k in np.arange(0, lon_len):
                if(fdir == 0):
                    if(mts[iyr,j,k]<= thold):
                        ts[iyr,j,k] = np.nan
                else:
                    if(mts[iyr,j,k]>= thold):
                        ts[iyr,j,k] = np.nan
    return ts    

def obtain_logging_impact(modname, modname_hrv, varlist, sel_yr, output_ratio = False):
    """
    Obtain the logging impact either in the original unit or in the percentage
    Note: Weight shall be multiplied with the corresponding variables before calling 
    and shall be divided by the return value (if in the original unit) after 
    calling the function. 
    # modname - list contains name of all model cases
    # modname_hrv - list contains name of all model cases with wood harvest activated 
    # varlist - list contains the data of corresponding variables
    # sel_yr - if scalar, pick a specific year of data; if 1d array, the first index is
    #        the begin year and the second index is the end year
    # output_ratio - true to output relative ratio (% change), false to output absolute change
    """
    case_len = len(modname)
    var_diff = []
    var_logging_global = []
    var_logging_global_ts = []
    var_diff_deno = []
    var_diff_all = []
    var_diff_all_deno = []
    var_logging_global_deno = []
    var_logging_global_ts_deno = []
    for i in np.arange(0, case_len):
        # Contour plot
        count = len(sel_yr)
        if (count > 1):
            begidx = sel_yr[0]
            endidx = sel_yr[count - 1]
            if(modname[i][4:5] == 'c' or modname[i][4:5] == 'a'):
                if(modname[i][4:5] == 'c'):
                    if(modname[i] != 'har_c'):
                        var_diff.append(np.nanmean(varlist[i][begidx:endidx,:,:] - varlist[i+1][begidx:endidx,:,:], 0))
                        var_diff_deno.append(np.nanmean(varlist[i][begidx:endidx,:,:], 0))
                        var_diff_all.append(varlist[i][begidx:endidx,:,:] - varlist[i+1][begidx:endidx,:,:])
                        var_diff_all_deno.append(varlist[i][begidx:endidx,:,:])
                    else:
                        var_diff.append(np.nanmean(varlist[i][begidx:endidx,:,:] - varlist[1][begidx:endidx,:,:], 0))
                        var_diff_deno.append(np.nanmean(varlist[i][begidx:endidx,:,:], 0))
                        var_diff_all.append(varlist[i][begidx:endidx,:,:] - varlist[1][begidx:endidx,:,:])
                        var_diff_all_deno.append(varlist[i][begidx:endidx,:,:])
                else:
                    if(modname[i][4:5] == 'a'):
                        var_diff.append(np.nanmean(varlist[i][begidx:endidx,:,:] - varlist[i-1][begidx:endidx,:,:], 0))
                        var_diff_deno.append(np.nanmean(varlist[i][begidx:endidx,:,:], 0))
                        var_diff_all.append(varlist[i][begidx:endidx,:,:] - varlist[i-1][begidx:endidx,:,:])
                        var_diff_all_deno.append(varlist[i][begidx:endidx,:,:])
        else:
            if(modname[i][4:5] == 'c' or modname[i][4:5] == 'a'):
                if(modname[i][4:5] == 'c'):
                    if(modname[i] != 'har_c'):
                        var_diff.append(varlist[i][sel_yr[0],:,:] - varlist[i+1][sel_yr[0],:,:])
                        var_diff_deno.append(varlist[i][sel_yr[0],:,:])
                    else:
                        var_diff.append(varlist[i][sel_yr[0],:,:] - varlist[1][sel_yr[0],:,:])
                        var_diff_deno.append(varlist[i][sel_yr[0],:,:])
                else:
                    if(modname[i][4:5] == 'a'):
                        var_diff.append(varlist[i][sel_yr[0],:,:] - varlist[i-1][sel_yr[0],:,:])
                        var_diff_deno.append(varlist[i][sel_yr[0],:,:])
    for i in np.arange(0, len(modname_hrv)):
        var_logging_global.append(np.nanmean(np.nanmean(var_diff[i], 0), 0))
        var_logging_global_deno.append(np.nanmean(np.nanmean(var_diff_deno[i], 0), 0))
        if (count > 1):
            var_logging_global_ts.append(np.nanmean(np.nanmean(var_diff_all[i], 1), 1))
            var_logging_global_ts_deno.append(np.nanmean(np.nanmean(var_diff_all_deno[i], 1), 1))
    if(output_ratio):
        global_total = 100.0 * np.array(var_logging_global)/np.array(var_logging_global_deno)
        if (count>1):
            global_std = np.std(100.0 * np.array(var_logging_global_ts)/np.array(var_logging_global_ts_deno), axis=1)
        else:
            global_std = 0.0
    else:
        global_total = np.array(var_logging_global)
        if (count>1):
            global_std = np.std(np.array(var_logging_global_ts), axis=1)
        else:
            global_std = 0.0
    return [var_diff, global_total, global_std]

def obtain_logging_impact_age(modname, varlist, sel_yr, age_bin = 8, output_ratio = False):
    """
    Obtain the logging impact either in the original unit or in the percentage
    Version to account for different age distribution
    Only difference is that varlist contain the age dimension (8 bins)
    # age_bin - number of age categories
    """
    var_diff = []
    var_logging_global = []
    var_diff_deno = []
    var_logging_global_deno = []
    for i in np.arange(0, case_len):
        for j in np.arange(0, age_bin):
            count = len(sel_yr)
            if (count > 1):
                begidx = sel_yr[0]
                endidx = sel_yr[count - 1]
                if(modname[i][4:5] == 'c' or modname[i][4:5] == 'a'):
                    if(modname[i][4:5] == 'c'):
                        if(modname[i] != 'har_c'):
                            var_diff.append(np.nanmean(varlist[i][begidx:endidx,j,:,:] - varlist[i+1][begidx:endidx,j,:,:], 0))
                            var_diff_deno.append(np.nanmean(varlist[i][begidx:endidx,j,:,:], 0))
                        else:
                            var_diff.append(np.nanmean(varlist[i][begidx:endidx,j,:,:] - varlist[1][begidx:endidx,j,:,:], 0))
                            var_diff_deno.append(np.nanmean(varlist[i][begidx:endidx,j,:,:], 0))
                    else:
                        if(modname[i][4:5] == 'a'):
                            var_diff.append(np.nanmean(varlist[i][begidx:endidx,j,:,:] - varlist[i-1][begidx:endidx,j,:,:], 0))
                            var_diff_deno.append(np.nanmean(varlist[i][begidx:endidx,j,:,:], 0))
            else:
                if(modname[i][4:5] == 'c' or modname[i][4:5] == 'a'):
                    if(modname[i][4:5] == 'c'):
                        if(modname[i] != 'har_c'):
                            var_diff.append(varlist[i][sel_yr[0],j,:,:] - varlist[i+1][sel_yr[0],j,:,:])
                            var_diff_deno.append(varlist[i][sel_yr[0],j,:,:])
                        else:
                            var_diff.append(varlist[i][sel_yr[0],j,:,:] - varlist[1][sel_yr[0],j,:,:])
                            var_diff_deno.append(varlist[i][sel_yr[0],j,:,:])
                    else:
                        if(modname[i][4:5] == 'a'):
                            var_diff.append(varlist[i][sel_yr[0],j,:,:] - varlist[i-1][sel_yr[0],j,:,:])
                            var_diff_deno.append(varlist[i][sel_yr[0],j,:,:])
    for i in np.arange(0, len(modname_hrv)):
        for j in np.arange(0, age_bin):
            var_logging_global.append(np.nanmean(np.nanmean(var_diff[i*age_bin+j], 0), 0))
            var_logging_global_deno.append(np.nanmean(np.nanmean(var_diff_deno[i*age_bin+j], 0), 0))
    if(output_ratio):
        global_total = 100.0 * np.array(var_logging_global)/np.array(var_logging_global_deno)
    else:
        global_total = np.array(var_logging_global)
    return [var_diff, global_total]

def obtain_logging_elm(modname, varlist, sel_yr, output_ratio = False):
    """
    Variant of obtain_logging_impact for ELM-only model
    All parameters are the same as the original version
    # modname - list contains name of all model cases
    # varlist - list contains the data of corresponding variables
    # sel_yr - if scalar, pick a specific year of data; if 1d array, the first index is
    #        the begin year and the second index is the end year
    # output_ratio - true to output relative ratio (% change), false to output absolute change
    """
    var_diff = []
    var_diff_ts = []
    var_logging_global = []
    var_logging_global_ts = []
    var_diff_deno = []
    var_diff_ts_deno = []
    var_logging_global_deno = []
    var_logging_global_ts_deno = []
    # Contour plot
    count = len(sel_yr)
    if (count > 1):
        begidx = sel_yr[0]
        endidx = sel_yr[count - 1]
        var_diff.append(np.nanmean(varlist[1][begidx:endidx,:,:] - varlist[0][begidx:endidx,:,:], 0))
        var_diff_deno.append(np.nanmean(varlist[0][begidx:endidx,:,:], 0))
        var_diff_ts.append(varlist[1][begidx:endidx,:,:] - varlist[0][begidx:endidx,:,:])
        var_diff_ts_deno.append(varlist[0][begidx:endidx,:,:])
    else:
        var_diff.append(varlist[1][sel_yr[0],:,:] - varlist[0][sel_yr[0],:,:])
        var_diff_deno.append(varlist[0][sel_yr[0],:,:])
    var_logging_global.append(np.nanmean(np.nanmean(var_diff[0], 0), 0))
    var_logging_global_deno.append(np.nanmean(np.nanmean(var_diff_deno[0], 0), 0))
    if (count > 1):
        var_logging_global_ts.append(np.nanmean(np.nanmean(var_diff_ts[0], 1), 1))
        var_logging_global_ts_deno.append(np.nanmean(np.nanmean(var_diff_ts_deno[0], 1), 1))
    if(output_ratio):
        global_total = 100.0 * np.array(var_logging_global)/np.array(var_logging_global_deno)
        if (count > 1):
            global_std = np.std(100.0 * np.array(var_logging_global_ts)/np.array(var_logging_global_ts_deno))
        else:
            global_std = 0.0
    else:
        global_total = np.array(var_logging_global)
        global_std = np.std(np.array(var_logging_global_ts))
    return [var_diff, global_total, global_std]

def obtain_logging_impact_ts(modname, modname_hrv, varlist, output_ratio = False):
    """
    Obtain the whole time series of the logging impact
    Use the January output from with and without case to perform diff
    # modname - list contains name of all model cases
    # varlist - list contains the data of corresponding variables
    # area_4x5 - grid area at the same resolution as input data
    # weighted - if true, the area will be applied as weight; otherwise we will directly use raw value
    # output_ratio - true to output relative ratio (% change), false to output absolute change
    """
    var_logging = []
    var_logging_global = []
    var_logging_deno = []
    var_logging_global_deno = []
    case_len = len(modname)
    idx = 0
    nyr = np.shape(varlist[0])[0]/12
    for i in np.arange(0, case_len):
        var_logging_deno.append(varlist[i][idx:12*nyr+idx:12,:,:])
        if(modname[i][4:5] == 'c'):
            if(modname[i] != 'har_c'):
                var_logging.append(varlist[i][idx:12*nyr+idx:12,:,:] - varlist[i+1][idx:12*nyr+idx:12,:,:])
            else:
                var_logging.append(varlist[i][idx:12*nyr+idx:12,:,:] - varlist[1][idx:12*nyr+idx:12,:,:])
        else:
            if(modname[i][4:5] == 'a'):
                var_logging.append(varlist[i][idx:12*nyr+idx:12,:,:] - varlist[i-1][idx:12*nyr+idx:12,:,:])
    for i in np.arange(0, len(modname_hrv)):
        var_logging[i][var_logging[i]<-2e2] = np.nan
        var_logging[i][var_logging[i]>2e2] = np.nan
        var_logging_deno[i][var_logging_deno[i]<-2e2] = np.nan
        var_logging_deno[i][var_logging_deno[i]>2e2] = np.nan
        var_logging_global.append(np.nanmean(np.nanmean(var_logging[i][:,:,:], 1), 1))
        var_logging_global_deno.append(np.nanmean(np.nanmean(var_logging_deno[i][:,:,:], 1), 1))
    if(output_ratio):
        var_logging_r = 100.0 * np.array(var_logging_global)/np.array(var_logging_global_deno)
    else:
        var_logging_r = 0.0
    return [var_logging, var_logging_global, var_logging_r]

def obtain_dir_logging_impact_ts(modname, modname_hrv, varlist, begyr, endyr, area_4x5, weighted = True, output_ratio = False):
    """
    Obtain the direct logging impact, direct impact is calculated by subtracting the variable values 
    after logging (Feb) from before the logging (Jan)
    # Excluded 1850 since it is the begin year
    # modname - list contains name of all model cases
    # varlist - list contains the data of corresponding variables
    # begyr - the begin year of the data used
    # endyr - the end year of the data used
    # area_4x5 - grid area at the same resolution as fpath
    # weighted - if true, the area will be applied as weight; otherwise we will directly use raw value
    # output_ratio - true to output relative ratio (% change), false to output absolute change
    """
    var_logging_loss = []
    var_logging_deno = []
    # Impact due to climate variation, unfortunately we cannot obtain this directly from model output
#     var_logging_var = []
    # Global aggregation
    var_logging_global = []
    var_logging_global_deno = []
    begidx_aft = 12*begyr+1
    endidx_aft = 12*endyr+1
    begidx_bef = 12*begyr
    endidx_bef = 12*endyr
    case_len = len(modname)
    # Weighted average require the pre-defined area_4x5.
    for i in np.arange(0,case_len):
        if(modname[i][4:5] == 'c' or modname[i][4:5] == 'a'):
            # Direct loss (Wt_Feb - Wt_Jan) and denominator
            if(weighted):
                temp_array = varlist[i][begidx_aft:endidx_aft:12,:,:] - varlist[i][begidx_bef:endidx_bef:12,:,:]
#                 temp_array[temp_array<-2e-1] = np.nan
#                 temp_array[temp_array>2e-1] = np.nan
                var_logging_loss.append(area_4x5[None,:,:] * temp_array[:,:,:])
                var_logging_deno.append(area_4x5[None,:,:] * varlist[i][begidx_aft:endidx_aft:12,:,:])
            else:
                temp_array = varlist[i][begidx_aft:endidx_aft:12,:,:] - varlist[i][begidx_bef:endidx_bef:12,:,:]
#                 temp_array[temp_array<-2e-1] = np.nan
#                 temp_array[temp_array>2e-1] = np.nan
                var_logging_loss.append(temp_array[:,:,:])
                var_logging_deno.append(varlist[i][begidx_aft:endidx_aft:12,:,:])
#            # Variation (Wot_Feb - Wot_Jan)
#             if(modname[i][4:5] == 'c'):
#                 if(modname[i] != 'har_c'):
#                     if(weighted):
#                         temp_array = varlist[i+1][begidx_aft:endidx_aft:12,:,:] - varlist[i+1][begidx_bef:endidx_bef:12,:,:]
#                     else:
#                         temp_array = varlist[i+1][begidx_aft:endidx_aft:12,:,:] - varlist[i+1][begidx_bef:endidx_bef:12,:,:]
#                 else:
#                     if(weighted):
#                         temp_array = varlist[1][begidx_aft:endidx_aft:12,:,:] - varlist[1][begidx_bef:endidx_bef:12,:,:]
#                     else:
#                         temp_array = varlist[1][begidx_aft:endidx_aft:12,:,:] - varlist[1][begidx_bef:endidx_bef:12,:,:]
#             else:
#                 if(modname[i][4:5] == 'a'):
#                     temp_array = varlist[i-1][begidx_aft:endidx_aft:12,:,:] - varlist[i-1][begidx_bef:endidx_bef:12,:,:]
#             temp_array[temp_array<-2e-1] = np.nan
#             temp_array[temp_array>2e-1] = np.nan
#             if(weighted):
#                 var_logging_var.append(area_4x5[None,:,:] * temp_array[:,:,:])
#             else:                        
#                 var_logging_var.append(temp_array[:,:,:])
    for i in np.arange(0,len(modname_hrv)):
        if(weighted):
            for iyr in np.arange(0,(endyr-begyr)):
                var_logging_loss[i][iyr,:,:] = var_logging_loss[i][iyr,:,:]/area_4x5[:,:]# - var_logging_var[i])
                var_logging_deno[i][iyr,:,:] = var_logging_deno[i][iyr,:,:]/area_4x5[:,:]
        var_logging_loss[i][var_logging_loss[i]<-2e2] = np.nan
        var_logging_loss[i][var_logging_loss[i]>2e2] = np.nan
        var_logging_deno[i][var_logging_deno[i]<-2e2] = np.nan
        var_logging_deno[i][var_logging_deno[i]>2e2] = np.nan
        var_logging_global.append(np.nanmean(np.nanmean(var_logging_loss[i], 1), 1))
        var_logging_global_deno.append(np.nanmean(np.nanmean(var_logging_deno[i], 1), 1))
    if(output_ratio):
        var_logging_ratio = 100.0 * np.array(var_logging_global)/np.array(var_logging_global_deno)# - var_logging_var[i])
    else:
        var_logging_ratio = 0.0

    return [var_logging_loss, var_logging_global, var_logging_ratio]

def obtain_regrowth_ts(modname, modname_hrv, varlist, begyr, endyr, months):
    """
    Change of variables due to natural + CO2 fertilization, time series
    Can accept either a consecutive sequence of months [begmonth, endmonth] or a single month
    # modname - list contains name of all model cases
    # varlist - list contains the data of corresponding variables
    # begyr - the begin year of the data used
    # endyr - the end year of the data used
    # months - determine how many months of regrowth accounted for, the first is start month and
    #         the second is end month
    """
    var_regrowth = []
    var_regrowth_global = []
    for i in np.arange(0,case_len):
        if(modname[i][4:5] == 'c' or modname[i][4:5] == 'a'):
            begidx_bef = 12*begyr+months[0]-1
            endidx_bef = 12*endyr+months[0]-1
            if(len(months) > 1):
                # Multiple months
                begidx_aft = 12*begyr+months[1]-1
                endidx_aft = 12*endyr+months[1]-1    
            else:
                # One month
                begidx_aft = 12*begyr+months[0]
                endidx_aft = 12*endyr+months[0]
            var_regrowth.append(varlist[i][begidx_aft:endidx_aft:12,:,:] - varlist[i][begidx_bef:endidx_bef:12,:,:])        
    for i in np.arange(0, len(modname_hrv)):
        var_regrowth_global.append(np.nanmean(np.nanmean(var_regrowth[i], 1), 1))
    return [var_regrowth, var_regrowth_global]

def obtain_regrowth_impact_ts(modname, varlist, begyr, endyr):
    """
    Approximate regrowth rate from natural revovery by excluding CO2 fertilization
    The difference of Next Jan - Current Feb between with logging and without logging
    This method is not perfect since wt and wot cases are not starting from the same point 
    in each year January, it is only an approximation.
    # modname - list contains name of all model cases
    # varlist - list contains the data of corresponding variables
    # begyr - the begin year of the data used
    # endyr - the end year of the data used
    """
    var_regrowth = []
    var_regrowth_global = []
    var_diff_dec_jan = []
    # First calculate diff between Feb and next Jan
    begidx_bef = 12*begyr+1
    endidx_bef = 12*(endyr-1)+1
    begidx_aft = 12*begyr+12
    endidx_aft = 12*endyr
    for i in np.arange(0,case_len):
        var_diff_dec_jan.append(varlist[i][begidx_aft:endidx_aft:12,:,:] - varlist[i][begidx_bef:endidx_bef:12,:,:])  
    # Difference between with and without logging
    for i in np.arange(0,case_len):
        if(modname[i][4:5] == 'c'):
            if(modname[i] != 'har_c'):
                var_regrowth.append(var_diff_dec_jan[i][:,:,:] - var_diff_dec_jan[i+1][:,:,:])
            else:
                var_regrowth.append(var_diff_dec_jan[i][:,:,:] - var_diff_dec_jan[1][:,:,:])
        else:
            if(modname[i][4:5] == 'a'):
                var_regrowth.append(var_diff_dec_jan[i][:,:,:] - var_diff_dec_jan[i-1][:,:,:])
    for i in np.arange(0, len(modname_hrv)):
        var_regrowth_global.append(np.nanmean(np.nanmean(var_regrowth[i], 1), 1))
        
    return [var_regrowth, var_regrowth_global]

def qa_fates(modname, cr_area_cl):
    """
    Quality assurance for ELM-FATES cases only.
    Check all gridcells with a canopy coverage decrease larger than 10% under no harvest simulation
    For other logging cases we use 50% as a threshold.
    These outliers shall be removed when calculating the BGP effect due to a certain issue with FATES global simualtions. 
    Issue: Sudden decrease of canopy coverage until reaching the minimum value (17%) is triggered at random time.
    # modname - list contains name of all model cases
    # cr_area_cl - Crown area from model output
    """
    cr_area_1850s = []
    cr_area_2000s = []
    outliers = []
    case_len = len(modname)
    nyr = np.shape(cr_area_cl)[1]/12
    cr_area_yr = np.zeros((nyr, 46, 72))
    for i in np.arange(0,case_len):
        for iyr in np.arange(0,nyr):
            # canopy coverage (time x cl x lat x lon) at top layer (cl = 0)
            cr_area_yr[iyr,:,:] = np.nanmean(cr_area_cl[i][12*iyr:12*iyr+12,0,:,:], 0)
        cr_area_1850s.append(np.nanmean(cr_area_yr[0:10,:,:], 0))
        cr_area_2000s.append(np.nanmean(cr_area_yr[nyr-10:nyr,:,:], 0))
    for k in np.arange(0, case_len):
        cr_area_diff = cr_area_1850s[k] - cr_area_2000s[k]
        if(modname[k][4:5] == 'c' or modname[i][4:5] == 'a'):
            outliers.append(np.where(cr_area_diff > 0.5))
        else:
            outliers.append(np.where(cr_area_diff > 0.1))

    return outliers

def qa_elm(modname, tlai_elm):
    """
    Check all gridcells with a LAI change larger than 1.2 (about 40% of global mean LAI) when comparing harvest to no harvest simulation
    For ELM cases only
    These outliers shall be removed when calculating the BGP effect due to certain issues with FATES global simualtions.
    # modname - list contains name of all model cases
    # tlai_elm - total leaf area index from 2 ELM cases
    """
    case_len = len(modname)
    nyr = np.shape(tlai_elm)[1]/12
    diff_lai_pt = np.zeros((nyr, 46, 72))
    for iyr in np.arange(0, nyr):
        if(modname[0] == 'elm_n'):
            diff_lai_pt[iyr,:,:] = np.nanmean(tlai_elm[1][12*iyr:12*iyr+12,:,:]-tlai_elm[0][12*iyr:12*iyr+12,:,:], 0) 
        else:
            diff_lai_pt[iyr,:,:] = np.nanmean(tlai_elm[0][12*iyr:12*iyr+12,:,:]-tlai_elm[1][12*iyr:12*iyr+12,:,:], 0) 
    outliers = np.ma.where(abs(diff_lai_pt) > 1.2)
    
    return outliers

def calc_rlength_fates(modname, cr_area_cl, hite, z0mr=0.055, rl_soil=0.01):
    """
    Notably it is a simplfied approach that certain understory grass may have a bit larger roughness length comparing
    to soil, however for simplcity we ignore their contribution
    # modname - list contains name of all model cases
    # cr_area_cl - crown area
    # hite - crown area weighted canopy height
    # z0mr - ratio of momentum roughness length to change of canopy height
    # rl_soil -  roughness length of bareground soil
    """
    # crown area-weighted mean height of canopy plants
    # Forest + (grass & crop) part
    case_len = len(modname)
    nyr = np.shape(hite)[1]/12
    rl_full = []
    for i in np.arange(0,case_len):
        rl_per_case = np.zeros((12*nyr, 46, 72))
        rl_per_case[:,:,:] = (z0mr * hite[i][:,:,:] * cr_area_cl[i][:,0,:,:])
        rl_full.append(rl_per_case)

    return rl_full
    
def calc_rlength_elm(modname, htop, z0mr=0.055, rl_soil=0.01):
    """
    Notably it is a simplfied approach that certain understory grass may have a bit larger roughness length comparing
    to soil, however for simplcity we ignore their contribution
    # modname - list contains name of all model cases
    # cr_area_cl - crown area
    # hite - crown area weighted canopy height
    # z0mr - ratio of momentum roughness length to change of canopy height
    # rl_soil -  roughness length of bareground soil
    """
    case_len = len(modname)
    nyr = np.shape(htop)[1]/12
    rl_full = []
    # For ELM-only model
    for i in np.arange(0,case_len):
        rl_per_case = np.zeros((12*nyr, 46, 72))
        rl_per_case[:,:,:] = z0mr * htop[i][:,:,:]
        rl_full.append(rl_per_case)

    return rl_full

def moving_average(a, n=10) :
    """
    Quick moving average function
    # a - input array
    # n - the window of moving average
    """
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n
