###########################################################################
# read_data module
#
#	In support of NSST work
#
###########################################################################
#
# Import python modules
#

import xarray as xr
import os
from os.path import exists
import numpy as np
import glob
import matplotlib.pyplot as plt
import pickle
import datetime as dt
import sys
from netCDF4 import Dataset

from daynight_terminator import terminator

from tools_analysis import get_verif_date_anl
from tools_analysis import get_verif_date_f06

###########################################################################
#
# FUNCTIONS for reading input datasets
#

#pi = 4.0 * np.arctan(1.0)

#===============================================================================================
# Get OSTIA observations on ocean model grid
#
def read_ostia(yyyymmdd):

        # NOTE: These files were interpolated by someone else. There is no error information.
    ostia_fname = '/scratch1/NCEPDEV/da/Katherine.Lukens/NSST/data/OSTIA_MOM6_grid/from_Shastri/'+yyyymmdd+'12_ostia_2_global_0.25deg.nc'
    #ostia_fname = '/scratch1/NCEPDEV/da/Katherine.Lukens/NSST/data/OSTIA_MOM6_grid/ostia_'+yyyymmdd+'12_global_0.25deg.nc4'

    ds = xr.open_dataset(ostia_fname)

    # Get the SST data variable
    sst = np.squeeze(ds['analysed_sst'][:])		# np.squeeze: removes all axes of length 1 from array
    #sst = np.squeeze(ds['sst'][:])
    lat = np.squeeze(ds['lat'][:])
    lon = np.squeeze(ds['lon'][:])
 
    # Close the NetCDF file
    ds.close()

    return sst.values, lat, lon

#===============================================================================================
# Get OSTIA data from original files
#
def read_ostia_orig_files(yyyymmdd):

    ostia_fname = '/scratch1/NCEPDEV/da/Katherine.Lukens/NSST/data/OSTIA_orig_files/'+yyyymmdd+'/'+yyyymmdd+'120000-UKMO-L4_GHRSST-SSTfnd-OSTIA-GLOB-v02.0-fv02.0.nc'

    # open dataset
    ds = Dataset(ostia_fname, mode='r', format='NetCDF4 Classic')

    # extract variables
            # _FillValue = -32768
    sst         = np.asarray(ds.variables['analysed_sst'])              # analysed sea surface temperature: K (add offset 273.15), scale factor 0.01
    sst_error   = np.asarray(ds.variables['analysis_error'])            # estimated error standard deviation of analysed_sst: K (no offset added), scale factor 0.01
            # _FillValue = -128
    seaice      = np.asarray(ds.variables['sea_ice_fraction'])          # sea ice area fraction: no units, scale factor 0.01
    mask        = np.asarray(ds.variables['mask'])                      # land sea ice lake bit mask: 1 = water, 2 = land, 4 = optional_lake_surface, 8 = sea_ice, 16 = optional_river_surface

    time        = np.asarray(ds.variables['time'])                      # reference time of sst field: seconds since 1981-01-01 00:00:00
    lat         = np.asarray(ds.variables['lat'])                       # latitude: degrees_north | Latitude geographical coordinates,WGS84 projectio
    lon         = np.asarray(ds.variables['lon'])                       # longitude: degrees_east | Longitude geographical coordinates,WGS84 projectio

    # get file attributes
    attrs_sst    = []
    attrs_error  = []
    attrs_seaice = []
    attrs_mask   = []
    for name, variable in ds.variables.items():
      if name=='analysed_sst':
        for attrname in variable.ncattrs():
            attrs_str = str(attrname)+": "+str(getattr(variable, attrname))
            attrs_sst.append(attrs_str)
            del attrs_str
      elif name=='analysis_error':
        for attrname in variable.ncattrs():
            attrs_str = str(attrname)+": "+str(getattr(variable, attrname))
            attrs_error.append(attrs_str)
            del attrs_str
      elif name=='sea_ice_fraction':
        for attrname in variable.ncattrs():
            attrs_str = str(attrname)+": "+str(getattr(variable, attrname))
            attrs_seaice.append(attrs_str)
            del attrs_str
      elif name=='mask':
        for attrname in variable.ncattrs():
            attrs_str = str(attrname)+": "+str(getattr(variable, attrname))
            attrs_mask.append(attrs_str)
            del attrs_str

    ds.close()

    return sst, sst_error, seaice, mask, time, lat, lon, ostia_fname, attrs_sst, attrs_error, attrs_seaice, attrs_mask

#===============================================================================================
# Get diag data from netcdf
#       Note: data has already been binned onto lat/lon grid
#
def read_diags_ocean_gridded(dpath, yyyy, mm, dd, hour, exptname, obstype, degree, opt):

    if int(opt)==0:
      meansumstr = "means"
    elif int(opt)==1:
      meansumstr = "sums"

    diag_fname = dpath+'DIAGS.'+exptname+'.'+str(yyyy)+str(mm)+str(dd)+str(hour)+'.'+str(obstype)+'_obs.global_'+str(degree)+'deg.'+str(meansumstr)+'.nc4'
    print("read_diags_gridded: input file = "+str(diag_fname))

    if exists(diag_fname)==False:
      omb = "NO FILES"
      return omb,omb,omb,omb,omb,omb,omb,omb,omb,omb,omb,omb,omb,omb,omb
    elif exists(diag_fname)==True:
      print("read_diags_gridded: file exists")

      ds = xr.open_dataset(diag_fname)

      # Get the diag data
      omb         = np.squeeze(ds['omb'][:,:])             # np.squeeze: removes all axes of length 1 from array
      oma         = np.squeeze(ds['oma'][:,:])
      obserr      = np.squeeze(ds['obserr_postQC'][:,:])
      obserr_orig = np.squeeze(ds['obserr_beforeQC'][:,:])
      obsval      = np.squeeze(ds['obsvalue'][:,:])
      bkg         = np.squeeze(ds['bkg'][:,:])
      anl         = np.squeeze(ds['anl'][:,:])

      SDomb         = np.squeeze(ds['SDomb'][:,:])             # np.squeeze: removes all axes of length 1 from array
      SDoma         = np.squeeze(ds['SDoma'][:,:])
      SDobserr      = np.squeeze(ds['SDobserr_postQC'][:,:])
      SDobserr_orig = np.squeeze(ds['SDobserr_beforeQC'][:,:])
      SDobsval      = np.squeeze(ds['SDobsvalue'][:,:])
      SDbkg         = np.squeeze(ds['SDbkg'][:,:])
      SDanl         = np.squeeze(ds['SDanl'][:,:])

      nobs          = np.squeeze(ds['nobs'][:,:])

      # Close the NetCDF file
      ds.close()

      return nobs, omb, SDomb, oma, SDoma, obserr, SDobserr, obserr_orig, SDobserr_orig, obsval, SDobsval, bkg, SDbkg, anl, SDanl

#------------------------------------------------
# Read ocean model background error (SD)
def read_bkgerr(exptpath, yyyymmdd, hh):

    f = exptpath+'gdas.'+str(yyyymmdd)+'/'+str(hh)+'/analysis/ocean/gdas.t'+str(hh)+'z.ocn.bkgerr_stddev.nc'

    ds = xr.open_dataset(f)

        # read surface T and convert to deg-C
    #tT = np.squeeze(ds['Temp'][0,0,:,:])      # top layer T
    tT = np.squeeze(ds['Temp'][0,:,:,:])      # all layers

    T = np.where(tT==0, np.nan, tT)

    #lat  = np.squeeze(ds['yaxis_1'].values)
    #lon  = np.squeeze(ds['xaxis_1'].values)

    ds.close()

    return T#, lat, lon

#------------------------------------------------
# Read horizontal and vertical decorrelation length scale file
def read_lengthscales(exptpath, yyyymmdd, hh):

    #f = '/scratch1/NCEPDEV/stmp2/Katherine.Lukens/RUNDIRS/converter_test2/gdasocnanal_06/ocn.cor_rh.incr.0001-01-01T00:00:00Z.nc'
    f = exptpath+'gdasocnanal_12/vt_scales.nc'

    ds = xr.open_dataset(f)

        # read surface T and convert to deg-C
    #tT = np.squeeze(ds['ave_ssh'][0,:,:])      # top layer T
    thz = np.squeeze(ds['hz'][:,:])
    tvt = np.squeeze(ds['vt'][:,:,:])

    #T = np.where(tT==0, np.nan, tT)
    hz = np.where(thz==0, np.nan, thz)

    shapevt = np.shape(tvt)
    vt = tvt
    for i in range(shapevt[0]):
      vt[i,:,:] = np.where(tvt[0,:,:]==0, np.nan, tvt[i,:,:])

    ds.close()
    del f

        # get lat and lon 1D arrays (for plotting)
    f = exptpath+'gdasocnanal_12/soca_gridspec.nc'
    ds = xr.open_dataset(f)
    lats = np.squeeze(ds['lath'][0,:]).values
    lons = np.squeeze(ds['lonh'][0,:]).values

    return hz, vt, lats, lons

#------------------------------------------------
# Read horizontal decorrelation length scale file
def read_lengthscale_horiz(exptpath, yyyymmdd, hh):

    f = '/scratch1/NCEPDEV/stmp2/Katherine.Lukens/RUNDIRS/converter_test2/gdasocnanal_06/ocn.cor_rh.incr.0001-01-01T00:00:00Z.nc'

    ds = xr.open_dataset(f)

        # read surface T and convert to deg-C
    tT = np.squeeze(ds['ave_ssh'][0,:,:])      # top layer T

    T = np.where(tT==0, np.nan, tT)

    #lat  = np.squeeze(ds['yaxis_1'].values)
    #lon  = np.squeeze(ds['xaxis_1'].values)

    ds.close()

    return T

#------------------------------------------------
# Read vertical decorrelation length scale file
def read_lengthscale_vert(exptpath, yyyymmdd, hh):

    f = '/scratch1/NCEPDEV/stmp2/Katherine.Lukens/RUNDIRS/converter_test2/gdasocnanal_06/ocn.cor_rv.incr.0001-01-01T00:00:00Z.nc'

    ds = xr.open_dataset(f)

        # read surface T and convert to deg-C
    tT = np.squeeze(ds['ave_ssh'][0,:,:])      # top layer T

    T = np.where(tT==0, np.nan, tT)

    #lat  = np.squeeze(ds['yaxis_1'].values)
    #lon  = np.squeeze(ds['xaxis_1'].values)

    ds.close()

    return T

#------------------------------------------------
# Extract skin T (from atmosphere analysis)
#       On atmosphere grid
def read_anl_skinT(exptpath, yyyymmdd, hh):

    fskin = exptpath+'gdas.'+str(yyyymmdd)+'/'+str(hh)+'/analysis/atmos/gdas.t'+str(hh)+'z.sfcanl.nc'

    ds = xr.open_dataset(fskin)

	# read surface T and convert to deg-C
    T = np.squeeze(ds['tmpsfc'][0,:,:].values) - 273.15

    lat  = np.squeeze(ds['grid_yt'].values)
    lon  = np.squeeze(ds['grid_xt'].values)

    ds.close()

    return T, lat, lon

#------------------------------------------------
# Extract skin T (from atmosphere surface backgrounds)
#       On atmosphere grid
def read_bkg_skinT(exptpath, yyyymmdd, hh):

    #fskin = exptpath+'gdas.'+str(yyyymmdd)+'/'+str(hh)+'/analysis/atmos/gdas.t'+str(hh)+'z.sfcanl.nc'
    fskin003 = exptpath+'gdas.t'+str(hh)+'z.sfcf003.nc'
    fskin006 = exptpath+'gdas.t'+str(hh)+'z.sfcf006.nc'
    fskin009 = exptpath+'gdas.t'+str(hh)+'z.sfcf009.nc'

    ds003 = xr.open_dataset(fskin003)
    ds006 = xr.open_dataset(fskin006)
    ds009 = xr.open_dataset(fskin009)

	# read surface T and convert to deg-C
    T003 = np.squeeze(ds003['tmpsfc'][0,:,:].values) - 273.15
    T006 = np.squeeze(ds006['tmpsfc'][0,:,:].values) - 273.15
    T009 = np.squeeze(ds009['tmpsfc'][0,:,:].values) - 273.15

    lat  = np.squeeze(ds006['grid_yt'].values)
    lon  = np.squeeze(ds006['grid_xt'].values)

    ds003.close()
    ds006.close()
    ds009.close()

    return T003, T006, T009, lat, lon

#===============================================================================================
# Get ocean model background SST (for time series plots)
#	Output daily mean SST
def read_ocnbkg_dailymean(exptpath, yyyy, mm, dd):

    #gdas_fnames = glob.glob(exptpath+'gdas.'+yyyy+mm+dd+'/??/model_data/ocean/history/gdas.t??z.ocnf*.nc')
    gdas_fnames = glob.glob(exptpath+'gdas.'+yyyy+mm+dd+'/??/model_data/ocean/history/gdas.t??z.ocnf006.nc')

    if np.size(gdas_fnames)==0: return "NO FILES"

    cnt = 0.0
    for gdas_fname in gdas_fnames:
      ds = xr.open_dataset(gdas_fname)
      if cnt == 0:
        sst = np.squeeze(ds['Temp'][0,0,:,:])
      else:
        sst += np.squeeze(ds['Temp'][0,0,:,:])

      cnt +=1

      ds.close()

    return sst.values/cnt

#===============================================================================================
# Get ocean analysis SST (for time series plots)
#       Output daily mean SST
def read_ocnanl_dailymean(exptpath, yyyy, mm, dd):

    gdas_fnames = glob.glob(exptpath+'gdas.'+yyyy+mm+dd+'/??/model_data/ocean/history/gdas.t??z.ocnana.nc')

    if np.size(gdas_fnames)==0: return "NO FILES"

    cnt = 0.0
    for gdas_fname in gdas_fnames:
      ds = xr.open_dataset(gdas_fname)
      if cnt == 0:
        sst = np.squeeze(ds['Temp'][0,0,:,:])
      else:
        sst += np.squeeze(ds['Temp'][0,0,:,:])

      cnt +=1

      ds.close()

    return sst.values/cnt

#===============================================================================================
# Get diag data from netcdf
#       Note: data are in obs space (each point corresponds to an observation 4D location)
#
#	opt: option to extract day, night, or all points
#               = -1 ... all points
#		=  0 ... day points
#		=  1 ... night points
def read_diags_ocean(exptpath, yyyy, mm, dd, hour, obstype, nlons, opt):

    if obstype=="all_sats":
        # all SST obs diags (from satellite observing instruments)
      fnames = 'gdas.'+yyyy+mm+dd+'/'+str(hour)+'/analysis/ocean/diags/sst*.'+yyyy+mm+dd+str(hour)+'.nc4'
      gdas_fnames = glob.glob(exptpath+fnames)
    elif obstype=="all_geo":
      fname0 = exptpath+'gdas.'+yyyy+mm+dd+'/'+str(hour)+'/analysis/ocean/diags/sst_abi_g16_l3c.'+yyyy+mm+dd+str(hour)+'.nc4'
      fname1 = exptpath+'gdas.'+yyyy+mm+dd+'/'+str(hour)+'/analysis/ocean/diags/sst_abi_g17_l3c.'+yyyy+mm+dd+str(hour)+'.nc4'
      fname2 = exptpath+'gdas.'+yyyy+mm+dd+'/'+str(hour)+'/analysis/ocean/diags/sst_ahi_h08_l3c.'+yyyy+mm+dd+str(hour)+'.nc4'
      gdas_fnames = [fname0]
      gdas_fnames.append(fname1)
      gdas_fnames.append(fname2)
    elif obstype=="all_leo":
      fname0 = exptpath+'gdas.'+yyyy+mm+dd+'/'+str(hour)+'/analysis/ocean/diags/sst_avhrr_ma_l3u.'+yyyy+mm+dd+str(hour)+'.nc4'
      fname1 = exptpath+'gdas.'+yyyy+mm+dd+'/'+str(hour)+'/analysis/ocean/diags/sst_avhrr_mb_l3u.'+yyyy+mm+dd+str(hour)+'.nc4'
      fname2 = exptpath+'gdas.'+yyyy+mm+dd+'/'+str(hour)+'/analysis/ocean/diags/sst_viirs_n20_l3u.'+yyyy+mm+dd+str(hour)+'.nc4'
      fname3 = exptpath+'gdas.'+yyyy+mm+dd+'/'+str(hour)+'/analysis/ocean/diags/sst_viirs_npp_l3u.'+yyyy+mm+dd+str(hour)+'.nc4'
      gdas_fnames = [fname0]
      gdas_fnames.append(fname1)
      gdas_fnames.append(fname2)
      gdas_fnames.append(fname3)
    else:
        # individual obs diags
      fnames = 'gdas.'+yyyy+mm+dd+'/'+str(hour)+'/analysis/ocean/diags/sst_'+str(obstype)+'*.'+yyyy+mm+dd+str(hour)+'.nc4'
      gdas_fnames = glob.glob(exptpath+fnames)

    print("get_omba: gdas_fnames = "+str(gdas_fnames)+" size = "+str(np.size(gdas_fnames)))

    if np.size(gdas_fnames)==0:
      omb = "NO FILES"
      return omb,omb,omb,omb,omb,omb,omb,omb,omb

    cnt = 0.0
    for gdas_fname in gdas_fnames:
     if exists(gdas_fname):
      ds = Dataset(gdas_fname)
      print("get_omba: gdas_fname = "+str(gdas_fname)+" | count = "+str(cnt))

      #---------------------------------------
      # Extract innovations and other variables
      
      omb    = np.squeeze(ds.groups['ombg'].variables['seaSurfaceTemperature'][:])
      oma    = np.squeeze(ds.groups['oman'].variables['seaSurfaceTemperature'][:])
      obserr = np.squeeze(ds.groups['EffectiveError1'].variables['seaSurfaceTemperature'][:])   # obs error after QC, inflation, BC, etc.
      obserr_orig = np.squeeze(ds.groups['ObsError'].variables['seaSurfaceTemperature'][:])     # original obs error (from data producers)
      obsval = np.squeeze(ds.groups['ObsValue'].variables['seaSurfaceTemperature'][:])
      qc     = np.squeeze(ds.groups['EffectiveQC1'].variables['seaSurfaceTemperature'][:])        # 0 = obs are accepted/assimilated
      time   = np.squeeze(ds.groups['MetaData'].variables['dateTime'][:])                 # seconds since 1970-01-01T00:00:00Z
      lats   = np.squeeze(ds.groups['MetaData'].variables['latitude'][:])
      lons   = np.squeeze(ds.groups['MetaData'].variables['longitude'][:])

      #---------------------------------------
      # Only include accepted (assimilated) obs
      
      idx = np.where(qc==0)

      qomb    = omb[idx]
      qoma    = oma[idx]
      qobserr = obserr[idx]
      qobserr_orig = obserr_orig[idx]
      qobsval = obsval[idx]
      qqc     = qc[idx]
      qtime   = time[idx]
      qlats   = lats[idx]
      qlons   = lons[idx]
      
      del idx

      #---------------------------------------
      # Get day and night lat/lon points
      # 	checks all the lats at longitude ilon and sets them to 1 if they are < termlats at longitude ilon: day=1, night=0
      #			ex) daynight[:,ilon] = np.where(lat[:,ilon] < termlats[ilon], 0, daynight[:,ilon])

      if opt>=0:
        lonmin = -180
        lonmax = 180
        termlats, termlons = terminator(yyyy, mm, dd, str(hour), "00", lonmin, lonmax, nlons)
	
        if opt==0:
	  # day points
          for ilon in range(nlons):
            qomb    = np.where(qlats<termlats[ilon], np.nan, qomb)
            qoma    = np.where(qlats<termlats[ilon], np.nan, qoma)
            qobserr = np.where(qlats<termlats[ilon], np.nan, qobserr)
            qobserr_orig = np.where(qlats<termlats[ilon], np.nan, qobserr_orig)
            qobsval = np.where(qlats<termlats[ilon], np.nan, qobsval)
            qqc     = np.where(qlats<termlats[ilon], np.nan, qqc)
            qtime   = np.where(qlats<termlats[ilon], np.nan, qtime)
            qlons   = np.where(qlats<termlats[ilon], np.nan, qlons)
            qlats   = np.where(qlats<termlats[ilon], np.nan, qlats)
            print("day size: "+str(np.size(qomb)))
                # keep only non-missing values
            qomb    = qomb[~np.isnan(qomb)]
            qoma    = qoma[~np.isnan(qoma)]
            qobserr = qobserr[~np.isnan(qobserr)]
            qobserr_orig = qobserr_orig[~np.isnan(qobserr_orig)]
            qobsval = qobsval[~np.isnan(qobsval)]
            qqc     = qqc[~np.isnan(qqc)]
            qtime   = qtime[~np.isnan(qtime)]
            qlats   = qlats[~np.isnan(qlats)]
            qlons   = qlons[~np.isnan(qlons)]
            print("day non-miss size: "+str(np.size(qomb)))

        elif opt==1:
	  # night points
          for ilon in range(nlons):
            qomb    = np.where(qlats>termlats[ilon], np.nan, qomb)
            qoma    = np.where(qlats>termlats[ilon], np.nan, qoma)
            qobserr = np.where(qlats>termlats[ilon], np.nan, qobserr)
            qobserr_orig = np.where(qlats>termlats[ilon], np.nan, qobserr_orig)
            qobsval = np.where(qlats>termlats[ilon], np.nan, qobsval)
            qqc     = np.where(qlats>termlats[ilon], np.nan, qqc)
            qtime   = np.where(qlats>termlats[ilon], np.nan, qtime)
            qlons   = np.where(qlats>termlats[ilon], np.nan, qlons)     
            qlats   = np.where(qlats>termlats[ilon], np.nan, qlats)
            print("night size: "+str(np.size(qomb)))
                # keep only non-missing values
            qomb    = qomb[~np.isnan(qomb)]
            qoma    = qoma[~np.isnan(qoma)]
            qobserr = qobserr[~np.isnan(qobserr)]
            qobserr_orig = qobserr_orig[~np.isnan(qobserr_orig)]
            qobsval = qobsval[~np.isnan(qobsval)]
            qqc     = qqc[~np.isnan(qqc)]
            qtime   = qtime[~np.isnan(qtime)]
            qlats   = qlats[~np.isnan(qlats)]
            qlons   = qlons[~np.isnan(qlons)]
            print("night non-miss size: "+str(np.size(qomb)))

        del termlats, termlons

      #---------------------------------------
      # Append to total arrays

      if cnt==0:
        omb_arr    = qomb
        oma_arr    = qoma
        obserr_arr = qobserr
        obserr_orig_arr = qobserr_orig
        obsval_arr = qobsval
        qc_arr     = qqc
        time_arr   = qtime
        lats_arr   = qlats
        lons_arr   = qlons
      else:
        omb_arr    = np.append(omb_arr, qomb, axis=0)
        oma_arr    = np.append(oma_arr, qoma, axis=0)
        obserr_arr = np.append(obserr_arr, qobserr, axis=0)
        obserr_orig_arr = np.append(obserr_orig_arr, qobserr_orig, axis=0)
        obsval_arr = np.append(obsval_arr, qobsval, axis=0)
        qc_arr     = np.append(qc_arr, qqc, axis=0)
        time_arr   = np.append(time_arr, qtime, axis=0)
        lats_arr   = np.append(lats_arr, qlats, axis=0)
        lons_arr   = np.append(lons_arr, qlons, axis=0)

      #---------------------------------------

      del omb,oma,obserr,obserr_orig,obsval,qc,time,lats,lons
      del qomb,qoma,qobserr,qobserr_orig,qobsval,qqc,qtime,qlats,qlons

      cnt += 1

      ds.close()

    if cnt==0:
      omb = "NO FILES"
      return omb,omb,omb,omb,omb,omb,omb,omb,omb
    else:
      return omb_arr, oma_arr, obserr_arr, obserr_orig_arr, obsval_arr, qc_arr, time_arr, lats_arr, lons_arr

#------------------------------------------------
# Extract ocean analysis increment
def read_ocninc(exptpath, yyyymmdd, hh):

    finc = exptpath+'gdas.'+str(yyyymmdd)+'/'+str(hh)+'/analysis/ocean/gdas.t'+str(hh)+'z.ocninc.nc'

    ds = xr.open_dataset(finc)

        # read surface T and convert to deg-C
    Tinc = np.squeeze(ds['Temp'][0,:,:].values)

    ds.close()

    return Tinc

#===============================================================================================
# Get ocean model background SST (f06 background forecast) - for time series (TS) plots
#
def read_ocnbkg_TS(exptpath, yyyy, mm, dd, filetype):

    # '/scratch2/NCEPDEV/ocean/Guillaume.Vernieres/runs/cp0/cp0.b/COMROT/cp0.b/gdas.'+yyyy+mm+dd+'/??/model_data/ocean/history/gdas.t??z.ocnf006.nc'
    if filetype=="all_fcsts":
      gdas_fnames = glob.glob(exptpath+'gdas.'+yyyy+mm+dd+'/??/model_data/ocean/history/gdas.ocean.t??z.inst.f00*.nc')
    elif filetype=="f06":
      gdas_fnames = glob.glob(exptpath+'gdas.'+yyyy+mm+dd+'/??/model_data/ocean/history/gdas.ocean.t??z.inst.f006.nc')
    print("get_bkg: gdas_fnames = "+str(gdas_fnames))

    if np.size(gdas_fnames)==0: return "NO FILES"

    cnt = 0.0
    for gdas_fname in gdas_fnames:
        ds = xr.open_dataset(gdas_fname)
        if cnt == 0:
            try:
                sst = np.squeeze(ds['Temp'][0,0,:,:])
                #sst = np.squeeze(ds['trefanlocn'][:,:]) - 273.15

                #sst = np.squeeze(ds['sst'][0,:,:])
            except:
                try:
                    sst = ds['TS_FOUND'][:,:] - 273.15
                except:
                    sst = ds['trefanlocn'][:,:] - 273.15
        else:
            try:
                sst += np.squeeze(ds['Temp'][0,0,:,:])
                #sst += np.squeeze(ds['trefanlocn'][:,:]) - 273.15

                #sst += np.squeeze(ds['sst'][0,:,:])
            except:
                try:
                    sst += ds['TS_FOUND'][:,:] - 273.15
                except:
                    sst += ds['trefanlocn'][:,:] - 273.15

#                sst += ds['TS_FOUND'][:,:] -273.15

        cnt +=1
        ds.close()
        #print(gdas_fname, ' ', cnt)
    return sst.values/cnt

#===============================================================================================
# Get ocean analysis SST - for time series plots
#
def read_ocnanl_TS(exptpath, yyyy, mm, dd):

    # '/scratch2/NCEPDEV/ocean/Guillaume.Vernieres/runs/cp0/cp0.b/COMROT/cp0.b/gdas.'+yyyy+mm+dd+'/??/analysis/ocean/gdas.t??z.ocnana.nc'
    gdas_fnames = glob.glob(exptpath+"gdas."+yyyy+mm+dd+"/??/analysis/ocean/gdas.t??z.ocnana.nc")
    print("get_anl: gdas_fnames = "+str(gdas_fnames))

    if np.size(gdas_fnames)==0: return "NO FILES"

    #tfilearr = []
    cnt = 0.0
    for gdas_fname in gdas_fnames:
      #tfile = exists(gdas_fname)
      #print("get_gsianl: tfile = "+str(tfile))
      #if tfile==True:
        ds = xr.open_dataset(gdas_fname)
        if cnt == 0:
            try:
                #sst = np.squeeze(ds[sfcTvar][:,:]) - 273.15
                sst = np.squeeze(ds['Temp'][0,0,:,:])
                #print("get_gsianl: get Temp")

            except:
            #    try:
                    sst = ds['TS_FOUND'][:,:] - 273.15
            #    except:
            #        sst = ds['trefanlocn'][:,:] - 273.15
        else:
            try:
                #sst += np.squeeze(ds[sfcTvar][:,:]) - 273.15
                sst += np.squeeze(ds['Temp'][0,0,:,:])
                #print("get_gsianl: get Temp")

            except:
            #    try:
                    sst += ds['TS_FOUND'][:,:] - 273.15
            #    except:
            #        sst += ds['trefanlocn'][:,:] - 273.15

#                sst += ds['TS_FOUND'][:,:] -273.15

        cnt +=1
        ds.close()

      #tfilearr.append(tfile)
      #if np.array(tfilearr).all()==False:
      #  fexist = False
      #else:
      #  fexist = True
      #print("get_gsianl: fexist = "+str(fexist))
      #del tfile

    return sst.values/cnt#, fexist

#===============================================================================================
# Compute daily mean and get T from each file called
#	Get SST or Skin T
#
#	Return: Daily mean, Data for each file opened, approx. sun position, lat and lon arrays
#
def get_dailymean_eachfile(dpath, yyyy, mm, dd, hours, filetype, var):

    if var=="temp": varname = 'Temp'            # ocean temperature, all levels
    if var=="sst": varname = 'Temp'             # coean temperature, top layer (SST)
    if var=="skint" or var=="atmosT": varname = 'tmpsfc'         # ocean temperature, surface (from atmos)
    
	# date before yyyymmdd
            # month
    if int(dd)==1:
      tmmB = int(mm) - 1
    else:
      tmmB = int(mm)
    if tmmB<10:
      mmB = "0"+str(tmmB)
    else:
      mmB = str(tmmB)
            # day
    if mmB==mm:
      tddB = int(dd) - 1
      if tddB<10:
        ddB = "0"+str(tddB)
      else:
        ddB = str(tddB)
    else:
      if mmB=="04" or mmB=="06" or mmB=="09" or mmB=="11":
        ddB = "30"
      elif mmB=="02":
        if int(yyyy)%4==0:
          ddB = "29"
        else:
          ddB = "28"
      else:
        ddB = "31"
            # hour
    HHarr = ['00','06','12','18']
    if hours=='00': hhB = '18'
    if hours=='06': hhB = '00'
    if hours=='12': hhB = '06'
    if hours=='18': hhB = '12'

	# get filenames
    if filetype=="all_fcsts":		# all background forecasts

      if dpath.find("role")!=-1 or dpath.find("marine_candidate")!=-1:
        if dpath.find("marine_candidate")!=-1:
          modelstr = "model"
        else:
          modelstr = "model_data"
        #dpathB = dpath.replace(yyyy+mm+dd+hours, yyyy+mmB+ddB+hhB)
        if var=="sst" or var=="temp":
          # before date
          Bfnames18_f006 = dpath+'/'+yyyy+mmB+ddB+'18/gdas.'+yyyy+mmB+ddB+'/18/'+modelstr+'/ocean/history/gdas.ocean.t18z.inst.f006.nc'
          Bfnames18_f009 = dpath+'/'+yyyy+mmB+ddB+'18/gdas.'+yyyy+mmB+ddB+'/18/'+modelstr+'/ocean/history/gdas.ocean.t18z.inst.f009.nc'
          # current date
          fnames00_f003 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/'+modelstr+'/ocean/history/gdas.ocean.t00z.inst.f003.nc'
          fnames00_f006 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/'+modelstr+'/ocean/history/gdas.ocean.t00z.inst.f006.nc'
          fnames00_f009 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/'+modelstr+'/ocean/history/gdas.ocean.t00z.inst.f009.nc'
          fnames06_f003 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/'+modelstr+'/ocean/history/gdas.ocean.t06z.inst.f003.nc'
          fnames06_f006 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/'+modelstr+'/ocean/history/gdas.ocean.t06z.inst.f006.nc'
          fnames06_f009 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/'+modelstr+'/ocean/history/gdas.ocean.t06z.inst.f009.nc'
          fnames12_f003 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/'+modelstr+'/ocean/history/gdas.ocean.t12z.inst.f003.nc'
          fnames12_f006 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/'+modelstr+'/ocean/history/gdas.ocean.t12z.inst.f006.nc'
          fnames12_f009 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/'+modelstr+'/ocean/history/gdas.ocean.t12z.inst.f009.nc'
          fnames18_f003 = dpath+'/'+yyyy+mm+dd+'18/gdas.'+yyyy+mm+dd+'/18/'+modelstr+'/ocean/history/gdas.ocean.t18z.inst.f003.nc'
        elif var=="skint" or var=="atmosT":
          # before date
          Bfnames18_f006 = dpath+'/'+yyyy+mmB+ddB+'18/gdas.'+yyyy+mmB+ddB+'/18/'+modelstr+'/atmos/history/gdas.t18z.sfcf006.nc'
          Bfnames18_f009 = dpath+'/'+yyyy+mmB+ddB+'18/gdas.'+yyyy+mmB+ddB+'/18/'+modelstr+'/atmos/history/gdas.t18z.sfcf009.nc'
          # current date
          fnames00_f003 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/'+modelstr+'/atmos/history/gdas.t00z.sfcf003.nc'
          fnames00_f006 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/'+modelstr+'/atmos/history/gdas.t00z.sfcf006.nc'
          fnames00_f009 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/'+modelstr+'/atmos/history/gdas.t00z.sfcf009.nc'
          fnames06_f003 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/'+modelstr+'/atmos/history/gdas.t06z.sfcf003.nc'
          fnames06_f006 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/'+modelstr+'/atmos/history/gdas.t06z.sfcf006.nc'
          fnames06_f009 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/'+modelstr+'/atmos/history/gdas.t06z.sfcf009.nc'
          fnames12_f003 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/'+modelstr+'/atmos/history/gdas.t12z.sfcf003.nc'
          fnames12_f006 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/'+modelstr+'/atmos/history/gdas.t12z.sfcf006.nc'
          fnames12_f009 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/'+modelstr+'/atmos/history/gdas.t12z.sfcf009.nc'
          fnames18_f003 = dpath+'/'+yyyy+mm+dd+'18/gdas.'+yyyy+mm+dd+'/18/'+modelstr+'/atmos/history/gdas.t18z.sfcf003.nc'
      else: 
       if var=="sst" or var=="temp":
          # before date
        Bfnames18_f006 = dpath+'gdas.'+yyyy+mmB+ddB+'/18/model_data/ocean/history/gdas.ocean.t18z.inst.f006.nc'
        Bfnames18_f009 = dpath+'gdas.'+yyyy+mmB+ddB+'/18/model_data/ocean/history/gdas.ocean.t18z.inst.f009.nc'
	  # current date
        fnames00_f003 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.ocean.t00z.inst.f003.nc'
        fnames00_f006 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.ocean.t00z.inst.f006.nc'
        fnames00_f009 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.ocean.t00z.inst.f009.nc'
        fnames06_f003 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.ocean.t06z.inst.f003.nc'
        fnames06_f006 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.ocean.t06z.inst.f006.nc'
        fnames06_f009 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.ocean.t06z.inst.f009.nc'
        fnames12_f003 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.ocean.t12z.inst.f003.nc'
        fnames12_f006 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.ocean.t12z.inst.f006.nc'
        fnames12_f009 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.ocean.t12z.inst.f009.nc'
        fnames18_f003 = dpath+'gdas.'+yyyy+mm+dd+'/18/model_data/ocean/history/gdas.ocean.t18z.inst.f003.nc'
       elif var=="skint" or var=="atmosT":
          # before date
        Bfnames18_f006 = dpath+'gdas.'+yyyy+mmB+ddB+'/18/model_data/atmos/history/gdas.t18z.sfcf006.nc'
        Bfnames18_f009 = dpath+'gdas.'+yyyy+mmB+ddB+'/18/model_data/atmos/history/gdas.t18z.sfcf009.nc'
          # current date
        fnames00_f003 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/atmos/history/gdas.t00z.sfcf003.nc'
        fnames00_f006 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/atmos/history/gdas.t00z.sfcf006.nc'
        fnames00_f009 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/atmos/history/gdas.t00z.sfcf009.nc'
        fnames06_f003 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/atmos/history/gdas.t06z.sfcf003.nc'
        fnames06_f006 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/atmos/history/gdas.t06z.sfcf006.nc'
        fnames06_f009 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/atmos/history/gdas.t06z.sfcf009.nc'
        fnames12_f003 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/atmos/history/gdas.t12z.sfcf003.nc'
        fnames12_f006 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/atmos/history/gdas.t12z.sfcf006.nc'
        fnames12_f009 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/atmos/history/gdas.t12z.sfcf009.nc'
        fnames18_f003 = dpath+'gdas.'+yyyy+mm+dd+'/18/model_data/atmos/history/gdas.t18z.sfcf003.nc'

      if np.size(hours)==1:
        if hours=='00':
          tfnames = [Bfnames18_f006]
          tfnames.append(Bfnames18_f009)
          tfnames.append(fnames00_f003)
        elif hours=='06':
          tfnames = [fnames00_f006]
          tfnames.append(fnames00_f009)
          tfnames.append(fnames06_f003)
        elif hours=='12':
          tfnames = [fnames06_f006]
          tfnames.append(fnames06_f009)
          tfnames.append(fnames12_f003)
        elif hours=='18':
          tfnames = [fnames12_f006]
          tfnames.append(fnames12_f009)
          tfnames.append(fnames18_f003)

      fnamesDAILY = [Bfnames18_f006]
      fnamesDAILY.append(Bfnames18_f009)
      fnamesDAILY.append(fnames00_f003)
      fnamesDAILY.append(fnames00_f006)
      fnamesDAILY.append(fnames00_f009)
      fnamesDAILY.append(fnames06_f003)
      fnamesDAILY.append(fnames06_f006)
      fnamesDAILY.append(fnames06_f009)
      fnamesDAILY.append(fnames12_f003)
      fnamesDAILY.append(fnames12_f006)
      fnamesDAILY.append(fnames12_f009)
      fnamesDAILY.append(fnames18_f003)
     
      del Bfnames18_f006, Bfnames18_f009
      del fnames00_f003, fnames00_f006, fnames00_f009
      del fnames06_f003, fnames06_f006, fnames06_f009
      del fnames12_f003, fnames12_f006, fnames12_f009
      del fnames18_f003

    elif filetype=="f06":		# only 6-h background forecasts
      if dpath.find("role")!=-1 or dpath.find("marine_candidate")!=-1:
       if dpath.find("marine_candidate")!=-1:
          modelstr = "model"
       else:
          modelstr = "model_data"
       if var=="sst" or var=="temp":
          # before date
        Bfnames18_f006 = dpath+'/'+yyyy+mmB+ddB+'18/gdas.'+yyyy+mmB+ddB+'/18/'+modelstr+'/ocean/history/gdas.ocean.t18z.inst.f006.nc'
	  # current date
        fnames00_f006 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/'+modelstr+'/ocean/history/gdas.ocean.t00z.inst.f006.nc'
        fnames06_f006 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/'+modelstr+'/ocean/history/gdas.ocean.t06z.inst.f006.nc'
        fnames12_f006 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/'+modelstr+'/ocean/history/gdas.ocean.t12z.inst.f006.nc'
       elif var=="skint" or var=="atmosT":
          # before date
        Bfnames18_f006 = dpath+'/'+yyyy+mmB+ddB+'18/gdas.'+yyyy+mmB+ddB+'/18/'+modelstr+'/atmos/history/gdas.t18z.sfcf006.nc'
          # current date
        fnames00_f006 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/'+modelstr+'/atmos/history/gdas.t00z.sfcf006.nc'
        fnames06_f006 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/'+modelstr+'/atmos/history/gdas.t06z.sfcf006.nc'
        fnames12_f006 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/'+modelstr+'/atmos/history/gdas.t12z.sfcf006.nc'
      else:
       if var=="sst" or var=="temp":
          # before date
        Bfnames18_f006 = dpath+'gdas.'+yyyy+mmB+ddB+'/18/model_data/ocean/history/gdas.ocean.t18z.inst.f006.nc'
	  # current date
        fnames00_f006 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.ocean.t00z.inst.f006.nc'
        fnames06_f006 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.ocean.t06z.inst.f006.nc'
        fnames12_f006 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.ocean.t12z.inst.f006.nc'
       elif var=="skint" or var=="atmosT":
          # before date
        Bfnames18_f006 = dpath+'gdas.'+yyyy+mmB+ddB+'/18/model_data/atmos/history/gdas.t18z.sfcf006.nc'
          # current date
        fnames00_f006 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/atmos/history/gdas.t00z.sfcf006.nc'
        fnames06_f006 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/atmos/history/gdas.t06z.sfcf006.nc'
        fnames12_f006 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/atmos/history/gdas.t12z.sfcf006.nc'

      if np.size(hours)==1:
        if hours=='00':
          tfnames = [Bfnames18_f006]
        elif hours=='06':
          tfnames = [fnames00_f006]
        elif hours=='12':
          tfnames = [fnames06_f006]
        elif hours=='18':
          tfnames = [fnames12_f006]

      fnamesDAILY = [Bfnames18_f006]
      fnamesDAILY.append(fnames00_f006)
      fnamesDAILY.append(fnames06_f006)
      fnamesDAILY.append(fnames12_f006)
     
    elif filetype=="anl":               # analysis files
      if dpath.find("role")!=-1:
       if var=="sst" or var=="temp":
        fnames00_ana = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/analysis/ocean/gdas.t00z.ocnana.nc'
        fnames06_ana = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/analysis/ocean/gdas.t06z.ocnana.nc'
        fnames12_ana = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/analysis/ocean/gdas.t12z.ocnana.nc'
        fnames18_ana = dpath+'/'+yyyy+mm+dd+'18/gdas.'+yyyy+mm+dd+'/18/analysis/ocean/gdas.t18z.ocnana.nc'
       elif var=="skint" or var=="atmosT":
        fnames00_ana = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/analysis/atmos/gdas.t00z.sfcanl.nc'
        fnames06_ana = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/analysis/atmos/gdas.t06z.sfcanl.nc'
        fnames12_ana = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/analysis/atmos/gdas.t12z.sfcanl.nc'
        fnames18_ana = dpath+'/'+yyyy+mm+dd+'18/gdas.'+yyyy+mm+dd+'/18/analysis/atmos/gdas.t18z.sfcanl.nc'
      else:
       if var=="sst" or var=="temp":
        fnames00_ana = dpath+'gdas.'+yyyy+mm+dd+'/00/analysis/ocean/gdas.t00z.ocnana.nc'
        fnames06_ana = dpath+'gdas.'+yyyy+mm+dd+'/06/analysis/ocean/gdas.t06z.ocnana.nc'
        fnames12_ana = dpath+'gdas.'+yyyy+mm+dd+'/12/analysis/ocean/gdas.t12z.ocnana.nc'
        fnames18_ana = dpath+'gdas.'+yyyy+mm+dd+'/18/analysis/ocean/gdas.t18z.ocnana.nc'
       elif var=="skint" or var=="atmosT":
        fnames00_ana = dpath+'gdas.'+yyyy+mm+dd+'/00/analysis/atmos/gdas.t00z.sfcanl.nc'
        fnames06_ana = dpath+'gdas.'+yyyy+mm+dd+'/06/analysis/atmos/gdas.t06z.sfcanl.nc'
        fnames12_ana = dpath+'gdas.'+yyyy+mm+dd+'/12/analysis/atmos/gdas.t12z.sfcanl.nc'
        fnames18_ana = dpath+'gdas.'+yyyy+mm+dd+'/18/analysis/atmos/gdas.t18z.sfcanl.nc'

      if np.size(hours)==1:
        if hours=='00':
          tfnames = [fnames00_ana]
        elif hours=='06':
          tfnames = [fnames06_ana]
        elif hours=='12':
          tfnames = [fnames12_ana]
        elif hours=='18':
          tfnames = [fnames18_ana]

      fnamesDAILY = [fnames00_ana]
      fnamesDAILY.append(fnames06_ana)
      fnamesDAILY.append(fnames12_ana)
      fnamesDAILY.append(fnames18_ana)
 
    elif filetype=="inc":		# analysis increments
      if dpath.find("role")!=-1:
       if var=='sst' or var=='temp':
        fnames00_inc = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/analysis/ocean/gdas.t00z.ocninc.nc'
        fnames06_inc = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/analysis/ocean/gdas.t06z.ocninc.nc'
        fnames12_inc = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/analysis/ocean/gdas.t12z.ocninc.nc'
        fnames18_inc = dpath+'/'+yyyy+mm+dd+'18/gdas.'+yyyy+mm+dd+'/18/analysis/ocean/gdas.t18z.ocninc.nc'
       elif var=="skint" or var=="atmosT":
        fnames00_inc = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/analysis/atmos/gdas.t00z.sfcanl.nc'
        fnames06_inc = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/analysis/atmos/gdas.t06z.sfcanl.nc'
        fnames12_inc = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/analysis/atmos/gdas.t12z.sfcanl.nc'
        fnames18_inc = dpath+'/'+yyyy+mm+dd+'18/gdas.'+yyyy+mm+dd+'/18/analysis/atmos/gdas.t18z.sfcanl.nc'
      else:
        fnames00_inc = dpath+'gdas.'+yyyy+mm+dd+'/00/analysis/ocean/gdas.t00z.ocninc.nc'
        fnames06_inc = dpath+'gdas.'+yyyy+mm+dd+'/06/analysis/ocean/gdas.t06z.ocninc.nc'
        fnames12_inc = dpath+'gdas.'+yyyy+mm+dd+'/12/analysis/ocean/gdas.t12z.ocninc.nc'
        fnames18_inc = dpath+'gdas.'+yyyy+mm+dd+'/18/analysis/ocean/gdas.t18z.ocninc.nc'

      if np.size(hours)==1:
        if hours=='00':
          tfnames = [fnames00_inc]
        elif hours=='06':
          tfnames = [fnames06_inc]
        elif hours=='12':
          tfnames = [fnames12_inc]
        elif hours=='18':
          tfnames = [fnames18_inc]

      fnamesDAILY = [fnames00_inc]
      fnamesDAILY.append(fnames06_inc)
      fnamesDAILY.append(fnames12_inc)
      fnamesDAILY.append(fnames18_inc)

    #print("get_dailymean: fnames = "+str(fnamesDAILY))

	# get number of files
    nfiles = np.size(tfnames)

    if nfiles==0:
      print("get_dailymean: NO FILES for "+yyyy+mm+dd+" at "+str(hours))
      sst_sum = []
      sst     = []
      del nfiles, nfilesDAILY, tfnames
      return sst_sum, sst
    del nfiles

	# compute daily mean
    sst_sum = []								# for computing daily mean
    cntDAILY  = 0.0	# daily mean file count
    for i in range(np.size(fnamesDAILY)):
     if exists(str(fnamesDAILY[i]))==True:

      ds = xr.open_dataset(str(fnamesDAILY[i]))

      if var=="sst" or var=="temp":
          t_sst = np.squeeze(ds[varname][0,0,:,:])             # units = degrees-C
      elif var=="skint" or var=="atmosT":
          t_sst = np.squeeze(ds[varname][0,:,:]) - 273.15      # units = degrees-C
          mask = np.squeeze(ds['land'][0,:,:])
          sstmask = np.where(mask>0, np.nan, t_sst)            # mask land and sea ice
          del mask, t_sst

      ds.close()

      if cntDAILY == 0:
        if var=="sst" or var=="temp":
          sst_sum = t_sst
          del t_sst
        elif var=="skint" or var=="atmosT":
          sst_sum = sstmask
          del sstmask
      else:
        if var=="sst" or var=="temp":
          sst_sum += t_sst
          del t_sst
        elif var=="skint" or var=="atmosT":
          sst_sum += sstmask
          del sstmask

      cntDAILY +=1

    if np.size(sst_sum)==0:
      daily_mean = []
    else:
      daily_mean = np.array(sst_sum)/cntDAILY
    del cntDAILY, sst_sum

	# extract SST for specified hour block
    cnt = 0	# sst hour file count
    ilat = 0
    for ih in range(np.size(hours)):			
      fnames = tfnames
      for fname in fnames:
        if exists(str(fname))==True:
          print("get_dailymean: fname = "+str(fname))

          ds = xr.open_dataset(fname)

          if var=="sst":
            tt_sst = np.squeeze(ds[varname][0,0,:,:])             # units = degrees-C
          elif var=="temp":
            tt_sst = np.squeeze(ds[varname][0,:,:,:])             # units = degrees-C
          elif var=="skint" or var=="atmosT":
            tt_sst = np.squeeze(ds[varname][0,:,:]) - 273.15      # units = degrees-C
            mask = np.squeeze(ds['land'][0,:,:])
            tsstmask = np.where(mask>0, np.nan, tt_sst)            # mask land and sea ice
          shape_t = np.shape(tt_sst)

          ds.close()

          print("get_dailymean: min/max w/o nan: sst = "+str(np.nanmin(tt_sst))+" "+str(np.nanmax(tt_sst)))

          if ilat==0:
            if dpath.find("role")!=-1 or dpath.find("marine_candidate")!=-1:
              if dpath.find("marine_candidate")!=-1:
                modelstr = "model"
              else:
                modelstr = "model_data"
              if var=="skint" or var=="atmosT":
                f006  = dpath+'/'+yyyy+mm+dd+hours+'/gdas.'+yyyy+mm+dd+'/'+hours+'/'+modelstr+'/atmos/history/gdas.t'+hours+'z.sfcf006.nc'
              else:
                f006  = dpath+'/'+yyyy+mm+dd+hours+'/gdas.'+yyyy+mm+dd+'/'+hours+'/'+modelstr+'/ocean/history/gdas.ocean.t'+hours+'z.inst.f000.nc'
            else:
              modelstr = "model_data"
              if var=="skint" or var=="atmosT":
                f006  = dpath+'/gdas.'+yyyy+mm+dd+'/'+hours+'/'+modelstr+'/atmos/history/gdas.t'+hours+'z.sfcf006.nc'
              else:
                f006  = dpath+'/gdas.'+yyyy+mm+dd+'/'+hours+'/'+modelstr+'/ocean/history/gdas.ocean.t'+hours+'z.inst.f000.nc'

            ds006 = xr.open_dataset(f006)

            if var=="skint" or var=="atmosT":
              lat = np.squeeze(ds006['grid_yt'][:])
              lon = np.squeeze(ds006['grid_xt'][:])
            else:
              lat = np.squeeze(ds006['yh'][:]) 
              lon = np.squeeze(ds006['xh'][:])

            ds006.close()
            del f006, ds006
          ilat += 1

          if cnt == 0:
            if var=="temp":
              #sst = np.empty([1,shape_t[0],shape_t[1],shape_t[2]],dtype=float)
              #sst[0,:,:,:] = np.array(tt_sst)
              sst = np.array(tt_sst)
            elif var=="sst":
              sst = np.empty([1,shape_t[0],shape_t[1]],dtype=float)
              sst[0,:,:] = tt_sst
            elif var=="skint" or var=="atmosT":
              sst = np.empty([1,shape_t[0],shape_t[1]],dtype=float)
              sst[0,:,:] = np.array(tsstmask)
              del mask, tsstmask
            del tt_sst
          else:
            if var=="temp":
              t_sst = np.array(tt_sst)
              sst = sst + t_sst           # sum instead of append. Average takes place outside of this loop.
              del t_sst
            elif var=="sst":
              t_sst = np.empty([1,shape_t[0],shape_t[1]],dtype=float)	# for looking at individual hour blocks
              t_sst[0,:,:] = tt_sst
              sst = sst + t_sst
              #sst = np.append(sst, np.array(t_sst), axis=0)
              ##sst[int(cnt),:,:] = np.array(t_sst)
              del t_sst
            elif var=="skint" or var=="atmosT":
              sstmask = np.empty([1,shape_t[0],shape_t[1]],dtype=float)     # for looking at individual hour blocks
              sstmask[0,:,:] = tsstmask
              sst = sst + np.array(sstmask)
              #sst = np.append(sst, np.array(sstmask), axis=0)
              ##sst[int(cnt),:,:] = np.array(sstmask)
              del mask, tsstmask, sstmask
            del tt_sst
          del shape_t

          cnt +=1

        elif exists(str(fname))==False:
          print("get_dailymean: "+fname+" does not exist!")
          continue

      # compute average sst for 3D variable (for cross-sections)
      if cnt>0: # and var=="temp": 
        sst = sst/cnt
      #  if filetype=="anl": 
          #f006   = dpath+'gdas.'+yyyy+mm+dd+'/'+hours+'/model_data/ocean/history/gdas.ocean.t'+hours+'z.inst.f006.nc'
          #ds006  = xr.open_dataset(f006)
          #sst006 = np.squeeze(ds[varname][0,:,:,:])
      #    sst = np.where(np.isnan(sst006)==True, np.nan, sst)
      #del sst006

    if cnt==0:
      sst = []
      lat = []
      lon = []

    print("get_dailymean: shape lat = "+str(np.shape(lat))+" lon = "+str(np.shape(lon)))

    return daily_mean, sst, lat, lon
    
#===============================================================================================
# Compute ocean T profiles from each file called
#
#	Return: Data for each file opened
#
def get_ocnprofiles(dpath, yyyy, mm, dd, hours, filetype, var, cgrd, grdstr):

    if var=="temp": varname = 'Temp'
    
	# date before yyyymmdd
            # month
    if int(dd)==1:
      tmmB = int(mm) - 1
    else:
      tmmB = int(mm)
    if tmmB<10:
      mmB = "0"+str(tmmB)
    else:
      mmB = str(tmmB)
            # day
    if mmB==mm:
      tddB = int(dd) - 1
      if tddB<10:
        ddB = "0"+str(tddB)
      else:
        ddB = str(tddB)
    else:
      if mmB=="04" or mmB=="06" or mmB=="09" or mmB=="11":
        ddB = "30"
      elif mmB=="02":
        if int(yyyy)%4==0:
          ddB = "29"
        else:
          ddB = "28"
      else:
        ddB = "31"
            # hour
    HHarr = ['00','06','12','18']
    if hours=='00': hhB = '18'
    if hours=='06': hhB = '00'
    if hours=='12': hhB = '06'
    if hours=='18': hhB = '12'

	# get filenames
    if filetype=="all_fcsts":		# all background forecasts

      if dpath.find("role")!=-1:
        #dpathB = dpath.replace(yyyy+mm+dd+hours, yyyy+mmB+ddB+hhB)
        if var=="sst" or var=="temp":
          # before date
          Bfnames18_f006 = dpath+'/'+yyyy+mmB+ddB+'18/gdas.'+yyyy+mmB+ddB+'/18/model_data/ocean/history/gdas.ocean.t18z.inst.f006.nc'
          Bfnames18_f009 = dpath+'/'+yyyy+mmB+ddB+'18/gdas.'+yyyy+mmB+ddB+'/18/model_data/ocean/history/gdas.ocean.t18z.inst.f009.nc'
          # current date
          fnames00_f003 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.ocean.t00z.inst.f003.nc'
          fnames00_f006 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.ocean.t00z.inst.f006.nc'
          fnames00_f009 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.ocean.t00z.inst.f009.nc'
          fnames06_f003 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.ocean.t06z.inst.f003.nc'
          fnames06_f006 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.ocean.t06z.inst.f006.nc'
          fnames06_f009 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.ocean.t06z.inst.f009.nc'
          fnames12_f003 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.ocean.t12z.inst.f003.nc'
          fnames12_f006 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.ocean.t12z.inst.f006.nc'
          fnames12_f009 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.ocean.t12z.inst.f009.nc'
          fnames18_f003 = dpath+'/'+yyyy+mm+dd+'18/gdas.'+yyyy+mm+dd+'/18/model_data/ocean/history/gdas.ocean.t18z.inst.f003.nc'
        elif var=="skint":
          # before date
          Bfnames18_f006 = dpath+'/'+yyyy+mmB+ddB+'18/gdas.'+yyyy+mmB+ddB+'/18/model_data/atmos/history/gdas.t18z.sfcf006.nc'
          Bfnames18_f009 = dpath+'/'+yyyy+mmB+ddB+'18/gdas.'+yyyy+mmB+ddB+'/18/model_data/atmos/history/gdas.t18z.sfcf009.nc'
          # current date
          fnames00_f003 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/model_data/atmos/history/gdas.t00z.sfcf003.nc'
          fnames00_f006 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/model_data/atmos/history/gdas.t00z.sfcf006.nc'
          fnames00_f009 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/model_data/atmos/history/gdas.t00z.sfcf009.nc'
          fnames06_f003 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/model_data/atmos/history/gdas.t06z.sfcf003.nc'
          fnames06_f006 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/model_data/atmos/history/gdas.t06z.sfcf006.nc'
          fnames06_f009 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/model_data/atmos/history/gdas.t06z.sfcf009.nc'
          fnames12_f003 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/model_data/atmos/history/gdas.t12z.sfcf003.nc'
          fnames12_f006 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/model_data/atmos/history/gdas.t12z.sfcf006.nc'
          fnames12_f009 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/model_data/atmos/history/gdas.t12z.sfcf009.nc'
          fnames18_f003 = dpath+'/'+yyyy+mm+dd+'18/gdas.'+yyyy+mm+dd+'/18/model_data/atmos/history/gdas.t18z.sfcf003.nc'
      else: 
       if var=="sst" or var=="temp":
          # before date
        Bfnames18_f006 = dpath+'gdas.'+yyyy+mmB+ddB+'/18/model_data/ocean/history/gdas.ocean.t18z.inst.f006.nc'
        Bfnames18_f009 = dpath+'gdas.'+yyyy+mmB+ddB+'/18/model_data/ocean/history/gdas.ocean.t18z.inst.f009.nc'
	  # current date
        fnames00_f003 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.ocean.t00z.inst.f003.nc'
        fnames00_f006 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.ocean.t00z.inst.f006.nc'
        fnames00_f009 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.ocean.t00z.inst.f009.nc'
        fnames06_f003 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.ocean.t06z.inst.f003.nc'
        fnames06_f006 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.ocean.t06z.inst.f006.nc'
        fnames06_f009 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.ocean.t06z.inst.f009.nc'
        fnames12_f003 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.ocean.t12z.inst.f003.nc'
        fnames12_f006 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.ocean.t12z.inst.f006.nc'
        fnames12_f009 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.ocean.t12z.inst.f009.nc'
        fnames18_f003 = dpath+'gdas.'+yyyy+mm+dd+'/18/model_data/ocean/history/gdas.ocean.t18z.inst.f003.nc'
       elif var=="skint":
          # before date
        Bfnames18_f006 = dpath+'gdas.'+yyyy+mmB+ddB+'/18/model_data/atmos/history/gdas.t18z.sfcf006.nc'
        Bfnames18_f009 = dpath+'gdas.'+yyyy+mmB+ddB+'/18/model_data/atmos/history/gdas.t18z.sfcf009.nc'
          # current date
        fnames00_f003 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/atmos/history/gdas.t00z.sfcf003.nc'
        fnames00_f006 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/atmos/history/gdas.t00z.sfcf006.nc'
        fnames00_f009 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/atmos/history/gdas.t00z.sfcf009.nc'
        fnames06_f003 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/atmos/history/gdas.t06z.sfcf003.nc'
        fnames06_f006 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/atmos/history/gdas.t06z.sfcf006.nc'
        fnames06_f009 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/atmos/history/gdas.t06z.sfcf009.nc'
        fnames12_f003 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/atmos/history/gdas.t12z.sfcf003.nc'
        fnames12_f006 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/atmos/history/gdas.t12z.sfcf006.nc'
        fnames12_f009 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/atmos/history/gdas.t12z.sfcf009.nc'
        fnames18_f003 = dpath+'gdas.'+yyyy+mm+dd+'/18/model_data/atmos/history/gdas.t18z.sfcf003.nc'

      if np.size(hours)==1:
        if hours=='00':
          tfnames = [Bfnames18_f006]
          tfnames.append(Bfnames18_f009)
          tfnames.append(fnames00_f003)
        elif hours=='06':
          tfnames = [fnames00_f006]
          tfnames.append(fnames00_f009)
          tfnames.append(fnames06_f003)
        elif hours=='12':
          tfnames = [fnames06_f006]
          tfnames.append(fnames06_f009)
          tfnames.append(fnames12_f003)
        elif hours=='18':
          tfnames = [fnames12_f006]
          tfnames.append(fnames12_f009)
          tfnames.append(fnames18_f003)

      del Bfnames18_f006, Bfnames18_f009
      del fnames00_f003, fnames00_f006, fnames00_f009
      del fnames06_f003, fnames06_f006, fnames06_f009
      del fnames12_f003, fnames12_f006, fnames12_f009
      del fnames18_f003

    elif filetype=="f06":		# only 6-h background forecasts
      if dpath.find("role")!=-1:
       if var=="sst" or var=="temp":
          # before date
        Bfnames18_f006 = dpath+'/'+yyyy+mmB+ddB+'18/gdas.'+yyyy+mmB+ddB+'/18/model_data/ocean/history/gdas.ocean.t18z.inst.f006.nc'
	  # current date
        fnames00_f006 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.ocean.t00z.inst.f006.nc'
        fnames06_f006 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.ocean.t06z.inst.f006.nc'
        fnames12_f006 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.ocean.t12z.inst.f006.nc'
       elif var=="skint":
          # before date
        Bfnames18_f006 = dpath+'/'+yyyy+mmB+ddB+'18/gdas.'+yyyy+mmB+ddB+'/18/model_data/atmos/history/gdas.t18z.sfcf006.nc'
          # current date
        fnames00_f006 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/model_data/atmos/history/gdas.t00z.sfcf006.nc'
        fnames06_f006 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/model_data/atmos/history/gdas.t06z.sfcf006.nc'
        fnames12_f006 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/model_data/atmos/history/gdas.t12z.sfcf006.nc'
      else:
       if var=="sst" or var=="temp":
          # before date
        Bfnames18_f006 = dpath+'gdas.'+yyyy+mmB+ddB+'/18/model_data/ocean/history/gdas.ocean.t18z.inst.f006.nc'
	  # current date
        fnames00_f006 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.ocean.t00z.inst.f006.nc'
        fnames06_f006 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.ocean.t06z.inst.f006.nc'
        fnames12_f006 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.ocean.t12z.inst.f006.nc'
       elif var=="skint":
          # before date
        Bfnames18_f006 = dpath+'gdas.'+yyyy+mmB+ddB+'/18/model_data/atmos/history/gdas.t18z.sfcf006.nc'
          # current date
        fnames00_f006 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/atmos/history/gdas.t00z.sfcf006.nc'
        fnames06_f006 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/atmos/history/gdas.t06z.sfcf006.nc'
        fnames12_f006 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/atmos/history/gdas.t12z.sfcf006.nc'

      if np.size(hours)==1:
        if hours=='00':
          tfnames = [Bfnames18_f006]
        elif hours=='06':
          tfnames = [fnames00_f006]
        elif hours=='12':
          tfnames = [fnames06_f006]
        elif hours=='18':
          tfnames = [fnames12_f006]

    elif filetype=="anl":               # analysis files
      if dpath.find("role")!=-1:
       if var=="sst" or var=="temp":
        fnames00_ana = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/analysis/ocean/gdas.t00z.ocnana.nc'
        fnames06_ana = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/analysis/ocean/gdas.t06z.ocnana.nc'
        fnames12_ana = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/analysis/ocean/gdas.t12z.ocnana.nc'
        fnames18_ana = dpath+'/'+yyyy+mm+dd+'18/gdas.'+yyyy+mm+dd+'/18/analysis/ocean/gdas.t18z.ocnana.nc'
       elif var=="skint":
        fnames00_ana = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/analysis/atmos/gdas.t00z.sfcanl.nc'
        fnames06_ana = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/analysis/atmos/gdas.t06z.sfcanl.nc'
        fnames12_ana = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/analysis/atmos/gdas.t12z.sfcanl.nc'
        fnames18_ana = dpath+'/'+yyyy+mm+dd+'18/gdas.'+yyyy+mm+dd+'/18/analysis/atmos/gdas.t18z.sfcanl.nc'
      else:
       if var=="sst" or var=="temp":
        fnames00_ana = dpath+'gdas.'+yyyy+mm+dd+'/00/analysis/ocean/gdas.t00z.ocnana.nc'
        fnames06_ana = dpath+'gdas.'+yyyy+mm+dd+'/06/analysis/ocean/gdas.t06z.ocnana.nc'
        fnames12_ana = dpath+'gdas.'+yyyy+mm+dd+'/12/analysis/ocean/gdas.t12z.ocnana.nc'
        fnames18_ana = dpath+'gdas.'+yyyy+mm+dd+'/18/analysis/ocean/gdas.t18z.ocnana.nc'
       elif var=="skint":
        fnames00_ana = dpath+'gdas.'+yyyy+mm+dd+'/00/analysis/atmos/gdas.t00z.sfcanl.nc'
        fnames06_ana = dpath+'gdas.'+yyyy+mm+dd+'/06/analysis/atmos/gdas.t06z.sfcanl.nc'
        fnames12_ana = dpath+'gdas.'+yyyy+mm+dd+'/12/analysis/atmos/gdas.t12z.sfcanl.nc'
        fnames18_ana = dpath+'gdas.'+yyyy+mm+dd+'/18/analysis/atmos/gdas.t18z.sfcanl.nc'

      if np.size(hours)==1:
        if hours=='00':
          tfnames = [fnames00_ana]
        elif hours=='06':
          tfnames = [fnames06_ana]
        elif hours=='12':
          tfnames = [fnames12_ana]
        elif hours=='18':
          tfnames = [fnames18_ana]

    elif filetype=="inc":		# analysis increments
      if dpath.find("role")!=-1:
        fnames00_inc = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/analysis/ocean/gdas.t00z.ocninc.nc'
        fnames06_inc = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/analysis/ocean/gdas.t06z.ocninc.nc'
        fnames12_inc = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/analysis/ocean/gdas.t12z.ocninc.nc'
        fnames18_inc = dpath+'/'+yyyy+mm+dd+'18/gdas.'+yyyy+mm+dd+'/18/analysis/ocean/gdas.t18z.ocninc.nc'
      else:
        fnames00_inc = dpath+'gdas.'+yyyy+mm+dd+'/00/analysis/ocean/gdas.t00z.ocninc.nc'
        fnames06_inc = dpath+'gdas.'+yyyy+mm+dd+'/06/analysis/ocean/gdas.t06z.ocninc.nc'
        fnames12_inc = dpath+'gdas.'+yyyy+mm+dd+'/12/analysis/ocean/gdas.t12z.ocninc.nc'
        fnames18_inc = dpath+'gdas.'+yyyy+mm+dd+'/18/analysis/ocean/gdas.t18z.ocninc.nc'

      if np.size(hours)==1:
        if hours=='00':
          tfnames = [fnames00_inc]
        elif hours=='06':
          tfnames = [fnames06_inc]
        elif hours=='12':
          tfnames = [fnames12_inc]
        elif hours=='18':
          tfnames = [fnames18_inc]

	# get number of files
    nfiles = np.size(tfnames)

    if nfiles==0:
      print("get_ocnprofiles: NO FILES for "+yyyy+mm+dd+" at "+str(hours))
      sst_sum = []
      sst     = []
      del nfiles, nfilesDAILY, tfnames
      return sst_sum, sst
    del nfiles

	# extract SST for specified hour block
    cnt = 0	# sst hour file count
    ilat = 0
    for ih in range(np.size(hours)):			
      fnames = tfnames
      for fname in fnames:
        if exists(str(fname))==True:
          print("get_ocnprofiles: fname = "+str(fname))

          ds = xr.open_dataset(fname)

          if var=="temp":
            tt_sst = np.squeeze(ds[varname][0,:,:,:])             # units = degrees-C
          #shape_t = np.shape(tt_sst)

          ds.close()

          print("get_ocnprofiles: min/max w/o nan: sst = "+str(np.nanmin(tt_sst))+" "+str(np.nanmax(tt_sst)))

          # get lat/lon arrays, and find array space that is closest to cross-section lat or lon
          if ilat==0:
            if dpath.find("role")!=-1:
              f006  = dpath+'/'+yyyy+mm+dd+hours+'/gdas.'+yyyy+mm+dd+'/'+hours+'/model_data/ocean/history/gdas.ocean.t'+hours+'z.inst.f000.nc'
            else:
              f006  = dpath+'/gdas.'+yyyy+mm+dd+'/'+hours+'/model_data/ocean/history/gdas.ocean.t'+hours+'z.inst.f000.nc'
            print("get_ocnprofiles: f006 for lat/lon = "+str(f006))
            ds006 = xr.open_dataset(f006)
            lat = np.squeeze(ds006['yh'][:]) 
            lon = np.squeeze(ds006['xh'][:])
            #sst006 = np.squeeze(ds[varname][0,:,:,:])
            ds006.close()
            del f006, ds006

	    # get array space equating to cross-section lat or lon
            if grdstr=="lat" or grdstr=="LAT":
              flats = np.array(lat)
              # find min latitude difference to crosslat
              tdiff = np.around(np.array(abs(flats - cgrd)), 3)
              tmin  = np.nanmin(tdiff)
              if np.size(tmin)>1:
                tminmin = tmin[0]
              else:
                tminmin = tmin
              for i in range(np.size(tdiff)):
                if tdiff[i]==tminmin: imin = i
              del tdiff, tmin, tminmin, flats
            elif grdstr=="lon" or grdstr=="LON":
              flons = np.array(lon)
              # find min longitude difference to crosslon
              tflons = np.where(flons<=-180, flons + 360, flons)        # convert ocean model lons [-300, 60] to [-180, 180]
              tdiff = np.around(np.array(abs(tflons - cgrd)), 3)
              tmin  = np.min(tdiff)
              if np.size(tmin)>1:
                tminmin = tmin[0]
              else:
                tminmin = tmin
              for i in range(np.size(tdiff)):
                if tdiff[i]==tminmin: imin = i
              del tdiff, tmin, tminmin, tflons, flons
          ilat += 1

          if grdstr=="lat" or grdstr=="LAT":
            vert0 = tt_sst[:,imin,:]
          elif grdstr=="lon" or grdstr=="LON":
            vert0 = tt_sst[:,:,imin]
          del tt_sst

	  # get T profiles at specified cross-section lat or lon
          if cnt == 0:
            if var=="temp":
              sst = np.array(vert0)
            del vert0
          else:
            if var=="temp":
              t_sst = np.array(vert0)
              sst = sst + t_sst           # sum instead of append. Average takes place outside of this loop.
              del t_sst
            del vert0

          cnt +=1

        elif exists(str(fname))==False:
          print("get_ocnprofiles: "+fname+" does not exist!")
          #continue

      # compute average sst for 3D variable (for cross-sections)
      if cnt>0: # and var=="temp": 
        sst = sst/cnt

    if cnt==0:
      sst = []
      lat = []
      lon = []
    print("get_ocnprofiles: cnt = "+str(cnt))

    #if grdstr=="lat" or grdstr=="LAT":
    #  grd = lat
    #elif grdstr=="lon" or grdstr=="LON":
    #  grd = lon

    return sst, lat, lon
    
#===============================================================================================
# Get ocean model background or analysis files, Compute daily mean centered around 12z
#       For comparisons with OSTIA (which is centered around 12z)
#
#	Only for Foundation T
#
#	Return: Daily mean
#
def get_dailymean_12zcenter(dpath, yyyy, mm, dd, filetype, fhour, var):

    if var=="temp": varname = 'Temp'            # ocean temperature, all levels
    if var=="sst": varname = 'Temp'             # coean temperature, top layer (SST)
    if var=="skint" or var=="atmosT": varname = 'tmpsfc'         # ocean temperature, surface (from atmos)

    hours = ['12']
    nhours = np.size(hours)

    hoursDAY = ['00','06','12','18']		# used for computing daily mean

	# date before yyyymmdd
            # month
    if int(dd)==1:
      tmmB = int(mm) - 1
    else:
      tmmB = int(mm)
    if tmmB<10:
      mmB = "0"+str(tmmB)
    else:
      mmB = str(tmmB)
            # day
    if mmB==mm:
      tddB = int(dd) - 1
      if tddB<10:
        ddB = "0"+str(tddB)
      else:
        ddB = str(tddB)
    else:
      if mmB=="04" or mmB=="06" or mmB=="09" or mmB=="11":
        ddB = "30"
      elif mmB=="02":
        if int(yyyy)%4==0:
          ddB = "29"
        else:
          ddB = "28"
      else:
        ddB = "31"

	# get filenames
    if filetype=="all_fcsts":		# all background forecasts
      fill = 0
      latname = "yh"
      lonname = "xh"
      if dpath.find("role")!=-1:
        # GFSv17 prototype runs
        if var=="sst" or var=="temp":
          # before date
          Bfnames18_f006 = dpath+'/'+yyyy+mmB+ddB+'18/gdas.'+yyyy+mmB+ddB+'/18/model_data/ocean/history/gdas.ocean.t18z.inst.f006.nc'
          Bfnames18_f009 = dpath+'/'+yyyy+mmB+ddB+'18/gdas.'+yyyy+mmB+ddB+'/18/model_data/ocean/history/gdas.ocean.t18z.inst.f009.nc'
          # current date
          fnames00_f003 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.ocean.t00z.inst.f003.nc'
          fnames00_f006 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.ocean.t00z.inst.f006.nc'
          fnames00_f009 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.ocean.t00z.inst.f009.nc'
          fnames06_f003 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.ocean.t06z.inst.f003.nc'
          fnames06_f006 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.ocean.t06z.inst.f006.nc'
          fnames06_f009 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.ocean.t06z.inst.f009.nc'
          fnames12_f003 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.ocean.t12z.inst.f003.nc'
          fnames12_f006 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.ocean.t12z.inst.f006.nc'
          fnames12_f009 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.ocean.t12z.inst.f009.nc'
          fnames18_f003 = dpath+'/'+yyyy+mm+dd+'18/gdas.'+yyyy+mm+dd+'/18/model_data/ocean/history/gdas.ocean.t18z.inst.f003.nc'
        elif var=="skint" or var=="atmosT":
          # before date
          Bfnames18_f006 = dpath+'/'+yyyy+mmB+ddB+'18/gdas.'+yyyy+mmB+ddB+'/18/model_data/atmos/history/gdas.t18z.sfcf006.nc'
          Bfnames18_f009 = dpath+'/'+yyyy+mmB+ddB+'18/gdas.'+yyyy+mmB+ddB+'/18/model_data/atmos/history/gdas.t18z.sfcf009.nc'
          # current date
          fnames00_f003 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/model_data/atmos/history/gdas.t00z.sfcf003.nc'
          fnames00_f006 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/model_data/atmos/history/gdas.t00z.sfcf006.nc'
          fnames00_f009 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/model_data/atmos/history/gdas.t00z.sfcf009.nc'
          fnames06_f003 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/model_data/atmos/history/gdas.t06z.sfcf003.nc'
          fnames06_f006 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/model_data/atmos/history/gdas.t06z.sfcf006.nc'
          fnames06_f009 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/model_data/atmos/history/gdas.t06z.sfcf009.nc'
          fnames12_f003 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/model_data/atmos/history/gdas.t12z.sfcf003.nc'
          fnames12_f006 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/model_data/atmos/history/gdas.t12z.sfcf006.nc'
          fnames12_f009 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/model_data/atmos/history/gdas.t12z.sfcf009.nc'
          fnames18_f003 = dpath+'/'+yyyy+mm+dd+'18/gdas.'+yyyy+mm+dd+'/18/model_data/atmos/history/gdas.t18z.sfcf003.nc'
      else:
        # Katie's new runs
                # before date
        Bfnames18_f006 = dpath+'gdas.'+yyyy+mmB+ddB+'/18/model_data/ocean/history/gdas.ocean.t18z.inst.f006.nc'
        Bfnames18_f009 = dpath+'gdas.'+yyyy+mmB+ddB+'/18/model_data/ocean/history/gdas.ocean.t18z.inst.f009.nc'
		# current date
        fnames00_f003 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.ocean.t00z.inst.f003.nc'
        fnames00_f006 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.ocean.t00z.inst.f006.nc'
        fnames00_f009 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.ocean.t00z.inst.f009.nc'
        fnames06_f003 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.ocean.t06z.inst.f003.nc'
        fnames06_f006 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.ocean.t06z.inst.f006.nc'
        fnames06_f009 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.ocean.t06z.inst.f009.nc'
        fnames12_f003 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.ocean.t12z.inst.f003.nc'
        fnames12_f006 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.ocean.t12z.inst.f006.nc'
        fnames12_f009 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.ocean.t12z.inst.f009.nc'
        fnames18_f003 = dpath+'gdas.'+yyyy+mm+dd+'/18/model_data/ocean/history/gdas.ocean.t18z.inst.f003.nc'

      if dpath.find("ocn-da")!=-1:
        # Katie's old runs
                # before date
        Bfnames18_f006 = dpath+'gdas.'+yyyy+mmB+ddB+'/18/model_data/ocean/history/gdas.t18z.ocnf006.nc'
        Bfnames18_f009 = dpath+'gdas.'+yyyy+mmB+ddB+'/18/model_data/ocean/history/gdas.t18z.ocnf009.nc'
                # current date
        fnames00_f003 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.t00z.ocnf003.nc'
        fnames00_f006 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.t00z.ocnf006.nc'
        fnames00_f009 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.t00z.ocnf009.nc'
        fnames06_f003 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.t06z.ocnf003.nc'
        fnames06_f006 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.t06z.ocnf006.nc'
        fnames06_f009 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.t06z.ocnf009.nc'
        fnames12_f003 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.t12z.ocnf003.nc'
        fnames12_f006 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.t12z.ocnf006.nc'
        fnames12_f009 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.t12z.ocnf009.nc'
        fnames18_f003 = dpath+'gdas.'+yyyy+mm+dd+'/18/model_data/ocean/history/gdas.t18z.ocnf003.nc'

      fnamesDAILY = [Bfnames18_f006]
      fnamesDAILY.append(Bfnames18_f009)
      fnamesDAILY.append(fnames00_f003)
      fnamesDAILY.append(fnames00_f006)
      fnamesDAILY.append(fnames00_f009)
      fnamesDAILY.append(fnames06_f003)
      fnamesDAILY.append(fnames06_f006)
      fnamesDAILY.append(fnames06_f009)
      fnamesDAILY.append(fnames12_f003)
      fnamesDAILY.append(fnames12_f006)
      fnamesDAILY.append(fnames12_f009)
      fnamesDAILY.append(fnames18_f003)
      
    elif filetype=="f06":		# 6-h background forecasts only
      fill = 0
      latname = "yh"
      lonname = "xh"
      if dpath.find("role")!=-1:
      		# before date
        Bfnames18_f006 = dpath+'/'+yyyy+mmB+ddB+'18/gdas.'+yyyy+mmB+ddB+'/18/model_data/ocean/history/gdas.ocean.t18z.inst.f006.nc'
		# current date
        fnames00_f006 = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.ocean.t00z.inst.f006.nc'
        fnames06_f006 = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.ocean.t06z.inst.f006.nc'
        fnames12_f006 = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.ocean.t12z.inst.f006.nc'
      else:
		# before date
        Bfnames18_f006 = dpath+'gdas.'+yyyy+mmB+ddB+'/18/model_data/ocean/history/gdas.ocean.t18z.inst.f006.nc'
		# current date
        fnames00_f006 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.ocean.t00z.inst.f006.nc'
        fnames06_f006 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.ocean.t06z.inst.f006.nc'
        fnames12_f006 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.ocean.t12z.inst.f006.nc'

      fnamesDAILY = [Bfnames18_f006]
      fnamesDAILY.append(fnames00_f006)
      fnamesDAILY.append(fnames06_f006)
      fnamesDAILY.append(fnames12_f006)

    elif filetype=="anl":               # analyses
      fill = 0
      latname = "yaxis_1"
      lonname = "xaxis_1"
      if dpath.find("role")!=-1:
       if var=="sst" or var=="temp":
        fnames00_ana = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/analysis/ocean/gdas.t00z.ocnana.nc'
        fnames06_ana = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/analysis/ocean/gdas.t06z.ocnana.nc'
        fnames12_ana = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/analysis/ocean/gdas.t12z.ocnana.nc'
        fnames18_ana = dpath+'/'+yyyy+mm+dd+'18/gdas.'+yyyy+mm+dd+'/18/analysis/ocean/gdas.t18z.ocnana.nc'
       elif var=="skint" or var=="atmosT":
        fnames00_ana = dpath+'/'+yyyy+mm+dd+'00/gdas.'+yyyy+mm+dd+'/00/analysis/atmos/gdas.t00z.sfcanl.nc'
        fnames06_ana = dpath+'/'+yyyy+mm+dd+'06/gdas.'+yyyy+mm+dd+'/06/analysis/atmos/gdas.t06z.sfcanl.nc'
        fnames12_ana = dpath+'/'+yyyy+mm+dd+'12/gdas.'+yyyy+mm+dd+'/12/analysis/atmos/gdas.t12z.sfcanl.nc'
        fnames18_ana = dpath+'/'+yyyy+mm+dd+'18/gdas.'+yyyy+mm+dd+'/18/analysis/atmos/gdas.t18z.sfcanl.nc'
      else:
        fnames00_ana = dpath+'gdas.'+yyyy+mm+dd+'/00/analysis/ocean/gdas.t00z.ocnana.nc'
        fnames06_ana = dpath+'gdas.'+yyyy+mm+dd+'/06/analysis/ocean/gdas.t06z.ocnana.nc'
        fnames12_ana = dpath+'gdas.'+yyyy+mm+dd+'/12/analysis/ocean/gdas.t12z.ocnana.nc'
        fnames18_ana = dpath+'gdas.'+yyyy+mm+dd+'/18/analysis/ocean/gdas.t18z.ocnana.nc'

      fnamesDAILY = [fnames00_ana]
      fnamesDAILY.append(fnames06_ana)
      fnamesDAILY.append(fnames12_ana)
      fnamesDAILY.append(fnames18_ana)

    elif filetype=="fcst_leads":               # forecast long-lead times (day 1, day 2, etc)
      fill = 0
      #varname = 'sst'      # should not use 'daily' because it only has SST, not temperature profile which we need here
      #fnames00_fA = dpath+'gfs.'+yyyy+mm+dd+'/00/model_data/ocean/history/gfs.ocean.t00z.daily.f'+str(fhour)+'.nc'
      #fnamesDAILY = [fnames00_fA]

      if dpath.find("atm_model")!=-1:
          # AtmModel AtmDA run
        varname = "tmpsfc"
        latname = "grid_yt"
        lonname = "grid_xt"
        if int(fhour)/24==1:  fA='006'; fB='012'; fC='018'; fD='024'; fE='003'; fF='009'; fG='015'; fH='021'
        if int(fhour)/24==2:  fA='030'; fB='036'; fC='042'; fD='048'; fE='027'; fF='033'; fG='039'; fH='045'
        if int(fhour)/24==3:  fA='054'; fB='060'; fC='066'; fD='072'; fE='051'; fF='057'; fG='063'; fH='069'
        if int(fhour)/24==4:  fA='078'; fB='084'; fC='090'; fD='096'; fE='075'; fF='081'; fG='087'; fH='093'
        if int(fhour)/24==5:  fA='102'; fB='108'; fC='114'; fD='120'; fE='099'; fF='105'; fG='111'; fH='117'
        if int(fhour)/24==6:  fA='126'; fB='132'; fC='138'; fD='144'; fE='123'; fF='129'; fG='135'; fH='141'
        if int(fhour)/24==7:  fA='150'; fB='156'; fC='162'; fD='168'; fE='147'; fF='153'; fG='159'; fH='165'
        if int(fhour)/24==8:  fA='174'; fB='180'; fC='186'; fD='192'; fE='171'; fF='177'; fG='183'; fH='189'
        if int(fhour)/24==9:  fA='198'; fB='204'; fC='210'; fD='216'; fE='195'; fF='201'; fG='207'; fH='213'
        if int(fhour)/24==10: fA='222'; fB='228'; fC='234'; fD='240'; fE='219'; fF='225'; fG='231'; fH='237'
        print("get_dailymean_12z: fhour = "+str(fhour)+": fA = "+str(fA)+" fB = "+str(fB)+" fC = "+str(fC)+" fD = "+str(fD)+" fE = "+str(fE)+" fF = "+str(fF)+" fG = "+str(fG)+" fH = "+str(fH))
	
        if dpath.find("role")!=-1:
          fnames_fA = dpath+'/'+yyyy+mm+dd+'00/gfs.'+yyyy+mm+dd+'/00/model_data/atmos/history/gfs.t00z.sfcf'+fA+'.nc'
          fnames_fB = dpath+'/'+yyyy+mm+dd+'00/gfs.'+yyyy+mm+dd+'/00/model_data/atmos/history/gfs.t00z.sfcf'+fB+'.nc'
          fnames_fC = dpath+'/'+yyyy+mm+dd+'00/gfs.'+yyyy+mm+dd+'/00/model_data/atmos/history/gfs.t00z.sfcf'+fC+'.nc'
          fnames_fD = dpath+'/'+yyyy+mm+dd+'00/gfs.'+yyyy+mm+dd+'/00/model_data/atmos/history/gfs.t00z.sfcf'+fD+'.nc'
          fnames_fE = dpath+'/'+yyyy+mm+dd+'00/gfs.'+yyyy+mm+dd+'/00/model_data/atmos/history/gfs.t00z.sfcf'+fE+'.nc'
          fnames_fF = dpath+'/'+yyyy+mm+dd+'00/gfs.'+yyyy+mm+dd+'/00/model_data/atmos/history/gfs.t00z.sfcf'+fF+'.nc'
          fnames_fG = dpath+'/'+yyyy+mm+dd+'00/gfs.'+yyyy+mm+dd+'/00/model_data/atmos/history/gfs.t00z.sfcf'+fG+'.nc'
          fnames_fH = dpath+'/'+yyyy+mm+dd+'00/gfs.'+yyyy+mm+dd+'/00/model_data/atmos/history/gfs.t00z.sfcf'+fH+'.nc'
        else:
          fnames_fA = dpath+'gfs.'+yyyy+mm+dd+'/00/model_data/atmos/history/gfs.t00z.sfcf'+fA+'.nc'
          fnames_fB = dpath+'gfs.'+yyyy+mm+dd+'/00/model_data/atmos/history/gfs.t00z.sfcf'+fB+'.nc'
          fnames_fC = dpath+'gfs.'+yyyy+mm+dd+'/00/model_data/atmos/history/gfs.t00z.sfcf'+fC+'.nc'
          fnames_fD = dpath+'gfs.'+yyyy+mm+dd+'/00/model_data/atmos/history/gfs.t00z.sfcf'+fD+'.nc'
          fnames_fE = dpath+'gfs.'+yyyy+mm+dd+'/00/model_data/atmos/history/gfs.t00z.sfcf'+fE+'.nc'
          fnames_fF = dpath+'gfs.'+yyyy+mm+dd+'/00/model_data/atmos/history/gfs.t00z.sfcf'+fF+'.nc'
          fnames_fG = dpath+'gfs.'+yyyy+mm+dd+'/00/model_data/atmos/history/gfs.t00z.sfcf'+fG+'.nc'
          fnames_fH = dpath+'gfs.'+yyyy+mm+dd+'/00/model_data/atmos/history/gfs.t00z.sfcf'+fH+'.nc'

        fnamesDAILY = [fnames_fA]
        fnamesDAILY.append(fnames_fB)
        fnamesDAILY.append(fnames_fC)
        fnamesDAILY.append(fnames_fD)
        fnamesDAILY.append(fnames_fE)
        fnamesDAILY.append(fnames_fF)
        fnamesDAILY.append(fnames_fG)
        fnamesDAILY.append(fnames_fH)

      else:
          # S2SModel runs
        varname = 'temp'
        latname = 'yh'
        lonname = 'xh'
        if int(fhour)/24==1:  fA='006'; fB='012'; fC='018'; fD='024'
        if int(fhour)/24==2:  fA='030'; fB='036'; fC='042'; fD='048'
        if int(fhour)/24==3:  fA='054'; fB='060'; fC='066'; fD='072'
        if int(fhour)/24==4:  fA='078'; fB='084'; fC='090'; fD='096'
        if int(fhour)/24==5:  fA='102'; fB='108'; fC='114'; fD='120'
        if int(fhour)/24==6:  fA='126'; fB='132'; fC='138'; fD='144'
        if int(fhour)/24==7:  fA='150'; fB='156'; fC='162'; fD='168'
        if int(fhour)/24==8:  fA='174'; fB='180'; fC='186'; fD='192'
        if int(fhour)/24==9:  fA='198'; fB='204'; fC='210'; fD='216'
        if int(fhour)/24==10: fA='222'; fB='228'; fC='234'; fD='240'
        #print("get_dailymean_12z: fhour = "+str(fhour)+": fA = "+str(fA)+" fB = "+str(fB)+" fC = "+str(fC)+" fD = "+str(fD))
	
        if dpath.find("role")!=-1:
          fnames_fA = dpath+'/'+yyyy+mm+dd+'00/gfs.'+yyyy+mm+dd+'/00/model_data/ocean/history/gfs.ocean.t00z.6hr_avg.f'+fA+'.nc'
          fnames_fB = dpath+'/'+yyyy+mm+dd+'00/gfs.'+yyyy+mm+dd+'/00/model_data/ocean/history/gfs.ocean.t00z.6hr_avg.f'+fB+'.nc'
          fnames_fC = dpath+'/'+yyyy+mm+dd+'00/gfs.'+yyyy+mm+dd+'/00/model_data/ocean/history/gfs.ocean.t00z.6hr_avg.f'+fC+'.nc'
          fnames_fD = dpath+'/'+yyyy+mm+dd+'00/gfs.'+yyyy+mm+dd+'/00/model_data/ocean/history/gfs.ocean.t00z.6hr_avg.f'+fD+'.nc'
        else:
          fnames_fA = dpath+'gfs.'+yyyy+mm+dd+'/00/model_data/ocean/history/gfs.ocean.t00z.6hr_avg.f'+fA+'.nc'
          fnames_fB = dpath+'gfs.'+yyyy+mm+dd+'/00/model_data/ocean/history/gfs.ocean.t00z.6hr_avg.f'+fB+'.nc'
          fnames_fC = dpath+'gfs.'+yyyy+mm+dd+'/00/model_data/ocean/history/gfs.ocean.t00z.6hr_avg.f'+fC+'.nc'
          fnames_fD = dpath+'gfs.'+yyyy+mm+dd+'/00/model_data/ocean/history/gfs.ocean.t00z.6hr_avg.f'+fD+'.nc'

        fnamesDAILY = [fnames_fA]
        fnamesDAILY.append(fnames_fB)
        fnamesDAILY.append(fnames_fC)
        fnamesDAILY.append(fnames_fD)

    #print("get_dailymean_12z: fnames = "+str(fnamesDAILY))

	# get number of files
    nfilesDAILY = np.size(fnamesDAILY)

    if nfilesDAILY==0:
      print("get_dailymean_12z: NO FILES for "+yyyy+mm+dd+" at "+str(hours))
      sst_sum = []
      sst     = []
      return sst_sum, sst

	# compute daily mean of SST (forecasts) centered around 12z (to match OSTIA)
    sst_sum = []								# for computing daily mean
    cntDAILY  = 0.0	# daily mean file count
    ilat=0
    for i in range(np.size(fnamesDAILY)):
      if exists(str(fnamesDAILY[i]))==True:
        #print("get_dailymean_12z: fname = "+str(fnamesDAILY[i]))
        ds = xr.open_dataset(str(fnamesDAILY[i]))

        if ilat==0:
          latf = np.squeeze(ds[latname][:])
          lonf = np.squeeze(ds[lonname][:])
          ilat = 1

        if dpath.find("atm_model")!=-1:
          fcst  = np.squeeze(ds[varname][0,:,:]) - 273.15
          fmask = np.squeeze(ds['land'][:,:])
          t_sst = np.where(fmask>0, np.nan, fcst)      # mask=0 ocean; mask=1 land; mask=2 sea ice
          del fmask, fcst
          s_sst = np.where(t_sst==fill, np.nan, t_sst)
        else:
#          t_sst = np.squeeze(ds[varname][0,1,:,:])             # 2nd layer of ocean. units = degrees-C
          if var=="sst" or var=="temp":
            t_sst = np.squeeze(ds[varname][0,1,:,:])             # units = degrees-C
            s_sst = np.where(t_sst==fill, np.nan, t_sst)
          elif var=="skint" or var=="atmosT":
            t_sst = np.squeeze(ds[varname][0,:,:]) - 273.15      # units = degrees-C
            mask = np.squeeze(ds['land'][0,:,:])
            sstmask = np.where(mask>0, np.nan, t_sst)            # mask land and sea ice
            del mask, t_sst
            s_sst = np.where(sstmask==fill, np.nan, sstmask)
            del sstmask
        
        if cntDAILY == 0:
          sst_sum = s_sst
        else:
          sst_sum += s_sst
        del t_sst, s_sst
        #print("get_dailymean_12z: sst_sum shape = "+str(np.shape(sst_sum)))

        cntDAILY += 1
        ds.close()
      else:
        continue

    if np.size(sst_sum)==0:
      daily_mean = []
      latf = []
      lonf = []
    else:
      daily_mean = np.array(sst_sum)/cntDAILY

    return daily_mean, latf, lonf
    
#===============================================================================================
# Extract gfs forecast and corresponding verification data
#
#	Return: Daily mean, Data for each file opened, approx. sun position, lat and lon arrays
#
def read_fcst_verif(dpath, yyyy, mm, dd, var, tfhour):

    print("read_fcst_verif: FORECAST DATE = "+str(tfhour)+"-hour forecast made on "+str(yyyy)+str(mm)+str(dd)+" at 00z")

    fhour = tfhour

    # forecast names
    if var=="sst": 
            # ocean sea surface temperature
      if dpath.find('s2s_da')!=-1:
        model = 'ocean'
        modelV = model
        varnamef = 'SST'
        latnamef = 'geolat'   # should be 2D
        lonnamef = 'geolon'   # should be 2D
        ffile = 'gfs.ocean.t00z.6hr_avg.f'+str(fhour)+'.nc'
        vfile = 'gdas.t00z.ocnana.nc'
        varnameV = 'Temp'
      elif dpath.find('atm_da')!=-1:
        model = 'ocean'
        modelV = 'atmos'
        varnamef = 'SST'
        latnamef = 'geolat'   # should be 2D
        lonnamef = 'geolon'   # should be 2D
        ffile = 'gfs.ocean.t00z.6hr_avg.f'+str(fhour)+'.nc'
        vfile = 'gdas.t00z.sfcanl.nc'
        varnameV = 'tmpsfc'
        #varnamef = 'tmp2m' 
      else:
        model = 'atmos'
        modelV = model
        varnamef = 'tmpsfc'
        #varnamef = 'tmp2m'
        latnamef = 'lat' # should be 2D
        lonnamef = 'lon' # shoudl be 2D
        ffile = 'gfs.t00z.sfcf'+str(fhour)+'.nc'
        vfile = 'gdas.t00z.sfcanl.nc'
        varnameV = varnamef
    elif var=="atmosT":
            # atmos surface temperature
      model = 'atmos'
      modelV = model
      varnamef = 'tmpsfc'  
      #varnamef = 'tmp2m'  
      latnamef = 'lat' # should be 2D 
      lonnamef = 'lon' # shoudl be 2D
      ffile = 'gfs.t00z.sfcf'+str(fhour)+'.nc'
      vfile = 'gdas.t00z.sfcanl.nc'
      varnameV = varnamef

    varstr = varnamef

    # forecast path, etc.
    if dpath.find("role")!=-1:
      fpath = dpath+"/"+yyyy+mm+dd+"00/gfs."+str(yyyy)+str(mm)+str(dd)+"/00/model_data/"+str(model)+"/history/"
    else:
      fpath = dpath+"/gfs."+str(yyyy)+str(mm)+str(dd)+"/00/model_data/"+str(model)+"/history/"
    fname = fpath+"/"+ffile

    # verification date
    yyyyV, mmV, ddV = get_verif_date_f06(yyyy, mm, dd, fhour)
    print("read_fcst_verif: VERIF DATE = "+str(yyyyV)+str(mmV)+str(ddV))

    # verification path, etc.
    if dpath.find("role")!=-1:
      vpath = dpath+"/"+yyyy+mm+dd+"00/gdas."+str(yyyyV)+str(mmV)+str(ddV)+"/00/analysis/"+str(modelV)+"/"
    else:
      vpath = dpath+"/gdas."+str(yyyyV)+str(mmV)+str(ddV)+"/00/analysis/"+str(modelV)+"/"
    vname = vpath+"/"+vfile
    vstr  = "ANL"

    print("read_fcst_verif: forecast file = "+fname)
    if exists(str(fname))==True:
      print("read_fcst_verif: forecast file exists")
      ds   = xr.open_dataset(str(fname))
      latf = np.squeeze(ds[latnamef][:,:])
      lonf = np.squeeze(ds[lonnamef][:,:])
      if var=="atmosT":
        fcst = np.squeeze(ds[varnamef][0,:,:]) - 273.15     # convert to degrees-C
        fmask = np.squeeze(ds['land'][:,:])                 # land: mask=1      
        fcstmask = np.where(fmask>0, np.nan, fcst)      # mask=0 ocean; mask=1 land; mask=2 sea ice
      elif var=="sst":
        fcst = np.squeeze(ds[varnamef][0,:,:])
        fcstmask = np.where(fcst<-999, np.nan, fcst)
      ds.close()
    else:
      fcstmask = "NO_FILE"
      latf = fcstmask
      lonf = fcstmask

    print("read_fcst_verif: verification file = "+vname)
    if exists(str(vname))==True:
      print("read_fcst_verif: verification file exists")
      ds    = xr.open_dataset(str(vname))
      if var=="atmosT":
        verif = np.squeeze(ds[varnameV][0,:,:]) - 273.15
        vmask = np.squeeze(ds['land'][:,:])               
        verifmask = np.where(vmask==1, np.nan, verif)
      elif var=="sst" and dpath.find('s2s_da')!=-1:
        verif = np.squeeze(ds[varnameV][0,0,:,:])
        verifmask = np.where(verif==0, np.nan, verif)
      else:
        verif = np.squeeze(ds[varnameV][0,:,:]) - 273.15
        vmask = np.squeeze(ds['land'][:,:])
        verifmask = np.where(vmask==1, np.nan, verif)
      ds.close()
    else:
      verifmask = "NO_FILE"

    return fcstmask, verifmask, latf, lonf, varstr, vstr #, lat2dV, lon2dV





####################################################################################################
####################################################################################################
####################################################################################################
# OLD

#===============================================================================================
# Get ocean model background or analysis or increment, one file.
#
#	Return: Data for file, approx. sun position
#
def get_file(dpath, yyyy, mm, dd, hour, filetype):
   
    varname = 'Temp'

    nhours = np.size(hour)
    print("get_file: nhours = "+str(nhours))

	# date before yyyymmdd
    		# below is for 2021-07 only
    tddB = int(dd) - 1
    if tddB<10: 
      ddB = "0"+str(tddB)
    else:
      ddB = str(tddB)

	# set Sun location
    if np.size(hour)==1:
      if hour[0]=='00':
        sunlat = 20
        sunlon = -170
      elif hour[0]=='06':
        sunlat = 20
        sunlon = 100
      elif hour[0]=='12':
        sunlat = 25
        sunlon = 5
      elif hour[0]=='18':
        sunlat = 25
        sunlon = -85

	# get filenames
    if filetype=="f06":		# only 6-h background forecast
      fname = dpath+'gdas.'+yyyy+mm+dd+'/'+str(hour[0])+'/model_data/ocean/history/gdas.t'+str(hour[0])+'z.ocnf006.nc'
    elif filetype=="anl":       # analysis file
      fname = dpath+'gdas.'+yyyy+mm+dd+'/'+str(hour[0])+'/analysis/ocean/gdas.t'+str(hour[0])+'z.ocnana.nc'
    elif filetype=="inc":	# analysis increment
      fname = dpath+'gdas.'+yyyy+mm+dd+'/'+str(hour[0])+'/analysis/ocean/gdas.t'+str(hour[0])+'z.ocninc.nc'
      
    tfnames = fname

	# get number of files
    nfiles = np.size(tfnames)
    print("get_file: nfiles = "+str(nfiles))

    if nfiles==0:
      print("get_file: NO FILES for "+yyyy+mm+dd+" at "+str(hour))
      sst = []
      return sst, sunlon, sunlat

	# get shape of data
		# file for setting continents to NaN
    testfile = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.t00z.ocnf006.nc'
    ds = xr.open_dataset(testfile)
    t = np.squeeze(ds[varname][0,0,:,:])
    shape_t = np.shape(t)
    #print("min/max t = "+str(np.min(t))+" "+str(np.max(t)))
    t = np.where(t==0, np.nan, t)
    #print("min/max t = "+str(np.min(t))+" "+str(np.max(t)))
    #del t
    ds.close()

	# extract SST for specified hour block
    sst = np.empty([nfiles,shape_t[0],shape_t[1]],dtype=float)	# for looking at individual hour blocks			
    fname = tfnames
    ds = xr.open_dataset(fname)
    print("get_file: file = "+str(fname))

		# SST by the hour
    t_sst = np.squeeze(ds[varname][0,0,:,:])             # units = degrees-C
    sst[0,:,:] = np.array(t_sst)
    del t_sst

    ds.close()
    print("get_file: min/max sst = "+str(np.nanmin(sst))+" "+str(np.nanmax(sst)))

    return sst, sunlon, sunlat


#===============================================================================================
# Get ocean SST analysis OmB, OmA from diags
#
def get_omba_diags(exptpath, yyyy, mm, dd, hour, absopt, obstype):

    if obstype=="all":
	# all SST obs diags
      if dpath.find("role")!=-1:
        gdas_fnames = glob.glob(exptpath+yyyy+mm+dd+hour+'/gdas.'+yyyy+mm+dd+'/'+str(hour)+'/analysis/ocean/diags/sst*.'+yyyy+mm+dd+str(hour)+'.nc4')
      else:
        gdas_fnames = glob.glob(exptpath+'gdas.'+yyyy+mm+dd+'/'+str(hour)+'/analysis/ocean/diags/sst*.'+yyyy+mm+dd+str(hour)+'.nc4')
    else:
	# individual obs diags
      if dpath.find("role")!=-1:
        gdas_fnames = glob.glob(exptpath+yyyy+mm+dd+hour+'/gdas.'+yyyy+mm+dd+'/'+str(hour)+'/analysis/ocean/diags/sst_'+str(obstype)+'*.'+yyyy+mm+dd+str(hour)+'.nc4')
      else:
        gdas_fnames = glob.glob(exptpath+'gdas.'+yyyy+mm+dd+'/'+str(hour)+'/analysis/ocean/diags/sst_'+str(obstype)+'*.'+yyyy+mm+dd+str(hour)+'.nc4')

    print("get_omba: gdas_fnames = "+str(gdas_fnames)+" size = "+str(np.size(gdas_fnames)))

    if np.size(gdas_fnames)==0: 
      omb = "NO FILES"
      return omb,omb,omb,omb,omb,omb,omb

    cnt = 0.0
    for gdas_fname in gdas_fnames:
      ds = Dataset(gdas_fname)

      omb  = np.squeeze(ds.groups['ombg'].variables['seaSurfaceTemperature'][:])
      oma  = np.squeeze(ds.groups['oman'].variables['seaSurfaceTemperature'][:])
      obserr = np.squeeze(ds.groups['ObsError'].variables['seaSurfaceTemperature'][:]) 
      qc   = np.squeeze(ds.groups['EffectiveQC1'].variables['seaSurfaceTemperature'][:])	# 0 = obs are accepted/assimilated
      time = np.squeeze(ds.groups['MetaData'].variables['dateTime'][:])			# seconds since 1970-01-01T00:00:00Z
      lats = np.squeeze(ds.groups['MetaData'].variables['latitude'][:])
      lons = np.squeeze(ds.groups['MetaData'].variables['longitude'][:])

                # only include accepted (assimilated) obs
      idx = np.where(qc==0)
      qomb    = omb[idx]
      qoma    = oma[idx]
      qobserr = obserr[idx]
      qqc     = qc[idx]
      qtime   = time[idx]
      qlats   = lats[idx]
      qlons   = lons[idx]
      del idx

      if cnt==0:
        omb_arr    = qomb
        oma_arr    = qoma
        obserr_arr = qobserr
        qc_arr     = qqc
        time_arr   = qtime
        lats_arr   = qlats
        lons_arr   = qlons
      else:
        omb_arr    = np.append(omb_arr, qomb, axis=0)
        oma_arr    = np.append(oma_arr, qoma, axis=0)
        obserr_arr = np.append(obserr_arr, qobserr, axis=0)
        qc_arr     = np.append(qc_arr, qqc, axis=0)
        time_arr   = np.append(time_arr, qtime, axis=0)
        lats_arr   = np.append(lats_arr, qlats, axis=0)
        lons_arr   = np.append(lons_arr, qlons, axis=0)

      del omb,oma,obserr,qc,time,lats,lons
      del qomb,qoma,qobserr,qqc,qtime,qlats,qlons
      cnt += 1
      ds.close()
        
    return omb_arr, oma_arr, obserr_arr, qc_arr, time_arr, lats_arr, lons_arr
    
#===============================================================================================
# Get ocean SST analysis OmB, OmA from diags
#
#	For B-matrix testing only
#
def get_omba_diags_Bmatrix(exptpath, yyyy, mm, dd, hour, absopt, obstype, matrix):

    if obstype=="all":
	# all SST obs diags
      gdas_fnames = glob.glob(exptpath+'diags_'+str(matrix)+'/sst*.'+yyyy+mm+dd+str(hour)+'.nc4')
    else:
	# individual obs diags
      gdas_fnames = glob.glob(exptpath+'diags_'+str(matrix)+'/sst_'+str(obstype)+'*.'+yyyy+mm+dd+str(hour)+'.nc4')

    print("get_omba: gdas_fnames = "+str(gdas_fnames)+" size = "+str(np.size(gdas_fnames)))

    if np.size(gdas_fnames)==0: 
      omb = "NO FILES"
      return omb,omb,omb,omb,omb,omb,omb

    cnt = 0.0
    for gdas_fname in gdas_fnames:
      ds = Dataset(gdas_fname)

      omb  = np.squeeze(ds.groups['ombg'].variables['seaSurfaceTemperature'][:])
      oma  = np.squeeze(ds.groups['oman'].variables['seaSurfaceTemperature'][:])
      obserr = np.squeeze(ds.groups['ObsError'].variables['seaSurfaceTemperature'][:]) 
      qc   = np.squeeze(ds.groups['EffectiveQC1'].variables['seaSurfaceTemperature'][:])	# 0 = obs are accepted/assimilated
      time = np.squeeze(ds.groups['MetaData'].variables['dateTime'][:])			# seconds since 1970-01-01T00:00:00Z
      lats = np.squeeze(ds.groups['MetaData'].variables['latitude'][:])
      lons = np.squeeze(ds.groups['MetaData'].variables['longitude'][:])

                # only include accepted (assimilated) obs
      idx = np.where(qc==0)
      qomb    = omb[idx]
      qoma    = oma[idx]
      qobserr = obserr[idx]
      qqc     = qc[idx]
      qtime   = time[idx]
      qlats   = lats[idx]
      qlons   = lons[idx]
      del idx

      if cnt==0:
        if absopt==0:
          omb_arr    = qomb
          oma_arr    = qoma
        elif absopt==1:
          omb_arr    = abs(qomb)
          oma_arr    = abs(qoma)
        obserr_arr = qobserr
        qc_arr     = qqc
        time_arr   = qtime
        lats_arr   = qlats
        lons_arr   = qlons
      else:
        if absopt==0:
          omb_arr    = np.append(omb_arr, qomb, axis=0)
          oma_arr    = np.append(oma_arr, qoma, axis=0)
        elif absopt==1:
          omb_arr    = np.append(omb_arr, abs(qomb), axis=0)
          oma_arr    = np.append(oma_arr, abs(qoma), axis=0)
        obserr_arr = np.append(obserr_arr, qobserr, axis=0)
        qc_arr     = np.append(qc_arr, qqc, axis=0)
        time_arr   = np.append(time_arr, qtime, axis=0)
        lats_arr   = np.append(lats_arr, qlats, axis=0)
        lons_arr   = np.append(lons_arr, qlons, axis=0)

      del omb,oma,obserr,qc,time,lats,lons
      del qomb,qoma,qobserr,qqc,qtime,qlats,qlons
      cnt += 1
      ds.close()
        
    return omb_arr, oma_arr, obserr_arr, qc_arr, time_arr, lats_arr, lons_arr
    
#===============================================================================================
# Get ocean SST analysis OmB, OmA from diags - for time series
#
def get_omba_diags_TS(exptpath, yyyy, mm, dd, hour, absopt):

    if np.size(hour)==1:

	# all SST diags
      gdas_fnames = glob.glob(exptpath+'gdas.'+yyyy+mm+dd+'/'+str(hour)+'/analysis/ocean/diags/sst*.'+yyyy+mm+dd+str(hour)+'.nc4')
	# VIIRS only
      #gdas_fnames = glob.glob(exptpath+'gdas.'+yyyy+mm+dd+'/'+str(hour)+'/analysis/ocean/diags/sst_viirs*.'+yyyy+mm+dd+str(hour)+'.nc4')
	# METOP* only
      #gdas_fnames = glob.glob(exptpath+'gdas.'+yyyy+mm+dd+'/'+str(hour)+'/analysis/ocean/diags/sst_metop*.'+yyyy+mm+dd+str(hour)+'.nc4')

    else:

        # all SST diags
      gdas_fnames = glob.glob(exptpath+'gdas.'+yyyy+mm+dd+'/??/analysis/ocean/diags/sst*.'+yyyy+mm+dd+'??.nc4')
        # VIIRS only
      #gdas_fnames = glob.glob(exptpath+'gdas.'+yyyy+mm+dd+'/??/analysis/ocean/diags/sst_viirs*.'+yyyy+mm+dd+'??.nc4')
        # METOP* only
      #gdas_fnames = glob.glob(exptpath+'gdas.'+yyyy+mm+dd+'/??/analysis/ocean/diags/sst_metop*.'+yyyy+mm+dd+'??.nc4')

    print("get_omba: gdas_fnames = "+str(gdas_fnames))

    if np.size(gdas_fnames)==0: return "NO FILES"

    cnt = 0.0
    for gdas_fname in gdas_fnames:
        #ds = xr.open_dataset(gdas_fname)
        ds = Dataset(gdas_fname)

        omb  = np.squeeze(ds.groups['ombg'].variables['seaSurfaceTemperature'][:])
        oma  = np.squeeze(ds.groups['oman'].variables['seaSurfaceTemperature'][:])
        obserr = np.squeeze(ds.groups['ObsError'].variables['seaSurfaceTemperature'][:]) 
        qc   = np.squeeze(ds.groups['EffectiveQC1'].variables['seaSurfaceTemperature'][:])        # 0 = obs are accepted/assimilated

		# only include accepted (assimilated) obs
        idx = np.where(qc==0)
        qomb    = omb[idx]
        qoma    = oma[idx]
        qobserr = obserr[idx]
        del idx

        if cnt == 0:
          if absopt==0:
            meanomb = np.mean(qomb)
            meanoma = np.mean(qoma)
            meanobserr = np.mean(qobserr)
          elif absopt==1: 
            meanomb = np.mean(abs(qomb))
            meanoma = np.mean(abs(qoma))
            meanobserr = np.mean(abs(qobserr))
          SDomb = np.std(qomb)
          SDoma = np.std(qoma)
          nobs = np.size(qomb)
        else:
          if absopt==0:
            meanomb += np.mean(qomb)
            meanoma += np.mean(qoma)
            meanobserr += np.mean(qobserr)
          elif absopt==1:
            meanomb += np.mean(abs(qomb))
            meanoma += np.mean(abs(qoma))
            meanobserr += np.mean(abs(qobserr))
          SDomb += np.std(qomb)
          SDoma += np.std(qoma)
          nobs += np.size(qomb)

        cnt +=1
        del omb,oma,obserr,qc,qomb,qoma,qobserr
        ds.close()
        
    return meanomb/cnt, meanoma/cnt, SDomb/cnt, SDoma/cnt, meanobserr/cnt, nobs
    
#===============================================================================================
# Get ocean model background, Compute daily mean and anomalies (deviations from daily mean)
#	(all background files for each day)
#
#	Return: Daily mean, Data for each file opened, approx. sun position
#
def get_anl_inc_Bmatrix(dpath, yyyy, mm, dd, hours, filetype, var):

    if var=="SST": varname = 'Temp'
    if var=="SSH": varname = 'ave_ssh'

    nhours = np.size(hours)
    print("get_bkg_mean_anom: hour = "+str(hours)+" nhours = "+str(nhours))

    hoursDAY = ['00','06','12','18']		# used for computing daily mean

	# date before yyyymmdd
    		# below is for 2021-07 only
    tddB = int(dd) - 1
    if tddB<10: 
      ddB = "0"+str(tddB)
    else:
      ddB = str(tddB)

	# set Sun location
    if np.size(hours)==1:
      if hours[0]=='00':
        sunlat = 20
        sunlon = -170
      elif hours[0]=='06':
        sunlat = 20
        sunlon = 100
      elif hours[0]=='12':
        sunlat = 25
        sunlon = 5
      elif hours[0]=='18':
        sunlat = 25
        sunlon = -85

	# get filenames
    if filetype=="anl":               # analysis files
      fnames00_ana = dpath+'ocn.3dvarfgat_pseudo.an.'+yyyy+'-'+mm+'-'+dd+'T00:00:00Z.nc'
      fnames06_ana = dpath+'ocn.3dvarfgat_pseudo.an.'+yyyy+'-'+mm+'-'+dd+'T06:00:00Z.nc'
      fnames12_ana = dpath+'ocn.3dvarfgat_pseudo.an.'+yyyy+'-'+mm+'-'+dd+'T12:00:00Z.nc'
      fnames18_ana = dpath+'ocn.3dvarfgat_pseudo.an.'+yyyy+'-'+mm+'-'+dd+'T18:00:00Z.nc'

      if np.size(hours)==1:
        if hours[0]=='00':
          tfnames = [fnames00_ana]
        elif hours[0]=='06':
          tfnames = [fnames06_ana]
        elif hours[0]=='12':
          tfnames = [fnames12_ana]
        elif hours[0]=='18':
          tfnames = [fnames18_ana]

    elif filetype=="inc":		# analysis increments
      fnames00_inc = dpath+'ocn.3dvarfgat_pseudo.incr.'+yyyy+'-'+mm+'-'+dd+'T00:00:00Z.nc'
      fnames06_inc = dpath+'ocn.3dvarfgat_pseudo.incr.'+yyyy+'-'+mm+'-'+dd+'T06:00:00Z.nc'
      fnames12_inc = dpath+'ocn.3dvarfgat_pseudo.incr.'+yyyy+'-'+mm+'-'+dd+'T12:00:00Z.nc'
      fnames18_inc = dpath+'ocn.3dvarfgat_pseudo.incr.'+yyyy+'-'+mm+'-'+dd+'T18:00:00Z.nc'

      if np.size(hours)==1:
        if hours[0]=='00':
          tfnames = [fnames00_inc]
        elif hours[0]=='06':
          tfnames = [fnames06_inc]
        elif hours[0]=='12':
          tfnames = [fnames12_inc]
        elif hours[0]=='18':
          tfnames = [fnames18_inc]


	# get number of files
    nfiles = np.size(tfnames)
    print("get_bkg_mean_anom: nfiles = "+str(nfiles))

    if nfiles==0:
      print("get_bkg_mean_anom: NO FILES for "+yyyy+mm+dd+" at "+str(hours))
      sst_sum = []
      sst     = []
      return sst_sum, sst

	# get shape of data
		# file for setting continents to NaN
    testfile = dpath+'ocn.3dvarfgat_pseudo.incr.'+yyyy+'-'+mm+'-'+dd+'T06:00:00Z.nc'
    ds = xr.open_dataset(testfile)
		# SST
    if var=="SST":
      tSST = np.squeeze(ds[varname][0,:,:,:])
      shape_tSST = np.shape(tSST)
      tSST = np.where(tSST==0, np.nan, tSST)
		# SSH
    if var=="SSH":
      tSST = np.squeeze(ds[varname][0,:,:])
      shape_tSST = np.shape(tSST)
      tSST = np.where(tSST==0, np.nan, tSST)
    ds.close()

	# extract SST for specified hour block
    if var=="SST": 
      sst = np.empty([nfiles,shape_tSST[0],shape_tSST[1],shape_tSST[2]],dtype=float)
      x   = np.empty([nfiles,shape_tSST[2]],dtype=float)		# x index
      y   = np.empty([nfiles,shape_tSST[1]],dtype=float)		# y index
      z   = np.empty([nfiles,shape_tSST[0]],dtype=float)		# z index
      h   = np.empty([nfiles,shape_tSST[0],shape_tSST[1],shape_tSST[2]],dtype=float)		# thickness of each ocean layer
    if var=="SSH": ssh = np.empty([nfiles,shape_tSSH[0],shape_tSSH[1]],dtype=float)
    cnt = 0.0	# sst hour file count
    for ih in range(np.size(hours)):			
      fnames = tfnames
      for fname in fnames:
        ds = xr.open_dataset(fname)
        print("get_bkg_mean_anom: file = "+str(fname))

		# SST by the hour
        if var=="SST":
          t_sst = np.squeeze(ds[varname][0,:,:,:])             # units = degrees C
          tx    = np.squeeze(ds["xaxis_1"][:])
          ty    = np.squeeze(ds["yaxis_1"][:])
          tz    = np.squeeze(ds["zaxis_1"][:])
          print("-----dpath = "+str(dpath))
          if dpath.find("updatedGDAS")==-1: 
            print("no updatedGDAS") 
            th = np.squeeze(ds["h"][0,:,:,:])
          if cnt == 0:
            sst[0,:,:,:] = np.array(t_sst)
            x[0,:]     = np.array(tx)
            y[0,:]     = np.array(ty)
            z[0,:]     = np.array(tz)
            if dpath.find("updatedGDAS")==-1: 
              print("no updatedGDAS") 
              h[0,:,:,:] = np.array(th) 
              del th
            del t_sst,tx,ty,tz
          else:
            sst[int(cnt),:,:,:] = np.array(t_sst)
            x[int(cnt),:]       = np.array(tx)
            y[int(cnt),:]       = np.array(ty)
            z[int(cnt),:]       = np.array(tz)
            if dpath.find("updatedGDAS")==-1: 
              print("no updatedGDAS") 
              h[int(cnt),:,:,:] = np.array(th) 
              del th
            del t_sst,tx,ty,tz

        if var=="SSH":
          t_sst = np.squeeze(ds[varname][0,:,:])             # units = m ??
          if cnt == 0:
            sst[0,:,:] = np.array(t_sst)
            del t_sst
          else:
            sst[int(cnt),:,:] = np.array(t_sst)
            del t_sst

        cnt +=1
        ds.close()
    print("get_Bmatrix: shape h = "+str(np.shape(h)))
    print("get_Bmatrix: min/max sst = "+str(np.nanmin(sst))+" "+str(np.nanmax(sst)))
    print("get_Bmatrix: min/max h = "+str(np.nanmin(h))+" "+str(np.nanmax(h)))

    return sst, h, z, y, x, sunlon, sunlat
    
#===============================================================================================
# Get grid (lat,lon) values for Bmatrix test plots 
#
#	Return: latitudes, longitudes
#
def get_grid_Bmatrix(dpath):

    fname = glob.glob(dpath+'soca_gridspec*')

    for iname in fname:		# should only be one file
      ds = xr.open_dataset(iname)

      tlat = np.squeeze(ds['lat'][0,:,:])
      tlon = np.squeeze(ds['lon'][0,:,:])
    
      ds.close()

    lat = np.asarray(tlat)
    lon = np.asarray(tlon)

    return lat, lon
    
#===============================================================================================
# Get ocean model background, Compute NIGHTTIME daily mean, and anomalies (deviations from daily mean)
#	(all background files for each day)
#
#	Need to make sure night "moves" with each forecast time
#	Need to compute the day/night terminator position
#
#	Return: Daily mean, Data for each file opened, approx. sun position
#
def get_timemean_anom_nightonly(yyyy, mm, dd, nlon):

    #dpath = '/scratch2/NCEPDEV/ocean/Guillaume.Vernieres/runs/cp0/cp0.b/COMROT/cp0.b/'
    dpath = '/scratch1/NCEPDEV/da/Katherine.Lukens/NSST/expts/wcda_expt.NO_marineDA/cp0.no-ocn-da/COMROT/cp0.no-ocn-da/'
    #dpath = '/scratch1/NCEPDEV/da/Katherine.Lukens/NSST/expts/wcda_expt.marineDA_sst/cp0.ocn-da.sst/COMROT/cp0.ocn-da.sst/'

    varname = 'Temp'

    #hours = ['00']
    #hours = ['06']
    #hours = ['12']
    hours = ['18']
    nhours = np.size(hours)

    hoursDAY = ['00','06','12','18']		# used for computing daily mean

	# date before yyyymmdd
		# below is for 2021-07 only
    tddB = int(dd) - 1
    if tddB<10:
      ddB = "0"+str(tddB)
    else:
      ddB = str(tddB)

	# set Sun location
    if np.size(hours)==1:
      if hours[0]=='00':
        sunlat = 20
        sunlon = -170
      elif hours[0]=='06':
        sunlat = 20
        sunlon = 100
      elif hours[0]=='12':
        sunlat = 25
        sunlon = 5
      elif hours[0]=='18':
        sunlat = 25
        sunlon = -85

	# get filenames
		# before date
    Bfnames18_f006 = dpath+'gdas.'+yyyy+mm+ddB+'/18/model_data/ocean/history/gdas.t18z.ocnf006.nc'
    Bfnames18_f009 = dpath+'gdas.'+yyyy+mm+ddB+'/18/model_data/ocean/history/gdas.t18z.ocnf009.nc'
		# current date
    fnames00_f003 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.t00z.ocnf003.nc'
    fnames00_f006 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.t00z.ocnf006.nc'
    fnames00_f009 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.t00z.ocnf009.nc'
    fnames06_f003 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.t06z.ocnf003.nc'
    fnames06_f006 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.t06z.ocnf006.nc'
    fnames06_f009 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.t06z.ocnf009.nc'
    fnames12_f003 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.t12z.ocnf003.nc'
    fnames12_f006 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.t12z.ocnf006.nc'
    fnames12_f009 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.t12z.ocnf009.nc'
    fnames18_f003 = dpath+'gdas.'+yyyy+mm+dd+'/18/model_data/ocean/history/gdas.t18z.ocnf003.nc'

    if np.size(hours)==1:
      if hours[0]=='00':
        tfnames = [Bfnames18_f006]
        tfnames.append(Bfnames18_f009)
        tfnames.append(fnames00_f003)
      elif hours[0]=='06':
        tfnames = [fnames00_f006]
        tfnames.append(fnames00_f009)
        tfnames.append(fnames06_f003)
      elif hours[0]=='12':
        tfnames = [fnames06_f006]
        tfnames.append(fnames06_f009)
        tfnames.append(fnames12_f003)
      elif hours[0]=='18':
        tfnames = [fnames12_f006]
        tfnames.append(fnames12_f009)
        tfnames.append(fnames18_f003)

    fnamesDAILY = [Bfnames18_f006]
    fnamesDAILY.append(Bfnames18_f009)
    fnamesDAILY.append(fnames00_f003)
    fnamesDAILY.append(fnames00_f006)
    fnamesDAILY.append(fnames00_f009)
    fnamesDAILY.append(fnames06_f003)
    fnamesDAILY.append(fnames06_f006)
    fnamesDAILY.append(fnames06_f009)
    fnamesDAILY.append(fnames12_f003)
    fnamesDAILY.append(fnames12_f006)
    fnamesDAILY.append(fnames12_f009)
    fnamesDAILY.append(fnames18_f003)

	# get number of files
    nfiles = np.size(tfnames)
    nfilesDAILY = np.size(fnamesDAILY)
    print("get_bkg_night: nfiles = "+str(nfiles))

    if nfiles==0:
      print("get_bkg_night: NO FILES for "+yyyy+mm+dd+" at "+str(hours))
      sst_sum = []
      sst     = []
      return sst_sum, sst

	# get shape of data
    ds = xr.open_dataset(Bfnames18_f006)
    t = np.squeeze(ds[varname][0,0,:,:])
    shape_t = np.shape(t)
    del t
    ds.close()

    lats, lons = get_grid()
    shapelat = np.shape(lats)

	# compute daily mean
    print("get_bkg_night: daily mean LOOP")
    sst_sum  = []								# for computing daily mean
    cntnight = 0.0
    cntDAILY   = 0.0			# daily mean file count
    for i in range(np.size(fnamesDAILY)):
      ds = xr.open_dataset(str(fnamesDAILY[i]))
      print("get_bkg_night: file = "+str(fnamesDAILY[i]))

        	# get day/night terminator and create mask
			# get hour for terminator computation
      if i==0: 
        hh = 0			# first file in day: f06 forecast from 18z on previous day
      elif i==1:
        hh = 3			# second file in day: f09 forecast from 18z on previous day
      elif str(fnamesDAILY[i]).find('t00z.ocnf003')!=-1:
        hh = 3 
      elif str(fnamesDAILY[i]).find('t00z.ocnf006')!=-1: 
        hh = 6
      elif str(fnamesDAILY[i]).find('t00z.ocnf009')!=-1: 
        hh = 9
      elif str(fnamesDAILY[i]).find('t06z.ocnf003')!=-1: 
        hh = 9
      elif str(fnamesDAILY[i]).find('t06z.ocnf006')!=-1:
        hh = 12
      elif str(fnamesDAILY[i]).find('t06z.ocnf009')!=-1: 
        hh = 15
      elif str(fnamesDAILY[i]).find('t12z.ocnf003')!=-1: 
        hh = 15
      elif str(fnamesDAILY[i]).find('t12z.ocnf006')!=-1:
        hh = 18
      elif str(fnamesDAILY[i]).find('t12z.ocnf009')!=-1: 
        hh = 21
      elif str(fnamesDAILY[i]).find('t18z.ocnf003')!=-1: 
        hh = 21
      print("get_bkg_night: hh = "+str(hh))

      print("get_bkg_night: date = "+str(yyyy)+str(mm)+str(dd)+str(hh))
      lonmin = -180
      lonmax = 180
      term_lats, term_lons = terminator(yyyy,mm,dd,hh,00, lonmin, lonmax, nlon)

                # daily mean
			# nighttime gridpoints (approx.) during July 2021
      flat = np.linspace(-90,90,shapelat[0])
      flon = np.linspace(-180,180,shapelat[1])
      meshlon, meshlat = np.meshgrid(flon,flat)
      inight = np.ones(np.shape(meshlon),int)
      print("shape inight = "+str(np.shape(inight)))
      for k in range(shapelat[1]):
        if int(mm) == 7:
          #print("term_lats = "+str(term_lats[k])+" term_lons = "+str(term_lons[k]))
          inight[:,k] = np.where(meshlat[:,k]>term_lats[k], 0, 1)	# night gridpoints = 1

      if cntDAILY == 0:
        t_sst = np.squeeze(ds[varname][0,0,:,:])             # units = degrees-C
        sst_sum = t_sst*inight #np.where(inight==1, t_sst, 0)
        #sst_sum = np.zeros(np.shape(t_sst),float)
        #for k in range(shapelat[1]):
        #  if int(mm) == 7:
        #    sst_sum[:,k] = np.where(meshlat[:,k]>term_lats[k], 0, t_sst[:,k])
        cntnight = inight
        del t_sst
      else:
        t_sst = np.squeeze(ds[varname][0,0,:,:])
        tsst_sum = t_sst*inight #np.where(inight==1, t_sst, 0)
        #tss = np.zeros(np.shape(t_sst),float)
        #for k in range(shapelat[1]):
        #  if int(mm) == 7:
        #    tss[:,k] = np.where(meshlat[:,k]>term_lats[k], 0, t_sst[:,k])
        #sst_sum += tss
        #del tss
        sst_sum += tsst_sum
        cntnight += inight
        del t_sst,tsst_sum

      cntDAILY +=1
      del hh,term_lats,term_lons,flat,flon,meshlon,meshlat
      del inight
      ds.close()

    print("cntnight = "+str(cntnight))

    if np.size(sst_sum)==0:
      daily_mean = []
    else:
      daily_mean = np.array(sst_sum)/cntnight #cntDAILY
    print("get_bkg_night: min/max daily_mean = "+str(np.nanmin(daily_mean))+" "+str(np.nanmax(daily_mean)))

	# extract SST for specified hour block
    print("get_bkg_night: SST hour LOOP")
    sst     = np.empty([nfiles,shape_t[0],shape_t[1]],dtype=float)	# for looking at individual hour blocks
    cnt     = 0.0	# sst hour file count
    for ih in range(np.size(hours)):			
      fnames = tfnames
      for fname in fnames:
        ds = xr.open_dataset(fname)
        print("get_bkg_night: file = "+str(fname))

		# SST by the hour
        if cnt == 0:
            t_sst = np.squeeze(ds[varname][0,0,:,:])             # units = degrees-C
            sst[0,:,:] = np.array(t_sst)
            del t_sst
        else:
            t_sst = np.squeeze(ds[varname][0,0,:,:])
            sst[int(cnt),:,:] = np.array(t_sst)
            del t_sst

        cnt +=1
        ds.close()
    print("get_bkg_night: count = "+str(cnt))
    print("get_bkg_night: min/max sst = "+str(np.nanmin(sst))+" "+str(np.nanmax(sst)))

    return daily_mean, sst, sunlon, sunlat
    
#===============================================================================================
# Compare ocean model background forecasts at "same time": f09 from previous cycle and f03 from current cycle
#
#	Return: Data for each file opened, approx. sun position
#
def compare_bkg_sametime(yyyy, mm, dd):

    #dpath = '/scratch2/NCEPDEV/ocean/Guillaume.Vernieres/runs/cp0/cp0.b/COMROT/cp0.b/'
    #dpath = '/scratch1/NCEPDEV/da/Katherine.Lukens/NSST/expts/wcda_expt.NO_marineDA/cp0.no-ocn-da/COMROT/cp0.no-ocn-da/'
    dpath = '/scratch1/NCEPDEV/da/Katherine.Lukens/NSST/expts/wcda_expt.marineDA_sst/cp0.ocn-da.sst/COMROT/cp0.ocn-da.sst/'
    
    varname = 'Temp'

    #hours = ['00']
    #hours = ['06']
    #hours = ['12']
    hours = ['18']
    nhours = np.size(hours)
    print("compare_bkg_sametime: nhours = "+str(nhours))

    hoursDAY = ['00','06','12','18']		# used for computing daily mean

	# date before yyyymmdd
    		# below is for 2021-07 only
    tddB = int(dd) - 1
    if tddB<10: 
      ddB = "0"+str(tddB)
    else:
      ddB = str(tddB)

	# set Sun location
    if np.size(hours)==1:
      if hours[0]=='00':
        sunlat = 20
        sunlon = -170
      elif hours[0]=='06':
        sunlat = 20
        sunlon = 100
      elif hours[0]=='12':
        sunlat = 25
        sunlon = 5
      elif hours[0]=='18':
        sunlat = 25
        sunlon = -85

	# get filenames
                # before date
    Bfnames18_f009 = dpath+'gdas.'+yyyy+mm+ddB+'/18/model_data/ocean/history/gdas.t18z.ocnf009.nc'
		# current date
    fnames00_f003 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.t00z.ocnf003.nc'
    fnames00_f009 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.t00z.ocnf009.nc'
    fnames06_f003 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.t06z.ocnf003.nc'
    fnames06_f009 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.t06z.ocnf009.nc'
    fnames12_f003 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.t12z.ocnf003.nc'
    fnames12_f009 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.t12z.ocnf009.nc'
    fnames18_f003 = dpath+'gdas.'+yyyy+mm+dd+'/18/model_data/ocean/history/gdas.t18z.ocnf003.nc'

    if np.size(hours)==1:
      if hours[0]=='00':
        tfnames = [Bfnames18_f009]
        tfnames.append(fnames00_f003)
      elif hours[0]=='06':
        tfnames = [fnames00_f009]
        tfnames.append(fnames06_f003)
      elif hours[0]=='12':
        tfnames = [fnames06_f009]
        tfnames.append(fnames12_f003)
      elif hours[0]=='18':
        tfnames = [fnames12_f009]
        tfnames.append(fnames18_f003)

	# get number of files
    nfiles = np.size(tfnames)
    print("compare_bkg_sametime: nfiles = "+str(nfiles))

    if nfiles==0:
      print("compare_bkg_sametime: NO FILES for "+yyyy+mm+dd+" at "+str(hours))
      sst     = []
      return sst_sum, sst

	# get shape of data
    ds = xr.open_dataset(Bfnames18_f009)
    t = np.squeeze(ds[varname][0,0,:,:])
    shape_t = np.shape(t)
    del t
    ds.close()

	# extract SST for specified hour block
    print("compare_bkg_sametime: SST hour LOOP")
    sstdiff = np.empty([1,shape_t[0],shape_t[1]],dtype=float)
			
    fnames = tfnames
    ds0 = xr.open_dataset(fnames[0])			  # f09 (previous cycle)
    ds1 = xr.open_dataset(fnames[1])			  # f03 (current cycle)
    print("compare_bkg_sametime: files = "+str(fnames))

    tsst0 = np.squeeze(ds0[varname][0,0,:,:])		  # units = degrees-C
    tsst1 = np.squeeze(ds1[varname][0,0,:,:])
    sst0 = np.array(tsst0)
    sst1 = np.array(tsst1)
    
    sstdiff[0,:,:] = sst1 - sst0			  # diff: f03(current) - f09(previous)
    
    del tsst0,tsst1,sst0,sst1

    ds0.close()
    ds1.close()
    print("compare_bkg_sametime: min/max sst = "+str(np.nanmin(sstdiff))+" "+str(np.nanmax(sstdiff)))

    return sstdiff, sunlon, sunlat
    
#===============================================================================================
# Decompose analysis difference between 2 adjacent times into model and DA components: 
#	Components: analysis diff, model background - analysis, increments
#
#	Return: Daily mean, Data for each file opened, approx. sun position
#
def decomp_anl_components(yyyy, mm, dd):

    #dpath = '/scratch2/NCEPDEV/ocean/Guillaume.Vernieres/runs/cp0/cp0.b/COMROT/cp0.b/'
    ###dpath = '/scratch1/NCEPDEV/da/Katherine.Lukens/NSST/expts/wcda_expt.NO_marineDA/cp0.no-ocn-da/COMROT/cp0.no-ocn-da/'
    dpath = '/scratch1/NCEPDEV/da/Katherine.Lukens/NSST/expts/wcda_expt.marineDA_sst/cp0.ocn-da.sst/COMROT/cp0.ocn-da.sst/saved/'
    
    varname = 'Temp'

    hours = ['00']
    #hours = ['06']
    #hours = ['12']
    #hours = ['18']
    nhours = np.size(hours)
    print("get_bkg_anl_inc: nhours = "+str(nhours))

    hoursDAY = ['00','06','12','18']		# used for computing daily mean

	# date before yyyymmdd
    		# below is for 2021-07 only
    tddB = int(dd) - 1
    if tddB<10: 
      ddB = "0"+str(tddB)
    else:
      ddB = str(tddB)

	# set Sun location
    if np.size(hours)==1:
      if hours[0]=='00':
        sunlat = 20
        sunlon = -170
      elif hours[0]=='06':
        sunlat = 20
        sunlon = 100
      elif hours[0]=='12':
        sunlat = 25
        sunlon = 5
      elif hours[0]=='18':
        sunlat = 25
        sunlon = -85

	# get filenames: anl, bkg
                # before date
    Bfnames18_anl = dpath+'gdas.'+yyyy+mm+ddB+'/18/analysis/ocean/gdas.t18z.ocnana.nc'
		# current date
    fnames00_f006 = dpath+'gdas.'+yyyy+mm+dd+'/00/model_data/ocean/history/gdas.t00z.ocnf006.nc'
    fnames00_anl  = dpath+'gdas.'+yyyy+mm+dd+'/00/analysis/ocean/gdas.t00z.ocnana.nc'
    fnames06_f006 = dpath+'gdas.'+yyyy+mm+dd+'/06/model_data/ocean/history/gdas.t06z.ocnf006.nc'
    fnames06_anl  = dpath+'gdas.'+yyyy+mm+dd+'/06/analysis/ocean/gdas.t06z.ocnana.nc'
    fnames12_f006 = dpath+'gdas.'+yyyy+mm+dd+'/12/model_data/ocean/history/gdas.t12z.ocnf006.nc'
    fnames12_anl  = dpath+'gdas.'+yyyy+mm+dd+'/12/analysis/ocean/gdas.t12z.ocnana.nc'
    fnames18_f006 = dpath+'gdas.'+yyyy+mm+dd+'/18/model_data/ocean/history/gdas.t18z.ocnf006.nc'
    fnames18_anl  = dpath+'gdas.'+yyyy+mm+dd+'/18/analysis/ocean/gdas.t18z.ocnana.nc'
    
    	# get filenames: increments
    fnames00_inc = dpath+'gdas.'+yyyy+mm+dd+'/00/analysis/ocean/gdas.t00z.ocninc.nc'
    fnames06_inc = dpath+'gdas.'+yyyy+mm+dd+'/06/analysis/ocean/gdas.t06z.ocninc.nc'
    fnames12_inc = dpath+'gdas.'+yyyy+mm+dd+'/12/analysis/ocean/gdas.t12z.ocninc.nc'
    fnames18_inc = dpath+'gdas.'+yyyy+mm+dd+'/18/analysis/ocean/gdas.t18z.ocninc.nc'

    if np.size(hours)==1:
      if hours[0]=='00':
        tfnames = [Bfnames18_anl]
        tfnames.append(fnames00_f006)
        tfnamesAA = [Bfnames18_anl]
        tfnamesAA.append(fnames00_anl)
        tfnamesinc = [fnames00_inc]
      elif hours[0]=='06':
        tfnames = [fnames00_anl]
        tfnames.append(fnames06_f006)
        tfnamesAA = [fnames00_anl]
        tfnamesAA.append(fnames06_anl)
        tfnamesinc = [fnames06_inc]
      elif hours[0]=='12':
        tfnames = [fnames06_anl]
        tfnames.append(fnames12_f006)
        tfnamesAA = [fnames06_anl]
        tfnamesAA.append(fnames12_anl)
        tfnamesinc = [fnames12_inc]
      elif hours[0]=='18':
        tfnames = [fnames12_anl]
        tfnames.append(fnames18_f006)
        tfnamesAA = [fnames12_anl]
        tfnamesAA.append(fnames18_anl)
        tfnamesinc = [fnames18_inc]

	# get number of files
    nfiles = np.size(tfnames)
    print("get_bkg_anl_inc: nfiles = "+str(nfiles))

    if nfiles==0:
      print("get_bkg_anl_inc: NO FILES for "+yyyy+mm+dd+" at "+str(hours))
      A2mA1 = []
      B2mA1 = []
      inc2  = []
      return A2mA1, B2mA1, inc2, sunlon, sunlat

	# get shape of data
    ds = xr.open_dataset(Bfnames18_anl)
    t = np.squeeze(ds[varname][0,0,:,:])
    shape_t = np.shape(t)
    del t
    ds.close()

	# extract SST for specified hour block
    print("get_bkg_anl_inc: SST hour LOOP")
    ds0   = xr.open_dataset(tfnames[0])
    print("get_bkg_anl_inc: files = "+str(tfnames))
    sst0   = np.squeeze(ds0[varname][0,0,:,:])           # units = degrees-C
    ds0.close()
    ds1   = xr.open_dataset(tfnames[1])
    print("get_bkg_anl_inc: files = "+str(tfnames))
    sst1   = np.squeeze(ds1[varname][0,0,:,:])
    ds1.close()
    print("get_bkg_anl_inc: shape sst0 = "+str(np.shape(sst0))+" sst1 = "+str(np.shape(sst1)))

              # SST by the hour
    tA2mA1 = sst1 - sst0
    del sst0,sst1
    A2mA1 = np.empty([1,shape_t[0],shape_t[1]],dtype=float)
    A2mA1[0,:,:] = np.array(tA2mA1)
    del tA2mA1

    dsAA0   = xr.open_dataset(tfnamesAA[0])
    dsAA1   = xr.open_dataset(tfnamesAA[1])
    print("get_bkg_anl_inc: files = "+str(tfnamesAA))
    sstAA0   = np.squeeze(dsAA0[varname][0,0,:,:])           # units = degrees-C
    sstAA1   = np.squeeze(dsAA1[varname][0,0,:,:])
    tB2mA1 = sstAA1 - sstAA0
    del sstAA0,sstAA1
              # SST by the hour
    B2mA1 = np.empty([1,shape_t[0],shape_t[1]],dtype=float)
    B2mA1[0,:,:] = np.array(tB2mA1)
    del tB2mA1
    dsAA0.close()
    dsAA1.close()

    inc2  = np.empty([1,shape_t[0],shape_t[1]],dtype=float)
    dsinc = xr.open_dataset(tfnamesinc[0])
    print("get_bkg_anl_inc: files = "+str(tfnamesinc))
    sstinc = np.squeeze(dsinc[varname][0,0,:,:])      
              # SST by the hour
    inc2[0,:,:]  = np.array(sstinc)
    del sstinc
    dsinc.close()

    return A2mA1, B2mA1, inc2, sunlon, sunlat

#===============================================================================================
# Get ocean model background, Compute daily mean and anomalies (deviations from daily mean)
#	(all background files for each day)
#
#	Return: Daily mean, Data for each file opened, approx. sun position
#
def get_skinT_timemean_daily_hour(dpath, yyyy, mm, dd, hours):

    maskval = 1		# value indicating what feature to mask: 0=ocean, 1=land, 2=sea ice

    varname = 'tmpsfc'

    nhours = np.size(hours)

    hoursDAY = ['00','06','12','18']		# used for computing daily mean

	# date before yyyymmdd
    		# below is for 2021-07 only
    tddB = int(dd) - 1
    if tddB<10: 
      ddB = "0"+str(tddB)
    else:
      ddB = str(tddB)

	# set Sun location
    if np.size(hours)==1:
      if hours[0]=='00':
        sunlat = 20
        sunlon = -170
      elif hours[0]=='06':
        sunlat = 20
        sunlon = 100
      elif hours[0]=='12':
        sunlat = 25
        sunlon = 5
      elif hours[0]=='18':
        sunlat = 25
        sunlon = -85

	# get filenames
    fnames00 = dpath+'gdas.'+yyyy+mm+dd+'/00/analysis/atmos/gdas.t00z.sfcanl.nc'
    fnames06 = dpath+'gdas.'+yyyy+mm+dd+'/06/analysis/atmos/gdas.t06z.sfcanl.nc'
    fnames12 = dpath+'gdas.'+yyyy+mm+dd+'/12/analysis/atmos/gdas.t12z.sfcanl.nc'
    fnames18 = dpath+'gdas.'+yyyy+mm+dd+'/18/analysis/atmos/gdas.t18z.sfcanl.nc'

	# for specified hour block
    if np.size(hours)==1:
      if hours[0]=='00':
        tfnames = [fnames00]
      elif hours[0]=='06':
        tfnames = [fnames06]
      elif hours[0]=='12':
        tfnames = [fnames12]
      elif hours[0]=='18':
        tfnames = [fnames18]

	# for daily mean
    fnamesDAILY = [fnames00]
    fnamesDAILY.append(fnames06)
    fnamesDAILY.append(fnames12)
    fnamesDAILY.append(fnames18)

	# get number of files
    nfiles = np.size(tfnames)
    nfilesDAILY = np.size(fnamesDAILY)
    print("get_skinT_mean_anom: nfiles = "+str(nfiles))

    if nfiles==0:
      print("get_skinT_mean_anom: NO FILES for "+yyyy+mm+dd+" at "+str(hours))
      sst_sum = []
      sst     = []
      return sst_sum, sst

	# get shape of data
    ds = xr.open_dataset(fnames00)
    t = np.squeeze(ds[varname][0,:,:])
    lon = np.squeeze(ds['lon'][:,:])
    lat = np.squeeze(ds['lat'][:,:])
    shape_t = np.shape(t)
    del t
    ds.close()

	# compute daily mean
    print("get_skinT_mean_anom: daily mean LOOP")
    sst_sum = []								# for computing daily mean
    cntDAILY  = 0.0	# daily mean file count
    for i in range(np.size(fnamesDAILY)):
      print("get_skinT_mean_anom: fname = "+str(fnamesDAILY[i]))
      ds = xr.open_dataset(str(fnamesDAILY[i]))

                # daily mean
      if cntDAILY == 0:
        mask = np.squeeze(ds['land'][0,:,:])
        t_sst = np.squeeze(ds[varname][0,:,:]) - 273.15           # units = degrees-C
        sstmask = np.where(mask==maskval, np.nan, t_sst)
        sst_sum = sstmask
        del t_sst,mask,sstmask
      else:
        mask = np.squeeze(ds['land'][0,:,:])
        t_sst = np.squeeze(ds[varname][0,:,:]) - 273.15
        sstmask = np.where(mask==maskval, np.nan, t_sst)
        sst_sum += sstmask
        del t_sst,mask,sstmask

      cntDAILY +=1
      ds.close()

    if np.size(sst_sum)==0:
      daily_mean = []
    else:
      daily_mean = np.array(sst_sum)/cntDAILY
    print("get_skinT_mean_anom: min/max daily_mean = "+str(np.nanmin(daily_mean))+" "+str(np.nanmax(daily_mean)))

	# extract SST for specified hour block
    print("get_skinT_mean_anom: SST hour LOOP")
    sst = np.empty([nfiles,shape_t[0],shape_t[1]],dtype=float)	# for looking at individual hour blocks
    cnt = 0.0	# sst hour file count
    for ih in range(np.size(hours)):			
      fnames = tfnames
      for fname in fnames:
        ds = xr.open_dataset(fname)
        print("get_skinT_mean_anom: file = "+str(fname))

		# SST by the hour
        if cnt == 0:
          mask = np.squeeze(ds['land'][0,:,:])
          t_sst = np.squeeze(ds[varname][0,:,:]) - 273.15                # units = degrees-C
          sstmask = np.where(mask==maskval, np.nan, t_sst)
          sst[0,:,:] = np.array(sstmask)
          del t_sst,mask,sstmask
        else:
          mask = np.squeeze(ds['land'][0,:,:])
          t_sst = np.squeeze(ds[varname][0,:,:]) - 273.15
          sstmask = np.where(mask==maskval, np.nan, t_sst)
          sst[int(cnt),:,:] = np.array(sstmask)
          del t_sst,mask,sstmask

        cnt +=1
        ds.close()
    print("get_skinT_mean_anom: count = "+str(cnt))
    print("get_skinT_mean_anom: min/max sst = "+str(np.nanmin(sst))+" "+str(np.nanmax(sst)))

    return daily_mean, sst, sunlon, sunlat, lon, lat



#===============================================================================================
