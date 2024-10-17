###################################################################################################################
###################################################################################################################
# Analysis and Plot of SST Diurnal Cycle
###################################################################################################################
###################################################################################################################
#
# Output: Figure 
#
# How to run: python NameOfScript.py yyyymmdd yyyymmdd hh
#	First yyyymmdd  = START DATE
#	Second yyyymmdd = END DATE
#       hh              = analysis hour (00, 06, 12, or 18)
#
###################################################################################################################
# Import python modules
###################################################################################################################
print("***** Import Python Modules *****")

#import xarray as xr
import os
import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import datetime as dt
import sys
from scipy.interpolate import interpn
from scipy.interpolate import RectBivariateSpline

from read_data import read_diags_ocean_gridded
from read_data import get_dailymean_eachfile

from tools_analysis import get_ocngrid
from tools_analysis import regrid

from tools_plotting import plot_map

###################################################################################################################
print("========== BEGIN MAIN PROGRAM ==========")

print("START: "+str(dt.datetime.now()))

#=============================================
# Read raw user input from command line

yyyymmddhhS = sys.argv[1]         # Start date: year month day
yyyymmddhhE = sys.argv[2]         # End date: year month day
outpath     = sys.argv[3]
exptname0   = sys.argv[4]
texptpath0  = sys.argv[5]
texptpath0_diags  = sys.argv[6]
exptname1   = sys.argv[7]
texptpath1  = sys.argv[8]
texptpath1_diags  = sys.argv[9]
var0        = sys.argv[10]
var1        = sys.argv[11]
degree      = sys.argv[12]
meansumopt  = sys.argv[13]
hour        = str(sys.argv[14])            # analysis hour (00, 06, 12, or 18 UTC)
obstype     = str(sys.argv[15])
filetype    = sys.argv[16]

print("obstype = "+str(obstype))
print("exptname0 = "+str(exptname0))
print("exptname1 = "+str(exptname1))
print("outpath = "+str(outpath))

#=============================================
# Set global parameters

fill = -999.0

if filetype=="all_fcsts":
    plotstr = "bkg"
elif filetype=="f06":
    plotstr = "bkgf06"
elif filetype=="anl":
    plotstr = "anl"

#gridsize = 0.5	# units: degrees
#gridsize = 1	# units: degrees
#gridsize = 1.25	# units: degrees
#gridsize = 2
gridsize = int(degree)
flats = list(np.arange(-90,91,gridsize))
flons = list(np.arange(-180,181,gridsize))
gridtype_out="1x1"

xnew = flons            # new lon to conform to
ynew = flats            # new lat to conform to

xstr = "Longitude"
ystr = "Latitude"

#=============================================

# start date
yyyyS = yyyymmddhhS[0:4]
mmS   = yyyymmddhhS[4:6]
ddS   = yyyymmddhhS[6:8]
hhS   = yyyymmddhhS[8:10]
# end date
yyyyE = yyyymmddhhE[0:4]
mmE   = yyyymmddhhE[4:6]
ddE   = yyyymmddhhE[6:8]
hhE   = yyyymmddhhE[8:10]

#start_date = dt.date(int(yyyyS), int(mmS), int(ddS))
#end_date = dt.date(int(yyyyE), int(mmE), int(ddE))
start_date = yyyymmddhhS
end_date = yyyymmddhhE
print("start_date = "+str(start_date))
print("end_date = "+str(end_date))

#---------------------------------------------------------------
# Get data and compute daily means, anomalies

lat, lon = get_ocngrid()
shape_lat = np.shape(lat)
nlat = shape_lat[0]
nlon = shape_lat[1]
lat1d = lat[:,0]
lon1d = lon[0,:]
print("shape lat = "+str(shape_lat)+" nlat = "+str(nlat)+" nlon = "+str(nlon))
#del lat,lon

avgopt = 0

lat = ynew
lon = xnew

    #-----------------------------------------------

mmARR = [ '01','02','03','04','05','06','07','08','09','10','11','12' ]
ddARR = [ "01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31" ]
hhARR = [ '00','06','12','18' ]

#=========================================================================================
# LOOP
#=========================================================================================
print("BEGIN LOOP")

idatef0=0
idatef1=0
idate0=0
idate1=0

yyyy = yyyyS
iyy  = int(yyyy)
iyyE = int(yyyyE)
while iyy <= int(yyyyE):
  if yyyyS==yyyyE:
    Simm = mmARR.index(mmS)
    Eimm = mmARR.index(mmE)
  elif yyyyS!=yyyyE:
    if iyy==int(yyyyS):
      Simm = mmARR.index(mmS)
      Eimm = 12-1
    elif iyy>int(yyyyS) and iyy<iyyE:
      Simm = 0
      Eimm = 12-1
    elif iyy==iyyE:
      Simm = 0
      Eimm = mmARR.index(mmE)

  imm=Simm
  while imm <= Eimm:
    mm = mmARR[imm]
    if mm=="02":
      if int(yyyy)%4==0:
        # leap year
        nDD = 29
      else:
        nDD = 28
    elif mm=="04" or mm=="06" or mm=="09" or mm=="11":
      nDD = 30
    else:
      nDD = 31

    if yyyyS==yyyyE:
      if mmS==mmE:
        Sidd = ddARR.index(ddS)
        Eidd = ddARR.index(ddE)
      else:
        if mm==mmS:
          Sidd = ddARR.index(ddS)
          Eidd = nDD-1
        elif mm==mmE:
          Sidd = 0
          Eidd = ddARR.index(ddE)
        else:
          Sidd = 0
          Eidd = nDD-1
    elif yyyyS!=yyyyE:
      if yyyy==yyyyS and mm==mmS:
        Sidd = ddARR.index(ddS)
        Eidd = nDD-1
      elif iyy==int(yyyyE) and mm==mmE:
        Sidd = 0
        Eidd = ddARR.index(ddE)
      else:
        Sidd = 0
        Eidd = nDD-1

    idd=Sidd
    while idd <= Eidd:
        dd = ddARR[idd]
        print("----- DATE: "+str(yyyy)+"-"+str(mm)+"-"+str(dd)+" -----")

        dd_str = str(dd).zfill(2)		# add zeroes to beginning of string
	
        opt = -1        # all points (do not stratify into day/night)

        #-------------------------------------------------------------------
	# EXPT 0
        print("..... EXPERIMENT 0 .....")
        #if texptpath0.find("role")!=-1: 
        #    exptpath0 = texptpath0+"/"+str(yyyy)+str(mm)+str(dd_str)+str(hour)+"/"
        #else:
        exptpath0 = texptpath0
        print("exptpath0 = "+str(exptpath0))

        #````````````````````````````````````````````
        # Compute diurnal cycle in bkg, anl, or inc

            # get var0 from expt 0
        tvar = var0
        if tvar=="sstavg":
          ocn_mean, ocn, ocnlat, ocnlon, skin_mean, skin, skinlat, skinlon = get_dailymean_eachfile(exptpath0, yyyy, mm, dd_str, hour, filetype, tvar)

          # regrid skin to ocean model grid
          lat_in  = skinlat
          lon_in  = skinlon
          lat_out = lat     # from get_ocngrd
          lon_out = lon     # from get_ocngrd
          print("shape ocnlat = "+str(np.shape(ocnlat))+" skinlat = "+str(np.shape(skinlat)))

          print("skin_mean")
          var_in = skin_mean
          gridtype_in = "atmos model"
          skin_mean_regrid = regrid(var_in, lat_in, lon_in, lat_out, lon_out, gridtype_in, gridtype_out)
          del var_in
          print("skin")
          var_in = skin[0,:,:]
          skin_regrid = regrid(var_in, lat_in, lon_in, lat_out, lon_out, gridtype_in, gridtype_out)
          del var_in

          # average ocean SST + skin T
          sst_mean0 = 0.5*(ocn_mean + skin_mean_regrid)
          sst0 = 0.5*(ocn + skin_regrid)

          lat_in_anom0 = lat1d
          lon_in_anom0 = lon1d

          del ocn_mean, ocn, ocnlat, ocnlon, skin_mean, skin, skinlat, skinlon
          del lat_in, lon_in, lat_out, lon_out
          del skin_mean_regrid, skin_regrid

        else:
          sst_mean0, sst0, slat, slon = get_dailymean_eachfile(exptpath0, yyyy, mm, dd_str, hour, filetype, tvar)
          #if tvar=="skint":
          #  del lat, lon
          lat_in_anom0 = slat
          lon_in_anom0 = slon
          del slat, slon

            # get var1 from expt 0
        tvar = var1
        if tvar=="sstavg":
          ocn_mean, ocn, ocnlat, ocnlon, skin_mean, skin, skinlat, skinlon = get_dailymean_eachfile(exptpath0, yyyy, mm, dd_str, hour, filetype, tvar)

          # regrid skin to ocean model grid
          lat_in  = skinlat
          lon_in  = skinlon
          lat_out = lat     # from get_ocngrd
          lon_out = lon     # from get_ocngrd
          print("shape ocnlat = "+str(np.shape(ocnlat))+" skinlat = "+str(np.shape(skinlat)))

          print("skin_mean")
          var_in = skin_mean
          gridtype_in = "atmos model"
          skin_mean_regrid = regrid(var_in, lat_in, lon_in, lat_out, lon_out, gridtype_in, gridtype_out)
          del var_in
          print("skin")
          var_in = skin[0,:,:]
          skin_regrid = regrid(var_in, lat_in, lon_in, lat_out, lon_out, gridtype_in, gridtype_out)
          del var_in

          # average ocean SST + skin T on ocean model (SST) grid
          sst_mean1 = 0.5*(ocn_mean + skin_mean_regrid)
          sst1 = 0.5*(ocn + skin_regrid)

          lat_in_anom1 = lat1d
          lon_in_anom1 = lon1d

          del ocn_mean, ocn, ocnlat, ocnlon, skin_mean, skin, skinlat, skinlon
          del lat_in, lon_in, lat_out, lon_out
          del skin_mean_regrid, skin_regrid

        else:
          sst_mean1, sst1, slat, slon = get_dailymean_eachfile(exptpath0, yyyy, mm, dd_str, hour, filetype, tvar)
          #if tvar=="skint":
          #  del lat, lon
          lat_in_anom1 = slat
          lon_in_anom1 = slon
          del slat, slon

        #print("shape sst0 = "+str(np.shape(sst0))+" sst_mean0 = "+str(np.shape(sst_mean0)))
        if np.size(sst0)>1 and np.size(sst1)>1:
        #  del sst_mean0, sst0, slat, slon
        #  idd += 1
        #  continue

          # SST anomalies (deviations from daily mean)
          sst_anom0 = sst0 - sst_mean0
          sst_anom1 = sst1 - sst_mean1

          # regrid to 1x1 degree grid: VAR0
          print("VAR 0 from expt 0")
          if var0=="sst" or var0=="sstavg":
              gridtype_in = "ocean model"
          elif var0=="skint":
              gridtype_in = "atmos model"
          lon_out_anom, lat_out_anom = np.meshgrid(flons,flats)
          sst_anom0_regrid = regrid(sst_anom0[0,:,:], lat_in_anom0, lon_in_anom0, lat_out_anom, lon_out_anom, gridtype_in, gridtype_out)
          del sst_anom0
          del lat_in_anom0, lon_in_anom0#, lat_out_anom, lon_out_anom
          print("shape sst_anom0_regrid = "+str(np.shape(sst_anom0_regrid)))

          tshape = np.shape(sst_anom0_regrid)
          sst_anom0_regrid3d = np.empty((1,tshape[0],tshape[1]),float)
          sst_anom0_regrid3d[0,:,:] = sst_anom0_regrid[:,:]
          del tshape

          if idatef0==0:
            sst_anom_arr0 = sst_anom0_regrid3d
          elif idatef0>0:
            sst_anom_arr0 = np.append(sst_anom_arr0, sst_anom0_regrid3d, axis=0)
          idatef0 += 1

          del sst_anom0_regrid, sst_anom0_regrid3d

          # regrid to 1x1 degree grid: VAR1
          print("VAR 1 from expt 0")
          if var1=="sst" or var1=="sstavg":
              gridtype_in = "ocean model"
          elif var1=="skint":
              gridtype_in = "atmos model"
          sst_anom1_regrid = regrid(sst_anom1[0,:,:], lat_in_anom1, lon_in_anom1, lat_out_anom, lon_out_anom, gridtype_in, gridtype_out)
          del sst_anom1
          del lat_in_anom1, lon_in_anom1
          del lat_out_anom, lon_out_anom
          print("shape sst_anom1_regrid = "+str(np.shape(sst_anom1_regrid)))

          tshape = np.shape(sst_anom1_regrid)
          sst_anom1_regrid3d = np.empty((1,tshape[0],tshape[1]),float)
          sst_anom1_regrid3d[0,:,:] = sst_anom1_regrid[:,:]
          del tshape

          if idatef1==0:
            sst_anom_arr1 = sst_anom1_regrid3d
          elif idatef1>0:
            sst_anom_arr1 = np.append(sst_anom_arr1, sst_anom1_regrid3d, axis=0)
          idatef1 += 1

          del sst_anom1_regrid, sst_anom1_regrid3d

        del sst_mean0, sst0
        del sst_mean1, sst1

        #-------------------------------------------------------------------
        # EXPT 1
        #if exptname1!=exptname0:
        print("..... EXPERIMENT 1 .....")
        exptpath1 = texptpath1
        print("exptpath1 = "+str(exptpath1))

          #```````````````````````````````````````````
          # Compute diurnal cycle in the obs/diags
        exptpath1 = texptpath1_diags
        print("exptpath1 = "+str(exptpath1))

        nobs, omb, SDomb, oma, SDoma, obserr, SDobserr, obserr_orig, SDobserr_orig, obsval, SDobsval, bkg, SDbkg, anl, SDanl = read_diags_ocean_gridded(exptpath1, yyyy, mm, dd_str, hour, exptname1, obstype, degree, meansumopt)
        ombHH = omb
        omaHH = oma
        obsHH = obsval
        bkgHH = bkg
        anlHH = anl
        nobsHH = nobs
        del omb, oma, obsval, bkg, anl, nobs
        del SDomb, SDoma, obserr, SDobserr, obserr_orig, SDobserr_orig, SDobsval, SDbkg, SDanl

        if np.size(ombHH)>1:

            # Compute daily mean
            print("compute daily means")
            print("... get 00z")
            nobs, omb, SDomb, oma, SDoma, obserr, SDobserr, obserr_orig, SDobserr_orig, obsval, SDobsval, bkg, SDbkg, anl, SDanl = read_diags_ocean_gridded(exptpath1, yyyy, mm, dd_str, "00", exptname1, obstype, degree, meansumopt)
            omb00 = omb
            oma00 = oma
            obs00 = obsval
            bkg00 = bkg
            anl00 = anl
            nobs00 = nobs
            del omb, oma, obsval, bkg, anl, nobs
            del SDomb, SDoma, obserr, SDobserr, obserr_orig, SDobserr_orig, SDobsval, SDbkg, SDanl

            print("... get 06z")
            nobs, omb, SDomb, oma, SDoma, obserr, SDobserr, obserr_orig, SDobserr_orig, obsval, SDobsval, bkg, SDbkg, anl, SDanl = read_diags_ocean_gridded(exptpath1, yyyy, mm, dd_str, "06", exptname1, obstype, degree, meansumopt)
            omb06 = omb
            oma06 = oma
            obs06 = obsval
            bkg06 = bkg
            anl06 = anl
            nobs06 = nobs
            del omb, oma, obsval, bkg, anl, nobs
            del SDomb, SDoma, obserr, SDobserr, obserr_orig, SDobserr_orig, SDobsval, SDbkg, SDanl

            print("... get 12z")
            nobs, omb, SDomb, oma, SDoma, obserr, SDobserr, obserr_orig, SDobserr_orig, obsval, SDobsval, bkg, SDbkg, anl, SDanl = read_diags_ocean_gridded(exptpath1, yyyy, mm, dd_str, "12", exptname1, obstype, degree, meansumopt)
            omb12 = omb
            oma12 = oma
            obs12 = obsval
            bkg12 = bkg
            anl12 = anl
            nobs12 = nobs
            del omb, oma, obsval, bkg, anl, nobs
            del SDomb, SDoma, obserr, SDobserr, obserr_orig, SDobserr_orig, SDobsval, SDbkg, SDanl

            print("... get 18z")
            nobs, omb, SDomb, oma, SDoma, obserr, SDobserr, obserr_orig, SDobserr_orig, obsval, SDobsval, bkg, SDbkg, anl, SDanl = read_diags_ocean_gridded(exptpath1, yyyy, mm, dd_str, "18", exptname1, obstype, degree, meansumopt)
            omb18 = omb
            oma18 = oma
            obs18 = obsval
            bkg18 = bkg
            anl18 = anl
            nobs18 = nobs
            del omb, oma, obsval, bkg, anl, nobs
            del SDomb, SDoma, obserr, SDobserr, obserr_orig, SDobserr_orig, SDobsval, SDbkg, SDanl

            ndaily=0
            if np.size(omb00)>1: ndaily += 1
            if np.size(omb06)>1: ndaily += 1
            if np.size(omb12)>1: ndaily += 1
            if np.size(omb18)>1: ndaily += 1
            print("ndaily = "+str(ndaily))

            # set missing points to NaN
            if np.size(omb00)>1:
              omb00 = np.where(omb00>999, np.nan, omb00)
              oma00 = np.where(oma00>999, np.nan, oma00)
              obs00 = np.where(obs00>999, np.nan, obs00)
              bkg00 = np.where(bkg00>999, np.nan, bkg00)
              anl00 = np.where(anl00>999, np.nan, anl00)
              nobs00 = np.where(nobs00>999, np.nan, nobs00)
            else:
              omb00 = 0
              oma00 = 0
              obs00 = 0
              bkg00 = 0
              anl00 = 0
              nobs00 = 0

            if np.size(omb06)>1:
              omb06 = np.where(omb06>999, np.nan, omb06)
              oma06 = np.where(oma06>999, np.nan, oma06)
              obs06 = np.where(obs06>999, np.nan, obs06)
              bkg06 = np.where(bkg06>999, np.nan, bkg06)
              anl06 = np.where(anl06>999, np.nan, anl06)
              nobs06 = np.where(nobs06>999, np.nan, nobs06)
            else:
              omb06 = 0
              oma06 = 0
              obs06 = 0
              bkg06 = 0
              anl06 = 0
              nobs06 = 0

            if np.size(omb12)>1:
              omb12 = np.where(omb12>999, np.nan, omb12)
              oma12 = np.where(oma12>999, np.nan, oma12)
              obs12 = np.where(obs12>999, np.nan, obs12)
              bkg12 = np.where(bkg12>999, np.nan, bkg12)
              anl12 = np.where(anl12>999, np.nan, anl12)
              nobs12 = np.where(nobs12>999, np.nan, nobs12)
            else:
              omb12 = 0
              oma12 = 0
              obs12 = 0
              bkg12 = 0
              anl12 = 0
              nobs12 = 0

            if np.size(omb18)>1:
              omb18 = np.where(omb18>999, np.nan, omb18)
              oma18 = np.where(oma18>999, np.nan, oma18)
              obs18 = np.where(obs18>999, np.nan, obs18)
              bkg18 = np.where(bkg18>999, np.nan, bkg18)
              anl18 = np.where(anl18>999, np.nan, anl18)
              nobs18 = np.where(nobs18>999, np.nan, nobs18)
            else:
              omb18 = 0
              oma18 = 0
              obs18 = 0
              bkg18 = 0
              anl18 = 0
              nobs18 = 0

            ombHH = np.where(ombHH>999, np.nan, ombHH)
            omaHH = np.where(omaHH>999, np.nan, omaHH)
            obsHH = np.where(obsHH>999, np.nan, obsHH)
            bkgHH = np.where(bkgHH>999, np.nan, bkgHH)
            anlHH = np.where(anlHH>999, np.nan, anlHH)

            # compute averages
            sum_nobs = nobs00 + nobs06 + nobs12 + nobs18
            del nobs00, nobs06, nobs12, nobs18

            if int(meansumopt)==0:
              # means
              mean_omb = (omb00 + omb06 + omb12 + omb18)/(1.0*ndaily)
              mean_oma = (oma00 + oma06 + oma12 + oma18)/(1.0*ndaily)
              mean_obs = (obs00 + obs06 + obs12 + obs18)/(1.0*ndaily)
              mean_bkg = (bkg00 + bkg06 + bkg12 + bkg18)/(1.0*ndaily)
              mean_anl = (anl00 + anl06 + anl12 + anl18)/(1.0*ndaily)
            elif int(meansumopt)==1:
              # sums
              mean_omb = (omb00 + omb06 + omb12 + omb18)/(1.0*sum_nobs) #float(ndaily)
              mean_oma = (oma00 + oma06 + oma12 + oma18)/(1.0*sum_nobs) #float(ndaily)
              mean_obs = (obs00 + obs06 + obs12 + obs18)/(1.0*sum_nobs) #float(ndaily)
              mean_bkg = (bkg00 + bkg06 + bkg12 + bkg18)/(1.0*sum_nobs) #float(ndaily)
              mean_anl = (anl00 + anl06 + anl12 + anl18)/(1.0*sum_nobs) #float(ndaily)
            del omb00, omb06, omb12, omb18
            del oma00, oma06, oma12, oma18
            del obs00, obs06, obs12, obs18
            del bkg00, bkg06, bkg12, bkg18
            del anl00, anl06, anl12, anl18
            del sum_nobs
            del ndaily

            # Compute anomaly (for diurnal cycle)
            print("anomalies")
            shapeomb = np.shape(ombHH)
            anom_omb        = np.empty((1,shapeomb[0],shapeomb[1]),float)
            anom_oma        = np.empty((1,shapeomb[0],shapeomb[1]),float)
            anom_obs        = np.empty((1,shapeomb[0],shapeomb[1]),float)
            anom_bkg        = np.empty((1,shapeomb[0],shapeomb[1]),float)
            anom_anl        = np.empty((1,shapeomb[0],shapeomb[1]),float)
            sum_nobs3d      = np.empty((1,shapeomb[0],shapeomb[1]),float)
            anom_omb[0,:,:] = ombHH - mean_omb
            anom_oma[0,:,:] = omaHH - mean_oma
            anom_obs[0,:,:] = obsHH - mean_obs
            anom_bkg[0,:,:] = bkgHH - mean_bkg
            anom_anl[0,:,:] = anlHH - mean_anl
            sum_nobs3d[0,:,:] = nobsHH
            del mean_omb, mean_oma, mean_obs, mean_bkg, mean_anl
            del ombHH, omaHH, obsHH, bkgHH, anlHH
            del shapeomb

            # Append to total arrays
            print("append to total arrays")
            if idate1==0:
              print("idate1 = 0")
              omb1_arr = anom_omb
              oma1_arr = anom_oma
              obs1_arr = anom_obs
              bkg1_arr = anom_bkg
              anl1_arr = anom_anl
              nobs1_arr = sum_nobs3d
              idate1 += 1
            elif idate1>0:
              print("idate1 > 0")
              print("shape omb1_arr = "+str(np.shape(omb1_arr))+" anom_omb = "+str(np.shape(anom_omb)))
              omb1_arr = np.append(omb1_arr, anom_omb, axis=0)
              oma1_arr = np.append(oma1_arr, anom_oma, axis=0)
              obs1_arr = np.append(obs1_arr, anom_obs, axis=0)
              bkg1_arr = np.append(bkg1_arr, anom_bkg, axis=0)
              anl1_arr = np.append(anl1_arr, anom_anl, axis=0)
              nobs1_arr = np.append(nobs1_arr, sum_nobs3d, axis=0)
              idate1 += 1
            del anom_omb, anom_oma, anom_obs, anom_bkg, anom_anl
            del sum_nobs3d
        else:
            del ombHH, omaHH, obsHH, bkgHH, anlHH

        #-------------------------------------------------------------------

        idd += 1
    imm += 1
  iyy += 1

print("END LOOP")
#=========================================================================================

print("idate0 = "+str(idate0)+" idate1 = "+str(idate1))
print("idatef0 = "+str(idatef0)+" idatef1 = "+str(idatef1))

        #---------------------------------
	# get 2D grid
del lat,lon
lon, lat = np.meshgrid(xnew, ynew)

	#---------------------------------
	# compute mean anomalies and their differences
if idatef0>0:
    # var0 from expt 0
  sst_anom_mean0 = np.nanmean(sst_anom_arr0, axis=0)
  #del sst_anom_arr0

if idatef1>0:
    # var1 from expt 0
  sst_anom_mean1 = np.nanmean(sst_anom_arr1, axis=0)
  #del sst_anom_arr1

if idate1>0:
      # mean anomaly
  obs_anom_mean1 = np.nanmean(obs1_arr, axis=0)
  nobs_total_sum1 = np.nansum(nobs1_arr, axis=0)

  #del omb1_arr, oma1_arr, obs1_arr, bkg1_arr, anl1_arr, nobs1_arr

# absolute differences of SST diurnal cycle
if idatef0>0 and idatef1>0 and idate1>0:
  sst0_obs1_diff    = sst_anom_mean0 - obs_anom_mean1
  sst0_obs1_absdiff = abs(sst_anom_mean0) - abs(obs_anom_mean1)
  sst1_obs1_diff    = sst_anom_mean1 - obs_anom_mean1
  sst1_obs1_absdiff = abs(sst_anom_mean1) - abs(obs_anom_mean1)


#*****************************************************************************************
#*****************************************************************************************
# PLOTS
#*****************************************************************************************
#*****************************************************************************************

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# Paths and Figure names

    # edit outpath
if yyyymmddhhS==yyyymmddhhE:
  daterange     = str(yyyyS)+str(mmS)+str(ddS)
else:
  daterange     = str(yyyyS)+str(mmS)+str(ddS)+"-"+str(yyyyE)+str(mmE)+str(ddE)     # study period plotted
datepath      = daterange+"/"
hourpath      = "hourly/"

varpath = "sst/" #var+"/"
if obstype=="all":
  geopath = "global/"
else:
  geopath = obstype+"/"

outpath = outpath+datepath+hourpath+varpath+geopath

os.system('mkdir -m 775 -p '+outpath)

    # figure name parts
figtype         = "DIURNAL_DIAGS_DIFF"                                # type of plot
#figsubtype      = "OMBA"                       # type of info plotted
figexpt0        = str(exptname0)                         # experiment name
figexpt1        = str(exptname1)                         # experiment name
figobs          = str(obstype)                          # observation instrument

#figfiletype = str(plotstr)

figext          = "png"                                 # figure extension: examples: .png, .jpg, .pdf, .gif

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

#===============================================================
# output stats in text file

statstr = [var0+'_Anomaly', var0+'_SD_Anomaly', var1+'_Anomaly', var1+'_SD_Anomaly', 'Obs_Anomaly', 'Obs_SD_Anomaly']

if idatef0>0 and idatef1>0:
  figexpt = figexpt0
  i=0
  for istat in statstr:
    if istat.find(var0)!=-1:
      stat = sst_anom_arr0
      obs  = obs1_arr
      var  = istat
    elif istat.find(var1)!=-1:
      stat = sst_anom_arr1
      obs  = obs1_arr
      var  = istat
    elif istat.find('Obs')!=-1:
      stat  = obs1_arr
      count = nobs_total_sum1
      var   = "sst_"+istat
      figexpt = figexpt1

    title = 'DIURNAL CYCLE STATS: '+str(var)+' ('+figexpt+'): '+str(hour)+'z '+daterange
    outname = outpath+"STATS_DIURNAL_CYCLE."+figexpt+"."+str(hour)+"z."+str(var.lower())

        # compute and print stats of differences
    avg  = np.nanmean(stat)                  # mean
    sd   = np.nanstd(stat)                   # standard deviation

    if istat.find(var0)!=-1 or istat.find(var1)!=-1:
      diff = stat - obs
      rmsd = np.sqrt(np.nanmean(diff**2))       # RMSD
      del diff

                # ADD TEXT INSIDE PANEL
                # add title for difference stats
    text = title
                # count
    if istat.find('Obs')!=-1:
      textc = "... Count = "+str(np.sum(np.where(~np.isnan(stat))))
                # add mean diff
    textm = "... Mean = "+str(np.round_(avg,decimals=3))
                # add SD of diff
    texts = "... StdDev = "+str(np.round_(sd,decimals=3))
                # add RMSD
    if istat.find(var0)!=-1 or istat.find(var1)!=-1:
      textrmsd = "... RMSD ("+var+" - Obs) = "+str(np.round_(rmsd,decimals=3))

                # OUTPUT STATS TO TEXT FILE
    os.system("rm -rf "+outname+".txt")
    file1 = open(outname+".txt","a+")            # open file for reading and writing. create file if it does not exist. data is appended to end of file.
    file1.write(text+"\n")
    if istat.find('Obs')!=-1: file1.write(textc+"\n")
    file1.write(textm+"\n")
    file1.write(texts+"\n")
    if istat.find(var0)!=-1 or istat.find(var1)!=-1: file1.write(textrmsd+"\n")

    file1.close()
    del stat,avg,sd,text,textm,texts
    del title, outname

    i += 1

#===============================================================

if hour=='00':
    sunlat = 20
    sunlon = -170
elif hour=='06':
    sunlat = 20
    sunlon = 100
elif hour=='12':
    sunlat = 25
    sunlon = 5
elif hour=='18':
    sunlat = 25
    sunlon = -85

zmean0 = "none"
#zmean0 = sst_mean_mean0
if exptname0!=exptname1:
  zmean1 = "none"
#  zmean1 = sst_mean_mean1

#----------------------------------------------
# Plot diurnal cycle differences (SST from bkg/anl/inc minus obs)

print("outpath = "+str(outpath))

pltmin = -0.5
pltmax = 0.5
plottype = "Diff"       # for colorbar

# plot diurnal cycle (daily mean anomaly)
print("plot daily mean anomalies")

if idatef0>0:
        # var0 from expt 0
  z0      = sst_anom_mean0
  var     = var0
  figexpt = figexpt0

    # compute stats for plotting
  figmean  = np.round_(np.nanmean(z0),decimals=3)
  figsd    = np.round_(np.nanstd(z0),decimals=3)

  title = 'Mean Anomaly of '+var.upper()+' '+plotstr.upper()+' ('+figexpt+'): \n'+str(hour)+'z '+daterange+'\n Mean = '+str(figmean)+' | Stddev = '+str(figsd)
  plot_map(lon, lat, z0, zmean0, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figexpt+"."+str(hour)+"z."+var+"_"+plotstr+".anom_mean."+figext
  plt.savefig(outname)
  del z0
  del figmean, figsd

if idatef1>0:
        # var1 from expt 0
  z0      = sst_anom_mean1
  var     = var1
  figexpt = figexpt0

    # compute stats for plotting
  figmean  = np.round_(np.nanmean(z0),decimals=3)
  figsd    = np.round_(np.nanstd(z0),decimals=3)

  title = 'Mean Anomaly of '+var.upper()+' '+plotstr.upper()+' ('+figexpt+'): \n'+str(hour)+'z '+daterange+'\n Mean = '+str(figmean)+' | Stddev = '+str(figsd)
  plot_map(lon, lat, z0, zmean0, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figexpt+"."+str(hour)+"z."+var+"_"+plotstr+".anom_mean."+figext
  plt.savefig(outname)
  del z0
  del figmean, figsd

if idate1>0:
        # obs from expt 1
  z0      = obs_anom_mean1
  z0count = nobs_total_sum1
  var     = "sst"
  figexpt = figexpt1

    # compute stats for plotting
  figmean  = np.round_(np.nanmean(z0),decimals=3)
  figsd    = np.round_(np.nanstd(z0),decimals=3)
  figcount = np.sum(np.where(~np.isnan(z0count)))

  title = 'Mean Anomaly of '+var.upper()+' Obs ('+figexpt+'): \n'+str(hour)+'z '+daterange+'\n Mean = '+str(figmean)+' | Stddev = '+str(figsd)+' | Count = '+str(figcount)
  plot_map(lon, lat, z0, zmean0, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figexpt+"."+str(hour)+"z."+var+"_obs.anom_mean."+figext
  plt.savefig(outname)
  del z0
  del figmean, figsd, figcount

# plot diurnal cycle DIFFERENCES
print("plot diffs in diurnal cycles")

if idatef0>0 and idate1>0:
  z0      = sst0_obs1_absdiff
  var     = var0
  figexpt = figexpt0

    # compute stats for plotting
  figmean = np.round_(np.nanmean(z0),decimals=3)
  figsd   = np.round_(np.nanstd(z0),decimals=3)
  diff = z0
  rmsd = np.sqrt(np.nanmean(diff**2))       # RMSD
  del diff
  figrmsd = np.round_(rmsd,decimals=3)
  del rmsd

  title = 'Abs. Diff in Mean Anomalies of '+var.upper()+' for '+figexpt+': \n|'+plotstr.upper()+' - Obs '+figobs+'| \n'+str(hour)+'z '+daterange+'\n Mean = '+str(figmean)+' | RMSD = '+str(figrmsd)
  plot_map(lon, lat, z0, zmean0, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+str(hour)+"z."+figexpt+"_"+var+"_anom_mean."+plotstr+"_MINUS_obs_"+figobs+".absdiff."+figext
  plt.savefig(outname)
  del z0
  del figmean, figsd, figrmsd

if idatef1>0 and idate1>0:
  z0      = sst1_obs1_absdiff
  var     = var1
  figexpt = figexpt0

    # compute stats for plotting
  figmean = np.round_(np.nanmean(z0),decimals=3)
  figsd   = np.round_(np.nanstd(z0),decimals=3)
  diff = z0
  rmsd = np.sqrt(np.nanmean(diff**2))       # RMSD
  del diff
  figrmsd = np.round_(rmsd,decimals=3)
  del rmsd

  title = 'Abs. Diff in Mean Anomalies of '+var.upper()+' for '+figexpt+': \n|'+plotstr.upper()+' - Obs '+figobs+'| \n'+str(hour)+'z '+daterange+'\n Mean = '+str(figmean)+' | RMSD = '+str(figrmsd)
  plot_map(lon, lat, z0, zmean0, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+str(hour)+"z."+figexpt+"_"+var+"_anom_mean."+plotstr+"_MINUS_obs_"+figobs+".absdiff."+figext
  plt.savefig(outname)
  del z0
  del figmean, figsd, figrmsd

if idatef0>0 and idatef1>0 and idate1>0:
  z0      = sst1_obs1_absdiff - sst0_obs1_absdiff
  figexpt = figexpt0

    # compute stats for plotting
  figmean = np.round_(np.nanmean(z0),decimals=3)
  figsd   = np.round_(np.nanstd(z0),decimals=3)
  diff = z0
  rmsd = np.sqrt(np.nanmean(diff**2))       # RMSD
  del diff
  figrmsd = np.round_(rmsd,decimals=3)
  del rmsd

  title = 'Diff in Abs. Diff of Mean Anomalies of '+var1.upper()+'-'+var0.upper()+'for '+figexpt+': \n|'+str(var1.upper())+' '+plotstr.upper()+' - Obs '+figobs+'| - |'+str(var0.upper())+' '+plotstr.upper()+' - Obs '+figobs+'| \n'+str(hour)+'z '+daterange+'\n Mean = '+str(figmean)+' | RMSD = '+str(figrmsd)
  plot_map(lon, lat, z0, zmean0, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+str(hour)+"z."+figexpt+"."+var1+"_anom_mean__MINUS__"+var0+"_anom_mean."+plotstr+"_MINUS_obs_"+figobs+".absdiff."+figext
  plt.savefig(outname)
  del z0
  del figmean, figsd, figrmsd

plt.close("all")

###################################################################################################################
print("END: "+str(dt.datetime.now()))
print("========== END MAIN PROGRAM ==========")
###################################################################################################################
###################################################################################################################
