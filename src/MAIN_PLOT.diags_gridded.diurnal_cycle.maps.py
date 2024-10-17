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

from tools_analysis import get_ocngrid

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
exptname1   = sys.argv[6]
texptpath1  = sys.argv[7]
var         = sys.argv[8]
degree      = sys.argv[9]
meansumopt  = sys.argv[10]
hour        = str(sys.argv[11])            # analysis hour (00, 06, 12, or 18 UTC)
obstype     = str(sys.argv[12])

print("obstype = "+str(obstype))
print("exptname0 = "+str(exptname0))
print("exptname1 = "+str(exptname1))
print("outpath = "+str(outpath))

#=============================================
# Set global parameters

fill = -999.0

#gridsize = 0.5	# units: degrees
#gridsize = 1	# units: degrees
#gridsize = 1.25	# units: degrees
#gridsize = 2
gridsize = int(degree)
flats = list(np.arange(-90,91,gridsize))
flons = list(np.arange(-180,181,gridsize))

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
        # Compute diurnal cycle in the obs and diags
        nobs, omb, SDomb, oma, SDoma, obserr, SDobserr, obserr_orig, SDobserr_orig, obsval, SDobsval, bkg, SDbkg, anl, SDanl = read_diags_ocean_gridded(exptpath0, yyyy, mm, dd_str, hour, exptname0, obstype, degree, meansumopt)
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
          nobs, omb, SDomb, oma, SDoma, obserr, SDobserr, obserr_orig, SDobserr_orig, obsval, SDobsval, bkg, SDbkg, anl, SDanl = read_diags_ocean_gridded(exptpath0, yyyy, mm, dd_str, "00", exptname0, obstype, degree, meansumopt)
          omb00 = omb
          oma00 = oma
          obs00 = obsval
          bkg00 = bkg
          anl00 = anl
          nobs00 = nobs
          del omb, oma, obsval, bkg, anl, nobs
          del SDomb, SDoma, obserr, SDobserr, obserr_orig, SDobserr_orig, SDobsval, SDbkg, SDanl

          print("... get 06z")
          nobs, omb, SDomb, oma, SDoma, obserr, SDobserr, obserr_orig, SDobserr_orig, obsval, SDobsval, bkg, SDbkg, anl, SDanl = read_diags_ocean_gridded(exptpath0, yyyy, mm, dd_str, "06", exptname0, obstype, degree, meansumopt)
          omb06 = omb
          oma06 = oma
          obs06 = obsval
          bkg06 = bkg
          anl06 = anl
          nobs06 = nobs
          del omb, oma, obsval, bkg, anl, nobs
          del SDomb, SDoma, obserr, SDobserr, obserr_orig, SDobserr_orig, SDobsval, SDbkg, SDanl

          print("... get 12z")
          nobs, omb, SDomb, oma, SDoma, obserr, SDobserr, obserr_orig, SDobserr_orig, obsval, SDobsval, bkg, SDbkg, anl, SDanl = read_diags_ocean_gridded(exptpath0, yyyy, mm, dd_str, "12", exptname0, obstype, degree, meansumopt)
          omb12 = omb
          oma12 = oma
          obs12 = obsval
          bkg12 = bkg
          anl12 = anl
          nobs12 = nobs
          del omb, oma, obsval, bkg, anl, nobs
          del SDomb, SDoma, obserr, SDobserr, obserr_orig, SDobserr_orig, SDobsval, SDbkg, SDanl

          print("... get 18z")
          nobs, omb, SDomb, oma, SDoma, obserr, SDobserr, obserr_orig, SDobserr_orig, obsval, SDobsval, bkg, SDbkg, anl, SDanl = read_diags_ocean_gridded(exptpath0, yyyy, mm, dd_str, "18", exptname0, obstype, degree, meansumopt)
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
          print("shape mean_omb = "+str(np.shape(mean_omb)))
          print("min/max mean_omb = "+str(np.nanmin(mean_omb))+" "+str(np.nanmax(mean_omb)))
          print("min/max mean_oma = "+str(np.nanmin(mean_oma))+" "+str(np.nanmax(mean_oma)))	

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
          print("shape anom_omb = "+str(np.shape(anom_omb)))
          print("min/max anom_omb = "+str(np.min(anom_omb))+" "+str(np.max(anom_omb)))
          print("min/max anom_oma = "+str(np.min(anom_oma))+" "+str(np.max(anom_oma)))
          print("min/max anom_obs = "+str(np.min(anom_obs))+" "+str(np.max(anom_obs)))
          print("min/max anom_bkg = "+str(np.min(anom_bkg))+" "+str(np.max(anom_bkg)))
          print("min/max anom_anl = "+str(np.min(anom_anl))+" "+str(np.max(anom_anl)))
          print("min/max sum_nobs = "+str(np.min(sum_nobs3d))+" "+str(np.max(sum_nobs3d)))
          print("min/max anom_omb no NaN = "+str(np.nanmin(anom_omb))+" "+str(np.nanmax(anom_omb)))
          print("min/max anom_oma no NaN = "+str(np.nanmin(anom_oma))+" "+str(np.nanmax(anom_oma)))
          print("min/max anom_obs no NaN = "+str(np.nanmin(anom_obs))+" "+str(np.nanmax(anom_obs)))
          print("min/max anom_bkg no NaN = "+str(np.nanmin(anom_bkg))+" "+str(np.nanmax(anom_bkg)))
          print("min/max anom_anl no NaN = "+str(np.nanmin(anom_anl))+" "+str(np.nanmax(anom_anl)))
          print("min/max sum_nobs no NaN = "+str(np.nanmin(sum_nobs3d))+" "+str(np.nanmax(sum_nobs3d)))

          # Append to total arrays
          print("append to total arrays")
          if idate0==0:
            print("idate0 = 0")
            omb0_arr = anom_omb
            oma0_arr = anom_oma
            obs0_arr = anom_obs
            bkg0_arr = anom_bkg
            anl0_arr = anom_anl
            nobs0_arr = sum_nobs3d
            idate0 += 1
          elif idate0>0:
            print("idate0 > 0")
            print("shape omb0_arr = "+str(np.shape(omb0_arr))+" anom_omb = "+str(np.shape(anom_omb)))
            omb0_arr = np.append(omb0_arr, anom_omb, axis=0)
            oma0_arr = np.append(oma0_arr, anom_oma, axis=0)
            obs0_arr = np.append(obs0_arr, anom_obs, axis=0)
            bkg0_arr = np.append(bkg0_arr, anom_bkg, axis=0)
            anl0_arr = np.append(anl0_arr, anom_anl, axis=0)
            nobs0_arr = np.append(nobs0_arr, sum_nobs3d, axis=0)
            #omb0_arr += anom_omb
            #oma0_arr += anom_oma
            #obs0_arr += anom_obs
            #bkg0_arr += anom_bkg
            #anl0_arr += anom_anl
            idate0 += 1
          del anom_omb, anom_oma, anom_obs, anom_bkg, anom_anl	
          del sum_nobs3d
        else:
          del ombHH, omaHH, obsHH, bkgHH, anlHH

        #-------------------------------------------------------------------
        # EXPT 1
        if exptname1!=exptname0:
          print("..... EXPERIMENT 1 .....")
          exptpath1 = texptpath1
          print("exptpath1 = "+str(exptpath1))

          #```````````````````````````````````````````
          # Compute diurnal cycle in the obs/diags
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

            # !!! ADD THE REST OF THIS FROM EXPT 0

        #-------------------------------------------------------------------

        idd += 1
    imm += 1
  iyy += 1

print("END LOOP")
#=========================================================================================

        #---------------------------------
	# get 2D grid
del lat,lon
lon, lat = np.meshgrid(xnew, ynew)

	#---------------------------------
	# compute mean anomalies and their differences
if idatef0>0:
  sst_anom_mean0 = np.nanmean(sst_anom_arr0, axis=0)
  del sst_anom_arr0

if idatef1>0:
  sst_anom_mean1 = np.nanmean(sst_anom_arr1, axis=0)
  del sst_anom_arr1

if idate0>0:
      # mean anomaly
  omb_anom_mean0 = np.nanmean(omb0_arr, axis=0)
  oma_anom_mean0 = np.nanmean(oma0_arr, axis=0)
  obs_anom_mean0 = np.nanmean(obs0_arr, axis=0)
  bkg_anom_mean0 = np.nanmean(bkg0_arr, axis=0)
  anl_anom_mean0 = np.nanmean(anl0_arr, axis=0)
  nobs_total_sum0 = np.nansum(nobs0_arr, axis=0)
#  omb_anom_mean0 = omb0_arr / idate0
#  oma_anom_mean0 = oma0_arr / idate0
#  obs_anom_mean0 = obs0_arr / idate0
#  bkg_anom_mean0 = bkg0_arr / idate0
#  anl_anom_mean0 = anl0_arr / idate0
  print("shape omb_anom_mean0 = "+str(np.shape(omb_anom_mean0)))
  print("min/max omb_anom_mean0 = "+str(np.nanmin(omb_anom_mean0))+" "+str(np.nanmax(omb_anom_mean0)))
  print("min/max oma_anom_mean0 = "+str(np.nanmin(oma_anom_mean0))+" "+str(np.nanmax(oma_anom_mean0)))

  del omb0_arr, oma0_arr, obs0_arr, bkg0_arr, anl0_arr, nobs0_arr

if idate1>0:
      # mean anomaly
  omb_anom_mean1 = np.nanmean(omb1_arr, axis=0)
  oma_anom_mean1 = np.nanmean(oma1_arr, axis=0)
  obs_anom_mean1 = np.nanmean(obs1_arr, axis=0)
  bkg_anom_mean1 = np.nanmean(bkg1_arr, axis=0)
  anl_anom_mean1 = np.nanmean(anl1_arr, axis=0)
  nobs_total_sum1 = np.nansum(nobs1_arr, axis=0)

  del omb1_arr, oma1_arr, obs1_arr, bkg1_arr, anl1_arr, nobs1_arr

# absolute differences of SST diurnal cycle
if idate0>0 and idate1>0:
  omb_anom_diff = abs(omb_anom_mean1) - abs(omb_anom_mean0)
  oma_anom_diff = abs(oma_anom_mean1) - abs(oma_anom_mean0)
  obs_anom_diff = abs(obs_anom_mean1) - abs(obs_anom_mean0)
  bkg_anom_diff = abs(bkg_anom_mean1) - abs(bkg_anom_mean0)
  anl_anom_diff = abs(anl_anom_mean1) - abs(anl_anom_mean0)


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

varpath = var+"/"
if obstype=="all":
  geopath = "global/"
else:
  geopath = obstype+"/"

outpath = outpath+datepath+hourpath+varpath+geopath

os.system('mkdir -m 775 -p '+outpath)

    # figure name parts
figtype         = "DIURNAL_DIAGS"                                # type of plot
#figsubtype      = "OMBA"                       # type of info plotted
figexpt0        = str(exptname0)                         # experiment name
if exptname0!=exptname1:
  figexpt1        = str(exptname1)                         # experiment name
figobs          = str(obstype)                          # observation instrument

#figfiletype = str(plotstr)

figext          = "png"                                 # figure extension: examples: .png, .jpg, .pdf, .gif

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

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
# Plots

    #``````````````````````````````````````````
    # plot number of obs
pltmin = 1
pltmax = 300
plottype = ""

if idate0>0:

    figexpt = figexpt0

    z0 = nobs_total_sum0
    title = 'Sum Total Number of '+var.upper()+' Obs '+figobs+' ('+figexpt+'): '+str(hour)+'z '+daterange
    plot_map(lon, lat, z0, zmean0, title, sunlat, sunlon, plottype, pltmin, pltmax)
    outname = outpath+figtype+"."+figexpt+"."+str(hour)+"z."+var+"_"+figobs+"_anom_mean.omb."+figext
    plt.savefig(outname)
    del z0

if idate1>0:

    figexpt = figexpt1

    z0 = nobs_total_sum1
    title = 'Sum Total Number of '+var.upper()+' Obs '+figobs+' ('+figexpt+'): '+str(hour)+'z '+daterange
    plot_map(lon, lat, z0, zmean0, title, sunlat, sunlon, plottype, pltmin, pltmax)
    outname = outpath+figtype+"."+figexpt+"."+str(hour)+"z."+var+"_"+figobs+"_anom_mean.omb."+figext
    plt.savefig(outname)
    del z0

    #``````````````````````````````````````````
    # plot diurnal cycle (daily mean anomaly)
print("plot daily mean anomalies")

pltmin = -0.5
pltmax = 0.5
plottype = "Diff"       # for colorbar

if idate0>0:

    figexpt = figexpt0

    z0 = omb_anom_mean0
    title = 'Mean Anomaly of '+var.upper()+' OMB '+figobs+' ('+figexpt+'): '+str(hour)+'z '+daterange
    plot_map(lon, lat, z0, zmean0, title, sunlat, sunlon, plottype, pltmin, pltmax)
    outname = outpath+figtype+"."+figexpt+"."+str(hour)+"z."+var+"_"+figobs+"_anom_mean.omb."+figext
    plt.savefig(outname)
    del z0

    z0 = oma_anom_mean0
    title = 'Mean Anomaly of '+var.upper()+' OMA '+figobs+' ('+figexpt+'): '+str(hour)+'z '+daterange
    plot_map(lon, lat, z0, zmean0, title, sunlat, sunlon, plottype, pltmin, pltmax)
    outname = outpath+figtype+"."+figexpt+"."+str(hour)+"z."+var+"_"+figobs+"_anom_mean.oma."+figext
    plt.savefig(outname)
    del z0

    z0 = obs_anom_mean0
    title = 'Mean Anomaly of '+var.upper()+' Obs '+figobs+' ('+figexpt+'): '+str(hour)+'z '+daterange
    plot_map(lon, lat, z0, zmean0, title, sunlat, sunlon, plottype, pltmin, pltmax)
    outname = outpath+figtype+"."+figexpt+"."+str(hour)+"z."+var+"_"+figobs+"_anom_mean.obs."+figext
    plt.savefig(outname)
    del z0

    z0 = bkg_anom_mean0
    title = 'Mean Anomaly of '+var.upper()+' BKG '+figobs+' ('+figexpt+'): '+str(hour)+'z '+daterange
    plot_map(lon, lat, z0, zmean0, title, sunlat, sunlon, plottype, pltmin, pltmax)
    outname = outpath+figtype+"."+figexpt+"."+str(hour)+"z."+var+"_"+figobs+"_anom_mean.bkg."+figext
    plt.savefig(outname)
    del z0

    z0 = anl_anom_mean0
    title = 'Mean Anomaly of '+var.upper()+' ANL '+figobs+' ('+figexpt+'): '+str(hour)+'z '+daterange
    plot_map(lon, lat, z0, zmean0, title, sunlat, sunlon, plottype, pltmin, pltmax)
    outname = outpath+figtype+"."+figexpt+"."+str(hour)+"z."+var+"_"+figobs+"_anom_mean.anl."+figext
    plt.savefig(outname)
    del z0

plt.close("all")

if idate1>0:

    figexpt = figexpt1

    z0 = omb_anom_mean1
    title = 'Mean Anomaly of '+var.upper()+' OMB '+figobs+' ('+figexpt+'): '+str(hour)+'z '+daterange
    plot_map(lon, lat, z0, zmean0, title, sunlat, sunlon, plottype, pltmin, pltmax)
    outname = outpath+figtype+"."+figexpt+"."+str(hour)+"z."+var+"_"+figobs+"_anom_mean.omb."+figext
    plt.savefig(outname)
    del z0

    z0 = oma_anom_mean1
    title = 'Mean Anomaly of '+var.upper()+' OMA '+figobs+' ('+figexpt+'): '+str(hour)+'z '+daterange
    plot_map(lon, lat, z0, zmean0, title, sunlat, sunlon, plottype, pltmin, pltmax)
    outname = outpath+figtype+"."+figexpt+"."+str(hour)+"z."+var+"_"+figobs+"_anom_mean.oma."+figext
    plt.savefig(outname)
    del z0

    z0 = obs_anom_mean1
    title = 'Mean Anomaly of '+var.upper()+' Obs '+figobs+' ('+figexpt+'): '+str(hour)+'z '+daterange
    plot_map(lon, lat, z0, zmean0, title, sunlat, sunlon, plottype, pltmin, pltmax)
    outname = outpath+figtype+"."+figexpt+"."+str(hour)+"z."+var+"_"+figobs+"_anom_mean.obs."+figext
    plt.savefig(outname)
    del z0

    z0 = bkg_anom_mean1
    title = 'Mean Anomaly of '+var.upper()+' BKG '+figobs+' ('+figexpt+'): '+str(hour)+'z '+daterange
    plot_map(lon, lat, z0, zmean0, title, sunlat, sunlon, plottype, pltmin, pltmax)
    outname = outpath+figtype+"."+figexpt+"."+str(hour)+"z."+var+"_"+figobs+"_anom_mean.bkg."+figext
    plt.savefig(outname)
    del z0

    z0 = anl_anom_mean1
    title = 'Mean Anomaly of '+var.upper()+' ANL '+figobs+' ('+figexpt+'): '+str(hour)+'z '+daterange
    plot_map(lon, lat, z0, zmean0, title, sunlat, sunlon, plottype, pltmin, pltmax)
    outname = outpath+figtype+"."+figexpt+"."+str(hour)+"z."+var+"_"+figobs+"_anom_mean.anl."+figext
    plt.savefig(outname)
    del z0

plt.close("all")

if idate0>0 and idate1>0:

    z0 = omb_anom_diff
    title = 'Abs. Diff of Mean '+var.upper()+' OMB ('+figexpt1+' - '+figexpt0+'): '+str(hour)+'z '+daterange
    plot_map(lon, lat, z0, zmean0, title, sunlat, sunlon, plottype, pltmin, pltmax)
    outname = outpath+figtype+"."+figexpt1+"_MINUS_"+figexpt0+"."+str(hour)+"z."+var+"_"+figobs+"_anom_mean.omb_absdiff."+figext
    plt.savefig(outname)
    del z0

plt.close("all")

###################################################################################################################
print("END: "+str(dt.datetime.now()))
print("========== END MAIN PROGRAM ==========")
###################################################################################################################
###################################################################################################################
