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
from tools_plotting import density_scatter
from tools_plotting import plot_map_points2d_ce

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
print("obstype = "+str(obstype))

if var=="sst": units = "deg-C"

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
del lat,lon

avgopt = 0

lat = ynew
lon = xnew

    #-----------------------------------------------

mmARR = [ '01','02','03','04','05','06','07','08','09','10','11','12' ]
ddARR = [ "01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31" ]
hhARR = [ '00','06','12','18' ]

#```````````````````````````````````````
# LOOP
#```````````````````````````````````````
print("BEGIN LOOP")

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
        if texptpath0.find("role")!=-1: 
            exptpath0 = texptpath0+"/"+str(yyyy)+str(mm)+str(dd_str)+str(hour)+"/"
        else:
            exptpath0 = texptpath0
        print("exptpath0 = "+str(exptpath0))

        tnobs, tomb, tSDomb, toma, tSDoma, tobserr, tSDobserr, tobserr_orig, tSDobserr_orig, tobsval, tSDobsval, tbkg, tSDbkg, tanl, tSDanl = read_diags_ocean_gridded(exptpath0, yyyy, mm, dd_str, hour, exptname0, obstype, degree, meansumopt)

        shapeomb = np.shape(tomb)
        print("shapeomb 0 = "+str(shapeomb))

        if np.size(tomb)>1:
          print("min/max omb = "+str(np.min(tomb))+" "+str(np.max(tomb)))
          print("min/max omb no NaN = "+str(np.nanmin(tomb))+" "+str(np.nanmax(tomb)))

          nobs        = np.empty((1,shapeomb[0],shapeomb[1]),int)
          omb         = np.empty((1,shapeomb[0],shapeomb[1]),float)
          oma         = np.empty((1,shapeomb[0],shapeomb[1]),float)
          obserr      = np.empty((1,shapeomb[0],shapeomb[1]),float)
          obserr_orig = np.empty((1,shapeomb[0],shapeomb[1]),float)
          obsval      = np.empty((1,shapeomb[0],shapeomb[1]),float)
          bkg         = np.empty((1,shapeomb[0],shapeomb[1]),float)
          anl         = np.empty((1,shapeomb[0],shapeomb[1]),float)
          nobs[0,:,:]        = tnobs
          omb[0,:,:]         = tomb
          oma[0,:,:]         = toma
          obserr[0,:,:]      = tobserr
          obserr_orig[0,:,:] = tobserr_orig
          obsval[0,:,:]      = tobsval
          bkg[0,:,:]         = tbkg
          anl[0,:,:]         = tanl

          SDomb         = np.empty((1,shapeomb[0],shapeomb[1]),float)
          SDoma         = np.empty((1,shapeomb[0],shapeomb[1]),float)
          SDobserr      = np.empty((1,shapeomb[0],shapeomb[1]),float)
          SDobserr_orig = np.empty((1,shapeomb[0],shapeomb[1]),float)
          SDobsval      = np.empty((1,shapeomb[0],shapeomb[1]),float)
          SDbkg         = np.empty((1,shapeomb[0],shapeomb[1]),float)
          SDanl         = np.empty((1,shapeomb[0],shapeomb[1]),float)
          SDomb[0,:,:]         = tSDomb
          SDoma[0,:,:]         = tSDoma
          SDobserr[0,:,:]      = tSDobserr
          SDobserr_orig[0,:,:] = tSDobserr_orig
          SDobsval[0,:,:]      = tSDobsval
          SDbkg[0,:,:]         = tSDbkg
          SDanl[0,:,:]         = tSDanl

            # append to total array
          if idate0==0:
            nobs0_arr        = nobs
            omb0_arr         = omb
            oma0_arr         = oma
            obserr0_arr      = obserr
            obserr_orig0_arr = obserr_orig
            obsval0_arr      = obsval
            bkg0_arr         = bkg
            anl0_arr         = anl

            SDomb0_arr         = SDomb
            SDoma0_arr         = SDoma
            SDobserr0_arr      = SDobserr
            SDobserr_orig0_arr = SDobserr_orig
            SDobsval0_arr      = SDobsval
            SDbkg0_arr         = SDbkg
            SDanl0_arr         = SDanl

            idate0 += 1
          elif idate0>0:
            nobs0_arr    = np.append(nobs0_arr, nobs, axis=0)
            omb0_arr    = np.append(omb0_arr, omb, axis=0)
            oma0_arr    = np.append(oma0_arr, oma, axis=0)
            obserr0_arr = np.append(obserr0_arr, obserr, axis=0)
            obserr_orig0_arr = np.append(obserr_orig0_arr, obserr_orig, axis=0)
            obsval0_arr = np.append(obsval0_arr, obsval, axis=0)
            bkg0_arr    = np.append(bkg0_arr, bkg, axis=0)
            anl0_arr    = np.append(anl0_arr, anl, axis=0)

            SDomb0_arr    = np.append(SDomb0_arr, SDomb, axis=0)
            SDoma0_arr    = np.append(SDoma0_arr, SDoma, axis=0)
            SDobserr0_arr = np.append(SDobserr0_arr, SDobserr, axis=0)
            SDobserr_orig0_arr = np.append(SDobserr_orig0_arr, SDobserr_orig, axis=0)
            SDobsval0_arr = np.append(SDobsval0_arr, SDobsval, axis=0)
            SDbkg0_arr    = np.append(SDbkg0_arr, SDbkg, axis=0)
            SDanl0_arr    = np.append(SDanl0_arr, SDanl, axis=0)

            idate0 += 1
          del nobs, omb, SDomb, oma, SDoma, obserr, SDobserr, obserr_orig, SDobserr_orig, obsval, SDobsval, bkg, SDbkg, anl, SDanl
        del shapeomb
        del tnobs, tomb, tSDomb, toma, tSDoma, tobserr, tSDobserr, tobserr_orig, tSDobserr_orig, tobsval, tSDobsval, tbkg, tSDbkg, tanl, tSDanl

        #-------------------------------------------------------------------
        # EXPT 1
        if exptname1!=exptname0:
          if texptpath1.find("role")!=-1: 
              exptpath1 = texptpath1+"/"+str(yyyy)+str(mm)+str(dd_str)+str(hour)+"/"
          else:
              exptpath1 = texptpath1
          print("exptpath1 = "+str(exptpath1))

          tnobs, tomb, tSDomb, toma, tSDoma, tobserr, tSDobserr, tobserr_orig, tSDobserr_orig, tobsval, tSDobsval, tbkg, tSDbkg, tanl, tSDanl = read_diags_ocean_gridded(exptpath1, yyyy, mm, dd_str, hour, exptname1, obstype, degree, meansumopt)
 
          shapeomb = np.shape(tomb)
          print("shapeomb 1 = "+str(shapeomb))

          if np.size(tomb)>1:
            nobs        = np.empty((1,shapeomb[0],shapeomb[1]),int)
            omb         = np.empty((1,shapeomb[0],shapeomb[1]),float)
            oma         = np.empty((1,shapeomb[0],shapeomb[1]),float)
            obserr      = np.empty((1,shapeomb[0],shapeomb[1]),float)
            obserr_orig = np.empty((1,shapeomb[0],shapeomb[1]),float)
            obsval      = np.empty((1,shapeomb[0],shapeomb[1]),float)
            bkg         = np.empty((1,shapeomb[0],shapeomb[1]),float)
            anl         = np.empty((1,shapeomb[0],shapeomb[1]),float)
            nobs[0,:,:]        = tnobs
            omb[0,:,:]         = tomb
            oma[0,:,:]         = toma
            obserr[0,:,:]      = tobserr
            obserr_orig[0,:,:] = tobserr_orig
            obsval[0,:,:]      = tobsval
            bkg[0,:,:]         = tbkg
            anl[0,:,:]         = tanl

            SDomb         = np.empty((1,shapeomb[0],shapeomb[1]),float)
            SDoma         = np.empty((1,shapeomb[0],shapeomb[1]),float)
            SDobserr      = np.empty((1,shapeomb[0],shapeomb[1]),float)
            SDobserr_orig = np.empty((1,shapeomb[0],shapeomb[1]),float)
            SDobsval      = np.empty((1,shapeomb[0],shapeomb[1]),float)
            SDbkg         = np.empty((1,shapeomb[0],shapeomb[1]),float)
            SDanl         = np.empty((1,shapeomb[0],shapeomb[1]),float)
            SDomb[0,:,:]         = tSDomb
            SDoma[0,:,:]         = tSDoma
            SDobserr[0,:,:]      = tSDobserr
            SDobserr_orig[0,:,:] = tSDobserr_orig
            SDobsval[0,:,:]      = tSDobsval
            SDbkg[0,:,:]         = tSDbkg
            SDanl[0,:,:]         = tSDanl

            # append to total array
            if idate1==0:
              nobs1_arr        = nobs
              omb1_arr         = omb
              oma1_arr         = oma
              obserr1_arr      = obserr
              obserr_orig1_arr = obserr_orig
              obsval1_arr      = obsval
              bkg1_arr         = bkg
              anl1_arr         = anl

              SDomb1_arr         = SDomb
              SDoma1_arr         = SDoma
              SDobserr1_arr      = SDobserr
              SDobserr_orig1_arr = SDobserr_orig
              SDobsval1_arr      = SDobsval
              SDbkg1_arr         = SDbkg
              SDanl1_arr         = SDanl

              idate1 += 1
            elif idate1>0:
              nobs1_arr   = np.append(nobs1_arr, nobs, axis=0)
              omb1_arr    = np.append(omb1_arr, omb, axis=0)
              oma1_arr    = np.append(oma1_arr, oma, axis=0)
              obserr1_arr = np.append(obserr1_arr, obserr, axis=0)
              obserr_orig1_arr = np.append(obserr_orig1_arr, obserr_orig, axis=0)
              obsval1_arr = np.append(obsval1_arr, obsval, axis=0)
              bkg1_arr    = np.append(bkg1_arr, bkg, axis=0)
              anl1_arr    = np.append(anl1_arr, anl, axis=0)

              SDomb1_arr    = np.append(SDomb1_arr, SDomb, axis=0)
              SDoma1_arr    = np.append(SDoma1_arr, SDoma, axis=0)
              SDobserr1_arr = np.append(SDobserr1_arr, SDobserr, axis=0)
              SDobserr_orig1_arr = np.append(SDobserr_orig1_arr, SDobserr_orig, axis=0)
              SDobsval1_arr = np.append(SDobsval1_arr, SDobsval, axis=0)
              SDbkg1_arr    = np.append(SDbkg1_arr, SDbkg, axis=0)
              SDanl1_arr    = np.append(SDanl1_arr, SDanl, axis=0)

              idate1 += 1
            del nobs, omb, SDomb, oma, SDoma, obserr, SDobserr, obserr_orig, SDobserr_orig, obsval, SDobsval, bkg, SDbkg, anl, SDanl
          del shapeomb
          del tnobs, tomb, tSDomb, toma, tSDoma, tobserr, tSDobserr, tobserr_orig, tSDobserr_orig, tobsval, tSDobsval, tbkg, tSDbkg, tanl, tSDanl

        #-------------------------------------------------------------------

        idd += 1
    imm += 1
  iyy += 1

print("END LOOP")

#sys.exit()

#===============================================================
#===============================================================
# PLOT prep
#===============================================================

    # edit outpath
if (str(yyyyS)+str(mmS)+str(ddS))==(str(yyyyE)+str(mmE)+str(ddE)):
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
figtype         = "MAP"                                 # type of plot
figsubtype      = "DIAGS"                                # plot subtype
figobs          = str(obstype)                          # observation instrument

#===============================================================
# output stats in text file

statstr = ['OMB', 'OMA', 'ObsError', 'ObsError_Orig']

if idate0>0:
  for istat in statstr:
    if istat=='OMB':
      stat = omb0_arr
    elif istat=='OMA':
      stat = oma0_arr
    elif istat=='ObsError':
      stat = obserr0_arr
    elif istat=='ObsError_Orig':
      stat = obserr_orig0_arr

    title = 'GRIDDED DIAG STATS: '+var.upper()+' '+str(istat)+' '+str(obstype)+' ('+str(exptname0)+') for '+str(hour)+' UTC: '+daterange
    outname = outpath+"STATS_DIAGS_GRIDDED."+exptname0+"."+str(hour)+"z."+var+"."+str(obstype)+"."+str(istat.lower())

        # compute and print stats of differences
    avg  = np.nanmean(stat)                  # mean
    sd   = np.nanstd(stat)                   # standard deviation
    #rmsd = np.sqrt(np.mean(diff**2))      # RMSD

                # ADD TEXT INSIDE PANEL
                # add title for difference stats
    text = title
                # count
    textc = "... Count = "+str(np.size(np.where(~np.isnan(stat))))
                # add mean diff
    textm = "... Mean = "+str(np.round_(avg,decimals=3))
                # add SD of diff
    texts = "... StdDev = "+str(np.round_(sd,decimals=3))
                # add RMSD
    #texts = "RMSD = "+str(np.round_(rmsd,decimals=2))

                # OUTPUT STATS TO TEXT FILE
    file1 = open(outname+".txt","a+")            # open file for reading and writing. create file if it does not exist. data is appended to end of file.
    file1.write(text+"\n")
    file1.write(textc+"\n")
    file1.write(textm+"\n")
    file1.write(texts+"\n")
    #textrmsd = "RMSD = "+str(np.round_(rmsd,decimals=2))
    file1.close()
    del stat,avg,sd,text,textm,texts
    del title, outname

if idate1>0:
  for istat in statstr:
    if istat=='OMB':
      stat = omb1_arr
    elif istat=='OMA':
      stat = oma1_arr
    elif istat=='ObsError':
      stat = obserr1_arr
    elif istat=='ObsError_Orig':
      stat = obserr_orig1_arr

    title = 'GRIDDED DIAG STATS: '+var.upper()+' '+str(istat)+' '+str(obstype)+' ('+str(exptname1)+') for '+str(hour)+' UTC: '+daterange
    outname = outpath+"STATS_DIAGS_GRIDDED."+exptname1+"."+str(hour)+"z."+var+"."+str(obstype)+"."+str(istat.lower())

        # compute and print stats of differences
    avg  = np.nanmean(stat)                  # mean
    sd   = np.nanstd(stat)                   # standard deviation
    #rmsd = np.sqrt(np.mean(diff**2))      # RMSD

                # ADD TEXT INSIDE PANEL
                # add title for difference stats
    text = title
                # count
    textc = "... Count = "+str(np.size(np.where(~np.isnan(stat))))
                # add mean diff
    textm = "... Mean = "+str(np.round_(avg,decimals=3))
                # add SD of diff
    texts = "... StdDev = "+str(np.round_(sd,decimals=3))
                # add RMSD
    #texts = "RMSD = "+str(np.round_(rmsd,decimals=2))

                # OUTPUT STATS TO TEXT FILE
    file1 = open(outname+".txt","a+")            # open file for reading and writing. create file if it does not exist. data is appended to end of file.
    file1.write(text+"\n")
    file1.write(textc+"\n")
    file1.write(textm+"\n")
    file1.write(texts+"\n")
    #textrmsd = "RMSD = "+str(np.round_(rmsd,decimals=2))
    file1.close()
    del stat,avg,sd,text,textm,texts
    del title, outname

#===============================================================

figext          = "png"                                 # figure extension: examples: .png, .jpg, .pdf, .gif

        #---------------------------------
        # Means of arrays

if idate0>0:
  omb2d0         = np.nanmean(omb0_arr, axis=0)
  oma2d0         = np.nanmean(oma0_arr, axis=0)
  obserr2d0      = np.nanmean(obserr0_arr, axis=0)
  obserr_orig2d0 = np.nanmean(obserr_orig0_arr, axis=0)
  obsval2d0      = np.nanmean(obsval0_arr, axis=0)
  bkg2d0         = np.nanmean(bkg0_arr, axis=0)
  anl2d0         = np.nanmean(anl0_arr, axis=0)

  SDomb2d0         = np.nanmean(SDomb0_arr, axis=0)
  SDoma2d0         = np.nanmean(SDoma0_arr, axis=0)
  SDobserr2d0      = np.nanmean(SDobserr0_arr, axis=0)
  SDobserr_orig2d0 = np.nanmean(SDobserr_orig0_arr, axis=0)
  SDobsval2d0      = np.nanmean(SDobsval0_arr, axis=0)
  SDbkg2d0         = np.nanmean(SDbkg0_arr, axis=0)
  SDanl2d0         = np.nanmean(SDanl0_arr, axis=0)

if idate1>0:
  omb2d1         = np.nanmean(omb1_arr, axis=0)
  oma2d1         = np.nanmean(oma1_arr, axis=0)
  obserr2d1      = np.nanmean(obserr1_arr, axis=0)
  obserr_orig2d1 = np.nanmean(obserr_orig1_arr, axis=0)
  obsval2d1      = np.nanmean(obsval1_arr, axis=0)
  bkg2d1         = np.nanmean(bkg1_arr, axis=0)
  anl2d1         = np.nanmean(anl1_arr, axis=0)

  SDomb2d1         = np.nanmean(SDomb1_arr, axis=0)
  SDoma2d1         = np.nanmean(SDoma1_arr, axis=0)
  SDobserr2d1      = np.nanmean(SDobserr1_arr, axis=0)
  SDobserr_orig2d1 = np.nanmean(SDobserr_orig1_arr, axis=0)
  SDobsval2d1      = np.nanmean(SDobsval1_arr, axis=0)
  SDbkg2d1         = np.nanmean(SDbkg1_arr, axis=0)
  SDanl2d1         = np.nanmean(SDanl1_arr, axis=0)

        #---------------------------------
        # Difference: Expt1 - Expt0

if idate0>0 and idate1>0:
  omb_diff = omb2d1 - omb2d0
  oma_diff = oma2d1 - oma2d0
  omb_absdiff = abs(omb2d1) - abs(omb2d0)
  oma_absdiff = abs(oma2d1) - abs(oma2d0)

            # diff between SST background and observations in observation space
  diff0_bkg_obs = bkg2d0 - obsval2d0
  diff1_bkg_obs = bkg2d1 - obsval2d1
            # diff between SST obs from each experiment
  diff_obs_obs = obsval2d1 - obsval2d0

        #---------------------------------
	# get 2D grid
del lat,lon
lon, lat = np.meshgrid(xnew, ynew)

#*****************************************************************************************
#*****************************************************************************************
# PLOTS
#*****************************************************************************************
#*****************************************************************************************

zmean = "none"

	# sun location (for July 1-15, 2021)
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

#------------------------------------
# Observations, Backgrounds, and Analysis Fields

pltmin = 0
pltmax = 30
plottype = ""

        #------------------------------------
        # EXPT 0
if idate0>0:
  figexpt = str(exptname0)

  z = obsval2d0
  title = var.upper()+' obs '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".obs."+figext
  plt.savefig(outname)
  del z

  z = bkg2d0
  title = var.upper()+' bkg '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".bkg."+figext
  plt.savefig(outname)
  del z

  z = anl2d0
  title = var.upper()+' anl '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".anl."+figext
  plt.savefig(outname)
  del z

        #------------------------------------
        # EXPT 1
if idate1>0:
  figexpt = str(exptname1)

  z = obsval2d1
  title = var.upper()+' obs '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".obs."+figext
  plt.savefig(outname)
  del z

  z = bkg2d1
  title = var.upper()+' bkg '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".bkg."+figext
  plt.savefig(outname)
  del z

  z = anl2d1
  title = var.upper()+' anl '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".anl."+figext
  plt.savefig(outname)
  del z

#------------------------------------
# Innovations and Observation Errors

        #------------------------------------
        # EXPT 0
if idate0>0:
  figexpt = exptname0

  z = omb2d0
  pltmin = -0.5
  pltmax = 0.5
  plottype = "Diff"
  title = var.upper()+' OMB '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".omb."+figext
  plt.savefig(outname)
  del z

  z = oma2d0
  pltmin = -0.5
  pltmax = 0.5
  plottype = "Diff"
  title = var.upper()+' OMA '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".oma."+figext
  plt.savefig(outname)
  del z

  z = obserr2d0
  pltmin = 0
  pltmax = 0.5
  plottype = ""
  title = var.upper()+' Obs Error '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".obs_error."+figext
  plt.savefig(outname)
  del z

        # OMB / obs error ratio
  z = omb2d0/obserr2d0
  pltmin = -0.5
  pltmax = 0.5
  plottype = "Diff"
  title = var.upper()+' OMB/ObsError '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".ratio_omb2obserror."+figext
  plt.savefig(outname)
  del z

  z = oma2d0/obserr2d0
  pltmin = -0.5
  pltmax = 0.5
  plottype = "Diff"
  title = var.upper()+' OMA/ObsError '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".ratio_oma2obserror."+figext
  plt.savefig(outname)
  del z

        #------------------------------------
        # EXPT 1
if idate1>0:
  figexpt = exptname1

  z = omb2d1
  pltmin = -0.5
  pltmax = 0.5
  plottype = "Diff"
  title = var.upper()+' OMB '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".omb."+figext
  plt.savefig(outname)
  del z

  z = oma2d1
  pltmin = -0.5
  pltmax = 0.5
  plottype = "Diff"
  title = var.upper()+' OMA '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".oma."+figext
  plt.savefig(outname)
  del z

  z = obserr2d1
  pltmin = 0
  pltmax = 0.5
  plottype = ""
  title = var.upper()+' Obs Error '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".obs_error."+figext
  plt.savefig(outname)
  del z

#------------------------------------
# Standard Deviations

pltmin = 0
pltmax = 0.5
plottype = "SD"

        # EXPT 0
if idate0>0:
  figexpt = exptname0

    # Innovations
  z = SDomb2d0
  title = 'SD of '+var.upper()+' OMB '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".omb_stddev."+figext
  plt.savefig(outname)
  del z

  z = SDoma2d0
  title = 'SD of '+var.upper()+' OMA '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".oma_stddev."+figext
  plt.savefig(outname)
  del z

    # Obs Errors
  z = SDobserr2d0
  title = 'SD of '+var.upper()+' Obs Error '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".obs_error_stddev."+figext
  plt.savefig(outname)
  del z

    # Obs, Backgrounds, Analysis fields
  z = SDobsval2d0
  title = 'SD of '+var.upper()+' obs '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".obs_stddev."+figext
  plt.savefig(outname)
  del z

  z = SDbkg2d0
  title = 'SD of '+var.upper()+' bkg '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".bkg_stddev."+figext
  plt.savefig(outname)
  del z

  z = SDanl2d0
  title = 'SD of '+var.upper()+' anl '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".anl_stddev."+figext
  plt.savefig(outname)
  del z

        # EXPT 1
if idate1>0:
  figexpt = exptname1

    # Innovations
  z = SDomb2d1
  title = 'SD of '+var.upper()+' OMB '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".omb_stddev."+figext
  plt.savefig(outname)
  del z

  z = SDoma2d1
  title = 'SD of '+var.upper()+' OMA '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".oma_stddev."+figext
  plt.savefig(outname)
  del z

    # Obs Errors
  z = SDobserr2d1
  title = 'SD of '+var.upper()+' Obs Error '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".obs_error_stddev."+figext
  plt.savefig(outname)
  del z

    # Obs, Backgrounds, Analysis fields
  z = SDobsval2d1
  title = 'SD of '+var.upper()+' obs '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".obs_stddev."+figext
  plt.savefig(outname)
  del z

  z = SDbkg2d1
  title = 'SD of '+var.upper()+' bkg '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".bkg_stddev."+figext
  plt.savefig(outname)
  del z

  z = SDanl2d1
  title = 'SD of '+var.upper()+' anl '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".anl_stddev."+figext
  plt.savefig(outname)
  del z

#------------------------------------
# Differences between experiments

if idate0>0 and idate1>0:

  figexpt0 = exptname0
  figexpt1 = exptname1

  plottype = "Diff"

  pltmin = -0.5
  pltmax = 0.5

        # diff in innovations between experiments
  z = omb_diff
  title = 'Diff in '+var.upper()+' OMB '+str(obstype)+' ('+str(figexpt1)+'-'+str(figexpt0)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt1+"_MINUS_"+figexpt0+"."+str(hour)+"z."+var+"."+figobs+".omb_diff."+figext
  plt.savefig(outname)
  del z

  z = omb_absdiff
  title = 'Abs Diff in '+var.upper()+' OMB '+str(obstype)+' ('+str(figexpt1)+'-'+str(figexpt0)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt1+"_MINUS_"+figexpt0+"."+str(hour)+"z."+var+"."+figobs+".omb_absdiff."+figext
  plt.savefig(outname)
  del z

  z = oma_diff
  title = 'Diff in '+var.upper()+' OMA '+str(obstype)+' ('+str(figexpt1)+'-'+str(figexpt0)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt1+"_MINUS_"+figexpt0+"."+str(hour)+"z."+var+"."+figobs+".oma_diff."+figext
  plt.savefig(outname)
  del z

  z = oma_absdiff
  title = 'Abs Diff in '+var.upper()+' OMA '+str(obstype)+' ('+str(figexpt1)+'-'+str(figexpt0)+') for '+str(hour)+' UTC: '+daterange
  plot_map(lon, lat, z, zmean, title, sunlat, sunlon, plottype, pltmin, pltmax)
  outname = outpath+figtype+"."+figsubtype+"."+figexpt1+"_MINUS_"+figexpt0+"."+str(hour)+"z."+var+"."+figobs+".oma_absdiff."+figext
  plt.savefig(outname)
  del z

###################################################################################################################
print("END: "+str(dt.datetime.now()))
print("========== END MAIN PROGRAM ==========")
###################################################################################################################
###################################################################################################################
