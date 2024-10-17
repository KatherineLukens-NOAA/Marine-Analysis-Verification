###################################################################################################################
###################################################################################################################
#
# Plot diag data as points on maps
#
###################################################################################################################
###################################################################################################################
#
# Output: Figures
#
###################################################################################################################
# Import python modules
###################################################################################################################
print("***** Import Python Modules *****")

import os
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import sys

from read_data import read_diags_ocean

from tools_plotting import plot_map_points2d_ce
from tools_plotting import density_scatter

###################################################################################################################
print("========== BEGIN MAIN PROGRAM ==========")

print("START: "+str(dt.datetime.now()))

#=============================================
# Read raw user input from command line

yyyymmddhhS = sys.argv[1]         # Start date: year month day
yyyymmddhhE = sys.argv[2]         # End date: year month day
outpath   = sys.argv[3]
exptname0 = sys.argv[4]
texptpath0 = sys.argv[5]
#exptname1 = sys.argv[6]
#texptpath1 = sys.argv[7]
var       = sys.argv[6]
hour      = str(sys.argv[7])            # analysis hour (00, 06, 12, or 18 UTC)
obstype   = str(sys.argv[8])

print("obstype = "+str(obstype))
print("exptname0 = "+str(exptname0))
#print("exptname1 = "+str(exptname1))

if var=="sst": units = "deg-C"

#=============================================
# Set global parameters

fill = -999.0

    # for plotting only
gridsize = 1	# units: degrees
flats = list(np.arange(-90,91,gridsize))
flons = list(np.arange(-180,181,gridsize))
nlon = np.size(flons)

xnew = flons            # new lon to conform to
ynew = flats            # new lat to conform to

lat = ynew
lon = xnew

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

avgopt = 0

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

        omb, oma, obserr, obserr_orig, obsval, qc, time, lats, lons = read_diags_ocean(exptpath0, yyyy, mm, dd_str, hour, obstype, nlon, opt)

        if np.size(omb)>1:

          bkg = obsval - omb
          anl = obsval - oma

            # append to total array
          if idate0==0:
            omb0_arr    = omb
            oma0_arr    = oma
            obserr0_arr = obserr
            obserr_orig0_arr = obserr_orig
            obsval0_arr = obsval
            qc0_arr     = qc
            time0_arr   = time
            lats0_arr   = lats
            lons0_arr   = lons
            bkg0_arr    = bkg
            anl0_arr    = anl
            idate0 += 1
          elif idate0>0:
            omb0_arr    = np.append(omb0_arr, omb, axis=0)
            oma0_arr    = np.append(oma0_arr, oma, axis=0)
            obserr0_arr = np.append(obserr0_arr, obserr, axis=0)
            obserr_orig0_arr = np.append(obserr_orig0_arr, obserr_orig, axis=0)
            obsval0_arr = np.append(obsval0_arr, obsval, axis=0)
            qc0_arr     = np.append(qc0_arr, qc, axis=0)
            time0_arr   = np.append(time0_arr, time, axis=0)
            lats0_arr   = np.append(lats0_arr, lats, axis=0)
            lons0_arr   = np.append(lons0_arr, lons, axis=0)
            bkg0_arr    = np.append(bkg0_arr, bkg, axis=0)
            anl0_arr    = np.append(anl0_arr, anl, axis=0)
            idate0 += 1
          del bkg, anl
        del omb, oma, obserr, obsval, qc, time, lats, lons#, bkg, anl

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

    title = 'DIAG STATS: '+var.upper()+' '+str(istat)+' '+str(obstype)+' ('+str(exptname0)+') for '+str(hour)+' UTC: '+daterange
    outname = outpath+"STATS_DIAGS."+exptname0+"."+str(hour)+"z."+var+"."+str(obstype)+"."+str(istat.lower())

        # compute and print stats of differences
    avg  = np.nanmean(stat)                  # mean
    sd   = np.nanstd(stat)                   # standard deviation
    if istat=='OMB' or istat=='OMA':
      rmsd = np.sqrt(np.mean(stat**2))      # RMSD

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
    if istat=='OMB' or istat=='OMA':
      textrmsd = "... RMSD = "+str(np.round_(rmsd,decimals=3))

                # OUTPUT STATS TO TEXT FILE
    os.system("rm "+outname+".txt")
    file1 = open(outname+".txt","a+")            # open file for reading and writing. create file if it does not exist. data is appended to end of file.
    file1.write(text+"\n")
    file1.write(textc+"\n")
    file1.write(textm+"\n")
    file1.write(texts+"\n")
    if istat=='OMB' or istat=='OMA':
      file1.write(textrmsd+"\n")
    file1.close()
    del stat,avg,sd,text,textm,texts
    del title, outname
    if istat=='OMB' or istat=='OMA': del textrmsd

#===============================================================

figext          = "png"                                 # figure extension: examples: .png, .jpg, .pdf, .gif

#************************************************
# Plot variables as points on map
#************************************************

alphaval = 0.25         # transparency factor (between 0 and 1, with 1=opaque)

opt=0
marksize = 2 #15

if var.find("sst")!=-1 or var.find("temp")!=-1: units = "deg-C"

if idate0>0:
  figexpt = str(exptname0)

  x = lons0_arr
  y = lats0_arr

  pltmin = 0
  pltmax = 30
  plottype = ""

  # plot obs value, bkg, anl on maps
  z = obsval0_arr
  title = var.upper()+' obs '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map_points2d_ce(var.upper(), x, y, z, lon, lat, title, marksize, alphaval, daterange, units, opt, plottype, pltmin, pltmax)
  outname = outpath+figtype+"_POINTS."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".obs."+figext
  plt.savefig(outname)
  del outname
  del z

  z = bkg0_arr
  title = var.upper()+' bkg '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map_points2d_ce(var.upper(), x, y, z, lon, lat, title, marksize, alphaval, daterange, units, opt, plottype, pltmin, pltmax)
  outname = outpath+figtype+"_POINTS."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".bkg."+figext
  plt.savefig(outname)
  del outname
  del z

  z = anl0_arr
  title = var.upper()+' anl '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map_points2d_ce(var.upper(), x, y, z, lon, lat, title, marksize, alphaval, daterange, units, opt, plottype, pltmin, pltmax)
  outname = outpath+figtype+"_POINTS."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".anl."+figext
  plt.savefig(outname)
  del outname
  del z

  # plot density scatterplots of bkg/anl vs obs
  regionstr = "Global"

  z0 = obsval0_arr
  z1 = bkg0_arr
  plotvar = 'bkg_vs_obs'
  title = var.upper()+' '+plotvar+' '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  toutname = outpath+"DENS_SCAT_POINTS."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+"."+plotvar
  density_scatter(z0, z1, 'OBS', 'BKG', toutname, plotvar, title, var, units, regionstr, pltmin, pltmax+10)
  outname = toutname+"."+figext
  plt.savefig(outname)
  del z0, z1

  z0 = obsval0_arr
  z1 = anl0_arr
  plotvar = 'anl_vs_obs'
  title = var.upper()+' '+plotvar+' '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  toutname = outpath+"DENS_SCAT_POINTS."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+"."+plotvar
  density_scatter(z0, z1, 'OBS', 'ANL', toutname, plotvar, title, var, units, regionstr, pltmin, pltmax+10)
  outname = toutname+"."+figext
  plt.savefig(outname)
  del z0, z1

  # plot OMB and OMA on maps
  pltmin = -0.5
  pltmax = 0.5
  plottype = "Diff"

  z = omb0_arr
  title = var.upper()+' OMB '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map_points2d_ce(var.upper(), x, y, z, lon, lat, title, marksize, alphaval, daterange, units, opt, plottype, pltmin, pltmax)
  outname = outpath+figtype+"_POINTS."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".omb."+figext
  plt.savefig(outname)
  del outname
  del z

  z = oma0_arr
  title = var.upper()+' OMA '+str(obstype)+' ('+str(figexpt)+') for '+str(hour)+' UTC: '+daterange
  plot_map_points2d_ce(var.upper(), x, y, z, lon, lat, title, marksize, alphaval, daterange, units, opt, plottype, pltmin, pltmax)
  outname = outpath+figtype+"_POINTS."+figsubtype+"."+figexpt+"."+str(hour)+"z."+var+"."+figobs+".oma."+figext
  plt.savefig(outname)
  del outname
  del z

  del x,y

###################################################################################################################
print("END: "+str(dt.datetime.now()))
print("========== END MAIN PROGRAM ==========")
###################################################################################################################
###################################################################################################################
