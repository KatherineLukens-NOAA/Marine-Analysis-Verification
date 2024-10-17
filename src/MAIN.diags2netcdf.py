###################################################################################################################
###################################################################################################################
#
# Output diag data (originally in obs space) onto 2D lat/lon grid
#
###################################################################################################################
###################################################################################################################
#
# Output: NetCDF-4 
#
###################################################################################################################
# Import python modules
###################################################################################################################
print("***** Import Python Modules *****")

import os
import numpy as np

import datetime as dt
import sys
from netCDF4 import Dataset

from read_data import read_diags_ocean

from tools_analysis import get_ocngrid
from tools_analysis import vars_to_2d

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
var         = sys.argv[6]
hour        = str(sys.argv[7])            # analysis hour (00, 06, 12, or 18 UTC)
obstype     = str(sys.argv[8])
degree      = sys.argv[9]
meansumopt  = sys.argv[10]

print("obstype = "+str(obstype))
print("exptname0 = "+str(exptname0))
print("meansumopt = "+str(meansumopt)+" type = "+str(type(meansumopt)))

if var=="sst": units = "deg-C"

#=============================================
# Set global parameters

fill = -999.0

gridsize = float(degree)	# units: degrees
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

#lat, lon = get_ocngrid()
#shape_lat = np.shape(lat)
#nlat = shape_lat[0]
#nlon = shape_lat[1]
#print("shape lat = "+str(shape_lat)+" nlat = "+str(nlat)+" nlon = "+str(nlon))
#del lat,lon

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
	
        cdate = str(yyyy)+str(mm)+str(dd_str)+str(hour)

        opt = -1        # all points (do not stratify into day/night)

        #-------------------------------------------------------------------
	# EXPT 0
        if texptpath0.find("role")!=-1 or texptpath0.find("marine_candidate")!=-1: 
            exptpath0 = texptpath0+"/"+str(yyyy)+str(mm)+str(dd_str)+str(hour)+"/"
        else:
            exptpath0 = texptpath0
        print("exptpath0 = "+str(exptpath0))

        nlon = 0    # only needed for day/night stratification
        omb, oma, obserr, obserr_orig, obsval, qc, time, lats, lons = read_diags_ocean(exptpath0, yyyy, mm, dd_str, hour, obstype, nlon, opt)

        if np.size(omb)>1:

          bkg = obsval - omb
          anl = obsval - oma
	  
	  #------------------------------------------------
	  # Conform 1d arrays to 2d (lat, lon)
	  #------------------------------------------------

		# output mean (first var) and SD (second var) per grid cell
          xold             = lons         # old lon
          yold             = lats         # old lats
          omb_expt         = omb
          oma_expt         = oma
          obserr_expt      = obserr
          obserr_orig_expt = obserr_orig
          obsval_expt      = obsval
          bkg_expt         = bkg
          anl_expt         = anl
          del bkg, anl

          num2d, omb2d, SDomb2d, oma2d, SDoma2d, obserr2d, SDobserr2d, obserr_orig2d, SDobserr_orig2d, obsval2d, SDobsval2d, bkg2d, SDbkg2d, anl2d, SDanl2d = vars_to_2d(xstr, ystr, xnew, ynew, xold, yold, meansumopt, omb_expt, oma_expt, obserr_expt, obserr_orig_expt, obsval_expt, bkg_expt, anl_expt)
          del xold, yold, omb_expt, oma_expt, obserr_expt, obserr_orig_expt, obsval_expt, bkg_expt, anl_expt

          print("shape omb2d = "+str(np.shape(omb2d)))
          print("size flats = "+str(np.size(flats))+" flons = "+str(np.size(flons)))
          print("max |omb2d| = "+str(np.max(abs(omb2d)))+" sans nan = "+str(np.nanmax(abs(omb2d))))
          print("min/max omb2d = "+str(np.min(omb2d))+" "+str(np.max(omb2d)))
          print("min/max omb2d sans nan = "+str(np.nanmin(omb2d))+" "+str(np.nanmax(omb2d)))

                # set missing value to zero
#          if np.isnan(np.max(abs(omb2d))) or np.nanmax(abs(omb2d))>999:
#            print("set missing value to zero")
#
#            omb2d           = np.where(np.isnan(omb2d)==True, 0, omb2d)
#            print("min/max omb2d after missing=0: "+str(np.min(omb2d))+" "+str(np.max(omb2d)))
#
#            num2d           = np.where(np.isnan(omb2d)==True or abs(omb2d)>999, 0, num2d)
#            omb2d           = np.where(np.isnan(omb2d)==True or abs(omb2d)>999, 0, omb2d)
#            oma2d           = np.where(np.isnan(omb2d)==True or abs(omb2d)>999, 0, oma2d)
#            obserr2d        = np.where(np.isnan(omb2d)==True or abs(omb2d)>999, 0, obserr2d)
#            obserr_orig2d   = np.where(np.isnan(omb2d)==True or abs(omb2d)>999, 0, obserr_orig2d)
#            obsval2d        = np.where(np.isnan(omb2d)==True or abs(omb2d)>999, 0, obsval2d)
#            bkg2d           = np.where(np.isnan(omb2d)==True or abs(omb2d)>999, 0, bkg2d)
#            anl2d           = np.where(np.isnan(omb2d)==True or abs(omb2d)>999, 0, anl2d)
#
#            SDomb2d           = np.where(np.isnan(omb2d)==True or abs(omb2d)>999, 0, SDomb2d)
#            SDoma2d           = np.where(np.isnan(omb2d)==True or abs(omb2d)>999, 0, SDoma2d)
#            SDobserr2d        = np.where(np.isnan(omb2d)==True or abs(omb2d)>999, 0, SDobserr2d)
#            SDobserr_orig2d   = np.where(np.isnan(omb2d)==True or abs(omb2d)>999, 0, SDobserr_orig2d)
#            SDobsval2d        = np.where(np.isnan(omb2d)==True or abs(omb2d)>999, 0, SDobsval2d)
#            SDbkg2d           = np.where(np.isnan(omb2d)==True or abs(omb2d)>999, 0, SDbkg2d)
#            SDanl2d           = np.where(np.isnan(omb2d)==True or abs(omb2d)>999, 0, SDanl2d)
#
#            print("min/max omb2d after missing=0: "+str(np.min(omb2d))+" "+str(np.max(omb2d)))

          #------------------------------------------------
          # OUTPUT NETCDF
          print("output netcdf")

          if int(meansumopt)==0:
            meansumstr = "means"
          elif int(meansumopt)==1:
            meansumstr = "sums"

          nc_out_filename = "DIAGS."+str(exptname0)+"."+str(yyyy)+str(mm)+str(dd_str)+str(hour)+"."+str(obstype)+"_obs.global_"+str(degree)+"deg."+str(meansumstr)+".nc4"

            # remove old file
          os.system("rm -rf "+outpath+nc_out_filename)

          nc_out = Dataset( #.................................................... Dataset object for output
                      outpath+nc_out_filename  , # Dataset input: Output file name
                      "w"              , # Dataset input: Make file write-able
                      format="NETCDF4" , # Dataset input: Set output format to netCDF4
                    )
          print("... nc_out = "+str(outpath)+str(nc_out_filename))

            # create dimensions for data variables
          print("... create lat/lon dims")
          nc_out.createDimension("lat", np.size(flats))   
          nc_out.createDimension("lon", np.size(flons)) 

            # create dimension variables
          print("... create lat/lon vars")
          y = nc_out.createVariable("lat", "f8", ("lat",))
          x = nc_out.createVariable("lon", "f8", ("lon",))

            # fill dimension variables
          print("... fill lat/lon vars")
          nc_out["lat"][:] = np.array(flats)
          nc_out["lon"][:] = np.array(flons)

            # create data variables
          print("... create vars")
          out_omb           = nc_out.createVariable("omb", "f8", ("lat","lon"))
          out_oma           = nc_out.createVariable("oma", "f8", ("lat","lon"))
          out_obserr        = nc_out.createVariable("obserr_postQC", "f8", ("lat","lon"))
          out_obserr_orig   = nc_out.createVariable("obserr_beforeQC", "f8", ("lat","lon"))
          out_obsval        = nc_out.createVariable("obsvalue", "f8", ("lat","lon"))
          out_bkg           = nc_out.createVariable("bkg", "f8", ("lat","lon"))
          out_anl           = nc_out.createVariable("anl", "f8", ("lat","lon"))

          out_SDomb         = nc_out.createVariable("SDomb", "f8", ("lat","lon"))
          out_SDoma         = nc_out.createVariable("SDoma", "f8", ("lat","lon"))
          out_SDobserr      = nc_out.createVariable("SDobserr_postQC", "f8", ("lat","lon"))
          out_SDobserr_orig = nc_out.createVariable("SDobserr_beforeQC", "f8", ("lat","lon"))
          out_SDobsval      = nc_out.createVariable("SDobsvalue", "f8", ("lat","lon"))
          out_SDbkg         = nc_out.createVariable("SDbkg", "f8", ("lat","lon"))
          out_SDanl         = nc_out.createVariable("SDanl", "f8", ("lat","lon"))

          out_num           = nc_out.createVariable("nobs", "i", ("lat","lon"))

            # fill data variables
          print("... fill output vars")
          out_omb[:,:]         = omb2d
          out_oma[:,:]         = oma2d
          out_obserr[:,:]      = obserr2d
          out_obserr_orig[:,:] = obserr_orig2d
          out_obsval[:,:]      = obsval2d
          out_bkg[:,:]         = bkg2d
          out_anl[:,:]         = anl2d

          out_SDomb[:,:]         = SDomb2d
          out_SDoma[:,:]         = SDoma2d
          out_SDobserr[:,:]      = SDobserr2d
          out_SDobserr_orig[:,:] = SDobserr_orig2d
          out_SDobsval[:,:]      = SDobsval2d
          out_SDbkg[:,:]         = SDbkg2d
          out_SDanl[:,:]         = SDanl2d

          out_num[:,:]           = num2d

            # add attributes to each variable
                # long names
          print("... add longnames to vars")
          out_omb.long_name = "Observation-minus-Background (OMB)"
          out_oma.long_name = "Observation-minus-Analysis (OMA)"
          out_obserr.long_name = "Observation error post quality-control and bias-correction, if applicable"
          out_obserr_orig.long_name = "Original observation error from data producers (before quality-control and bias-correction)"
          out_obsval.long_name = "Observation value"
          out_bkg.long_name = "Background value"
          out_anl.long_name = "Analysis value"

          out_SDomb.long_name = "Standard deviation of "+str(out_omb.long_name)
          out_SDoma.long_name = "Standard deviation of "+str(out_oma.long_name)
          out_SDobserr.long_name = "Standard deviation of "+str(out_obserr.long_name)
          out_SDobserr_orig.long_name = "Standard deviation of "+str(out_obserr_orig.long_name)
          out_SDobsval.long_name = "Standard deviation of "+str(out_obsval.long_name)
          out_SDbkg.long_name = "Standard deviation of "+str(out_bkg.long_name)
          out_SDanl.long_name = "Standard deviation of "+str(out_anl.long_name)

          out_num.long_name = "Number of observations at each grid point"

          y.long_name = "Latitude"
          x.long_name = "Longitude"

                # units
          print("... add units to vars")
          out_omb.units = units
          out_oma.units = units
          out_obserr.units = units
          out_obserr_orig.units = units
          out_obsval.units = units
          out_bkg.units = units
          out_anl.units = units

          out_SDomb.units = units
          out_SDoma.units = units
          out_SDobserr.units = units
          out_SDobserr_orig.units = units
          out_SDobsval.units = units
          out_SDbkg.units = units
          out_SDanl.units = units

          y.units = "degrees-north"
          x.units = "degrees-east"

            # add global attributes
          print("... add global atts")
          nc_out.title = "Diag data converted to ["+str(degree)+"-degree x "+str(degree)+"-degree] grid"
          nc_out.notes = "The mean (or stddev) of each variable is computed for each grid point, and is then assigned to the corresponding grid point."
          nc_out.experiment_path = texptpath0
          nc_out.cycle = cdate

          nc_out.creation_date = str(dt.datetime.now())

          nc_out.close()

        #------------------------------------------

          del num2d, omb2d, SDomb2d, oma2d, SDoma2d, obserr2d, SDobserr2d, obserr_orig2d, SDobserr_orig2d, obsval2d, SDobsval2d, bkg2d, SDbkg2d, anl2d, SDanl2d

        del omb, oma, obserr, obserr_orig, obsval, qc, time, lats, lons#, bkg, anl

        idd += 1
    imm += 1
  iyy += 1

print("END LOOP")

###################################################################################################################
print("END: "+str(dt.datetime.now()))
print("========== END MAIN PROGRAM ==========")
###################################################################################################################
###################################################################################################################
