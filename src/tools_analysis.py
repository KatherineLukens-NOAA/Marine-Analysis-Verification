###########################################################################
# general_tools module
#
#	In support of NSST work
#
###########################################################################
#
# Import python modules
#

import xarray as xr
import os
import numpy as np
import glob
#import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import pickle
import datetime as dt
import sys
import warnings
import time

from scipy import asarray as ar,exp
from scipy.stats import norm 
from scipy.interpolate import griddata

###########################################################################
#
# FUNCTIONS 
#

pi = 4.0 * np.arctan(1.0)

#------------------------------------------------
# Get ocean model grid (lats, lons)
def get_ocngrid():

    # Open the NetCDF file using xarray
    #grid_fname = os.path.join('/scratch2/NCEPDEV/ocean/Guillaume.Vernieres/runs/low-res/soca_gridspec.nc')
    grid_fname = os.path.join('/scratch1/NCEPDEV/da/Katherine.Lukens/NSST/data/soca_gridspec.nc')

    ds = xr.open_dataset(grid_fname)

    lat = np.squeeze(ds['lat'][:])
    lon = np.squeeze(ds['lon'][:])

    ds.close()

    return lat.values, lon.values

#------------------------------------------------
# Regrid some variable ("in" grid) to ocean model grid ("out" grid)
#
# INPUT
#       var_in ....................... input variable to regrid to OUTPUT ocean model grid (lon_out, lat_out)
#       lat_in ....................... latitude array for var_in (1D)
#       lon_in ....................... longitude array for var_in (1D)
#       lat_out ...................... latitude array of OUTPUT ocean model grid (2D)
#       lon_out ...................... longitude array of OUTPUT ocean model grid (2D)
#
# OUTPUT
#       var_in_regrid ................ input variable that has been regridded to OUTPUT ocean model grid
#
def regrid(var_in, lat_in, lon_in, lat_out, lon_out, gridtype_in, gridtype_out):
    #clock1 = time.time()        # starting-time of task

    print("regrid: gridtype_in = "+str(gridtype_in)+" gridtype_out = "+str(gridtype_out))

    interp_method = 'linear' #'nearest'

    shape_var_in = np.shape(var_in)
    ndims_var_in = len(shape_var_in)
    print("regrid: var_in shape = "+str(shape_var_in))
    print("regrid: lat_in shape = "+str(np.shape(lat_in))+" min = "+str(np.min(lat_in))+" max = "+str(np.max(lat_in)))
    print("regrid: lon_in shape = "+str(np.shape(lon_in))+" min = "+str(np.min(lon_in))+" max = "+str(np.max(lon_in)))
    print("regrid: lat_out shape = "+str(np.shape(lat_out))+" min = "+str(np.min(lat_out))+" max = "+str(np.max(lat_out)))
    print("regrid: lon_out shape = "+str(np.shape(lon_out))+" min = "+str(np.min(lon_out))+" max = "+str(np.max(lon_out)))

    shape_lat = np.shape(lat_out)
    shape_lon = np.shape(lon_out)
    nlat_out = shape_lat[0]
    nlon_out = shape_lat[1]
    print("regrid: nlat_out, nlon_out = "+str(nlat_out)+" "+str(nlon_out))

    lon_in_max  = np.max(lon_in)
    lon_out_max = np.max(lon_out)

    lon_in = np.array(lon_in)

    if gridtype_in=="ocean model": #and gridtype_out!="atmos model":
        # get ocean model grid and rotate to match rectangular output grid [-180,180]
      lon_in[lon_in < -180] = lon_in[lon_in < -180] + 360
    elif gridtype_in=="atmos model" and gridtype_out=="ocean model" and lon_in_max > lon_out_max:
        # get atmos grid and rotate to match ocean model grid
      lon_in[lon_in >= lon_out_max] = lon_in[lon_in >= lon_out_max] - 360
    elif gridtype_in=="atmos model" and gridtype_out!="ocean model":
        # get atmos grid and rotate to match rectangular output grid [-180,180]
      lon_in[lon_in > 180] = lon_in[lon_in > 180] - 360
    print("regrid: lon_in rotated min = "+str(np.min(lon_in))+" max = "+str(np.max(lon_in)))

	# create 2D meshgrids
    lon_in_grid, lat_in_grid   = np.meshgrid(lon_in, lat_in)
    #lon_out_grid, lat_out_grid = np.meshgrid(lon_out, lat_out)
    print("regrid: lat_in_grid shape = "+str(np.shape(lat_in_grid)))
    print("regrid: lon_in_grid shape = "+str(np.shape(lon_in_grid)))

        # input for griddata function
    points_in = np.array([lon_in_grid.flatten(), lat_in_grid.flatten()]).T
    values_in = var_in.flatten()
    print("regrid: points_in size = "+str(np.size(points_in))+" shape = "+str(np.shape(points_in)))
    print("regrid: values_in size = "+str(np.size(values_in))+" shape = "+str(np.shape(values_in)))
    del lon_in_grid, lat_in_grid

    #lon_out_grid, lat_out_grid = np.meshgrid(lon_out, lat_out)
    lat_out_grid = lat_out
    lon_out_grid = lon_out
    gridpoints_out = np.array([lon_out_grid.flatten(), lat_out_grid.flatten()]).T
    shape_lon_out_grid = lon_out_grid.shape
    del lon_out_grid, lat_out_grid

        # regrid input variable to ocean grid
    if ndims_var_in==2:
            # 2D var_in
      #var_in_regrid = griddata((lon_in_grid.reshape(-1), lat_in_grid.reshape(-1)), np.squeeze(var_in.reshape(-1)), (lon_out, lat_out), method='nearest')
      tvar_in_regrid = griddata(points_in, values_in, gridpoints_out, method=interp_method)
      var_in_regrid = tvar_in_regrid.reshape(shape_lon_out_grid) #lon_out_grid.shape)
      del tvar_in_regrid

    del points_in, values_in, gridpoints_out, shape_lon_out_grid

    #elif ndims_var_in==3:
    #        # 3D var_in
    #  var_in_regrid = np.empty((shape_var_in[0], nlat_out, nlon_out), float)
    #  print("var_in_regrid shape = "+str(np.shape(var_in_regrid)))
    #  for i in range(shape_var_in[0]):
    #    tvar_in_regrid = griddata((lon_in_grid.reshape(-1), lat_in_grid.reshape(-1)), np.squeeze(var_in[i,:,:].reshape(-1)), (lon_out, lat_out), method='nearest')
    #    var_in_regrid[i,:,:] = tvar_in_regrid[:,:]
    #
    #elif ndims_var_in==4:
    #        # 4D var_in
    #  var_in_regrid = np.empty((shape_var_in[0], shape_var_in[1], nlat_out, nlon_out), float)
    #  for i in range(shape_var_in[0]):
    #    for ip in range(shape_var_in[1]):
    #        var_in_regrid[i,ip,:,:] = griddata((lon_in_grid.reshape(-1), lat_in_grid.reshape(-1)), np.squeeze(var_in[i,ip,:,:].reshape(-1)), (lon_out, lat_out), method='nearest')

    #clock2 = time.time()        # ending-time of task
    #clock = clock2 - clock1     # time elapsed (sec)
    #print('regrid: time elapsed = '+str(clock)+' seconds')

    return var_in_regrid

#===============================================================================================
# Get verification date (for analysis field) corresponding to given forecast (fhour)
#
# INPUT: Forecast date
# 
# OUTPUT: Verification date
#
def get_verif_date_anl(yyyy, mm, dd, fhour):

    # number of days in future that forecast hour (fhour) corresponds to
    #   e.g., if fhour==72, then nfdays==3 --> forecast fhour corresponds to 3 days in the future
    nfdays = int(int(fhour) / 24)
    #print("get_verif_date: nfdays = "+str(nfdays))

    if nfdays==0:
      yyyyV = yyyy
      mmV   = mm
      ddV   = dd
    else:
      # !!! Below is ONLY for experiments spanning 2021-07 to 2021-08
      yyyyV = yyyy

      tdd = int(dd) + nfdays
      if tdd>31:
        newdd = tdd - 31
        if int(mm)==7: mmV = "08"
        if int(mm)==8: mmV = "09"
      else:
        mmV = mm
        newdd = tdd
      #print("get_verif_date: newdd = "+str(newdd))

      if newdd<10:
        ddV = "0"+str(newdd)
      else:
        ddV = str(newdd)

    return yyyyV, mmV, ddV
    
#===============================================================================================
# Get verification date (for 6-hr background forecast field) corresponding to given forecast (fhour)
#
# INPUT: Forecast date
# 
# OUTPUT: Verification date
#
def get_verif_date_f06(yyyy, mm, dd, fhour):

    # number of days in future that forecast hour (fhour) corresponds to
    #   e.g., if fhour==72, then nfdays==3 --> forecast fhour corresponds to 3 days in the future
    nfdays = int(int(fhour) / 24)
    #print("get_verif_date: nfdays = "+str(nfdays))

    if nfdays==0:
      yyyyV = yyyy
      mmV   = mm
      ddV   = dd
    else:
      # !!! Below is ONLY for experiments spanning 2021-07 to 2021-08
      yyyyV = yyyy

      tdd = int(dd) + nfdays - 1		# when using 18z f06 as verification
      if tdd>31:
        newdd = tdd - 31
        if int(mm)==7: mmV = "08"
        if int(mm)==8: mmV = "09"
      else:
        mmV = mm
        newdd = tdd
      #print("get_verif_date: newdd = "+str(newdd))

      if newdd<10:
        ddV = "0"+str(newdd)
      else:
        ddV = str(newdd)

    return yyyyV, mmV, ddV

#===============================================================================================
# Compute Fit Function
#       Use given Gaussian parameters and variable x
#
def fit_func(x,a,mu,sigma,c):

    #"""gaussian function used for the fit"""

    return a * norm.pdf(x,loc=mu,scale=sigma) + c

#===============================================================================================
# Compute Guassian PDF
#	Use given Gaussian parameters and variable x
#
def gaus(x,c,x_mean,x_sd):

    G = c*exp(-(x - x_mean)**2 / (2*x_sd**2))

    return G

#===============================================================================================
# Correlate two arrays (Pearson correlation) 
#
# Source: https://www.geeksforgeeks.org/exploring-correlation-in-python/
#
def pearson_corr(X,Y):

    shapeX = np.shape(X); shapeY = np.shape(Y)
    print("pearson_corr: shape x = "+str(shapeX)+" y = "+str(shapeY))
    #if len(X)==len(Y):
    if shapeX==shapeY:
      if shapeX[0]>1:
        meanX = np.mean(X, axis=0)
        meanY = np.mean(Y, axis=0)
      #elif shapeX[0]==1:
      #  meanX = np.mean(X)
      #  meanY = np.mean(Y)
      Sum_xy = np.sum((X - meanX)*(Y - meanY), axis=0)
      #Sum_xy = np.sum((X[:,:,:] - meanX[:,:])*(Y[:,:,:] - meanY[:,:]), axis=0)

      Sum_x_squared = np.sum((X - meanX)**2, axis=0)
      Sum_y_squared = np.sum((Y - meanY)**2, axis=0)  
     
      sqrt_xy = np.sqrt(Sum_x_squared * Sum_y_squared)
      #print("pearson_corr: shape sqrt_xy = "+str(np.shape(sqrt_xy)))
      #print("pearson_corr: min/max sqrt_xy = "+str(np.nanmin(sqrt_xy))+" "+str(np.nanmax(sqrt_xy)))

      corr = Sum_xy / sqrt_xy
      print("pearson_corr: shape corr = "+str(np.shape(corr))+" X = "+str(np.shape(X)))

      return corr

    else: 
      print("No correlation performed. Exit program")
      sys.exit()

#===============================================================================================
# Interpolate from one rectangular grid to another rectangular grid
#	Like 'linint2' function from NCL
#	Source code from NCL Pivot to Python geocat.comp: https://github.com/NCAR/geocat-comp/blob/main/geocat/comp/interpolation.py
#
def interp_grid2grid(data_in, lat_in, lon_in, lat_out, lon_out, method: str="linear", fill_value: np.number=np.nan):

    data_in = xr.DataArray(data_in, dims=['lat','lon'], coords={'lat': lat_in, 'lon': lon_in})

    #output_coords = {'lat': lat_out, 'lon': lon_out}
    output_coords = {data_in.dims[-1]: lon_out, data_in.dims[-2]: lat_out,}

    data_out = data_in.interp(output_coords, method=method, kwargs={'fill_value': fill_value})

    return data_out

#===============================================================================================
# Interpolate atmos model grid to (MOM6) ocean model grid
#       Source code: https://github.com/NOAA-EMC/GDASApp/blob/a3c3c10a62ac54dde502598129f8b543197325a6/ush/socaincr2mom6.py#L77
#
# INPUT:
#	data_in = path to input file with variable to be regridded
#	invarname = variable name to be regridded, as specified in data_in
#	inlatname = latitude dim name to be regridded, as specified in data_in
#	inlonname = longitude dim name to be regridded, as specified in data_in
#
# OUTPUT:
#	outvar = input variable regridded to output (ocean) grid
#
def interp_to_ocngrid(data_in, invarname, inlatname, inlonname, ocngridfile, fill_value: np.number=np.nan):

	# get data to be regridded
    infile = xr.open_dataset(data_in)
    invar  = infile[invarname].values
    
	# get ocean grid
    ocngrid = xr.open_dataset(ocngridfile)
    ocnlats = ocngrid['lat'].values.squeeze()
    ocnlons = ocngrid['lon'].values.squeeze()

	# get atmos Gaussian grid and rotate to match ocean grid
    atmlats = infile[inlatname].values
    atmlons = infile[inlonname].values

    ocnlon_max = np.max(ocnlons)

    atmlons[atmlons >= ocnlon_max] = atmlons[atmlons >= ocnlon_max] - 360
    atmlongrid, atmlatgrid = np.meshgrid(atmlons, atmlats)

    outvar = griddata((atmlongrid.reshape(-1), atmlatgrid.reshape(-1)), np.squeeze(invar).reshape(-1), (atmlons, atmlats), method='nearest')

    return outvar

# -------------------------------------------------------------------------
# Find observation number density for a single dataset within grid cells as defined by (x,y), and conform resulting array to 2 dimensions
#
def var_to_2d_counts_1dset(xstr,ystr,tx,ty,txvar0,tyvar0):

        # convert lists to arrays
    x     = np.asarray(tx)
    y     = np.asarray(ty)
    xvar0 = np.asarray(txvar0)
    yvar0 = np.asarray(tyvar0)

        # loop to compute sums per grid cell (determined by x,y arrays)
                # 'nsum' contains the number of obs within each grid cell, as determined by
                # the conditions within np.where() in the loop.
    nsum0 = np.zeros([len(y),len(x)], dtype=float)
    for j in range(len(y)-1):
      if ystr.find("Pressure")!=-1 or ystr.find("Height")!=-1:
        ydiff = abs(y[j+1] - y[j])/2.0                  # half the difference between each gridded z value
        ymin  = y[j]-ydiff
        ymax  = y[j]+ydiff
        del ydiff
      else:
        ymin  = y[j]
        ymax  = y[j+1]

                # 'xidx' is an array containing the number of obs that satisfy the conditions within np.where()
      xidx0 = np.zeros([len(x)], dtype=int)
      if xstr.find("Time")!=-1:
        xidx0 = [np.size(np.where((xvar0>=x[k])*(yvar0>=ymin)*(yvar0<ymax))) for k in range(len(x))]
      else:
        xidx0 = [np.size(np.where((xvar0>=x[k])*(xvar0<x[k+1])*(yvar0>=ymin)*(yvar0<ymax))) for k in range(len(x)-1)]
                        # compute density for cyclic grid point
        xidx0_last = [np.size(np.where((xvar0>=x[len(x)-1])*(yvar0>=ymin)*(yvar0<ymax)))]
        xidx0 = np.append(xidx0,xidx0_last,axis=0)
        del xidx0_last

      sizeidx0 = np.size(xidx0)                 # size of xidx
      nsum0[j,0:sizeidx0] = xidx0               # write 'xidx' to each j dimension of 'nsum'

      del xidx0,sizeidx0
      del ymin,ymax

        # mean per grid cell
    warnings.filterwarnings('ignore', category=RuntimeWarning)

        # fill grid cells where nsum=0 with nan
    nsum_nan0 = np.where(nsum0==0,np.nan,nsum0)

    return nsum_nan0            # return counts per grid cell

# -------------------------------------------------------------------------
# Conform 1 variable to 2 dimensions
#
# Input:
#       xstr .......................... name of new x dimension
#       ystr .......................... name of new y dimension
#       x ............................. new x dimension
#       y ............................. new y dimension
#       xdim .......................... 1d x dimension variable to be transformed
#       ydim ........................., 1d y dimension variable to be transformed
#       var0 .......................... variable from dataset to be transformed to 2d grid (x,y)
#
# Output:
#       mean0 ......................... mean of var0 on 2d grid
#       SDx2d ......................... SD of var0 on 2d grid
#
def var1_to_2d(xstr,ystr,x,y,xdim,ydim,var0):

        # sort y to be in ascending order
    y.sort()

        # loop to compute sums per grid cell (determined by x,y arrays)
    sumx2d = np.zeros([len(y),len(x)], dtype=float)
    SDx2d  = np.zeros([len(y),len(x)], dtype=float)
    nsum   = np.zeros([len(y),len(x)], dtype=float)
    for j in range(len(y)-1):
      print("var_to_2d: j = "+str(j)+"/"+str(len(y)-1))
      if ystr.find("Pressure")!=-1 or ystr.find("Height")!=-1:
        ydiff = abs(y[j+1] - y[j])/2.0                  # half the difference between each gridded z value
        ymin  = y[j]-ydiff
        ymax  = y[j]+ydiff
        del ydiff
      else:
        ymin  = y[j]
        ymax  = y[j+1]

      for k in range(len(x)-1):
        if xstr.find("Time")!=-1:
          xidx = np.where((xdim==x[k])*(ydim>=ymin)*(ydim<ymax))
        else:
          xidx = np.where((xdim>=x[k])*(xdim<x[k+1])*(ydim>=ymin)*(ydim<ymax))
        sumx2d[j,k] = np.sum(var0[xidx])
        SDx2d[j,k]  = np.std(var0[xidx])
        nsum[j,k]   = np.size(xidx)
        del xidx

		# compute density for cyclic grid point
      xidx = np.where((xdim>=x[len(x)-1])*(ydim>=ymin)*(ydim<ymax))
      sumx2d[j,k] += np.sum(var0[xidx])
      SDx2d[j,k]  += np.std(var0[xidx])
      nsum[j,k]   += np.size(xidx)
      del xidx
      del ymin,ymax

        # mean per grid cell
    warnings.filterwarnings('ignore', category=RuntimeWarning)
    mean0    = sumx2d / nsum

    return mean0, SDx2d
    
# -------------------------------------------------------------------------
# Conform 2 variables to 2 dimensions simultaneously
#
# Input:
#       xstr .......................... name of new x dimension
#       ystr .......................... name of new y dimension
#       x ............................. new x dimension
#       y ............................. new y dimension
#       xdim .......................... 1d x dimension variable to be transformed
#       ydim ........................., 1d y dimension variable to be transformed
#       var0 .......................... variable from dataset 0 to be transformed to 2d grid (x,y)
#       var1 .......................... variable from dataset 1 to be transformed to 2d grid (x,y)
#
# Output:
#       mean0 ......................... mean of var0 on 2d grid
#       SDx2d0 ........................ SD of var0 on 2d grid
#       mean1 ......................... mean of var1 on 2d grid
#       SDx2d1 ........................ SD of var1 on 2d grid
#
# NOTE: var0 and var1 must be the same size!
#
def var2_to_2d(xstr,ystr,x,y,xdim,ydim,var0,var1):

        # sort y to be in ascending order
    y.sort()

        # loop to compute sums per grid cell (determined by x,y arrays)
    sumx2d0 = np.zeros([len(y),len(x)], dtype=float)
    SDx2d0  = np.zeros([len(y),len(x)], dtype=float)
    sumx2d1 = np.zeros([len(y),len(x)], dtype=float)
    SDx2d1  = np.zeros([len(y),len(x)], dtype=float)
    nsum    = np.zeros([len(y),len(x)], dtype=float)
    for j in range(len(y)-1):
      print("var2_to_2d: j = "+str(j)+"/"+str(len(y)-1))
      if ystr.find("Pressure")!=-1 or ystr.find("Height")!=-1:
        ydiff = abs(y[j+1] - y[j])/2.0                  # half the difference between each gridded z value
        ymin  = y[j]-ydiff
        ymax  = y[j]+ydiff
        del ydiff
      else:
        ymin  = y[j]
        ymax  = y[j+1]

      for k in range(len(x)-1):
        if xstr.find("Time")!=-1:
          xidx = np.where((xdim==x[k])*(ydim>=ymin)*(ydim<ymax))
        else:
          xidx = np.where((xdim>=x[k])*(xdim<x[k+1])*(ydim>=ymin)*(ydim<ymax))

        sumx2d0[j,k] = np.sum(var0[xidx])
        SDx2d0[j,k]  = np.std(var0[xidx])
        sumx2d1[j,k] = np.sum(var1[xidx])
        SDx2d1[j,k]  = np.std(var1[xidx])

        nsum[j,k]   = np.size(xidx)
        del xidx

                # compute density for cyclic grid point
      xidx = np.where((xdim>=x[len(x)-1])*(ydim>=ymin)*(ydim<ymax))

      sumx2d0[j,k] += np.sum(var0[xidx])
      SDx2d0[j,k]  += np.std(var0[xidx])
      sumx2d1[j,k] += np.sum(var1[xidx])
      SDx2d1[j,k]  += np.std(var1[xidx])

      nsum[j,k]   += np.size(xidx)
      del xidx
      del ymin,ymax

        # mean per grid cell
    warnings.filterwarnings('ignore', category=RuntimeWarning)
    mean0    = sumx2d0 / nsum
    mean1    = sumx2d1 / nsum

    return mean0, SDx2d0, mean1, SDx2d1

# -------------------------------------------------------------------------
# Conform 3+ variables to 2 dimensions simultaneously
#
# Input:
#       xstr .......................... name of new x dimension
#       ystr .......................... name of new y dimension
#       x ............................. new x dimension
#       y ............................. new y dimension
#       xdim .......................... 1d x dimension variable to be transformed
#       ydim ........................., 1d y dimension variable to be transformed
#       var0 .......................... variable from dataset to be transformed to 2d grid (x,y)
#       var1 .......................... variable from dataset to be transformed to 2d grid (x,y)
#       var2 .......................... variable from dataset to be transformed to 2d grid (x,y)
#       var3 .......................... variable from dataset to be transformed to 2d grid (x,y)
#
# Output:
#       mean0 ......................... mean of var0 on 2d grid
#       SDx2d0 ........................ SD of var0 on 2d grid
#       mean1 ......................... mean of var1 on 2d grid
#       SDx2d1 ........................ SD of var1 on 2d grid
#       mean2 ......................... mean of var2 on 2d grid
#       SDx2d2 ........................ SD of var2 on 2d grid
#
# NOTE: var0 and var1 must be the same size!
#
def vars_to_2d(xstr,ystr,x,y,xdim,ydim,topt,var0,var1,var2,var3,var4,var5,var6):

    opt = int(topt)

        # sort y to be in ascending order
    y.sort()

        # loop to compute sums per grid cell (determined by x,y arrays)
    sumx2d0 = np.zeros([len(y),len(x)], dtype=float)
    sumx2d1 = np.zeros([len(y),len(x)], dtype=float)
    sumx2d2 = np.zeros([len(y),len(x)], dtype=float)
    sumx2d3 = np.zeros([len(y),len(x)], dtype=float)
    sumx2d4 = np.zeros([len(y),len(x)], dtype=float)
    sumx2d5 = np.zeros([len(y),len(x)], dtype=float)
    sumx2d6 = np.zeros([len(y),len(x)], dtype=float)

    SDx2d0  = np.zeros([len(y),len(x)], dtype=float)
    SDx2d1  = np.zeros([len(y),len(x)], dtype=float)
    SDx2d2  = np.zeros([len(y),len(x)], dtype=float)
    SDx2d3  = np.zeros([len(y),len(x)], dtype=float)
    SDx2d4  = np.zeros([len(y),len(x)], dtype=float)
    SDx2d5  = np.zeros([len(y),len(x)], dtype=float)
    SDx2d6  = np.zeros([len(y),len(x)], dtype=float)
   
    nsum    = np.zeros([len(y),len(x)], dtype=float)
    
    for j in range(len(y)-1):
      print("vars_to_2d: j = "+str(j)+"/"+str(len(y)-1))
      if ystr.find("Pressure")!=-1 or ystr.find("Height")!=-1:
        ydiff = abs(y[j+1] - y[j])/2.0                  # half the difference between each gridded z value
        ymin  = y[j]-ydiff
        ymax  = y[j]+ydiff
        del ydiff
      else:
        ymin  = y[j]
        ymax  = y[j+1]

      for k in range(len(x)-1):
        if xstr.find("Time")!=-1:
          xidx = np.where((xdim==x[k])*(ydim>=ymin)*(ydim<ymax))
        else:
          xidx = np.where((xdim>=x[k])*(xdim<x[k+1])*(ydim>=ymin)*(ydim<ymax))

        sumx2d0[j,k] = np.sum(var0[xidx])
        sumx2d1[j,k] = np.sum(var1[xidx])
        sumx2d2[j,k] = np.sum(var2[xidx])
        sumx2d3[j,k] = np.sum(var3[xidx])
        sumx2d4[j,k] = np.sum(var4[xidx])
        sumx2d5[j,k] = np.sum(var5[xidx])
        sumx2d6[j,k] = np.sum(var6[xidx])

        SDx2d0[j,k]  = np.std(var0[xidx])
        SDx2d1[j,k]  = np.std(var1[xidx])
        SDx2d2[j,k]  = np.std(var2[xidx])
        SDx2d3[j,k]  = np.std(var3[xidx])
        SDx2d4[j,k]  = np.std(var4[xidx])
        SDx2d5[j,k]  = np.std(var5[xidx])
        SDx2d6[j,k]  = np.std(var6[xidx])

        nsum[j,k]    = np.size(xidx)

        del xidx

                # compute density for cyclic grid point
      xidx = np.where((xdim>=x[len(x)-1])*(ydim>=ymin)*(ydim<ymax))

      sumx2d0[j,k] += np.sum(var0[xidx])
      sumx2d1[j,k] += np.sum(var1[xidx])
      sumx2d2[j,k] += np.sum(var2[xidx])
      sumx2d3[j,k] += np.sum(var3[xidx])
      sumx2d4[j,k] += np.sum(var4[xidx])
      sumx2d5[j,k] += np.sum(var5[xidx])
      sumx2d6[j,k] += np.sum(var6[xidx])

      SDx2d0[j,k]  += np.std(var0[xidx])
      SDx2d1[j,k]  += np.std(var1[xidx])
      SDx2d2[j,k]  += np.std(var2[xidx])
      SDx2d3[j,k]  += np.std(var3[xidx])
      SDx2d4[j,k]  += np.std(var4[xidx])
      SDx2d5[j,k]  += np.std(var5[xidx])
      SDx2d6[j,k]  += np.std(var6[xidx])

      nsum[j,k]    += np.size(xidx)

      del xidx
      del ymin,ymax

        # mean per grid cell
    warnings.filterwarnings('ignore', category=RuntimeWarning)
    mean0    = sumx2d0 / nsum
    mean1    = sumx2d1 / nsum
    mean2    = sumx2d2 / nsum
    mean3    = sumx2d3 / nsum
    mean4    = sumx2d4 / nsum
    mean5    = sumx2d5 / nsum
    mean6    = sumx2d6 / nsum

    if opt==0:
        # return means and SDs
      print("vars_to_2d: return means")
      return nsum, mean0, SDx2d0, mean1, SDx2d1, mean2, SDx2d2, mean3, SDx2d3, mean4, SDx2d4, mean5, SDx2d5, mean6, SDx2d6
    elif opt==1:
        # return sums and SDs
      print("vars_to_2d: return sums")
      return nsum, sumx2d0, SDx2d0, sumx2d1, SDx2d1, sumx2d2, SDx2d2, sumx2d3, SDx2d3, sumx2d4, SDx2d4, sumx2d5, SDx2d5, sumx2d6, SDx2d6

#===============================================================================================
# Calculate depth of each ocean layer
#
def get_ocean_depths(th):

    h = th[0,:,:,:]
    #print("shape h = "+str(np.shape(h)))

	# total depth of ocean at each lat,lon point
	#	Result: 2D array (lat, lon)
    total_depth = np.sum(h, axis=0)
    print("total depth = "+str(total_depth))
    print("shape total depth = "+str(np.shape(total_depth)))    
    
    shapeh = np.shape(h)
    
    depths = np.empty(shapeh,dtype=float)
    for ih in range(shapeh[0]):
      if ih==0:
        depths[ih,:,:] = h[ih,:,:] / 2.0
      else:
        depths[ih,:,:] = depths[ih-1,:,:] + (h[ih,:,:] / 2.0)
    print("depths shape = "+str(np.shape(depths)))
    print("depths min/max = "+str(np.min(depths))+" "+str(np.max(depths)))

    return depths

#===============================================================================================
# Calculate anomaly correlation
#
#   Source: $NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl --> pattern_cor2
#
def anom_corr(x, y):

    # make x and y 1D arrays
    xarr = np.asarray(x); del x
    yarr = np.asarray(y); del y
    x1d = xarr.flatten()
    y1d = yarr.flatten()

    xsum = np.sum(x1d)
    ysum = np.sum(y1d)

    xanom = x1d - xsum
    yanom = y1d - ysum

    xycov = np.sum(xanom * yanom)

    xanom2 = np.sum(xanom**2)
    yanom2 = np.sum(yanom**2)

    if xanom2>0 and yanom2>0:
      r = xycov / (np.sqrt(xanom2) * np.sqrt(yanom2))
    else:
      r = -999

    del x1d, y1d, xsum, ysum, xanom, yanom, xycov, xanom2, yanom2

    return r




#===============================================================================================
