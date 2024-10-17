###########################################################################
# daynight_terminator module
#
#	In support of NSST work
#
# Source: https://github.com/JoGall/terminator/blob/master/terminator.R
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

import math

###########################################################################
#
# FUNCTIONS 
#

pi = 4.0 * np.arctan(1.0)

#===============================================================================================
# Convert datetime to seconds since epoch (1 Jan 1970)
#
#	time_delta ... units = days
#
def time_since_epoch(yyyy, mm, dd, hh, mn):

    t_dt = dt.datetime(int(yyyy),int(mm),int(dd),int(hh),int(mn))
    epoch_time = dt.datetime(1970,1,1)

    time_delta = t_dt - epoch_time
    time_delta_sec = time_delta.total_seconds()		# units = seconds
    print("time since epoch = "+str(time_delta_sec))

    del t_dt, epoch_time, time_delta

    return time_delta_sec

#===============================================================================================
# Convert radians to degrees
#
def rad2deg(rad):

    deg = (rad * 180.0) / pi

    return deg

#===============================================================================================
# Convert degrees to radians
#
def deg2rad(deg):

    rad = (deg * pi) / 180.0

    return rad

#===============================================================================================
# Get Julian day
#	Continuous count of days since the beginning of the Julian period, i.e., number of days since
#		noon on January 1, 4713 BC; 2440587.5 is number of days between Julian epoch and UNIX epoch
#
#	time ... units = seconds since 1 Jan 1970
#	jday ... units = days
#
def get_julian_day(time):

    jday = (int(time) / 86400) + 2440587.5
    #print("jday = "+str(jday))

    return jday

#===============================================================================================
# Calculate Greenwich Mean Sidereal Time (GMST): hour angle of the average position of the vernal equinox
#
def get_gmst(jday):

    t_gmst = jday - 2451545.0
    gmst = (18.697374558 + 24.06570982441908 * t_gmst) % 24
    #print("gmst = "+str(gmst))

    del t_gmst

    return gmst

#===============================================================================================
# Compute the position of Sun in ecliptic coordinates
#
def sun_ecliptic_pos(jday):

	# days since start of J2000.0
    n = jday - 2451545.0

	# mean longitude of sun
    t_mlon = 280.460 + (0.9856474 * n)
    mlon = t_mlon % 360
    del t_mlon

	# mean anomlay of sun
    t_asun = 357.528 + (0.9856003 * n)
    asun = t_asun % 360
    del t_asun

	# ecliptic longitude of sun
    elon = mlon + (1.915 * np.sin(deg2rad(asun))) + (0.02 * np.sin(2.0 * deg2rad(asun)))

	# distance from sun in AU
    xsun = 1.00014 - (0.01671 * np.cos(deg2rad(asun))) - (0.0014 * np.cos(2.0 * deg2rad(asun)))
 
    return elon, xsun

#===============================================================================================
# Compute ecliptic obliquity
#
def sun_ecliptic_obliq(jday):

    	# days since start of J2000.0
    n = jday - 2451545.0

    	# Julian centuries since J2000.0
    jcent = n / 36525

	# compute epsilon
    epsilon = 23.43929111 - jcent * ((46.836769 / 3600) - jcent * ((0.0001831 / 3600) + jcent * ((0.00200340 / 3600) - jcent * ((0.576e-6 / 3600) - jcent * (4.34e-8 / 3600)))))

    return epsilon

#===============================================================================================
# Compute the Sun's equatorial position from its ecliptic position
#
def sun_equatorial_pos(seclip_pos, eclip_obliq):

    alpha = rad2deg(np.arctan(np.cos(deg2rad(eclip_obliq)) * np.tan(deg2rad(seclip_pos))))
    delta = rad2deg(np.arcsin(np.sin(deg2rad(eclip_obliq)) * np.sin(deg2rad(seclip_pos))))

    tlquad = seclip_pos/90
    trquad = alpha/90

    lquad = math.floor(tlquad) * 90
    rquad = math.floor(trquad) * 90

    alpha = alpha + (lquad - rquad)

    return alpha, delta

#===============================================================================================
# Calculate hour angle of Sun for a longitude on Earth 
#
def hour_angle(lon, salpha, gmst, mm):

    t_angle = gmst + (lon / 15)

	# boreal winter
    angle = (t_angle * 15) - salpha
	# boreal summer
    ##angle = (t_angle * 15) + salpha

    return angle

#===============================================================================================
# Compute latitude of the day/night terminator for a given hour angle and sun position
#
def lat_terminator(angle, sdelta):

    term_lat = rad2deg(np.arctan(-np.cos(deg2rad(angle)) / np.tan(deg2rad(sdelta))))

    return term_lat

#===============================================================================================
# Calculate latitude and longitude of the day/night terminator with specified range using time (datetime in seconds since 1 Jan 1970)
#
#	lonmin = -180, lonmax = 180
#
#	time ... units = days
#
def terminator(yyyy,mm,dd,hh,mn, lonmin, lonmax, nlons):

    time_seconds = time_since_epoch(yyyy,mm,dd,hh,mn)

    jday = get_julian_day(time_seconds)		# time units = seconds since epoch
    
    gmst = get_gmst(jday)

    seclip_pos, dist_sun = sun_ecliptic_pos(jday)
    eclip_obliq          = sun_ecliptic_obliq(jday)
    salpha, sdelta       = sun_equatorial_pos(seclip_pos, eclip_obliq)

    flons = np.linspace(lonmin, lonmax, nlons)

    flats = np.nan*np.ones_like(flons)
    for k in range(np.size(flons)):
      ha = hour_angle(flons[k], salpha, gmst, mm)
      print("hour angle = "+str(ha))
      ##longitude = flons[k] + ha
      ##tlat = lat_terminator(longitude, sdelta)
      tlat = lat_terminator(ha, sdelta)
      flats[k] = tlat
      del ha,tlat#,longitude

    del time_seconds, jday, gmst
    del seclip_pos, dist_sun, eclip_obliq, salpha, sdelta

    return flats, flons
    


#===============================================================================================
