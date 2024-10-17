###########################################################################
# plotting_tools module
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
import warnings

import cartopy
cartopy.config['data_dir'] = '/home/Katherine.Lukens/.local/share/cartopy'
import cartopy.crs as ccrs

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
from matplotlib.cm import get_cmap

import scipy.stats as ss
from scipy.stats import norm 
from scipy.stats import levy_stable
from scipy.optimize import curve_fit
from scipy.stats import ttest_ind

import pickle
import datetime as dt
import sys

from tools_analysis import gaus
from tools_analysis import fit_func
from tools_analysis import var_to_2d_counts_1dset

###########################################################################
#
# FUNCTIONS 
#

#````````````````````````````````````````
# Constants

pi = 4.0 * np.arctan(1.0)

#````````````````````````````````````````
# Plotting specs

fonttitle  = 18
fontaxis = 10
fontlegend = 12

#===============================================================================================
# Plot map of differences (e.g., anomalies)
#
#	Optional Arguments: pltmin, pltmax
#		... uses default unless specified in function call
#
def plot_map(x, y, z, zmean, title, sunlat, sunlon, plottypeflag, pltmin=-999, pltmax=999):

    print("plot_map: shape x = "+str(np.shape(x))+" y = "+str(np.shape(y))+" z = "+str(np.shape(z)))

    	# create plot
    #print("create plot")
    fig = plt.figure(figsize=(10,7.5))
	# create map
    proj  = ccrs.PlateCarree()
    axmap = plt.axes(projection=proj)

	# set colormap and color range values
    varmin = np.nanmin(z)
    varmax = np.nanmax(z)
    if plottypeflag.find("diff")!=-1 or plottypeflag.find("Diff")!=-1 or plottypeflag.find("Sign")!=-1 or plottypeflag.find("sign")!=-1:
      ppcmap = "seismic"
      #ppcmap = "RdBu_r"
      if pltmin==-999:
        if varmax!=np.nan and varmax<1:   pltmax=0.5
        if varmax!=np.nan and varmax>=1:  pltmax=1 
        if varmax!=np.nan and varmax>=5:  pltmax=5 
        if varmax!=np.nan and varmax>=10: pltmax=10
        if varmax!=np.nan and varmax>=15: pltmax=15
        if varmax!=np.nan and varmax>=20: pltmax=25
        if varmax!=np.nan and varmax>=30: pltmax=50
        pltmin = -pltmax
    elif plottypeflag.find("SD")!=-1:
      ppcmap = "YlOrRd" #"hot"
      if pltmin==-999:
        if varmax!=np.nan and varmax>=1:  pltmax=1  ; pltmin = 0
        if varmax!=np.nan and varmax>=5:  pltmax=10 ; pltmin = 1
        if varmax!=np.nan and varmax>=10: pltmax=15 ; pltmin = 1
        if varmax!=np.nan and varmax>=20: pltmax=25 ; pltmin = 1
        if varmax!=np.nan and varmax>=30: pltmax=50 ; pltmin = 1
    else:
      #ppcmap = "jet"
      ppcmap = "nipy_spectral"
      if pltmin==-999 and varmin<0:
        if varmax!=np.nan and varmax<1:   pltmax=0.5
        if varmax!=np.nan and varmax>=1:  pltmax=1 
        if varmax!=np.nan and varmax>=5:  pltmax=5 
        if varmax!=np.nan and varmax>=10: pltmax=10
        if varmax!=np.nan and varmax>=15: pltmax=15
        if varmax!=np.nan and varmax>=20: pltmax=25
        if varmax!=np.nan and varmax>=30: pltmax=50
        pltmin = -pltmax
      elif pltmin==-999 and varmin>=0:
        if varmax!=np.nan and varmax>=1:  pltmax=1  ; pltmin = 0
        if varmax!=np.nan and varmax>=5:  pltmax=10 ; pltmin = 1
        if varmax!=np.nan and varmax>=10: pltmax=15 ; pltmin = 1
        if varmax!=np.nan and varmax>=20: pltmax=25 ; pltmin = 1
        if varmax!=np.nan and varmax>=30: pltmax=50 ; pltmin = 1
    #print("plot_map: ppcmap = "+str(ppcmap))

	# plot data on map
    #axmap = plt.axes(projection=proj)
    pp = plt.pcolormesh(x,y,z,cmap=ppcmap,vmin=pltmin,vmax=pltmax,transform=proj)

	# add ccastlines
    axmap.coastlines()

	# add grid lines
    gl = axmap.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.7, linestyle='--')
    gl.xlines = False            # plot longitudes
    gl.ylines = False            # plot latitudes
    #gl.xlines = True
    #gl.ylines = True
        	# set GLOBAL domain limits
    latmin = -90.0
    latmax = 90.0
    lonmin = -180.0
    lonmax = 180.0
    clat   = (latmax + latmin) / 2   
    clon   = (lonmax + lonmin) / 2
    axmap.set_extent([lonmin,lonmax,latmin,latmax], crs=proj)
                	# define x-axis tickmarks
    gl.xlocator = mticker.FixedLocator([-180,-120,-60,0,60,120,180])
    #gl.xlocator = mticker.FixedLocator([-160,-140,-120,-100,-80,-60])
                	# define y-axis tickmarks
    gl.ylocator = mticker.FixedLocator([-90,-60,-30,0,30,60,90])
    #gl.ylocator = mticker.FixedLocator([-45,-30,-15,0,15,30,45])

    gl.xlabel_style = {'size': fontaxis, 'color': 'black'}#'gray'}
    gl.ylabel_style = {'size': fontaxis, 'color': 'black'}#'gray'}

    gl.xlabels_top=False
    gl.ylabels_right=False
    #gl.top_labels=None
    #gl.right_labels=None

	# overlay daily mean contours
    if zmean != "none":
	# for anomalies
      linemin = pltmin #-50
      linemax = pltmax #50
      if abs(pltmax)<1:
        #pltlevels = [-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5]
        pltlevels = [-1, 0, 1]
      elif abs(pltmax)<10:
        linestride = 1
      elif abs(pltmax)<20:
        linestride = 2
      elif abs(pltmax)<50:
        linestride = 5
      else:
        linestride = 10
      #cc = plt.contour(x,y,zmean,colors="navy",levels=range(linemin,linemax,linestride),linewidth=1)
      cc = plt.contour(x,y,zmean,colors="navy",levels=range(-1,1,1),linewidth=1)
      #cc = plt.contour(x,y,zmean,colors="navy",levels=np.linspace(linemin,linemax,num=linestride*10+1),linewidth=1)
	# for increments
      #cc = plt.contour(x,y,zmean,colors="navy",levels=[-0.1,0.1],linewidth=1)

      axmap.clabel(cc,inline=1,fontsize=fontlegend,fmt='%1.0f')

	# add location of Sun
    if sunlon != -999 and sunlat != -999:
      axmap.scatter(sunlon, sunlat, c="black", s=150, marker="*", alpha=0.75, label="Sun Location")
      leg = plt.legend(loc='lower right', prop={'size': 8})

	# add colorbar to plot
    if plottypeflag.find("Sign")!=-1 or plottypeflag.find("sign")!=-1:		# IS sign change plot
     redlabel = "OmA>0 and OmB<0"
     bluelabel = "OmA<0 and OmB>0"
     custom_markers = [Line2D([0],[0],color="red",label=redlabel,linestyle="None",marker="s"), Line2D([0],[0],color="blue",label=bluelabel,linestyle="None",marker="s")]
     legP1 = axmap.legend(handles=custom_markers, loc='lower center', bbox_to_anchor=(0.5, 0.), prop={'size': fontlegend})
     for lhP1 in legP1.legendHandles:
       lhP1.set_alpha(1)
     legendmarkersize = 30
     for lhP1 in range(np.size(legP1.legendHandles)):
       legP1.legendHandles[lhP1]._sizes = [legendmarkersize]
    elif plottypeflag.find("Sign")==-1 or plottypeflag.find("sign")==-1:		# IS NOT sign change plot
      cbarstr = 'SST (deg-C)'
      cbar = fig.colorbar(pp, label=cbarstr, orientation='horizontal', pad=0.05)#0.125)
      cbar.ax.tick_params(labelsize=fontlegend)

	# add title
    axmap.set_title(title)

    return fig
    
#===============================================================================================
# Plot 2D cross-section
#
#       Optional Arguments: pltmin, pltmax
#               ... uses default unless specified in function call
#
def plot_cross_section(x, y, z, zmean, xvals, yvals, xlabel, ylabel, title, plottypeflag, pltmin=-999, pltmax=999):

        # create plot
    print("create plot")
    fig = plt.figure(figsize=(10,7.5))

    ax = fig.add_subplot()

        # set axes labels
    ax.set_xlabel(xlabel, fontsize=fontaxis)
    ax.set_ylabel(ylabel, fontsize=fontaxis)

            # x-axis
    if (abs(np.nanmin(xvals))+abs(np.nanmax(xvals))) < 200:
      denom = 180        # xvals is LATITUDE
    else:
      denom = 360       # xvals is LONGITUDE

    nspacing = int(30 * (np.size(xvals)/denom))
    txlabels = xvals[::nspacing]
    del nspacing
                    # print minus sign (not hyphen) in front of negative numbers on axis
    nspacing = int(30 * (np.size(x)/denom))
    xvalues = np.arange(0,np.size(x),nspacing)
    xlabels = [[] for i in range(np.size(x))]
    i=0
    j=0
    while (i < np.size(x)):
      xlabels[i] = str(txlabels[j])
      i += nspacing
      j += 1
    ax.set_xticks(xvalues)
    ax.set_xticklabels(xlabels[::nspacing], fontsize=fontaxis-2)
    del xvalues,txlabels,xlabels,nspacing

            # y-axis
    nspacing = 5
    yvalues = np.arange(0,np.size(yvals),nspacing)
    ylabels = yvals[::nspacing]
    ax.set_yticks(yvalues)
    ax.set_yticklabels(ylabels, fontsize=fontaxis-2)
    del yvalues,ylabels,nspacing

        # set colormap and color range values
    varmin = np.nanmin(z)
    varmax = np.nanmax(z)
    if plottypeflag.find("diff")!=-1 or plottypeflag.find("Diff")!=-1 or plottypeflag.find("Sign")!=-1 or plottypeflag.find("sign")!=-1:
      ppcmap = "seismic"
      if pltmin==-999:
        if varmax!=np.nan and varmax<1:   pltmax=0.5
        if varmax!=np.nan and varmax>=1:  pltmax=1
        if varmax!=np.nan and varmax>=5:  pltmax=5
        if varmax!=np.nan and varmax>=10: pltmax=10
        if varmax!=np.nan and varmax>=15: pltmax=15
        if varmax!=np.nan and varmax>=20: pltmax=25
        if varmax!=np.nan and varmax>=30: pltmax=50
        pltmin = -pltmax
    elif plottypeflag.find("SD")!=-1:
      ppcmap = "YlOrRd" #"hot"
      if pltmin==-999:
        if varmax!=np.nan and varmax>=1:  pltmax=1  ; pltmin = 0
        if varmax!=np.nan and varmax>=5:  pltmax=10 ; pltmin = 1
        if varmax!=np.nan and varmax>=10: pltmax=15 ; pltmin = 1
        if varmax!=np.nan and varmax>=20: pltmax=25 ; pltmin = 1
        if varmax!=np.nan and varmax>=30: pltmax=50 ; pltmin = 1
    else:
      #if plottypeflag.find("bkgerr")!=-1:
      #  ppcmap = "jet"
      #else:
      ppcmap = "nipy_spectral"
      if pltmin==-999 and varmin<0:
        if varmax!=np.nan and varmax<1:   pltmax=0.5
        if varmax!=np.nan and varmax>=1:  pltmax=1
        if varmax!=np.nan and varmax>=5:  pltmax=5
        if varmax!=np.nan and varmax>=10: pltmax=10
        if varmax!=np.nan and varmax>=15: pltmax=15
        if varmax!=np.nan and varmax>=20: pltmax=25
        if varmax!=np.nan and varmax>=30: pltmax=50
        pltmin = -pltmax
      elif pltmin==-999 and varmin>=0:
        if varmax!=np.nan and varmax>=1:  pltmax=1  ; pltmin = 0
        if varmax!=np.nan and varmax>=5:  pltmax=10 ; pltmin = 1
        if varmax!=np.nan and varmax>=10: pltmax=15 ; pltmin = 1
        if varmax!=np.nan and varmax>=20: pltmax=25 ; pltmin = 1
        if varmax!=np.nan and varmax>=30: pltmax=50 ; pltmin = 1

        # plot data on map
    originstr = "lower"
    #pp = plt.imshow(z, cmap=ppcmap, vmin=pltmin, vmax=pltmax, origin=originstr, aspect='auto')
    pp = plt.imshow(z, cmap=ppcmap, vmin=pltmin, vmax=pltmax, origin=originstr, aspect='auto')

        # add colorbar to plot
    if plottypeflag.find("Sign")!=-1 or plottypeflag.find("sign")!=-1:          # IS sign change plot
      redlabel = "OmA>0 and OmB<0"
      bluelabel = "OmA<0 and OmB>0"
      custom_markers = [Line2D([0],[0],color="red",label=redlabel,linestyle="None",marker="s"), Line2D([0],[0],color="blue",label=bluelabel,linestyle="None",marker="s")]
      legP1 = ax.legend(handles=custom_markers, loc='lower center', bbox_to_anchor=(0.5, 0.), prop={'size': fontlegend})
      for lhP1 in legP1.legendHandles:
        lhP1.set_alpha(1)
      legendmarkersize = 30
      for lhP1 in range(np.size(legP1.legendHandles)):
        legP1.legendHandles[lhP1]._sizes = [legendmarkersize]

    elif plottypeflag.find("Sign")==-1 or plottypeflag.find("sign")==-1:                # IS NOT sign change plot
      #cbarstr = 'SST'
      #cbar = fig.colorbar(pp, label=cbarstr, orientation='horizontal', pad=0.05)#0.125)
      cbar = fig.colorbar(pp, orientation='horizontal', pad=0.1)
      cbar.ax.tick_params(labelsize=fontlegend)

        # invert y-axis for so that depth level 1 is at top
    plt.gca().invert_yaxis()

        # add title
    ax.set_title(title)

    return fig





#----- BELOW is OLD -----

#===============================================================================================
# Plot map of differences (e.g., anomalies)
#
def plot_map_diff(x, y, z, zmean, title, sunlat, sunlon):

    	# create plot
    print("create plot")
    fig = plt.figure(figsize=(10,7.5))
	# create map
    print("create map")
    proj  = ccrs.PlateCarree()
    axmap = plt.axes(projection=proj)

	# colormap
    ppcmap = "seismic"
    pltmin = -20
    pltmax = 20
    varmax = np.nanmax(z)  
    if varmax!=np.nan and varmax<1:  pltmax=0.5 ; pltmin = -0.5
    if varmax!=np.nan and varmax>=1:  pltmax=1  ; pltmin = -1
    if varmax!=np.nan and varmax>=5:  pltmax=5  ; pltmin = -5
    if varmax!=np.nan and varmax>=10: pltmax=10 ; pltmin = -10
    if varmax!=np.nan and varmax>=15: pltmax=15 ; pltmin = -15
    if varmax!=np.nan and varmax>=20: pltmax=25 ; pltmin = -25
    if varmax!=np.nan and varmax>=30: pltmax=50 ; pltmin = -50
    pltmin = -0.5
    pltmax = 0.5
    #pltmin = -1
    #pltmax = 1
    #pltmin = -5.
    #pltmax = 5. 

	# plot data on map
    print("plot data on map")
    print("shape x = "+str(np.shape(x))+" y = "+str(np.shape(y))+" z = "+str(np.shape(z)))
    pp = plt.pcolormesh(x,y,z,cmap=ppcmap,vmin=pltmin,vmax=pltmax,transform=proj)

	# add ccastlines
    print("add coastlines")
    axmap.coastlines()

	# add grid lines
    print("add gridlines")
    gl = axmap.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.7, linestyle='--')
    gl.xlines = False#True            # plot longitudes
    gl.ylines = False#True            # plot latitudes
        	# set GLOBAL domain limits
    latmin = -90.0
    latmax = 90.0
    lonmin = -180.0
    lonmax = 180.0
    clat   = (latmax + latmin) / 2   
    clon   = (lonmax + lonmin) / 2
    axmap.set_extent([lonmin,lonmax,latmin,latmax], crs=proj)
                	# define x-axis tickmarks
    gl.xlocator = mticker.FixedLocator([-180,-120,-60,0,60,120,180])
                	# define y-axis tickmarks
    gl.ylocator = mticker.FixedLocator([-90,-60,-30,0,30,60,90])

    gl.xlabel_style = {'size': fontaxis, 'color': 'black'}#'gray'}
    gl.ylabel_style = {'size': fontaxis, 'color': 'black'}#'gray'}

    gl.xlabels_top=False
    gl.ylabels_right=False
    #gl.top_labels=None
    #gl.right_labels=None

	# overlay daily mean contours
    if zmean != "none":
	# for anomalies
      #linemin = -50
      #linemax = 50
      #linestride = 10
      #cc = plt.contour(x,y,zmean,colors="navy",levels=range(linemin,linemax,linestride),linewidth=1)
	# for increments
      cc = plt.contour(x,y,zmean,colors="navy",levels=[-0.1,0.1],linewidth=1)

      axmap.clabel(cc,inline=1,fontsize=10,fmt='%1.0f')

	# add location of Sun
    if sunlon != -999 and sunlat != -999:
      axmap.scatter(sunlon, sunlat, c="black", s=150, marker="*", alpha=0.75, label="Sun Location")
      leg = plt.legend(loc='lower right', prop={'size': 8})

	# add colorbar to plot
    print("add colorbar")
    cbarstr = 'SST (deg-C)'

    cbar = fig.colorbar(pp, label=cbarstr, orientation='horizontal', pad=0.05)#0.125)
    cbar.ax.tick_params(labelsize=fontlegend)

	# add title
    print("add title")
    axmap.set_title(title)

    return fig

#===============================================================================================
# Plot map of SD of differences (e.g., anomalies)
#
def plot_map_SDdiff(x, y, z, zmean, title, sunlat, sunlon):

	# create plot
    print("create plot")
    fig = plt.figure(figsize=(10,7.5))
	# create map
    print("create map")
    proj  = ccrs.PlateCarree()
    axmap = plt.subplot(1,1,1,projection=proj)

	# colormap
    ppcmap = "YlOrRd" #"hot"
    pltmin = 1
    pltmax = 20
    varmax = np.nanmax(z)
    if varmax!=np.nan and varmax>=1:  pltmax=1  ; pltmin = 0
    if varmax!=np.nan and varmax>=5:  pltmax=10 ; pltmin = 1
    if varmax!=np.nan and varmax>=10: pltmax=15 ; pltmin = 1
    if varmax!=np.nan and varmax>=20: pltmax=25 ; pltmin = 1
    if varmax!=np.nan and varmax>=30: pltmax=50 ; pltmin = 1
    pltmin = 0
    pltmax = 0.5

	# plot data on map
    print("plot data on map")
    pp = plt.pcolormesh(x,y,z,cmap=ppcmap,vmin=pltmin,vmax=pltmax,transform=proj)

	# add coastlines
    print("add coastlines")
    axmap.coastlines()

	# add grid lines
    print("add gridlines")
    gl = axmap.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.7, linestyle='--')
    gl.xlines = False#True            # plot longitudes
    gl.ylines = False#True            # plot latitudes
		# set GLOBAL domain limits
    latmin = -90.0
    latmax = 90.0
    lonmin = -180.0
    lonmax = 180.0
    clat   = (latmax + latmin) / 2
    clon   = (lonmax + lonmin) / 2
    axmap.set_extent([lonmin,lonmax,latmin,latmax], crs=proj)
                	# define x-axis tickmarks
    gl.xlocator = mticker.FixedLocator([-180,-120,-60,0,60,120,180])
                	# define y-axis tickmarks
    gl.ylocator = mticker.FixedLocator([-90,-60,-30,0,30,60,90])

    gl.xlabel_style = {'size': fontaxis, 'color': 'black'}#'gray'}
    gl.ylabel_style = {'size': fontaxis, 'color': 'black'}#'gray'}

    gl.xlabels_top=False
    gl.ylabels_right=False
    #gl.top_labels=None
    #gl.right_labels=None

	# overlay daily mean contours
    if zmean != "none":
      #cc = plt.contour(x,y,zmean,colors="navy",levels=range(0,50,10),linewidth=1)
      cc = plt.contour(x,y,zmean,colors="navy",levels=[0.1],linewidth=1)
      axmap.clabel(cc,inline=1,fontsize=10,fmt='%1.0f')

	# add location of Sun
    if sunlon != -999 and sunlat != -999:
      axmap.scatter(sunlon, sunlat, c="black", s=150, marker="*", alpha=0.75, label="Sun Location")
      leg = plt.legend(loc='lower right', prop={'size': 8})

	# add colorbar to plot
    print("add colorbar")
    cbarstr = 'SD (deg-C)'

    cbar = fig.colorbar(pp, label=cbarstr, orientation='horizontal', pad=0.05)#0.125)
    cbar.ax.tick_params(labelsize=fontlegend)

	# add title
    print("add title")
    axmap.set_title(title)

    return fig

#===============================================================================================
# Plot histogram
#
def plot_hist(y, y_names, alphaval, colors, units, startdate, enddate, hour, obstype):

    fontaxis = 11
    fontlegend = 10

        # create plot
    fig = plt.figure(figsize=(9,8.5))
    ax1 = plt.subplot()

        # histogram specs
    label = "SST "
    for ilabel in y_names:
      label += ilabel+" "               # axis label, i.e., what's being plotting

    hist_dir = "vertical"
    histmin  = -1
    histmax  = 1
    binsize  = 0.025
    bins     = int((2*int(histmax))/binsize)

    xpos     = histmin
    posdiff  = binsize

    if hist_dir == "vertical":
      ylabel = label+" ("+units+")"
      xlabel = "Obs Count"
    elif hist_dir == "horizontal":
      xlabel = label+" ("+units+")"
      ylabel = "Obs Count"

    trange = (histmin,histmax)
    w      = 1 #0.9                        # fraction of bin width = width of bars

        # count number of obs
    count = np.size(y)

        # get counts per dependent dataset for legend labels
    legendlabel = np.nan * np.ones(np.size(y_names), dtype=list)
    for j in range(np.size(y_names)):
      avg  = np.mean(y[j])                      # mean
      sd   = np.std(y[j])                       # standard deviation
      var  = np.var(y[j])                       # variance
      text = "Mean="+str(np.round_(avg,decimals=2))+" SD="+str(np.round_(sd,decimals=2))#+" RMSD="+str(np.round_(rmsd,decimals=2))
      legendlabel[j] = str(y_names[j])+" | count="+str(np.size(y[j]))+" "+text

      if j==0:
        ymean = avg
        ysd   = np.sqrt(var)
      else:
        ymean = np.append(ymean,avg)
        ysd   = np.append(ysd,np.sqrt(var))

      del avg,sd,var

        # plot
    for j in range(np.size(y_names)):
      print("----- plot_hist: y_name = "+str(y_names[j])+" -----")

        # create histogram from y
        #       counts = bin heights
      counts, bin_edges, _ = plt.hist(y[j], bins, trange, histtype='stepfilled', edgecolor=colors[j], facecolor=colors[j], rwidth=w, alpha=0.2, label=legendlabel[j], orientation=hist_dir, density=False)

      if y_names[j].find("OmB")!=-1 or y_names[j].find("OmA")!=-1:

        bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2
        xx = np.linspace(bin_centers[0], bin_centers[-1], bins)
        mean = np.mean(xx)
        std  = np.std(xx)
        print("plot_hist: mean = "+str(mean)+" std = "+str(std))

        # Gaussian distribution (normal, not heavy-tailed)
        #       p = best fit values for parameters of givel model (gaus)
        p, _ = curve_fit(gaus, bin_centers, counts, p0=[np.max(counts), 0., 1.], maxfev=1000)
        ax1.plot(xx, gaus(xx, *p), c=colors[j], linestyle='dashed')#, label='fit')

        # heavy-tailed distribution (Gaussian-esque)
        #       alpha, beta are shape parameters
        alpha=1.8; beta=-0.5
        ax1.plot(xx, counts+levy_stable.pdf(xx, alpha, beta), c=colors[j])

      del counts,bin_edges

    ax1.set_ylabel(xlabel, fontsize=fontaxis)
    ax1.set_xlabel(ylabel,  fontsize=fontaxis)

    ax1.axvline(x=0,color="black")             # vertical line

        # plot title
    if startdate==enddate:
        plt.title("Histogram of "+str(obstype.upper())+" SST OmB,OmA for "+str(hour)+"z: "+str(startdate), fontsize=14)
    else:
        plt.title("Histogram of "+str(obstype.upper())+" SST OmB,OmA for "+str(hour)+"z: "+str(startdate)+" to "+str(enddate), fontsize=14)

        # add legend to plot
        #       'bbox_to_anchor' = location of legend box. (x,y) = (left=0 and right=1, bottom=0 and top=1)
        #       'center'         = center of bounding box at coords 'bbox_to_anchor'
        #       'centerleft '    = center of left edge of bounding box at coords 'bbox_to_anchor'
        #       'center right'   = center of right edge of bounding box at coords 'bbox_to_anchor'

    legendloc = "upper left"
    leg = plt.legend(loc=legendloc, prop={'size': 12})

    for lh in leg.legendHandles:
      lh.set_alpha(1)

    legendmarkersize = 30
    for lh in range(np.size(leg.legendHandles)):
      leg.legendHandles[lh]._sizes = [legendmarkersize]

    # custom legend for line styles
    ax1b = ax1.twinx()
    ax1b.axes.get_yaxis().set_visible(False)    # hide axis ticks and labels
    custom_markers = [Line2D([0],[0],color="black",label="Levy, \u03B1="+str(alpha),linestyle="solid"), Line2D([0],[0],color="black",label="Gaussian, \u03B1=2",linestyle="dashed")]
    legP1 = ax1b.legend(handles=custom_markers, loc='upper right', prop={'size': fontlegend})
    for lhP1 in range(np.size(legP1.legendHandles)):
      legP1.legendHandles[lhP1]._sizes = [legendmarkersize]

    return fig

#===============================================================================================
# Plot time series
#
def plot_time_series(nDuniq_list,date_uniq,Dyyyymmdd,Dlat,Dlon,Dvert,tx,ty,x_name,tname,match_str,acolors,units,regionstr,latmin,latmax,outname,**kwargs):

        # set font sizes for plotting
    ftsz = 9
    legendmarkersize = ftsz

        # create plot
    fig = plt.figure(figsize=(12,12))

        # create grid of panels
    gs = fig.add_gridspec(nrows=4, ncols=1)

        # create each panel
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0])
    ax3 = fig.add_subplot(gs[2,0])
    ax4 = fig.add_subplot(gs[3,0])

        # panel specs
    if match_str=="HLOS":
      label = "HLOS Wind Velocity"
      labelyaxis = "Wind"
      axismin = -5.0
      axismax = 5.0
      diffpos = 0.2
    elif match_str=="Wind Speed":
      label = match_str
      labelyaxis = "Wind"
      axismin = -10.0
      axismax = 10.0
      diffpos = 5
    elif match_str=="Pressure":
      label = match_str
      labelyaxis = label
      axismin = -15.0
      axismax = 15.0
      diffpos = -5
    elif match_str=="Height":
      label = match_str
      labelyaxis = label
      axismin = -2.0
      axismax = 2.0
      diffpos = 0.5

        # set initial positions of plot text boxes
    if match_str=="Pressure":
      xpos = 0
      ypos = axismin
    elif match_str=="Height":
      xpos = 0
      ypos = axismin
    else:
      xpos = 0
      ypos = axismin

    warnings.filterwarnings('ignore', category=FutureWarning)           # ignore warnings
    datestr  = date_uniq
    ndatestr = len(datestr)
    datelabels = datestr

        # set xaxis
    warnings.filterwarnings('ignore', category=RuntimeWarning)          # ignore warnings
    idx = np.where(Dyyyymmdd==datestr[ndatestr-1])
    if np.size(idx)<=0:
      ilast = ndatestr-1
    else:
      ilast = ndatestr
    del idx
                # x-axis specs
    timestepstr = "days"              # includes 00-23 hours
              # x-axis tickmark values and labels
    if ndatestr<10:
      stride  = 1
    elif ndatestr>=10 and ndatestr<30:
      stride = 5
    elif ndatestr>=30 and ndatestr<93:
      stride = 10
    elif ndatestr>=93:
      stride = 30
    xaxis   = np.arange(0, ilast, 1)
    xvalues = np.arange(0, ilast, stride)
    xlabels = datelabels[0:ilast:stride]

        # OUTPUT STATS TO TEXT FILE
    file1 = open(outname+".txt","a+")            # open file for reading and writing. create file if it does not exist. data is appended to end of file.

        # compute stats per DEPENDENT dataset
    xstat = []
    ystat = []
    nobsmax = 0
    for i in range(np.size(acolors)):
        # variable to plot from DRIVER (x) and DEPENDENT (y)
      txx = tx[i]
      tyy = ty[i]
      ttmp_Dyyyymmdd = Dyyyymmdd[i]
      tlat = Dlat[i]
      idx = np.where((tlat>=latmin)*(tlat<latmax))
      x = txx[idx]
      y = tyy[idx]
      tmp_Dyyyymmdd = ttmp_Dyyyymmdd[idx]
      del txx,tyy,tlat,idx,ttmp_Dyyyymmdd
        # adding text inside the plot
      if match_str=="Pressure":
        xpos = xpos
        ypos = ypos + 50
      elif match_str=="Height":
        xpos = xpos
        ypos = ypos - 1
      else:
        xpos = xpos                                             # value on x-axis where text will begin
        ypos = ypos                                             # value on y-axis where text will begin

        # loop thru all obs for dataset
      diff = []
      corr = []
      sd   = []
      rmsd = []
      nobs = []
      difftot = []
      sdtot   = []
      rmsdtot = []
      nobstot = []
      for j in range(ilast):
        warnings.filterwarnings('ignore', category=RuntimeWarning)              # ignore warnings
                # extract vars based on time
        idx = np.where((x!=np.nan)*(y!=np.nan)*(tmp_Dyyyymmdd==datestr[j]))
        if np.size(idx)>0:
          varxdate = x[idx]
          varydate = y[idx]
          tnobs = np.size(varxdate)
          nobs.append(tnobs)
                # find absolute max of nobs
          if tnobs>nobsmax:
            nobsmax = tnobs
                #compute stats per timestep
          tdiff  = varydate - varxdate                          # y-x = DEPENDENT minus DRIVER
          ttcorr = np.corrcoef(varxdate,varydate)               # correlation
          tcorr  = ttcorr[0,1]
          tavg   = np.mean(tdiff)                               # mean
          tsd    = np.std(tdiff)                                # standard deviation
          trmsd  = np.sqrt(np.mean(tdiff**2))                   # RMSD
                #append stats per timestep
          diff.append(tavg)
          corr.append(tcorr)
          sd.append(tsd)
          rmsd.append(trmsd)
                # stats for total means
          difftot.append(tavg)
          sdtot.append(tsd)
          rmsdtot.append(trmsd)

          del varxdate,varydate,tdiff,ttcorr,tcorr,tavg,tsd,trmsd,tnobs
        else:
          diff.append(np.nan)
          corr.append(np.nan)
          sd.append(np.nan)
          rmsd.append(np.nan)
          nobs.append(np.nan)
        del idx

        # compute stats
      diffarr = np.asarray(diff)
      sdarr   = np.asarray(sd)
      rmsdarr = np.asarray(rmsd)

      sdifftot = np.asarray(difftot)
      ssdtot   = np.asarray(sdtot)
      srmsdtot = np.asarray(rmsdtot)

      xarr     = np.asarray(x)
      yarr     = np.asarray(y)
      tcorr    = np.corrcoef(xarr,yarr)           # correlation
      corr_all = tcorr[0,1]

        # average stats over entire time period
      avg_all  = np.mean(sdifftot)
      sd_all   = np.mean(ssdtot)
      rmsd_all = np.mean(srmsdtot)

      legendlabel1text = str(tname[i])+" | Mean_Diff="+str(np.round_(avg_all,decimals=2))
      legendlabel1 = str(tname[i])#+" | Mean_Diff="+str(np.round_(avg_all,decimals=2))
      del tcorr,xarr,yarr

        # plot 1: time series: difference
      ax1.plot(xaxis, diff, '-', color=acolors[i], linewidth=1, label=legendlabel1)

        # plot 2: correlation coefficient
      legendlabel2text = str(tname[i])+" | r="+str(np.round_(corr_all,decimals=2))
      legendlabel2 = str(tname[i])#+" | r="+str(np.round_(corr_all,decimals=2))
      #ax2.plot(xaxis, corr, color=acolors[i], label=legendlabel2)

        # plot 3: standard deviation, RMSD
                # color lines legend
      legendlabel3text = str(tname[i])+" | SD_Diff="+str(np.round_(sd_all,decimals=2))+" RMSD="+str(np.round_(rmsd_all,decimals=2))
      legendlabel3 = str(tname[i])#+" | SD_Diff="+str(np.round_(sd_all,decimals=2))+" RMSD="+str(np.round_(rmsd_all,decimals=2))

                # black lines legend
      custom_lines = [Line2D([0],[0],color="black",label="SD_Diff",linestyle="-"), Line2D([0],[0],color="black",label="RMSD",linestyle="--")]

      ax3.plot(xaxis, sd, color=acolors[i], label=legendlabel3)
      ax3.plot(xaxis, rmsd, color=acolors[i], linestyle='dashed')

        # plot 4: number of observations
      legendlabel4text = str(tname[i])+" | count="+str(np.size(x))
      legendlabel4 = str(tname[i])#+" | count="+str(np.size(x))
      #ax4.plot(xaxis, nobs, color=acolors[i], label=legendlabel4)              # add legend to bottom plot

      del x,y,sd,rmsd,tmp_Dyyyymmdd

        # OUTPUT STATS TO TEXT FILE
      file1.write(legendlabel1text+"\n")
      file1.write(legendlabel2text+"\n")
      file1.write(legendlabel3text+"\n")
      file1.write(legendlabel4text+"\n")

        # set tick marks for all panels
    ax1.set_xticks(xvalues)
    ax2.set_xticks(xvalues)
    ax3.set_xticks(xvalues)
    ax4.set_xticks(xvalues)
        # set tick mark labels for bottom panel only
    ax1.set_xticklabels([])
    ax2.set_xticklabels([])
    ax3.set_xticklabels([])
    ndatestr = 5
    ax4.set_xticklabels(xlabels, rotation=45)

        # add DRIVER obs count to plot 1 legend
    Dlegendlabel = "Total Unique DRIVER Count = "+str(nDuniq_list)
    file1.write(Dlegendlabel+"\n")

        # axis limits
    nobsarr = np.asarray(nobs)
    print("time_series: nobsarr max = "+str(np.max(nobsarr)))

    ax1.set_xlim([min(xaxis),max(xaxis)])
    ax1.set_ylim([axismin,axismax])
    ax2.set_xlim([min(xaxis),max(xaxis)])
    ax2.set_ylim([0,1])
    ax3.set_xlim([min(xaxis),max(xaxis)])
    ax3.set_ylim([0,axismax*2])

    ax4.set_xlim([min(xaxis),max(xaxis)])
    y4values = [10,100,1000,10000,25000]
    y4labels = ["10","100","1000","10000","25000"]
    ax4.set_yticks(y4values)
    ax4.set_yticklabels(y4labels)#, fontsize=fontaxis+2)
    ax4.set_yscale('log')

        # plot 1 labels
    ax1.axhline(y=0,color="black")                              # horizontal line
    ax1.axvline(x=0,color="black")                              # vertical line
    if x_name.find("Aeolus")!=-1: x_name = "Aeolus"
    ax1.set_title("Daily Mean Time Series: "+regionstr+" "+label+" Difference [DEPENDENT (colors) - "+x_name+"]", fontsize=fonttitle-3)
    ax1.set_ylabel(labelyaxis+" Diff. ("+units+")", fontsize=10)

        # legend
    leg1 = ax1.legend(loc='lower right', prop={'size': ftsz})
    for lh in leg1.legendHandles:
      lh.set_alpha(1)
    for lh in range(np.size(leg1.legendHandles)):
      leg1.legendHandles[lh]._sizes = [legendmarkersize]

        # plot 2 labels
    ax2.axhline(y=0,color="black")                              # horizontal line
    ax2.axvline(x=0,color="black")                              # vertical line
    ax2.set_ylabel("Corr. Coeff. (r)", fontsize=10)

                # legend
    #leg2 = ax2.legend(loc='lower right', prop={'size': ftsz})
    #for lh in leg2.legendHandles:
    #  lh.set_alpha(1)
    #for lh in range(np.size(leg2.legendHandles)):
    #  leg2.legendHandles[lh]._sizes = [legendmarkersize]

        # plot 3 labels
    ax3.axhline(y=0,color="black")                              # horizontal line
    ax3.axvline(x=0,color="black")                              # vertical line
    ax3.set_ylabel("RMSD, SD of Diff ("+units+")", fontsize=10)

                # color lines legend
    #leg3a = ax3.legend(loc='lower right', prop={'size': ftsz})
    #for lh in leg3a.legendHandles:
    #  lh.set_alpha(1)
    #for lh in range(np.size(leg3a.legendHandles)):
    #  leg3a.legendHandles[lh]._sizes = [legendmarkersize]

                # black lines legend
    ax3b = ax3.twinx()
    leg3b = ax3b.legend(handles=custom_lines, loc='lower left', prop={'size': ftsz})
    for lh in leg3b.legendHandles:
      lh.set_alpha(1)
    for lh in range(np.size(leg3b.legendHandles)):
      leg3b.legendHandles[lh]._sizes = [legendmarkersize]
    ax3b.tick_params(right=False,labelright=False)

        # plot 4 labels
    ax4.axhline(y=0,color="black")                              # horizontal line
    ax4.axvline(x=0,color="black")                              # vertical line
    ax4.set_ylabel("Count", fontsize=10)

                # legend
    #leg4 = ax4.legend(loc='lower right', prop={'size': ftsz})
    #for lh in leg4.legendHandles:
    #  lh.set_alpha(1)
    #for lh in range(np.size(leg4.legendHandles)):
    #  leg4.legendHandles[lh]._sizes = [legendmarkersize]

                # DRIVER count legend
    #ax4D = ax4.twinx()
    #Dlegend = ax4D.plot([],[],' ',label=Dlegendlabel)
    #leg4D = ax4D.legend(handles=Dlegend, loc='upper left', prop={'size': ftsz})
    #for lh in leg4D.legendHandles:
    #  lh.set_alpha(1)
    #for lh in range(np.size(leg4D.legendHandles)):
    #  leg4D.legendHandles[lh]._sizes = [legendmarkersize]
    #ax4D.tick_params(right=False,labelright=False)

    del diff,nobs,diffarr,sdarr

    file1.close()       # CLOSE TEXT FILE

    return fig

#===============================================================================================
# Plot 2D cross-section
#
#	Optional Arguments: pltmin, pltmax
#		... uses default unless specified in function call
#
def plot_cross_section(x, y, z, zmean, xvals, yvals, xlabel, ylabel, title, plottypeflag, pltmin=-999, pltmax=999):

    	# create plot
    print("create plot")
    fig = plt.figure(figsize=(10,7.5))

    ax = fig.add_subplot()


        # set axes labels
    ax.set_xlabel(xlabel, fontsize=fontaxis)
    ax.set_ylabel(ylabel, fontsize=fontaxis)
            # x-axis
    nspacing = int(30 * (np.size(xvals)/360))
    print("plot_cross: nspacing = "+str(nspacing)+" size xvals = "+str(np.size(xvals)))
    txlabels = xvals[::nspacing]
    del nspacing
                    # print minus sign (not hyphen) in front of negative numbers on axis
    nspacing = int(30 * (np.size(x)/360))
    xvalues = np.arange(0,np.size(x),nspacing)
    xlabels = [[] for i in range(np.size(x))]
    print("size x = "+str(np.size(x))+" xvals = "+str(np.size(xvals)))
    i=0
    j=0
    while (i < np.size(x)):
      xlabels[i] = str(txlabels[j])
      i += nspacing
      j += 1
    ax.set_xticks(xvalues)
    ax.set_xticklabels(xlabels[::nspacing], fontsize=fontaxis-2)
    del xvalues,txlabels,xlabels,nspacing
            # y-axis
    nspacing = 5
    yvalues = np.arange(0,np.size(yvals),nspacing)
    ylabels = yvals[::nspacing]
    ax.set_yticks(yvalues)
    ax.set_yticklabels(ylabels, fontsize=fontaxis-2)
    del yvalues,ylabels,nspacing

	# set colormap and color range values
    varmin = np.nanmin(z)
    varmax = np.nanmax(z)
    if plottypeflag.find("diff")!=-1 or plottypeflag.find("Diff")!=-1 or plottypeflag.find("Sign")!=-1 or plottypeflag.find("sign")!=-1:
      ppcmap = "seismic"
      if pltmin==-999:
        if varmax!=np.nan and varmax<1:   pltmax=0.5
        if varmax!=np.nan and varmax>=1:  pltmax=1 
        if varmax!=np.nan and varmax>=5:  pltmax=5 
        if varmax!=np.nan and varmax>=10: pltmax=10
        if varmax!=np.nan and varmax>=15: pltmax=15
        if varmax!=np.nan and varmax>=20: pltmax=25
        if varmax!=np.nan and varmax>=30: pltmax=50
        pltmin = -pltmax
    elif plottypeflag.find("SD")!=-1:
      ppcmap = "YlOrRd" #"hot"
      if pltmin==-999:
        if varmax!=np.nan and varmax>=1:  pltmax=1  ; pltmin = 0
        if varmax!=np.nan and varmax>=5:  pltmax=10 ; pltmin = 1
        if varmax!=np.nan and varmax>=10: pltmax=15 ; pltmin = 1
        if varmax!=np.nan and varmax>=20: pltmax=25 ; pltmin = 1
        if varmax!=np.nan and varmax>=30: pltmax=50 ; pltmin = 1
    else:
      #if plottypeflag.find("bkgerr")!=-1:
      #  ppcmap = "jet"
      #else:
      ppcmap = "nipy_spectral"
      if pltmin==-999 and varmin<0:
        if varmax!=np.nan and varmax<1:   pltmax=0.5
        if varmax!=np.nan and varmax>=1:  pltmax=1 
        if varmax!=np.nan and varmax>=5:  pltmax=5 
        if varmax!=np.nan and varmax>=10: pltmax=10
        if varmax!=np.nan and varmax>=15: pltmax=15
        if varmax!=np.nan and varmax>=20: pltmax=25
        if varmax!=np.nan and varmax>=30: pltmax=50
        pltmin = -pltmax
      elif pltmin==-999 and varmin>=0:
        if varmax!=np.nan and varmax>=1:  pltmax=1  ; pltmin = 0
        if varmax!=np.nan and varmax>=5:  pltmax=10 ; pltmin = 1
        if varmax!=np.nan and varmax>=10: pltmax=15 ; pltmin = 1
        if varmax!=np.nan and varmax>=20: pltmax=25 ; pltmin = 1
        if varmax!=np.nan and varmax>=30: pltmax=50 ; pltmin = 1

	# plot data on map
    originstr = "lower"
    #pp = plt.imshow(z, cmap=ppcmap, vmin=pltmin, vmax=pltmax, origin=originstr, aspect='auto')
    pp = plt.imshow(z, cmap=ppcmap, vmin=pltmin, vmax=pltmax, origin=originstr, aspect='auto')
    #pp = plt.pcolormesh(x, y, z, cmap=ppcmap, vmin=pltmin, vmax=pltmax, transform=proj)

    #ax.set_xlabel(xlabel)
    #ax.set_xticks(x)
    ##ax.set_xticklabels(x, fontsize=8)
    #ax.set_xticklabels(xlabels, fontsize=8)
    #ax.set_ylabel(ylabel)
    #ax.set_yticks(y)
    ##ax.set_yticklabels(y, fontsize=8)
    #ax.set_yticklabels(ylabels, fontsize=8)
 
	# overlay daily mean contours
#    if zmean != "none":
#	# for anomalies
#      #linemin = -50
#      #linemax = 50
#      #linestride = 10
#      #cc = plt.contour(x,y,zmean,colors="navy",levels=range(linemin,linemax,linestride),linewidth=1)
#	# for increments
#      cc = plt.contour(x,y,zmean,colors="navy",levels=[-0.1,0.1],linewidth=1)
#      ax.clabel(cc,inline=1,fontsize=10,fmt='%1.0f')

	# add colorbar to plot
    if plottypeflag.find("Sign")!=-1 or plottypeflag.find("sign")!=-1:		# IS sign change plot
      redlabel = "OmA>0 and OmB<0"
      bluelabel = "OmA<0 and OmB>0"
      custom_markers = [Line2D([0],[0],color="red",label=redlabel,linestyle="None",marker="s"), Line2D([0],[0],color="blue",label=bluelabel,linestyle="None",marker="s")]
      legP1 = ax.legend(handles=custom_markers, loc='lower center', bbox_to_anchor=(0.5, 0.), prop={'size': fontlegend})
      for lhP1 in legP1.legendHandles:
        lhP1.set_alpha(1)
      legendmarkersize = 30
      for lhP1 in range(np.size(legP1.legendHandles)):
        legP1.legendHandles[lhP1]._sizes = [legendmarkersize]

    elif plottypeflag.find("Sign")==-1 or plottypeflag.find("sign")==-1:		# IS NOT sign change plot
      #cbarstr = 'SST'
      #cbar = fig.colorbar(pp, label=cbarstr, orientation='horizontal', pad=0.05)#0.125)
      cbar = fig.colorbar(pp, orientation='horizontal', pad=0.1)
      cbar.ax.tick_params(labelsize=fontlegend)

	# invert y-axis for so that depth level 1 is at top
    plt.gca().invert_yaxis()

	# add title
    ax.set_title(title)

    return fig

# -------------------------------------------------------------------------
# Density Scatter Plot
#
def density_scatter(tx, ty, x_name, y_name, outname, plotvar, title, var, units, regionstr, pltmin=-999, pltmax=999, ax=None,**kwargs):

        # create plot
    fig = plt.figure(figsize=(10,8.5))

    if ax is None:
      ax = plt.subplot()

    x = tx
    y = ty

        # get colormap 'cmap'
    ppcmap = plt.cm.get_cmap("jet")

        # set plot specs
    if var=="sst":
      label = regionstr+" SST"
      axismin = pltmin
      axismax = pltmax
      txpos   = pltmin+1
      typos   = pltmax
      if pltmax<=1:
        txpos = -0.1
        diffpos = 0.1
      else:
        diffpos = 5

        # conform x,y to 2D array
    if pltmax<=1:
      gridcellsize = 0.05
    else:
      gridcellsize = 0.25 #1
    xaxis = list(np.arange(axismin,axismax,gridcellsize))
    yaxis = xaxis
    var2d = var_to_2d_counts_1dset(x_name,y_name,xaxis,yaxis,x,y)

        # set range of plotted colored contours
    warnings.filterwarnings('ignore', category=RuntimeWarning)          # ignore warnings
    varmax = np.nanmax(var2d)
    varmin = np.nanmin(var2d)
        # reassign color range max
    pltmin = 1
    pltmax = 100
    if varmax!=np.nan and varmax>=100000  : pltmax=1000000
    if varmax!=np.nan and varmax<100000: pltmax=10000
    if varmax!=np.nan and varmax<10000: pltmax=1000
    if varmax!=np.nan and varmax<1000: pltmax=100
    if varmax!=np.nan and varmax<100: pltmax=10
    var2d  = np.where(var2d==np.nan,0,var2d)
    del varmax,varmin

        # plot
    originstr = 'lower'
    if pltmax==10:
      pp = plt.imshow(var2d,cmap=ppcmap,vmin=pltmin,vmax=pltmax, extent=[axismin,axismax,axismin,axismax], origin=originstr, aspect='auto')
    else:
      pp = plt.imshow(var2d,cmap=ppcmap,norm=LogNorm(vmin=pltmin,vmax=pltmax), extent=[axismin,axismax,axismin,axismax], origin=originstr, aspect='auto')

        # statistical significance
    signif = 95.0                                      # statistical signifiance in percent (%). Type: float

    tsigmax = (100.0-signif)/100.0                     # max p-value allowed for statistical signiificance
    tsig = ttest_ind(x,y).pvalue                       # get p-value from Student's t-test
    print("dens_scat: tsig = "+str(tsig))
    if tsig<=tsigmax:                                  # if p-value < tsigmax (difference is statistically significant at 95% level)
      statsigstr = "Diffs are signif. at "+str(int(signif))+"%"

        # set axis limits
    ax.set_xlim([axismin,axismax])
    ax.set_ylim([axismin,axismax])
    ax.tick_params(labelsize=fontaxis+2)

        # draw specific lines on plot
    ax.axhline(y=0,color="black")             # horizontal line
    ax.axvline(x=0,color="black")             # vertical line
    #ax.axline((0,0),slope=1.0,color="black")  # one-to-one line
            # one-to-one line
    xx = np.linspace(-100, 100, 500)
    yy = xx
    plt.plot(xx, yy, color='black')
    del xx, yy

        # add text inside the plot
    xpos = txpos                                      # value on x-axis where text will begin
    ypos = typos                                      # value on y-axis where text will begin
    ftsz = 12                                         # font size of text_* (see below)

        # ADD TEXT INSIDE PANEL
        # add obs count to dataset name for legend
    xlegendlabel = str(x_name)+" count = "+str(np.size(x))
    ylegendlabel = str(y_name)+" count = "+str(np.size(y))
    ypos = ypos - diffpos
    #plt.text(xpos, ypos, legendlabel, fontsize = ftsz)

        # compute and print stats of differences
    diff  = y - x
    tcorr = np.corrcoef(x,y)               # correlation
    corr  = tcorr[0,1]
    avg   = np.mean(diff)                  # mean
    sd    = np.std(diff)                   # standard deviation
    rmsd  = np.sqrt(np.mean(diff**2))      # RMSD

                # ADD TEXT INSIDE PANEL
                # add title for difference stats
    text = "----- Difference Stats ("+str(y_name)+" - "+str(x_name)+") -----"
    #ypos = ypos - diffpos
    #plt.text(xpos, ypos, text, fontsize = ftsz)
                # add correlation
    textr = "r = "+str(corr) #str(np.round_(corr,decimals=3))
    ypos = ypos - diffpos/2
    plt.text(xpos, ypos, textr, fontsize = ftsz)
                # add mean diff
    textm = "Mean_Diff = "+str(np.round_(avg,decimals=3))
    ypos = ypos - diffpos/2
    plt.text(xpos, ypos, textm, fontsize = ftsz)
                # add SD of diff
    texts = "SD_Diff = "+str(np.round_(sd,decimals=3))
    ypos = ypos - diffpos/2
    plt.text(xpos, ypos, texts, fontsize = ftsz)
                # add RMSD
    #texts = "RMSD = "+str(np.round_(rmsd,decimals=2))
    #ypos = ypos - diffpos/2
    #plt.text(xpos, ypos, texts, fontsize = ftsz)
    if tsig<=tsigmax:
                # add space
      texts = " "
      #ypos = ypos - diffpos/2
      #plt.text(xpos,ypos,texts,fontsize=ftsz)
                # state whether differences are statistcially significant
      texts = statsigstr
      ypos = ypos - diffpos/2
      plt.text(xpos,ypos,texts,fontsize=ftsz)
      #del tsig,tsigmax,signif
    #del diff,tcorr,corr,avg,sd,rmsd,text,textr,textm,texts

                # OUTPUT STATS TO TEXT FILE
    file1 = open(outname+".txt","a+")            # open file for reading and writing. create file if it does not exist. data is appended to end of file.
    file1.write(text+"\n")
    file1.write(xlegendlabel+"\n")
    file1.write(ylegendlabel+"\n")
    file1.write(textr+"\n")
    file1.write(textm+"\n")
    file1.write(texts+"\n")
    textrmsd = "RMSD = "+str(np.round_(rmsd,decimals=2))
    file1.write(textrmsd+"\n")
    if tsig<=tsigmax: file1.write(statsigstr+"\n")
    file1.close()
    del diff,tcorr,corr,avg,sd,rmsd,text,textr,textm,texts
    del tsig,tsigmax,signif

      # plot title and axis labels
    ax.set_title(title)#, fontsize=fonttitle)
    ax.set_xlabel(var.upper()+" "+x_name+" ("+str(units)+")", fontsize=fontaxis+2)
    ax.set_ylabel(var.upper()+" "+y_name+" ("+str(units)+")", fontsize=fontaxis+2)

        # Add a colorbar to the bottom of the plot.
    cbarstr = "Number Density"
    frac = 0.046
    cbar = fig.colorbar(pp, label=cbarstr, orientation='vertical', fraction=frac, pad=0.04)
    #cbar = fig.colorbar(pp, orientation='vertical', fraction=frac, pad=0.04)
    cbar.ax.tick_params(labelsize=fontlegend+2)
    #cbar.ax.set_label(cbarstr)#, size="medium")

    return fig

# -------------------------------------------------------------------------
# Map of Locations of Matched (Collocated) Observations
#
def plot_map_points2d_ce(varstr, tx, ty, ss, xaxis, yaxis, title, marksize, alphaval, datein, units, opt, plottypeflag, pltmin=-999, pltmax=999, **kwargs):

        # create plot
    fig = plt.figure(figsize=(10,7.5))

        # create map
    proj = ccrs.PlateCarree()
    #ax = plt.subplot(1,1,1,projection=proj)
    ax = plt.axes(projection=proj)

    ax.coastlines()

    gl = ax.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.7, linestyle='--')
    gl.xlines = False            # plot longitudes
    gl.ylines = False            # plot latitudes
    #gl.xlines = True
    #gl.ylines = True

    #gl = ax.gridlines(draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
    #gl.xlines = True            # plot longitudes
    #gl.ylines = True            # plot latitudes
    #gl.top_labels=None
    #gl.right_labels=None

    region_flag=0
    if region_flag==0:
        # set GLOBAL domain limits
      latmin = -90.0
      latmax = 90.0
      lonmin = -180.0
      lonmax = 180.0
      clat   = (latmax + latmin) / 2
      clon   = (lonmax + lonmin) / 2
      ax.set_extent([lonmin,lonmax,latmin,latmax], crs=proj)
        # define x-axis tickmarks
      gl.xlocator = mticker.FixedLocator([-180,-120,-60,0,60,120,180])
        # define x-axis tickmarks
      gl.ylocator = mticker.FixedLocator([-90,-60,-30,0,30,60,90])

    #gl.xlabel_style = {'size': 15, 'color': 'gray'}
    #gl.ylabel_style = {'size': 15, 'color': 'gray'}

    #gl.top_labels=None
    #gl.right_labels=None

    gl.xlabel_style = {'size': fontaxis, 'color': 'black'}#'gray'}
    gl.ylabel_style = {'size': fontaxis, 'color': 'black'}#'gray'}

    gl.xlabels_top=False
    gl.ylabels_right=False
    #gl.top_labels=None
    #gl.right_labels=None

    warnings.filterwarnings('ignore', category=RuntimeWarning)          # ignore warnings

        # set colormap and color range values
    varmin = np.nanmin(ss)
    varmax = np.nanmax(abs(ss))
    if plottypeflag.find("diff")!=-1 or plottypeflag.find("Diff")!=-1 or plottypeflag.find("Sign")!=-1 or plottypeflag.find("sign")!=-1:
      ppcmap = "seismic"
      if pltmin==-999:
        if varmax!=np.nan and varmax<1:   pltmax=0.5
        if varmax!=np.nan and varmax>=1:  pltmax=1
        if varmax!=np.nan and varmax>=5:  pltmax=5
        if varmax!=np.nan and varmax>=10: pltmax=10
        if varmax!=np.nan and varmax>=15: pltmax=15
        if varmax!=np.nan and varmax>=20: pltmax=25
        if varmax!=np.nan and varmax>=30: pltmax=50
        pltmin = -pltmax
    elif plottypeflag.find("SD")!=-1:
      ppcmap = "YlOrRd" #"hot"
      if pltmin==-999:
        if varmax!=np.nan and varmax>=1:  pltmax=1  ; pltmin = 0
        if varmax!=np.nan and varmax>=5:  pltmax=10 ; pltmin = 1
        if varmax!=np.nan and varmax>=10: pltmax=15 ; pltmin = 1
        if varmax!=np.nan and varmax>=20: pltmax=25 ; pltmin = 1
        if varmax!=np.nan and varmax>=30: pltmax=50 ; pltmin = 1
    else:
      #ppcmap = "jet"
      ppcmap = "nipy_spectral"
      if pltmin==-999 and varmin<0:
        if varmax!=np.nan and varmax<1:   pltmax=0.5
        if varmax!=np.nan and varmax>=1:  pltmax=1
        if varmax!=np.nan and varmax>=5:  pltmax=5
        if varmax!=np.nan and varmax>=10: pltmax=10
        if varmax!=np.nan and varmax>=15: pltmax=15
        if varmax!=np.nan and varmax>=20: pltmax=25
        if varmax!=np.nan and varmax>=30: pltmax=50
        pltmin = -pltmax
      elif pltmin==-999 and varmin>=0:
        if varmax!=np.nan and varmax>=1:  pltmax=1  ; pltmin = 0
        if varmax!=np.nan and varmax>=5:  pltmax=10 ; pltmin = 1
        if varmax!=np.nan and varmax>=10: pltmax=15 ; pltmin = 1
        if varmax!=np.nan and varmax>=20: pltmax=25 ; pltmin = 1
        if varmax!=np.nan and varmax>=30: pltmax=50 ; pltmin = 1
    del varmax,varmin

        # add obs counts to dataset name for legend
    #Dlegendlabel = str(x_name)

        # plot DRIVER points
    Dlons  = tx
    Dlats  = ty

        # plot DRIVER on map
    #sc = ax.scatter(Dlons, Dlats, c=ss, cmap=ppcmap, vmin=pltmin, vmax=pltmax, s=marksize, marker="o", label=Dlegendlabel, edgecolors='none')
    sc = ax.scatter(Dlons, Dlats, c=ss, cmap=ppcmap, vmin=pltmin, vmax=pltmax, s=marksize, marker="o", edgecolors='none')

    plt.title(title, fontsize=9)

        # Add a colorbar to the bottom of the plot.
    cbar = plt.colorbar(sc, orientation='horizontal', pad=0.1)
    cbar.set_label(str(units))

    return fig








#===============================================================================================
