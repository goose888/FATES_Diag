#!/usr/bin/env python

"""
This package contains function for plotting FATES global output.
Figures are for the manuscript: 
"""
#%%

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib import ticker, cm, colors
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import cartopy.mpl.ticker as cticker
from cartopy.util import add_cyclic_point

def fast_plot_case(lonnc, latnc, plotvar, set_levels, set_ticks, save_path, savefig = False):
    """
    Make contour plot for only one case
    # lonnc - longitude of the variable
    # latnc - latitude of the variable
    # plotvar - 2-D array of data for plot
    # set_levels - array for setting colorbar levels
    # set_ticks - array for setting colorbar ticks
    # save_path - text, path the figure will save to
    # savefig - turn on or off to save figure
    This function does not have a return value
    """
    fig = plt.figure(figsize=(10,10), dpi=300)
    ax = plt.axes(projection=ccrs.PlateCarree())
    cs = plt.contourf(lonnc, latnc, plotvar[:,:], levels=set_levels, \
                transform=ccrs.PlateCarree(), cmap=plt.cm.jet, extend='both')
    ax.coastlines()

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.4, color='gray', alpha=0.5, linestyle='--')
    gl.xformatter = LONGITUDE_FORMATTER

    gl.yformatter = LATITUDE_FORMATTER
    cbar = fig.colorbar(cs, fraction=0.1, orientation='horizontal', ticks=set_ticks )
    plt.show()
    if(savefig):
        plt.savefig(save_path)   

def plot_all_cases(modname, sel_yr, lonnc, latnc, varlist, set_levels, set_ticks, set_unit_tag, save_path, savefig = False):
    """
    Make contour plot for all cases
    # modname - list contains name of all model cases
    # sel_yr - selected frame of the variable
    # lonnc - longitude of the variable
    # latnc - latitude of the variable
    # varlist - list contains the data of corresponding variables
    # set_levels - array for setting colorbar levels
    # set_ticks - array for setting colorbar ticks
    # set_unit_tag - text, Unit of the variable
    # save_path - text, path the figure will save to
    # savefig - turn on or off to save figure
    """
    # Plot 4 x 3
    # Define the figure and each axis for the 4 rows and 3 columns
    fig, axs = plt.subplots(nrows=4,ncols=3,
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        figsize=(24,16))

    # axs is a 2 dimensional array of `GeoAxes`.  We will flatten it into a 1-D array
    axs=axs.flatten()

    # Subplot of all cases 
    for i in np.arange(0, len(modname)):
        # Contour plot
        cs=axs[i].contourf(lonnc, latnc, varlist[i][:,:], levels=set_levels, \
                       transform=ccrs.PlateCarree(), cmap=plt.cm.jet, extend='both')

        # Title each subplot with the name of the model
        axs[i].set_title(modname[i], fontsize=32)

        # Draw the coastines for each subplot
        axs[i].coastlines()
    
        gl = axs[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=0.4, color='gray', alpha=0.5, linestyle='--')
        gl.xlocator = ticker.FixedLocator([-180, -120, -60, 0, 60, 120, 180])
        gl.ylocator = ticker.FixedLocator([-60, -30, 0, 30, 60, 90])
    
        # Define the xticks for longitude
        axs[i].set_xticks(np.arange(-180,181,60), crs=ccrs.PlateCarree())
        lon_formatter = cticker.LongitudeFormatter()
        axs[i].xaxis.set_major_formatter(lon_formatter)

        # Define the yticks for latitude
        axs[i].set_yticks(np.arange(-60,91,30), crs=ccrs.PlateCarree())
        lat_formatter = cticker.LatitudeFormatter()
        axs[i].yaxis.set_major_formatter(lat_formatter)

    # Add a colorbar axis at the bottom of the graph
    cbar_ax = fig.add_axes([0.2, 0.1, 0.6, 0.02])

    # Draw the colorbar
    cbar=fig.colorbar(cs, cax=cbar_ax, orientation='horizontal', spacing='proportional', 
                  extend='both', ticks = set_ticks)
    cbar.ax.tick_params(labelsize=20)
    fig.text(0.05, 0.1, set_unit_tag, fontsize=28)
    # plt.tight_layout()
    if(savefig):
        plt.savefig(save_path, dpi=300)
    else:
        plt.show()

def plot_logging_cases(modname, sel_yr, lonnc, latnc, varlist, set_levels, set_ticks, set_unit_tag, save_path, savefig = False, inc_elm = False):
    """
    Make contour plot for all logging cases
    Selected variable for plot shall include the results from all cases with logging on
    # modname - list contains name of all model cases
    # sel_yr - selected frame of the variable
    # lonnc - longitude of the variable
    # latnc - latitude of the variable
    # varlist - list contains the data of corresponding variables
    # set_levels - array for setting colorbar levels
    # set_ticks - array for setting colorbar ticks
    # set_unit_tag - text, Unit of the variable
    # save_path - text, path the figure will save to
    # savefig - turn on or off to save figure
    # inc_elm - if the figure contains ELM case output or not
    """
    # Plot 2 x 4
    # Define the figure and each axis for the 4 rows and 3 columns
    fig, axs = plt.subplots(nrows=4,ncols=2,
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        figsize=(16,20))

    # axs is a 2 dimensional array of `GeoAxes`.  We will flatten it into a 1-D array
    axs=axs.flatten()

    # Subplot of all cases 
    for i in np.arange(0, len(modname)):
        # Contour plot
        cs=axs[i].contourf(lonnc, latnc, varlist[i][:,:], levels=set_levels, \
                       transform=ccrs.PlateCarree(), cmap=plt.cm.jet, extend='both')

        # Title each subplot with the name of the model
        axs[i].set_title(modname[i], fontsize=32)

        # Draw the coastines for each subplot
        axs[i].coastlines()
    
        gl = axs[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=0.4, color='gray', alpha=0.5, linestyle='--')
        gl.xlocator = ticker.FixedLocator([-180, -120, -60, 0, 60, 120, 180])
        gl.ylocator = ticker.FixedLocator([-60, -30, 0, 30, 60, 90])
    
        # Define the xticks for longitude
        axs[i].set_xticks(np.arange(-180,181,60), crs=ccrs.PlateCarree())
        lon_formatter = cticker.LongitudeFormatter()
        axs[i].xaxis.set_major_formatter(lon_formatter)

        # Define the yticks for latitude
        axs[i].set_yticks(np.arange(-60,91,30), crs=ccrs.PlateCarree())
        lat_formatter = cticker.LatitudeFormatter()
        axs[i].yaxis.set_major_formatter(lat_formatter)
        
    if(inc_elm):
        cs=axs[len(modname)].contourf(lonnc, latnc, varlist[len(modname)][:,:], levels=set_levels, \
                       transform=ccrs.PlateCarree(), cmap=plt.cm.jet, extend='both')

        # Title each subplot with the name of the model
        axs[len(modname)].set_title('elm_a', fontsize=32)

        # Draw the coastines for each subplot
        axs[len(modname)].coastlines()
    
        gl = axs[len(modname)].gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=0.4, color='gray', alpha=0.5, linestyle='--')
        gl.xlocator = ticker.FixedLocator([-180, -120, -60, 0, 60, 120, 180])
        gl.ylocator = ticker.FixedLocator([-60, -30, 0, 30, 60, 90])
    
        # Define the xticks for longitude
        axs[len(modname)].set_xticks(np.arange(-180,181,60), crs=ccrs.PlateCarree())
        lon_formatter = cticker.LongitudeFormatter()
        axs[len(modname)].xaxis.set_major_formatter(lon_formatter)

        # Define the yticks for latitude
        axs[len(modname)].set_yticks(np.arange(-60,91,30), crs=ccrs.PlateCarree())
        lat_formatter = cticker.LatitudeFormatter()
        axs[len(modname)].yaxis.set_major_formatter(lat_formatter)

    # Add a colorbar axis at the bottom of the graph
    cbar_ax = fig.add_axes([0.2, 0.1, 0.6, 0.02])

    # Draw the colorbar
    cbar=fig.colorbar(cs, cax=cbar_ax, orientation='horizontal', spacing='proportional', 
                  extend='both', ticks = set_ticks)
#     cbar.set_ticklabels(['-1%', '-0.75%', '-0.5%', '-0.25%', '0', '0.25%', '0.5%', '0.75%', '1%'])
    cbar.ax.tick_params(labelsize=20)
    fig.text(0.05, 0.1, set_unit_tag, fontsize=28)
    # plt.tight_layout()
    if(savefig):
        plt.savefig(save_path, dpi=300)
    else:
        plt.show()
    return [fig, axs]
