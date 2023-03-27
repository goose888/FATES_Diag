# remap_harvest_rate.r

# this script reads the GCAM Agro-eco region map, LUH2 (or data from other 
# sources with the same data format) wood C harvest rate and ELM-FATES (or 
# alternative vegetation dynamic model output) forest carbon inventory to remap 
# LUH2 biomass harvest rates following the biomass density map
# from FATES calculation and conserve the total harvested biomass.

# required variables to set below (lines after the grid mapping function):
#	luh2_file: luh2 transitions file - source of harvest data (in kg C per year) (1/4 degree res)
#	grid_file: elm griddata file - source of new grid; needs to match the land_file grid
#	land_file: elm landuse timeseries file - the original land surface file to copy and modify
#	out_land_file: new elm landuse file - biomass harvest data are now in the original variables
#	start_year_in: historical file first record is 850, future file fist record is 2015
#					note that luh2 harvest data do not exist for the final year in the file
#                   e.g., there are no 2015 harvest data in the historical files that spans 850-2015

# output:
#	new landuse timeseries file with biomass harvest in kg C per year
#	note that this is harvest removed, after 30% loss of estimated total carbon affected
#	the original fraction data are moved to different variable names
#	the new biomass harvest data are added using the original variable names
#	this is because the model code looks specifically for the original variable names
#	existing out_land_file will be overwritten

# notes:
#	1) include the full directory path with each file, so that they can be wherever
#	2) using the griddata file ensures the correct mapping the lat-lon info in 
# the land use file is not sufficient to determine the grid
# 3) current version does not include the remap based on gcam regions
# only diagnosis based on gcam region is available option 
#	4) currently the harvest data for year are stored in year+1 in the land surface file
#	so duplicate the first year data to store in the first year slot 
# (even though it is not used). This should be changed in e3sm phase 3 restructuring
# 5) also: the diagnostic time series data have not been shifted and so are 
# plotted accordingly
# 6) This script maybe highly inefficient handling high-res input data. I only 
# tested the inputs with a global 4x5 degree resolution

library(raster)
library(rasterVis)
library(RColorBrewer)
library(terra)
library(ncdf4)
library(parallel)
library(plyr)

################# calculate max mode function [non area-weighted]
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

################# adjust search area if touching boundary
# note: This logic is not suitable along latitude. Think if we touch the 
# northern boundary we shall search for grids to another hemisphere. However,
# in this case the search will definetely outside the continent for sure.
# Thus right now I ignore this case and the search will shift from northern to
# southern hemisphere, which likely be bareground and not producing weird 
# cross-continental remapping of the harvest rate.
adjust_bound <- function(gid, shift, bound) {
  if(gid+shift <= 0) {
    shift = bound + shift
  }
  if(gid+shift > bound) {
    shift = - bound + shift
  }
  return(shift)
}

################# remap the harvest rate based on certain rules
# This function require the following inputs:
#  -> diff map of forest C density between total harvest demand and available 
#     forest C as base map
#  -> Forest C density from basemap
#  -> Accumulated harvest rate over 1850 - 2010
# 1) Gridcells with negative value from the diff map contain harvest debt.
# 2) For each gridcell with debt, we define a radius starts from 1 and search 
#    the gridcells falling in this search area.
# 3) Calculate the total available C in this area. If the total available C is 
#    higher than harvest rate, we distribute the harvest debt based on the relative
#    fraction of forest C density with the search area. If not, continue step 2)
#    by expanding radius until the debt is clear.
# 4) Update the harvest rate time series.
#
remap_harv_rate <- function(nonforest, nlon_out, nlat_out, nyr, forest_c, harv_c, harv_ts)	{
  # Arguments:
  # nonforest - logical argument to decide if non-forest harvest rate is
  #             accounted in the input
  # nlon_out - length of longitude
  # nlat_out - length of latitude
  # nyr - length of year
  # forest_c - Forest C density (kgC m-2 grid area)
  # harv_c - accumulated harvested C density (kgC m-2 grid area)
  # harv_ts - Time series of wood harvest rate (kgC per grid)

  surplus_map = (forest_c - harv_c) *1e6*area_grid
  surplus_map[which(surplus_map <= 0)] <- NA
  harv_debt = (forest_c - harv_c) *1e6*area_grid
  harv_debt[which(harv_debt >= 0)] <- NA
  # Longitude and latitude loop
  for (ln in 1:nlon_out){
    for (lt in 1:nlat_out){
      # If current grid has harvest debt
      cur_harv_debt = harv_debt[ln,lt] 
      if(!is.na(cur_harv_debt) && cur_harv_debt < 0.0){
        # Reset output
        redis_factor = array(0, dim = c(nlon_out,nlat_out))
        # Set search radius
        i = 0
        while (cur_harv_debt < 0.0) {
          hotspot_for_c = 0.0
          i = i + 1
          # Sum the available forest C within the range
          for (t in -i:i) {
            for (n in -i:i) {
              n = adjust_bound(ln, n, nlon_out)
              t = adjust_bound(lt, t, nlat_out)
              # Set NA value to 0.0
              if(!is.na(surplus_map[ln+n,lt+t])){
                hotspot_for_c = hotspot_for_c + surplus_map[ln+n,lt+t]
              }
            }
          }
          cur_harv_debt = harv_debt[ln,lt] + hotspot_for_c
        } # End of radius search
        # Calculate the adjustment factor of harvest rate to redistribute 
        # harvest debt over the region based on the proportion of
        # surplus forest C
        redis_factor[ln,lt] = harv_debt[ln,lt] / (harv_c[ln,lt]*1e6*area_grid[ln,lt])
        for (t in -i:i) {
          for (n in -i:i) {
            n = adjust_bound(ln, n, nlon_out)
            t = adjust_bound(lt, t, nlat_out)
            if(!is.na(surplus_map[ln+n,lt+t])){
              redis_factor[ln+n,lt+t] = - redis_factor[ln,lt] * surplus_map[ln+n,lt+t] / hotspot_for_c
              # Update surplus map
              surplus_map[ln+n,lt+t] = surplus_map[ln+n,lt+t] - redis_factor[ln+n,lt+t] * (harv_c[ln,lt]*1e6*area_grid[ln,lt])
            }
          }
        }
        # Update harvest debt map
        harv_debt[ln,lt] = 0.0
        # Adjust the Time Series of harvest rate
        for (yr in 1:nyr){
          harv_ts[1,,,yr] = harv_ts[1,,,yr] + harv_ts[1,ln,lt,yr] * redis_factor[,]
          if(nonforest){
            harv_ts[2,,,yr] = harv_ts[2,,,yr] + harv_ts[2,ln,lt,yr] * redis_factor[,]
          }
        }
      } # End of harvest debt check
    }
  } # End of cell loop

  return(harv_ts)
  
} # end of harvest rate remap function

################# remap the harvest rate for a specific year
# A duplication of "remap_harv_rate" for secondary forest.
# There are two major difference in this function
# 1) The calculated revision factors will only be applied to the harvest rate in current year
# 2) After the revision of harvest rate, forest density will be updated to exclude
# the harvested amount in all the following years after the current year
remap_harv_rate_per_yr <- function(nonforest, nlon_out, nlat_out, nyear_out, cyr, forest_c, harv_c, harv_c_allcat, area_grid, out_list)	{
  # Arguments:
  # nonforest - logical argument to determine if the non-forest harvest rate is 
  # accounted for. Note the logic here is different from the remap function of
  # primary forest
  # nlon_out - length of longitude
  # nlat_out - length of latitude
  # cyr - current year
  # forest_c - forest C density of current year
  # harv_c - harvested C density of current year
  # harv_c_allcat - harvested C density of current year from all categories
  # area_grid - grid area
  # Return:
  # out_list - a list contains two outputs:
  # 1) harv_ts - wood harvest rate of the current year
  # 2) new_forest_c - forest C density time series after subtracting the 
  # revised harvested C 
  # 
  
  # Initialization
  small_tol = 0.0
  
  surplus_map = (forest_c - harv_c) * 1e6 * area_grid
  surplus_map[which(surplus_map <= 0)] <- NA
  harv_debt = (forest_c - harv_c) * 1e6 * area_grid
  harv_debt[which(harv_debt >= 0)] <- NA
  harv_ts <- out_list[[1]]
  
  # Longitude and latitude loop
  for (ln in 1:nlon_out){
    for (lt in 1:nlat_out){
      # If current grid has harvest debt
      cur_harv_debt = harv_debt[ln,lt] 
      if(!is.na(cur_harv_debt) && cur_harv_debt < 0.0){
        # Reset output
        redis_factor = array(0, dim = c(nlon_out,nlat_out))
        # Set search radius
        i = 0
        while (cur_harv_debt < 0.0) {
          hotspot_for_c = 0.0
          i = i + 1
          if(i>30){
              cat("Too large radius for point: ", ln, lt, "\n")
              cat("Current harvest debt?", cur_harv_debt, "\n")
              
              cat("Current year?", cyr, "\n")
              cat("non-forest?", nonforest, "\n")
              readline(prompt="Pause, press [enter] to continue")
          }
          # Sum the available forest C within the range
          for (t in -i:i) {
            for (n in -i:i) {
              n = adjust_bound(ln, n, nlon_out)
              t = adjust_bound(lt, t, nlat_out)
              # Exclude NA value
              if(!is.na(surplus_map[ln+n,lt+t]) && surplus_map[ln+n,lt+t] > small_tol){
                hotspot_for_c = hotspot_for_c + surplus_map[ln+n,lt+t]
              }
            }
          }
          if(i>30){
            cat("hotspot_for_c?", hotspot_for_c, "\n")
            readline(prompt="Pause, press [enter] to continue")
          }
          cur_harv_debt = harv_debt[ln,lt] + hotspot_for_c
        } # End of search

        # Calculate the adjustment factor of harvest rate to redistribute 
        # harvest debt over the region based on the proportion of
        # surplus forest C
        redis_factor[ln,lt] = harv_debt[ln,lt] / (harv_c[ln,lt] * 1e6 * area_grid[ln,lt])
        if(nonforest){
          harv_ts[4,ln,lt] = harv_ts[4,ln,lt] + (redis_factor[ln,lt] * harv_c_allcat[4,ln,lt])*1e6*area_grid[ln,lt]
          harv_ts[5,ln,lt] = harv_ts[5,ln,lt] + (redis_factor[ln,lt] * harv_c_allcat[5,ln,lt])*1e6*area_grid[ln,lt]
        } else {
          harv_ts[3,ln,lt] = harv_ts[3,ln,lt] + (redis_factor[ln,lt] * harv_c[ln,lt])*1e6*area_grid[ln,lt]
        }
        for (t in -i:i) {
          for (n in -i:i) {
            n = adjust_bound(ln, n, nlon_out)
            t = adjust_bound(lt, t, nlat_out)
            if(!is.na(surplus_map[ln+n,lt+t])){
              redis_factor[ln+n,lt+t] = - redis_factor[ln,lt] * surplus_map[ln+n,lt+t] / hotspot_for_c
              # if((ln+n) == 26 && (lt+t) == 34){
              #   cat("Year?", cyr, " \n")
              #   cat("Forest C of the point.", forest_c[ln+n,lt+t], " \n")
              #   cat("Point found before.", surplus_map[ln+n,lt+t], " \n")
              # }
              
              # Update surplus map
              surplus_map[ln+n,lt+t] = surplus_map[ln+n,lt+t] - redis_factor[ln+n,lt+t] * harv_c[ln,lt] * 1e6 * area_grid[ln,lt]
              # Adjust the harvest rate of current year
              if(nonforest){
                harv_ts[4,ln+n,lt+t] = harv_ts[4,ln+n,lt+t] + (redis_factor[ln+n,lt+t] * harv_c_allcat[4,ln,lt])*1e6*area_grid[ln,lt]
                harv_ts[5,ln+n,lt+t] = harv_ts[5,ln+n,lt+t] + (redis_factor[ln+n,lt+t] * harv_c_allcat[5,ln,lt])*1e6*area_grid[ln,lt]
              } else {
                harv_ts[3,ln+n,lt+t] = harv_ts[3,ln+n,lt+t] + (redis_factor[ln+n,lt+t] * harv_c[ln,lt])*1e6*area_grid[ln,lt]
              }
              # if((ln+n) == 26 && (lt+t) == 34){
              #   cat("Point found after.", surplus_map[ln+n,lt+t], " \n")
              #   cat("Redis factor.", redis_factor[ln+n,lt+t] * harv_c[ln,lt], " \n")
              #   cat("Source harv c.", harv_c[ln,lt], " \n")
              #   cat("The grid id, ", ln+n, lt+t, "\n")
              #   cat("The source grid id, ", ln, lt, "\n")
              #   cat("The total available C for source is, ", hotspot_for_c, "\n")
              #   readline(prompt="Pause, press [enter] to continue")
              # }
            }
          }
        }
        # Update harvest debt map
        harv_debt[ln,lt] = 0.0
        # Adjust the harvest rate of current year
        # harv_ts[,] = harv_ts[,] + harv_ts[ln,lt]*redis_factor[,]
        out_list[[1]] <- harv_ts
      } # End of harvest debt check
    }
  } # End of cell loop
  
  # Update harvest rate
  out_list[[1]] <- harv_ts
  # Adjust the Time Series of forest C
  new_forest_c <- out_list[[2]]
  for (i in cyr:nyear_out) {
    if(nonforest){
      new_forest_c[,,i] = new_forest_c[,,i] - (harv_ts[4,,]+harv_ts[5,,])*1e-6 / area_grid
    } else {
      new_forest_c[,,i] = new_forest_c[,,i] - harv_ts[3,,]*1e-6 / area_grid
    }
  }

  # tmp_test = harv_ts*1e-6 / area_grid
  # cat("harv_ts?", tmp_test[2,38], "\n")
  # cat("new_forest_c?", new_forest_c[2,38,1:10], "\n")
  # readline(prompt="Press [enter] to continue")
  
  # Check if any NAN exists in harv_ts
  # tmp_array = is.nan(harv_ts)
  # if(length(harv_ts[tmp_array]) != 0){
  #   readline(prompt="Press [enter] to continue")
  #   cat("In which year we have NAN?", cyr, "\n")
  #   cat("How many NANs do we have?", length(harv_ts[tmp_array]), "\n")
  #   tmp_out = new_forest_c[,,cyr]
  #   cat("New_forest_c:", forest_c[tmp_array], "\n")
  # }
  # Do a negative check on new_forest_c
  for (ln in 1:nlon_out){
    for (lt in 1:nlat_out){
      if(!is.na(new_forest_c[ln,lt,cyr]) && new_forest_c[ln,lt,cyr] < -1e-16){
        cat("In which year we have negative?", cyr, "\n")
        cat("Location?", ln, lt, "\n")
        cat("New_forest_c:", new_forest_c[ln,lt,cyr], "\n")
        cat("Harv_c:", harv_c[ln,lt], "\n")
        readline(prompt="Press [enter] to continue")
      }
    }
  }
  # Tiny negative forest C due to less vegc in the future year is possible, thus force 
  # those negative values to zero and those gridcells will not contribute in future years
  new_forest_c[new_forest_c<0] = 0.0
  out_list[[2]] <- new_forest_c
  
  return(out_list)
  
} # end of harvest rate remap function

################# Upscale categorical raster data by finding the max-mode
# Since the input here is categorical map, we pick the max-mode instead of 
# performing numerical average
upscale_cell <- function(cell, lonW, lonE, latS, latN, p4string, raster_stack)	{
  # get the current out cell bounds
  # this assumes that cell starts at 0 to eliminate a calc inside
  lt = trunc(cell / nlon_out) + 1
  ln = cell %% nlon_out + 1
  lnmin = lonW[ln,lt]
  lnmax = lonE[ln,lt]
  ltmin = latS[lt,lt]
  ltmax = latN[lt,lt]
  
  # need to shift the longitude values so that they range from -180 to 180 rather than -xmin to xmax
  # deal with the lon boundary straddle also
  # lat values are already within -90 to 90
  if (lnmin < 180 & lnmax > 180) {
    # need two spatial objects
    # need each corner of each cell, clockwise, and first and last are the same
    lnmax_tmp = lnmax
    lnmax = 180
    corners = cbind(c(lnmin, lnmax, lnmax, lnmin, lnmin), c(ltmax, ltmax, ltmin, ltmin, ltmax))
    lnmin = -180
    lnmax =  -180 + (lnmax_tmp - 180)
    corners2 = cbind(c(lnmin, lnmax, lnmax, lnmin, lnmin), c(ltmax, ltmax, ltmin, ltmin, ltmax))
    p1 = Polygons( list(Polygon(coords = corners, hole=FALSE)) , 1)
    p2 = Polygons( list(Polygon(coords = corners2, hole=FALSE)) , 2)
    spobj = SpatialPolygons( list(p1, p2) , proj4string=CRS(p4string))
  } else { # only one spatial object needed
    if (lnmin >= 180 & lnmax > 180) {
      lnmin = -180 + (lnmin - 180)
      lnmax = -180 + (lnmax - 180)
    }
    # need each conrer of the cell, clockwise, and first and last are the same
    corners = cbind(c(lnmin, lnmax, lnmax, lnmin, lnmin), c(ltmax, ltmax, ltmin, ltmin, ltmax))
    spobj = SpatialPolygons( list(Polygons( list(Polygon(coords = corners, hole=FALSE)) , 1)) , proj4string=CRS(PROJ4_STRING))
  } # end if-else for lon adjustment
  
  # Get the max-mode
  # put the values in a single row of a data frame with the cell index for id
  harv_cell_df = extract(raster_stack, spobj, weights=TRUE, normalizeWeights=FALSE, cellnumbers=TRUE, df=TRUE, small=TRUE, exact=TRUE, na.rm=TRUE)
  out_df = harv_cell_df[1,c("cell", "layer")]
  out_df$cell = cell
  out_df[1, "layer"] = calculate_mode(harv_cell_df[,"layer"])

  return(out_df)
  
} # end single-cell mapping function


############ change this to reflect correct directory structure and files ! ###################
gcam_reg_file = "C:/Users/sshu3/anaconda_wkspace/FATES/AEZ_orig.grd"
grid_file = "C:/Users/sshu3/anaconda_wkspace/FATES/griddata_4x5_060404.nc"
land_file = "C:/Users/sshu3/anaconda_wkspace/FATES/landuse.timeseries_4x5_HIST_simyr1700-2015.biomass_harvest.nc"
# land_file = "C:/Users/sshu3/anaconda_wkspace/FATES/landuse.timeseries_4x5_hist_simyr1850-2015_200311_biomass_harvest.nc"
biom_den_file = "C:/Users/sshu3/anaconda_wkspace/FATES/fates.forest.vegc.avg_1700_1707.nc"
secb_den_file = "C:/Users/sshu3/anaconda_wkspace/FATES/fates.sec.forest.vegc.ts_1700_2015.nc"
out_land_file = "C:/Users/sshu3/anaconda_wkspace/FATES/landuse.timeseries_4x5_hist_harmonized_simyr1700-2015.biomass_harvest.nc"

# Options
# Stage 1 --- Harmonize primary forest
# Stage 2 --- Harmonize secondary forest
stage = 2
# Set true to include non-forest harvest rate
inc_non_forest = TRUE
start_year_in = 1800
stop_year_in = 2015
# Validate the regional total based on GCAM regions
gcam_constraint = FALSE
# Check the diff map between harvest rate and harvestable forest C
diag_plots = FALSE
num_cores = detectCores()
# kilograms to Petagrams C
k2P = 1/1000000000000

# Global management parameters to select part of the harvestable forest C
# percentage of trunk in the total harvestable forest C
pct_trunk = 0.67
# Direct logging ratio, usually 1.0 for global run
logging_direct = 1.0
# The percentage that will be exported from on-site to local product pool
# usually 1.0 for global run
logging_export = 1.0

cat("Start remap_harvest_rate at", date(), "\n")

# 1. Read in LUH2 dataset and calculate the required accumulated 
#    harvested C density [kgC m-2, grid-level], 1850 - 2010 

# the number of harvest variables and their names
# the names are matched by order for the in and out files
# primary forest, primary non-forest, secondary mature forest, secondary young forest, secondary non-forest
num_harvest = 5
luh2_bioharv_names = c("primf_bioh", "primn_bioh", "secmf_bioh", "secyf_bioh", "secnf_bioh")
ls_harv_names = c("HARVEST_VH1", "HARVEST_VH2", "HARVEST_SH1", "HARVEST_SH2", "HARVEST_SH3")
ls_harv_names_frac = c("HARVEST_VH1_LUH2", "HARVEST_VH2_LUH2", "HARVEST_SH1_LUH2", "HARVEST_SH2_LUH2", "HARVEST_SH3_LUH2")

# open the land surface file and get years present and lat-lon
lsid = nc_open(land_file)
nlon_out = lsid$dim$lsmlon$len
nlat_out = lsid$dim$lsmlat$len
# Test
nyears_out = lsid$dim$time$len
years_out = ncvar_get(lsid, varid = "YEAR", start=c(1), count=c(nyears_out))
latixy_out = ncvar_get(lsid, varid = "LATIXY", start=c(1,1), count=c(nlon_out,nlat_out))
longxy_out = ncvar_get(lsid, varid = "LONGXY", start=c(1,1), count=c(nlon_out,nlat_out))
ncells_out = nlon_out * nlat_out
harv_rate_in = array(dim=c(num_harvest, nlon_out, nlat_out, nyears_out))
for (h in 1:num_harvest) {
  harv_rate_in[h,,,] = ncvar_get(lsid, varid = ls_harv_names[h], start=c(1,1,1), count=c(nlon_out,nlat_out,nyears_out))
}
nc_close(lsid)

# open the new grid file and get the new cell boundaries and lat-lon
gdid = nc_open(grid_file)
nlon_grid = gdid$dim$lsmlon$len
nlat_grid = gdid$dim$lsmlat$len
lonE = ncvar_get(gdid,varid="LONE",start=c(1,1), count=c(nlon_grid, nlat_grid))
lonW = ncvar_get(gdid,varid="LONW",start=c(1,1), count=c(nlon_grid, nlat_grid))
latN = ncvar_get(gdid,varid="LATN",start=c(1,1), count=c(nlon_grid, nlat_grid))
latS = ncvar_get(gdid,varid="LATS",start=c(1,1), count=c(nlon_grid, nlat_grid))
latixy_grid = ncvar_get(gdid, varid = "LATIXY", start=c(1,1), count=c(nlon_grid,nlat_grid))
longxy_grid = ncvar_get(gdid, varid = "LONGXY", start=c(1,1), count=c(nlon_grid,nlat_grid))
area_grid = ncvar_get(gdid, varid = "AREA", start=c(1,1), count=c(nlon_grid,nlat_grid))
nc_close(gdid)

# Sum harvested C from 1700 to 2015 and calculate the density
harv_rate_acc = array(dim=c(num_harvest, nlon_out, nlat_out))
harv_rate_acc[,,]= 0
harv_rate_den = array(dim=c(num_harvest, nlon_out, nlat_out))
harv_rate_den[,,]= 0
harv_rate_den_ts = array(dim=c(num_harvest, nlon_out, nlat_out, nyears_out))
harv_rate_den_ts[,,,]= 0
for (h in 1:num_harvest) {
  harv_rate_acc[h,,] = rowSums(harv_rate_in[h,,,], dims = 2)
  harv_rate_den[h,,] = harv_rate_acc[h,,]*1e-6 / area_grid
  for (i in 1:nyears_out) {
    harv_rate_den_ts[h,,,i] = harv_rate_in[h,,,i]*1e-6 / area_grid
  }
}

if(gcam_constraint){
   # Read in GCAM agro-ecological zone (AEZ) file
   gcam_reg <- brick(gcam_reg_file)
   hdr(gcam_reg, format = 'ENVI')
   PROJ4_STRING = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
   # PROJ4_STRING=crs(gcam_reg)

   # Upscale to 4 by 5 degs
   cat("\nstarting cell-by-cell mapping for GCAM regional mask \n")
   mcout = mclapply(c(0:(ncells_out-1)), function(i) upscale_cell(cell = i, lonW, lonE, latS, latN, PROJ4_STRING, gcam_reg), mc.cores = 1)
   cat("finishing cell-by-cell mapping. \n")

   reg2elm_out = rbind.fill(mcout)
   reg2elm_out = luh2elm_out[order(luh2elm_out$cell),]
   gcam_out = array(dim=c(nlon_out, nlat_out))
   gcam_out[,] = reg2elm_out[, "layer"]
}
# 2. Read in biomass density from ELM-FATES
if(stage == 1) {
   bioid = nc_open(biom_den_file)
   # Potential forest is calculated through FATES spinup under pre-industrial 
   # CO2 level
   forestc_pri_den = ncvar_get(bioid, varid = "FATES_PRI_VEGC_FOREST", start=c(1,1), count=c(nlon_out,nlat_out))
   forestc_pri_den = forestc_pri_den*pct_trunk*logging_direct*logging_export
   nc_close(bioid)
} else if(stage == 2) {
   # Potential secondary forest is calculated through FATES via
   # assuming no secondary forest harvest but do apply primary forest harvest
   # with the pre-industrial level (or historical?) CO2
   bioid = nc_open(secb_den_file)
   nyears_cal = stop_year_in - start_year_in
   forestc_sec_den = ncvar_get(bioid, varid = "FATES_SEC_MATURE_VEGC_FOREST", start=c(1,1,1), count=c(nlon_out,nlat_out,nyears_out))
   forestc_sec_den = forestc_sec_den*pct_trunk*logging_direct*logging_export
   # Young secondary forest biomass is all the available biomass 
   # with age less than 94 years
   # Mature secondary forest contains all the available biomass after 
   # subtracting young forest from total secondary forest
   forestc_secy_den = ncvar_get(bioid, varid = "FATES_SEC_YOUNG_VEGC_FOREST", start=c(1,1,1), count=c(nlon_out,nlat_out,nyears_out))
   forestc_secy_den = forestc_secy_den*pct_trunk*logging_direct*logging_export
   nc_close(bioid)
}

# Quick diagnosis of the harvest rate diff
if(diag_plots) {
  diff_den = harv_rate_den[1,,] - forestc_pri_den
  pdf(file=paste0(outpath, "check_harvrate_luh2_minus_fatesvegc.pdf"), width=8.5, height=4)
  cols <- brewer.pal(3, "RdYlGn")
  pal <- colorRampPalette(cols)
  print(levelplot(diff_den, margin=FALSE, at=seq(-10, 10, 1), col.regions = pal, colorkey = list(title = list("[LUH2 - FATES, Primary]", cex = 1, fontface = "bold", col = "red"))))
  dev.off()
}

# 3. Harmonize LUH2 harvest rate with FATES forest C distribution
if(stage == 1) {
  if(inc_non_forest){
     tmp_harv_rate_den = harv_rate_den[1,,] + harv_rate_den[2,,] 
     new_harv_ts = remap_harv_rate(inc_non_forest, nlon_out, nlat_out, nyears_out, forestc_pri_den, tmp_harv_rate_den, harv_rate_in[,,,])
  } else {
     new_harv_ts = remap_harv_rate(inc_non_forest, nlon_out, nlat_out, nyears_out, forestc_pri_den, harv_rate_den[1,,], harv_rate_in[,,,])
  }
} else if(stage == 2) {
  if(inc_non_forest){
    new_harv_ts = array(dim=c(num_harvest, nlon_out, nlat_out, nyears_out))
    new_harv_ts[,,,]= 0
    tmp_ls <- vector(mode='list',length=2)
    # Secondary mature
    for (i in (nyears_out-nyears_cal+1):nyears_out){
      tmp_ls[[1]] <- harv_rate_in[,,,i]
      tmp_ls[[2]] <- forestc_sec_den[,,]
      tmp_ls = remap_harv_rate_per_yr(FALSE, nlon_out, nlat_out, nyears_out, i, forestc_sec_den[,,i], harv_rate_den_ts[3,,,i], harv_rate_den_ts[,,,i], area_grid, tmp_ls)
      # update forestc_sec_den
      new_harv_ts[,,,i] <- tmp_ls[[1]]
      forestc_sec_den <- tmp_ls[[2]]
    }
    # Store the harmonized mature secondary forest harvest rate in temporary container
    tmp_harv_ts = new_harv_ts[3,,,]
    # Secondary young
    for (i in (nyears_out-nyears_cal+1):nyears_out){
      tmp_ls[[1]] <- harv_rate_in[,,,i]
      tmp_ls[[2]] <- forestc_secy_den[,,]
      tmp_ls = remap_harv_rate_per_yr(TRUE, nlon_out, nlat_out, nyears_out, i, forestc_secy_den[,,i], (harv_rate_den_ts[4,,,i]+harv_rate_den_ts[5,,,i]), harv_rate_den_ts[,,,i], area_grid, tmp_ls)
      # update forestc_sec_den
      new_harv_ts[,,,i] <- tmp_ls[[1]]
      forestc_secy_den <- tmp_ls[[2]]
    }
    new_harv_ts[3,,,(nyears_out-nyears_cal+1):nyears_out] = tmp_harv_ts[,,(nyears_out-nyears_cal+1):nyears_out]
    # Secondary forest harvest rate before the calculation year shall be 0
    new_harv_ts[3:5,,,1:(nyears_out-nyears_cal)] = 0.
  } else {
    new_harv_ts = array(dim=c(nlon_out, nlat_out, nyears_out))
    new_harv_ts[,,]= 0
    tmp_ls <- vector(mode='list',length=2)
    for (i in 1:nyears_out){
      tmp_ls[[1]] <- harv_rate_in[3,,,i]
      tmp_ls[[2]] <- forestc_sec_den[,,]
      tmp_ls = remap_harv_rate_per_yr(FALSE, nlon_out, nlat_out, nyears_out, i, forestc_sec_den[,,i], harv_rate_den_ts[3,,,i], harv_rate_den_ts[,,,i], area_grid, tmp_ls)
      # update forestc_sec_den
      new_harv_ts[,,i] <- tmp_ls[[1]]
      forestc_sec_den <- tmp_ls[[2]]
    }
  }
}

# 4. Validate the new harvest rate time series vs. before
harvest_in_globe_series = array(dim=c(nyears_out))
harvest_in_globe_series[] = 0
harvest_out_globe_series = array(dim=c(nyears_out))
harvest_out_globe_series[] = 0
if(stage == 1) {
   for (yind_out in 1:(nyears_out-1)) {
      # match the in year to the out year
      harvest_in_globe_series[yind_out] = sum(harv_rate_in[1,,,yind_out], na.rm=TRUE) * k2P
      harvest_out_globe_series[yind_out] = sum(new_harv_ts[1,,,yind_out], na.rm=TRUE) * k2P
   }
   plot(years_out[1:(nyears_out-1)], harvest_in_globe_series[1:(nyears_out-1)], type = "l", lty = 1, main = "Annual total global biomass harvest",
        xlab = "Year", ylab = "Harvested carbon in Pg C")
   lines(years_out[1:(nyears_out-1)], harvest_out_globe_series[1:(nyears_out-1)], lty=2, col = "red")
   if(inc_non_forest){
     for (yind_out in 1:(nyears_out-1)) {
       # match the in year to the out year
       harvest_in_globe_series[yind_out] = sum(harv_rate_in[2,,,yind_out], na.rm=TRUE) * k2P
       harvest_out_globe_series[yind_out] = sum(new_harv_ts[2,,,yind_out], na.rm=TRUE) * k2P
     }
     lines(years_out[1:(nyears_out-1)], harvest_in_globe_series[1:(nyears_out-1)], lty=1, col = "blue")
     lines(years_out[1:(nyears_out-1)], harvest_out_globe_series[1:(nyears_out-1)], lty=2, col = "green")
     legend("topleft", legend = c("LUH2 Primary Forest Harvest", "Harmonized Primary Forest Harvest", "LUH2 Primary Non-forest Harvest", "Harmonized Primary Non-forest Harvest"), lty=c(1,2,1,2), col = c("black", "red", "blue", "green"), cex = 0.8)
   } else {
     legend("topleft", legend = c("LUH2 Primary Forest Harvest", "Harmonized Primary Forest Harvest"), lty=c(1,2), col = c("black", "red"), cex = 0.8)
   }
} else if(stage == 2){
   for (yind_out in 1:(nyears_out-1)) {
      # match the in year to the out year
      harvest_in_globe_series[yind_out] = sum(harv_rate_in[3,,,yind_out], na.rm=TRUE) * k2P
      harvest_out_globe_series[yind_out] = sum(new_harv_ts[3,,,yind_out], na.rm=TRUE) * k2P
   }
   plot(years_out[(nyears_out-nyears_cal+1):(nyears_out-1)], harvest_in_globe_series[(nyears_out-nyears_cal+1):(nyears_out-1)], type = "l", lty = 1, main = "Annual total global biomass harvest",
        xlab = "Year", ylab = "Harvested carbon in Pg C")
   lines(years_out[(nyears_out-nyears_cal+1):(nyears_out-1)], harvest_out_globe_series[(nyears_out-nyears_cal+1):(nyears_out-1)], lty=2, col = "red")
   if(inc_non_forest){
     for (yind_out in 1:(nyears_out-1)) {
       # match the in year to the out year
       harvest_in_globe_series[yind_out] = sum(harv_rate_in[4,,,yind_out], na.rm=TRUE) * k2P
       harvest_out_globe_series[yind_out] = sum(new_harv_ts[4,,,yind_out], na.rm=TRUE) * k2P
     }
     lines(years_out[(nyears_out-nyears_cal+1):(nyears_out-1)], harvest_in_globe_series[(nyears_out-nyears_cal+1):(nyears_out-1)], lty=1, col = "blue")
     lines(years_out[(nyears_out-nyears_cal+1):(nyears_out-1)], harvest_out_globe_series[(nyears_out-nyears_cal+1):(nyears_out-1)], lty=2, col = "green")
     for (yind_out in 1:(nyears_out-1)) {
       # match the in year to the out year
       harvest_in_globe_series[yind_out] = sum(harv_rate_in[5,,,yind_out], na.rm=TRUE) * k2P
       harvest_out_globe_series[yind_out] = sum(new_harv_ts[5,,,yind_out], na.rm=TRUE) * k2P
     }
     lines(years_out[(nyears_out-nyears_cal+1):(nyears_out-1)], harvest_in_globe_series[(nyears_out-nyears_cal+1):(nyears_out-1)], lty=1, col = "gold")
     lines(years_out[(nyears_out-nyears_cal+1):(nyears_out-1)], harvest_out_globe_series[(nyears_out-nyears_cal+1):(nyears_out-1)], lty=2, col = "cyan")
     legend("topleft", legend = c("LUH2 Secondary Mature Forest Harvest", "Harmonized Secondary Mature Forest Harvest", 
                                  "LUH2 Secondary Young Forest Harvest", "Harmonized Secondary Young Forest Harvest",
                                  "LUH2 Secondary Non-forest Harvest", "Harmonized Secondary Non-forest Harvest"), 
            lty=c(1,2,1,2,1,2), col = c("black", "red", "blue", "green", "gold", "cyan"), cex = 0.8)
   } else {
     legend("topleft", legend = c("LUH2 Secondary Forest Harvest", "Harmonized Secondary Forest Harvest"), lty=c(1,2), col = c("black", "red"), cex = 0.8)
   }
}

 # 5. Overwrite the wood harvest TS
 if(stage == 1){
   # clear all secondary forest
   new_harv_ts[3:5,,,] = 0.
 }

 # this overwrites an existing out_land_file
 file.copy(land_file, out_land_file, overwrite=TRUE)

 cat("Before writing new file\n")

 # rename original fraction data and add biomass harvest data with original names
 lsid = nc_open(out_land_file, write=TRUE)

 # get the dimension objects
 lon_dim = lsid$dim[['lsmlon']]
 lat_dim = lsid$dim[['lsmlat']]
 time_dim = lsid$dim[['time']]

 for (h in 1:num_harvest) {
   # rename fraction data
   lsid <- ncvar_rename(lsid, ls_harv_names[h], ls_harv_names_frac[h])

   # define and add new variables; missing value is 0
   new_var <- ncvar_def(ls_harv_names[h], "kg C per year", list(lon_dim, lat_dim, time_dim), prec="double")
   lsid <- ncvar_add(lsid, new_var)

   # add biomass harvest data
   ncvar_put(lsid, varid = ls_harv_names[h], vals = new_harv_ts[h,,,], start = c(1,1,1), count = c(nlon_out, nlat_out, nyears_out))
 } # end for loop over h for updating output file

 nc_close(lsid)

 cat("Finish remap_harvest_rate at", date(), "\n")

# # 6. Post-processing and generate plots
# LVP<-vector(mode='list',length=2)
# 
# r1 <-c(rast(new_harv_ts[,,1]), rast(new_harv_ts[,,80]), rast(new_harv_ts[,,160]))
# r2 <-c(rast(harv_rate_in[1,,,1]), rast(harv_rate_in[1,,,80]), rast(harv_rate_in[1,,,160]))
# 
# crs(r1) <- PROJ4_STRING
# crs(r2) <- PROJ4_STRING
# 
# LVP[[1]]<- levelplot(r1, layout=c(1,r1@ptr$nlyr()))
# LVP[[2]]<- levelplot(r2, layout=c(1,r1@ptr$nlyr()))
# 
# do.call(grid.arrange,c(LVP,ncol=2))


