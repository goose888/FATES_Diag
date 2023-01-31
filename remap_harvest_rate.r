# remap_harvest_rate.r

# this script reads the GCAM Agro-eco region map, LUH2 biomass harvest and ELM-FATES 
# forest carbon to remap LUH2 biomass harvest rates following the biomass density map
# from FATES calculation and conserve the total harvested biomass for each region.
#

# required variables to set below (lines after the grid mapping function):
#	luh2_file: luh2 transitions file - source of harvest data (in kg C per year) (1/4 degree res)
#	grid_file: elm griddata file - source of new grid; needs to match the land_file grid
#	land_file: elm landuse timeseries file - the original land surface file to copy and modify
#	out_land_file: new elm landuse file - biomass harvest data are now in the original variables
#	start_year_in: historical file first record is 850, future file fist record is 2015
#					note that luh2 harvest data do not exist for the final year in the file
#                   e.g., there are no 2015 harvest data in the historical files that spans 850-2015

# note:
#	include the full directory path with each file, so that they can be wherever

# note:
#	using the griddata file ensures the correct mapping
#	the lat-lon info in the land use file is not sufficient to determine the grid

# note:
#	this has been tested only for regridding to coarser resolution (> 1/4 degree)
#   it should work going to finer resolution (if i understand the extract function), but maybe less efficiently

# output:
#	new landuse timeseries file with biomass harvest in kg C per year
#	note that this is harvest removed, after 30% loss of estimated total carbon affected
#	the original fraction data are moved to different variable names
#	the new biomass harvest data are added using the original variable names
#	this is because the model code looks specifically for the original variable names
#	existing out_land_file will be overwritten

# note:
#	currently the harvest data for year are stored in year+1 in the land surface file
#	so duplicate the first year data to store in the first year slot (even though it is not used)
#	this should be changed in e3sm phase 3 restructuring
# also: the diagnostic time series data have not been shifted and so are plotted accordingly

library(raster)
library(rasterVis)
library(RColorBrewer)
library(ncdf4)
library(parallel)
library(plyr)

##### the single cell function for parallelizing must be defined first !!!!!!!

################# calculate max mode function
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

################# adjust search area if touching boundary
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
remap_harv_rate <- function(nlon_out, nlat_out, nyr, forest_c, harv_c, harv_ts)	{
  # Arguments:
  # nlon_out - length of longitude
  # nlat_out - length of latitude
  # nyr - length of year
  # forest_c - Forest C density
  # harv_c - accumulated harvested C density
  # harv_ts - Time series of wood harvest rate

  surplus_map = forest_c - harv_c
  surplus_map[which(surplus_map <= 0)] <- NA
  harv_debt = forest_c - harv_c
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
        redis_factor[ln,lt] = harv_debt[ln,lt] / harv_c[ln,lt]  
        for (t in -i:i) {
          for (n in -i:i) {
            n = adjust_bound(ln, n, nlon_out)
            t = adjust_bound(lt, t, nlat_out)
            if(!is.na(surplus_map[ln+n,lt+t])){
              redis_factor[ln+n,lt+t] = - redis_factor[ln,lt] * 
                surplus_map[ln+n,lt+t] / hotspot_for_c
              # Update surplus map
              surplus_map[ln+n,lt+t] = forest_c[ln+n,lt+t] -
                harv_c[ln+n,lt+t] + redis_factor[ln+n,lt+t] * harv_c[ln,lt]
            }
          }
        }
        # Update harvest debt map
        harv_debt[ln,lt] = 0.0
        # Adjust the Time Series of harvest rate
        for (yr in 1:nyr){
          harv_ts[,,yr] = harv_ts[,,yr] + harv_ts[ln,lt,yr] * redis_factor[,]
        }
      } # End of harvest debt check
    }
  } # End of cell loop

  return(harv_ts)
  
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
land_file = "C:/Users/sshu3/anaconda_wkspace/FATES/landuse.timeseries_4x5_hist_simyr1850-2015_200311_biomass_harvest.nc"
biom_den_file = "C:/Users/sshu3/anaconda_wkspace/FATES/fates.forest.vegc.avg_2000_2009.nc"
out_land_file = "C:/Users/sshu3/anaconda_wkspace/Fates/landuse.timeseries_4x5_hist_harmonized_simyr1850-2015.biomass_harvest.nc"
start_year_in = 1850
diag_plots = FALSE
num_cores = detectCores()
# kilograms to Petagrams C
k2P = 1/1000000000000

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

# Sum harvested C from 1850 to 2015 and calculate the density
harv_rate_acc = array(dim=c(num_harvest, nlon_out, nlat_out))
harv_rate_acc[,,]= 0
harv_rate_den = array(dim=c(num_harvest, nlon_out, nlat_out))
harv_rate_den[,,]= 0
for (h in 1:num_harvest) {
  harv_rate_acc[h,,] = rowSums(harv_rate_in[h,,,], dims = 2)
  harv_rate_den[h,,] = harv_rate_acc[h,,]*1e-6 / area_grid 
}

# Read in GCAM agro-ecological zone (AEZ) file
gcam_reg <- brick(gcam_reg_file)
hdr(gcam_reg, format = 'ENVI')
PROJ4_STRING = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
# PROJ4_STRING=crs(gcam_reg)

# # Upscale to 4 by 5 degs
# cat("\nstarting cell-by-cell mapping for GCAM regional mask \n")
# mcout = mclapply(c(0:(ncells_out-1)), function(i) upscale_cell(cell = i, lonW, lonE, latS, latN, PROJ4_STRING, gcam_reg), mc.cores = 1)
# cat("finishing cell-by-cell mapping. \n")
# 
# luh2elm_out = rbind.fill(mcout)
# luh2elm_out = luh2elm_out[order(luh2elm_out$cell),]
# gcam_out = array(dim=c(nlon_out, nlat_out))
# gcam_out[,] = luh2elm_out[, "layer"]



# 2. Read in biomass density from ELM-FATES
bioid = nc_open(biom_den_file)
forestc_den = ncvar_get(bioid, varid = "FATES_VEGC_FOREST", start=c(1,1), count=c(nlon_out,nlat_out))
nc_close(bioid)

# Quick diagnosis of the harvest rate diff
if(diag_plots) {
  diff_den = harv_rate_den[1,,] - forestc_den
  pdf(file=paste0(outpath, "check_harvrate_luh2_minus_fatesvegc.pdf"), width=8.5, height=4)
  cols <- brewer.pal(3, "RdYlGn")
  pal <- colorRampPalette(cols)
  print(levelplot(diff_den, margin=FALSE, at=seq(-10, 10, 1), col.regions = pal, colorkey = list(title = list("[LUH2 - FATES, Primary]", cex = 1, fontface = "bold", col = "red"))))
  dev.off()
}

# 3. Harmonize LUH2 harvest rate with FATES forest C distribution
new_harv_ts = remap_harv_rate(nlon_out, nlat_out, nyears_out, forestc_den, harv_rate_den[1,,], harv_rate_in[1,,,])

# 4. Validate the new harvest rate time series vs. old one.
harvest_in_globe_series = array(dim=c(nyears_out))
harvest_in_globe_series[] = 0
harvest_out_globe_series = array(dim=c(nyears_out))
harvest_out_globe_series[] = 0
for (yind_out in 1:(nyears_out-1)) {
  # match the in year to the out year
  harvest_in_globe_series[yind_out] = sum(harv_rate_in[1,,,yind_out], na.rm=TRUE) * k2P
  harvest_out_globe_series[yind_out] = sum(new_harv_ts[,,yind_out], na.rm=TRUE) * k2P
}
plot(years_out[1:(nyears_out-1)], harvest_in_globe_series[1:(nyears_out-1)], type = "l", lty = 1, main = "Annual total global biomass harvest",
     xlab = "Year", ylab = "Harvested carbon in kg C")
lines(years_out[1:(nyears_out-1)], harvest_out_globe_series[1:(nyears_out-1)], lty=2, col = "red")
legend("topleft", legend = c("LUH2 Primary Forest Harvest", "Harmonized Primary Forest Harvest"), lty=c(1,2), col = c("black", "red"), cex = 0.8)

# 5. Overwrite the wood harvest TS
# Shijie: These code need to be revised later
harv_rate_in[1,,,] = new_harv_ts[,,]
for (h in 2:num_harvest) {
  harv_rate_in[h,,,] = 0.0
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
  ncvar_put(lsid, varid = ls_harv_names[h], vals = harv_rate_in[h,,,], start = c(1,1,1), count = c(nlon_out, nlat_out, nyears_out))
} # end for loop over h for updating output file	

nc_close(lsid)

cat("Finish remap_harvest_rate at", date(), "\n")


# # projection info for all files, data, and rasters
# PROJ4_STRING = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
# 
# # precision tolerance for error checking
# tol = 1e-4
# 
# 
# # check that griddata file and land use file have same lat-lon cell centers
# if (nlon_out != nlon_grid | nlat_out != nlat_grid) {
#   stop(paste("Error: grid file", grid_file, "lat", nlat_grid, " and lon", nlon_grid, "dims do not match land surface file",
#              land_file, "lat", nlat_out, " and lon", nlon_out, "dims"))
# } else {
#   d = latixy_grid - latixy_out
#   if (sum(d) != 0) {
#     stop(paste("Error: grid latixy do not match land surface latixy"))
#   }
#   d = longxy_grid - longxy_out
#   if (sum(d) != 0) {
#     stop(paste("Error: grid longxy do not match land surface longxy"))
#   }
# } # end if-else check on grid specs
# 
# # open the luh2 file and get lat-lon bounds and number of years
# # there is no year variable, the first time record depends on whether it is historical or future
# # the data start in the upper left corner (-180,90 corner), with 0.25 degree resolution
# luh2id = nc_open(luh2_file)
# nlon_in = luh2id$dim$lon$len
# nlat_in = luh2id$dim$lat$len
# nyears_in = luh2id$dim$time$len
# nbounds_in = luh2id$dim$bounds$len
# lat_bounds_in = ncvar_get(luh2id, varid = "lat_bounds", start=c(1,1), count=c(nbounds_in, nlat_in))
# lon_bounds_in = ncvar_get(luh2id, varid = "lon_bounds", start=c(1,1), count=c(nbounds_in, nlon_in))
# 
# # set up some arrays
# harvest_in = array(dim=c(num_harvest, nlon_in, nlat_in))
# harvest_in_globe = array(dim=c(num_harvest))
# harvest_in_globe_series = array(dim=c(nyears_out))
# harvest_in_globe_series[] = 0
# harvest_out = array(dim=c(num_harvest, nlon_out, nlat_out, nyears_out))
# harvest_out_globe = array(dim=c(num_harvest))
# harvest_out_globe_series = array(dim=c(nyears_out))
# harvest_out[,,,]= 0
# harvest_out_globe_series[] = 0
# 
# # loop over years to read in and map the data
# # collect the output data in an array for later transfer to file
# for (yind_out in 1:(nyears_out-1)) {
#   # match the in year to the out year
#   year_out = years_out[yind_out]
#   yind_in = year_out - start_year_in + 1
#   
#   # read this year's in data and convert to rasters and put into a stack
#   # need to flip the lat so that the matrix origin is at the lower left corner
#   harvest_in[,,] = 0
#   harvest_in_globe[] = 0
#   harvest_out_globe[] = 0
#   harvest_in_rast_stack=stack()
#   for (h in 1:num_harvest) {
#     harvest_in[h,,] = ncvar_get(luh2id, varid=luh2_bioharv_names[h], start=c(1, 1, yind_in), count=c(nlon_in, nlat_in, 1))
#     harvest_in_t = t(harvest_in[h,,])
#     harvest_rast = raster(x=harvest_in_t, xmn=-180, ymn=-90, xmx=180, ymx=90, crs=PROJ4_STRING)
#     harvest_in_rast_stack = addLayer(harvest_in_rast_stack , harvest_rast)
#     # sum global harvest for comparison with out
#     harvest_in_globe[h] = sum(harvest_in[h,,], na.rm=TRUE)
#   } # end for loop over h for in data
#   # set the raster layer names
#   names(harvest_in_rast_stack) <- luh2_bioharv_names
#   
#   # process each out cell in parallel
#   # cell index starts at zero to avoid an extra calculation inside the function
#   # mclapply forks the calling environment and makes the data and code available to each process
#   cat("\nstarting cell-by-cell mapping for year", year_out, "at", date(), "\n")
#   mcout = mclapply(c(0:(ncells_out-1)), function(i) upscale_cell(cell = i, lonW, lonE, latS, latN, PROJ4_STRING, luh2_bioharv_names, harvest_in_rast_stack), mc.cores = 1)
#   cat("finishing cell-by-cell mapping (before sorting and writing) for year", year_out, "at", date(), "\n")
#   
#   # merge the out data and order by cell index
#   # each mcout is a one-record df with the cell index from above and
#   #     the weighted sum of the biomass harvest in kg C for each type
#   #     with the type names from luh2_bioharv_names
#   luh2elm_out = rbind.fill(mcout)
#   luh2elm_out = luh2elm_out[order(luh2elm_out$cell),]
#   
#   #cat("Before storage loop", year_out, "\n")
#   
#   # store the out data and sum this year for globe
#   # NAs have been removed in the function above
#   # compare in and out global harvest in kg C for this year
#   # need to store this is in the year+1 index
#   for (h in 1:num_harvest) {
#     harvest_out[h,,,yind_out+1] = luh2elm_out[, luh2_bioharv_names[h]]
#     harvest_out_globe[h] = sum(harvest_out[h,,,yind_out+1])
#     
#     # duplicate the first record in the output data
#     if (yind_out == 1) {
#       harvest_out[h,,,yind_out] = harvest_out[h,,,yind_out+1]
#     }
#     
#     if (harvest_out_globe[h] != harvest_in_globe[h]) {
#       # check for precision tolerance
#       # if ( abs(harvest_out_globe[h] - harvest_in_globe[h]) >= tol) {
#       # 	stop(paste("Error:", year_out, "global out harvest", harvest_out_globe[h],
#       # 			"kg C does not equal global in harvest", harvest_in_globe[h],
#       # 			"for harvest type", luh2_bioharv_names[h]))
#       # }
#     }
#     
#   } # end for loop over h for out data
#   
#   #cat("Before annual tracking", year_out, "\n")
#   
#   # keep track of and check the annual sum across harvest types in Pg C
#   # note that these data have not been year shifted
#   harvest_in_globe_series[yind_out] = sum(harvest_in_globe) * k2P
#   harvest_out_globe_series[yind_out] = sum(harvest_out_globe) * k2P
#   if (harvest_out_globe_series[yind_out] != harvest_in_globe_series[yind_out]) {
#     # Bypassed # check for precision tolerance
#     # if ( abs(harvest_out_globe_series[yind_out] - harvest_in_globe_series[yind_out]) >= tol) {
#     # 	stop(paste("Error:", year_out, "total global out harvest", harvest_out_globe_series[yind_out],
#     # 				"Pg C does not equal total global in harvest", harvest_in_globe_series[yind_out]))
#     # }
#   }
#   
#   cat("Finish year", year_out, "at", date(), "\n")
#   
# } # end for loop over out years (yind_out)
# 
# nc_close(luh2id)
# 
# # plot the global total in and out time series in Pg C
# # these data have not been year shifted, so drop the last element (2015 has no data)
# plot(years_out[1:(nyears_out-1)], harvest_in_globe_series[1:(nyears_out-1)], type = "l", lty = 1, main = "Annual total global biomass harvest",
#      xlab = "Year", ylab = "Harvested carbon in kg C")
# lines(years_out[1:(nyears_out-1)], harvest_out_globe_series[1:(nyears_out-1)], lty=2, col = "red")
# legend("topleft", legend = c("LUH2 file", "ELM file"), lty=c(1,2), col = c("black", "red"), cex = 0.8)
# 
# # this overwrites an existing out_land_file
# file.copy(land_file, out_land_file, overwrite=TRUE)
# 
# cat("Before writing new file\n")
# 
# # rename original fraction data and add biomass harvest data with original names
# lsid = nc_open(out_land_file, write=TRUE)
# 
# # get the dimension objects
# lon_dim = lsid$dim[['lsmlon']]
# lat_dim = lsid$dim[['lsmlat']]
# time_dim = lsid$dim[['time']]
# 
# for (h in 1:num_harvest) {
#   # rename fraction data
#   lsid <- ncvar_rename(lsid, ls_harv_names[h], ls_harv_names_frac[h])
#   
#   # define and add new variables; missing value is 0
#   new_var <- ncvar_def(ls_harv_names[h], "kg C per year", list(lon_dim, lat_dim, time_dim), prec="double")
#   lsid <- ncvar_add(lsid, new_var)
#   
#   # add biomass harvest data
#   ncvar_put(lsid, varid = ls_harv_names[h], vals = harvest_out[h,,,], start = c(1,1,1), count = c(nlon_out, nlat_out, nyears_out))
#   
# } # end for loop over h for updating output file	
# 
# nc_close(lsid)
# 
# cat("Finish proc_luh2_biomass_harv at", date(), "\n")
# 


