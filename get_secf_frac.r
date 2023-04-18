# get_secf_frac.r

# this script overlap the harvest rate before 1700 using secondary forest area 
# fraction and rewrite the harvest rate at year 1700.
# For area-based harvest calculation only.

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

library(raster)
library(rasterVis)
library(RColorBrewer)
library(terra)
library(ncdf4)
library(parallel)
library(plyr)

##### the single cell function for parallelizing must be defined first !!!!!!!

################# luh2elm single-cell grid mapping function

# assume that mclapply is being used, which utilizes fork()
# this means that all of the calling environment is forked to each child
# such that the memory is shared and available without exporting or passing
# copies are made only when variables are modified
luh2elm_cell <- function(cell, lonW, lonE, latS, latN, p4string, type_names, raster_stack)	{
  # get the current out cell bounds
  # this assumes that cell starts at 0 to eliminate a calc inside
  lt = trunc(cell / nlon_out) + 1
  ln = cell %% nlon_out + 1
  # if (ln > 1 & lt > 1) {
  lnmin = lonW[ln,lt]
  lnmax = lonE[ln,lt]
  ltmin = latS[lt,lt]
  ltmax = latN[lt,lt]
  #     } else {
  #         lnmin = lonW
  # 	    lnmax = lonE
  # 	    ltmin = latS
  # 	    ltmax = latN
  #     }
  
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
  
  # calculate the weighted sum for each harvest type (in case there are fractions of input cells in the output cell)
  # put the values in a single row of a data frame with the cell index for id
  harv_cell_df = extract(raster_stack, spobj, weights=TRUE, normalizeWeights=FALSE, cellnumbers=TRUE, df=TRUE, small=TRUE, exact=TRUE, na.rm=TRUE)
  out_df = harv_cell_df[1,c("cell", type_names)]
  out_df$cell = cell
  for (h in 1:num_secf) {
    out_df[1,type_names[h]] = mean(harv_cell_df[,type_names[h]] * harv_cell_df$weight, na.rm=TRUE)
  }
  
  return(out_df)
  
} # end single-cell mapping function

############ change this to reflect correct directory structure and files ! ###################
luh2_file = "D:/LUH2/states_hist_hi.nc"
grid_file = "C:/Users/sshu3/anaconda_wkspace/griddata_4x5_060404.nc"
land_file = "C:/Users/sshu3/anaconda_wkspace/landuse.timeseries_4x5_HIST_hi_simyr1700-2015_c03242023.nc"
out_land_file = "C:/Users/sshu3/anaconda_wkspace/Fates/landuse.timeseries_4x5_HIST_hi_simyr1700.accumulated.nc"
start_year_in = 850

cat("Start get_secf_frac at", date(), "\n")

# the number of secondary forest and non-forest area fraction
num_secf = 2
luh2_secf_names = c("secdf", "secdn")

ls_harv_names = c("HARVEST_VH1", "HARVEST_VH2")

# projection info for all files, data, and rasters
PROJ4_STRING = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

num_cores = detectCores()

# precision tolerance for error checking
tol = 1e-4

# open the land surface file and get years present and lat-lon
lsid = nc_open(land_file)
nlon_out = lsid$dim$lsmlon$len
nlat_out = lsid$dim$lsmlat$len
nyears_out = lsid$dim$time$len
years_out = ncvar_get(lsid, varid = "YEAR", start=c(1), count=c(nyears_out))
latixy_out = ncvar_get(lsid, varid = "LATIXY", start=c(1,1), count=c(nlon_out,nlat_out))
longxy_out = ncvar_get(lsid, varid = "LONGXY", start=c(1,1), count=c(nlon_out,nlat_out))
ncells_out = nlon_out * nlat_out
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
nc_close(gdid)

# check that griddata file and land use file have same lat-lon cell centers
if (nlon_out != nlon_grid | nlat_out != nlat_grid) {
  stop(paste("Error: grid file", grid_file, "lat", nlat_grid, " and lon", nlon_grid, "dims do not match land surface file",
             land_file, "lat", nlat_out, " and lon", nlon_out, "dims"))
} else {
  d = latixy_grid - latixy_out
  if (sum(d) != 0) {
    stop(paste("Error: grid latixy do not match land surface latixy"))
  }
  d = longxy_grid - longxy_out
  if (sum(d) != 0) {
    stop(paste("Error: grid longxy do not match land surface longxy"))
  }
} # end if-else check on grid specs

# open the luh2 file and get lat-lon bounds and number of years
# there is no year variable, the first time record depends on whether it is historical or future
# the data start in the upper left corner (-180,90 corner), with 0.25 degree resolution
luh2id = nc_open(luh2_file)
nlon_in = luh2id$dim$lon$len
nlat_in = luh2id$dim$lat$len
nyears_in = luh2id$dim$time$len
nbounds_in = luh2id$dim$bounds$len
lat_bounds_in = ncvar_get(luh2id, varid = "lat_bounds", start=c(1,1), count=c(nbounds_in, nlat_in))
lon_bounds_in = ncvar_get(luh2id, varid = "lon_bounds", start=c(1,1), count=c(nbounds_in, nlon_in))

# Select the year period from year 1 to the previous year of the begin year
# from land use timeseries file (land_file)
year_beg = 1
year_end = years_out[1]- start_year_in
y_len = year_end - year_beg + 1

# set up some arrays
harvest_in_tmp = array(dim=c(num_secf, nlon_in, nlat_in))
harvest_in = array(dim=c(num_secf, nlon_in, nlat_in))
harvest_in_globe = array(dim=c(num_secf))
harvest_in_globe_series = array(dim=c(nyears_out))
harvest_in_globe_series[] = 0
harvest_out = array(dim=c(num_secf, nlon_out, nlat_out, nyears_out))
harvest_out_globe = array(dim=c(num_secf))
harvest_out_globe_series = array(dim=c(nyears_out))
harvest_out[,,,]= 0
harvest_out_globe_series[] = 0

# read these years' in data, get primary forest and non-forest fraction
# as harvest rate and convert to rasters and put into a stack
# need to flip the lat so that the matrix origin is at the lower left corner
harvest_in_tmp[,,] = 0.
harvest_in[,,] = 0
harvest_in_globe[] = 0
harvest_out_globe[] = 0
harvest_in_rast_stack=stack()
iyr = years_out[1]- start_year_in
for (h in 1:num_secf) {
  harvest_in_tmp[h,,] = ncvar_get(luh2id, varid=luh2_secf_names[h], start=c(1, 1, iyr), count=c(nlon_in, nlat_in, 1))
  harvest_in[h,,] = harvest_in[h,,] + harvest_in_tmp[h,,]
  harvest_in <- ifelse(is.na(harvest_in), 0.0, harvest_in)
}

for (h in 1:num_secf) {
  harvest_in_t = t(harvest_in[h,,])
  harvest_rast = raster(x=harvest_in_t, xmn=-180, ymn=-90, xmx=180, ymx=90, crs=PROJ4_STRING)
  harvest_in_rast_stack = addLayer(harvest_in_rast_stack , harvest_rast)
  # sum global harvest for comparison with out
  harvest_in_globe[h] = sum(harvest_in[h,,], na.rm=TRUE)
} # end for loop over h for in data
# set the raster layer names
names(harvest_in_rast_stack) <- luh2_secf_names

# process each out cell in parallel
# cell index starts at zero to avoid an extra calculation inside the function
# mclapply forks the calling environment and makes the data and code available to each process
cat("\nstarting cell-by-cell mapping for the historical accumulation at ", date(), "\n")
mcout = mclapply(c(0:(ncells_out-1)), function(i) luh2elm_cell(cell = i, lonW, lonE, latS, latN, PROJ4_STRING, luh2_secf_names, harvest_in_rast_stack), mc.cores = 1)
cat("finishing cell-by-cell mapping (before sorting and writing) for the historical accumulation at", date(), "\n")

# merge the out data and order by cell index
# each mcout is a one-record df with the cell index from above and
#     the weighted sum of the biomass harvest in kg C for each type
#     with the type names from luh2_secf_names
luh2elm_out = rbind.fill(mcout)
luh2elm_out = luh2elm_out[order(luh2elm_out$cell),]

#cat("Before storage loop", year_out, "\n")

# store the out data and sum this year for globe
# NAs have been removed in the function above
# compare in and out global harvest in kg C for this year
# need to store this is in the year+1 index
harvest_out_tmp = array(dim=c(num_secf, nlon_out, nlat_out))
harvest_out_tmp[,,] = 0.
for (h in 1:num_secf) {
  harvest_out_tmp[h,,] = luh2elm_out[, luh2_secf_names[h]]
  harvest_out_globe[h] = sum(harvest_out_tmp[h,,])
  
  # if (harvest_out_globe[h] != harvest_in_globe[h]) {
  #   # check for precision tolerance
  #   # if ( abs(harvest_out_globe[h] - harvest_in_globe[h]) >= tol) {
  #   # 	stop(paste("Error:", year_out, "global out harvest", harvest_out_globe[h],
  #   # 			"kg C does not equal global in harvest", harvest_in_globe[h],
  #   # 			"for harvest type", luh2_secf_names[h]))
  #   # }
  # }
  
} # end for loop over h for out data

#cat("Before annual tracking", year_out, "\n")

cat("Finish aggregating accumulated forest harvest rate at", date(), "\n")

file.copy(land_file, out_land_file, overwrite=TRUE)
harvest_out[,,,2] = harvest_out_tmp
cat("Before writing new file\n")

# rename original fraction data and add biomass harvest data with original names
lsid = nc_open(out_land_file, write=TRUE)

# get the dimension objects
lon_dim = lsid$dim[['lsmlon']]
lat_dim = lsid$dim[['lsmlat']]
time_dim = lsid$dim[['time']]

for (h in 1:num_secf) {
  
  # add biomass harvest data
  ncvar_put(lsid, varid = ls_harv_names[h], vals = harvest_out[h,,,], start = c(1,1,1), count = c(nlon_out, nlat_out, nyears_out))
  
} # end for loop over h for updating output file	

nc_close(lsid)

cat("Finish writing new file\n")
