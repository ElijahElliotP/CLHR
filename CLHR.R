# CLHR R code, version 4

# Change History ############################################

# 1. JC saved file with new name and gave it version 2.0
# 2. Eli removed the need for GRASS with the help of Jon Clayden, writes CSVs V3.3
# 3. Write GeoTiffs V3.4
# 4. Loop structure in place
# 5. Some error handling added v3.6.3
# 6. Reads .asc or .txt, writes .asc, writes a Summary Sesult CSV, bugfixes,
#   separate preferences and functions file v3.7
# 7. The entire script has been broken into major pieces. 3 files.

# To-do: JPEG export
# To-do: vector ingestion?

# File structure input constraints.
# Region folder names may have prefixes or suffixes to their resolutions.
# Resolution numbers must be 1) integers 2) unique for each folder.
# Observations or "Individuals" folder names may have prefixes or suffixes.
# Observation prefixes are used to id different individuals or species
# if a sequence of them is to be run. Also useful in avoiding overwriting work.
# Region and Individuals folder names must have corresponding resolution numbers.
# Filenames must match folder names.
# All individuals subfolders must be in the same location.

# R-environment prep ####
library(raster)
library(rgdal)
library(data.table)
library(igraph)
library(gtools)
library(stringr)
library(mmand)
rm(list = ls())
timeStart <- proc.time()

# FUNCTIONS ####
source(file = 'CLHR-functions.R')

# Configuration ####
source(file = 'CLHR-config.txt')

allRegionFolders    <- mixedsort(dir(regionFolder, all.files = F, full.names = F))
allIndivFolders    <- mixedsort(dir(individualsFolder, all.files = F, full.names = F))

##### THE LARGE LOOP #####
for (pRS in seq_along(resolutions)) {

  RS <- resolutions[pRS]
  dir.create(file.path(CLHRDirectory, RS), showWarnings = F, recursive = T)
  regionMSTFile <- str_c(file.path(CLHRDirectory, RS), '/', RS, "-Region-MSTs.csv", sep = '')
    if (!file.exists(regionMSTFile)) {
    sink(file = regionMSTFile)
    cat('Individual,MST_Lengths_(px),MST_Distance_(m)\n')
    sink()
    } else {
      file.remove(regionMSTFile)
      sink(file = regionMSTFile)
      cat('Individual,MST_Lengths_(px),MST_Distance_(m)\n')
      sink()
    }
  animsByRes <- grep(str_c("\\D?", RS, "(?![0-9])", sep = ''), allIndivFolders, value = T, perl = T)
  currentFolder <- grep(str_c("\\D?", RS, "(?![0-9])", sep = ''), allRegionFolders, value = T, perl = T)
  if (length(currentFolder) > 1) {
    message("The current resolution specified matches more than one region folder.", appendLF = T)
    message("Please name folder inputs so they contain a unique integer.", appendLF = T)
    message("e.g. '5', '10', '50', '100', '365'.", appendLF = T)
    message("They must have this number in common with the 'individuals' files names.", appendLF = T)
    stop(currentFolder)
  }
  locationFilepath <- file.path(regionFolder, currentFolder)
  locationFilename <- str_c(locationFilepath, '/', currentFolder, ascExtension, sep = '')
  cat('Starting:', locationFilename)
  plot(1, main = str_c('Region Resolution:', RS, sep = " "), type = "n", axes = F, xlab = "", ylab = "")

  # Import location ----------------------------------------------------------------
  regionTile <- tryCatch({
    prepASCtoMatrix.fn(locationFilename) # this function is for an ascii file
    },
    error = function(err) {
      print("The script couldn't find a region folder or file," )
      print("or the region folders and files do not have the same name." )
      print("e.g. A folder 'Boston50' must contain the ASCII 'Boston50.txt' OR 'Boston50.asc'.")
      stop(err)
    })

  m <- dim(regionTile)[1]
  n <- dim(regionTile)[2]

  {
    if (f.drawMaps | f.writeGeoTIF | f.writeASC) { # Rasterize location (regionTile) and retrieve its extents & crs
      rastertile <- copySourceGeoSpatToRas.fn(regionTile, locationFilename)
    }
    if (f.writeGeoTIF | f.writeCSV | f.writeASC) { # create the output directory
      theDir <- file.path(CLHRDirectory, RS, '1.Region')
      dir.create(theDir, showWarnings = F, recursive = T)
    }
    if (f.drawMaps) { # plot the location
      plot(rastertile, col = c('green', 'blue'), main = 'Region - Full', xlab = "Easting", ylab = "Northing", legend = F)
    }
    if (f.writeGeoTIF) { # Write location to gtif
      rasterToGTIFF.fn(rastertile, theDir, "Region", 'full')
    }
    if (f.writeCSV) { # Write location to csv
      matrixToSimpleCSV.fn(regionTile, theDir, "Region", 'full')
    }
    if (f.writeASC) {
      writeOutAsc.fn(rastertile, theDir, "Region", 'full')
    }
  }

  # Eroding location matrix to majorSkeleton ----------------------------------------------
  majorSkeleton <- extractMajorSkeleton.fn(regionTile)

  {
    if (f.drawMaps | f.writeGeoTIF | f.writeASC) { # Rasterize majorSkeleton and apply the extents & crs
      skelRas <- copySourceGeoSpatToRas.fn(majorSkeleton, locationFilename)
    }
    if (f.writeGeoTIF | f.writeCSV | f.writeASC) { # create the output directory
      theDir <- file.path(CLHRDirectory, RS, '2.Skeleton')
      dir.create(theDir, showWarnings = F, recursive = T)
    }
    if (f.drawMaps) { # plot the majorSkeleton
      plot(skelRas, col = c('green', 'blue'), main = 'Region - Skeleton', xlab = "Easting", ylab = "Northing", legend = F)
    }
    if (f.writeGeoTIF) { # Write majorSkeleton to gtif
      rasterToGTIFF.fn(skelRas, theDir, "Region", 'skeleton')
    }
    if (f.writeCSV) { # Write majorSkeleton to csv
      matrixToSimpleCSV.fn(majorSkeleton, theDir, "Region", 'skeleton')
    }
    if (f.writeASC) {
      writeOutAsc.fn(skelRas, theDir, "Region", 'skeleton')
    }
  }

  baseConnectivity <- connectivityArray.fn(majorSkeleton)

  #### OBSERVATIONS SUBLOOP ####
  for (pAS in seq_along(animsByRes)) {
    AS <- animsByRes[pAS]
    truncAni <- str_locate(pattern = as.character(RS), string = AS)
    indivNickname <- (str_sub(AS, start = (truncAni[1] - truncAni[1] + 1), end = truncAni[2]))

    dataFilepath <- file.path(individualsFolder, AS) # species location info
    dataFilename <- str_c(basename(dataFilepath), ascExtension)

    cat('\nStarting:', str_c(dataFilepath, ascExtension), indivNickname, '\n')
    indivTileFile <- file.path(dataFilepath, dataFilename)
    if (file.exists(file.path(CLHRDirectory, RS, indivNickname, '/.done')) == TRUE) {
      cat(indivNickname, 'has already run. \nSkipping. If this is in error please delete the ".done" file.\n')
      next
    }
    # Import observations -------------------------------------------------
    sdata = tryCatch({
      # as.matrix(read.table(file.path(dataFilepath, dataFilename), header = F, skip = 6, sep = " "))
      prepASCtoMatrix.fn(indivTileFile)
    },
    error = function(err) {
      print("The name of the individual should preceed the resolution and" )
      print("the individuals folders and files must have the same name." )
      print("e.g. the folder 'SnapTurt50' should contain the ascii 'SnapTurt50.txt'.")
      stop(err)
    })

    numOfIndivs <- sum(sdata)

    {
      if (f.drawMaps | f.writeGeoTIF | f.writeASC) { # Rasterize the observation points (sdata) and apply the extents & crs
        sdataRas <- copySourceGeoSpatToRas.fn(sdata, locationFilename)
      }
      if (f.writeGeoTIF | f.writeCSV | f.writeASC) { # create the output directory
        theDir <- file.path(CLHRDirectory, RS, indivNickname, '1.IndividualLocations')
        dir.create(theDir, showWarnings = F, recursive = T)
      }
      if (f.drawMaps) { # plot the the observation points
        plot(sdataRas, main = str_c(indivNickname,'Locations', sep = ' '),
             mtext(str_c(numOfIndivs, "individuals", sep = ' ')),
             xlab = "Easting", ylab = "Northing", col = c('white', 'red'), legend = F)
      }
      if (f.writeGeoTIF) { # Write the observation points to gtif
        rasterToGTIFF.fn(sdataRas, theDir, indivNickname, 'ObsPtsOnly')
      }
      if (f.writeCSV) { # Write the observation points to csv
        matrixToSimpleCSV.fn(sdata, theDir, indivNickname, 'ObsPtsOnly')
      }
      if (f.writeASC) { # Write the observation points to asc
        writeOutAsc.fn(sdataRas, theDir, indivNickname, 'ObsPtsOnly')
      }
    }

    # Observations on FULL LANDSCAPE
    mapIndivsOnLandFull = regionTile
    mapIndivsOnLandFull[sdata == 1] = 10  # where there is an individual
    mapIndivsOnLandFull[mapIndivsOnLandFull == 1] = 11 # where there is "habitat"
    mapIndivsOnLandFull[mapIndivsOnLandFull == 0] = 12 # where there is "no habitat"

    {
      if (f.drawMaps | f.writeGeoTIF | f.writeASC) {
        mapRas <- copySourceGeoSpatToRas.fn(mapIndivsOnLandFull, locationFilename)
      }
      if (f.writeGeoTIF | f.writeCSV | f.writeASC) { # create the output directory
        theDir <- file.path(CLHRDirectory, RS, indivNickname, '1.IndividualLocations')
        dir.create(theDir, showWarnings = F, recursive = T)
      }
      if (f.drawMaps) { # plot the Obs Points on the Landscape
        plot(mapRas, main = str_c(indivNickname, "Locations in Region", sep = ' '),
             mtext(str_c(numOfIndivs, "individuals", sep = ' ')),
             xlab = "Easting", ylab = "Northing", col = c('red', 'blue', 'green'), legend = F)
      }
      if (f.writeGeoTIF) { # Write Obs-Points-on-Landscape to gtif
        rasterToGTIFF.fn(mapRas, theDir, indivNickname, 'IndividualsInRegion')
      }
      if (f.writeCSV) { # Write Obs-Points-on-Landscape to csv
        matrixToSimpleCSV.fn(mapIndivsOnLandFull, theDir, indivNickname, 'IndividualsInRegion')
      }
      if (f.writeASC) { # Write Obs-Points-on-Landscape to asc
        writeOutAsc.fn(mapRas, theDir, indivNickname, 'IndividualsInRegion')
      }
    }

    # observations on the skeleton
    mapIndivsOnLandSkel <- shiftDataToSkeleton.fn(sdata, majorSkeleton)

    {
      if (f.drawMaps | f.writeGeoTIF | f.writeASC) {
        map2Ras <- copySourceGeoSpatToRas.fn(mapIndivsOnLandSkel, locationFilename)
      }
      if (f.writeGeoTIF | f.writeCSV | f.writeASC) { # create the output directory
        theDir <- file.path(CLHRDirectory, RS, indivNickname, '2.LocationsOnSkeleton')
        dir.create(theDir, showWarnings = F, recursive = T)
      }
      if (f.drawMaps) { # plot the Obs Points on the Skeleton
        plot(map2Ras, main = str_c(indivNickname,'Locations on Region Skeleton', sep = ' '),
             mtext(str_c(numOfIndivs, "individuals", sep = ' ')),
             xlab = "Easting", ylab = "Northing", col = c('red', 'blue', 'green'), legend = F)
      }
      if (f.writeGeoTIF) { # Write Obs Points on the Skeleton to gtif
        rasterToGTIFF.fn(map2Ras, theDir, indivNickname, 'LocationsOnSkeleton')
      }
      if (f.writeCSV) { # Write Obs Points on the Skeleton to csv
        matrixToSimpleCSV.fn(mapIndivsOnLandSkel, theDir, indivNickname, 'LocationsOnSkeleton')
      }
      if (f.writeASC) { # Write Obs-Points-on-Landscape to asc
        writeOutAsc.fn(map2Ras, theDir, indivNickname, 'LocationsOnSkeleton')
      }
    }

    # Counting points ---------------------------------------------------------

    # Numbering data points and joints -----------------------------------------------

    serializedJointsNPointsList <- makeJointsNPointsSequential.fn(baseConnectivity, shiftedPoints)

    joints            <- serializedJointsNPointsList$aa
    seqDataNum        <- serializedJointsNPointsList$bb
    distanceSize      <- serializedJointsNPointsList$cc
    seqShiftedPts     <- serializedJointsNPointsList$dd

    # the largest part of the operation
    mapIndivsOnMST <- calculateTheMinimumSpanningTree.fn()

    {
      if (f.drawMaps | f.writeGeoTIF | f.writeASC) {
        map3Ras <- copySourceGeoSpatToRas.fn(mapIndivsOnMST, locationFilename)
      }
      if (f.writeGeoTIF | f.writeCSV | f.writeASC) { # create the output directory
        theDir <- file.path(CLHRDirectory, RS, indivNickname, '3.MinimumSpanningTree')
        dir.create(theDir, showWarnings = F, recursive = T)
      }
      if (f.drawMaps) { # plot the MST and data points on skeleton
        plot(map3Ras, main = str_c(indivNickname, 'Minimum Spanning Tree', sep = ' '),
             mtext(str_c(numOfIndivs, "individuals", sep = ' ')),
             xlab = "Easting", ylab = "Northing", col = c('red', 'pink', 'blue', 'green'), legend = F)
      }
      if (f.writeGeoTIF) { # Write MST and data points on skeleton to gtif
        rasterToGTIFF.fn(map3Ras, theDir, indivNickname, 'MST')
      }
      if (f.writeCSV) { # Write MST and data points on skeleton to csv
        matrixToSimpleCSV.fn(mapIndivsOnMST, theDir, indivNickname, 'MST')
      }
      if (f.writeASC) {
        writeOutAsc.fn(map3Ras, theDir, indivNickname, 'MST')
      }

      # Output distance of MST ####
      if (f.writeResultText) {
      dir.create(file.path(CLHRDirectory, RS, indivNickname), showWarnings = F, recursive = T)
      sink(str_c(file.path(CLHRDirectory, RS, indivNickname), '/', indivNickname, "-MST.txt"))
      cat("The Minimum Spanning Tree (px) is:", mstdistance, "# of Indivs:", numOfIndivs, '\n')
      sink(file = str_c(file.path(CLHRDirectory, RS), '/', RS, "-Region-MSTs.csv"), append = T, split = T)
      cat(c(indivNickname, ",", mstdistance, ",", mstdistance*as.numeric(RS), "\n"))
      } else {
        sink(file = str_c(file.path(CLHRDirectory, RS), '/', RS, "-Region-MSTs.csv"), append = T)
        cat(c(indivNickname, ",", mstdistance, ",", mstdistance*as.numeric(RS), "\n"))
      }
      sink.reset()
    }
    cat('The MST distance for', indivNickname, 'is:', mstdistance*as.numeric(RS), "m. # of Indivs: ", numOfIndivs, '\n')
    invisible(file.create(file.path(CLHRDirectory, RS, indivNickname, '/.done')))
  }
}
timeEnd <- proc.time()
runtime <- timeEnd - timeStart
