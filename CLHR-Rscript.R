# CLHR R code, version 3.5

# Change History ############################################

# 1. JC saved file with new name and gave it version 2.0
# 2. Eli removed the need for GRASS with the help of Jon Clayden, writes CSVs V3.3
# 3. Write GeoTiffs V3.4
# 4. Loop structure in place

# To-do: correct distance discrepency
# To-do: vector-to-raster?
# To-do: how to ingest ascii files that do not have a .txt extension

# Setting working directory:
# Personal preference but the old script had the data files at the same folder
# level as the code folder and I've changed this so all data and ouputs are within
# the script folder. - Eli

# Filenames should not have the units of measure ('m' for example) in the names.

# R-environment prep ####
library(raster)
library(rgdal)
library(data.table)
library(igraph) # masks 'union' from :raster and :base and both 'decompose' and 'spectrum' from :stats
library(gtools)
library(stringr)
library(maptools)
# temp load dev version of mmand
# devtools::install_github("jonclayden/mmand")
library(mmand)
# dev.off() # or:
# frame() # or
rm(list = ls())

# FUNCTIONS ####
{
  skeletonise3_2d.fn <- function(x) {
    k1 <- matrix(c(0, NA, 1, 0, 1, 1, 0, NA, 1), 3, 3)
    k2 <- matrix(c(NA, 1, NA, 0, 1, 1, 0, 0, NA), 3, 3)
    rot.fn <- function(x) {t(apply(x, 2, rev))}
    hom <- function(x, k) morph(x, k, operator = "==", merge = "all", value = 1)

    repeat
    {
      previous <- x
      for (i in 1:4)
      {
        x <- x & !hom(x, k1)
        storage.mode(x) <- "integer"
        x <- x & !hom(x, k2)
        storage.mode(x) <- "integer"
        k1 <- rot.fn(k1)
        k2 <- rot.fn(k2)
      }
      if (isTRUE(all.equal(x, previous)))
        return(x)
    }
  }

  extractMajorSkeleton.fn <- function(x, kernel) {
    chunks <- components(x, kernel)
    sort(table(chunks), decreasing = T)
    chunksMaxVal <- as.integer(row.names(as.matrix(sort(table(chunks), decreasing = T))))[1]
    chunksMajComponent <- (chunks == chunksMaxVal)*1
    chunksMajComponent[is.na(chunksMajComponent)] <- 0
    skeletonise3_2d.fn(chunksMajComponent)
  }

  prepASCtoMatrix.fn <- function(x) {
    x[x == -9999] = 0
    x <- as.data.table(x) # recommended to reduce system load for large data frames
    x <- x[, which(unlist(lapply(x, function(x)!all(is.na(x))))), with = F]
    x <- as.matrix(x)
  }

  rasterToGTIFF.fn <- function(x, directiory, animal, res, suffix) {
    SPDF <- as(x, 'SpatialPixelsDataFrame')
    writeGDAL(SPDF, fname = paste(directiory, '/', animal, '-', res, '-', suffix, '.tif', sep = ''), drivername = 'GTiff', type = 'Byte', mvFlag = 0, options = "TFW=YES")
  }

  matrixToSimpleCSV.fn <- function(x, directory, animal, res, suffix) {
    write.table(x, file = paste(directory, '/', animal, '-', res, '-', suffix, '.csv', sep = ''), append = F, sep = ',', row.names = F, col.names = F)
  }

  sink.reset <- function() { # courtsey of Dason from Wesley Chapel, Florida via Stackoverflow.com
    for (i in seq_len(sink.number())) {
      sink(NULL)
    }
  }
}

# Configuration ####
# a separate config file is probably a good idea at this point.
{
  .f.writeCSVs        <- T # slows script
  .f.writeGeoTIF      <- F # slows script
  .f.drawMaps         <- T # slows script
  # .f.writeJPGs        <- T # slows script... no doubt
  .f.writeResultText  <- T
  regionFolder        <- '1.Shore' # region
  individualsFolder   <- '2.Turtles' # individuals
  CLHRDirectory       <- '3.CLHRs' # the files can go anywhere, why not here?
  indivPrefix         <- '' # planned
  indivSuffix         <- 'output'
  regionPrefix        <- '' # planned
  regionSuffix        <- '' # planned
  resolutions         <- c(20)
  nresolutions        <- length(resolutions)
}

plot(1, main = 'Restart', type = "n", axes = F, xlab = "", ylab = "")
allRegionFolders    <- mixedsort(dir(regionFolder, all.files = F, full.names = F))
allAnimalFolders    <- mixedsort(dir(individualsFolder, all.files = F, full.names = F))

##### THE LARGE LOOP #####
for (pRS in seq_along(resolutions)) {
  # pRS <- 1
  RS <- resolutions[pRS]
  dir.create(file.path(CLHRDirectory, RS), showWarnings = F, recursive = T)
  regionMSTFile <- paste(file.path(CLHRDirectory, RS), '/', RS, "-Region-MSTs.csv", sep = '')
    if (!file.exists(regionMSTFile)) {
    sink(file = regionMSTFile)
    cat('Individual,MST_PixelLengths,MST_Distance_(m)\n')
    sink()
    } else {
      file.remove(regionMSTFile)
      sink(file = regionMSTFile)
      cat('Individual,MST_PixelLengths,MST_Distance_(m)\n')
      sink()
    }
  animsByRes <- grep(str_c('\\S', RS, '\\D', sep = ''), allAnimalFolders, value = T, perl = T)
  currentFolder <- grep(str_c('\\S', RS, '\\D', sep = ''), allRegionFolders, value = T, perl = T)
  locationFilepath <- file.path(regionFolder, currentFolder)
  locationFilename <- str_c(locationFilepath, '/', currentFolder, '.txt', sep = '')
  plot(1, main = paste('Region Resolution:', RS, sep = " "), type = "n", axes = F, xlab = "", ylab = "")

  # Import location ----------------------------------------------------------------
  tile = as.matrix(read.table(locationFilename, header = F, skip = 6, sep = " "))
  matrixtile <- prepASCtoMatrix.fn(tile) # this function is for an ascii tile in .txt format
  # this expects non-occurance points = -9999, occurance = 1 and NA values and return a binary matrix
  m <- dim(matrixtile)[1]
  n <- dim(matrixtile)[2]

  {
    if (.f.drawMaps | .f.writeGeoTIF) { # Rasterize location (matrixtile) and retrieve its extents & crs
      rastertile <- raster(matrixtile)
      extent(rastertile) <- readGDAL(locationFilename)
      crs(rastertile) <- crs(readGDAL(locationFilename)) # readGDAL reads the accompanying .prj file
    }

    if (.f.writeGeoTIF | .f.writeCSVs) { # create the output directory
      theDir <- file.path(CLHRDirectory, RS, '1.Region')
      dir.create(theDir, showWarnings = F, recursive = T)
    }

    if (.f.drawMaps) { # plot the location
      plot(rastertile, col = c('green', 'blue'), main = 'Region - Full', xlab = "Easting", ylab = "Northing", legend = F)
    }

    if (.f.writeGeoTIF) { # Write location to gtif
      rasterToGTIFF.fn(rastertile, theDir, "Region", RS, 'full')
    }

    if (.f.writeCSVs) { # Write location to csv
      matrixToSimpleCSV.fn(matrixtile, theDir, "Region", RS, 'full')
    }
  }

  # Eroding location matrix to skeleton ----------------------------------------------
  # rast.skeleton <- skeletonise3_2d.fn(matrixtile)
  rast.skeleton <- extractMajorSkeleton.fn(matrixtile, shapeKernel(c(3,3), type = "box"))

  {
    if (.f.drawMaps | .f.writeGeoTIF) { # Rasterize skeleton (rast.skeleton) and apply the extents & crs
      skelRas <- raster(rast.skeleton)
      extent(skelRas) <- extent(rastertile)
      crs(skelRas) <- crs(rastertile)
    }

    if (.f.writeGeoTIF | .f.writeCSVs) { # create the output directory
      theDir <- file.path(CLHRDirectory, RS, '2.Skeleton')
      dir.create(theDir, showWarnings = F, recursive = T)
    }

    if (.f.drawMaps) { # plot the skeleton
      plot(skelRas, col = c('green', 'blue'), main = 'Region - Skeleton', xlab = "Easting", ylab = "Northing", legend = F)
    }

    if (.f.writeGeoTIF) { # Write skeleton to gtif
      rasterToGTIFF.fn(skelRas, theDir, "Region", RS, 'skeleton')
    }

    if (.f.writeCSVs) { # Write skeleton to csv
      matrixToSimpleCSV.fn(rast.skeleton, theDir, "Region", RS, 'skeleton')
    }
  }
  # increment skeleton values if they are non-zero, based on the number of non-zero 8-way neighbours
  # this creates a map of 'connectedness'
  skeleton <- rast.skeleton
  for (i in 1:m) {
    for (j in 1:n) {
      if (skeleton[i, j] > 0) { # if the current points is not zero
        joint = 0
        if ((i > 1) & (j > 1)) { # any point excluding the first row and column
          if (skeleton[i - 1, j - 1] > 0) { # compare with up-left
            joint = joint + 1
          }
        }
        if (i > 1) { # any column excluding the first row
          if (skeleton[i - 1, j] > 0) { # compare with up
            joint = joint + 1
          }
        }
        if ((i > 1) & (j < n)) { # any point excluding the first row and excluding the last column
          if (skeleton[i - 1, j + 1] > 0) { # compare with up-right
            joint = joint + 1
          }
        }
        if (j < n) { # any row excluding the last column
          if (skeleton[i, j + 1] > 0) { # compare with right
            joint = joint + 1
          }
        }
        if ((i < m) & (j < n)) { # any point excluding the last row and column
          if (skeleton[i + 1, j + 1] > 0) { # compare with down-right
            joint = joint + 1
          }
        }
        if (i < m) { # any column excluding the last
          if (skeleton[i + 1, j] > 0) { # compare with down
            joint = joint + 1
          }
        }
        if ((i < m) & (j > 1)) { # any point excluding the last row and any point excluding the first column
          if (skeleton[i + 1, j - 1] > 0) { # compare with down-left
            joint = joint + 1
          }
        }
        if (j > 1) { # any column excluding the first
          if (skeleton[i, j - 1] > 0) { # compare with left
            joint = joint + 1
          }
        }
        if (joint > 0) {
          skeleton[i, j] = joint # increment the current point by 1 for every non-zero 8-way neighbour
        }
      }
    }
  }

  # OBSERVATIONS SUBLOOP
  for (pAS in seq_along(animsByRes)) {
    # pAS <- 6
    AS <- animsByRes[pAS]
    truncAni <- str_locate(pattern = as.character(RS), string = AS)
    animalNickname <- (str_sub(AS, start = (truncAni[1] - truncAni[1] + 1), end = truncAni[2]))
    dataFPPrep <- paste(animalNickname, indivSuffix, sep = '')
    dataFilepath <- file.path(individualsFolder, dataFPPrep) # species location info
    dataFilename <- paste(basename(dataFilepath), '.txt', sep = '')

    cat(locationFilename, animalNickname, dataFilepath, dataFilename, '\n')
    plot(1, main = paste('Individuals in Set:', animalNickname, sep = " "), type = "n", axes = F, xlab = "", ylab = "")

    # Just in case, clear previous variable assingments
    # rm(list = c("data", "dist", "dist2", "distances", "distances2", "distances3", "dpairs", "dpairsleft", "joints", "map", "map2", "map3", "minimum", "minpathways", "minpathways2", "minpathways3", "minrow", "mst", "paths", "check", "col", "D", "d", "datanum", "dgraph", "distsize", "i", "j", "jointnum", "jsize", "k", "l", "mapRas", "map2Ras", "map3Ras", "minD", "minx", "miny", "mstdistance", "new", "new.dist.row", "new.minpath", "new.paths", "path", "row", "sdata", "sdataRas", "segment", "sjd", "skel2", "skelseg", "t", "used", "x", "y", "z"))


    # Import observations -------------------------------------------------
    sdata <- as.matrix(read.table(file.path(dataFilepath, dataFilename), header = F, skip = 6, sep = " "))
    sdata <- prepASCtoMatrix.fn(sdata)
    {
      if (.f.drawMaps | .f.writeGeoTIF) { # Rasterize the observation points (sdata) and apply the extents & crs
        sdataRas <- raster(sdata)
        extent(sdataRas) <- extent(skelRas)
        crs(sdataRas) <- crs(skelRas)
      }

      if (.f.writeGeoTIF | .f.writeCSVs) { # create the output directory
        theDir <- file.path(CLHRDirectory, RS, animalNickname, '1.IndividualLocations')
        dir.create(theDir, showWarnings = F, recursive = T)
      }

      if (.f.drawMaps) { # plot the the observation points
        plot(sdataRas, main = paste(animalNickname,'Locations', sep = " "), xlab = "Easting", ylab = "Northing",
             col = c('white', 'red'), legend = F)
      }

      if (.f.writeGeoTIF) { # Write the observation points to gtif
        rasterToGTIFF.fn(sdataRas, theDir, animalNickname, RS, 'ObsPtsOnly')
      }

      if (.f.writeCSVs) { # Write the observation points to csv
        matrixToSimpleCSV.fn(sdata, theDir, animalNickname, RS, 'ObsPtsOnly')
      }
    }

    # Observations on FULL LANDSCAPE
    map = matrixtile
    map[sdata == 1] = 10  # where there is an individual
    map[map == 1] = 11    # where there is "habitat"
    map[map == 0] = 12    # where there is "no habitat"

    {
      if (.f.drawMaps | .f.writeGeoTIF) { # Rasterize Obs Points on the landscape (map) and apply the extents & crs
        mapRas <- raster(map)
        extent(mapRas) <- extent(sdataRas)
        crs(mapRas) <- crs(sdataRas)
      }

      if (.f.writeGeoTIF | .f.writeCSVs) { # create the output directory
        theDir <- file.path(CLHRDirectory, RS, animalNickname, '1.IndividualLocations')
        dir.create(theDir, showWarnings = F, recursive = T)
      }

      if (.f.drawMaps) { # plot the Obs Points on the Landscape
        plot(mapRas, main = paste(animalNickname,'Locations in Region', sep = " "),
             xlab = "Easting", ylab = "Northing", col = c('red', 'blue', 'green'), legend = F)
      }

      if (.f.writeGeoTIF) { # Write Obs-Points-on-Landscape to gtif
        rasterToGTIFF.fn(mapRas, theDir, animalNickname, RS, 'LocationsInRegion')
      }

      if (.f.writeCSVs) { # Write Obs-Points-on-Landscape to csv
        matrixToSimpleCSV.fn(map, theDir, animalNickname, RS, 'LocationsInRegion')
      }
    }

    # Moving data points to skeleton ####
    data = matrix(0, m, n) # a blank matrix with the same dimensions as everything else
    for (i in 1:m) {
      for (j in 1:n) { # move down each row and across every column
        if (sdata[i, j] == 1) {
          minD = Inf
          for (k in 1:m) {
            for (l in 1:n) {
              if (rast.skeleton[k, l] == 1) {
                D = (i - k)^2 + (j - l)^2 # should this be the square root of the result instead?
                if (D < minD) {
                  minD = D
                  miny = k
                  minx = l
                }
              }
            }
          }
          data[miny, minx] = 1
        }
      }
    } # small wait time

    # Map of data points on SKELETON
    map2 = rast.skeleton
    map2[data == 1] = 10
    map2[map2 == 1] = 11
    map2[map2 == 0] = 12

    {
      if (.f.drawMaps | .f.writeGeoTIF) { # Rasterize Obs Points on the Skeleton (map2) and apply the extents & crs
        map2Ras <- raster(map2)
        extent(map2Ras) <- extent(mapRas)
        crs(map2Ras) <- crs(mapRas)
      }

      if (.f.writeGeoTIF | .f.writeCSVs) { # create the output directory
        theDir <- file.path(CLHRDirectory, RS, animalNickname, '2.LocationsOnSkeleton')
        dir.create(theDir, showWarnings = F, recursive = T)
      }

      if (.f.drawMaps) { # plot the Obs Points on the Skeleton
        plot(map2Ras, main = paste(animalNickname,'Locations On Region Skeleton', sep = " "),
             xlab = "Easting", ylab = "Northing", col = c('red', 'blue', 'green'), legend = F)
      }

      if (.f.writeGeoTIF) { # Write Obs Points on the Skeleton to gtif
        rasterToGTIFF.fn(map2Ras, theDir, animalNickname, RS, 'LocationsOnSkeleton')
      }

      if (.f.writeCSVs) { # Write Obs Points on the Skeleton to csv
        matrixToSimpleCSV.fn(map2, theDir, animalNickname, RS, 'LocationsOnSkeleton')
      }
    }

    # Counting points ---------------------------------------------------------

    # Joint matrix and numbering data points -----------------------------------------------
    {
      joints = matrix(0, m, n)
      jointnum = 1
      datanum = 1
      for (i in 1:m) {
        for (j in 1:n) {
          if (skeleton[i, j] > 2) {
            joints[i, j] = jointnum
            jointnum = jointnum + 1
          }
          if (data[i, j] > 0) {
            data[i, j] = datanum
            datanum = datanum + 1
          }
        }
      }
      jointnum = jointnum - 1
      datanum = datanum - 1
      distsize = jointnum + datanum
    }
    # Adding numbered data points and joints to skeleton -----------------------------------------------
    {
      sjd = data
      k = 1
      for (i in 1:m) {
        for (j in 1:n) {
          if (data[i, j] == 0 & joints[i, j] != 0) {
            sjd[i, j] = datanum + k
            k = k + 1
          }
        }
      }
      distsize = max(sjd)
      for (i in 1:m) {
        for (j in 1:n) {
          if (data[i, j] == 0 & joints[i, j] == 0 & skeleton[i, j] != 0) {
            sjd[i, j] = distsize + 1
          }
        }
      }
    }
    # Distance between all points -----------------------------------------------
    {
      travel <- function(skel2, i, j, d) {
        x = 0
        while (x == 0) {
          if ((i > 1) & (j > 1)) {
            if (skel2[i - 1, j - 1] > 0) {
              x = 1
              i = i - 1
              j = j - 1
              d = d + sqrt(2)
              break
            }
          }
          if (i > 1) {
            if (skel2[i - 1, j] > 0) {
              x = 1
              i = i - 1
              d = d + 1
              break
            }
          }
          if ((i > 1) & (j < n)) {
            if (skel2[i - 1, j + 1] > 0) {
              x = 1
              i = i - 1
              j = j + 1
              d = d + sqrt(2)
              break
            }
          }
          if (j < n) {
            if (skel2[i, j + 1] > 0) {
              x = 1
              i = i
              j = j + 1
              d = d + 1
              break
            }
          }
          if ((i < m) & (j < n)) {
            if (skel2[i + 1, j + 1] > 0) {
              x = 1
              i = i + 1
              j = j + 1
              d = d + sqrt(2)
              break
            }
          }
          if (i < m) {
            if (skel2[i + 1, j] > 0) {
              x = 1
              i = i + 1
              j = j
              d = d + 1
              break
            }
          }
          if ((i < m) & (j > 1)) {
            if (skel2[i + 1, j - 1] > 0) {
              x = 1
              i = i + 1
              j = j - 1
              d = d + sqrt(2)
              break
            }
          }
          if (j > 1) {
            if (skel2[i, j - 1] > 0) {
              x = 1
              i = i
              j = j - 1
              d = d + 1
              break
            }
          }
          x = 2
        }
        result = c(i, j, d, x)
        return(result)
      }
    }

    # Calculate distance between pairs of data points or joints
    {
      dist = {}
      skelseg = sjd
      row = 1
      segment = distsize + 2
      for (i in 1:distsize) {
        j = which((sjd == i), arr.ind = TRUE)
        x = j[1]
        y = j[2]
        jsize = skeleton[x, y]
        skel2 = skeleton
        for (t in 1:jsize) {
          check = 0
          d = 0
          x = j[1]
          y = j[2]
          skel2[x, y] = 0
          while (check == 0) { # if any of the following are true, check != 1, advance the loop
            new = travel(skel2, x, y, d)
            x = new[1]
            y = new[2]
            d = new[3]
            z = new[4]
            skel2[x, y] = 0
            if (sjd[x, y] != (distsize + 1)) {
              if (row == 1) {
                dist = c(i, sjd[x, y], d, segment)
                row = row + 1
                segment = segment + 1
                check = 1
              } else if (row == 2) {
                if (dist[1] == i & dist[2] == sjd[x, y]) {
                  if (dist[3] > d) {
                    dist[3] = d
                    row = row + 1
                    segment = segment + 1
                    check = 1
                  } else {
                    segment = segment + 1
                    check = 1
                  }
                } else {
                  new.dist.row = c(i, sjd[x, y], d, segment)
                  dist = rbind(dist, new.dist.row)
                  row = row + 1
                  segment = segment + 1
                  check = 1
                }
              } else if (any(dist[, 1] == i & dist[, 2] == sjd[x, y])) {
                oldrow = which((dist[, 1] == i & dist[, 2] == sjd[x, y]), arr.ind = TRUE)
                old = dist[oldrow, 3]
                if (old > d) {
                  new.dist.row = c(i, sjd[x, y], d, segment)
                  dist = rbind(dist, new.dist.row)
                  row = row + 1
                  segment = segment + 1
                  check = 1
                } else {
                  segment = segment + 1
                  check = 1
                }
              } else {
                new.dist.row = c(i, sjd[x, y], d, segment)
                dist = rbind(dist, new.dist.row)
                row = row + 1
                segment = segment + 1
                check = 1
              }
            } else {
              skelseg[x, y] = segment
            }
            if (z == 2) {
              check = 1
              segment = segment + 1
            }
          }
        }
      }
    }

    # Shortest paths between data points -----------------------------------------------
    dist2 = as.data.frame(dist[, 1:3])

    # Make column names -----------------------------------------------
    names(dist2) = c("from", "to", "distance")

    # Transform into graph -----------------------------------------------
    dgraph = graph.data.frame(dist2, directed = FALSE)

    # Finding shortest paths between data points -----------------------------------------------
    distances = matrix(0, datanum, datanum)
    paths = {}
    for (i in 1:(datanum - 1)) {
      for (j in (i + 1):datanum) {
        path = get.shortest.paths(dgraph, i, j, weights = E(dgraph)$distance)
        new.paths = list(c(i, j), c(unlist(path[1][[1]][1])))
        paths = rbind(paths, new.paths)
        distances[i, j] = shortest.paths(dgraph, i, j, weights = E(dgraph)$distance)
        distances[j, i] = distances[i, j]
      }
    }
    distances2 = distances

    # Minimum spanning tree - Has issues -----------------------------------------------

    # Change distances matrix for Prim's algorithm
    {
      distances3 = distances2
      for (i in 1:(datanum - 1)) {
        for (j in (i + 1):datanum) {
          for (k in 1:dim(paths)[1] ) {
            if (paths[, 1][[k]][1] == i & paths[, 1][[k]][2] == j) {
              for (l in 2:(length(paths[, 2][[k]]) - 1)) {
                if (paths[, 2][[k]][l] <= datanum) {
                  distances3[i, j] = Inf
                  distances3[j, i] = Inf
                }
              }
              break
            }
          }
        }
      }
      dpairs = matrix(0, datanum * (datanum - 1)/2, 3)
      row = 1
      for (i in 1:(datanum - 1)) {
        for (j in (i + 1):datanum) {
          dpairs[row, 1] = i
          dpairs[row, 2] = j
          dpairs[row, 3] = distances2[i, j]
          row = row + 1
        }
      }
    }

    # Finding minimum distances using Prim's algorithm
    dpairsleft = dpairs
    minimum = matrix(0, datanum - 1, 3)
    used = 1
    l = 1
    check = datanum - 1

    # PROBLEM ####
    # PROBLEM ####
    # PROBLEM ####
    tryCatch(while (check != 0) {

      minrow = matrix(NA, dim(dpairsleft)[1], 3)
      col = 1
      for (k in 1:l) {
        for (i in 1:dim(dpairsleft)[1]) {
          for (j in 1:2) {
            if (dpairsleft[i, j] == used[k]) {
              minrow[col, ] = dpairsleft[i, ]
              col = col + 1
            }
          }
        }
      }
      minimum[l, 3] = min(minrow[, 3], na.rm = TRUE)
      minimum[l, 1] = minrow[which(minrow[, 3] == minimum[l, 3])[1], 1]
      minimum[l, 2] = minrow[which(minrow[, 3] == minimum[l, 3])[1], 2]
      l = l + 1
      used[l] = minimum[l - 1, 1]
      if (any(duplicated(used))) {
        used[l] = minimum[l - 1, 2]
      }
      for (i in seq(dim(dpairsleft)[1], 1, by = -1)) {
        if (any(dpairsleft[i, 1] == used) & any(dpairsleft[i, 2] == used)) {
          dpairsleft = dpairsleft[-i, ]
        }
      }
      check = check - 1
    },
    # error = function(e){cat("ERROR :",conditionMessage(e), "\n")}) # Error in dpairsleft[i, 1] : incorrect number of dimensions
    error = function(e){}) # Error in dpairsleft[i, 1] : incorrect number of dimensions

    minpathways = {}
    for (i in 1:(datanum - 1)) {
      for (j in 1:dim(paths)[1]) {
        if (paths[, 1][[j]][1] == minimum[i, 1] & paths[, 1][[j]][2] == minimum[i, 2]) {
          row = j
        }
      }
      for (k in 1:(length(paths[, 2][[row]]) - 1)) {
        new.minpath = c(paths[, 2][[row]][k], paths[, 2][[row]][k + 1])
        minpathways = rbind(minpathways, new.minpath)
      }
    }

    # Getting rid of duplicated rows -----------------------------------------------
    # maybe simplify() from igraph will do this
    minpathways2 = minpathways
    for (i in 1:dim(minpathways)[1]) {
      if (minpathways[i, 1] > minpathways[i, 2]) {
        minpathways2[i, 1] = minpathways[i, 2]
        minpathways2[i, 2] = minpathways[i, 1]
      }
    }
    minpathways3 = minpathways2[!duplicated(minpathways2), ]

    minpathways3 = cbind(minpathways3, matrix(0, dim(minpathways3)[1], 2))
    for (i in 1:dim(minpathways3)[1]) {
      row = which(minpathways3[i, 1] == dist[, 2] & minpathways3[i, 2] == dist[, 1])
      minpathways3[i, 3] = dist[row, 3]
      minpathways3[i, 4] = dist[row, 4]
    }

    mstdistance = sum(minpathways3[, 3])

    mst = matrix(0, m, n)
    for (i in 1:m) {
      for (j in 1:n) {
        if (sjd[i, j] <= datanum & sjd[i, j] != 0) {
          mst[i, j] = 1
        } else {
          for (k in 1:dim(minpathways3)[1]) {
            if (skelseg[i, j] == minpathways3[k, 4]) {
              mst[i, j] = 1
            } else if (sjd[i, j] == minpathways3[k, 1]) {
              mst[i, j] = 1
            } else if (sjd[i, j] == minpathways3[k, 2]) {
              mst[i, j] = 1
            }
          }
        }
      }
    } # medium wait time

    # Create map of MST and data points on skeleton -----------------------------------------------
    map3 = rast.skeleton
    map3[mst == 1] = 7
    map3[data != 0] = 5
    map3[map3 == 1] = 10
    map3[map3 == 0] = 12

    {
      if (.f.drawMaps | .f.writeGeoTIF) { # Rasterize MST and data points on skeleton (map3) and apply the extents & crs
        map3Ras <- raster(map3)
        extent(map3Ras) <- extent(map2Ras)
        crs(map3Ras) <- crs(map2Ras)
      }

      if (.f.writeGeoTIF | .f.writeCSVs) { # create the output directory
        theDir <- file.path(CLHRDirectory, RS, animalNickname, '3.MinimumSpanningTree')
        dir.create(theDir, showWarnings = F, recursive = T)
      }

      if (.f.drawMaps) { # plot the MST and data points on skeleton
        plot(map3Ras, main = paste(animalNickname, 'Minimum Spanning Tree', sep = " "),
             xlab = "Easting", ylab = "Northing", col = c('red', 'pink', 'blue', 'green'), legend = F)
      }

      if (.f.writeGeoTIF) { # Write MST and data points on skeleton to gtif
        rasterToGTIFF.fn(map3Ras, theDir, animalNickname, RS, 'MST')
      }

      if (.f.writeCSVs) { # Write MST and data points on skeleton to csv
        matrixToSimpleCSV.fn(map3, theDir, animalNickname, RS, 'MST')
      }
      # Output distance of MST
      if (.f.writeResultText) {
      dir.create(file.path(CLHRDirectory, RS, animalNickname), showWarnings = F, recursive = T)
      sink(paste(file.path(CLHRDirectory, RS, animalNickname), '/', animalNickname, "-MST.txt", sep = ''))
      cat('The Minimum Spanning Tree is:', mstdistance, '\n')
      sink(file = paste(file.path(CLHRDirectory, RS), '/', RS, "-Region-MSTs.csv", sep = ''), append = T, split = T)
      cat(c(animalNickname, ",", mstdistance, ",",mstdistance*as.numeric(RS),"\n"))
      } else {
        sink(file = paste(file.path(CLHRDirectory, RS), '/', RS, "-Region-MSTs.csv", sep = ''), append = T)
        cat(c(animalNickname, ",", mstdistance, ",",mstdistance*as.numeric(RS),"\n"))
      }
      sink.reset()
      cat('The MST distance for', animalNickname, 'is:', mstdistance*as.numeric(RS), 'm\n')
    }
  }
}

