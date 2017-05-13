# CLHR R code, version 3.5

# Change History ############################################

# 1. JC saved file with new name and gave it version 2.0
# 2. Eli removed the need for GRASS with the help of Jon Clayden, writes CSVs V3.3
# 3. Write GeoTiffs V3.4
# 4. Loop structure in place

# To-do: REMOVE ALL LAKES
# mmand connected components
# feed in a subset of landscape and points to see if it runs without error
# if it doesn't find a problem in the subset then it's a connectivity thing

# To-do: correct distance discrepency
# To-do: multiple individuals / landscapes
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
library(igraph)
library(gtools)
library(stringr)
library(maptools)
# temp load dev version of mmand
# devtools::install_github("jonclayden/mmand")
library(mmand)
dev.off()
rm(list = ls())

# Configuration ####
{
.f.writeCSVs        <- T
.f.writeGeoTIF      <- F
.f.drawMaps         <- T
.f.writeResultText  <- T
shoreFolder         <- '1.Shore'
observationFolder   <- '2.Turtles'
outputDirectory     <- '3.CLHRs'
allRegionFolders    <- mixedsort(dir(shoreFolder, all.files = F, full.names = F))
allAnimalFolders    <- mixedsort(dir(observationFolder, all.files = F, full.names = F))
animalPrefix        <- ''
animalSuffix        <- 'output'
locationPrefix      <- ''
locationSuffix      <- ''
resolutions         <- c(100)
nresolutions        <- length(resolutions)
}

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
}

##### THE LARGE LOOP #####
for (pRS in seq_along(resolutions)) {
  pRS <- 1
  RS <- resolutions[pRS]
  animsByRes <- grep(str_c('\\S', RS, '\\D', sep = ''), allAnimalFolders, value = T, perl = T)
  currentFolder <- grep(str_c('\\S', RS, '\\D', sep = ''), allRegionFolders, value = T, perl = T)
  locationFilepath <- file.path(shoreFolder, currentFolder)
  locationFilename <- str_c(locationFilepath, '/', currentFolder, '.txt', sep = '')

  # cat(RS, 'pass', currentFolder,'\n')

  for (pAS in seq_along(animsByRes)) {
    pAS <- 5
    AS <- animsByRes[pAS]
    truncAni <- str_locate(pattern = as.character(RS), string = AS)
    animalNickname <- (str_sub(AS, start = (truncAni[1] - truncAni[1] + 1), end = truncAni[2]))
    dataFPPrep <- paste(animalNickname, animalSuffix, sep = '')
    dataFilepath <- file.path(observationFolder, dataFPPrep) # species location info
    dataFilename <- paste(basename(dataFilepath), '.txt', sep = '')

    cat(locationFilepath, locationFilename, animalNickname, dataFilepath, dataFilename, '\n')

    # (dist, dist2, distances, distances2, distances3, dpairs, dpairsleft, joints, map, map2, map3, minimum, minpathways, minpathways2, minpathways3, minrow, mst, paths, sjd, skel2, skeleton, skelseg, check, col, D, d, datanum, dgraph, distsize, i, j, joint, jointnum, jsize, k, l, m, minD, minx, miny, mstdistance, n, new, new.dist.row, new.minpath, new.paths, path, row, segment, t, x, y, z)

    # Import location ----------------------------------------------------------------
    tile = as.matrix(read.table(locationFilename, header = F, skip = 6, sep = " "))

    matrixtile <- prepASCtoMatrix.fn(tile) # this function is for the niagara turtles

    m <- dim(matrixtile)[1]
    n <- dim(matrixtile)[2]

    # Rasterize location (matrixtile) and retrieve its extents & crs
    if (.f.drawMaps | .f.writeGeoTIF) {
      rastertile <- raster(matrixtile)
      extent(rastertile) <- readGDAL(locationFilename)
      crs(rastertile) <- crs(readGDAL(locationFilename)) # readGDAL reads the accompanying .prj file
    }

    # create the output directory
    if (.f.writeGeoTIF | .f.writeCSVs) {
      theDir <- file.path(outputDirectory, RS, animalNickname, shoreFolder)
      dir.create(theDir, showWarnings = F, recursive = T)
    }

    # plot the location
    if (.f.drawMaps) {
      plot(rastertile, col = c('green', 'blue'), main = 'Landscape', xlab = "Easting", ylab = "Northing", legend = F)
    }

    # Write location to gtif
    if (.f.writeGeoTIF) {
      rasterToGTIFF.fn(rastertile, theDir, animalNickname, RS, 'locale')
    }

    # Write location to csv
    if (.f.writeCSVs) {
      matrixToSimpleCSV.fn(matrixtile, theDir, animalNickname, RS, 'locale')
    }

    # Eroding location matrix to skeleton ----------------------------------------------
    # rast.skeleton <- skeletonise3_2d.fn(matrixtile)
    p.test <- components(matrixtile, shapeKernel(c(3,3), type = "box"))
    sort(table(p.test), decreasing = T)
    p.testMaxVal <- as.integer(row.names(as.matrix(sort(table(p.test), decreasing = T))))[1]
    p.testLargestComponent <- (p.test == p.testMaxVal)*1
    p.testLargestComponent[is.na(p.testLargestComponent)] <- 0
    rast.skeleton <- skeletonise3_2d.fn(p.testLargestComponent)
    # perhaps igraph here can give us a list of connected elements and we then discard the shortest ones and continue on.
    # raster::clump might help
    # patchify might help
    # rasterToPolygons(..., dissolve=T, fun=function(x){x>6}) -> disaggregate -> slotNames
      # p.test <- components(matrixtile, shapeKernel(c(3,3), type="box"))
      # p.testR <- raster(p.test)
      # p.testP <- rasterToPolygons(p.testR, fun = function(x){x>0}, dissolve = T)
      # crs(p.testP) <- crs(readGDAL(locationFilename))
      # p.disagg <- disaggregate(p.testP)
      # p.disagg$layer <- factor(seq_len(length(p.disagg)))
    ###########################################
    # Rasterize skeleton (rast.skeleton) and apply the extents & crs
    if (.f.drawMaps | .f.writeGeoTIF) {
      skelRas <- raster(rast.skeleton)
      extent(skelRas) <- extent(rastertile)
      crs(skelRas) <- crs(rastertile)
    }

    # create the output directory
    if (.f.writeGeoTIF | .f.writeCSVs) {
      theDir <- file.path(outputDirectory, RS, animalNickname, '3.Skeleton')
      dir.create(theDir, showWarnings = F, recursive = T)
    }

    # plot the skeleton
    if (.f.drawMaps) {
      plot(skelRas, col = c('green', 'blue'), main = 'Landscape Skeleton', xlab = "Easting", ylab = "Northing", legend = F)
    }

    # Write skeleton to gtif
    if (.f.writeGeoTIF) {
      rasterToGTIFF.fn(skelRas, theDir, animalNickname, RS, 'skeleton')
    }

    # Write skeleton to csv
    if (.f.writeCSVs) {
      matrixToSimpleCSV.fn(rast.skeleton, theDir, animalNickname, RS, 'skeleton')
    }

    # Counting joints and counting skeleton pixels
    skeleton = rast.skeleton
    for (i in 1:m) {
      for (j in 1:n) {
        if (skeleton[i, j] > 0) {
          joint = 0
          if ((i > 1) & (j > 1)) {
            if (skeleton[i - 1, j - 1] > 0) {
              joint = joint + 1
            }
          }
          if (i > 1) {
            if (skeleton[i - 1, j] > 0) {
              joint = joint + 1
            }
          }
          if ((i > 1) & (j < n)) {
            if (skeleton[i - 1, j + 1] > 0) {
              joint = joint + 1
            }
          }
          if (j < n) {
            if (skeleton[i, j + 1] > 0) {
              joint = joint + 1
            }
          }
          if ((i < m) & (j < n)) {
            if (skeleton[i + 1, j + 1] > 0) {
              joint = joint + 1
            }
          }
          if (i < m) {
            if (skeleton[i + 1, j] > 0) {
              joint = joint + 1
            }
          }
          if ((i < m) & (j > 1)) {
            if (skeleton[i + 1, j - 1] > 0) {
              joint = joint + 1
            }

          }
          if (j > 1) {
            if (skeleton[i, j - 1] > 0) {
              joint = joint + 1
            }
          }
          if (joint > 0) {
            skeleton[i, j] = joint
          }
        }
      }
    }


    # Import observations -------------------------------------------------
    sdata <- as.matrix(read.table(file.path(dataFilepath, dataFilename), header = F, skip = 6, sep = " "))
    sdata <- prepASCtoMatrix.fn(sdata)

    # Rasterize the observation points (sdata) and apply the extents & crs
    if (.f.drawMaps | .f.writeGeoTIF) {
      sdataRas <- raster(sdata)
      extent(sdataRas) <- extent(skelRas)
      crs(sdataRas) <- crs(skelRas)
    }

    # create the output directory
    if (.f.writeGeoTIF | .f.writeCSVs) {
      theDir <- file.path(outputDirectory, RS, animalNickname, '2.ObservedAnimalPoints')
      dir.create(theDir, showWarnings = F, recursive = T)
    }

    # plot the the observation points
    if (.f.drawMaps) {
      plot(sdataRas, main = 'Observation Points', xlab = "Easting", ylab = "Northing", col = c('white', 'red'), legend = F)
    }

    # Write the observation points to gtif
    if (.f.writeGeoTIF) {
      rasterToGTIFF.fn(sdataRas, theDir, animalNickname, RS, 'ObsPtsOnly')
    }

    # Write the observation points to csv
    if (.f.writeCSVs) {
      matrixToSimpleCSV.fn(sdata, theDir, animalNickname, RS, 'ObsPtsOnly')
    }

    # Observations on landscape Map to be exported
    map = matrixtile
    map[sdata == 1] = 10
    map[map == 1] = 11
    map[map == 0] = 12

    # Rasterize Obs Points on the landscape (map) and apply the extents & crs
    if (.f.drawMaps | .f.writeGeoTIF) {
      mapRas <- raster(map)
      extent(mapRas) <- extent(sdataRas)
      crs(mapRas) <- crs(sdataRas)
    }

    # create the output directory
    if (.f.writeGeoTIF | .f.writeCSVs) {
      theDir <- file.path(outputDirectory, RS, animalNickname, '2.ObservedAnimalPoints')
      dir.create(theDir, showWarnings = F, recursive = T)
    }

    # plot the Obs Points on the Landscape
    if (.f.drawMaps) {
      plot(mapRas, main = 'Observation Points On Landscape', xlab = "Easting", ylab = "Northing", col = c('red', 'blue', 'green'), legend = F)
    }

    # Write Obs Points on the Landscape to gtif
    if (.f.writeGeoTIF) {
      rasterToGTIFF.fn(mapRas, theDir, animalNickname, RS, 'ObsPointsOnLandscape')
    }

    # Write Obs Points on the Landscape to csv
    if (.f.writeCSVs) {
      matrixToSimpleCSV.fn(map, theDir, animalNickname, RS, 'ObsPointsOnLandscape')
    }

    #### PROBLEM ####
    #### PROBLEM ####
    #### PROBLEM ####
    # Moving data points to skeleton ####
    data = matrix(0, m, n)
    for (i in 1:m) {
      for (j in 1:n) {
        if (sdata[i, j] == 1) {
          minD = Inf
          for (k in 1:m) {
            for (l in 1:n) {
              if (rast.skeleton[k, l] == 1) {
                D = (i - k)^2 + (j - l)^2
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
    # Can shift some points to bodies of water that are not connected to the skeleton.

    # Map of data points on skeleton
    map2 = rast.skeleton
    map2[data == 1] = 10
    map2[map2 == 1] = 11
    map2[map2 == 0] = 12

    # Rasterize Obs Points on the Skeleton (map2) and apply the extents & crs
    if (.f.drawMaps | .f.writeGeoTIF) {
      map2Ras <- raster(map2)
      extent(map2Ras) <- extent(mapRas)
      crs(map2Ras) <- crs(mapRas)
    }

    # create the output directory
    if (.f.writeGeoTIF | .f.writeCSVs) {
      theDir <- file.path(outputDirectory, RS, animalNickname, '4.PointsOnSkeleton')
      dir.create(theDir, showWarnings = F, recursive = T)
    }

    # plot the Obs Points on the Skeleton
    if (.f.drawMaps) {
      plot(map2Ras, main = 'Observation Points On Skeleton', xlab = "Easting", ylab = "Northing", col = c('red', 'blue', 'green'), legend = F)
    }

    # Write Obs Points on the Skeleton to gtif
    if (.f.writeGeoTIF) {
      rasterToGTIFF.fn(map2Ras, theDir, animalNickname, RS, 'PointsOnSkel')
    }

    # Write Obs Points on the Skeleton to csv
    if (.f.writeCSVs) {
      matrixToSimpleCSV.fn(map2, theDir, animalNickname, RS, 'PointsOnSkel')
    }

    # Counting points ---------------------------------------------------------

    # Joint matrix and numbering data points -----------------------------------------------
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

    # Adding numbered data points and joints to skeleton -----------------------------------------------
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

    # Distance between all points -----------------------------------------------
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
        if ((i<m) & (j < n)) {
          if (skel2[i + 1, j + 1] > 0) {
            x = 1
            i = i + 1
            j = j + 1
            d = d + sqrt(2)
            break
          }
        }
        if (i<m) {
          if (skel2[i + 1, j] > 0) {
            x = 1
            i = i + 1
            j = j
            d = d + 1
            break
          }
        }
        if ((i<m)&(j > 1)) {
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

    # Calculate distance between pairs of data points or joints
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
        while (check == 0) {
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
    } # for at least HK13 some observations move to isolated portions of water and cannot be connected
    distances2 = distances

    # Minimum spanning tree - Has issues -----------------------------------------------

    # Change distances matrix for Prim's algorithm
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
    error = function(e){cat("ERROR :",conditionMessage(e), "\n")}) # Error in dpairsleft[i, 1] : incorrect number of dimensions




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

    # Create map of minimum spanning tree and data points on skeleton -----------------------------------------------
    map3 = rast.skeleton
    map3[mst == 1] = 7
    map3[data != 0] = 5
    map3[map3 == 1] = 10
    map3[map3 == 0] = 12

    # Rasterize minimum spanning tree and data points on skeleton (map3) and apply the extents & crs
    if (.f.drawMaps | .f.writeGeoTIF) {
      map3Ras <- raster(map3)
      extent(map3Ras) <- extent(map2Ras)
      crs(map3Ras) <- crs(map2Ras)
    }

    # create the output directory
    if (.f.writeGeoTIF | .f.writeCSVs) {
      theDir <- file.path(outputDirectory, RS, animalNickname, '5.MinimumSpanningTree')
      dir.create(theDir, showWarnings = F, recursive = T)
    }

    # plot the minimum spanning tree and data points on skeleton
    if (.f.drawMaps) {
      plot(map3Ras, main = 'Minimum Spanning Tree', xlab = "Easting", ylab = "Northing", col = c('red', 'pink', 'blue', 'green'), legend = F)
    }

    # Write minimum spanning tree and data points on skeleton to gtif
    if (.f.writeGeoTIF) {
      rasterToGTIFF.fn(map3Ras, theDir, animalNickname, RS, 'MST')
    }

    # Write minimum spanning tree and data points on skeleton to csv
    if (.f.writeCSVs) {
      matrixToSimpleCSV.fn(map3, theDir, animalNickname, RS, 'MST')
    }

    # Output distance of minimum spanning tree
    if (.f.writeResultText) {
      dir.create(file.path(outputDirectory, RS, animalNickname), showWarnings = F, recursive = T)
      sink(paste(file.path(outputDirectory, RS, animalNickname), '/', animalNickname, "-MST-dist.txt", sep = ''))
      cat('The minimum spanning tree is:\n')
      cat(c(animalNickname, mstdistance))
      sink()
    }
    cat('The minimum spanning tree is:', mstdistance)



  }
}

