prepASCtoMatrix.fn <- function(tileFile) {
  # x <- locationFilename
  x <- as.matrix(read.table(file = tileFile, header = F, skip = 6, sep = ' '))
  x[x == -9999] = 0
  x <- as.data.table(x) # was recommended to reduce system load for large data frames
  x <- x[, which(unlist(lapply(x, function(x)!all(is.na(x))))), with = F]
  x <- as.matrix(x)
  binarise(x)
}

extractMajorSkeleton.fn <- function(x, kernel = shapeKernel(c(3,3), type = 'box')) {
  chunks <- components(x, kernel)
  # sort(table(chunks), decreasing = T)
  chunksMaxVal <- as.integer(row.names(as.matrix(sort(table(chunks), decreasing = T))))[1]
  chunksMajComponent <- (chunks == chunksMaxVal)*1
  chunksMajComponent[is.na(chunksMajComponent)] <- 0
  skeletonise(x = chunksMajComponent, method = "hitormiss")
}

copySourceGeoSpatToRas.fn <- function(targetMat, source) {
  rasTarget <- raster(targetMat)
  extent(rasTarget) <- readGDAL(source, silent = T)
  crs(rasTarget) <- crs(readGDAL(source, silent = T))
  return(rasTarget)
}

rasterToGTIFF.fn <- function(x, dir = getwd(), indiv = NULL, suffix = NULL) {
  SPDF <- as(x, 'SpatialPixelsDataFrame')
  writeGDAL(SPDF, fname = str_c(dir, '/', indiv, '-', suffix, '.tif', sep = ''), drivername = 'GTiff', type = 'Byte', mvFlag = 0, options = 'TFW=YES')
}

matrixToSimpleCSV.fn <- function(x, dir = getwd(), indiv = NULL, suffix = NULL) {
  write.table(x, file = str_c(dir, '/', indiv, '-', suffix, '.csv', sep = ''), append = F, sep = ',', row.names = F, col.names = F)
}

writeOutAsc.fn <- function(x, dir = getwd(), indiv = NULL, suffix = NULL) {
  writeRaster(x, filename = str_c(dir, '/', indiv, '-', suffix, '.asc', collapse = T), format = 'ascii', datatype = 'INT1U', overwrite = T, prj = T)
}

sink.reset <- function() { # courtsey of Dason from Wesley Chapel, Florida via Stackoverflow.com
  for (i in seq_len(sink.number())) {
    sink(NULL)
  }
}

connectivityArray.fn <- function(x) {
  m <- dim(x)[1] # increment majorSkeleton values if they are non-zero, based on the number of non-zero 8-way neighbours
  n <- dim(x)[2]   # this creates a map of 'connectedness'
  for (i in 1:m) {
    for (j in 1:n) {
      if (x[i, j] > 0) { # if the current points is not zero
        joint = 0
        if ((i > 1) & (j > 1)) { # any point excluding the first row and column
          if (x[i - 1, j - 1] > 0) { # compare with up-left
            joint = joint + 1
          }
        }
        if (i > 1) { # any column excluding the first row
          if (x[i - 1, j] > 0) { # compare with up
            joint = joint + 1
          }
        }
        if ((i > 1) & (j < n)) { # any point excluding the first row and excluding the last column
          if (x[i - 1, j + 1] > 0) { # compare with up-right
            joint = joint + 1
          }
        }
        if (j < n) { # any row excluding the last column
          if (x[i, j + 1] > 0) { # compare with right
            joint = joint + 1
          }
        }
        if ((i < m) & (j < n)) { # any point excluding the last row and column
          if (x[i + 1, j + 1] > 0) { # compare with down-right
            joint = joint + 1
          }
        }
        if (i < m) { # any column excluding the last
          if (x[i + 1, j] > 0) { # compare with down
            joint = joint + 1
          }
        }
        if ((i < m) & (j > 1)) { # any point excluding the last row and any point excluding the first column
          if (x[i + 1, j - 1] > 0) { # compare with down-left
            joint = joint + 1
          }
        }
        if (j > 1) { # any column excluding the first
          if (x[i, j - 1] > 0) { # compare with left
            joint = joint + 1
          }
        }
        if (joint > 0) {
          x[i, j] = joint # increment the current point by 1 for every non-zero 8-way neighbour
        }
      }
    }
  }
  return(x)
}

shiftDataToSkeleton.fn <- function(points, skeleton) {
  m = dim(points)[1]
  n = dim(points)[2]
  k = dim(skeleton)[1]
  l = dim(skeleton)[2]
  data = matrix(0, m, n)
  # shift = data.frame(matrix(0, m, n), fix.empty.names = T)
  for (i in 1:m) {
    for (j in 1:n) { # move down each row
      if (points[i, j] == 1) { # move across each column in points
        minDist = Inf
        for (k in 1:m) { # move down each row
          for (l in 1:n) { # move across each column in majorSheleton
            if (skeleton[k, l] == 1) {
              dist = (i - k)^2 + (j - l)^2 # starting with a value at infinity compare the new value to the old, if smaller
              if (dist < minDist) { # take the smaller value as minDist
                minDist = dist
                miny = k
                minx = l
              }
            }
          }
        }
        data[miny, minx] = 1
        # shift[i, j] = list(str_c(miny, minx, minDist, sep = ","))
      }
    }
  }
  # result = c(data, shift)
  shiftedPoints <<- data
  mapIndivsOnLandSkel = skeleton
  mapIndivsOnLandSkel[data == 1] = 10
  mapIndivsOnLandSkel[mapIndivsOnLandSkel == 1] = 11
  mapIndivsOnLandSkel[mapIndivsOnLandSkel == 0] = 12
  # return(data)
  return(mapIndivsOnLandSkel)
} # The cells on the data map get values that increase from top left to bottom right

makeJointsNPointsSequential.fn <- function(connectivityMap, pointsOnSkeletonMap) {
  m = dim(connectivityMap)[1]
  n = dim(connectivityMap)[2]

  lastJointNum = 1
  lastDataNum = 1
  joints = matrix(0, m, n)
  for (i in 1:m) { # the nested for loop in each of these steps creates a consistent orientation to the generation
    for (j in 1:n) { # of the points top left to botom right
      if (connectivityMap[i, j] > 2) { # if the connectivity map has a value >2 it's a joint location
        joints[i, j] = lastJointNum # mark that value into the the "joints" map at that location
        lastJointNum = lastJointNum + 1 # and increment the lastJointNum by one (these are jointIDs)
      }
      if (pointsOnSkeletonMap[i, j] > 0) {
        pointsOnSkeletonMap[i, j] = lastDataNum
        lastDataNum = lastDataNum + 1
      }
    }
  }
  lastJointNum = lastJointNum - 1 # reduce by one so lastJointNum has the value of the number of joints
  lastDataNum = lastDataNum - 1
  distSize = lastJointNum + lastDataNum # distSize is all of the joints plus all of the points...
  results = list(aa = joints, bb = lastDataNum, cc = distSize, dd = pointsOnSkeletonMap)
  return(results)
}

travel.fn <- function(con, i, j, d) { # Distance between all points
  x = 0
  while (x == 0) {
    if ((i > 1) & (j > 1)) {
      if (con[i - 1, j - 1] > 0) {
        x = 1
        i = i - 1
        j = j - 1
        d = d + sqrt(2)
        break
      }
    }
    if (i > 1) {
      if (con[i - 1, j] > 0) {
        x = 1
        i = i - 1
        d = d + 1
        break
      }
    }
    if ((i > 1) & (j < n)) {
      if (con[i - 1, j + 1] > 0) {
        x = 1
        i = i - 1
        j = j + 1
        d = d + sqrt(2)
        break
      }
    }
    if (j < n) {
      if (con[i, j + 1] > 0) {
        x = 1
        i = i
        j = j + 1
        d = d + 1
        break
      }
    }
    if ((i < m) & (j < n)) {
      if (con[i + 1, j + 1] > 0) {
        x = 1
        i = i + 1
        j = j + 1
        d = d + sqrt(2)
        break
      }
    }
    if (i < m) {
      if (con[i + 1, j] > 0) {
        x = 1
        i = i + 1
        j = j
        d = d + 1
        break
      }
    }
    if ((i < m) & (j > 1)) {
      if (con[i + 1, j - 1] > 0) {
        x = 1
        i = i + 1
        j = j - 1
        d = d + sqrt(2)
        break
      }
    }
    if (j > 1) {
      if (con[i, j - 1] > 0) {
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

# the monster ####

calculateTheMinimumSpanningTree.fn <- function(sPoints = shiftedPoints, manyJoints = joints,
                                               datanum = seqDataNum, distsize = distanceSize,
                                               seqShiftedPoints = seqShiftedPts,
                                               connectivityLayout = baseConnectivity) {
  m = dim(sPoints)[1]
  n = dim(sPoints)[2]

  sjd = seqShiftedPoints
  k = 1
  for (i in 1:m) {
    for (j in 1:n) {
      if (seqShiftedPoints[i, j] == 0 & manyJoints[i, j] != 0) {
        sjd[i, j] = datanum + k
        k = k + 1
      }
    }
  }

  distsize = max(sjd) # the largest value from "where there is not indiv but there is a joint"
  for (i in 1:m) {
    for (j in 1:n) {
      # if it's not a data point (coincident indiv and skel) AND it's not a joint AND it's on the skeleton
      if (seqShiftedPoints[i, j] == 0 & manyJoints[i, j] == 0 & connectivityLayout[i, j] != 0) {
        sjd[i, j] = distsize + 1 # mark these points with that value + 1
      } # if a point overlaps a joint, then that spot is not counted among the joints.
    } # all non-point, non-joint positions on the skeleton now have the same value.
  }

  # Calculate distance between pairs of data points or joints
  {
    dist = {}
    skeletonSegment = sjd # new layer
    row = 1
    segment = distsize + 2 # no segment is larger than this
    for (i in 1:distsize) {
      j = which((sjd == i), arr.ind = TRUE)
      x = j[1]
      y = j[2]
      jsize = connectivityLayout[x, y]
      skel2 = connectivityLayout
      for (t in 1:jsize) {
        check = 0
        d = 0
        x = j[1]
        y = j[2]
        skel2[x, y] = 0
        while (check == 0) { # if any of the following are true, check != 1
          new = travel.fn(skel2, x, y, d)
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
            skeletonSegment[x, y] = segment
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
  dist2 = as.data.frame(dist[, 1:3]) # subset of distances between points

  # Make column names -----------------------------------------------
  names(dist2) = c("from", "to", "distance")

  # Transform into graph -----------------------------------------------
  dgraph = graph.data.frame(dist2, directed = FALSE)

  # Finding shortest paths between data points -----------------------------------------------
  distances = matrix(0, datanum, datanum)
  pathBucket = {}
  for (i in 1:(datanum - 1)) {
    for (j in (i + 1):datanum) {
      path = get.shortest.paths(dgraph, i, j, weights = E(dgraph)$distance)
      new.paths = list(c(i, j), c(unlist(path[1][[1]][1])))
      pathBucket = rbind(pathBucket, new.paths)
      distances[i, j] = shortest.paths(dgraph, i, j, weights = E(dgraph)$distance)
      distances[j, i] = distances[i, j]
    }
  }
  distances2 = distances

  # Minimum spanning tree -----------------------------------------------

  # Change distances matrix for Prim's algorithm
  {
    distances3 = distances2
    for (i in 1:(datanum - 1)) {
      for (j in (i + 1):datanum) {
        for (k in 1:dim(pathBucket)[1] ) {
          if (pathBucket[, 1][[k]][1] == i & pathBucket[, 1][[k]][2] == j) {
            for (l in 2:(length(pathBucket[, 2][[k]]) - 1)) {
              if (pathBucket[, 2][[k]][l] <= datanum) {
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
    for (i in seq(dim(dpairsleft)[1], 1, by = -1)) { # unsure of this
      if (any(dpairsleft[i, 1] == used) & any(dpairsleft[i, 2] == used)) {
        dpairsleft = dpairsleft[-i, ]
      }
    }
    check = check - 1
  }, error = function(e){}) # silenced # Error in dpairsleft[i, 1] : incorrect number of dimensions

  minpathways = {}
  for (i in 1:(datanum - 1)) {
    for (j in 1:dim(pathBucket)[1]) {
      if (pathBucket[, 1][[j]][1] == minimum[i, 1] & pathBucket[, 1][[j]][2] == minimum[i, 2]) {
        row = j
      }
    }
    for (k in 1:(length(pathBucket[, 2][[row]]) - 1)) {
      new.minpath = c(pathBucket[, 2][[row]][k], pathBucket[, 2][[row]][k + 1])
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

  mstdistance <<- sum(minpathways3[, 3]) # The MST !!!

  mst = matrix(0, m, n) # The MST illustrated !!!
  for (i in 1:m) {
    for (j in 1:n) {
      if (sjd[i, j] <= datanum & sjd[i, j] != 0) {
        mst[i, j] = 1
      } else {
        for (k in 1:dim(minpathways3)[1]) {
          if (skeletonSegment[i, j] == minpathways3[k, 4]) {
            mst[i, j] = 1
          } else if (sjd[i, j] == minpathways3[k, 1]) {
            mst[i, j] = 1
          } else if (sjd[i, j] == minpathways3[k, 2]) {
            mst[i, j] = 1
          }
        }
      }
    }
  } #$$$ 9 Essential Arguments; 37 Variables; 3 Outputs

  mstMap <<- mst # Just the MST line

  # Create mapIndivsOnLandFull of MST and data points on skeleton ####
  indivsOnMST = majorSkeleton
  indivsOnMST[mst == 1] = 7
  indivsOnMST[seqShiftedPoints != 0] = 5
  indivsOnMST[indivsOnMST == 1] = 10
  indivsOnMST[indivsOnMST == 0] = 12
  return(indivsOnMST) # Returns the Image
}
