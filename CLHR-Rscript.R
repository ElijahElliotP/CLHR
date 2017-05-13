tempAscii <- read.asciigrid(paste(filepathLocation,'.txt', sep = ''),as.image = T)
tempAscii <- read.asciigrid(paste(filepathLocation,'.txt', sep = ''),as.image = T)
tempAscii <- read.asciigrid(paste(filepathLocation,'.txt', sep = ''),as.image = T)
# Loading location files --------------------------------------------------------

# CLHR R code, version 3.0



#############################################################
# Change History:
# 1. JC saved file with new name and gave it version 2.0
# 2.
## END Change History
#############################################################

# Setting working directory
# set the working directory to the same directory where this script is found.
# All paths in the parameter file should then be relative to that location. For example,
# "../1.Shore/" accessed that folder.

library(mmand)
rm(list = ls())

CLHRparameters = read.table("CLHR-parameters.txt", as.is = T)
CLHRparameters
landscapefile.fname = CLHRparameters[1,1]
landscapefile.fname
animalnickname = CLHRparameters[2,1]
animalnickname
observationpoints.fname = CLHRparameters[3,1]
observationpoints.fname

outputdirectory = "../3.CLHRs/"
shorelinefolder = "1.Shoreline"
observationfolder = "2.ObservedAnimalPoints/"

afn.makedirectorystructure = function(outputdirectory,animalnickname, shorelinefolder)
{
  animalfullpath = paste(outputdirectory, animalnickname, sep = "")
  animalfullpath
  # thecall = paste("mkdir", animalfullpath)
  thecall = animalfullpath
  thecall
  dir.create(thecall, recursive = T)

  shorelinefullpath = paste(animalfullpath,"/",shorelinefolder,sep = "")
  shorelinefullpath
  # thecall = paste("mkdir", shorelinefullpath)
  thecall = shorelinefullpath
  thecall
  dir.create(thecall, recursive = T)

}

afn.makedirectorystructure(outputdirectory,animalnickname,shorelinefolder)


## Location
# Enter file name of location
# filepath.location = readline("Please enter the file name of the location:") # have them enter the full path
filepath.location = file.path('..',shorelinefolder,'5.Niagara-100m-output','niagara100moutput.txt')
#niagara100moutput.txt
# Enter number of rows in file
# numrowstoread = readline("Please enter the number of rows in the location file:") # instead, read this from the initial lines of the file that they specified
numrowstoread = 239
#239
# Uploading file
# tile = as.matrix(read.table(filepath.location, header=F, skip=6, sep=" ", nrows=as.integer(numrowstoread)))
tile = as.matrix(read.table(filepath.location, header = F, skip = 6, sep = " "))
tile[tile == -9999] = 0
tile[is.na(tile)] = 0


# to reorient properly
# tile2 = t(apply(tile,2,rev))
# tile3 = tile2[-c(2241),]  # hack that should be removed most likely

# Find dimensions
# m = dim(tile3)[1]
# n = dim(tile3)[2]
m <- dim(tile)[1]
n <- dim(tile)[2]
diff <- m - n


# Making a square matrix, because GRASS, which is used later, requires that it operate on a square matrix
# tile4 = tile3
# # diff = m-n
# # x = diff
# if (diff>0){
#   while (x!=0){
#     tile4 = cbind(tile4,0)
#     x = x-1
#   }
# } else if (diff<0){
#   while (x!=0){
#     tile4 = rbind(tile4,0)
#     x = x+1
#   }
# }

## Rasterize tile4
#install.packages("raster")
library(raster)

# rastertile4 = raster(tile4)
# Write raster to file
# filename.tile4 = readline("Please enter the desired file name of the output for Step 2:") # Remove this choice, and just write it out to 1.whichturtle/1.shoreline
# filename.tile4 = file.path(outputdirectory,animalnickname,shorelinefolder,'tile4.asc')
#tile4.asc
# rf = writeRaster(rastertile4, filename=filename.tile4, overwrite=TRUE)

# OLD Eroding raster to skeleton ----------------------------------------------

#install.packages("spgrass6")
### library(spgrass6)

### # Location of GRASS installation:
### ## JC says our ideal is " running R stand-alone and creating a throw-away GRASS environment from within R"
### loc = initGRASS("C:\\Program Files (x86)\\GRASS GIS 6.4.4", home=tempdir(), override = TRUE)
### loc
###
###
### # Import the file to GRASS:
### parseGRASS("r.in.gdal") # commmand description
### execGRASS("r.in.gdal", flags="o", parameters=list(input=filename.tile4, output="tile4"))
### execGRASS("g.region", parameters=list(rast="tile4"))
### gmeta6()
###
###
###
### # Thin the raster map so it can be converted to vectors:
### execGRASS("r.thin", parameters=list(input="tile4", output="tile4thin"))
### tile4thin = readRAST6("tile4thin")
### square.rast.skeleton = as.matrix(tile4thin)
### square.rast.skeleton[is.na(square.rast.skeleton)] = 0
###
### # Reorient
### square.rast.skeleton = t(square.rast.skeleton)
###
### # Return to original dimensions
### rast.skeleton = square.rast.skeleton
### if (diff>0){
###   rast.skeleton = rast.skeleton[,-(n+1:m)]
### } else if (diff<0){
###   rast.skeleton = rast.skeleton[-(m+1:n),]
### }

# NEW Eroding raster to skeleton ----------------------------------------------

skeletonise3_2d.fn <- function(x)
{
  k1 <- matrix(c(0,NA,1,0,1,1,0,NA,1), 3, 3)
  k2 <- matrix(c(NA,1,NA,0,1,1,0,0,NA), 3, 3)
  rot.fn <- function(x) {t(apply(x, 2, rev))}
  hom <- function(x,k) morph(x, k, operator = "==", merge = "all", value = 1)

  repeat
  {
    previous <- x
    for (i in 1:4)
    {
      x <- x & !hom(x,k1)
      storage.mode(x) <- "integer"
      x <- x & !hom(x,k2)
      storage.mode(x) <- "integer"
      k1 <- rot.fn(k1)
      k2 <- rot.fn(k2)
    }

    if (isTRUE(all.equal(x, previous)))
      return(x)
  }
}

tileRas <- raster(tile)
plot(tileRas, main = 'original')

rast.skeleton <- skeletonise3_2d.fn(tile)
skelRas <- raster(rast.skeleton)
plot(skelRas, main = 'skelaton_3_2d')


# Counting joints and counting skeleton pixels
skeleton = rast.skeleton
for (i in 1:m){
  for (j in 1:n){
    if (skeleton[i,j]>0){
      joint = 0
      if ((i>1)&(j>1)){
        if (skeleton[i-1,j-1]>0){
          joint = joint+1
        }
      }
      if (i>1){
        if (skeleton[i-1,j]>0){
          joint = joint+1
        }
      }
      if ((i>1)&(j<n)){
        if (skeleton[i-1,j+1]>0){
          joint = joint+1
        }
      }
      if (j<n){
        if (skeleton[i,j+1]>0){
          joint = joint+1
        }
      }
      if ((i<m)&(j<n)){
        if (skeleton[i+1,j+1]>0){
          joint = joint+1
        }
      }
      if (i<m){
        if (skeleton[i+1,j]>0){
          joint = joint+1
        }
      }
      if ((i<m)&(j>1)){
        if (skeleton[i+1,j-1]>0){
          joint = joint+1
        }

      }
      if (j>1){
        if (skeleton[i,j-1]>0){
          joint = joint+1
        }
      }
      if (joint>0){
        skeleton[i,j] = joint
      }
    }
  }
}

# Data points to skeleton -------------------------------------------------

## Sample data
# Enter file name of data points
# filepath.data = readline("Please enter the file name of the data points:")
filepath.data = file.path('2.Turtles','cq100output','cq100moutput.txt')
#cq100moutput.txt
#cq13100moutput.txt
#cw13100moutput.txt
#hk100moutput.txt
#hk13100moutput.txt
#hm100moutput.txt
#hn100moutput.txt
#hp100moutput.txt
#hq100moutput.txt
#hq13100moutput.txt
#hv100moutput.txt
#hw100moutput.txt
# Enter number of rows in file
# numrowstoread = readline("Please enter the number of rows in the data points file:")
numrowstoread = dim(skeleton)[1]
#239
# Upload file
sdata = as.matrix(read.table(filepath.data, header=F, skip=6, sep=" ", nrows=as.integer(numrowstoread)))
sdata[sdata==-9999] = 0
# to reorient properly
# sdata2 = t(apply(sdata,2,rev))
# sdata3 = sdata2[-c(2241),]

# Map to be exported
map = tile
# map[sdata3==1] = 10
map[sdata==1] = 10
map[map==1] = 11
map[map==0] = 12

# pdf of original with data points
pdf1.name = readline("Please enter the desired file name of the pdf file
                     with the original map and data points:")
#original2.pdf
pdf(pdf1.name)
image(map,main="Example 2: Original location with data points of hw",
      xlab="Easting",ylab="Northing")
dev.off()

# Moving data points to skeleton
data = matrix(0, m, n)
for (i in 1:m){
  for (j in 1:n){
    # if (sdata3[i,j]==1){
    if (sdata[i,j]==1){
      minD = Inf
      for (k in 1:m){
        for (l in 1:n){
          if (rast.skeleton[k,l]==1){
            D = (i-k)^2+(j-l)^2
            if (D<minD){
              minD = D
              miny = k
              minx = l
            }
          }
        }
      }
      data[miny,minx] = 1
    }
  }
}

# Map of data points on skeleton
map2 = rast.skeleton
map2[data==1] = 10
map2[map2==1] = 11
map2[map2==0] = 12

# pdf of data points on skeleton
pdf2.name = readline("Please enter the desired file name of the pdf file
                     with the skeleton and associated data points:")
#datapoints_skeleton2.pdf
pdf(pdf2.name)
image(map2,main="Example 2: Data points for hk13 associated to skeleton",xlab="Easting",ylab="Northing")
dev.off()

# Counting points ---------------------------------------------------------

# Joint matrix and numbering data points
joints = matrix(0,m,n)
jointnum = 1
datanum = 1
for (i in 1:m){
  for (j in 1:n){
    if (skeleton[i,j]>2){
      joints[i,j] = jointnum
      jointnum = jointnum+1
    }
    if (data[i,j]>0){
      data[i,j] = datanum
      datanum = datanum+1
    }
  }
}
jointnum = jointnum-1
datanum = datanum-1
distsize = jointnum+datanum

# Adding numbered data points and joints to skeleton
sjd = data
k = 1
for (i in 1:m){
  for (j in 1:n){
    if (data[i,j]==0 & joints[i,j]!=0){
      sjd[i,j] = datanum+k
      k = k+1
    }
  }
}
distsize = max(sjd)
for (i in 1:m){
  for (j in 1:n){
    if (data[i,j]==0 & joints[i,j]==0 & skeleton[i,j]!=0){
      sjd[i,j] = distsize+1
    }
  }
}

# Distance between all points ---------------------------------------------

travel <- function(skel2,i,j,d){
  x = 0
  while (x==0){
    if ((i>1)&(j>1)){
      if (skel2[i-1,j-1]>0){
        x = 1
        i = i-1
        j = j-1
        d = d+sqrt(2)
        break
      }
    }
    if (i>1){
      if (skel2[i-1,j]>0){
        x = 1
        i = i-1
        d = d+1
        break
      }
    }
    if ((i>1)&(j<n)){
      if (skel2[i-1,j+1]>0){
        x = 1
        i = i-1
        j = j+1
        d = d+sqrt(2)
        break
      }
    }
    if (j<n){
      if (skel2[i,j+1]>0){
        x = 1
        i = i
        j = j+1
        d = d+1
        break
      }
    }
    if ((i<m)&(j<n)){
      if (skel2[i+1,j+1]>0){
        x = 1
        i = i+1
        j = j+1
        d = d+sqrt(2)
        break
      }
    }
    if (i<m){
      if (skel2[i+1,j]>0){
        x = 1
        i = i+1
        j = j
        d = d+1
        break
      }
    }
    if ((i<m)&(j>1)){
      if (skel2[i+1,j-1]>0){
        x = 1
        i = i+1
        j = j-1
        d = d+sqrt(2)
        break
      }
    }
    if (j>1){
      if (skel2[i,j-1]>0){
        x = 1
        i = i
        j = j-1
        d = d+1
        break
      }
    }
    x = 2
  }
  result = c(i,j,d,x)
  return(result)
}

# Calculate distance between pairs of data points or joints
dist = {}
skelseg = sjd
row = 1
segment = distsize+2
for (i in 1:distsize){
  j = which((sjd==i), arr.ind=TRUE)
  x = j[1]
  y = j[2]
  jsize = skeleton[x,y]
  skel2 = skeleton
  for (t in 1:jsize){
    check = 0
    d = 0
    x = j[1]
    y = j[2]
    skel2[x,y] = 0
    while (check==0){
      new = travel(skel2,x,y,d)
      x = new[1]
      y = new[2]
      d = new[3]
      z = new[4]
      skel2[x,y] = 0
      if (sjd[x,y]!=(distsize+1)){
        if (row==1){
          dist = c(i,sjd[x,y],d,segment)
          row = row+1
          segment = segment+1
          check = 1
        } else if (row==2){
          if (dist[1]==i&dist[2]==sjd[x,y]){
            if (dist[3]>d){
              dist[3] = d
              row = row+1
              segment = segment+1
              check = 1
            } else{
              segment = segment+1
              check = 1
            }
          } else{
            new.dist.row = c(i,sjd[x,y],d,segment)
            dist = rbind(dist, new.dist.row)
            row = row+1
            segment = segment+1
            check = 1
          }
        } else if (any(dist[,1]==i&dist[,2]==sjd[x,y])){
          oldrow = which((dist[,1]==i&dist[,2]==sjd[x,y]),arr.ind=TRUE)
          old = dist[oldrow,3]
          if (old>d){
            new.dist.row = c(i,sjd[x,y],d,segment)
            dist = rbind(dist, new.dist.row)
            row = row+1
            segment = segment+1
            check = 1
          } else{
            segment = segment+1
            check = 1
          }
        } else{
          new.dist.row = c(i,sjd[x,y],d,segment)
          dist = rbind(dist, new.dist.row)
          row = row+1
          segment = segment+1
          check = 1
        }
      } else{
        skelseg[x,y] = segment
      }
      if (z==2){
        check = 1
        segment = segment+1
      }
    }
  }
}

# Shortest paths between data points  -------------------------------------------------------------------

#install.packages("igraph")
library(igraph)

dist2 = as.data.frame(dist[,1:3])

# Make column names
names(dist2) = c("from","to","distance")

# Transform into graph
dgraph = graph.data.frame(dist2,directed=FALSE)

# Finding shortest paths between data points
distances = matrix(0,datanum,datanum)
paths = {}
for (i in 1:(datanum-1)){
  for (j in (i+1):datanum){
    path = get.shortest.paths(dgraph,i,j,weight=E(dgraph)$distance)
    new.paths = list(c(i,j), c(unlist(path[1][[1]][1])))
    paths = rbind(paths, new.paths)
    distances[i,j] = shortest.paths(dgraph,i,j,weight=E(dgraph)$distance)
    distances[j,i] = distances[i,j]
  }
}
distances2 = distances

# Minimum spanning tree ---------------------------------------------------

# Change distances matrix for Prim's algorithm
distances3 = distances2
for (i in 1:(datanum-1)){
  for (j in (i+1):datanum){
    for (k in 1:dim(paths)[1]){
      if (paths[,1][[k]][1]==i&paths[,1][[k]][2]==j){
        for (l in 2:(length(paths[,2][[k]])-1)){
          if (paths[,2][[k]][l]<=datanum){
            distances3[i,j] = Inf
            distances3[j,i] = Inf
          }
        }
        break
      }
    }
  }
}
dpairs = matrix(0,datanum*(datanum-1)/2,3)
row=1
for (i in 1:(datanum-1)){
  for (j in (i+1):datanum){
    dpairs[row,1] = i
    dpairs[row,2] = j
    dpairs[row,3] = distances2[i,j]
    row = row+1
  }
}

# Finding minimum distances using Prim's algorithm
dpairsleft = dpairs
minimum = matrix(0,datanum-1,3)
used = 1
l = 1
check = datanum-1
while (check!=0){
  minrow = matrix(NA,dim(dpairsleft)[1],3)
  col = 1
  for (k in 1:l){
    for (i in 1:dim(dpairsleft)[1]){
      for (j in 1:2){
        if (dpairsleft[i,j]==used[k]){
          minrow[col,] = dpairsleft[i,]
          col = col+1
        }
      }
    }
  }
  minimum[l,3] = min(minrow[,3],na.rm=TRUE)
  minimum[l,1] = minrow[which(minrow[,3]==minimum[l,3])[1],1]
  minimum[l,2] = minrow[which(minrow[,3]==minimum[l,3])[1],2]
  l = l+1
  used[l] = minimum[l-1,1]
  if (any(duplicated(used))){
    used[l] = minimum[l-1,2]
  }
  for (i in seq(dim(dpairsleft)[1],1,by=-1)){
    if (any(dpairsleft[i,1]==used)&any(dpairsleft[i,2]==used)){
      dpairsleft = dpairsleft[-i,]
    }
  }
  check = check-1
}

minpathways = {}
for (i in 1:(datanum-1)){
  for (j in 1:dim(paths)[1]){
    if (paths[,1][[j]][1]==minimum[i,1]&paths[,1][[j]][2]==minimum[i,2]){
      row = j
    }
  }
  for (k in 1:(length(paths[,2][[row]])-1)){
    new.minpath = c(paths[,2][[row]][k],paths[,2][[row]][k+1])
    minpathways = rbind(minpathways,new.minpath)
  }
}

# Getting rid of duplicated rows
minpathways2 = minpathways
for (i in 1:dim(minpathways)[1]){
  if (minpathways[i,1]>minpathways[i,2]){
    minpathways2[i,1] = minpathways[i,2]
    minpathways2[i,2] = minpathways[i,1]
  }
}
minpathways3 = minpathways2[!duplicated(minpathways2),]

minpathways3 = cbind(minpathways3,matrix(0,dim(minpathways3)[1],2))
for (i in 1:dim(minpathways3)[1]){
  row = which(minpathways3[i,1]==dist[,2]&minpathways3[i,2]==dist[,1])
  minpathways3[i,3] = dist[row,3]
  minpathways3[i,4] = dist[row,4]
}

mstdistance=sum(minpathways3[,3])

mst=matrix(0,m,n)
for (i in 1:m){
  for (j in 1:n){
    if (sjd[i,j]<=datanum&sjd[i,j]!=0){
      mst[i,j] = 1
    } else {
      for (k in 1:dim(minpathways3)[1]){
        if (skelseg[i,j]==minpathways3[k,4]){
          mst[i,j] = 1
        } else if (sjd[i,j]==minpathways3[k,1]){
          mst[i,j] = 1
        } else if (sjd[i,j]==minpathways3[k,2]){
          mst[i,j] = 1
        }
      }
    }
  }
}

# Map of data points on skeleton
map3 = rast.skeleton
map3[mst==1] = 7
map3[data!=0] = 5
map3[map3==1] = 10
map3[map3==0] = 12

# pdf of minimum spanning tree
pdf3.name = readline("Please enter the desired file name of the pdf file
                     with the minimum spanning tree:")
#minspantree2.pdf
pdf(pdf3.name)
image(map3,main="Example 2: Minimum Spanning Tree for hw",xlab="Easting",ylab="Northing")
dev.off()

# Output distance of minimum spanning tree
write.table(mstdistance, "mstdistancehw.txt", sep="\t")
