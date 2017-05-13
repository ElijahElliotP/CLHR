# Complex Linear Home Range Estimator
The `clhr` package uses kernel-based image processing and weighted, undirected graphs to establish a consistent and repeatable method of estimating complex linear home ranges for species living in rivers and around island systems. This process was first automated with R by Nicole McLaren in 2015 under the supervision of Jeff Cardille based on the paper by Ouellette and Cardille from 2011, _The Complex Linear Home Range Estimator: Representing the Home Range of River Turtles Moving in Multiple Channels_. More generally this package allows the user to compare the lengths of sets of points in network systems.

New in this release is the switch away from `GRASS` (Geographic Resources Analysis Support System) to [`mmand`][mmand] (Mathematical Morphology in Any Number of Dimensions) by Jon Clayden, to perform the required skeletonisation.

### Basic Operation
The basic operation is as follows: set-up your files (see [Input Specifications](#input-specifications)), configure your preferences text (an example is included in the text itself), and then run (source) the `CLHR-Rscript.R` script.

## Contents
- [Example Data](#example-Data)
- [ASCII File to Binary Matrix](#ascii-file-to-binary-matrix)
- [Extract Major Skeleton](#extract-major-skeleton)
- [Shift Data Points to a Region Skeleton](#shift-data-points-to-a-region-skeleton)
- [Connectivity Matrix](#connectivity-matrix)
- [Input Specifications](#input-specifications)
- [Output Files](#output-files)
- [Known Issues and Limitations](#known-issues-and-limitations)

## Example Data
A sample data set of observations of Map Turtles and Red-Eared Sliders in the Niagara River system is provided by Bryan Haas. The river system delineated here has an area of approximately 53,000 km<sup>2</sup>. This set includes 12 different sets of observations at 3 different resolutions, 10 m, 20 m and 100 m and a region file for each resolution. They have already been converted from `.shp` files into in ASCII file format with both `.txt` and `.asc` extensions for demonstration.
![Niagara Region Rivers](http://i.imgur.com/BeyDkW8.jpg)
Above: The Region of Interest courtesy of Google Earth and DigitalGlobe.

## ASCII File to Binary Matrix
`prepASCtoMatrix.fn` serves to ‘clean up’ and convert `nodata` and `NA` values the ascii file might have, as well as convert any non-zero presence values to 1, giving a binary matrix as the name implies.

```R
regionTile <- prepASCtoMatrix.fn('2.Turtle/hk20output/hk20output.asc')
rastertile <- copySourceGeoSpatToRas.fn(regionTile, '2.Turtle/hk20output/hk20output.asc')
plot(rastertile, col = c('white', 'black'))
```
![The River System in ASCII](http://i.imgur.com/jPoPGDQ.png)

The river system of the Niagara region at 20 m resolution, from the provided `.asc` file. At this resolution the entire dataset (all 12 observation sets) takes a couple of hours to run.

## Extract Major Skeleton
In addition to extracting the centreline or skeleton for the image, we need that skeleton to be a single unbroken line. Instead of throwing an error in the case of non-contiguous regions, the workaround is to keep the largest portion of the potential home-range and collapse all the individual location points to that line. This tool is partially powered by the [`components()`](https://github.com/jonclayden/mmand#connected-components) function in the [mmand package by Jon Clayden][mmand] which can identify and group parts of an image based on their adjacency.
See [Known Issues and Limitations](#known-issues-and-limitations) for more on this.

### Skeletonisation with mmand
Switching to [`skeletonise()`](https://github.com/jonclayden/mmand#skeletonisation) in [mmand][mmand], away from the GRASS tool `r.thin`, is the principal change to the way this script functions. This allows the user to stay within R and use easily installable libraries. Briefly, this function (as used here) erodes binary matrices using [hit-or-miss transform][homt] and a pair of 3x3 kernels until only a centreline remains.

```R
majorSkeleton <- extractMajorSkeleton.fn(regionTile)
skelRas <- copySourceGeoSpatToRas.fn(majorSkeleton, locationFilename)
plot(skelRas, col = c('white', 'black')
```
![The Major Skeleton](http://i.imgur.com/Kmw1dmC.png)

The Extracted Major Skeleton: Because of the relatively fine lines produced, the skeleton for the 20 m resolution would appear broken at this scale.  The skeleton is overlaid on the original region and the depicted line in the 100 m resolution for clarity. Often branches or “spurs” can occur during skeletonisation. These are exaggerations of features, typically protrusions, from the original shape and can be expected; they are a consequence of the method of skeletonisation. By its nature a Minimum Spanning Tree (MST), the line connecting all of our points of interest, ignores these features.
    
## Shift Data Points to a Region Skeleton
Given a binary presence map and a binary skeleton, it is then necessary to associate each presence point with the nearest location on the skeleton of the river system. This is done by comparing the location of every individual with every point on the skeleton in a ‘moving window’ fashion. The result is a new ‘presence’ layer that resembles the skeleton populated with observations and a new object in the global environment called shiftedPoints that is to be passed to subsequent CLHR operations.
    
![Individual Initial Positions](http://i.imgur.com/mgl42nb.png)
Above, the observations in red in their initial locations along a portion of the skeleton.

```R
mapIndivsOnLandSkel <- shiftDataToSkeleton.fn(sdata, majorSkeleton)
map2Ras <- copySourceGeoSpatToRas.fn(mapIndivsOnLandSkel, locationFilename) 
plot(map2Ras, col = c('red', 'black', ‘white’) 
```
![Individual Final Positions](http://i.imgur.com/aGoywmP.png)
Now in orange, the observation points have been shifted onto the nearest points on the skeleton.

## From Start to Finish
The final operations are still grouped together inside `calculateTheMinimumSpanningTree.fn()` and neither they nor it should be called directly (it's still quite messy). Instead, set-up your files, configure your preferences text, and then run (source) the script so that all the bookkeeping is done for you. Eventually all intermediate steps will be compartmentalised to facility other uses.

![The MST](http://i.imgur.com/mkVGEQh.png)

The observation points, still in orange are all connected by the green line: The Minimum Spanning Tree. The script also creates a matrix of the MST alone and gives the final figures for the total length.

## Connectivity Matrix
`connectivityArray.fn()` outputs a matrix illustrating the connectedness of every non-zero point of the input matrix in terms of the numbers of neighbours it has. This is a ‘moving window’ that travels through every column and every row and looks for 8-connectedness (both a ‘cross’ shape and the diagonals).

## Joint and Point ID Assignment
In order to work with graphs using [`igraph`](https://github.com/igraph/igraph), we need a map of joints or junctions and a map of observation points that will become the vertices for those graphs. `makeJointsNPointsSequential.fn()` assigns IDs sequentially to all non-zero points on the observation layer and all points greater than 2 on the connectivity layer. Points with values of 1 or 2 on that layer represent line ends and midpoints respectively.

## Input Specifications
-	‘Region’ and ‘Individual’ folders should be in the same folder as the script files.
-	Images must be projected in UTM.
-	ASCII files, like those from ArcGIS, should be accompanied by a `.prj` file
    -	They may have either a `.asc` or `.txt` extension, but do not mix them.
-	Region folders and their files must have the same name.
-	‘Individual' folders and files must have the same name.
-	The region folders must share a unique number in common with the individuals' folders.
    - This does not necessarily have to be a resolution, only a unique ID. In this case, you must multiply the MST pixel values by the desired resolution.
-	All ‘individual' folders go together in the user-specified containing folder
    - This allows iteration among different individuals in the same region.

### Example Input Folder Structure
![Example Input Folders](http://i.imgur.com/9w6rIyY.png)

Please excuse the 'output' in the folder names, that is leftover from converting the files into ASCII format.

## Output Files
The script provides automatic generation of `.asc` files and geoTiff files with the appropriate spatial information (the `.asc` files are accompanied by a `.prj`). As well it can be helpful to look at the matrix at many different stages, so a simple `.csv` can be output for the regions and the individuals. One text file per individual is written out with the MST and a `.csv` file is made in the root folder to have all the MSTs tabulated. All of the output files will be found in the output folder specified in the configuration file, inside of the folder with the corresponding resolution number. There is also going to be a `.done` file in each individual’s folder. This is to facilitate faster recovery time if something goes wrong. To redo an individual, you must delete that this.

### Example Output Folder Structure
![Example Output Folders](http://i.imgur.com/SrVlTSh.png)

## Known Issues and Limitations
This segment of the river illustrates the need for the careful creation of the region file.

![Poor Use-Case Before](http://i.imgur.com/nFAk3QQ.png?1)

A landscape with discontinuities like this is not a good use-case for this CLHR transformation / estimation. There are several disconnected pieces of the river system here which would compromise the skeleton and prevent the creation of a Minimum Spanning Tree. To work around this issue all but the largest (by area) contiguous region of the river is removed, and only then do we apply the skeletonisation.

![Poor Use-Case After](http://i.imgur.com/L1NdiiH.png?1)

The result of the workaround is clear. Many points have been collapsed together and the MST is a small portion of what it would otherwise be. Again, this allows the script to function and return an MST but can potentially lead to significant shifts in the number of observation points and an underestimation of total range length. It is advisable to subset the region and run that segment separately, or connect the waterways and re-run. In this instance the script runs but all of the points that were near the eliminated segment shift west.

Less severe, but worth noting, is that two or more points equidistant from a single point on the skeleton would also be collapsed to one individual point. As this script works on binary rasters it can cause the number of total observation points to decrease as any point on the resulting skeleton raster becomes a single presence regardless of how many individuals are recorded there. To resolve any two points the maximum cell size (resolution) must be smaller than their corresponding least-distance. If this is the case then consider a higher image resolution over a subset region to keep computation time down. High resolution sets can be quite slow.
 
[mmand]: https://github.com/jonclayden/mmand
[homt]: http://homepages.inf.ed.ac.uk/rbf/HIPR2/hitmiss.htm
