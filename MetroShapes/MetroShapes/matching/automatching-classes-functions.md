# Functions Automatching 

`MetroSahpes/core/MetroShapes.cpp` is the primary class; calling the match, smooth, mixed layout functions. 



`MetroShapes/matching/Stationmatching` contains my old Manuel matching approach `



`MetroShapes/matching/Automatching` contains the new approaches to automatically find a shape inside the metro system. Containing following functions: 

+ `findPath()` Starts the search for paths. In this function, for each vertex a new spaceDijgast run is started, and the best found path stored. 
+ `spaceDijgstar()` is the main algorithm. This part is more or less equal to the compass path approach. A pseudocode of the Djigstar like approach can be found in her phd thesis https://elib.uni-stuttgart.de/handle/11682/3036
+ `computeDFD()` initialised the computation of the Frechet distance.
  + `computeSegemntDFD()` computes the distance between two curves. The basis of the recursive algorithm is taken from here http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.90.937&rep=rep1&type=pdf
  + `costSegmentDFD()` returns the actual cost of two sections of the two polygons. Because we compare directions rather than spatial distances, the messurement is taken form here: https://www.win.tue.nl/~mdberg/Papers/2011/bc-gfdfd.pdf#; I also found there old java implementation, but I have not used it: https://u.pcloud.link/publink/show?code=kZGSTeXZ0ULcWVPG7gp0YWLEU8i0AjaauyDX

+ All functions including the `...IPR()` are functions which I have created during the experiments with the inflection point representation of the paths. 