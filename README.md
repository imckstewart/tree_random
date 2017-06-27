Project tree_random
===================

The Treegrid algorithm attempts to generate random points within an orthotopic working space which follow a distribution given by some user-supplied number-density function. (An orthotope is a generalization to N dimensions of the concept of a rectangle.) The algorithm is designed to avoid the lengthy run times which can plague the (otherwise robust) rejection method when there are significant variations in the value of the point number density function within the working space. The algorithm makes use of a 2^D tree approach (where D is the number of spatial dimensions) to split the working space into smaller and smaller segments until a space is reached within which the variation of the number density function is within acceptable bounds.

The algorithm is performed in two sections: the first (implemented in function _constructTheTree) decides on the structure of the tree; the second (implemented in function _fillTheTree) decides on the number of points required for each sub-cell, then generates these points via the standard rejection method.

The rejection method takes a starting set of points which have a constant (i.e. flat) probability distribution, then tests each one of these, rejecting points based on the ratio between the value of the desired number density at that point and its maximum value in the relevant subcell. The user has the choice to use either completely random points for the starting set, or points which have a flat probability distribution but which are chosen from a quasi-random sequence, the default being the Halton sequence:

  Halton, J. (1964), Algorithm 247: Radical-inverse quasi-random point sequence, ACM, p. 701

One reason for choosing a quasi-random sequence is to avoid pairs of points which are close together. This is useful if the points are used for the vertices of a finite-element mesh, because it is inefficient to solve these sorts of problems if the FE cells are very divergent in size and shape.

Points chosen from a quasi-random sequence require no seed and thus are always the same set. This has obvious implications for aliasing-type problems if the locations of the starting points at least have the same fractional positions in each of the sub-cells. Two modifications have been made to tree_random to ameliorate such problems, as given below:

  - The division between subcells is not at exactly 1/2 the distance along each dimension axis of the parent cell: some dither is introduced.

  - The axes of each new sub-cell are randomly permuted, and the signs of displacement along each axis are chosen randomly.

One final adjustment which was needed for good results with quasi-random sequences was to introduce a small buffer zone near sub-cell boundaries within which points are not allowed. Without this, points in abutting sub-cells could approach each other arbitrarily closely, even though points within the same sub-cell cannot.

Any C program which wishes to use this code should look broadly like the following:

```
#include "tree_random.h"

double exampleNumDensyFunc(double *r){
  < Appropriate code to return the desired relative number density of the random points >
}

int main() {
  treeRandConstType rinc;
  double (*outRandLocations)[N_DIMS]=NULL;
  double *outRandDensities=NULL;

  setConstDefaults(&rinc);

  < Set non-default rinc values >

  treeGenerateRandoms(&rinc, exampleNumDensyFunc, outRandLocations, outRandDensities);

  < Make use of the output values >

  free(outRandLocations);
  free(outRandDensities);

  return 0;
}
```

