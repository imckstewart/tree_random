#include <stdlib.h>
#include <cpgplot.h>
#include "dims.h"

typedef struct {
  double radius,(*gridDensMaxLoc)[N_DIMS],*gridDensMaxValues;
  int numPoints,numGridDensMaxima;
} configInfo;

double exampleNumDensyFunc(double *r);
void makeAndPlotPoints(configInfo *par);
