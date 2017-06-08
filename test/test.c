#include "test.h"

int main () {
  configInfo par;

  par.radius = 1.2;
  par.numPoints = 10000;
  par.numGridDensMaxima = 0;
  par.gridDensMaxValues = NULL;
  par.gridDensMaxLoc    = NULL;

  makeAndPlotPoints(&par);

  free(par.gridDensMaxLoc);
  free(par.gridDensMaxValues);

  return 0;
}

