#include "test.h"
#include "../src/tree_random.h"

double exampleNumDensyFunc0(double *r){
  /* Implements an offset Lorenzian profile with a centred circular cutoff. */

  double lorCentre[2] = {0.32,0.76}, lorHalfWidth=0.01, lorMaxValue=2.8;
  double lorHalfWidthSquared = lorHalfWidth*lorHalfWidth;
  double lorNumerator = lorMaxValue*lorHalfWidthSquared;
  double truncCentre[2] = {0.15,0.05}, truncRadius=0.85;
  double truncRadiusSquared = truncRadius*truncRadius;
  double ringCentre[2] = {0.47,-0.06}, ringInnerR=0.32, ringOuterR=0.35, ringExtraGridNumDensy=0.01;
  double rectX[2]={-0.3,0.2},rectY[2]={-0.3,0.1},rectAngle=0.2,rectExtraGridNumDensy=0.005;
  double diff,diffSquared,cosA,sinA,x1,y1;
  int i;
  double unnormalizedGridNumDens=0.0;

  const int addRing=1,addRect=1;

  if(N_DIMS!=2){
    printf("Error! par->numDims should ==2 but it ==%d\n", N_DIMS);
    exit(1);
  }

  /* Exclude points outside of the truncation radius. */
  diffSquared = 0.0;
  for (i=0;i<N_DIMS;i++){
    diff = r[i] - truncCentre[i];
    diffSquared += diff*diff;
  }
  if (diffSquared>truncRadiusSquared)
    return 0.0;
 
  /* Add Lorentzian. */
  diffSquared = 0.0;
  for (i=0;i<N_DIMS;i++){
    diff = r[i] - lorCentre[i];
    diffSquared += diff*diff;
  }
  unnormalizedGridNumDens = lorNumerator/(diffSquared + lorHalfWidthSquared);

  if(addRing){
    diffSquared = 0.0;
    for (i=0;i<N_DIMS;i++){
      diff = r[i] - ringCentre[i];
      diffSquared += diff*diff;
    }
    if(diffSquared>=ringInnerR*ringInnerR && diffSquared<ringOuterR*ringOuterR)
      unnormalizedGridNumDens += ringExtraGridNumDensy;
  }

  if(addRect){
    cosA = cos(rectAngle);
    sinA = sin(rectAngle);
    x1 = r[0]*cosA - r[1]*sinA;
    y1 = r[1]*cosA + r[0]*sinA;
    if(x1>=rectX[0] && x1<rectX[1] && y1>=rectY[0] && y1<rectY[1])
      unnormalizedGridNumDens += rectExtraGridNumDensy;
  }

  return unnormalizedGridNumDens;
}

double exampleNumDensyFunc1(double *r){
  /* This density function consists of several overlaid (and summing) rectangles, each of constant number density. The whole is truncated within a circular area. */

  double truncCentre[2] = {0.15,0.05}, truncRadius=0.95;
  double truncRadiusSquared = truncRadius*truncRadius;
  double rect0X[2]={-0.5, 0.3 },rect0Y[2]={-0.5, 0.1},rect0ExtraGridNumDensy=0.005;
  double rect1X[2]={ 0.0, 0.07},rect1Y[2]={-0.7, 0.7},rect1ExtraGridNumDensy=0.01;
  double rect2X[2]={-0.15,0.15},rect2Y[2]={-0.05,0.3},rect2ExtraGridNumDensy=0.05;
  double diff,diffSquared;
  int i;
  double unnormalizedGridNumDens=0.0005;

  const int addRing=1,addRect=1;

  if(N_DIMS!=2){
    printf("Error! par->numDims should ==2 but it ==%d\n", N_DIMS);
    exit(1);
  }

  /* Exclude points outside of the truncation radius. */
  diffSquared = 0.0;
  for (i=0;i<N_DIMS;i++){
    diff = r[i] - truncCentre[i];
    diffSquared += diff*diff;
  }
  if (diffSquared>truncRadiusSquared)
    return 0.0;
 
  if(r[0]>=rect0X[0] && r[0]<rect0X[1] && r[1]>=rect0Y[0] && r[1]<rect0Y[1])
    unnormalizedGridNumDens += rect0ExtraGridNumDensy;
  if(r[0]>=rect1X[0] && r[0]<rect1X[1] && r[1]>=rect1Y[0] && r[1]<rect1Y[1])
    unnormalizedGridNumDens += rect1ExtraGridNumDensy;
  if(r[0]>=rect2X[0] && r[0]<rect2X[1] && r[1]>=rect2Y[0] && r[1]<rect2Y[1])
    unnormalizedGridNumDens += rect2ExtraGridNumDensy;

  return unnormalizedGridNumDens;
}

/*....................................................................*/
void _prepareForMonitor(double *worldLo, double *worldWidth){
  int devId;

  /* Do setups for pgplot calls in testMonitoFunc.
  */
  devId = cpgopen("/XWIN");
  cpgsvp(0.01,0.99,0.03,0.99);
  cpgwnad(worldLo[0],worldLo[0]+worldWidth[0],worldLo[1],worldLo[1]+worldWidth[1]);
  cpgsfs(2);
}

/*....................................................................*/
void _testMonitorFunc2D(const int cellI, double fieldOrigin[N_DIMS]\
    , double fieldWidth[N_DIMS], double (*outRandLocations)[N_DIMS]\
    , unsigned int firstPointI, unsigned int numPoints){

  float *xs, *ys, xLo, xHi, yLo, yHi, xMid, yMid;
  unsigned int i_u,j_u;
  char numberStr[6];

  xLo = fieldOrigin[0];
  yLo = fieldOrigin[1];
  xHi = xLo + fieldWidth[0];
  yHi = yLo + fieldWidth[1];
  xMid = 0.5*(xLo + xHi);
  yMid = 0.5*(yLo + yHi);

  xs = malloc(sizeof(double)*numPoints);
  ys = malloc(sizeof(double)*numPoints);

  for(i_u=0;i_u<numPoints;i_u++){
    j_u = i_u + firstPointI;
    xs[i_u] = (float)outRandLocations[j_u][0];
    ys[i_u] = (float)outRandLocations[j_u][1];
    if (xs[i_u]<xLo || xs[i_u]>xHi)
      printf("Bad X value %f at %d\n", xs[i_u], (int)i_u);
    if (ys[i_u]<yLo || ys[i_u]>yHi)
      printf("Bad Y value %f at %d\n", ys[i_u], (int)i_u);
  }

  if (numPoints>0){
    cpgsci(2);
    cpgpt(numPoints, xs, ys, -1);
  }

  if(cellI>=0){
    sprintf(numberStr, "%d", cellI);
    cpgsci(3);
    cpgtext(xMid, yMid, numberStr);
  }

  cpgsci(1);
  cpgrect(xLo, xHi, yLo, yHi);

  free(ys);
  free(xs);

}

/*....................................................................*/
void _exampleMonitorFunc(const int numDims, const int cellI, double fieldOrigin[N_DIMS]\
    , double fieldWidth[N_DIMS], unsigned int desiredNumPoints, double (*outRandLocations)[N_DIMS]\
    , unsigned int firstPointI, unsigned int actualNumPoints){

  if(numDims!=2)
exit(1);

  _testMonitorFunc2D(cellI, fieldOrigin, fieldWidth, outRandLocations, firstPointI, actualNumPoints);
}

/*....................................................................*/
void _cleanUpAfterMonitor(){
  cpgend();
}

/*....................................................................*/
void _plotPoints(treeRandConstType *rinc, double (*outRandLocations)[N_DIMS], double *outRandDensities){
  _prepareForMonitor(rinc->wholeFieldOrigin, rinc->wholeFieldWidth);
  _testMonitorFunc2D(-1, rinc->wholeFieldOrigin, rinc->wholeFieldWidth, outRandLocations, 0, rinc->desiredNumPoints);
  _cleanUpAfterMonitor();
}

/*....................................................................*/
void makeAndPlotPoints(configInfo *par){
  double *outRandDensities=NULL;
  double (*outRandLocations)[N_DIMS]=NULL;
  treeRandConstType rinc;
  int levelI=0,di,i;
  unsigned int startI=0;

  setConstDefaults(&rinc);

  rinc.numDims = 2;
  rinc.desiredNumPoints = (unsigned int)par->numPoints;
  rinc.wholeFieldOrigin[0] = -par->radius;
  rinc.wholeFieldOrigin[1] = -par->radius;
  rinc.wholeFieldWidth[0] = 2.0*par->radius;
  rinc.wholeFieldWidth[1] = 2.0*par->radius;
  rinc.verbosity = 1;
  rinc.monitorFunc = NULL;/*_exampleMonitorFunc;*/

  rinc.totalNumHighPoints = par->numGridDensMaxima;

  if(rinc.totalNumHighPoints>0){
    rinc.allHighPointLoc   = malloc(sizeof(*(rinc.allHighPointLoc  ))*rinc.totalNumHighPoints);
    rinc.allHighPointDensy = malloc(sizeof(*(rinc.allHighPointDensy))*rinc.totalNumHighPoints);
    for(i=0;i<rinc.totalNumHighPoints;i++){
      for(di=0;di<rinc.numDims;di++){
        rinc.allHighPointLoc[i][di] = par->gridDensMaxLoc[i][di];
      }
      rinc.allHighPointDensy[i] = par->gridDensMaxValues[i];
    }
  }else{
    rinc.allHighPointLoc = NULL;
    rinc.allHighPointDensy = NULL;
  }

  outRandDensities = malloc(sizeof(*outRandDensities)*par->numPoints); /* Not used at present; and in fact they are not useful outside this routine, because they are not the values of the physical density at that point, just what exampleNumDensyFunc() returns, which is not necessarily the same thing. */
  outRandLocations = malloc(sizeof(*outRandLocations)*par->numPoints);

  if(rinc.monitorFunc!=NULL)
    _prepareForMonitor(rinc.wholeFieldOrigin, rinc.wholeFieldWidth);

  treeGenerateRandoms(&rinc, exampleNumDensyFunc0, outRandLocations, outRandDensities);

  if(rinc.monitorFunc==NULL)
    _plotPoints(&rinc, outRandLocations, outRandDensities);
  else
    _cleanUpAfterMonitor();

  free(outRandLocations);
  free(outRandDensities);
}







