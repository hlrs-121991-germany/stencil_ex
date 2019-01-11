/*
 * Simple Stencil example
 * patterns for the stencil define STENCIL_TYPE to choose
 *
 * Brian J Gravelle
 * gravelle@cs.uoregon.edu
 *
 */


#include "mesh.h"

#ifdef USE_CALI
#include <caliper/cali.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#ifndef STENCIL_TYPE
#define STENCIL_TYPE 0
#endif

#if STENCIL_TYPE == 0

// for point (x,y) in mesh get the list of neighbors 
// returns the number of neighbors
int get_neighbors(int x_size, int y_size, int x, int y, int neighbors[9][2]) {
#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  int _x, _y, n;

  n = 0;
  for (_x = -1; _x <= 1; _x++) {
    for (_y = -1; _y <= 1; _y++) {

      neighbors[n][X] = x + _x;
      neighbors[n][Y] = y + _y;

      if((neighbors[n][X] < 0) || (neighbors[n][X] >= x_size) )
        neighbors[n][X] = x;      
      if((neighbors[n][Y] < 0) || (neighbors[n][Y] >= y_size) )
        neighbors[n][Y] = y;

      n++;
    }
  }

  return n;
}


#elif STENCIL_TYPE == 1
// for point (x,y) in mesh get the list of neighbors 
// returns the number of neighbors
int get_neighbors(int x_size, int y_size, int x, int y, int neighbors[9][2]) {
#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  int _x, _y, n=0;

  for (_x = -1; _x <= 1; _x++) {

    neighbors[n][X] = x + _x;
    if((neighbors[n][X] < 0) || (neighbors[n][X] >= x_size) )
      neighbors[n][X] = x;     

    n++;
  }

  for (_y = -1; _y <= 1; _y++) {

    neighbors[n][Y] = y + _y;
    if((neighbors[n][Y] < 0) || (neighbors[n][Y] >= y_size) )
      neighbors[n][Y] = y;

    n++;
  }


  return n;
}

#endif