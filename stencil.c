/*
 * Simple Stencil example
 * Main program example
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

#define TRUE  1
#define FALSE 0

#define X 0
#define Y 1

// create and fill the mesh with starting values
int init_mesh(struct Mesh ***mesh, int x_size, int y_size) {
#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  struct Mesh **_mesh;

  int err = FALSE;
  int i, j;

  double H = 100;
  double V = 1000;

  // allocate memory and verify that it worked
  _mesh = (struct Mesh**) malloc(x_size * sizeof(struct Mesh*));
  if(_mesh == NULL) err = TRUE;
  for (i = 0; i < x_size; ++i) {
    _mesh[i] = (struct Mesh*) malloc(y_size * sizeof(struct Mesh));
    if(_mesh[i] == NULL) err = TRUE;
  }

  // define starting values
  for (i = 0; i < x_size; i++) {
    V = 1000;
    for (j = 0; j < y_size; j++) {
      _mesh[i][j].heat   = H;
      _mesh[i][j].volume = V;
      _mesh[i][j].jff    = H + V;
      V += 1000;
    }
    H += 100;
  }

  *mesh = _mesh;

  return err;
}

// for point (x,y) in mesh get the list of neighbors 
void get_neighbors(int x_size, int y_size, int x, int y, int neighbors[9][2]) {
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

}

// perform one iteration of the timestep
void do_timestep(struct Mesh **mesh, struct Mesh **new_mesh, int x_size, int y_size, double time) {
#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  int neighbors[9][2];

  int _x, _y, n;

  for (_x = 0; _x < x_size; _x++) {
    for (_y = 0; _y < y_size; _y++) {
      new_mesh[_x][_y].heat   = 0;
      new_mesh[_x][_y].volume = 0;
      new_mesh[_x][_y].jff    = 0;
    }
  }

  for (_x = 0; _x < x_size; _x++) {
    for (_y = 0; _y < y_size; _y++) {

      get_neighbors(x_size, y_size, _x, _y, neighbors);

      for(n = 0; n < 9; n++) {
        new_mesh[_x][_y].heat += mesh[neighbors[n][X]][neighbors[n][Y]].heat;
      }
      new_mesh[_x][_y].heat /= 9;

      for(n = 0; n < 9; n++) {
        new_mesh[neighbors[n][X]][neighbors[n][Y]].volume += mesh[_x][_y].volume/9;
      }

    }
  }

  for (_x = 0; _x < x_size; _x++)
    for (_y = 0; _y < y_size; _y++)  
      new_mesh[_x][_y].jff = new_mesh[_x][_y].heat + new_mesh[_x][_y].volume;

}

// print the mesh
void print_mesh(struct Mesh **mesh, int x_size, int y_size) {
#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  int i, j;

  for (i = 0; i < x_size; i++) {

    for (j = 0; j < y_size; j++) {
      printf("%10.2f", mesh[i][j].heat);
    }
    printf("\n");

    for (j = 0; j < y_size; j++) {
      printf("%10.2f", mesh[i][j].volume);
    }
    printf("\n");

    for (j = 0; j < y_size; j++) {
      printf("%10.2f", mesh[i][j].jff);
    }
    printf("\n\n");
  }

}

// liberate the memory
void free_mesh(struct Mesh **mesh, int x_size, int y_size) {
#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  int i;

  for (i = 0; i < x_size; ++i) {
    free(mesh[i]);
  }
  free(mesh);

}


// main function
int main(int argc, char **argv) {

#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
cali_id_t thread_attr = cali_create_attribute("thread_id", CALI_TYPE_INT, CALI_ATTR_ASVALUE | CALI_ATTR_SKIP_EVENTS);
#pragma omp parallel
{
cali_set_int(thread_attr, omp_get_thread_num());
}
#endif

  int err = FALSE;

  struct Mesh **mesh_1 = NULL;
  struct Mesh **mesh_2 = NULL;
  int x_size = 5;
  int y_size = 10;
  double time = 0.0;

  printf("init_mesh...\n");
  err = err | init_mesh(&mesh_1, x_size, y_size);
  err = err | init_mesh(&mesh_2, x_size, y_size);
  if(mesh_1 == NULL) return 1;
  if(mesh_2 == NULL) return 1;
  printf("print_mesh...\n");
  print_mesh(mesh_1, x_size, y_size);
  printf("do_timestep...\n");
  do_timestep(mesh_1, mesh_2, x_size, y_size, time);
  printf("print_mesh...\n");
  print_mesh(mesh_2, x_size, y_size);
  printf("free_mesh...\n");
  free_mesh(mesh_1, x_size, y_size);
  free_mesh(mesh_2, x_size, y_size);

  return err;
}




