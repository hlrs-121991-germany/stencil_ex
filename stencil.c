/*
 * Simple Stencil example
 * Main program example with openmp thread parallelism
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

#ifndef X_SIZE
#define X_SIZE 10000
#endif

#ifndef Y_SIZE
#define Y_SIZE 20000
#endif

#ifndef TIME
#define TIME 0.0
#endif

#ifndef STEP
#define STEP 1.0
#endif

#ifndef TIME_STOP
#define TIME_STOP 10.0
#endif



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
      _mesh[i][j].fancy  = H + V;
      V += 1000;
    }
    H += 100;
  }

  *mesh = _mesh;

  return err;
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

double pythag(double x1, double y1, double x2, double y2) {
  return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
}

// perform one iteration of the timestep
void do_timestep(struct Mesh **mesh, struct Mesh **new_mesh, int x_size, int y_size, double time, double dt) {
#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  int neighbors[9][2];
  double dt2 = dt*dt;
  double C = 0.25;

  int _x, _y, n;

  #pragma simd
  for (_x = 0; _x < x_size; _x++) {
    for (_y = 0; _y < y_size; _y++) {
      new_mesh[_x][_y].heat   = 0;
      new_mesh[_x][_y].volume = 0;
      new_mesh[_x][_y].fancy  = -2*dt2 * mesh[_x][_y].fancy * C;
    }
  }

  #pragma omp parallel for private(_y, n)
  for (_x = 0; _x < x_size; _x++) {
    for (_y = 0; _y < y_size; _y++) {

      int neighbors[9][2];

      get_neighbors(x_size, y_size, _x, _y, neighbors);

      for(n = 0; n < 9; n++) {
        new_mesh[_x][_y].heat += mesh[neighbors[n][X]][neighbors[n][Y]].heat;
      }
      new_mesh[_x][_y].heat /= 9;

      #pragma simd
      for(n = 0; n < 9; n++) {
        new_mesh[neighbors[n][X]][neighbors[n][Y]].volume += mesh[_x][_y].volume/9;
      }

      for(n = 0; n < 9; n++){
        double dist2 = pythag(_x, _y, neighbors[n][X], neighbors[n][Y]); // dx^2
        new_mesh[_x][_y].fancy += -2*dt2 * mesh[neighbors[n][X]][neighbors[n][Y]].fancy / (dist2*C);
      }

    }
  }

} // do time step

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
      printf("%10.2f", mesh[i][j].fancy);
    }
    printf("\n\n");
  }

}

int test_small_mesh() {
#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
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
  do_timestep(mesh_1, mesh_2, x_size, y_size, time, 1.0);
  printf("print_mesh...\n");
  print_mesh(mesh_2, x_size, y_size);
  printf("free_mesh...\n");
  free_mesh(mesh_1, x_size, y_size);
  free_mesh(mesh_2, x_size, y_size);

  return err;
}


int run_custom_mesh(int x_size, int y_size, double time, double step, double time_stop) {
#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  int err = FALSE;

  struct Mesh **mesh_1 = NULL;
  struct Mesh **mesh_2 = NULL;

  double wall_tot_start, wall_tot_end;
  double wall_init_start, wall_init_end;
  double wall_step_start, wall_step_end;
  double wall_free_start, wall_free_end;


  printf("\n\nRunning new Stencil with \n\
    x_size     = %d \n\
    y_size     = %d \n\
    start time = %f \n\
    time step  = %f \n\
    end time   = %f \n\n",
    x_size, y_size, time, step, time_stop);

  wall_tot_start = omp_get_wtime();
  wall_init_start = omp_get_wtime();
  printf("init_mesh......."); fflush(stdout);
  err = err | init_mesh(&mesh_1, x_size, y_size);
  err = err | init_mesh(&mesh_2, x_size, y_size);
  if(mesh_1 == NULL) return 1;
  if(mesh_2 == NULL) return 1;
  wall_init_end = omp_get_wtime();
  printf("%fs\n", (wall_init_end - wall_init_start));

  while(time < time_stop) {

    printf("timestep %.2f...", time); fflush(stdout);
    wall_step_start = omp_get_wtime();
    do_timestep(mesh_1, mesh_2, x_size, y_size, time, step);
    time += step;
    wall_step_end = omp_get_wtime();
    printf("%fs\n", (wall_step_end - wall_step_start));

    printf("timestep %.2f...", time); fflush(stdout);
    wall_step_start = omp_get_wtime();
    do_timestep(mesh_2, mesh_1, x_size, y_size, time, step);
    time += step;
    wall_step_end = omp_get_wtime();
    printf("%fs\n", (wall_step_end - wall_step_start));

  }


  printf("free_mesh.......\n"); fflush(stdout);
  wall_free_start = omp_get_wtime();
  free_mesh(mesh_1, x_size, y_size);
  free_mesh(mesh_2, x_size, y_size);
  wall_free_end = omp_get_wtime();
  printf("%fs\n", (wall_free_end - wall_free_start));

  wall_tot_end = omp_get_wtime();
  printf("\n total time: %fs\n", (wall_tot_end - wall_tot_start));

  return err;
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

  // err = err | test_small_mesh();
  // printf("\n\n");
  err = err | run_custom_mesh(X_SIZE, Y_SIZE, TIME, STEP, TIME_STOP);

  return err;
}




