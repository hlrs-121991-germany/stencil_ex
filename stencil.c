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

#ifndef FILE_NAME
#define FILE_NAME "openmp_results.txt"
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
      _mesh[i][j].avg = H;
      _mesh[i][j].sum = V;
      _mesh[i][j].pde = i*j;
      _mesh[i][j].dep = H+V;
      V += 1000;
    }
    H += 100;
  }

  *mesh = _mesh;

  return err;
}

double pythag(double x1, double y1, double x2, double y2) {
  return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
}

// perform one iteration of the timestep
void do_timestep(struct Mesh **mesh, struct Mesh **_temp_mesh, int x_size, int y_size, double time, double dt) {
#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  // OpenMP specific variables
  char row_status[x_size];
  int thr_id;
  struct Mesh **temp_mesh;

  double dt2 = dt*dt;
  double C = 0.25, dx2 = 1.0;

  int _x, _y, n, i, j;

  for (i = 0; i < x_size; i++) {
    row_status[i] = 0;
  }
  // goal set the last couple rows to start at a higher count since they have fewer dependancies
  for (i = x_size - TEMP_ROWS; i < x_size; i++) {
    row_status[i] = TEMP_ROWS + (i - x_size);
  }

  // looping over all rows of the matrix
  // main source of paralleism 
  #pragma omp parallel for private(_y, n, thr_id, temp_mesh, dx2)
  for (_x = 0; _x < x_size; _x++) {

    int neighbors[NUM_NEIGHBORS][2];

    // establish temporary mesh for this thread
    thr_id = omp_get_thread_num();
    temp_mesh = &_temp_mesh[thr_id*TEMP_ROWS];

    // set list of temporary indices to untaken (-1)
    int temp_inds[TEMP_ROWS];
    for (int t = 0; t < TEMP_ROWS; t++) temp_inds[t] = -1;

    // get next available temp index
    int temp_x = -1;
    int t;
    while (temp_x < 0) {

      // check temp rows to see if the original os done being used
      for (t = 0; t < TEMP_ROWS; t++) {
        if(temp_inds[t] >= 0){ 
          if(row_status[temp_inds[t]] == (TEMP_ROWS-1)){
            for (_y = 0; _y < y_size; _y++) {
              mesh[temp_inds[t]][_y].avg = temp_mesh[t][_y].avg;
              mesh[temp_inds[t]][_y].sum = temp_mesh[t][_y].sum;
              mesh[temp_inds[t]][_y].pde = temp_mesh[t][_y].pde;
              mesh[temp_inds[t]][_y].dep = temp_mesh[t][_y].dep;
            }
            temp_inds[t] = -1;
          }
        }
      }

      // choose next temporary row to use
      for (t = 0; t < TEMP_ROWS; t++) {
        if(temp_inds[t] < 0){
          temp_x = t;
          temp_inds[t] = _x;
          break;
        }

      }

    } // while

    // fill next temp row with old values
    for (_y = 0; _y < y_size; _y++) {
      temp_mesh[temp_x][_y].avg = 0;
      temp_mesh[temp_x][_y].sum = 0;
      temp_mesh[temp_x][_y].pde = -2*dt2 * mesh[_x][_y].pde * C;
      temp_mesh[temp_x][_y].dep = -2*dt2 * mesh[_x][_y].dep * C;
    }

    // actually do some computation
    for (_y = 0; _y < y_size; _y++) {

      get_neighbors(x_size, y_size, _x, _y, neighbors);

      for(n = 0; n < NUM_NEIGHBORS; n++) {
        temp_mesh[temp_x][_y].avg += mesh[neighbors[n][X]][neighbors[n][Y]].avg;
      }
      temp_mesh[temp_x][_y].avg /= NUM_NEIGHBORS;

      for(n = 0; n < NUM_NEIGHBORS; n++) {
        temp_mesh[temp_x][_y].sum += mesh[neighbors[n][X]][neighbors[n][Y]].sum/NUM_NEIGHBORS;
      }

      for(n = 0; n < NUM_NEIGHBORS; n++){
        dx2 = pythag(_x, _y, neighbors[n][X], neighbors[n][Y]); // dx^2
        temp_mesh[temp_x][_y].pde += (-2*dt2 * mesh[neighbors[n][X]][neighbors[n][Y]].pde) / ((dx2 + 1.0) * C);
      }

      for(n = 0; n < NUM_NEIGHBORS; n++){
        dx2 = pythag(_x, _y, neighbors[n][X], neighbors[n][Y]); // dx^2
        temp_mesh[temp_x][_y].dep += (mesh[neighbors[n][X]][neighbors[n][Y]].avg*dt2 * \
                                      mesh[neighbors[n][X]][neighbors[n][Y]].dep) / \
                                      ((dx2 + mesh[neighbors[n][X]][neighbors[n][Y]].sum) * C);
      }
    } // _y loop

    // set row satus
    #pragma omp critical
    for (t = 1; t < TEMP_ROWS-1; t++) {
      if((_x-t) >= 0) 
        row_status[_x-t] += 1;
    }

  } // _x loop

} // do time step

// print the mesh
void print_mesh(struct Mesh **mesh, int x_size, int y_size) {
#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  int i, j;

  for (i = 0; i < x_size; i++) {
    printf("x = %d\n", i);
    for (j = 0; j < y_size; j++) {
      printf("%10.2e ", mesh[i][j].avg);
    }
    printf("\n");

    for (j = 0; j < y_size; j++) {
      printf("%10.2e ", mesh[i][j].sum);
    }
    printf("\n");

    for (j = 0; j < y_size; j++) {
      printf("%10.2e ", mesh[i][j].pde);
    }
    printf("\n");

    for (j = 0; j < y_size; j++) {
      printf("%10.2e ", mesh[i][j].dep);
    }
    printf("\n\n");
  }

}

// print the mesh to file
void output_mesh(FILE* file, struct Mesh **mesh, int x_size, int y_size) {
#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  int i, j;

  for (i = 0; i < x_size; i++) {
    fprintf(file, "x = %d\n", i);
    for (j = 0; j < y_size; j++) {
      fprintf(file, "%10.2e ", mesh[i][j].avg);
    }
    fprintf(file, "\n");

    for (j = 0; j < y_size; j++) {
      fprintf(file, "%10.2e ", mesh[i][j].sum);
    }
    fprintf(file, "\n");

    for (j = 0; j < y_size; j++) {
      fprintf(file, "%10.2e ", mesh[i][j].pde);
    }
    fprintf(file, "\n");

    for (j = 0; j < y_size; j++) {
      fprintf(file, "%10.2e ", mesh[i][j].dep);
    }
    fprintf(file, "\n\n");
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
  err = err | init_mesh(&mesh_2, TEMP_ROWS, y_size);
  if(mesh_1 == NULL) return 1;
  if(mesh_2 == NULL) return 1;
  printf("print_mesh...\n");
  print_mesh(mesh_1, x_size, y_size);
  printf("do_timestep...\n");
  do_timestep(mesh_1, mesh_2, x_size, y_size, time, 1.0);
  printf("print_mesh...\n");
  print_mesh(mesh_1, x_size, y_size);
  printf("free_mesh...\n");
  free_mesh(mesh_1, x_size, y_size);
  free_mesh(mesh_2, TEMP_ROWS, y_size);

  return err;
}


int run_custom_mesh(int x_size, int y_size, double time, double step, double time_stop) {
#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  int err = FALSE;

  struct Mesh **main_mesh = NULL;
  struct Mesh **temp_mesh = NULL;

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
  err = err | init_mesh(&main_mesh, x_size, y_size);
  err = err | init_mesh(&temp_mesh, TEMP_ROWS, y_size);
  if(main_mesh == NULL) return 1;
  if(temp_mesh == NULL) return 1;
  wall_init_end = omp_get_wtime();
  printf("%fs\n", (wall_init_end - wall_init_start));

#ifdef DO_IO
  printf("output to file....."); fflush(stdout);
  double io_start = omp_get_wtime();
  FILE* file = fopen(FILE_NAME, "w+");
  fprintf(file, "\n\nRunning new Stencil with \n\
    x_size     = %d \n\
    y_size     = %d \n\
    start time = %f \n\
    time step  = %f \n\
    end time   = %f \n\n",
    x_size, y_size, time, step, time_stop);
  output_mesh(file, main_mesh, x_size, y_size);
  printf("%fs\n", (omp_get_wtime() - io_start));
#endif

  while(time < time_stop) {

    printf("timestep %.2f...", time); fflush(stdout);
    wall_step_start = omp_get_wtime();
    do_timestep(main_mesh, temp_mesh, x_size, y_size, time, step);
    time += step;
    wall_step_end = omp_get_wtime();
    printf("%fs\n", (wall_step_end - wall_step_start));

  }

#ifdef DO_IO
  io_start = omp_get_wtime();
  printf("output to file....."); fflush(stdout);
  fprintf(file, "\n\n");
  output_mesh(file, main_mesh, x_size, y_size);
  fprintf(file, "\n\n");
  fclose(file);
  printf("%fs\n", (omp_get_wtime() - io_start));
#endif

  printf("free_mesh.......\n"); fflush(stdout);
  wall_free_start = omp_get_wtime();
  free_mesh(main_mesh, x_size, y_size);
  free_mesh(temp_mesh, TEMP_ROWS, y_size);
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




