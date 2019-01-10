/*
 * Simple Stencil example
 *
 *
 * Brian J Gravelle
 * gravelle@cs.uoregon.edu
 *
 */

// point for each mesh location
// contains three "physical values"
//   - desciptions explain behavior on each time step
//   - 3 doubles are 192 bits which don't fit nicely into 
//     cache lines or vector registers
//   - note that 4 is 256
struct Mesh
{

  double heat;    // update to the average to itself and the neighbors
  double volume;  // divide the vaule amoung itself and neighbors
  double fancy;   // just for fun add the heat and volume
  
};

// create and fill the mesh with starting values
int init_mesh(struct Mesh ***mesh, int x_size, int y_size);

// for point (x,y) in mesh get the list of neighbors 
void get_neighbors(int x_size, int y_size, int x, int y, int neighbors[9][2]);

// perform one iteration of the timestep
void do_timestep(struct Mesh **mesh, struct Mesh **new_mesh, int x_size, int y_size, double time, double dt);

// print a grid of the mesh
void print_mesh(struct Mesh **mesh, int x_size, int y_size);

// liberate the memory
void free_mesh(struct Mesh **mesh, int x_size, int y_size);





