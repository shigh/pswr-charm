#include <vector>

#pragma once

// I have never had good luck with enums
typedef int boundary_t;
#define WEST  (boundary_t)0
#define EAST  (boundary_t)1
#define NORTH (boundary_t)2
#define SOUTH (boundary_t)3

class Region
{

private:

	// dt at each time step
	std::vector<double> dt_vals;
	// N where the chunk starts
	std::vector<int> chunk_start;
	// Number of steps in a chunk
	std::vector<int> chunk_size;
	// Initial values
	std::vector<double> x0;
	// Boundaries
	std::vector<double> west, east, north, west;

	int nt, nT, ny, nx, t;
	double dy, dx;

	// Initialize the PETc solver
	void build_solver();

	// Apply solver to current t
	void apply_solver();

public:

	Region(int nT_, int nt_, int ny_, double dy_, int nx_, double dx_):
		nT(nT_), nt(nt_), ny(ny_), dy(dy_), nx(nx_), dx(dx_)
	{
		dt_vals = std::vector<double>(nt,    0);
		x0      = std::vector<double>(nx*ny, 0);
		west    = std::vector<double>(ny*nt, 0);
		east    = std::vector<double>(ny*nt, 0);
		north   = std::vector<double>(nx*nt, 0);
		south   = std::vector<double>(nx*nt, 0);
		t = 0;
		build_solver();
	}

	/*! Advance one time step
	 */
	void time_step();

	/*! Advance N time steps
	 */
	void time_step(int N);
	
	/*! Advance each step in chunck N
	 */
	void time_step_chunk(int N);

	/*! Set dt for chunk NT
	 */
	void set_dt(int N, double dt);

	/*! Set dt for at all t
	 */
	void set_dt(double dt);

	/*! dt at step N
	 */
	double get_dt(int N);

	/*! Update boundary chunk N
	 */
	void update_boundary(boundary_t bndy, const double* vals, int N);

	/*! Update the entire boundary
	 */
	void update_boundary(boundary_t bndy, const double* vals);

	// std::vector<double> get_boundary(...);
	// msg* get_boundary(...);

};
