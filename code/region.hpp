/* Local Variables: */
/* c-set-offset: 'substatement-open 0 */
/* c++-tab-always-indent: t */
/* c-basic-offset: 4 */
/* c-indent-level: 4 */
/* tab-stop-list: '(4 8 12 16 20 24 28 32 36 40 44 48 52 56 60) */
/* tab-width: 4 */
/* indent-tabs-mode: t */
/* End: */

#include <vector>

#pragma once

#include "utils.hpp"
#include "solver.hpp"

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
	std::vector<double> west, east, north, south;

	// Number of chunks
	int K;
	int nt, ny, nx, t;
	double dy, dx;

	// Initialize the PETc solver
	void build_solver();

	// Apply solver to current t
	void apply_solver();

	std::vector<double>& get_boundary_vector(boundary_t bndy);

public:

	Region(int K_, int nt_, int ny_, double dy_, int nx_, double dx_);

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

	/*! dt at chunk N
	 */
	double get_dt(int N);

	/*! Update boundary chunk N
	 */
	void update_boundary(boundary_t bndy, const double* vals, int N);

	/*! Get boundary chunck N
	 */
	std::vector<double> get_boundary(boundary_t bndy, int N);

	// msg* get_boundary(...);

};
