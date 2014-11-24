/* Local Variables: */
/* c-set-offset: 'substatement-open 0 */
/* c++-tab-always-indent: t */
/* c-basic-offset: 4 */
/* c-indent-level: 4 */
/* tab-stop-list: '(4 8 12 16 20 24 28 32 36 40 44 48 52 56 60) */
/* tab-width: 4 */
/* indent-tabs-mode: t */
/* End: */

#pragma once

#include <vector>
#include <memory>
#include "utils.hpp"
#include "solver.hpp"

class Region
{

private:

	std::shared_ptr<Solver> solver;

	// dt at each time step
	std::vector<double> dt_vals;
	// N where the chunk starts
	std::vector<int> chunk_start;
	// Number of steps in a chunk
	std::vector<int> chunk_size;
	// Initial values
	std::vector<double> x0;
	// Temp vector
	std::vector<double> x;
	// Boundaries
	std::vector<double> west, east, north, south;

	// Number of chunks
	int K;
	int curr_chunk, curr_chunk_ind, curr_ind;

	// Model params
	int nt, ny, nx;
	double dy, dx;

	// Overlap logic
	int overlap;

	void update_boundary_arrays(const std::vector<double>& vec, int chunk, int chunk_ind);

	std::vector<double>& get_boundary_vector(boundary_t bndy);

	// ---These functions are usefull for indexing into the bndy arrays
	// Total number of elements in a chunk
	int get_chunk_n_elems(boundary_t bndy, int N);
	// The element index where the chunk starts
	int get_chunk_start_index(boundary_t bndy, int N);
	// Get the index of current chunk start point
	int get_curr_start_index(boundary_t bndy);
	// First index of a boundary in chunk
	int get_start_index(boundary_t bndy, int chunk, int chunk_ind);
	// ---

	// Advance one time step
	void time_step();

public:

	Region(int K_, int overlap, int nt_, int ny_, double dy_, int nx_, double dx_,
		   std::vector<double> x0, std::shared_ptr<Solver> solver);


	/*! Advance each step in chunk N
	 */
	void time_step_chunk();

	/*! Reset to start new iteration
	 */
	void reset();

	/*! Set dt for chunk NT
	 */
	void set_dt(double dt, int N);

	/*! dt at chunk N
	 */
	double get_dt(int N);

	/*! Number of steps in chunk N
	 */
	double get_chunk_size(int N);

	/*! Set boundary chunk N
	 */
	void set_boundary(boundary_t bndy, const double* vals, int N);

	/*! Get boundary chunck N
	 */
	std::vector<double> get_boundary(boundary_t bndy, int N);

	// msg* get_boundary(...);

};
