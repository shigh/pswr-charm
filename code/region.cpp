
#include "region.hpp"
#include "solver.hpp"

Region::Region(int K_, int overlap_, int nt_, int ny_,
			   double dy_, int nx_, double dx_,
			   std::vector<double> x0_, std::shared_ptr<Solver> solver_):
	K(K_), overlap(overlap_), nt(nt_), ny(ny_), dy(dy_), nx(nx_), dx(dx_),
	x0(x0_), solver(solver_)
{

	dt_vals = std::vector<double>(K,     0);
	x = x0; // Copy x0 into x
	west    = std::vector<double>(ny*nt, 0);
	east    = std::vector<double>(ny*nt, 0);
	north   = std::vector<double>(nx*nt, 0);
	south   = std::vector<double>(nx*nt, 0);

	// Setup chunk logic
	// TODO Handle case where nt not a multiple of K
	int cs = (int)nt/K; // chunk size
	chunk_start = std::vector<int>(K, 0);
	chunk_size  = std::vector<int>(K, cs);
	for(int i=0; i<K; i++)
		chunk_start[i] = cs*i;
	chunk_size[K-1] = nt-cs*(K-1);

	// Set the time t==0 boundary array values to
	// the intial values
	curr_chunk     = 0;
	curr_chunk_ind = 0;
	curr_ind       = 0;
	update_boundary_arrays();

	reset();

}

void Region::reset()
{

	x = x0; // Copy x0 into x
	// Set region to initial state
	curr_chunk     = 0;
	curr_chunk_ind = 1;
	curr_ind       = 1;
	solver->set_dt(dt_vals[curr_chunk]);

}

void Region::time_step()
{

	// Update solver boundarys
	double* west  = &west[get_curr_start_index(WEST)];
	double* east  = &east[get_curr_start_index(EAST)];
	double* north = &north[get_curr_start_index(NORTH)];
	double* south = &south[get_curr_start_index(SOUTH)];

	solver->set_rhs(x, west, east, north, south);

	solver->solve(x);

	update_boundary_arrays();

	++curr_chunk_ind;
	++curr_ind;
	if(curr_chunk_ind == chunk_size[curr_chunk])
	{
		++curr_chunk;
		curr_chunk_ind = 0;
		curr_ind = chunk_start[curr_chunk];
		solver->set_dt(dt_vals[curr_chunk]);
	}
	
}

void Region::update_boundary_arrays()
{

	int start;
	start = get_curr_start_index(EAST);
	for(int i=0; i<ny; i++)
		east[start+i] = x[i*nx+overlap];

	start = get_curr_start_index(WEST);
	for(int i=0; i<ny; i++)
		west[start+i] = x[(i+1)*nx-1-overlap];

	start = get_curr_start_index(NORTH);
	for(int i=0; i<nx; i++)
		north[start+i] = x[(ny-1)*nx+i-overlap*nx];

	start = get_curr_start_index(SOUTH);
	for(int i=0; i<nx; i++)
		south[start+i] = x[i+overlap*nx];

}

void Region::time_step_chunk()
{

	int n_steps = chunk_size[curr_chunk];
	// Handle the known initial values
	if(curr_chunk==0) --n_steps;

	for(int i=0; i<n_steps; i++)
		time_step();

}

void Region::set_dt(double dt, int N)
{
	dt_vals[N] = dt;
	if(N==curr_chunk)
		solver->set_dt(dt_vals[curr_chunk]);
}

double Region::get_dt(int N)
{
	return dt_vals[N];
}

double Region::get_chunk_size(int N)
{
	return chunk_size[N];
}

std::vector<double>& Region::get_boundary_vector(boundary_t bndy)
{

	if(bndy==NORTH)
		return north;
	else if(bndy==SOUTH)
		return south;
	else if(bndy==WEST)
		return west;
	else if(bndy==EAST)
		return east;
	assert(0);

}

int Region::get_chunk_n_elems(boundary_t bndy, int N)
{

	int n_set = -1;
	if(bndy==WEST || bndy==EAST)
		n_set = ny*chunk_size[N];
	else if(bndy==NORTH || bndy==SOUTH)
		n_set = nx*chunk_size[N];
	assert(n_set!=-1);

	return n_set;

}

int Region::get_chunk_start_index(boundary_t bndy, int N)
{

	int start = -1;
	if(bndy==WEST || bndy==EAST)
		start = ny*chunk_start[N];
	else if(bndy==NORTH || bndy==SOUTH)
		start = nx*chunk_start[N];
	assert(start!=-1);

	return start;

}

int Region::get_curr_start_index(boundary_t bndy)
{

	int start = -1;
	if(bndy==WEST || bndy==EAST)
		start = ny*chunk_start[curr_chunk]+ny*curr_chunk_ind;
	else if(bndy==NORTH || bndy==SOUTH)
		start = nx*chunk_start[curr_chunk]+nx*curr_chunk_ind;
	assert(start!=-1);

	return start;

}


void Region::update_boundary(boundary_t bndy, const double* vals, int N)
{

	int n_set = get_chunk_n_elems(bndy, N);
	int start = get_chunk_start_index(bndy, N);

	std::vector<double>& vec = get_boundary_vector(bndy);
	for(int i=chunk_start[N];
		i<chunk_size[N]; i++)
		for(int j=0; j<n_set; j++)
			vec[start+j] = vals[j];

}

std::vector<double> Region::get_boundary(boundary_t bndy, int N)
{

	int n_set = get_chunk_n_elems(bndy, N);
	int start = get_chunk_start_index(bndy, N);

	std::vector<double>& vec = get_boundary_vector(bndy);
	std::vector<double> vals = std::vector<double>(n_set, 0);

	for(int i=0; i<n_set; i++)
		vals[i] = vec[i+start];

	return vals;

}
