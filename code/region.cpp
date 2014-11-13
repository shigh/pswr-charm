
#include "region.hpp"

Region::Region(int K_, int nt_, int ny_,
			   double dy_, int nx_, double dx_):
	K(K_), nt(nt_), ny(ny_), dy(dy_), nx(nx_), dx(dx_)
{

	dt_vals = std::vector<double>(K,     0);
	x0      = std::vector<double>(nx*ny, 0);
	west    = std::vector<double>(ny*nt, 0);
	east    = std::vector<double>(ny*nt, 0);
	north   = std::vector<double>(nx*nt, 0);
	south   = std::vector<double>(nx*nt, 0);
	t = 0;
	build_solver();

	// Setup chunk logic
	int cs = (int)nt/K; // chunk size
	chunk_start = std::vector<int>(K, 0);
	chunk_size  = std::vector<int>(K, cs);
	for(int i=0; i<K; i++)
		chunk_start[i] = cs*i;
	
}

void Region::build_solver()
{

}

void Region::apply_solver()
{

}

void Region::time_step()
{

}

void Region::time_step(int N)
{

}

void Region::time_step_chunk(int N)
{

}

void Region::set_dt(double dt, int N)
{

}

void Region::set_dt(double dt)
{

}

double Region::get_dt(int N)
{
	return 0;
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

int Region::get_chunk_elem_start(boundary_t bndy, int N)
{

	int start = -1;
	if(bndy==WEST || bndy==EAST)
		start = ny*chunk_start[N];
	else if(bndy==NORTH || bndy==SOUTH)
		start = nx*chunk_start[N];
	assert(start!=-1);

	return start;

}

void Region::update_boundary(boundary_t bndy, const double* vals, int N)
{

	int n_set = get_chunk_n_elems(bndy, N);
	int start = get_chunk_elem_start(bndy, N);

	std::vector<double>& vec = get_boundary_vector(bndy);
	for(int i=chunk_start[N];
		i<chunk_size[N]; i++)
		for(int j=0; j<n_set; j++)
			vec[start+j] = vals[j];

}

std::vector<double> Region::get_boundary(boundary_t bndy, int N)
{

	int n_set = get_chunk_n_elems(bndy, N);
	int start = get_chunk_elem_start(bndy, N);

	std::vector<double>& vec = get_boundary_vector(bndy);
	std::vector<double> vals = std::vector<double>(n_set, 0);
	for(int i=chunk_start[N];
		i<chunk_size[N]; i++)
		for(int j=0; j<n_set; j++)
			vals[j] = vec[start+j];

	return vals;

}
