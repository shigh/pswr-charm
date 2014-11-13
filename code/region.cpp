
#include "region.hpp"

Region::Region(int K_, int nt_, int ny_,
			   double dy_, int nx_, double dx_):
	K(K_), nt(nt_), ny(ny_), dy(dy_), nx(nx_), dx(dx_)
{
	dt_vals = std::vector<double>(K,    0);
	x0      = std::vector<double>(nx*ny, 0);
	west    = std::vector<double>(ny*nt, 0);
	east    = std::vector<double>(ny*nt, 0);
	north   = std::vector<double>(nx*nt, 0);
	south   = std::vector<double>(nx*nt, 0);
	t = 0;
	build_solver();
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

void Region::set_dt(int N, double dt)
{

}

void Region::set_dt(double dt)
{

}

double Region::get_dt(int N)
{
	return 0;
}

void Region::update_boundary(boundary_t bndy, const double* vals, int N)
{

}

void Region::update_boundary(boundary_t bndy, const double* vals)
{

}

std::vector<double> Region::get_boundary(boundary_t bndy, int N)
{
	
}
