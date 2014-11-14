
#include "solver.hpp"

Solver::Solver(int ny_, double dy_, int nx_, double dx_):
	ny(ny_), dy(dy_), nx(nx_), dx(dx_)
{
	rhs = std::vector<double>(nx*ny, 0);
}

/*! Solve for x in Ax=b
 */
void Solver::solve(std::vector<double>& x)
{

}

// It would probablly be better to do this with iterators. We can
// look into that later.
void Solver::set_rhs(const std::vector<double>& b,
		double* west, double* east,
		double* north, double* south)
{

}

void Solver::set_dt(double dt_)
{

}

double Solver::get_dt()
{

}

