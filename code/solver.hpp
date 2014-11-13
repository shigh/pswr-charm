/* Local Variables: */
/* c-set-offset: 'substatement-open 0 */
/* c++-tab-always-indent: t */
/* c-basic-offset: 4 */
/* c-indent-level: 4 */
/* tab-stop-list: '(4 8 12 16 20 24 28 32 36 40 44 48 52 56 60)) */
/* tab-width: 4 */
/* indent-tabs-mode: t */
/* End: */



#include <vector>

#pragma once

/*
  Usage:
  Solver s = Solver(ny, dy, nx, dx);
  s.set_rhs(b, west, east, north, south, dt);
  s.solve(x);
*/

// Depending on how PETc does its solves we may need to change some of
// the data layout.  I am not sure yet what the best signature for the
// solve functions is, if you build it using these we should be able
// to modifify it relatively easily to do what we want.
class Solver
{

private:

	double dt;
	int ny, nx;
	double dy, dx;
	std::vector<double> rhs;

public:

	Solver(int ny_, double dy_, int nx_, double dx_):
		ny(ny_), dy(dy_), nx(nx_), dx(dx_)
	{
		rhs = std::vector<double>(nx*ny, 0);
	}

	/*! Solve for x in Ax=b
	 */
	void solve(std::vector<double>& x);

	// It would probablly be better to do this with iterators. We can
	// look into that later.
	void set_rhs(const std::vector<double>& b,
				 double* west, double* east,
				 double* north, double* south,
				 double dt);

};
