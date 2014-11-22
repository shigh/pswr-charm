
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <math.h>
#include <algorithm>
#include <cstddef>
#include <vector>
#include <iostream>
#include <fstream>
#include "../utils.hpp"
#include "../solver.hpp"

BOOST_AUTO_TEST_SUITE( solver_tests )

BOOST_AUTO_TEST_CASE( heatbtcs_zero_boundary_convg )
{
	// Using a general solution of the form:
	// u(x, y, t) = sin(kx)sin(ky)exp(-2k^2 t)
	// on a square domain
	// This is not a well resolved problem so we will use
	// a relatively higher tolerance
	double tol = 0.001;

	int N = 200;
	int ind;
	double dx, dt, L, k;
	double xj, yi;

	L  = M_PI;
	k  = 1.;
	dt = .01;
	dx = L/((double)(N-1));
	
	HeatSolverBTCS solver = HeatSolverBTCS(N, dx, N, dx);
	solver.set_dt(dt);

	auto b        = std::vector<double>(N*N, 0);
	auto x        = std::vector<double>(N*N, 0);
	auto expected = std::vector<double>(N*N, 0);
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
		{
			ind = j + i*N;
			xj = j*dx;
			yi = i*dx;
			b[ind]        = sin(k*xj)*sin(k*yi);
			expected[ind] = sin(k*xj)*sin(k*yi)*exp(-(2*k*k)*dt);
		}

	auto west  = std::vector<double>(N, 0);
	auto east  = std::vector<double>(N, 0);
	auto north = std::vector<double>(N, 0);
	auto south = std::vector<double>(N, 0);

	solver.set_rhs(b, &west[0], &east[0], &north[0], &south[0]);
	solver.solve(x);

	double error = 0.;
	for(int i=0; i<N*N; i++)
		error = std::max(error, std::abs(expected[i] - x[i]));

	BOOST_CHECK_LT( error, tol );

}


BOOST_AUTO_TEST_SUITE_END()
