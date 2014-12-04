
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <math.h>
#include <algorithm>
#include <cstddef>
#include <vector>
#include <iostream>
#include <memory>
#include <fstream>

#include "../utils.hpp"
#include "../solver.hpp"
#include "../region.hpp"


BOOST_AUTO_TEST_SUITE( region_tests )

BOOST_AUTO_TEST_CASE( region_time_step )
{

	int N = 10;
	int K = 5;
	double d = .1;
	// Do not set nt as a multiple of K to make sure we check
	// that region correctly handles this case
	int nt = 21; 

	// Set x0 to ones so the chunk 0 boundaries will all be one
	// (The first N values will be whatever is in x0)
	auto x0 = std::vector<double>(N*N, 1.);
	std::shared_ptr<Solver> solver = std::make_shared<DummySolver>(N, d, N, d);
	Region region = Region(K, 0, nt, N, d, N, d, x0, solver, nt);

	for(int i=0; i<K; i++)
		region.set_dt(i+1, i);
		
	for(int i=0; i<K; i++)
		region.time_step_chunk();

	for(int i=0; i<K; i++)
	{
		auto expected = std::vector<double>(region.get_chunk_size(i)*N, i+1);

		auto bndy = region.get_boundary(EAST, i);
		BOOST_CHECK(std::equal(expected.begin(), expected.end(), bndy.begin()));

		bndy = region.get_boundary(WEST, i);
		BOOST_CHECK(std::equal(expected.begin(), expected.end(), bndy.begin()));

		bndy = region.get_boundary(NORTH, i);
		BOOST_CHECK(std::equal(expected.begin(), expected.end(), bndy.begin()));

		bndy = region.get_boundary(SOUTH, i);
		BOOST_CHECK(std::equal(expected.begin(), expected.end(), bndy.begin()));

	}

	// Reset and run for another iteration
	region.reset();

	for(int i=0; i<K; i++)
		region.set_dt(i+K+1, i);
		
	for(int i=0; i<K; i++)
		region.time_step_chunk();

	for(int i=0; i<K; i++)
	{
		auto expected = std::vector<double>(region.get_chunk_size(i)*N, i+K+1);
		// x0 does not change, so chunk 0 will have its first N values still
		// set to 1
		if(i==0)
			for(int j=0; j<N; j++)
				expected[j] = 1;

		auto bndy = region.get_boundary(EAST, i);
		BOOST_CHECK(std::equal(expected.begin(), expected.end(), bndy.begin()));

		bndy = region.get_boundary(WEST, i);
		BOOST_CHECK(std::equal(expected.begin(), expected.end(), bndy.begin()));

		bndy = region.get_boundary(NORTH, i);
		BOOST_CHECK(std::equal(expected.begin(), expected.end(), bndy.begin()));

		bndy = region.get_boundary(SOUTH, i);
		BOOST_CHECK(std::equal(expected.begin(), expected.end(), bndy.begin()));

	}


}

BOOST_AUTO_TEST_CASE( region_get_set_boundary )
{
	// Make sure get/set boundary functions are working correctly
	int N = 10;
	int K = 5;
	int nt = 21;
	double d = .1;

	auto x0 = std::vector<double>(N*N, 1);
	std::shared_ptr<Solver> solver = std::make_shared<DummySolver>(N, d, N, d);
	Region region = Region(K, 0, nt, N, d, N, d, x0, solver, nt);

	auto set = std::vector<double>(region.get_chunk_size(1)*N, 1);
	region.set_boundary(EAST, &set[0], 1);
	auto get = region.get_boundary(EAST, 1);
	BOOST_CHECK(std::equal(get.begin(), get.end(), set.begin()));
	BOOST_CHECK_EQUAL(get.size(), set.size());

	set = std::vector<double>(region.get_chunk_size(2)*N, 2);
	region.set_boundary(WEST, &set[0], 2);
	get = region.get_boundary(WEST, 2);
	BOOST_CHECK(std::equal(get.begin(), get.end(), set.begin()));
	BOOST_CHECK_EQUAL(get.size(), set.size());

	set = std::vector<double>(region.get_chunk_size(3)*N, 3);
	region.set_boundary(NORTH, &set[0], 3);
	get = region.get_boundary(NORTH, 3);
	BOOST_CHECK(std::equal(get.begin(), get.end(), set.begin()));
	BOOST_CHECK_EQUAL(get.size(), set.size());

	set = std::vector<double>(region.get_chunk_size(4)*N, 4);
	region.set_boundary(SOUTH, &set[0], 4);
	get = region.get_boundary(SOUTH, 4);
	BOOST_CHECK(std::equal(get.begin(), get.end(), set.begin()));
	BOOST_CHECK_EQUAL(get.size(), set.size());

}

BOOST_AUTO_TEST_CASE( region_zero_boundary_convg )
{
	// Using a general solution of the form:
	// u(x, y, t) = sin(kx)sin(ky)exp(-2k^2 t)
	// on a square domain
	// Make sure that the region class correctly handles
	// the solver over many time steps

	double tol = 0.01; // Good enough for government work...

	int N = 200;
	int K = 5;
	int nt = 20;
	int ind;
	double dx, dt, L, k;
	double xj, yi;

	L  = M_PI;
	k  = 1.;
	dt = .01;
	dx = L/((double)(N-1));
	
	auto x0       = std::vector<double>(N*N, 0);
	auto expected = std::vector<double>(N*N, 0);
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
		{
			ind = j + i*N;
			xj = j*dx;
			yi = i*dx;
			x0[ind]       = sin(k*xj)*sin(k*yi);
			expected[ind] = sin(k*xj)*sin(k*yi)*exp(-(2*k*k)*dt*(nt-1));
		}

	std::shared_ptr<Solver> solver = std::make_shared<HeatSolverBTCS>(N, dx, N, dx);
	Region region = Region(K, 0, nt, N, dx, N, dx, x0, solver);

	for(int i=0; i<K; i++)
		region.set_dt(dt, i);
		
	for(int i=0; i<K; i++)
		region.time_step_chunk();

	auto x = region.get_x();

	double error = 0.;
	for(int i=0; i<N*N; i++)
		error = std::max(error, std::abs(expected[i] - x[i]));

	BOOST_CHECK_LT( error, tol );

}

BOOST_AUTO_TEST_CASE( region_four_domain_iteration )
{
	// Using a general solution of the form:
	// u(x, y, t) = sin(kx)sin(ky)exp(-2k^2 t)
	// on a square domain
	//
	// This is a full SWR test
	// 4 subdomains, 2 in each dimension.
	// This checks that everything works together.
	// This should give you a good idea how the charm
	// implementation will work

	double tol = 0.01; // Good enough for government work...

	int N = 200;
	int K = 1;
	int overlap = 20;
	int nt = 10;
	int ind;
	double dx, dt, L, k;
	double xj, yi;

	L  = 2*M_PI;
	k  = 1.;
	dt = .01;
	dx = L/((double)(N-1));

	auto start = std::vector<int>(2, 0);
	auto end   = std::vector<int>(2, 0);
	partition_domain(start, end, N, 2, overlap);

	// Full domain
	auto x0       = std::vector<double>(N*N, 0);
	auto expected = std::vector<double>(N*N, 0);
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
		{
			ind = j + i*N;
			xj = j*dx;
			yi = i*dx;
			x0[ind]       = sin(k*xj)*sin(k*yi);
			expected[ind] = sin(k*xj)*sin(k*yi)*exp(-(2*k*k)*dt*(nt-1));
		}
	std::shared_ptr<Solver> solver = std::make_shared<HeatSolverBTCS>(N, dx, N, dx);
	Region region = Region(K, 0, nt, N, dx, N, dx, x0, solver);

	// Build the four domains
	// West
	const int nxw = end[0]-start[0];
	const int nyw = end[0]-start[0];
	auto x0w       = std::vector<double>(nyw*nxw, 0);
	auto expectedw = std::vector<double>(nyw*nxw, 0);
	for(int i=0; i<nyw; i++)
		for(int j=0; j<nxw; j++)
		{
			ind = j + i*nxw;
			xj = j*dx;
			yi = i*dx;
			x0w[ind]       = sin(k*xj)*sin(k*yi);
			expectedw[ind] = sin(k*xj)*sin(k*yi)*exp(-(2*k*k)*dt*(nt-1));
		}
	std::shared_ptr<Solver> solverw = std::make_shared<HeatSolverBTCS>(nyw, dx, nxw, dx);
	Region regionw = Region(K, overlap, nt, nyw, dx, nxw, dx, x0w, solverw);
	regionw.hold_constant(WEST|SOUTH);

	// East
	const int nxe = end[1]-start[1];
	const int nye = end[0]-start[0];
	auto x0e       = std::vector<double>(nye*nxe, 0);
	auto expectede = std::vector<double>(nye*nxe, 0);
	for(int i=0; i<nye; i++)
		for(int j=0; j<nxe; j++)
		{
			ind = j + i*nxe;
			xj = (j+start[1])*dx;
			yi = i*dx;
			x0e[ind]       = sin(k*xj)*sin(k*yi);
			expectede[ind] = sin(k*xj)*sin(k*yi)*exp(-(2*k*k)*dt*(nt-1));
		}
	std::shared_ptr<Solver> solvere = std::make_shared<HeatSolverBTCS>(nye, dx, nxe, dx);
	Region regione = Region(K, overlap, nt, nye, dx, nxe, dx, x0e, solvere);
	regione.hold_constant(EAST|SOUTH);

	// North west
	const int nxnw = end[0]-start[0];
	const int nynw = end[1]-start[1];
	auto x0nw       = std::vector<double>(nxnw*nynw, 0);
	auto expectednw = std::vector<double>(nxnw*nynw, 0);
	for(int i=0; i<nynw; i++)
		for(int j=0; j<nxnw; j++)
		{
			ind = j + i*nxnw;
			xj = j*dx;
			yi = (i+start[1])*dx;
			x0nw[ind]       = sin(k*xj)*sin(k*yi);
			expectednw[ind] = sin(k*xj)*sin(k*yi)*exp(-(2*k*k)*dt*(nt-1));
		}
	std::shared_ptr<Solver> solvernw = std::make_shared<HeatSolverBTCS>(nynw, dx, nxnw, dx);
	Region regionnw = Region(K, overlap, nt, nynw, dx, nxnw, dx, x0nw, solvernw);
	regionnw.hold_constant(WEST|NORTH);

	// North east
	const int nxne = end[1]-start[1];
	const int nyne = end[1]-start[1];
	auto x0ne       = std::vector<double>(nxne*nyne, 0);
	auto expectedne = std::vector<double>(nxne*nyne, 0);
	for(int i=0; i<nyne; i++)
		for(int j=0; j<nxne; j++)
		{
			ind = j + i*nxne;
			xj = (j+start[1])*dx;
			yi = (i+start[1])*dx;
			x0ne[ind]       = sin(k*xj)*sin(k*yi);
			expectedne[ind] = sin(k*xj)*sin(k*yi)*exp(-(2*k*k)*dt*(nt-1));
		}
	std::shared_ptr<Solver> solverne = std::make_shared<HeatSolverBTCS>(nyne, dx, nxne, dx);
	Region regionne = Region(K, overlap, nt, nyne, dx, nxne, dx, x0ne, solverne);
	regionne.hold_constant(EAST|NORTH);


	// Solve whole domain expected
	region.set_dt(dt, 0);
	region.time_step_chunk();
	auto x = region.get_x();
	
	// Iterate over 4 domains
	regionw.set_dt(dt, 0);
	regione.set_dt(dt, 0);
	regionnw.set_dt(dt, 0);
	regionne.set_dt(dt, 0);

	double error, errorw, errore, errornw, errorne;
	std::vector<double> send1, send2;
	std::vector<double> xw, xe, xnw, xne;
	for(int i=0; i<5; i++)
	{
		// Time step
		regionw.time_step_chunk();
		regione.time_step_chunk();
		regionnw.time_step_chunk();
		regionne.time_step_chunk();

		// Exchange boundarys
		send1 = regionw.get_boundary(EAST, 0);
		send2 = regione.get_boundary(WEST, 0);
		regionw.set_boundary(EAST, &send2[0], 0);
		regione.set_boundary(WEST, &send1[0], 0);

		send1 = regionw.get_boundary(NORTH, 0);
		send2 = regionnw.get_boundary(SOUTH, 0);
		regionw.set_boundary(NORTH, &send2[0], 0);
		regionnw.set_boundary(SOUTH, &send1[0], 0);

		send1 = regione.get_boundary(NORTH, 0);
		send2 = regionne.get_boundary(SOUTH, 0);
		regione.set_boundary(NORTH, &send2[0], 0);
		regionne.set_boundary(SOUTH, &send1[0], 0);

		send1 = regionnw.get_boundary(EAST, 0);
		send2 = regionne.get_boundary(WEST, 0);
		regionnw.set_boundary(EAST, &send2[0], 0);
		regionne.set_boundary(WEST, &send1[0], 0);

		// Check error
		// I left these in the loop just in case I
		// need to come back and debug again
		xw = regionw.get_x();
		xe = regione.get_x();
		xnw = regionnw.get_x();
		xne = regionne.get_x();

		errorw = 0.;
		errore = 0.;
		errornw = 0.;
		errorne = 0.;
		for(int i=0; i<nyw*nxw; i++)
			errorw = std::max(errorw, std::abs(expectedw[i] - xw[i]));
		for(int i=0; i<nye*nxe; i++)
			errore = std::max(errore, std::abs(expectede[i] - xe[i]));
		for(int i=0; i<nynw*nxnw; i++)
			errornw = std::max(errornw, std::abs(expectednw[i] - xnw[i]));
		for(int i=0; i<nyne*nxne; i++)
			errorne = std::max(errorne, std::abs(expectedne[i] - xne[i]));

		// Reset for next iteration
		regionw.reset();
		regione.reset();
		regionnw.reset();
		regionne.reset();

 	}

	error = 0.;
	for(int i=0; i<N*N; i++)
		error = std::max(error, std::abs(expected[i] - x[i]));

	BOOST_CHECK_LT( error,   tol );
	BOOST_CHECK_LT( errorw,  tol );
	BOOST_CHECK_LT( errore,  tol );
	BOOST_CHECK_LT( errornw, tol );
	BOOST_CHECK_LT( errorne, tol );

}
BOOST_AUTO_TEST_SUITE_END()
