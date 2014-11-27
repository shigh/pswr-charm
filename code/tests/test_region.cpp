
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

BOOST_AUTO_TEST_CASE( region_constructor )
{
	// Check for explosions
	int N = 10;
	int K = 5;
	int nt = 20;
	double d = .1;

	auto x0 = std::vector<double>(N*N, 0);
	std::shared_ptr<Solver> solver = std::make_shared<DummySolver>(N, d, N, d);

	Region region = Region(K, 0, nt, N, d, N, d, x0, solver);

}

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
	Region region = Region(K, 0, nt, N, d, N, d, x0, solver);

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
	Region region = Region(K, 0, nt, N, d, N, d, x0, solver);

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

BOOST_AUTO_TEST_CASE( region_iteration )
{
	// Using a general solution of the form:
	// u(x, y, t) = sin(kx)sin(ky)exp(-2k^2 t)
	// on a square domain
	// Make sure that the region class correctly handles
	// the solver over many time steps

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
	std::cout << start[0] << ' ' << end[0] << ' '
			  << start[1] << ' ' << end[1] << std::endl;

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

	// West/East domains
	const int nxw = end[0] - start[0];
	auto x0w       = std::vector<double>(N*nxw, 0);
	auto expectedw = std::vector<double>(N*nxw, 0);
	for(int i=0; i<N; i++)
		for(int j=0; j<end[0]; j++)
		{
			ind = j + i*nxw;
			xj = j*dx;
			yi = i*dx;
			x0w[ind]       = sin(k*xj)*sin(k*yi);
			expectedw[ind] = sin(k*xj)*sin(k*yi)*exp(-(2*k*k)*dt*(nt-1));
		}

	std::shared_ptr<Solver> solverw = std::make_shared<HeatSolverBTCS>(N, dx, nxw, dx);
	Region regionw = Region(K, overlap, nt, N, dx, nxw, dx, x0w, solverw);
	regionw.hold_constant(WEST|NORTH|SOUTH);
	
	const int nxe = end[1] - start[1];
	auto x0e       = std::vector<double>(N*nxe, 0);
	auto expectede = std::vector<double>(N*nxe, 0);
	for(int i=0; i<N; i++)
		for(int j=0; j<nxe; j++)
		{
			ind = j + i*nxe;
			xj = (j+start[1])*dx;
			yi = i*dx;
			x0e[ind]       = sin(k*xj)*sin(k*yi);
			expectede[ind] = sin(k*xj)*sin(k*yi)*exp(-(2*k*k)*dt*(nt-1));
		}

	std::shared_ptr<Solver> solvere = std::make_shared<HeatSolverBTCS>(N, dx, nxe, dx);
	Region regione = Region(K, overlap, nt, N, dx, nxe, dx, x0e, solvere);
	regione.hold_constant(EAST|NORTH|SOUTH);


	// Solve expected
	region.set_dt(dt, 0);
	region.time_step_chunk();
	auto x = region.get_x();
	

	regionw.set_dt(dt, 0);
	regione.set_dt(dt, 0);

	double error  = 0.;
	double errorw = 0.;
	double errore = 0.;
	std::vector<double> sendr, sendl;
	std::vector<double> xw, xe;
	for(int i=0; i<5; i++)
	{
		regionw.time_step_chunk();
		regione.time_step_chunk();

		sendr = regionw.get_boundary(EAST, 0);
		sendl = regione.get_boundary(WEST, 0);

		regionw.set_boundary(EAST, &sendl[0], 0);
		regione.set_boundary(WEST, &sendr[0], 0);

		xw = regionw.get_x();
		xe = regione.get_x();

		errorw = 0.;
		errore = 0.;
		for(int i=0; i<N*nxw; i++)
			errorw = std::max(errorw, std::abs(expectedw[i] - xw[i]));
		for(int i=0; i<N*nxe; i++)
			errore = std::max(errore, std::abs(expectede[i] - xe[i]));

		std::cout << errorw << ' ' << errore << std::endl;

		regionw.reset();
		regione.reset();

 	}
	
	for(int i=0; i<N*N; i++)
		error = std::max(error, std::abs(expected[i] - x[i]));
	std::cout << error << std::endl;

	std::ofstream ol;
	ol.open("ol.txt");
	for(int i=0; i<N*nxw; i++)
		ol << xw[i] << ' ';
	ol << std::endl;
	ol.close();

	std::ofstream orr;
	orr.open("or.txt");
	for(int i=0; i<N*nxe; i++)
		orr << xe[i] << ' ';
	orr << std::endl;
	orr.close();

	BOOST_CHECK_LT( error, tol );

}

BOOST_AUTO_TEST_SUITE_END()
