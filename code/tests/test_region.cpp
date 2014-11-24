
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
	int nt = 20;
	double d = .1;

	auto x0 = std::vector<double>(N*N, 1.);
	std::shared_ptr<Solver> solver = std::make_shared<DummySolver>(N, d, N, d);
	Region region = Region(K, 0, nt, N, d, N, d, x0, solver);


	for(int i=0; i<K; i++)
		region.set_dt(i+1, i);
		
	for(int i=0; i<K; i++)
		region.time_step_chunk();

	for(int i=0; i<K; i++)
	{
		auto bndy = region.get_boundary(EAST, i);
		for(auto i:bndy)
			std::cout << i << ' ';
		std::cout << std::endl;
		bndy = region.get_boundary(WEST, i);
		for(auto i:bndy)
			std::cout << i << ' ';
		std::cout << std::endl;
		bndy = region.get_boundary(NORTH, i);
		for(auto i:bndy)
			std::cout << i << ' ';
		std::cout << std::endl;
		bndy = region.get_boundary(SOUTH, i);
		for(auto i:bndy)
			std::cout << i << ' ';
		std::cout << std::endl;

	}

	for(int i=0; i<K; i++)
		std::cout << region.get_dt(i) << ' ';
	std::cout << std::endl;

}


BOOST_AUTO_TEST_SUITE_END()
