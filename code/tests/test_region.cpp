
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

BOOST_AUTO_TEST_SUITE_END()
