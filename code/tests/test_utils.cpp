
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <math.h>
#include <algorithm>
#include <cstddef>
#include <vector>
#include "../utils.hpp"


BOOST_AUTO_TEST_SUITE( util_tests )

// Check for explosions
BOOST_AUTO_TEST_CASE( partition_grid_no_overlap )
{

	const int k = 0;
	const int N = 10;
	const int n = 3;

	std::vector<int> start;
	std::vector<int> end;
	partition_domain(start, end, N, n, k);

	std::vector<int> start_e(n);
	start_e[0] = 0;	start_e[1] = 3;	start_e[2] = 6;
	std::vector<int> end_e(n);
	end_e[0] = 4;	end_e[1] = 7;	end_e[2] = 10;

	BOOST_CHECK( std::equal(start.begin(), start.end(), start_e.begin()) &&
				 std::equal(end.begin(), end.end(), end_e.begin()) );

}

BOOST_AUTO_TEST_CASE( partition_grid_overlap )
{

	const int k = 1;
	const int N = 10;
	const int n = 3;

	std::vector<int> start;
	std::vector<int> end;
	partition_domain(start, end, N, n, k);

	std::vector<int> start_e(n);
	start_e[0] = 0;	start_e[1] = 3;	start_e[2] = 6;
	std::vector<int> end_e(n);
	end_e[0] = 5;	end_e[1] = 8;	end_e[2] = 10;

	BOOST_CHECK( std::equal(start.begin(), start.end(), start_e.begin()) &&
				 std::equal(end.begin(), end.end(), end_e.begin()) );

}

BOOST_AUTO_TEST_SUITE_END()
