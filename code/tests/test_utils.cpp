
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

	const std::size_t k = 0;
	const std::size_t N = 10;
	const std::size_t n = 3;

	std::vector<std::size_t> start;
	std::vector<std::size_t> end;
	partition_domain(start, end, N, n, k);

	std::vector<std::size_t> start_e(n);
	start_e[0] = 0;	start_e[1] = 3;	start_e[2] = 6;
	std::vector<std::size_t> end_e(n);
	end_e[0] = 4;	end_e[1] = 7;	end_e[2] = 10;

	BOOST_CHECK( std::equal(start.begin(), start.end(), start_e.begin()) &&
				 std::equal(end.begin(), end.end(), end_e.begin()) );

}

BOOST_AUTO_TEST_CASE( partition_grid_overlap )
{

	const std::size_t k = 1;
	const std::size_t N = 10;
	const std::size_t n = 3;

	std::vector<std::size_t> start;
	std::vector<std::size_t> end;
	partition_domain(start, end, N, n, k);

	std::vector<std::size_t> start_e(n);
	start_e[0] = 0;	start_e[1] = 3;	start_e[2] = 6;
	std::vector<std::size_t> end_e(n);
	end_e[0] = 5;	end_e[1] = 8;	end_e[2] = 10;

	BOOST_CHECK( std::equal(start.begin(), start.end(), start_e.begin()) &&
				 std::equal(end.begin(), end.end(), end_e.begin()) );

}


BOOST_AUTO_TEST_SUITE_END()
