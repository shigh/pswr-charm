
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
#include "../interpolator.hpp"


BOOST_AUTO_TEST_SUITE( interpolator_tests )

BOOST_AUTO_TEST_CASE( interpolator_expand )
{
  const int nx_src = 5;
  const int nx_dst = 8;
  const int nt_src = 7;
  const int nt_dst = 23;
  
  vector<double> src(nx_src * nt_src);
  vector<double> dst(nx_dst * nt_dst);
  
  for (int i = 0; i < nt_src; i++) {
    for (int j = 0; j < nx_src; j++) {
      src.at(i * nx_src + j) = (double)j/(nx_src - 1) + (double)i/(nt_src - 1);
    }
  }
  
  interpolate(nt_dst, nx_dst, dst, nt_src, nx_src, src);



  for (int i = 0; i < nt_dst; i++) {   
    for (int j = 0; j < nx_dst; j++) {
      if(fabs(dst.at(i * nx_dst + j) - ((double)j/(nx_dst - 1) + (double)i/(nt_dst - 1))) > 1e-5) {
	cout << dst.at(i * nx_dst + j) << " != " << (double)j/(nx_dst - 1) + (double)i/(nt_dst - 1)<< endl;
	cout << dst.at(i * nx_dst + j)  - ((double)j/(nx_dst - 1) + (double)i/(nt_dst - 1))<< endl;
	assert(0);
      }
    }
  }
}
BOOST_AUTO_TEST_SUITE_END()
