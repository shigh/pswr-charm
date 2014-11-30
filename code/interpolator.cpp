/* Local Variables: */
/* c-set-offset: 'substatement-open 0 */
/* c++-tab-always-indent: t */
/* c-basic-offset: 4 */
/* c-indent-level: 4 */
/* tab-stop-list: '(4 8 12 16 20 24 28 32 36 40 44 48 52 56 60) */
/* tab-width: 4 */
/* indent-tabs-mode: t */
/* End: */


#include "interpolator.hpp"
#include <iostream>
#include <cmath>

void interpolate(int nt1, int nx1, vector<double>& vec1, int nt2, int nx2, vector<double>& vec2) {  
  for (int t = 0; t < nt1; t++) {

    double tloc = (double)t / (double)(nt1 - 1)   * (double)(nt2 - 1);
    double t1 = floor(tloc);
    double t2 = ceil(tloc);

    for ( int x = 0; x < nx1; x++) {
      double xloc = (double)x  / (double)(nx1 - 1) * (double)(nx2 - 1);
      double x1 = floor(xloc);
      double x2 = ceil(xloc);

      double q11, q12, q21, q22;

      q11 = vec2.at(nx2 * t1 + x1);     
      
      q21 = vec2.at(nx2 * t1 + x2);
      
      q12 = vec2.at(nx2 * t2 + x1);
      
      q22 = vec2.at(nx2 * t2 + x2);
      
      vec1.at(x + nx1 * t) = 
	q11 * ((x1 + 1) - xloc) * ((t1 + 1) - tloc) 
	+ q21 * (xloc - x1) * ((t1 + 1) - tloc)
	+ q12 * ((x1 + 1) - xloc) * (tloc - t1)
	+ q12 * (xloc - x1) * (tloc - t1);
    }
  }
}


