/* Local Variables: */
/* c-set-offset: 'substatement-open 0 */
/* c++-tab-always-indent: t */
/* c-basic-offset: 4 */
/* c-indent-level: 4 */
/* tab-stop-list: '(4 8 12 16 20 24 28 32 36 40 44 48 52 56 60) */
/* tab-width: 4 */
/* indent-tabs-mode: t */
/* End: */


#pragma once

#include <vector>
#include <memory>
using namespace std;


//interpolates from vec2 into vec1
void interpolate(int nt1, int nx1, vector<double>& vec1, int nt2, int nx2, vector<double>& vec2);
