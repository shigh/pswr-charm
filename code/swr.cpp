#include <stdlib.h>
#include <vector>
#include <string>
#include "pup_stl.h"
#include "swr.hpp"
#include "region.hpp"
#include "solver.hpp"
#include "SWRCharm.decl.h"

/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ CProxy_SWRDomain domainProxy;

class Main: public CBase_Main
{
private:

	int N, count;

public:

    Main(CkArgMsg* m)
	{

		 mainProxy = thisProxy;

		 count = 0;
		 N = 2;

		 int nx = 200;
		 int ny = 200;

		 int K = 1;
		 int overlap = 20;
		 int nt = 10;
		 int ind;
		 double dx, dt, L, k;

		 L  = 2*M_PI;
		 k  = 1.;
		 dt = .01;
		 dx = L/((double)(N-1));

		 domainProxy = CProxy_SWRDomain::ckNew(K, overlap, nt, nx, dx, nx, dx, N, N, N, N);

    }

	void callback()
	{
		if(++count==N*N)
			CkExit();
	}
	

};

class SWRDomain: public CBase_SWRDomain
{

private:

	std::shared_ptr<Region> region;


	// Global ny/nx
	int gny, gnx;
	// Number of domains
	int GNx, GNy;

	// K = number of chunks
	int K, overlap;
	int nt, ny, nx;
	double dy, dx;

	int west, east, north, south;
	int xstart, xend, ystart, yend;

public:

	SWRDomain(CkMigrateMessage* M) {}

    SWRDomain(int K_, int overlap_, int nt_, int gny_, double dy_, int gnx_, double dx_,
			  int GNx_, int GNy_):
		K(K_), overlap(overlap_), nt(nt_), gny(gny_), dy(dy_), gnx(gnx_), dx(dx_),
		GNx(GNx_), GNy(GNy_)
	{

		auto start = std::vector<int>(GNx, 0);
		auto end   = std::vector<int>(GNy, 0);

		partition_domain(start, end, gnx, GNx, overlap);
		xstart = start[thisIndex.x];
		xend   = end[thisIndex.x];

		partition_domain(start, end, gny, GNy, overlap);
		ystart = start[thisIndex.y];
		yend   = end[thisIndex.y];

		nx = xend-xstart;
		ny = yend-ystart;

		auto x0 = std::vector<double>(ny*nx, 0);
		std::shared_ptr<Solver> solver = std::make_shared<DummySolver>(ny, dy, nx, dx);
 		region = std::make_shared<Region>(K, overlap, nt, ny, dy, nx, dx, x0, solver);

		west  = thisIndex.x>0     ? thisIndex.x-1:-1;
		east  = thisIndex.x<GNx-1 ? thisIndex.x+1:-1;
		north = thisIndex.y<GNy-1 ? thisIndex.y+1:-1;
		south = thisIndex.y>0     ? thisIndex.y-1:-1;

		run_simulation();

	}

	void run_simulation()
	{

		region->set_dt(.1, 0);
		region->time_step_chunk();

		mainProxy.callback();

	}

};

#include "SWRCharm.def.h"
