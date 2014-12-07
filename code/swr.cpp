
#include "swr.hpp"

/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ CProxy_SWRDomain domainProxy;

void init_node()
{
	PetscInitializeNoArguments();
}


Main::Main(CkArgMsg* m)
{

	if(m->argc < 5) CkAbort("Expects 4 args: N domains (each dimension), L, n_iter, lb_freq");

	N      = atoi(m->argv[1]);
	Lpi    = atoi(m->argv[2]);
	n_iter = atoi(m->argv[3]);
	lb_freq = atoi(m->argv[4]);
	delete m;

	mainProxy = thisProxy;

	// The last iteration to report its error
	iter_err = 0;
	// The last error reported
	err = -1;

	count = 0;
	int nx = 100*N;
	int ny = 100*N;

	int K = 1;
	int overlap = 20;
	int nt = 10;
	int ind;
	double dx, dt, L, k;

	L  = Lpi*M_PI;
	k  = 1.;
	dt = .01;
	dx = L/((double)(nx-1));

	domainProxy = CProxy_SWRDomain::ckNew(K, overlap, nt, dt, nx, dx, nx, dx, N, N, n_iter, lb_freq, N,  N);
	domainProxy.run_simulation();

}

void Main::maxErrorReduction(double error)
{
	err = error;
	ckout << "Iteration " << iter_err++ << ": " << error << endl;
}

void Main::done(int idx_x, int idx_y, int iteration)
{

	if(++count==N*N)
	{
		ckout << "Done at iteration " << iteration
			  << " - Error: " << err << endl;
		CkExit();
	}
}

SWRDomain::SWRDomain(int K_, int overlap_, int nt_, double dt_,
					 int gny_, double dy_, int gnx_, double dx_,
					 int GNx_, int GNy_, int n_iter_, int lb_freq_):
	K(K_), overlap(overlap_), nt(nt_), dt(dt_), gny(gny_), dy(dy_), gnx(gnx_), dx(dx_),
	GNx(GNx_), GNy(GNy_), n_iter(n_iter_), lb_freq(lb_freq_)
{
	usesAtSync = true;

	// Find this domains location in the global grid
	std::vector<int> start(GNx, 0);
	std::vector<int> end(GNx, 0);

	partition_domain(start, end, gnx, GNx, overlap);
	xstart = start[thisIndex.x];
	xend   = end[thisIndex.x];

	partition_domain(start, end, gny, GNy, overlap);
	ystart = start[thisIndex.y];
	yend   = end[thisIndex.y];

	nx = xend-xstart;
	ny = yend-ystart;

	// Setup test problem
	build_x0_expected();

	// Build data structures
	double scale = 1.;
	if(thisIndex.x<std::floor(GNx/2))
		scale = 2.;
	nt = nt*scale;
	dt = dt/scale;
	interp = std::vector<double>(nt*std::max(nx, ny), 0);

	Solver *solver = new HeatSolverBTCS(ny, dy, nx, dx);
	region = new Region(K, overlap, nt, ny, dy, nx, dx, x0, solver);
	region->set_dt(dt, 0);

	// Setup boundary comm logic
	n_recv = 0;
	int hold_constant = 0;
	if(thisIndex.x>0)
	{
		west = thisIndex.x-1;
		comm_west = true;
		++n_recv;
	}
	else
	{
		west = -1;
		comm_west = false;
		hold_constant |= WEST;
	}

	if(thisIndex.x<GNx-1)
	{
		east = thisIndex.x+1;
		comm_east = true;
		++n_recv;
	}
	else
	{
		east = -1;
		comm_east = false;
		hold_constant |= EAST;
	}

	if(thisIndex.y>0)
	{
		south = thisIndex.y-1;
		comm_south = true;
		++n_recv;
	}
	else
	{
		south = -1;
		comm_south = false;
		hold_constant |= SOUTH;
	}

	if(thisIndex.y<GNy-1)
	{
		north = thisIndex.y+1;
		comm_north = true;
		++n_recv;
	}
	else
	{
		north = -1;
		comm_north = false;
		hold_constant |= NORTH;
	}

	region->hold_constant(hold_constant);

}

SWRDomain::~SWRDomain()
{
}
void SWRDomain::ResumeFromSync() {
	run_iter();
}
void SWRDomain::pup(PUP::er& p){
	CBase_SWRDomain::pup(p);
	p | nt; p | gny; p | gnx;
	p | dt; p | dy; p | dx;
	p | ny; p | nx; p | overlap;
	p | GNx; p | GNy;
	p | K;
	p | west; p | east; p | north; p | south;
	p | comm_west; p | comm_east; p | comm_north; p | comm_south;
	p | xstart; p | xend; p | ystart; p | yend;
	p | x0; p | expected;
	p | iteration; p | n_recv; p | recv;
	p | n_iter; p | lb_freq;
	if (p.isUnpacking()) {
		region = new Region();
		p | (*region);
		region->reset();
		interp = std::vector<double>(nt*std::max(nx, ny), 0);
	}
	else {
		p | *region;
	}

}

void SWRDomain::build_x0_expected()
{

	x0       = std::vector<double>(ny*nx, 0);
	expected = std::vector<double>(ny*nx, 0);

	double k = 1.;
	double xj, yi;
	int ind;
	for(int i=0; i<ny; i++)
		for(int j=0; j<nx; j++)
		{
			ind = j + i*nx;
			xj  = (j+xstart)*dx;
			yi  = (i+ystart)*dy;
			x0[ind]       = sin(k*xj)*sin(k*yi);
			expected[ind] = sin(k*xj)*sin(k*yi)*exp(-(2*k*k)*dt*(nt-1));
		}

}



#include "SWRCharm.def.h"
