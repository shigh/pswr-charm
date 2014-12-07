
#pragma once

#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include "pup_stl.h"
#include "region.hpp"
#include "solver.hpp"
#include "SWRCharm.decl.h"


class Main: public CBase_Main
{
private:

	int N, Lpi, n_iter, count, lb_freq;

public:

    Main(CkArgMsg* m);

	void maxErrorReduction(double error);

	// Termination callback
	void done();

};


class SWRDomain: public CBase_SWRDomain
{

private:

	SWRDomain_SDAG_CODE

	Region *region;

	// Global ny/ny/nx
	int nt, gny, gnx, n_iter;
	double dt, dy, dx;
	// Local ny/nx
	int ny, nx, overlap;
	// Interpolation temporary
	std::vector<double> interp;
	// Number of domains
	int GNx, GNy;
	// K = number of chunks
	int K;

	// Domain indexes
	int west, east, north, south;
	// Domains to communicate with
	bool comm_west, comm_east, comm_north, comm_south;

	// This domains start and end grid points
	int xstart, xend, ystart, yend;

	int lb_freq;
	
	// Internal tmp vectors
	std::vector<double> x0, expected;

	// Charm variables for iteration logic
	int iteration, n_recv, recv;

public:

	SWRDomain(CkMigrateMessage* m) {}
	void pup(PUP::er& p);       
	SWRDomain(int K_, int overlap_, int nt_, double dt_,
			  int gny_, double dy_, int gnx_, double dx_,
			  int GNx_, int GNy_, int n_iter_, int lb_freq_);

	// Setup the test problem
	void build_x0_expected();
	virtual void ResumeFromSync();
	~SWRDomain();

};


