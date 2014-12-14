
#pragma once

#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include "pup_stl.h"
#include "region.hpp"
#include "solver.hpp"
#include "PSWRCharm.decl.h"


class Main: public CBase_Main
{
private:

	int N, Lpi, n_iter, count, lb_freq, iter_err;
	double err;

public:

    Main(CkArgMsg* m);

	void maxErrorReduction(double error);

	// Termination callback
	void done(int idx_x, int idx_y, int iteration, double err);

};


class PSWRDomain: public CBase_PSWRDomain
{

private:

	PSWRDomain_SDAG_CODE

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

	int next_itr, prev_itr, curr_chunk;

	double err;

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

	PSWRDomain(CkMigrateMessage* m) {}
	void pup(PUP::er& p);       
	PSWRDomain(int K_, int overlap_, int nt_, double dt_,
			  int gny_, double dy_, int gnx_, double dx_,
			  int GNx_, int GNy_, int n_iter_, int lb_freq_);

	// Setup the test problem
	void build_x0_expected();
	virtual void ResumeFromSync();
	~PSWRDomain();

};


