/* Local Variables: */
/* c-set-offset: 'substatement-open 0 */
/* c++-tab-always-indent: t */
/* c-basic-offset: 4 */
/* c-indent-level: 4 */
/* tab-stop-list: '(4 8 12 16 20 24 28 32 36 40 44 48 52 56 60)) */
/* tab-width: 4 */
/* indent-tabs-mode: t */
/* End: */

#include <vector>
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>
#include "pup.h"
#pragma once

#include "utils.hpp"

/*
  Usage:
  Solver s = Solver(ny, dy, nx, dx);
  s.set_rhs(b, west, east, north, south, dt);
  s.solve(x);
*/

// Depending on how PETc does its solves we may need to change some of
// the data layout.  I am not sure yet what the best signature for the
// solve functions is, if you build it using these we should be able
// to modifify it relatively easily to do what we want.
class Solver: public PUP::able
{

protected:

	int 	ny, nx;
	double 	dy, dx;
	double 	dt;

public:
	PUPable_decl(Solver);
	Solver() {}
	Solver(CkMigrateMessage* msg) : PUP::able(msg) {}
	Solver(int ny_, double dy_, int nx_, double dx_):
		ny(ny_), dy(dy_), nx(nx_), dx(dx_) {};

	virtual ~Solver() {};

	/*! Solve for x in Ax=b
	 */
	virtual void solve(std::vector<double>& x) {};

	// It would probablly be better to do this with iterators. We can
	// look into that later.
	virtual void set_rhs(const std::vector<double>& b,
						 double* west, double* east,
						 double* north, double* south) {};

	virtual void set_dt(double dt_) {};

	double get_dt();

	virtual void pup(PUP::er &p);

};


class HeatSolverBTCS: public Solver
{

private:

	Vec rhs;	// right hand side
	Vec temp;	// PETSc Vec to contain result before copying to c++ vector
	Mat A;
	PC pc;		// preconditioner 
	KSP ksp;	// linear solver context
	KSPConvergedReason reason;

public:
	PUPable_decl(HeatSolverBTCS);
	HeatSolverBTCS() {}
	HeatSolverBTCS(int ny_, double dy_, int nx_, double dx_);
	HeatSolverBTCS(CkMigrateMessage* msg) : Solver(msg) {}
	void solve(std::vector<double>& x);

	void set_rhs(const std::vector<double>& b,
				 double* west, double* east,
				 double* north, double* south);

	void set_dt(double dt_);
	virtual void pup(PUP::er &p);
	~HeatSolverBTCS();
	
};

// Used for testing
// Sets x to dt in solve
class DummySolver: public Solver
{
public:
	PUPable_decl(DummySolver);
	DummySolver() {}
	DummySolver(int ny_, double dy_, int nx_, double dx_);
	DummySolver(CkMigrateMessage* msg) : Solver(msg) {}
	
	void solve(std::vector<double>& x);

	void set_rhs(const std::vector<double>& b,
				 double* west, double* east,
				 double* north, double* south);

	void set_dt(double dt_);

	~DummySolver();
	
};
