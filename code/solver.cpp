
#include "solver.hpp"
#include <iostream>

#ifdef __CHARMC__
void operator| (PUP::er& p, Vec& v) {
  PetscInt sz;
  if (!p.isUnpacking()) {
    VecGetSize(v, &sz);
  }
  p | sz;
  if(p.isUnpacking()) {    
	  VecCreateSeq(PETSC_COMM_WORLD, sz,  &v);
	  VecSetFromOptions(v);
	  VecGetSize(v, &sz);

	  for (int i = 0; i < sz; i++) {
		  PetscScalar d;
		  p | d;
		  VecSetValue(v, i, d, INSERT_VALUES);
	  }
	  VecAssemblyBegin(v);
	  VecAssemblyEnd(v);
  }
  else {
	  for (int i = 0; i < sz; i++) {
		  PetscScalar d;
		  VecGetValues(v, 1, &i, &d);
		  p | d;
	  }
  }
}
#endif


double Solver::get_dt()
{
	return dt;
}

#ifdef __CHARMC__
void Solver::pup(PUP::er &p) 
{
	PUP::able::pup(p);
	p|nx; p|ny; p|dy; p|dx; p|dt;
}
#endif

void HeatSolverBTCS::set_dt(double dt_)
{
	dt = dt_;
	two_d_heat_BTCS_T_up(A, ny, dy, nx, dx, dt);
}

HeatSolverBTCS::HeatSolverBTCS(int ny_, double dy_, int nx_, double dx_):
	Solver(ny_, dy_, nx_, dx_)
{
	// Create a sparse matrix A
	two_d_heat_BTCS(A, ny, dy, nx, dx);
	dt = 1.0;
	// set up linear solver context (ksp) and preconditioner (pc)
	KSPCreate(PETSC_COMM_WORLD,&ksp);
	KSPSetOperators(ksp,A,A);
	KSPGetPC(ksp,&pc);
	PCSetType(pc,PCLU);
	KSPSetFromOptions(ksp);

	// create rhs
	VecCreateSeq(PETSC_COMM_SELF, nx*ny, &rhs);
	VecSetFromOptions(rhs);
	VecSet(rhs,0.0);

	// prepare temp
	VecCreateSeq(PETSC_COMM_SELF, nx*ny, &temp);
	VecSetFromOptions(temp);
}

HeatSolverBTCS::~HeatSolverBTCS()
{
	VecDestroy(&rhs);
	MatDestroy(&A);
}

/*! Solve for x in Ax=b
 */
void HeatSolverBTCS::solve(std::vector<double>& x)
{
	PetscInt i;
	PetscScalar s;
	KSPSolve(ksp, rhs, temp);
	KSPGetConvergedReason(ksp, &reason);
	for (i = 0; i < nx*ny; i++) {
		VecGetValues(temp, 1, &i, &s);
		x[i] = s;
	}
}

// It would probablly be better to do this with iterators. We can
// look into that later.
void HeatSolverBTCS::set_rhs(const std::vector<double>& b,
							 double* west, double* east,
							 double* north, double* south)
{
	// This method assumes that rhs is ordered so that it begins  
	// with the south border proceeding from west to east
	PetscScalar val;
	int j;
	PetscInt sz;
	VecGetSize(rhs, &sz);
	// Set the south boundary values.
	for (int i = 0; i < nx; i++)
	{
		val = south[i];
		VecSetValues(rhs, 1, &i, &val, INSERT_VALUES);
	} 

	// Set interior boundary values
	for (int i=nx; i<nx*(ny-1); i++)
	{
		val = b[i];
		VecSetValues(rhs, 1, &i, &val, INSERT_VALUES);
	}

	// Set north boundary values
	for(int i=0; i<nx; i++)
	{
		j = i+nx*(ny-1);
		val = north[i];
		VecSetValues(rhs, 1, &j, &val, INSERT_VALUES);
	}

	// Set the east and west boundary values
	for (int i=0; i<ny; i++)
	{
		j = i*nx;
		val = west[i];
		VecSetValues(rhs, 1, &j, &val, INSERT_VALUES);

		j = (i+1)*nx-1;
		val = east[i];
		VecSetValues(rhs, 1, &j, &val, INSERT_VALUES);
	}

}

#ifdef __CHARMC__
void HeatSolverBTCS::pup(PUP::er &p) 
{
  Solver::pup(p);  
  p|rhs; p|temp;
  if (p.isUnpacking()) {
    //unpack A
    two_d_heat_BTCS(A, ny, dy, nx, dx); 
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCLU);
    KSPSetFromOptions(ksp);
  } else {
    //pack A
  }
}
#endif

void DummySolver::set_dt(double dt_)
{
	dt = dt_;
}

DummySolver::DummySolver(int ny_, double dy_, int nx_, double dx_):
	Solver(ny_, dy_, nx_, dx_)
{
}

DummySolver::~DummySolver()
{
}

/*! Solve for x in Ax=b
 */
void DummySolver::solve(std::vector<double>& x)
{
	for(int i=0; i<x.size(); i++)
		x[i] = dt;
}

void DummySolver::set_rhs(const std::vector<double>& b,
						  double* west, double* east,
						  double* north, double* south)
{
	// Nothing to see here...
}

