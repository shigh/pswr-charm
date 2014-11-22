
#include "solver.hpp"
#include <iostream>

Solver::Solver(int ny_, double dy_, int nx_, double dx_):
	ny(ny_), nx(nx_), dy(dy_), dx(dx_)
{
	// Create a sparse matrix A
	two_d_heat_BTCS(A, 1.0, ny, dy, nx, dx, true);

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

Solver::~Solver()
{
	VecDestroy(&rhs);
	MatDestroy(&A);
}

/*! Solve for x in Ax=b
 */
void Solver::solve(std::vector<double>& x)
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
void Solver::set_rhs(const std::vector<double>& b,
					 double* west, double* east,
					 double* north, double* south)
{
	// This method assumes that rhs is ordered so that it begins  
	// with the south border proceeding from west to east
	PetscScalar val;
	int j;

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

void Solver::set_dt(double dt_)
{
	// Recreating A just to update dt is inefficient. We will do this differently later.
	two_d_heat_BTCS(A, dt_, ny, dy, nx, dx, false);
}

double Solver::get_dt()
{
	return dt;
}

