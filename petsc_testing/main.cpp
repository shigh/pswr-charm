
#include "OneAndTwoDHeatBTCS.hpp"
#include "solver.hpp"
#include <vector>

int main( int argc, char *argv[] ) {

	PetscReal dt,dx,dy;
	PetscInt nt,nx,ny;

	// Testing PetscInitialize:
	nt = 5; // nt = 100;
	nx = 5; // nx = 100;
	ny = 5; // ny = 100;
	dt = 1.0/(nt-1);
	dx = 2.0*PETSC_PI/(nx-1);
	dy = 2.0*PETSC_PI/(ny-1);

	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
	
	Mat Atest;
	two_d_heat_BTCS(Atest,dt,ny,dy,nx,dx,true);
	MatView(Atest,PETSC_VIEWER_STDOUT_WORLD);  


	// Example solver use:

	Solver *solver = new Solver(ny, dy, nx, dx);
	std::vector<double> b(nx*ny,1);
	std::vector<double> x(nx*ny,0);
	double testArrayW[5] = {2,2,2,2,2};
	double testArrayE[5] = {3,3,3,3,3};
	double testArrayN[5] = {4,4,4,4,4};
	double testArrayS[5] = {5,5,5,5,5};
	double *west = testArrayW;
	double *east = testArrayE;
	double *north = testArrayN;
	double *south = testArrayS;

	solver->set_rhs(b, west, east, north, south);
	solver->set_dt(dt);
	solver->solve(x);

	delete solver;
	MatDestroy(&Atest); 
	PetscFinalize();
	return 0;
}
