
#include<petscmat.h>

Mat one_d_heat_BTCS(PetscInt n, PetscReal dx, PetscReal dt);
Mat two_d_heat_BTCS(PetscReal dt, PetscInt ny, PetscReal dy, PetscInt nx, PetscReal dx);

int main( int argc, char *argv[] ) {

	PetscReal dt,dx,dy;
	PetscInt nt,nx,ny;

	// values set small for testing...
	nt = 5; // nt = 100;
	nx = 5; // nx = 100;
	ny = 5; // ny = 100;
	dt = 1.0/(nt-1);
	dx = 2.0*PETSC_PI/(nx-1);
	dy = 2.0*PETSC_PI/(ny-1);

	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

//	Mat A = two_d_heat_BTCS(dt,ny,dy,nx,dx);
	Mat A = one_d_heat_BTCS(nx,dx,dt);

	MatView(A,PETSC_VIEWER_STDOUT_WORLD);  

	MatDestroy(&A); 
	PetscFinalize();
	return 0;
}

Mat one_d_heat_BTCS(PetscInt n, PetscReal dx, PetscReal dt) {

	Mat A;
	PetscReal beta = dt/(dx*dx);
	PetscInt i;
	PetscInt non_zeros_per_row = 3; //TODO Scott, is this always 3?
	PetscReal stencil[] = {-beta, 1.0+2.0*beta, -beta};

	// Create a sparse matrix:
	MatCreateSeqAIJ(PETSC_COMM_WORLD,n,n,non_zeros_per_row,PETSC_NULL,&A);
	MatSetFromOptions(A);  //TODO figure out why MatSetFromOptions is recommended

	for (i = 1; i < n-1; i++) {
	
		// set stencil destination
		PetscInt j_to_change[] = {i-1, i, i+1};
		PetscInt i_to_change[] = {i};

		// set values
		MatSetValues(A,1,i_to_change,3,j_to_change,stencil,INSERT_VALUES);
	}
	MatSetValue(A,0,0,1,INSERT_VALUES);
	MatSetValue(A,n-1,n-1,1,INSERT_VALUES);

	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);  //TODO learn more about these 
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);  //TODO learn more about these

	return A;
}


Mat two_d_heat_BTCS(PetscReal dt, PetscInt ny, PetscReal dy, PetscInt nx, PetscReal dx) {

	Mat A;
	PetscReal kx,ky;
	PetscInt N,i,j,m;
	PetscInt non_zeros_per_row = 5; //TODO Scott, is this always 5?

	kx = dt/(dx*dx);
	ky = dt/(dy*dy);
	N = nx*ny;

	// Create a sparse matrix:
	MatCreateSeqAIJ(PETSC_COMM_WORLD,N,N,non_zeros_per_row,PETSC_NULL,&A);
	MatSetFromOptions(A);

	for (j = 0; j < ny; j++) {
		for (i = 0; i < nx; i++) {
			
			// row/col for point
			m = i+j*nx;

			// Either the boundary or the I part of the equation
			MatSetValue(A,m,m,1,INSERT_VALUES);

			// Differentials
			if (i!=0 && j!=0 && i!=ny-1 && j!=nx-1) {

				// x-direction
				MatSetValue(A,m,m,1+2*kx+2*ky,INSERT_VALUES);
				MatSetValue(A,m,m+1,-1*kx,INSERT_VALUES);
				MatSetValue(A,m,m-1,-1*kx,INSERT_VALUES);

				// y-direction
				MatSetValue(A,m,m+nx,-1*ky,INSERT_VALUES);
				MatSetValue(A,m,m-nx,-1*ky,INSERT_VALUES);
			}
		}
	}

	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); 
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	 
	return A;
}


// TODO: should be checking error values from PETSC calls
// e.g.  ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr); 



