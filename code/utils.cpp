
#include "utils.hpp"

void one_d_heat_BTCS(Mat &A, PetscInt n, PetscReal dx, PetscReal dt) {

	PetscReal beta = dt/(dx*dx);
	PetscInt i;
	PetscInt non_zeros_per_row = 3;
	PetscReal stencil[] = {-beta, 1.0+2.0*beta, -beta};

	// Create a sparse matrix:
	MatCreateSeqAIJ(PETSC_COMM_WORLD,n,n,non_zeros_per_row,PETSC_NULL,&A);
	MatSetFromOptions(A);

	for (i = 1; i < n-1; i++) {
	
		// set stencil destination
		PetscInt j_to_change[] = {i-1, i, i+1};
		PetscInt i_to_change[] = {i};

		// set values
		MatSetValues(A,1,i_to_change,3,j_to_change,stencil,INSERT_VALUES);
	}
	MatSetValue(A,0,0,1,INSERT_VALUES);
	MatSetValue(A,n-1,n-1,1,INSERT_VALUES);

	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); 
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
}

void two_d_heat_BTCS(Mat &A, PetscInt ny, PetscReal dy, PetscInt nx, PetscReal dx) {

	PetscReal kx,ky;
	PetscReal dt = 1.0;
	PetscInt N,i,j,m;
	PetscInt non_zeros_per_row = 5;

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
			if (i!=0 && j!=0 && i!=nx-1 && j!=ny-1) {

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
}

void two_d_heat_BTCS_T_up(Mat &A, PetscInt ny, PetscReal dy, PetscInt nx, PetscReal dx, PetscReal dt) {
	
	PetscReal kx,ky;
	PetscInt i,j,m;

	kx = dt/(dx*dx);
	ky = dt/(dy*dy);

	PetscReal stencil[] = {-1*ky, -1*kx, 1+2*kx+2*ky, -1*kx, -1*ky};

	// iterate only over rows that must be changed, apply stencil to each row
	for (j = 1; j < ny-1; j++) {		// executes ny-2 times
		for (i = 1; i < nx-1; i++) {	// executes nx-2 times
		
			// indices of diagonal
			m = i+j*nx;
	
			// set stencil destination
			PetscInt cols_to_change[] = {m-nx, m-1, m, m+1, m+nx};
			PetscInt rows_to_change[] = {m};

			// set values
			MatSetValues(A,1,rows_to_change,5,cols_to_change,stencil,INSERT_VALUES);
		}
	}

	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); 
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
}

void partition_domain(std::vector<int>& start, std::vector<int>& end,
					  const int N, const int n, const int k)
{

	start.resize(n);
	end.resize(n);

	const int domain_size = round(((double)(N+(n-1)*(1+k)))/((double)n));

	int last_end = k;
	for(int i=0; i<n; i++)
	{
		start[i] = last_end-k;
		end[i]   = last_end-k+domain_size;
		last_end = last_end-k+domain_size-1;
	}

	end[n-1] = N;

}
