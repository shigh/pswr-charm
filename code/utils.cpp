
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

void two_d_heat_BTCS(Mat &A, PetscReal dt, PetscInt ny, PetscReal dy, PetscInt nx, PetscReal dx, bool init) {

	PetscReal kx,ky;
	PetscInt N,i,j,m;
	PetscInt non_zeros_per_row = 5;

	kx = dt/(dx*dx);
	ky = dt/(dy*dy);
	N = nx*ny;

	// This will be removed after we fix our dt update
	if (init) {
		// Create a sparse matrix:
		MatCreateSeqAIJ(PETSC_COMM_WORLD,N,N,non_zeros_per_row,PETSC_NULL,&A);
		MatSetFromOptions(A);
	}

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
}

void partition_domain(std::vector<std::size_t>& start, std::vector<std::size_t>& end,
					  const std::size_t N, const std::size_t n, const std::size_t k)
{

	start.resize(n);
	end.resize(n);

	const std::size_t domain_size = round(((double)(N+(n-1)*(1+k)))/((double)n));

	std::size_t last_end = k;
	for(std::size_t i=0; i<n; i++)
	{
		start[i] = last_end-k;
		end[i]   = last_end-k+domain_size;
		last_end = last_end-k+domain_size-1;
	}

	end[n-1] = N;

}
