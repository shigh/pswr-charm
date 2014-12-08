
#pragma once

#include <vector>
#include <assert.h>
#include <petscmat.h>

// I have never had good luck with enums
typedef int boundary_t;

#define WEST  (boundary_t)0x01 // 0000 0001
#define EAST  (boundary_t)0x02 // 0000 0010
#define NORTH (boundary_t)0x04 // 0000 0100
#define SOUTH (boundary_t)0x08 // 0000 1000

void one_d_heat_BTCS(Mat &A, PetscInt n, PetscReal dx, PetscReal dt);
void two_d_heat_BTCS(Mat &A, PetscInt ny, PetscReal dy, PetscInt nx, PetscReal dx);
void two_d_heat_BTCS_T_up(Mat &A, PetscInt ny, PetscReal dy, PetscInt nx, PetscReal dx, PetscReal dt);
void partition_domain(std::vector<int>& start, std::vector<int>& end,
					  const int N, const int n, const int k);

