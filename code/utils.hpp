
#pragma once

#include <vector>
#include <assert.h>
#include <petscmat.h>

// I have never had good luck with enums
typedef int boundary_t;

#define WEST  (boundary_t)0
#define EAST  (boundary_t)1
#define NORTH (boundary_t)2
#define SOUTH (boundary_t)3

void one_d_heat_BTCS(Mat &A, PetscInt n, PetscReal dx, PetscReal dt);
void two_d_heat_BTCS(Mat &A, PetscReal dt, PetscInt ny, PetscReal dy, PetscInt nx, PetscReal dx, bool init);
void partition_domain(std::vector<std::size_t>& start, std::vector<std::size_t>& end,
					  const std::size_t N, const std::size_t n, const std::size_t k);

