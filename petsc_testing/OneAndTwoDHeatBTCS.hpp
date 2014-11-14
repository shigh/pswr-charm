
#pragma once

#include<petscmat.h>

void one_d_heat_BTCS(Mat &A, PetscInt n, PetscReal dx, PetscReal dt);
void two_d_heat_BTCS(Mat &A, PetscReal dt, PetscInt ny, PetscReal dy, PetscInt nx, PetscReal dx, bool init);




