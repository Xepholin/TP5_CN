/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

void set_GB_operator_colMajor_poisson1D(double *AB, int *lab, int *la, int *kv) {
	int index = 0;

	for (int i = 0; i < (*la); ++i) {
		if (i == 0) {
			AB[index + (*kv)] = 0.0;
			AB[index + (*kv) + 1] = 2.0;
			AB[index + (*kv) + 2] = -1.0;
		} else if (i == (*la) - 1) {
			AB[index + (*kv)] = -1.0;
			AB[index + (*kv) + 1] = 2.0;
			AB[index + (*kv) + 2] = 0.0;
		} else {
			AB[index + (*kv)] = -1.0;
			AB[index + (*kv) + 1] = 2.0;
			AB[index + (*kv) + 2] = -1.0;
		}

		index += (*lab);
	}
}

void set_GB_operator_colMajor_poisson1D_Id(double *AB, int *lab, int *la, int *kv) {
	int index = 0;

	for (int i = 0; i < (*la); ++i) {
		AB[index + (*kv) + 1] = 1.0;

		index += (*lab);
	}
}

void set_dense_RHS_DBC_1D(double *RHS, int *la, double *BC0, double *BC1) {
	RHS[0] = *BC0;

	for (int i = 1; i < (*la) - 1; ++i) {
		RHS[i] = 0.0;
	}

	RHS[*la - 1] = *BC1;
}

void set_analytical_solution_DBC_1D(double *EX_SOL, double *X, int *la, double *BC0, double *BC1) {
	for (int i = 0; i < (*la); ++i) {
		EX_SOL[i] = (*BC0) + X[i] * ((*BC1) - (*BC0));
	}
}

void set_grid_points_1D(double *x, int *la) {
	double h = 1.0 / ((*la) + 1.0);
	x[0] = h;
	for (int i = 1; i < (*la); ++i) {
		x[i] = x[i - 1] + h;
	}
}

double relative_forward_error(double *x, double *y, int *la) {
	cblas_daxpy((*la), -1.0, y, 1, x, 1);

	double abs_err = cblas_dnrm2((*la), x, 1);
	double y_norm = cblas_dnrm2((*la), y, 1);

	return abs_err / y_norm;
}

int indexABCol(int i, int j, int *lab) {
	return j * (*lab) + i;
}

int dgbtrftridiag(int *la, int *n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info) {
	ipiv[0] = 1;
	*info = 0;

	for (int i = 1; i < (*la); i++) {
		if (AB[((*lab) * i) - 2] == 0.0) {
			*info = i;
			return *info;
		}

		AB[((*lab) * i) - 1] = AB[((*lab) * i) - 1] / AB[((*lab) * i) - 2];
		AB[((*lab) * (i + 1)) - 2] = AB[((*lab) * (i + 1)) - 2] - AB[((*lab) * i) - 1] * AB[((*lab) * (i + 1)) - 3];

		ipiv[i] = i + 1;
	}

	return *info;
}
