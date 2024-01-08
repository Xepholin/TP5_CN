/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void eig_poisson1D(double *eigval, int *la) {
	for (int i = 0; i < (*la); ++i) {
		eigval[i] = 2.0 - 2.0 * cos((M_PI * (i + 1.0)) / ((*la) + 1.0));
	}
}

double eigmax_poisson1D(int *la) {
	return 2.0 - 2.0 * cos((M_PI * (*la)) / ((*la) + 1.0));
}

double eigmin_poisson1D(int *la) {
	return 2.0 - 2.0 * cos((M_PI * 1.0) / ((*la) + 1));
}

double richardson_alpha_opt(int *la) {
	return 2 / (eigmax_poisson1D(la) + eigmin_poisson1D(la));
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
	double *b = calloc((*la), sizeof(double));
	double norm_RHS = cblas_dnrm2((*la), RHS, 1);

	for ((*nbite) = 0; (*nbite) < (*maxit); ++(*nbite)) {
		cblas_dcopy((*la), RHS, 1, b, 1);
		cblas_dgbmv(CblasColMajor, CblasNoTrans, (*la), (*la), (*kl), (*ku), -1.0, AB, (*lab), X, 1.0, 1.0, b, 1.0);

		resvec[*nbite] = cblas_dnrm2((*la), b, 1) / norm_RHS;  // residu
		cblas_daxpy((*la), (*alpha_rich), b, 1, X, 1);

		if ((*tol) >= resvec[(*nbite)]) {
			break;
		}
	}

	free(b);

	printf("\nNombre d'itérations : %d\n", *nbite);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv) {
	int index = 0;

	for (int i = 0; i < (*la); ++i) {
		MB[index + 1] = AB[index + 1];
		index += (*lab);
	}
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv) {
	int index = 0;

	for (int i = 0; i < (*la); ++i) {
		MB[index + 1] = AB[index + 1];
		MB[index + 2] = AB[index + 2];
		index += (*lab);
	}
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
	int ku_minus = *ku - 1;
	int info = 0;
	int NRHS = 1;

	int *ipiv = calloc((*la), sizeof(int));
	double *B = calloc((*la), sizeof(double));
	double norm_RHS = cblas_dnrm2((*la), RHS, 1);

	dgbtrf_(la, la, kl, &ku_minus, MB, lab, ipiv, &info);

	for ((*nbite) = 0; (*nbite) < *maxit; ++(*nbite)) {
		cblas_dcopy(*la, RHS, 1, B, 1);
		cblas_dgbmv(CblasColMajor, CblasNoTrans, (*la), (*la), (*kl), (*ku), -1.0, AB, (*lab), X, 1, 1.0, B, 1);

		resvec[*nbite] = cblas_dnrm2(*la, B, 1) / norm_RHS;
		dgbtrs_("N", la, kl, &ku_minus, &NRHS, MB, lab, ipiv, B, la, &info);
		cblas_daxpy((*la), 1, B, 1, X, 1);

		if ((*tol) >= resvec[(*nbite)]) {
			break;
		}
	}
	free(B);
	free(ipiv);

	printf("\nNombre d'itérations : %d\n", *nbite);
}