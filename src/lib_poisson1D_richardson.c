/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void eig_poisson1D(double *eigval, int *la) {
	for (int i = 0; i <= (*la); ++i) {
		eigval[i] = 2 * (1 - cos(M_PI * (i + 1) / (*la) + 1));
	}
}

double eigmax_poisson1D(int *la) {
	return 2 * (1 - cos(M_PI * (*la) / (*la) + 1));
}

double eigmin_poisson1D(int *la) {
	return 2 * (1 - cos(M_PI * 1 / (*la) + 1));
}

double richardson_alpha_opt(int *la) {
	return 2 / (eigmax_poisson1D(la) + eigmin_poisson1D(la));
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
	double *b = calloc((*la), sizeof(double));
	double norm_RHS = cblas_dnrm2((*la), RHS, 1);

	for ((*nbite) = 0; (*nbite) < (*maxit); ++(*nbite)) {
		cblas_daxpy((*la), (*alpha_rich), b, 1, X, 1);	// x(k+1)
		cblas_dcopy((*la), RHS, 1, b, 1);
		cblas_dgbmv(CblasColMajor, CblasNoTrans, (*la), (*la), (*kl), (*ku), -1.0, AB, (*lab), X, 1.0, 1.0, b, 1.0);	// b = b - Ax

		resvec[*nbite] = cblas_dnrm2((*la), b, 1) / norm_RHS;	 // residu

		if ((*tol) >= resvec[(*nbite)]) {
			break;
		}
	}

	free(b);
}

// void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
// 	double *b = malloc((*la) * sizeof(double));
// 	double norm_RHS = cblas_dnrm2((*la), RHS, 1);

// 	for ((*nbite); (*nbite) < (*maxit); ++(*nbite)) {
// 		cblas_dcopy((*la), RHS, 1, b, 1);
// 		cblas_dgbmv(CblasColMajor, CblasNoTrans, (*la), (*la), (*kl), (*ku), (*alpha_rich), AB, (*lab), X, 1, 1.0, b, 1);	// b = b - Ax

// 		resvec[(*nbite)] = cblas_dnrm2((*la), b, 1) / norm_RHS;	 // residu
// 		cblas_daxpy((*la), (*alpha_rich), b, 1, X, 1);			 // x(k+1)

// 		if ((*tol) < resvec[(*nbite)]) {
// 			break;
// 		}
// 	}

// 	free(b);
// }

// void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite)
// {
//   double* Y = (double *) calloc(*la, sizeof(double));
//   const double d_norm_b = (1 / cblas_dnrm2(*la, RHS, 1)); // division de norme de b

//   // copy de b dans y
//   cblas_dcopy(*la, RHS, 1.0, Y, 1.0);
//   // y = y - Ax
//   cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1.0, 1.0, Y,  1.0);

//   // calcul residu
//   double residu = cblas_dnrm2(*la, Y, 1) * d_norm_b;
//   resvec[*nbite] = residu;
//   while (residu > *tol && *maxit > *nbite)
//   {
//     // calcule de x(k+1) -> X = X + alphaY
//     cblas_daxpy(*la, *alpha_rich, Y, 1.0, X, 1.0);
//     // résidu
//     cblas_dcopy(*la, RHS, 1.0, Y, 1.0);
//     cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1.0, 1.0, Y, 1.0);
//     residu = cblas_dnrm2(*la, Y, 1) * d_norm_b;
//     ++*nbite;
//     resvec[*nbite] = residu;
//     printf("resvec %d %f\n", *nbite, resvec[*nbite]);
//   }
// }

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv) {
	int index = 0;

	for (int i = 0; i < (*la); ++i) {
		MB[index + (*kv) + 1] = AB[index + (*kv) + 1];
		index += (*lab);
	}
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv) {
	int index = 0;

	for (int i = 0; i < (*la); ++i) {
		MB[index + (*kv) + 1] = AB[index + (*kv) + 1];
		MB[index + (*kv) + 2] = -(AB[index + (*kv) + 2]);
		index += (*lab);
	}
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
	// double *b = malloc((*la) * sizeof(double));
	// double norm_RHS = cblas_dnrm2((*la), RHS, 1);

	// for ((*nbite); (*nbite) < (*maxit); ++(*nbite)) {
	// 	cblas_dcopy((*la), RHS, 1, b, 1);
	// 	cblas_dgbmv(CblasColMajor, CblasNoTrans, (*la), (*lab), (*kl), (*ku), (*alpha_rich), AB, (*lab), X, 1, 1.0, b, 1);	// b = b - Ax

	// 	resvec[(*nbite)] = cblas_dnrm2((*la), b, 1) / norm_RHS;	 // residu
	// 	cblas_daxpy((*la), (*alpha_rich), b, 1, X, 1);			 // x(k+1)

	// 	if ((*tol) < resvec[(*nbite)]) {
	// 		break;
	// 	}
	// }

	// free(b);
}
