
#ifndef R_BLAS_H
#define R_BLAS_H

#include <R_ext/RS.h>		/* for F77_... */
#include <R_ext/Complex.h>	/* for Rcomplex */

#ifndef BLAS_extern
#define BLAS_extern extern
#endif

BLAS_extern void
F77_NAME(_dgemm)(const char *transa, const char *transb, const int *m,
		const int *n, const int *k, const double *alpha,
		const double *a, const int *lda,
		const double *b, const int *ldb,
		const double *beta, double *c, const int *ldc);
		
BLAS_extern void
F77_NAME(_dgemv)(const char *trans, const int *m, const int *n,
		const double *alpha, const double *a, const int *lda,
		const double *x, const int *incx, const double *beta,
		double *y, const int *incy);
#endif