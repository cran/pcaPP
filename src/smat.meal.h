/*
    SMat - Simple Matrix Classes v0.1beta
    Copyright (C) 2011 by Heinrich Fritz (heinrich_fritz@hotmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

//	smat.eal.h
//	Environment Abstraction Layer

	void meal_dgeev(const char* jobvl, const char* jobvr, const int* n, double* a, const int* lda, double* wr, double* wi, double* vl, const int* ldvl, double* vr, const int* ldvr, double* work, const int* lwork, int* info) ;

//matmult
	void meal_dgemm (const char *transa, const char *transb, const int *m, const int *n, const int *k, const double *alpha, const double *a, const int *lda, const double *b, const int *ldb, const double *beta, double *c, const int *ldc) ;

//svd
	void meal_dgesv (const int* n, const int* nrhs, double* a, const int* lda, int* ipiv, double* b, const int* ldb, int* info) ;
	void meal_dgesvd (const char* jobu, const char* jobvt, const int* m, const int* n, double* a, const int* lda, double* s, double* u, const int* ldu, double* vt, const int* ldvt, double* work, const int* lwork, int* info) ;

//	sorting
	void meal_sort (double *d, int l) ;
	void meal_sort_order (double *, int *, int) ;
	void meal_sort_order_rev (double *d, int *o, int l) ;

//	Random stuff
	void meal_PutRNGstate () ;
	void meal_GetRNGstate () ;

	double meal_unif_rand () ;
	double meal_norm_rand ();
	double meal_exp_rand  ();

//	void meal_runif (double *d, int l) ;
//	void meal_runif (double *d, int l, double dL, double dU) ;
//	void meal_runif_r (double *d, int l) ;
//	void meal_SampleNoReplace(int k, int n, int *y, int *x) ;

////////////////////////////////////
//	special values amd constants  //
////////////////////////////////////

	double	meal_NaN		() ;
	double	meal_PosInf	() ;
	double	meal_NegInf	() ;
	double	meal_NaReal	() ;
	int		meal_NaInt	() ;

	double  meal_PI () ;

//////////////////////////
//	printing functions  //
//////////////////////////

	void meal_printf (const char *, ...) ;
	void meal_warning (const char *) ;
	void meal_error (const char *) ;

//////////////////
//	Exceptions  //
//////////////////

	void meal_OnException (const char * szDate, const char * szFile, int nLine) ;
	void meal_OnUException () ;
