#ifndef LOGLIK_AFT_H
#define LOGLIK_AFT_H

void loglik_aft(int *dis,
	       int *mb, 
	       double *b, 
	       double *alpha,
	       double *gamma,
	       int *nn, 
	       double *z, 
	       double *time0, 
	       double *time, 
	       int *ind, 
	       double *offset,
	       double *f);

void d_loglik_aft(int *dis,
		 int *mb, 
		 double *b,
		 double *alpha,
		 double *gamma,
		 int *nn, 
		 double *z, 
		 double *time0, 
		 double *time, 
		 int *ind, 
		 double *offset, 
		 double *fp);

void d2_loglik_aft(int *dis,
		  int *mb, 
		  double *b,
		  double *alpha,
		  double *gamma,
		  int *nn, 
		  double *z, 
		  double *time0, 
		  double *time, 
		  int *ind, 
		  double *offset, 
		  double *fpp);

#endif
