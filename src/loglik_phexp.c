#define USE_FC_LEN_T
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif
#include <math.h>
#include <Rmath.h>
#include <R.h>

#include "phfun.h"
#include "loglik_phexp.h"

extern int dist;
extern ph0S_fun *S0;
extern ph0_fun *f0;
extern ph0_fun *h0;
extern ph0_fun *f0_t;
extern ph0_fun *h0_t;
extern ph0_fun *h0_tt;

void loglik_phexp(int * dis,
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
		  double *f){

/* ********************************************************************* */
/*    Calculates minus log likelihood (in one stratum), returns in 'f'.  */
/*    gamma is fixed here (pfix == 1). Doesn't mean anythng here, but... */

    int one = 1; 
    int i;
    double zb;

    int log_p;

    double res1, res2;

    dist = *dis;

    if (dist == 0){
	S0 = &S0_weibull;
	f0 = &f0_weibull;      
	h0 = &h0_weibull;      
	f0_t = &f0_t_weibull;     
	h0_t = &h0_t_weibull;     
	h0_tt = &h0_tt_weibull;    
    }else if (dist == 1){
	S0 = &S0_loglogistic;
	f0 = &f0_loglogistic;      
	h0 = &h0_loglogistic;      
	f0_t = &f0_t_loglogistic;     
	h0_t = &h0_t_loglogistic;     
	h0_tt = &h0_tt_loglogistic;
    }else if (dist == 2){    
	S0 = &S0_lognormal;
	f0 = &f0_lognormal;      
	h0 = &h0_lognormal;      
	f0_t = &f0_t_lognormal;     
	h0_t = &h0_t_lognormal;     
	h0_tt = &h0_tt_lognormal;   
    }else if ((dist == 3) || (dist == 4)){ /* 4 = Gomperts-Makeham */    
	S0 = &S0_ev;
	f0 = &f0_ev;      
	h0 = &h0_ev;      
	f0_t = &f0_t_ev;     
	h0_t = &h0_t_ev;     
	h0_tt = &h0_tt_ev;   
    }else{
	error("Unknown distribution");
    }

    res1 = 0.0;
    res2 = 0.0;

    log_p = 1;
    for (i = 0; i < *nn; i++){
      if (*mb){
	zb = offset[i] + F77_CALL(ddot)(mb, (z + i * (*mb)), &one, b, &one);
      }else{
	zb = offset[i];
      }
	if (ind[i]){
	    res1 += log(h(time[i], (*gamma), (*alpha))) + zb;
	}
	res2 += exp(zb) * (S(time0[i], *gamma, *alpha, log_p) - 
			   S(time[i], *gamma, *alpha, log_p));
    }

    *f = -(res1 - res2);
/* NOTE: The negative of the log likelihood (min is looked for) */
}

void d_loglik_phexp(int *dis,
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
		    double *fp){
    
/* *************************************************************** */
/*    Calculates minus score (in one stratum), returns in 'fp'.     */

/* Note that here the 'dimension' is (mb + 1),                     */
/* i.e., no strata and 'ipfixed'. 'bdim' is local here and      */
/* not (necessarily) the same as the "global" bdim!                */

    int one = 1; 
    int i, j;
    int log_p = 0;
    double ezb;
    double res1, res2;
    int bdim;

    double x, x0;

    if (dist == 0){
	S0 = &S0_weibull;
	f0 = &f0_weibull;      
	h0 = &h0_weibull;      
	f0_t = &f0_t_weibull;     
	h0_t = &h0_t_weibull;     
	h0_tt = &h0_tt_weibull;    
    }else if (dist == 1){
	S0 = &S0_loglogistic;
	f0 = &f0_loglogistic;      
	h0 = &h0_loglogistic;      
	f0_t = &f0_t_loglogistic;     
	h0_t = &h0_t_loglogistic;     
	h0_tt = &h0_tt_loglogistic;
    }else if (dist == 2){    
	S0 = &S0_lognormal;
	f0 = &f0_lognormal;      
	h0 = &h0_lognormal;      
	f0_t = &f0_t_lognormal;     
	h0_t = &h0_t_lognormal;     
	h0_tt = &h0_tt_lognormal;   
    }else if ((dist == 3) || (dist == 4)){    
	S0 = &S0_ev;
	f0 = &f0_ev;      
	h0 = &h0_ev;      
	f0_t = &f0_t_ev;     
	h0_t = &h0_t_ev;     
	h0_tt = &h0_tt_ev;   
    }else{
	error("Unknown distribution");
    }

    bdim = *mb + 1;

    log_p = 1;
    for (j = 0; j < *mb; j++){
	res1 = 0.0;
	res2 = 0.0;
        for (i = 0; i < *nn; i++){
/*            if (ind[i]) res1 += z[j * mb + i]; */
	    if (ind[i]) res1 += z[i * (*mb) + j];
            ezb = exp(F77_CALL(ddot)(&(*mb), (z + i * (*mb)), &one, b, &one));
	    ezb *= exp(offset[i]);
	    res2 += z[j + i * (*mb)] * ezb *  
                (S(time0[i], *gamma, *alpha, log_p) - 
		 S(time[i], *gamma, *alpha, log_p));
	}
	fp[j] = res1 - res2;
    }

    /* Scale, lambda = exp((*alpha)) */

    res1 = 0.0;
    res2 = 0.0;

    log_p = 0;
    for (i = 0; i < *nn; i++){
	x = g(time[i], *gamma, *alpha);
	x0 = g(time0[i], *gamma, *alpha); 
	if (*mb){
	  ezb = exp(F77_CALL(ddot)(&(*mb), (z + i * (*mb)), &one, b, &one));
	}else{
	  ezb = 1.0;
	}
	res1 += ind[i] * h_alpha(time[i], (*gamma), (*alpha)) / 
	    h(time[i], (*gamma), (*alpha));
	res2 += ezb * /*
			(S_alpha(time0[i], *gamma, *alpha, log_p) / 
			S(time0[i], *gamma, *alpha, log_p) -
			S_alpha(time[i], (*gamma), (*alpha)) / 
			S(time[i], (*gamma), (*alpha)));
		      */
	    (-g_alpha(time0[i], *gamma, *alpha) *  h0(x0) +
	     g_alpha(time[i], *gamma, *alpha) *  h0(x));
    }
    
    fp[(*mb)] = res1 - res2;

    for (i = 0; i < bdim; i++){
	fp[i] = -fp[i]; /* NOTE! Minimization!! */
    }
}
	
void d2_loglik_phexp(int *dis,
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
		     double *fpp){

/* *************************************************************** */
/*    Calculates minus hessian (in one stratum), returns in 'fpp'. */
/*                                                                 */
/* Note that here the 'dimension' is ((*mb) + 1),                     */
/* i.e., no strata and 'ipfixed'. 'bdim' is local here and      */
/* not (necessarily) the same as the "global" bdim!                */

/* fpp is "((*mb) + 2) x ((*mb) + 2) = bdim x bdim)".                    */

    int i, j, m;
    double *ezb;

    double alf = 1.0; 
    double beta = 0.0;
    int one = 1;
    char trans = 'T';
    int bdim;

    int log_p;
    double tmp, s1, s2, s3;

    double x, x0;

    dist = *dis;

    if (dist == 0){
	S0 = &S0_weibull;
	f0 = &f0_weibull;      
	h0 = &h0_weibull;      
	f0_t = &f0_t_weibull;     
	h0_t = &h0_t_weibull;     
	h0_tt = &h0_tt_weibull;    
    }else if (dist == 1){
	S0 = &S0_loglogistic;
	f0 = &f0_loglogistic;      
	h0 = &h0_loglogistic;      
	f0_t = &f0_t_loglogistic;     
	h0_t = &h0_t_loglogistic;     
	h0_tt = &h0_tt_loglogistic;
    }else if (dist == 2){    
	S0 = &S0_lognormal;
	f0 = &f0_lognormal;      
	h0 = &h0_lognormal;      
	f0_t = &f0_t_lognormal;     
	h0_t = &h0_t_lognormal;     
	h0_tt = &h0_tt_lognormal;   
    }else if ((dist == 3) || (dist == 4)){    
	S0 = &S0_ev;
	f0 = &f0_ev;      
	h0 = &h0_ev;      
	f0_t = &f0_t_ev;     
	h0_t = &h0_t_ev;     
	h0_tt = &h0_tt_ev;   
    }else{
	error("Unknown distribution");
    }

    bdim = (*mb) + 1;

    ezb = R_Calloc(*nn, double);

    /* ezb[i] = exp(sum(b_j * z_ij)), i = 0, (*nn - 1). (hopefully;-) */ 
    if (*mb){
      F77_CALL(dgemv)(&trans, mb, nn, &alf, z, mb,  
		      b, &one, &beta, ezb, &one FCONE);
      for (i = 0; i < *nn; i++) ezb[i] = exp(ezb[i] + offset[i]);
    }else{
	for (i = 0; i < *nn; i++) ezb[i] = exp(offset[i]);
    }
      /*
    Rprintf("b = %f\n", b[0]);
    Rprintf("z = %f, %f, %f, %f\n", z[0], z[1], z[2], z[3]);
    Rprintf("ezb = %f, %f, %f, %f\n", ezb[0], ezb[1], ezb[2], ezb[3]);
*/
    if ((*mb) >= 1){
/* beta and beta: */
	log_p = 1;
        for (j = 0; j < (*mb); j++){
            for (m = 0; m <= j; m++){
		tmp = 0.0;
		for (i = 0; i < *nn; i++){
/*
		    ezb = exp(F77_CALL(ddot)(&(*mb), (z + i * (*mb)),
					     &one, b, &one));
*/
/*
		    Rprintf("[d2_loglik] m = %d, j = %d, i = %d, mb = %d\n",
			    m, j, i, (*mb));
		    Rprintf("[d2_loglik] ezb[%d] = %f\n", i, ezb[i]);
		    Rprintf("[d2_loglik] time0[%d] = %f\n", i, time0[i]);
		    Rprintf("[d2_loglik] time[%d] = %f\n", i, time[i]);
		    Rprintf("[d2_loglik] gamma = %f\n", i, *gamma);
		    Rprintf("[d2_loglik] alpha = %f\n", *alpha);
		    Rprintf("[d2_loglik] S0[%d] = %f\n", i, 
			    S(time0[i], *gamma, *alpha));
		    Rprintf("[d2_loglik] S[%d] = %f\n", i, 
			    S(time[i], *gamma, *alpha));

		    Rprintf("z[%d] = %f\n", z[m + i * (*mb)]);
*/
		    tmp += z[m + i * (*mb)] * z[j + i * (*mb)] * ezb[i] *
			(S(time0[i], *gamma, *alpha, log_p) -
			 S(time[i], *gamma, *alpha, log_p));
		}
		fpp[j + bdim * m] = tmp;
		fpp[m + bdim * j] = tmp;
	    }
        }


/* beta and alpha (log(scale)): */

	log_p = 0; /* From now on! */
	for (j = 0; j < (*mb); j++){
	    tmp = 0.0;
	    for (i = 0; i < *nn; i++){
		x0 = g(time0[i], *gamma, *alpha);
		x = g(time[i], *gamma, *alpha);
		s1 = /* S_alpha(time0[i], (*gamma), (*alpha)) / 
		    S(time0[i], (*gamma), (*alpha)) -
		    S_alpha(time[i], (*gamma), (*alpha)) / 
		    S(time[i], (*gamma), (*alpha));
		     */
		    -g_alpha(time0[i], *gamma, *alpha) *  h0(x0) +
		    g_alpha(time[i], *gamma, *alpha) *  h0(x);


		tmp += z[i * (*mb) + j] * ezb[i] * s1;
	    }
	    fpp[(*mb) + bdim * j] = tmp;
	    fpp[j + bdim * (*mb)] = tmp; /* Symmetry! */
	}
	
    } /* END if ((*mb) >= 1) */


/* Scale -- scale; hess[k][k]: */

    tmp = 0.0;
    for (i = 0; i < *nn; i++){
	if (ind[i]){
	    s1 =  h_alpha2(time[i], (*gamma), (*alpha)) / 
		h(time[i], (*gamma), (*alpha)) -
		R_pow_di(h_alpha(time[i], (*gamma), (*alpha)) / 
			 h(time[i], (*gamma), (*alpha)), 2);
	    tmp -= s1;
/*	    Rprintf("time = %f\n", time[i]);

	    Rprintf("h_alpha = %f\n", h_alpha(time[i], (*gamma), (*alpha)));
	    Rprintf("h_alpha2 = %f\n", h_alpha2(time[i], (*gamma), (*alpha)));
	    Rprintf("h = %f\n", h(time[i], (*gamma), (*alpha)));
	    Rprintf("s1 = %f\n\n", s1);
*/
	}
	s2 = S_alpha2(time0[i], (*gamma), (*alpha)) / 
	    S(time0[i], (*gamma), (*alpha), log_p) -
	    R_pow_di(S_alpha(time0[i], (*gamma), (*alpha)) / 
		     S(time0[i], (*gamma), (*alpha), log_p), 2);
/*	Rprintf("time0 = %f\n", time0[i]);
	Rprintf("S_alpha = %f\n", S_alpha(time0[i], (*gamma), (*alpha)));
	Rprintf("S_alpha2 = %f\n", S_alpha2(time0[i], (*gamma), (*alpha)));
	Rprintf("S = %f\n", S(time0[i], (*gamma), (*alpha), log_p));
	    Rprintf("s2 = %f\n\n", s2);
*/	
	s3 = S_alpha2(time[i], (*gamma), (*alpha)) / 
	    S(time[i], (*gamma), (*alpha), log_p) -
	    R_pow_di(S_alpha(time[i], (*gamma), (*alpha)) / 
		     S(time[i], (*gamma), (*alpha), log_p), 2);
/*	Rprintf("time = %f\n", time[i]);
	Rprintf("gamma, alpha = %f  %f\n", *gamma, *alpha);
	Rprintf("S_alpha = %f\n", S_alpha(time[i], (*gamma), (*alpha)));
	Rprintf("S_alpha2 = %f\n", S_alpha2(time[i], (*gamma), (*alpha)));
	Rprintf("S = %f\n", S(time[i], (*gamma), (*alpha), log_p));
	    Rprintf("s3 = %f\n\n\n", s3);
*/
	tmp += ezb[i] * (s2 - s3);
    }
    fpp[(*mb) + bdim * (*mb)] = tmp;
    
    
    R_Free(ezb);
}
