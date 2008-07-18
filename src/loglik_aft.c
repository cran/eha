#include <math.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R.h>

#include "phfun.h"
#include "loglik_aft.h"

extern int dist;
extern ph0S_fun *S0;
extern ph0_fun *f0;
extern ph0_fun *h0;
extern ph0_fun *f0_t;
extern ph0_fun *h0_t;
extern ph0_fun *h0_tt;

void loglik_aft(int * dis,
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


    int one = 1; 
    int i;
    double zb, alphamzb;

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
    }else if ((dist == 3) || (dist == 4)){ /* 4 = Gomperts */    
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
	    zb = F77_CALL(ddot)(mb, (z + i * (*mb)), &one, b, &one);
	}else{
	    zb = 0.0;
	}
	alphamzb = *alpha - zb;
	if (ind[i]){
	    res1 += log(h(time[i], *gamma, alphamzb));
	}
	res2 += S(time0[i], *gamma, alphamzb, log_p) - 
	    S(time[i], *gamma, alphamzb, log_p);
    }

    *f = -(res1 - res2);
/* NOTE: The negative of the log likelihood (min is looked for) */
}

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
		  double *fp){

/* *************************************************************** */
/*    Calculates minus score (in one stratum), returns in 'fp'.     */

/* Note that here the 'dimension' is (mb + 2),                     */
/* i.e., no strata and no 'ipfixed'. 'bdim' is local here and      */
/* not (necessarily) the same as the "global" bdim!                */

    int one = 1; 
    int i, j;
    int log_p = 0;
    double zb, alphamzb;
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

    bdim = *mb + 2;

    log_p = 0;
    for (j = 0; j < *mb; j++){
	res1 = 0.0;
	res2 = 0.0;
        for (i = 0; i < *nn; i++){
	    if (*mb){
		zb = F77_CALL(ddot)(&(*mb), (z + i * (*mb)), &one, b, &one);
	    }else{
		zb = 0.0;
	    }
	    alphamzb = *alpha - zb;
	    if (ind[i]) res1 -= z[i * (*mb) + j] * /* NOTE: -= */
			    h_alpha(time[i], *gamma, alphamzb) / 
			    h(time[i], *gamma, alphamzb);
	    res2 += z[j + i * (*mb)] *   
                (S_alpha(time0[i], *gamma, alphamzb) /
		 S(time0[i], *gamma, alphamzb, log_p) - 
		 S_alpha(time[i], *gamma, alphamzb) / 
		 S(time[i], *gamma, alphamzb, log_p));
	}
	fp[j] = res1 + res2;
    }

    /* Scale, lambda = exp((*alpha)) */

    res1 = 0.0;
    res2 = 0.0;

    log_p = 0;
    for (i = 0; i < *nn; i++){
	if (*mb){
	  zb = F77_CALL(ddot)(&(*mb), (z + i * (*mb)), &one, b, &one);
	}else{
	  zb = 0.0;
	}
	alphamzb = *alpha - zb;
	x = g(time[i], *gamma, alphamzb);
	x0 = g(time0[i], *gamma, alphamzb); 
	x = g(time[i], *gamma, alphamzb);
	x0 = g(time0[i], *gamma, alphamzb);
	if (ind[i])
	    res1 += h_alpha(time[i], *gamma, alphamzb) / 
		h(time[i], (*gamma), alphamzb);
	res2 +=  S_alpha(time0[i], *gamma, alphamzb) / 
	    S(time0[i], *gamma, alphamzb, log_p) -
	    S_alpha(time[i], *gamma, alphamzb) / 
	    S(time[i], *gamma, alphamzb, log_p); 
/* CHECK this! */
		 
	    /*
	      (-g_alpha(time0[i], *gamma, alphamzb) *  h0(x0) +
	     g_alpha(time[i], *gamma, alphamzb) *  h0(x));
	    */
		  }
    
    fp[(*mb)] = res1 - res2;

    /* Shape, p = exp((*gamma)) */



    res1 = 0.0;
    res2 = 0.0;
    
    for (i = 0; i < *nn; i++){
	zb = F77_CALL(ddot)(&(*mb), (z + i * (*mb)), &one, b, &one);
	alphamzb = *alpha - zb;
	x = g(time[i], *gamma, alphamzb);
	x0 = g(time0[i], *gamma, alphamzb);
	if (ind[i])
	    res1 += h_gamma(time[i], *gamma, alphamzb) / 
		h(time[i], *gamma, alphamzb);
	res2 += S_gamma(time0[i], *gamma, alphamzb) / 
	    S(time0[i], *gamma, alphamzb, log_p) -
	    S_gamma(time[i], *gamma, alphamzb) / 
	    S(time[i], *gamma, alphamzb, log_p); 
/*KOLLA!*/
	/*
	    (-g_gamma(time0[i], *gamma, alphamzb) *  h0(x0) +
	     g_gamma(time[i], *gamma, alphamzb) *  h0(x));
	*/
    }
    
    fp[(*mb) + 1] = res1 - res2;


    for (i = 0; i < bdim; i++){
	fp[i] = -fp[i]; /* NOTE! Minimization!! */
    }
}
	
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
		  double *fpp){

/* *************************************************************** */
/*    Calculates minus hessian (in one stratum), returns in 'fpp'. */
/*                                                                 */
/* Note that here the 'dimension' is ((*mb) + 2),                     */
/* i.e., no strata and no 'ipfixed'. 'bdim' is local here and      */
/* not (necessarily) the same as the "global" bdim!                */

/* fpp is "((*mb) + 2) x ((*mb) + 2) = bdim x bdim)".                    */

    int i, j, m;
    double *zb;
    double alphamzb;

    double alf = 1.0; 
    double beta = 1.0;
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

    bdim = (*mb) + 2;

    zb = Calloc(*nn, double);

    /* ezb[i] = exp(sum(b_j * z_ij)), i = 0, (*nn - 1). (hopefully;-) */ 
    F77_CALL(dcopy)(nn, offset, &one, zb, &one); 
    if (*mb){
	F77_CALL(dgemv)(&trans, mb, nn, &alf, z, mb,  
			b, &one, &beta, zb, &one);
    }

      /*
    Rprintf("b = %f\n", b[0]);
    Rprintf("z = %f, %f, %f, %f\n", z[0], z[1], z[2], z[3]);
    Rprintf("ezb = %f, %f, %f, %f\n", ezb[0], ezb[1], ezb[2], ezb[3]);
*/
    if ((*mb) >= 1){
/* beta and beta: */
	log_p = 0;
        for (j = 0; j < (*mb); j++){
            for (m = 0; m <= j; m++){
		tmp = 0.0;
		for (i = 0; i < *nn; i++){
		    alphamzb = *alpha - zb[i];
		    if (ind[i]){
			s1 = z[j + i * (*mb)] * z[m + i * (*mb)] * (  
			    h_alpha2(time[i], *gamma, alphamzb) / 
			    h(time[i], *gamma, alphamzb) -
			    R_pow_di(h_alpha(time[i], *gamma, alphamzb) / 
				     h(time[i], *gamma, alphamzb), 2));
			tmp -= s1; /* NOTE! */
		    }
		    s2 = S_alpha2(time0[i], *gamma, alphamzb) / 
			S(time0[i], *gamma, alphamzb, log_p) -
			R_pow_di(S_gamma(time0[i], *gamma, alphamzb) / 
				 S(time0[i], *gamma, alphamzb, log_p), 2);
		    
		    s3 = S_alpha2(time[i], *gamma, alphamzb) / 
			S(time[i], *gamma, alphamzb, log_p) -
			R_pow_di(S_alpha(time[i], *gamma, alphamzb) / 
				 S(time[i], *gamma, alphamzb, log_p), 2);
		    tmp += z[j + i * (*mb)] * z[m + i * (*mb)] * 
			(s2 - s3); /* NOTE: += */

		}
		fpp[j + bdim * m] = tmp;
		fpp[m + bdim * j] = tmp;
	    }
        }

/* beta and gamma (log(shape)): */
	log_p = 0;
	for (j = 0; j < (*mb); j++){
	    tmp = 0.0;
	    for (i = 0; i < *nn; i++){
		alphamzb = *alpha - zb[i];
		if (ind[i]){
		    s1 = z[j + i * (*mb)] * (  
			h_gamma_alpha(time[i], *gamma, alphamzb) / 
			h(time[i], *gamma, alphamzb) -
			h_alpha(time[i], *gamma, alphamzb) * 
			h_gamma(time[i], *gamma, alphamzb) / 
			R_pow_di(h(time[i], *gamma, alphamzb), 2));
		    tmp += s1;
		}
		s2 = S_gamma_alpha(time0[i], *gamma, alphamzb) / 
		    S(time0[i], *gamma, alphamzb, log_p) -
		    S_alpha(time0[i], *gamma, alphamzb) * 
		    S_gamma(time0[i], *gamma, alphamzb) / 
		    R_pow_di(S(time0[i], *gamma, alphamzb, log_p), 2);
		
		s3 = S_gamma_alpha(time[i], *gamma, alphamzb) / 
		    S(time[i], *gamma, alphamzb, log_p) -
		    S_alpha(time[i], *gamma, alphamzb) * 
		    S_gamma(time[i], *gamma, alphamzb) / 
		    R_pow_di(S(time[i], *gamma, alphamzb, log_p), 2);
		tmp -= z[j + i * (*mb)] * 
		    (s2 - s3); /* NOTE: -= */
	    }
	    fpp[((*mb) + 1) + bdim * j] = tmp;
	    fpp[j + bdim * ((*mb) + 1)] = tmp; /* Symmetry! */
	}
    

/* beta and alpha (log(scale)): */
    
	for (j = 0; j < (*mb); j++){
	    tmp = 0.0;
	    for (i = 0; i < *nn; i++){
		alphamzb = *alpha - zb[i];
		if (ind[i]){
		    s1 = z[j + i * (*mb)] * (  
			h_alpha2(time[i], *gamma, alphamzb) / 
			h(time[i], *gamma, alphamzb) -
			R_pow_di(h_alpha(time[i], *gamma, alphamzb) / 
				 h(time[i], *gamma, alphamzb), 2));
		    tmp += s1;
		}
		s2 = S_alpha2(time0[i], *gamma, alphamzb) / 
		    S(time0[i], *gamma, alphamzb, log_p) -
		    R_pow_di(S_gamma(time0[i], *gamma, alphamzb) / 
			     S(time0[i], *gamma, alphamzb, log_p), 2);
		
		s3 = S_alpha2(time[i], *gamma, alphamzb) / 
		    S(time[i], *gamma, alphamzb, log_p) -
		    R_pow_di(S_alpha(time[i], *gamma, alphamzb) / 
			     S(time[i], *gamma, alphamzb, log_p), 2);
		tmp -= z[j + i * (*mb)] * 
		    (s2 - s3); /* NOTE: -= */
	    }
	    fpp[(*mb) + bdim * j] = tmp;
	    fpp[j + bdim * (*mb)] = tmp;
	}
    } /* END if ((*mb) >= 1) */

/* Shape -- shape; hess[k+1][k+1]: */
    tmp = 0.0;
    for (i = 0; i < *nn; i++){
	alphamzb = *alpha - zb[i];
	if (ind[i]){
	    s1 =  h_gamma2(time[i], *gamma, alphamzb) / 
		h(time[i], *gamma, alphamzb) -
		R_pow_di(h_gamma(time[i], *gamma, alphamzb) / 
			 h(time[i], *gamma, alphamzb), 2);
	    tmp -= s1;
	}
	s2 = S_gamma2(time0[i], *gamma, alphamzb) / 
	    S(time0[i], *gamma, alphamzb, log_p) -
	    R_pow_di(S_gamma(time0[i], *gamma, alphamzb) / 
		     S(time0[i], *gamma, alphamzb, log_p), 2);
	
	s3 = S_gamma2(time[i], *gamma, alphamzb) / 
	    S(time[i], *gamma, alphamzb, log_p) -
	    R_pow_di(S_gamma(time[i], *gamma, alphamzb) / 
		     S(time[i], *gamma, alphamzb, log_p), 2);
	/* Rprintf("[d2_loglik_aft]: s3 = %f\n", s3); */
	tmp += (s2 - s3);
    }
    fpp[(*mb) + 1 + ((*mb) + 1) * bdim] = tmp;
    
/* Shape -- scale hess[k+1][k], hess[k][k+1] */
    tmp = 0.0;
    for (i = 0; i < *nn; i++){
	alphamzb = *alpha - zb[i];
	if (ind[i]){
	    s1 = h_gamma_alpha(time[i], *gamma, alphamzb) / 
		h(time[i], *gamma, alphamzb) -
		h_alpha(time[i], *gamma, alphamzb) * 
		h_gamma(time[i], *gamma, alphamzb) /
		R_pow_di(h(time[i], *gamma, alphamzb), 2);
	    tmp -= s1;
	}
	s2 = S_gamma_alpha(time0[i], *gamma, alphamzb) / 
	    S(time0[i], *gamma, alphamzb, log_p) -
	    S_alpha(time0[i], *gamma, alphamzb) * 
	    S_gamma(time0[i], *gamma, alphamzb) /
	    R_pow_di(S(time0[i], *gamma, alphamzb, log_p), 2);
	s3 = S_gamma_alpha(time[i], *gamma, alphamzb) / 
	    S(time[i], *gamma, alphamzb, log_p) -
	    S_alpha(time[i], *gamma, alphamzb) * 
	    S_gamma(time[i], *gamma, alphamzb) /
	    R_pow_di(S(time[i], *gamma, alphamzb, log_p), 2);
	tmp += (s2 - s3);
    }
    fpp[(*mb) + 1 + bdim * (*mb)] = tmp;
    fpp[(*mb) + bdim * ((*mb)+1)] = tmp;

/* Scale -- scale; hess[k][k]: */

    tmp = 0.0;
    for (i = 0; i < *nn; i++){
	alphamzb = *alpha - zb[i];
	if (ind[i]){
	    s1 =  h_alpha2(time[i], *gamma, alphamzb) / 
		h(time[i], *gamma, alphamzb) -
		R_pow_di(h_alpha(time[i], *gamma, alphamzb) / 
			 h(time[i], *gamma, alphamzb), 2);
	    tmp -= s1;
	}
	s2 = S_alpha2(time0[i], *gamma, alphamzb) / 
	    S(time0[i], *gamma, alphamzb, log_p) -
	    R_pow_di(S_alpha(time0[i], *gamma, alphamzb) / 
		     S(time0[i], *gamma, alphamzb, log_p), 2);
	
	s3 = S_alpha2(time[i], *gamma, alphamzb) / 
	    S(time[i], *gamma, alphamzb, log_p) -
	    R_pow_di(S_alpha(time[i], *gamma, alphamzb) / 
		     S(time[i], *gamma, alphamzb, log_p), 2);
	tmp += (s2 - s3); /* NOTE: -= */
    }
    fpp[(*mb) + bdim * (*mb)] = tmp;
    
    Free(zb);
}
