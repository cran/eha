#include <R.h>
#include <R_ext/BLAS.h>
#include <Rmath.h>
#include "sup.h"
#include "eha_zeroin.h"
#include "coxfun.h"

extern RS_fun *eha_rs;

static double gam1_fun(double gam, void *info){

    RiskSet *risk;
    double dg, s, egam;
    int i, who;

    risk = info;

    dg = 0.0;
    for (i = 0; i < risk->size; i++){
	who = risk->riskset[i];
	s = score[who]; /* 'score' is global. */
	dg += s;
    }

    egam = exp(gam);
    for (i = 0; i < risk->antevents; i++){
	who = risk->riskset[i];
	s = score[who];
	dg += s / expm1(-s * egam);
    }


    return(-dg);
}

/* static double get1_gam(RiskSet *risk){ An error?! */
void get1_gam(RiskSet *risk){

    /* Binomial, cloglog link */
    int i, itmax;

    double gam, eps;

    double gmin, gmax;
    double ax, bx;

    int who;
    double what;


    if (risk->size == risk->antevents) {
	if (!(risk->out))
	    warning("[get1_gam] gamma positive infinite");
	risk->gamma = 1000.0;
	risk->hazard = 1.0;
	return;
    }else if (risk->size == 1){
	if (!(risk->out))
	    warning("[get1_gam] gamma negative infinite");
	risk->gamma = -1000;
	risk->hazard = 0.0;
	return;
    }

    itmax = 25;
    eps = 0.000001;

    risk->tot_score = 0.0;
    who = risk->riskset[0];
    risk->tot_score += score[who];
    gmin = lin[who];
    gmax = gmin;
    for (i = 1; i < risk->size; i++){
	who = risk->riskset[i];
	risk->tot_score += score[who];
	what = lin[who];
	if (what < gmin){ 
	    gmin = what;
	}else{
	    if (what > gmax) gmax = what;
	}
    }
    if (risk->antevents == 1){
	who = risk->eventset[0];
	/*
	  gam = 1.0 - 
	    R_pow_di(1.0 - score[who] / risk->tot_score, 1.0 / score[who]);
	*/
	gam = log(-log1p(-score[who] / risk->tot_score) / score[who]);
    }else{
	gam = 1.0 - (double)(risk->antevents) / risk->tot_score;
	gam = log(-log1p(-(double)(risk->antevents) / 
			 (double)(risk->size)) ); /* start value */
	ax = gam - gmax;
	bx = gam - gmin;
	if (abs (ax - bx) < eps){
	    risk->gamma = (ax + bx) / 2.0;
	}else{
	    if (gam1_fun(ax, risk) * gam1_fun(bx, risk) > 0.0){
		Rprintf("f(%f) = %f, f(%f) = %f\n", 
			ax, gam1_fun(ax, risk), bx, gam1_fun(bx, risk));
		Rprintf("antevents = %f\n", risk->antevents); 
		Rprintf("size = %f\n", risk->size); 
		error("\nWrong interval in [get0_gam]");
	    }
	    gam = eha_zeroin(ax, bx, &gam1_fun, risk, &eps, &itmax);
	}
    }
    risk->gamma = gam;
    risk->hazard = 1.0 - exp(-exp(gam));
}

void ml_rs(int what, RiskSet *risk,
	   double *b, double e_frac,
	   double *loglik, double *dloglik, double *d2loglik){

    double h1, h11;
    double *h21;

    int i, who;
    double egam, hil, ehil, bil, gil;

    double zero = 0.0;
    double one = 1.0;

    int izero = 0;
    int ione = 1;

    char up = 'U';

    if (risk->out) return;    if (risk->antevents == risk->size) return;

    /* get "gamma[rs]"; in risks[rs]: */
    get1_gam(risk);

/*  Initialize:  */
    h21 = sumdscore;
    h1 = zero;
    h11 = zero;
    F77_CALL(dcopy)(&p, &zero, &izero, h21, &ione);

    egam = exp(risk->gamma);

/* Events : */

    for (i = 0; i < risk->antevents; i++){
	who = risk->eventset[i];
	hil = egam * score[who];
	ehil = exp(-hil);
	/* loglik += log(one - ehil) + hil; */
	*loglik += log1p(-ehil) + hil;

	if (what >= 1){  /* First derivatives */
	    bil = hil / (one - ehil);
	    /* bil = -hil / expm1(-hil); */
            h1 += bil;
            /* Update h2[j] "+= x[j, who] * bil": */
	    F77_CALL(daxpy)(&p, &bil, (x + p * who), &ione, 
			    dloglik, &ione);
	    
	    if (what >= 2){ /* Second derivatives */
		gil = bil * (ehil + hil * ehil - one) / (one - ehil); 
		    /* gil = bil * (bil * ehil - one); */
		h11 += gil;
		/* Update h21: */
		F77_CALL(daxpy)(&p, &gil, (x + p * who), &ione, 
				h21, &ione);
		/* Update h22 (upper triangle): */
		F77_CALL(dsyr)(&up, &p, &gil, (x + p * who), &ione, 
			       d2loglik, &p);
	    }
	}
    }

/* All in riskset: */
    for (i = 0; i < risk->size; i++){
	who = risk->riskset[i];
	hil = -egam * score[who]; /* Note 1! */
	ehil = exp(hil);
	*loglik += hil;
	if (what >= 1){ /* First derivatives */
	    h1 += hil;
	    F77_CALL(daxpy)(&p, &hil, (x + p * who), &ione,
			    dloglik, &ione);

	    if (what >= 2){ /* Second derivatives: */
		hil = -hil; /* Note 2! */
		h11 += hil; /* Added 0.99-15 */
		F77_CALL(daxpy)(&p, &hil, (x + p * who), &ione, 
				h21, &ione);
		F77_CALL(dsyr)(&up, &p, &hil, (x + p * who), &ione,
			       d2loglik, &p);
	    }
	}
    }

    if (what >= 2){
	h11 = -one / h11;
	F77_CALL(dsyr)(&up, &p, &h11, h21, &ione, d2loglik, &p);
    }
}

static void cox_obs_rs(int what, RiskSet *risk,
		       double *b,
		       double *loglik, double *dloglik){

    int i, who;
    double one = 1.0;
    int ione = 1;

    if (risk->out) return;

    for (i = 0; i < risk->antevents; i++){
	who = risk->eventset[i];
	*loglik += lin[who];
	if (what >= 1) F77_CALL(daxpy)(&p, &one, (x + p * who), &ione,
				       dloglik, &ione);
    }
}
    
void breslow_rs(int what, RiskSet *risk, 
		double *b, double e_frac,
		double *loglik, double *dloglik, 
		double *d2loglik){

/* Calculates the "breslow" variation of the partial likelihood, */
/* and eventually its first and second partial derivatives.      */

    int i, who;
    char up = 'U';

    double zero = 0.0;
    int izero = 0;
    int ione = 1;
    double sumscore;
    int p2 = p * p;
    double alpha;
    
    if (risk->out) return;

    /* First the "observed" part (common to the 'efron' method): */
    
    cox_obs_rs(what, risk,
	       b,
	       loglik, dloglik);
    
    /* Then the "expected" part: */

    /* Initialize: */

    sumscore = 0.0;
    if (what >= 1){
	F77_CALL(dcopy)(&p, &zero, &izero, sumdscore, &ione);
	if (what >= 2){
	    F77_CALL(dcopy)(&p2, &zero, &izero, sumd2score, &ione);
	}
    }

    /* Go thru riskset: */

    for (i = 0; i < risk->size; i++){
	who = risk->riskset[i];
	sumscore += score[who];
	if (what >= 1){ /* First derivatives: */
	    F77_CALL(daxpy)(&p, (score + who), (x + p * who), &ione, 
			    sumdscore, &ione);

	    if (what >= 2){ /* Second derivatives: */
		F77_CALL(dsyr)(&up, &p, (score + who), 
			       (x + p * who), &ione,
			       sumd2score, &p);
	    }
	}
    }

    /* Add in: */

    *loglik -= risk->antevents * log(sumscore);
    if (what >= 1){
	alpha = -(double)(risk->antevents) / sumscore;
	F77_CALL(daxpy)(&p, &alpha, sumdscore, &ione, 
			dloglik, &ione);
	if (what >= 2){
	    alpha = -alpha;
	    F77_CALL(daxpy)(&p2, &alpha, sumd2score, &ione,
			    d2loglik, &ione);
	    alpha = -alpha / sumscore;
	    F77_CALL(dsyr)(&up, &p, &alpha, sumdscore, &ione,
			   d2loglik, &p);
	}
    }

}

void efron_rs(int what, RiskSet *risk, 
	      double *b, double e_frac,
			 double *loglik, double *dloglik, 
			 double *d2loglik){
  
    /*************************************************************
C     This subroutine calculates the 'efron' variation of the
C     log likelihood function, and its first and second order
C     partial derivatives.
    *************************************************************/

    int i, r, who;
    double sumscore;

/*************************************************************
C     +++
C     Local (note the deviation from strict standard here!):
**************************************************************/

    double escore;
    double *edscore, *ed2score;
    double w, ws;

    char up = 'U';

    double *temp; /*(antcov)*/
    int p2 = p * p;

    int izero = 0;
    int ione = 1;

    double zero = 0.0;
    double one = 1.0;

    double alpha;

    if (risk->out) return;

    /* First the "observed" part (common to the 'breslow' method): */

    cox_obs_rs(what, risk,
	       b,
	       loglik, dloglik);

    /* Then the "expected" part: */

    edscore = Calloc(p, double);
    ed2score = Calloc(p2, double);

    temp = Calloc(p, double);
		     

/*     Reset to zero: */
    sumscore = zero;
    escore = zero;
    if (what >= 1){
	F77_CALL(dcopy)(&p, &zero, &izero, sumdscore, &ione);
        F77_CALL(dcopy)(&p, &zero, &izero, edscore, &ione);
	if (what >= 2){
	    F77_CALL(dcopy)(&p2, &zero, &izero, sumd2score, &ione);
	    F77_CALL(dcopy)(&p2, &zero, &izero, ed2score, &ione);
	}
    }

/*     Go thru riskset(rs, j): */
    for (i = 0; i < risk->size; i++){
	who = risk->riskset[i];
	sumscore += score[who];
	if (what >= 1){
	    F77_CALL(daxpy)(&p, (score + who), (x + p * who),
			    &ione, sumdscore, &ione);
	    if (what >= 2){
		F77_CALL(dsyr)(&up, &p, (score + who),
                               (x + p * who), &ione, 
                               sumd2score, &p);
/*		F77_CALL(dger)(&p, &p, (score + who), (x + p * who), &ione,
		(x + p * who), &ione, sumd2score, &p); */
		
	    }
	}
    }
               
    if (risk->antevents == 1){ /* No ties */
/*     Add into loglik: */
	*loglik -= log(sumscore);
	if (what >= 1){
/*     Add into dloglik: */
	    alpha = -one / sumscore;
	    F77_CALL(daxpy)(&p, &alpha, 
			    sumdscore, &ione, dloglik, &ione);
	    if (what >= 2){
/*     Add into d2loglik: */
		alpha = (double)risk->antevents / sumscore; 
		F77_CALL(daxpy)(&p2, &alpha, sumd2score, &ione, 
				d2loglik, &ione);
		alpha = -(double)(risk->antevents) / 
		    (sumscore * sumscore);
		F77_CALL(dsyr)(&up, &p, &alpha, sumdscore, &ione, 
		d2loglik, &p); 
/*		F77_CALL(dger)(&p, &p, &alpha, sumdscore, &ione, 
			       sumdscore, &ione, d2loglik, &p); */
	    }
	}
    }else{
/*     +++ IF TIES: */
	
/* Go thru events and create escore, edscore and ed2score: */
	for (i = 0; i < risk->antevents; i++){
	    who = risk->eventset[i];
	    escore = escore + score[who];
	    if (what >= 1){ /* first derivatives */
		F77_CALL(daxpy)(&p, (score + who), (x + p * who), &ione,
				edscore, &ione);
		if (what >= 2){ /* second derivatives */
		    F77_CALL(dsyr)(&up, &p, (score + who), 
				   (x + p * who), &ione,
				   ed2score, &p); 
/*	       F77_CALL(dger)(&p, &p, (score + who), (x + p * who), &ione,
	       (x + p * who), &ione, ed2score, &p); */
		}
	    }
	}


	for (r = 0; r < risk->antevents; r++){
/* WRONG: (Fortran)   w = (double)(r - ione) / (double)(risk->antevents);*/
	    w = (double)r / (double)(risk->antevents);
	    ws = w * escore;
/*     Add into loglik: */
	    *loglik -= log(sumscore - ws);
	    if (what >= 1){
/*     Add into dloglik: */
		F77_CALL(dcopy)(&p, sumdscore, &ione, temp, &ione);
		alpha = -w;
		F77_CALL(daxpy)(&p, &alpha, edscore, &ione, temp, &ione);
		alpha = one / (sumscore - ws);
		F77_CALL(dscal)(&p, &alpha, temp, &ione);
		alpha = -one;
		F77_CALL(daxpy)(&p, &alpha, temp, &ione, dloglik, &ione);
		if (what >= 2){
/*     Add into d2loglik: */
		    alpha = one / (sumscore - ws);
		    F77_CALL(daxpy)(&p2, &alpha, sumd2score, &ione,
				    d2loglik, &ione);
		    alpha = -w / (sumscore - ws);
		    F77_CALL(daxpy)(&p2, &alpha, ed2score, &ione,
				    d2loglik, &ione);
		    alpha = -one;
		    F77_CALL(dsyr)(&up, &p, &alpha, temp, &ione,
				   d2loglik, &p);
		}
	    }
	}
    }
    Free(temp);
    Free(ed2score);
    Free(ed2score);
}


void mppl_rs(int what, RiskSet *risk,
	     double *b, double e_frac,
	     double *loglik, double *dloglik, double *d2loglik){

    if (risk->antevents == risk->size) return;
    if (risk->out) return;
    if (risk->antevents == 1){
	breslow_rs(what, risk,
		   b, e_frac,
		       loglik, dloglik, d2loglik);
    }else if (risk->antevents <= e_frac * risk->size){
	efron_rs(what, risk, 
		 b, e_frac,
		 loglik, dloglik, d2loglik);
    }else{
	ml_rs(what, risk,
	      b, e_frac,
	      loglik, dloglik, d2loglik);
    }
}

void coxfun(int what, int totrs, RiskSet *risks, 
/*	    int ml, int method, */ double e_frac, 
	    double *b, 
	    double *loglik, double *dloglik, double *d2loglik){

    /* Note that the following is deprecated:      */
    /* ml : 0 = "Cox", 1 = "Maximum likelihood"    */
    /* method: 0 = "efron",  if ml = 0,            */ 
    /*         1 = "breslow, if ml = 0             */
    /*                                             */
    /* method: 0 = "mppl",   if ml = 1,            */ 
    /*         1 = "ml,      if ml = 1             */
    /*                                             */


    int izero = 0;
    int ione = 1;
    double zero = 0.0;
    double one = 1.0;
    char trans = 'T';

    int p2 = p * p;
    int s, m, rs;

    /* Initialize: */

    *loglik = 0.0;
    if (what < 0) return;

    if (what >= 1){
	F77_CALL(dcopy)(&p, &zero, &izero, dloglik, &ione);
	if (what >= 2){
	    F77_CALL(dcopy)(&p2, &zero, &izero, d2loglik, &ione);  
	}
    }

    /* Calculate 'lin' and 'score' */

    F77_CALL(dcopy)(&nn, offset, &ione, lin, &ione);
    F77_CALL(dgemv)(&trans, &p, &nn, &one, x, 
		    &p, b, &ione, &one, lin, &ione); 
    for (m = 0; m < nn; m++) score[m] = exp(lin[m]);

    /* Start walking thru risksets: */
    
    for (rs = 0; rs < totrs; rs++){
	eha_rs(what, (risks + rs),
	       b, e_frac, 
	       loglik, dloglik, d2loglik);
    }
	/* Fill in the lower part of 'd2loglik': */

    if (what >= 2){
	for (s = 0; s < p; s++){
	    for (m = 0; m < s; m++){
		d2loglik[s + m * p] = d2loglik[m + s * p];
	    }
	}
    }
}

