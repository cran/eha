      subroutine sup(method, iter, eps, prl, 
     &     totevent, totrs, ns, 
     &     antrs, antevents, size,
     &     totsize, eventset, riskset, 
     &     nn, antcov, covar, offset,
     &     startbeta, beta,
     &     loglik, dloglik, d2loglik, sctest,
     &     score, sumdscore, sumd2score, 
     &     conver, f_conver, fail)

C +++ 
C     method   : 0 = breslow,
C                1 = efron.
C     iter     : On input = maxiter; on output = actual No. of iterations.
C     eps      : Convergence criterion; L2 < eps ===> convergence
C     prl      : Print level; 0 = nothing, 1 = more.
C     totevent : Total number of events.
C     totrs    : Total number of risk sets.
C     ns       : Number of strata.
C
C     antrs     : antrs(i) = No. of risk sets in stratum i, i = 1, ns.
C     antevents : number of events in each riskset.
C     size      : Size of each risk set.
C
C     totsize  : Sum of the risk set sizes.
C     eventset : pointers to events in risk sets (length totevent).
C     riskset  : pointers to members of risk sets (length totsize).
C
c     nn     : No. of spells.
C     antcov : No. of covariates.
C     covar  : matrix of covariates (nn x antcov).
C     offset : Vector of offsets for each spell (nn).
C     beta   : Vector of coefficients.
C
C     loglik   : return value.
C     dloglik  : return value.
C     d2loglik : return value.
C
C     score, sumdscore, sumd2score: 'Work areas', avoiding local
C                                   dynamic memory allocation.
C                                   Used by 'coxfun'.
C     conver : 1 if convergence, 
C              0 otherwise.
C     fail   : 1 if failure (i.e., linear dependency among covariates
C                                    or singular hessian).
C              0 if success
C +++ 

      implicit none

      integer method, iter, prl
      double precision eps
      integer totevent, totrs, ns, totsize, nn, antcov
      integer antrs(ns), antevents(totrs), size(totrs)
      integer eventset(totevent), riskset(totsize)
      double precision covar(nn, antcov)
      double precision offset(nn)

      double precision startbeta(antcov), beta(antcov)

      double precision loglik(2), dloglik(antcov)
      double precision d2loglik(antcov, antcov), sctest

      double precision score(nn)
      double precision sumdscore(antcov) 
      double precision sumd2score(antcov, antcov)
      integer f_conver
      integer conver, fail
C +++
C     Local:
C
      integer what
      double precision zero, one
      parameter (zero = 0.d0, one = 1.d0) 
      integer ione
      parameter (ione = 1)
      character*1 trans
      parameter (trans = 'N')

      double precision ll

      integer info, i, j
      double precision L2
      double precision db(antcov)

C +++ Called BLAS functions:
      double precision dnrm2, ddot

      integer itmax

      integer job
      double precision det(2)

C +++ For dpodi to calculate the inverse only ('det' is a dummy):

      job = 1
C
C +++

      itmax = iter

      what = 2

      call dcopy(antcov, startbeta, ione, beta, ione)

      call coxfun(what, method,
     &     totevent, totrs, ns, 
     &     antrs, antevents, size,
     &     totsize, eventset, riskset, 
     &     nn, antcov, covar, offset,
     &     beta,
     &     ll, dloglik, d2loglik,
     &     score, sumdscore, sumd2score)

C *** return value for loglik(1):
      loglik(1) = ll
C *** for convergence check purposes:
      loglik(2) = ll

      iter = 0
      conver = 0 
      fail = 0
      f_conver = 0

      do while ( (iter .lt. itmax) .and. (conver .eq. 0) )

         call dcopy(antcov, dloglik, ione, db, ione) 
         call dpofa(d2loglik, antcov, antcov, info)
         if (info .eq. 0) then
            call dposl(d2loglik, antcov, antcov, db)
         else
C            call intpr("fail in [dpofa]; info =", 23, info, 1) 
C            fail = 1
            fail = info
            return
         endif
         
C +++ Score test when iter .eq. 0:
         if (iter .eq. 0) 
     &        sctest = ddot(antcov, dloglik, ione, db, ione)
            
         L2 = dnrm2(antcov, db, ione)
         if (L2 .lt. eps) conver = 1
         if (abs(one - ll / loglik(2)) .lt. eps) f_conver = 1
         loglik(2) = ll
         if (prl .eq. 1) then
            call intpr(" ", 1, iter, 0)
            call intpr('*** Iteration ', 14, iter, 1)
            call dblepr('L2 = ', 5, L2, 1)
            call dblepr('loglik = ', 9, ll, 1)
         endif

C         if (conver .eq. 0) then
C +++ Update beta:
            call daxpy(antcov, one, db, ione, beta, ione)
            
            call coxfun(what, method,
     &           totevent, totrs, ns, 
     &           antrs, antevents, size,
     &           totsize, eventset, riskset, 
     &           nn, antcov, covar, offset,
     &           beta,
     &           ll, dloglik, d2loglik,
     &           score, sumdscore, sumd2score)


C +++ New iteration
            iter = iter + 1
C         endif

      end do

C +++ Done! The afterwork:

C     The inverse of the hessian:
C      call dcopy(antcov, dloglik, ione, db, ione) 
      call dpofa(d2loglik, antcov, antcov, info)
      if (info .eq. 0) then
         call dpodi(d2loglik, antcov, antcov, det, job)
C --- Fill in the 'lower half' of d2loglik:
         do i = 2, antcov
            do j = 1, i - 1
               d2loglik(i, j) = d2loglik(j, i)
            enddo
         enddo
      else
         fail = info
         return
      endif


      if (prl .eq. 1) then
         call intpr(" ", 1, iter, 0)
         call intpr('*** Iteration ', 14, iter, 1)
         if (conver .eq. 1) then
            call intpr('Convergence', 11, iter, 0)
         else
            call intpr('NOTE: No convergence!', 21, iter, 0)
         endif
         call dblepr('loglik = ', 9, ll, 1)
      endif

      loglik(2) = ll
         
      return
      end
