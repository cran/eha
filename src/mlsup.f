C ***
C
      subroutine gamma_iter(nn, antevents,
     &              size, eventset, riskset,
     &              score, gamma)

      implicit none

      integer nn, antevents, size
      integer eventset(antevents), riskset(size)
      double precision score(nn), gamma

C +++ Local variables:
      integer who, i
      double precision totscore, lscore(antevents), dg, d2g, upd
      logical conver
      integer iter, itmax
      parameter (itmax = 10)
      double precision eps, egam, tmp, ehil
      parameter (eps = 1.d-8)

      do i = 1, antevents
         who = eventset(i)
         lscore(i) = score(who)
      enddo

      totscore = 0.d0
      do i = 1, size
         who = riskset(i)
         totscore = totscore + score(who)
      enddo

      iter = 0
      conver = .false.
C +++ Do some Newton-Raphson iterations:
      do while ( (.not. conver) .and. (iter < itmax) ) 
         egam = exp(gamma)
         dg = -totscore
         d2g = 0.d0
         do i = 1, antevents
            ehil = exp(-egam * lscore(i))
            tmp = lscore(i) / (1.d0 - ehil)
            dg = dg + tmp
            d2g = d2g + (tmp**2) * egam * ehil
         enddo
         upd = dg / d2g
         conver = (abs(upd) .le. eps)
         gamma = gamma + upd
         iter = iter + 1
      enddo
      
      if (iter .ge. itmax) 
     &     call intpr('No convergence in [gamma_iter]', 30, iter, 1) 

      return
      end

C ***
C
      subroutine init_haz(ns, antrs, totrs, antevents, size, 
     &     hazard, nsk2, gamma, iterate, nn, score,
     &     totevent, eventset, totsize, riskset)

      implicit none

      integer ns, antrs(ns), totrs, nsk2, antevents(totrs), size(totrs)
      double precision hazard(totrs), gamma(nsk2)
      logical iterate
      integer nn, totevent, totsize
      double precision score(nn)
      integer eventset(totevent), riskset(totsize)
      
      integer rs, j, gindx, rsindx, eindx, sindx
      double precision lambda
C +++
C     Start values:
      eindx = 1
      sindx = 1
      gindx = 0
      rsindx = 0
      do rs = 1, ns
         do j = 1, antrs(rs)
            rsindx = rsindx + 1
            lambda = dble(antevents(rsindx)) / dble(size(rsindx))
            hazard(rsindx) = lambda
            if ( (antevents(rsindx) .lt. size(rsindx)) .and.
     &           (antevents(rsindx) .ge. 2) ) then
               gindx = gindx + 1
               gamma(gindx) = log(-log(1.d0 - lambda))
               if (iterate) call gamma_iter(nn, antevents(rsindx),
     &              size(rsindx), eventset(eindx), riskset(sindx),
     &              score, gamma(gindx))
            endif
            eindx = eindx + antevents(rsindx)
            sindx = sindx + size(rsindx)
         enddo
      enddo

      return
      end

C ***
C
      subroutine fill_haz(ns, antrs, totrs, antevents, size, 
     &     hazard, nsk2, gamma, nn, antcov, covar, score, 
     &     totevent, eventset, totsize, riskset)

      implicit none

      integer ns, antrs(ns), totrs, nsk2, antevents(totrs), size(totrs)
      double precision hazard(totrs), gamma(nsk2)
      integer nn, antcov
      double precision covar(nn, antcov), score(nn)
      integer totevent, totsize, eventset(totevent), riskset(totsize)

      integer rs, j, m, eindx, sindx, gindx, rsindx, who
      double precision totscore, p1, r1
C +++
C     Start values:
      eindx = 1
      sindx = 1
      gindx = 0
      rsindx = 0
      do rs = 1, ns
         do j = 1, antrs(rs)
            rsindx = rsindx + 1
            if ( (antevents(rsindx) .lt. size(rsindx)) ) then
               if (antevents(rsindx) .eq. 1) then
                  who = eventset(eindx)
                  r1 = score(who)
                  totscore = 0.d0
                  do m = 1, size(rsindx)
                     who = riskset(sindx + m - 1)
                     totscore = totscore + score(who)
                  enddo
                  p1 = r1 / totscore
                  hazard(rsindx) = 1.d0 - (1.d0 - p1)**(1.d0 / r1) 
               else
                  gindx = gindx + 1
                  hazard(rsindx) = 1.d0 - exp(-exp(gamma(gindx)))
               endif
            endif
            eindx = eindx + antevents(rsindx)
            sindx = sindx + size(rsindx)
         enddo
      enddo

      return
      end
C ***
C
      subroutine inv_hess(antcov, nsk2, h1, h2, h11, h21, h22, 
     &     f, fail)

      implicit none

      integer antcov, nsk2, fail
      double precision h1(nsk2), h2(antcov), h11(nsk2)
      double precision h21(antcov, nsk2), h22(antcov, antcov)
      double precision f(nsk2, antcov)

      integer i, j, job
      character*1 transa, transb
      double precision one, zero
      parameter (one = 1.d0)
      parameter (zero = 0.d0)
      double precision det(2)

      transa = 'N'
      transb = 'T'

      do i = 1, nsk2
         do j = 1, antcov
            f(i, j) = h21(j, i) / h11(i)
         enddo
      enddo
      
C *** Store 'J = h22 - h21%*%f' in h22:

      if (nsk2 .ge. 1)
     &     call dgemm(transa, transa, antcov, antcov, nsk2, -one, 
     &     h21, antcov, f, nsk2, one, h22, antcov)

C *** Invert J, AKA h22: 

      call dpofa(h22, antcov, antcov, fail)
      if (fail .eq. 0) then
         job = 01
         call dpodi(h22, antcov, antcov, det, job)
         do i = 2, antcov
            do j = 1, i - 1
               h22(i, j) = h22(j, i)
            enddo
         enddo
      else
         return
      endif

C +++ Store -J**(-1) %*% t(F), AKA -h22 %*% t(f) in h21:

      if (nsk2 .ge. 1)
     &     call dgemm(transa, transb, antcov, nsk2, antcov, -one, 
     &     h22, antcov, f, nsk2, zero, h21, antcov)

      return
      end
C ***
C
      subroutine next_step(nsk2, antcov, h1, h2, h11, h21, h22, f, 
     &     dg, db)

      implicit none

      integer nsk2, antcov
      double precision h1(nsk2), h2(antcov), h11(nsk2)
      double precision h21(antcov, nsk2), h22(antcov, antcov)
      double precision f(nsk2, antcov), dg(nsk2), db(antcov)

      integer i, j, m
      double precision tmp, zero
      parameter (zero = 0.d0)

C     First dg:
      do i = 1, nsk2
         dg(i) = h1(i) / h11(i)
         do j = 1, nsk2
            tmp = zero
            do m = 1, antcov
               tmp = tmp + f(i, m) * h21(m, j)
            enddo
            dg(i) = dg(i) - h1(j) * tmp
         enddo
         do j = 1, antcov
            dg(i) = dg(i) + h2(j) * h21(j, i)
         enddo
      enddo
C     Then db:
      do i = 1, antcov
         db(i) = zero
         do j = 1, nsk2
            db(i) = db(i) + h1(j) * h21(i, j)
         enddo
         do j = 1, antcov
            db(i) = db(i) + h2(j) * h22(i, j)
         enddo
      enddo

      return
      end

C************************************************************
C *** The 'MAIN' subroutine:
C
      subroutine mlsup(method, iter, eps, prl, 
     &     totevent, totrs, ns, 
     &     antrs, antevents, size,
     &     totsize, eventset, riskset, 
     &     nn, nsk2, antcov, covar, offset,
     &     startbeta, beta,
     &     loglik, h2, h22, sctest,
     &     hazard,
     &     score, sumdscore, sumd2score, 
     &     conver, f_conver, fail)


C +++ 
C     method   : 0 = ML,
C                1 = MPPL, hybrid.
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
C     offset : Vector of offsets (nn).

C     startbeta : Start values for beta (antcov). 
C     beta      : Vector of coefficients (antcov); return value.
C
C     loglik         : return value.
C     h2 (dloglik)   : return value.
C     h22 (d2loglik) : return value.
C     sctest         : score test statistic, return value.
C
C     hazard : estimated hazard atoms (totrs).
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
      integer totevent, totrs, ns, totsize, nn, antcov, nsk2
      integer antrs(ns), antevents(totrs), size(totrs)
      integer eventset(totevent), riskset(totsize)
      double precision covar(nn, antcov)
      double precision  offset(nn)

      double precision startbeta(antcov), beta(antcov)

      double precision loglik(2)
      double precision sctest
      double precision hazard(totrs)
      double precision score(nn)
      double precision sumdscore(antcov) 
      double precision sumd2score(antcov, antcov)
      integer conver, f_conver, fail
C ************************************************************
C +++
C     variables needed for ML:
      double precision ll, h1(nsk2), h2(antcov), h11(nsk2) 
      double precision h21(antcov, nsk2), h22(antcov, antcov)
      double precision gamma(nsk2)
      double precision dg(nsk2), db(antcov)

      double precision f(nsk2, antcov)
      integer what

      integer itmax, i

      logical iterate 

      character*1 transa, transb
      double precision one
      parameter (one = 1.d0)
      double precision zero
      parameter (zero = 0.d0)
      integer ione
      parameter (ione = 1)
      character*1 trans
      parameter (trans = 'N')

      double precision ddot, dnrm2, L2

      transa = 'N'
      transb = 'T'

C +++ Get initial values for gamma:

      if (dnrm2(antcov, startbeta, ione) .gt. eps) then
         iterate = .TRUE.
         call dcopy(nn, offset, ione, score, ione)
         call dgemv(trans, nn, antcov, one, covar, nn, startbeta, ione, 
     &     one, score, ione)
      
         do i = 1, nn
            score(i) = exp(score(i))
         enddo
      else
         iterate = .FALSE.
      endif

      call init_haz(ns, antrs, totrs, antevents, size, 
     &     hazard, nsk2, gamma, iterate, nn, score,
     &     totevent, eventset, totsize, riskset)
C --- done!

      itmax = iter

      what = 2

      call dcopy(antcov, startbeta, ione, beta, ione)

      call mlfun(what, method,
     &     totevent, totrs, nsk2, ns, 
     &     antrs, antevents, size,
     &     totsize, eventset, riskset, 
     &     nn, antcov, covar, offset,
     &     beta, gamma,
     &     ll, h1, h2, h11, h21, h22,
     &     score)

      loglik(1) = ll
      loglik(2) = ll

      iter = 0
      conver = 0 
      f_conver = 0
      fail = 0

      do while ( (iter .lt. itmax) .and. (conver .eq. 0) )
C         iter = iter + 1

         call inv_hess(antcov, nsk2, h1, h2, h11, h21, h22, 
     &        f, fail)

            
         if (fail .ne. 0) then
C            call intpr('Info from [inv_hess] = ', -1, fail, 1)
            return
         endif

C +++ Calculate 'next step':

         call next_step(nsk2, antcov, h1, h2, h11, h21, h22, f, 
     &        dg, db)
            
C +++ The score test statistic:
         if (iter .eq. 0) then
            sctest = ddot(antcov, db, ione, h2, ione)
         endif

         L2 = dnrm2(antcov, db, ione) + dnrm2(nsk2, dg, ione)
         if (L2 .le. eps) conver = 1
         if (abs(one - ll / loglik(2)) .le. eps) f_conver = 1 
         if (prl .eq. 1) then
            call intpr(" ", 1, iter, 0)
            call intpr('*** Iteration ', 14, iter, 1)
            call dblepr('L2 = ', 5, L2, 1)
            call dblepr('loglik = ', 9, ll, 1)
         endif

C         if (conver .eq. 0) then
C +++ Update gamma, beta:
     
         call daxpy(antcov, one, db, ione, beta, ione)
         call daxpy(nsk2, one, dg, ione, gamma, ione)

         call mlfun(what, method,
     &        totevent, totrs, nsk2, ns, 
     &        antrs, antevents, size,
     &        totsize, eventset, riskset, 
     &        nn, antcov, covar, offset,
     &        beta, gamma,
     &        ll, h1, h2, h11, h21, h22,
     &        score)

         iter = iter + 1
      enddo
C +++ Done iterating!
C
C +++ The 'afterwork':
         
      if (prl .eq. 1) then
         
         call intpr(" ", 1, iter, 0)
         call intpr('*** Iteration ', 14, iter, 1)
         if (conver .eq. 1) then
            call intpr('Convergence', -1, iter, 0)
         else
            call intpr('NOTE: No Convergence!', -1, iter, 0)
         endif
         call dblepr('loglik = ', 9, ll, 1)
      endif
      
      loglik(2) = ll

C +++ Get the variance(antcov, antcov) matrix in h22:
      call inv_hess(antcov, nsk2, h1, h2, h11, h21, h22, f, fail)

      if (fail .ne. 0) then
C         call intpr('Last Info from [inv_hess] = ', -1, fail, 1)
         return
      endif

C +++ Fill in 'hazard':
      call fill_haz(ns, antrs, totrs, antevents, size, 
     &     hazard, nsk2, gamma, nn, antcov, covar, score, 
     &     totevent, eventset, totsize, riskset)
      
      return
      end
