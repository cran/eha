C (C) Göran Broström (2003).

C ***
C
      subroutine get_prob(nn, size, riskset, score, prob)

      implicit none

      integer nn, size, riskset(size)
      double precision score(nn), prob(size)

      double precision totsize
      integer i, who

      totsize = 0.d0
      do i = 1, size
         who = riskset(i)
         prob(i) = score(who)
         totsize = totsize + prob(i)
      enddo
      do i = 1, size
         prob(i) = prob(i) / totsize
      enddo

C +++ Sort p in decreasing order, apply same permutation to riskset:
      call rvsort(prob, riskset, size)
      
      return
      end

C ***
C
      subroutine sample_events(antevents, size, riskset, eventset,
     &     prob)
C +++
C     This is an adaption and translation to Fortran of 
C     the  C  function ProbSampleNoReplace found in random.c 
C     in the  R  sources.
 
      implicit none

      integer antevents, size, riskset(size), eventset(antevents)
      double precision prob(size)

      integer perm(size), k, n, i, j
      double precision z, p(size), rt, mass, totalmass

      if (antevents .ge. size) then
         do i = 1, size
            eventset(i) = riskset(i)
         enddo
C +++ This should _never_ happen:
         if (antevents .gt. size) 
     &        call rexit("Error in [sample_events]. Mail a bug report!")
      endif

      n = size - 1

C *** Take copies to protect prob and riskset:
      do i = 1, size
        perm(i) = riskset(i)
        p(i) = prob(i)
      enddo

      totalmass = 1.d0

      do i = 1, antevents
         call ranf(z)
         rt = z * totalmass
         mass = 0.d0
         j = 0
         do while ( (mass .lt. rt) .and. (j .le. n) )
            j = j + 1
            mass = mass + p(j)
         enddo
         eventset(i) = perm(j)
         totalmass = totalmass - p(j)
         do k = j, n
            p(k) = p(k + 1)
            perm(k) = perm(k + 1)
         enddo
         n = n - 1
      enddo
   
      return
      end
C ***
C
      subroutine bootcox(who, boot, sample, sd_sample,
     &     method, iter, eps, prl, 
     &     totevent, totrs, ns, 
     &     antrs, antevents, size,
     &     totsize, eventset, riskset, 
     &     nn, antcov, covar, offset,
     &     startbeta, beta,
     &     loglik, dloglik, d2loglik, sctest,
     &     score, sumdscore, sumd2score, 
     &     conver, fail)

C +++ 
C     who      : 1 = coxreg, i.e., call 'sup'
C                2 = mlreg,  i.e., call 'mlsup'
C     boot     : No. of bootstrap samples.
C     sample   : The return: a matrix of bootstrapped parameter values.
C  var_sample  : Its estimated sd.
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
C     offset : Vector of offsets (nn).
C     beta   : Vector of coefficients (antcov).
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

      integer who, boot, method, iter, prl
      double precision eps
      integer totevent, totrs, ns, totsize, nn, antcov
      integer antrs(ns), antevents(totrs), size(totrs)
      integer eventset(totevent), riskset(totsize)
C     Note dimensions of 'sample'!
      double precision sample(antcov, boot), sd_sample(antcov, boot)
      double precision covar(nn, antcov)
      double precision offset(nn)

      double precision startbeta(antcov), beta(antcov)

      double precision loglik(2), dloglik(antcov)
      double precision d2loglik(antcov, antcov), sctest

      double precision score(nn)
      double precision sumdscore(antcov) 
      double precision sumd2score(antcov, antcov)
      integer conver, fail
C +++
C     Local:
C
      double precision zero, one
      parameter (zero = 0.d0, one = 1.d0) 
      integer ione
      parameter (ione = 1)
      character*1 trans
      parameter (trans = 'N')

      double precision hazard(totrs)
      integer nsk2
      integer i, j, rs, rep, eindx, sindx, rsindx

      double precision prob(totsize)

      integer maxit
      integer f_conver
      
      maxit = iter

C +++ Calculate score(i), i = 1, nn:
C     Only needed for calculation of selection probabilities.
      call dcopy(nn, offset, ione, score, ione)
      call dgemv(trans, nn, antcov, one, covar, nn, beta, ione, one,  
     &     score, ione)

      do i = 1, nn
         score(i) = exp(score(i))
      enddo

      nsk2 = 0
      if (who .eq. 2) then
         rsindx = 0
         do rs = 1, ns
            do j = 1, antrs(rs)
               rsindx = rsindx + 1
               if ( (antevents(rsindx) .ge. 2) .and.
     &              (antevents(rsindx) .lt. size(rsindx)) )
     &              nsk2 = nsk2 + 1
            enddo
         enddo
      endif

C +++ get the (sorted) sample probabilities.
C     sort 'riskset' accordingly. 
C     Note the consequences in calling function!!

      sindx = 1
      rsindx = 0
      do rs = 1, ns
         do j = 1, antrs(rs)
            rsindx = rsindx + 1
            call get_prob(nn, size(rsindx), riskset(sindx), 
     &           score, prob(sindx))
            sindx = sindx + size(rsindx)
         enddo
      enddo

      call randstart()


      do rep = 1, boot
         if (prl .eq. 1) then
            call intpr('rep = ', -1, rep, 1)
         endif
         sindx = 1
         eindx = 1
         rsindx = 0
C +++ Sample riskset, put in eventset: 
         do rs = 1, ns
            do j = 1, antrs(rs)
               rsindx = rsindx + 1
               call sample_events(antevents(rsindx), size(rsindx),
     &              riskset(sindx), eventset(eindx), prob(sindx))
               sindx = sindx + size(rsindx)
               eindx = eindx + antevents(rsindx)
            enddo
         enddo
C +++ Get bootstrap sample:
         iter = maxit
         if (who .eq. 1) then

            call sup(method, iter, eps, prl, 
     &           totevent, totrs, ns, 
     &           antrs, antevents, size,
     &           totsize, eventset, riskset, 
     &           nn, antcov, covar, offset,
     &           startbeta, beta,
     &           loglik, dloglik, d2loglik, sctest,
C     &           hazard, Note: should be changed in sup. See 'mlsup'! 
     &           score, sumdscore, sumd2score, 
     &           conver, f_conver, fail)
         else
            call mlsup(method, iter, eps, prl, 
     &           totevent, totrs, ns, 
     &           antrs, antevents, size,
     &           totsize, eventset, riskset, 
     &           nn, nsk2, antcov, covar, offset,
     &           startbeta, beta,
     &           loglik, dloglik, d2loglik, sctest,
     &           hazard,
     &           score, sumdscore, sumd2score, 
     &           conver, f_conver, fail)
         endif
         call dcopy(antcov, beta, ione, sample(1, rep), ione)
         do j = 1, antcov
            sd_sample(j, rep) = sqrt(d2loglik(j, j))
         enddo
      enddo 

      call randend()

      return
      end
