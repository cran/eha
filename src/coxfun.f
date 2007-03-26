C ***
C
      subroutine init_lik(n, loglik, dloglik, d2loglik)

      implicit none

      integer n
      double precision loglik, dloglik(n), d2loglik(n, n)

      loglik = 0.d0
      call dcopy(n, 0.d0, 0, dloglik, 1)
      call dcopy(n * n, 0.d0, 0, d2loglik, 1)

      return
      end

C ***
C
      subroutine obs_lik(what, ns, antrs, totrs, antevents, size,
     &     totevent, eventset, nn, score, 
     &     antcov, covar,
     &     loglik, dloglik)

C +++
C     This subroutine calculates the 'observed' part of the
C     log likelihood function, and its first order
C     partial derivatives.
C +++

      implicit none

      integer what, ns, totrs, totevent, nn, antcov
      integer antrs(ns), antevents(totrs), size(totrs)
      integer eventset(totevent)
      double precision score(nn), covar(nn, antcov)
      double precision loglik, dloglik(antcov)

      integer rsindx, indx, rs, j, i, s, who

      indx = 0
      rsindx = 0
      do rs = 1, ns
         do j = 1, antrs(rs)
            rsindx = rsindx + 1
C +++ The following is not correct for variance estimates!!
C     So we change; "atomic" risksets can be renoverd.
C            if (antevents(rsindx) < size(rsindx)) then
            if (size(rsindx) .ge. 2) then
               do i = 1, antevents(rsindx)
                  indx = indx + 1
                  who = eventset(indx)
                  loglik = loglik + score(who)
                  if (what .ge. 1) then
                     call daxpy(antcov, 1.d0, covar(who, 1), nn,
     &                    dloglik, 1)
                  endif
               enddo
            else
               indx = indx + antevents(rsindx)
            endif
C +++ Next risk set!
         enddo
      enddo
      
      return
      end

C ***
C
      subroutine exp_lik(what, ns, antrs, totrs, antevents, size, 
     &     totsize, riskset, nn, score, 
     &     antcov, covar,
     &     sumdscore, sumd2score,
     &     loglik, dloglik, d2loglik)

C +++
C     This subroutine calculates the 'expected' part of the
C     log likelihood function, and its first and second order
C     partial derivatives.
C ***

      implicit none

      integer what, ns, totrs, totsize, nn, antcov
      integer antrs(ns), antevents(totrs), size(totrs), riskset(totsize)
      double precision score(nn), covar(nn, antcov)

      double precision sumdscore(antcov), sumd2score(antcov, antcov)
      double precision loglik, dloglik(antcov)
      double precision d2loglik(antcov, antcov)

      integer rsindx, indx, rs, j, i, s, m, who
      double precision sumscore

      character*1 UPLO
      parameter (UPLO = 'L')

      rsindx = 0
      indx = 0
      do rs = 1, ns
         do j = 1, antrs(rs)
            rsindx = rsindx + 1
C     Reset to zero:
            sumscore = 0.d0
            if (what .ge. 1) call dcopy(antcov, 0.d0, 0, sumdscore, 1)
            if (what .eq. 2) call dcopy(antcov * antcov, 0.d0, 0, 
     &           sumd2score, 1)
C     Go thru riskset(rs, j):
C +++ Not Removed! We don't go thru 'single' risksets too
C            if (antevents(rsindx) < size(rsindx)) then
            if (size(rsindx) .ge. 2) then
               do i = 1, size(rsindx)
                  indx = indx +  1
                  who = riskset(indx)
                  sumscore = sumscore + score(who)
                  if (what .ge. 1) call daxpy(antcov, score(who), 
     &                 covar(who, 1), nn, sumdscore, 1)
                  if (what .eq. 2) then
                     call dsyr(UPLO, antcov, score(who), 
     &                    covar(who, 1), nn, sumd2score, antcov)
C                     call dger(antcov, antcov, score(who), 
C     &                    covar(who, 1), nn, covar(who, 1), nn, 
C     &                    sumd2score, antcov)
                  endif
               enddo
C     Add into loglik:
               loglik = loglik - antevents(rsindx) * log(sumscore)
               if (what .ge. 1) then
C     Add into dloglik:
                  call daxpy(antcov, -antevents(rsindx) / sumscore,
     &                 sumdscore, 1, dloglik, 1)
               endif
               if (what .eq. 2) then
                  call daxpy(antcov * antcov, 
     &                 antevents(rsindx) / sumscore, sumd2score, 1,
     &                 d2loglik, 1)
                  call dsyr(UPLO, antcov, 
     &                 -antevents(rsindx) / sumscore**2, sumdscore, 1,
     &                 d2loglik, antcov)
 
C                  call dger(antcov, antcov, 
C     &                 -antevents(rsindx) / sumscore**2, sumdscore, 1,
C     &                 sumdscore, 1, d2loglik, antcov)
               endif
            else
               indx = indx + size(rsindx)
            endif
C     Done adding in; next riskset.
         enddo
      enddo
      
      return
      end

C     ***
C     
      subroutine efron_lik(what, ns, antrs, totrs, antevents, size, 
     &     totevent, totsize, eventset, riskset, nn, score, 
     &     antcov, covar,
     &     sumdscore, sumd2score,
     &     loglik, dloglik, d2loglik)

C     +++
C     This subroutine calculates the 'expected' part of the
C     log likelihood function, and its first and second order
C     partial derivatives.
C     ***

      implicit none

      integer what, ns, totrs, totevent, totsize, nn, antcov
      integer antrs(ns), antevents(totrs), size(totrs)
      integer eventset(totevent), riskset(totsize)
      double precision score(nn), covar(nn, antcov)

      double precision sumdscore(antcov), sumd2score(antcov, antcov)
      double precision loglik, dloglik(antcov)
      double precision d2loglik(antcov, antcov)

      integer rsindx, indx, eindx, rs, j, i, r, s, m, who
      double precision sumscore

C     +++
C     Local (note the deviation from strict standard here!):

      double precision escore, edscore(antcov), ed2score(antcov, antcov)
      double precision w, ws

      character*1 UPLO
      parameter (UPLO = 'L')

      double precision temp(antcov)

      rsindx = 0
      indx = 0
      eindx = 0

      do 1000, rs = 1, ns
         do 900, j = 1, antrs(rs)
            rsindx = rsindx + 1
C     Reset to zero:
            sumscore = 0.d0
            escore = 0.d0
            if (what .ge. 1) then
               call dcopy(antcov, 0.d0, 0, sumdscore, 1)
               call dcopy(antcov, 0.d0, 0, edscore, 1)
               if (what .eq. 2) then
                  call dcopy(antcov * antcov, 0.d0, 0, sumd2score, 1)
                  call dcopy(antcov * antcov, 0.d0, 0, ed2score, 1)
               endif
            endif
C     Go thru eventset(rs, j):
C            if (antevents(rsindx) < size(rsindx)) then
            if (size(rsindx) .ge. 2) then
               do 600, i = 1, antevents(rsindx)
                  eindx = eindx +  1
                  who = eventset(eindx)
                  escore = escore + score(who)
                  if (what .ge. 1) then
                     call daxpy(antcov, score(who), covar(who, 1), nn,
     &                    edscore, 1)
                     if (what .eq. 2) then
C                     call dsyr(UPLO, antcov, score(who), 
C     &                    covar(who, 1), 1, ed2score, antcov)

                        call dger(antcov, antcov, score(who), 
     &                       covar(who, 1), nn, covar(who, 1), nn,
     &                       ed2score, antcov)
                     endif
                  endif
 600           continue
C     Go thru riskset(rs, j):
               do 700, i = 1, size(rsindx)
                  indx = indx +  1
C     +++ Removed! We go thru 'single' risksets too
C     if (antevents(rsindx) < size(rsindx)) then
C     Added: skip "atomic" risksets:
                  if (size(rsindx) .ge. 2) then
                     who = riskset(indx)
                     sumscore = sumscore + score(who)
                     if (what .ge. 1) then
                        call daxpy(antcov, score(who), covar(who, 1),
     &                       nn, sumdscore, 1)
                        if (what .eq. 2) then
C                           call dsyr(UPLO, antcov, score(who),
C     &                          covar(who, 1), nn, sumd2score, antcov)
                           call dger(antcov, antcov, score(who),
     &                          covar(who, 1), nn, covar(who, 1), nn,
     &                          sumd2score, antcov)
                        endif
                     endif
                  endif
 700           continue
               
               if (antevents(rsindx) .eq. 1) then
C     Add into loglik:
                  loglik = loglik - log(sumscore)
                  if (what .ge. 1) then
C     Add into dloglik:
                     call daxpy(antcov, -1.d0 / sumscore, 
     &                    sumdscore, 1, dloglik, 1)
                     if (what .eq. 2) then
C     Add into d2loglik:
                        call daxpy(antcov * antcov, 
     &                       antevents(rsindx) / sumscore, 
     &                       sumd2score, 1, d2loglik, 1)
C                        call dsyr(UPLO, antcov, 
C     &                       -antevents(rsindx) / sumscore**2,
C     &                       sumdscore, 1, d2loglik, antcov)
                        call dger(antcov, antcov, 
     &                       -antevents(rsindx) / sumscore**2, 
     &                       sumdscore, 1, sumdscore, 1, 
     &                       d2loglik, antcov)
                        endif
                  endif
C     Done adding in; next riskset.
               else

C     +++ IF TIES:

                  do 800, r = 1, antevents(rsindx)
                     w = dble(r - 1) / dble(antevents(rsindx))
                     ws = w * escore
C     Add into loglik:
                     loglik = loglik - log(sumscore - ws)
                     if (what .ge. 1) then
C     Add into dloglik:
                        call dcopy(antcov, sumdscore, 1, temp, 1)
                        call daxpy(antcov, -w, edscore, 1, temp, 1)
                        call dscal(antcov, 1.d0 / (sumscore -ws), 
     &                       temp, 1)
                        call daxpy(antcov, -1.d0,
     &                       temp, 1, dloglik, 1)
                        if (what .eq. 2) then
C     Add into d2loglik:
C +++ New, next line / Remove again!!:
C                           call dscal(antcov, 1.d0 / (sumscore -ws), 
C     &                          temp, 1)
                           call daxpy(antcov * antcov, 
     &                          1.d0/(sumscore - ws), sumd2score, 1,
     &                          d2loglik, 1)
                           call daxpy(antcov * antcov, 
     &                          -w/(sumscore - ws), ed2score, 1,
     &                          d2loglik, 1)
                           call dsyr(UPLO, antcov, -1.d0, temp, 1,
     &                          d2loglik, antcov)
                        endif
                     endif
 800              continue
C     Done adding in; next riskset.
               endif
            else
               eindx = eindx + antevents(rsindx)
               indx = indx + size(rsindx)
            endif
 900     continue
 1000 continue

      return
      end

C ***
C
      subroutine coxfun(what, method,
     &     totevent, totrs, ns, 
     &     antrs, antevents, size,
     &     totsize, eventset, riskset, 
     &     nn, antcov, covar, offset,
     &     beta,
     &     loglik, dloglik, d2loglik,
     &     score, sumdscore, sumd2score)

C +++ 
C     what     : 0 = Only loglihood.
C                1 = Loglihood and first derivatives.
C                2 = Loglihood, first derivatives, and the negative hessian.
C            other = Nothing.
C     method   : 0 = breslow,
C                1 = efron.
C     totevent : Total number of events.
C     totrs    : Total number of risk sets.
C     ns       : Number of strata.
C
C     antrs     : antrs(i) = No. of risk sets in stratum i, i = 1, ns.
C     antevents : number of events in each riskset.
C     size      : Size of each risk set.
C
C     totsize  : Sum of the risk set sizes.
C     eventset : pointers to events of risk sets (length totevents).
C     riskset  : pointers to members of risk sets (length totsize).
C
c     nn     : No. of spells.
C     antcov : No. of covariates.
C     covar  : matrix of covariates (nn x antcov).
C     offset : Vector of offsets (nn).
C     beta   : Vector of coefficients (antciv).
C
C     loglik   : return value.
C     dloglik  : return value.
C     d2loglik : return value.
C
C     score, sumdscore, sumd2score: 'Work areas', avoiding local
C                                   dynamic memory allocation.
C +++

      implicit none

      integer what, method
      integer totevent, totrs, ns, totsize, nn, antcov
      integer antrs(ns), antevents(totrs), size(totrs)
      integer eventset(totevent), riskset(totsize)
      double precision covar(nn, antcov)
      double precision  offset(nn)

      double precision beta(antcov)

      double precision loglik, dloglik(antcov)
      double precision d2loglik(antcov, antcov)

C +++ Work areas:
      double precision score(nn)
      double precision sumdscore(antcov) 
      double precision sumd2score(antcov, antcov)
C ************************************************************
C     Local:
C
      double precision zero, one
      parameter (zero = 0.d0, one = 1.d0) 
      integer ione
      parameter (ione = 1)
      character*1 trans
      parameter (trans = 'N')

      integer i, s, m

C      integer indx, rsindx

C      double precision sumscore 
C
C *************************************************************

      call init_lik(antcov, loglik, dloglik, d2loglik)

      if ( (what .lt. 0) .or. (what .gt. 2) ) return

C +++ Calculate score(i), i = 1, nn:
      call dcopy(nn, offset, ione, score, ione)
      call dgemv(trans, nn, antcov, one, covar, nn, beta, ione, one,  
     &     score, ione)
      
C +++ 
C     First,add in the 'observed' parts of loglik and dloglik
C     No calculations for d2loglik here (is analytically zero).
C

      call obs_lik(what, ns, antrs, totrs, antevents, size, 
     &     totevent, eventset, nn, score, 
     &     antcov, covar,
     &     loglik, dloglik)

C +++
C     Then, exponentiate 'score':
C
      do i = 1, nn
         score(i) = exp(score(i))
      enddo

C +++
C     and add in the 'expected' parts of loglik, dloglik, and
C     d2loglik.
C
      if (method .eq. 0) then
         call exp_lik(what, ns, antrs, totrs, antevents, size, 
     &        totsize, riskset, nn, score, 
     &        antcov, covar,
     &        sumdscore, sumd2score,
     &        loglik, dloglik, d2loglik)
      else
         call efron_lik(what, ns, antrs, totrs, antevents, size, 
     &        totevent, totsize, eventset, riskset, nn, score, 
     &        antcov, covar,
     &        sumdscore, sumd2score,
     &        loglik, dloglik, d2loglik)
      endif
C +++
C     Fill in the upper triangle part of d2loglik (symmetry!):
C     already done above! (Maybe shouldnt...) Change Back:

      if (what .eq. 2) then

         do s = 1, antcov
            do m = s + 1, antcov
               d2loglik(s, m) = d2loglik(m, s)
            enddo
         enddo

      endif

      
      return
      end
