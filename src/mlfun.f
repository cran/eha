C +++
C
      subroutine init_ml(nsk2, antcov, loglik,
     &     h1, h2, h11, h21, h22)

      implicit none

      integer nsk2, antcov 
      double precision loglik, h1(nsk2), h2(antcov)
      double precision h11(nsk2)
      double precision h21(antcov, nsk2), h22(antcov, antcov)

      integer i, j
      double precision zero
      parameter (zero = 0.d0)

      loglik = zero

      do i = 1, nsk2
         h1(i) = zero
         h11(i) = zero
         do j = 1, antcov
            h21(j, i) = zero
         enddo
      enddo

      do j = 1, antcov
         h2(j) = zero
         do i = 1, antcov
            h22(j, i) = zero
         enddo
      enddo

      return
      end

      subroutine ml_rs(what, antevents, size,
     &                 eventset, riskset, nn, score, 
     &                 antcov, covar, gamma,
     &                 loglik, h1, h2, h11, h21, h22)

C +++ Calculates contributions from one risk set. Note that gamma
C     then is a scalar (a fixed element corresponding to the risk set).
C
C     what = 0: Only 'loglihood'
C            1: Loglihood and first derivatives.
C            2: Loglihood, first derivatives, and the negative hessian.

      implicit none

      integer what
      integer antevents, size
      integer eventset(antevents), riskset(size)
      integer nn, antcov
      double precision covar(nn, antcov), score(nn)
      double precision gamma
      double precision loglik, h1, h2(antcov), h11
      double precision h21(antcov), h22(antcov, antcov)
      
      integer j, m, i, who
      double precision egam, hil, ehil, bil, gil

      double precision one
      parameter (one = 1.d0)

      egam = exp(gamma)

C +++ Events:
      do i = 1, antevents
         who = eventset(i)
         hil = egam * score(who)
         ehil = exp(-hil)
         loglik = loglik + log(one - ehil) + hil
         if (what .ge. 1) then

C +++ First derivatives:
            bil = hil / (one - ehil)
            h1 = h1 + bil
            do j = 1, antcov
               h2(j) = h2(j) + covar(who, j) * bil
            enddo
            if (what .eq. 2) then

C +++ Second derivatives:
               gil = bil * (ehil + hil * ehil - one) / (one - ehil)
               h11 = h11 + gil
               do j = 1, antcov
                  h21(j) = h21(j) + covar(who, j) * gil
                  do m = 1, j
                     h22(j, m) = h22(j, m) + 
     &                    covar(who, j) * covar(who, m) * gil
                  enddo
               enddo
            endif
         endif
      enddo

C +++ All in riskset:
      do i = 1, size
         who = riskset(i)
         hil = egam * score(who)
         ehil = exp(-hil)
         loglik = loglik - hil
         if (what .ge. 1) then

C +++ First derivatives:
            h1 = h1 - hil
            do j = 1, antcov
               h2(j) = h2(j) - covar(who, j) * hil
            enddo
            if (what .eq. 2) then

C +++ Second derivatives:
               h11 = h11 + hil
               do j = 1, antcov
                  h21(j) = h21(j) + covar(who, j) * hil
                  do m = 1, j
                     h22(j, m) = h22(j, m) + 
     &                    covar(who, j) * covar(who, m) * hil
                  enddo
               enddo
            endif
         endif
      enddo

      return
      end
C ***
C
      subroutine prof_rs(what, antevents, size, 
     &     eventset, riskset, nn, score, 
     &     antcov, covar,
     &     loglik, h2, h22)

      implicit none

      integer what, antevents, size
      integer eventset, riskset(size)
      integer antcov, nn
      double precision score(nn), covar(nn, antcov)
      double precision loglik, h2(antcov), h22(antcov, antcov)

C +++ Local:
      double precision p1, q1, lqdp, totscore, ev 
      integer i, j, m, who
      double precision cbar(antcov), c1(antcov), u, ccov
      double precision one
      parameter (one = 1.d0)

C +++ Brutal, but a brutal error!
      if (antevents .ne. 1) stop 'Error in antevents in [profile]!'

      totscore = 0.d0

      do i = 1, size
         who = riskset(i)
         totscore = totscore + score(who)
      enddo
      p1 = score(eventset) / totscore
      q1 = one - p1
      lqdp = log(q1) / p1

      loglik = loglik + q1 * lqdp + log(p1)

      if (what .ge. 1) then
         do j = 1, antcov
            c1(j) = covar(eventset, j)
            ev = 0.d0
            do i = 1, size
               who = riskset(i)
               ev = ev + score(who) * covar(who, j)
            enddo
            cbar(j) = ev / totscore
            h2(j) = h2(j) - 
     &           lqdp * (c1(j) - cbar(j))
         enddo

         if (what .eq. 2) then
            u = one + q1 * lqdp

            do j = 1, antcov
               do m = 1, j
                  ccov = 0.d0
                  do i = 1, size
                     who = riskset(i)
                     ccov = ccov + covar(who, j) * covar(who, m) * 
     &                    score(who)
                  enddo
                  ccov = ccov / totscore
                  h22(j, m) = h22(j, m) -
     &                 (u / q1) * (c1(j) - cbar(j)) * (c1(m) - cbar(m))
     &                 - lqdp * (ccov - cbar(j) * cbar(m))
               enddo
            enddo
         endif

      endif

      return
      end
      
C ***
C
      subroutine cox_rs(what, ns, 
     &     antevents, size, 
     &     totevent, totsize, eventset, 
     &     riskset, nn, 
     &     score, 
     &     antcov, covar,
     &     loglik, dloglik, d2loglik)

      implicit none

      integer what, ns, antevents, size, totevent, totsize
      integer eventset(antevents), riskset(size)
      integer nn, antcov
      double precision score(nn), covar(nn, antcov)
      double precision loglik, dloglik(antcov), d2loglik(antcov, antcov)

      integer i, m, s, who

      double precision sumscore, sumdscore(antcov)
      double precision sumd2score(antcov, antcov)


C +++ "Observed" part:

C     Go thru eventset:

      do i = 1, antevents
         who = eventset(i)
         loglik = loglik + log(score(who))
         if (what .ge. 1) then
            do s = 1, antcov
               dloglik(s) = dloglik(s) + covar(who, s)
            enddo
         endif
      enddo

C +++ "Expected" part:

C     Reset to zero:
      sumscore = 0.d0
      if (what .ge. 1) then
         do s = 1, antcov
            sumdscore(s) = 0.d0
            if (what .eq. 2) then
               do m = 1, s
                  sumd2score(s, m) = 0.d0
               enddo
            endif
         enddo
      endif

C     Go thru riskset:

      do i = 1, size
         who = riskset(i)
         sumscore = sumscore + score(who)
         if (what .ge. 1) then
            do s = 1, antcov
               sumdscore(s) = sumdscore(s) + 
     &              covar(who, s) * score(who)
               if (what .eq. 2) then
                  do m = 1, s
C     Note the sign here! It is a minus in coxfun! (Correctly)
C     Should be unified.
                     sumd2score(s, m) = sumd2score(s, m) +
     &                    covar(who, s) * covar(who, m) * 
     &                    score(who)
                  enddo
               endif
            enddo 
         endif
         
      enddo
C     Add into loglik:
      loglik = loglik - antevents * log(sumscore)
      if (what .ge. 1) then
C     Add into dloglik:
         do s = 1, antcov
            dloglik(s) = dloglik(s) - antevents * 
     &           sumdscore(s) / sumscore
            if (what .eq. 2) then
C     Add into d2loglik:
               do m = 1, s
                  d2loglik(s, m) = d2loglik(s, m) + 
     &                 antevents * 
     &                 ( sumd2score(s, m) / sumscore -
     &                 sumdscore(s) * sumdscore(m) /
     &                 sumscore**2 )
               enddo
            endif
         enddo
      endif

      return
      end
C ***
C
      subroutine mlfun(what, method,
     &     totevent, totrs, nsk2, ns, 
     &     antrs, antevents, size,
     &     totsize, eventset, riskset, 
     &     nn, antcov, covar, offset,
     &     beta, gamma,
     &     loglik, h1, h2, h11, h21, h22,
     &     score)

C     +++ 
C     what     : 0 = Only log likelihood.
C                1 = log likelihood and first derivatives.
C                2 = Everything.
C            other = Nothing.
C     method   : 0 = ML (pure ML),
C                1 = MPPL, the 'hybrid' between ML and MPL.
C     totevent : Total number of events.
C     totrs    : Total number of risk sets.
C     nsk2     : Total number of risk sets with 2 <= entevents < size.
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
C     beta   : Vector of coefficients (antcov).
C     gamma  : Vector of hazard atoms (nsk2). 
C     
C     loglik,
C     h1, h2,
C     h11, h21, h22 : return values.
C     
C     score  : exp(covar %*% beta)
C     dynamic memory allocation.
C     +++

      implicit none

      integer what, method
      integer totevent, totrs, ns, nsk2, totsize, nn, antcov
      integer antrs(ns), antevents(totrs), size(totrs)
      integer eventset(totevent), riskset(totsize)
      double precision covar(nn, antcov), offset(nn)

      double precision beta(antcov), gamma(nsk2)

      double precision loglik

C     +++
C     variables needed for ML:

      double precision h1(nsk2), h2(antcov)
      double precision h11(nsk2) 
      double precision h21(antcov, nsk2), h22(antcov, antcov)

C     +++ Work areas:
      double precision score(nn)

C     ************************************************************
C     Local:
C     
      double precision zero, one
      parameter (zero = 0.d0, one = 1.d0) 
      integer ione
      parameter (ione = 1)
      character*1 trans
      parameter (trans = 'N')

      integer i, j, s, m, rs
      integer aindx, eindx, rsindx, sindx


C     *************************************************************

      call init_ml(nsk2, antcov, loglik, 
     &     h1, h2, h11, h21, h22)


      if ( (what .lt. 0) .or. (what .gt. 2) ) return

C     +++ Calculate score(i), i = 1, nn:
C
      call dcopy(nn, offset, ione, score, ione)
      call dgemv(trans, nn, antcov, one, covar, nn, beta, ione, one,  
     &     score, ione)
      
      do i = 1, nn
         score(i) = exp(score(i))
      enddo
C     
C     ---

      aindx = 0
      sindx = 1
      eindx = 1
      rsindx = 0
      do rs = 1, ns
         do j = 1, antrs(rs)
            rsindx = rsindx + 1
            if (antevents(rsindx) .lt. size(rsindx))then
               if (antevents(rsindx) .ge. 2) then
C     +++ Here ML applies; ties present:
                  aindx = aindx + 1
                  call ml_rs(what, antevents(rsindx), size(rsindx), 
     &                 eventset(eindx), riskset(sindx), nn, score, 
     &                 antcov, covar, gamma(aindx),
     &                 loglik, h1(aindx), h2, h11(aindx), 
     &                 h21(1, aindx), h22)
               else
C     +++ No ties; (0) profile out gamma (ML); or (1) MPPL (Cox):
                  if (method .eq. 0) then
C     +++ Ordinary ML:
                     call prof_rs(what, antevents(rsindx), size(rsindx), 
     &                    eventset(eindx), riskset(sindx), nn, score, 
     &                    antcov, covar,
     &                    loglik, h2, h22)
                  else

C     +++ The 'hybrid', MPPL:
                     call cox_rs(what, ns, 
     &                    antevents(rsindx), size(rsindx), 
     &                    totevent, totsize, eventset(eindx), 
     &                    riskset(sindx), nn, 
     &                    score, 
     &                    antcov, covar,
     &                    loglik, h2, h22)
                  endif
               endif
            endif
            sindx = sindx + size(rsindx)
            eindx = eindx + antevents(rsindx)
C     +++ Next risk set!
         enddo
      enddo

C     Fill in the upper triangle part of h22 (symmetry!):
C     

      if (what .eq. 2) then
         do s = 1, antcov
            do m = s + 1, antcov
               h22(s, m) = h22(m, s)
            enddo
         enddo

      endif
      
      return
      end
