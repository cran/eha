C ***
C
      subroutine update_null(ord1, ord2, wind, wtime, 
     &     woffset, pfixed, p, alfa, bdim, b, 
     &     s, sy, syy)


      implicit none

      logical ord1, ord2
      integer wind, bdim
      double precision b(bdim), woffset
      logical pfixed
      double precision p, alfa
      double precision wtime, s, sy, syy

C +++ Local:
      double precision zero
      parameter (zero = 0.d0)
      double precision y, epy
      double precision tmp1

      y = log(wtime)
      
      epy = exp(p * (y + alfa))

      if (wind .eq. 2) epy = -epy
              
      s = s + epy
              
      if (ord1) then
         if (.not. pfixed) then
            tmp1 = y * epy
            sy = sy + tmp1
         endif
         
         if (ord2) then
            if (.not. pfixed) then
               syy = syy + tmp1 * y
            endif
         endif
c     +++          order .ge. 2
      endif
c     +++       order .ge. 1
              
      return
      end

C ***
C
      subroutine getsums_null(ord1, ord2, bdim, b, alfa, p, pfixed,
     +     nn, time, time0, ind, offset,
     +     s, sy, syy)

      implicit none

      logical ord1, ord2, pfixed

      integer bdim, nn, ind(nn)
      double precision time(nn), time0(nn), offset(nn)
      double precision b(bdim), alfa, p, s, sy, syy

C +++ Local:
      integer i

      double precision wtime
      integer wind

      double precision zero
      parameter (zero = 0.d0)

      s = zero
      sy = zero
      syy = zero

C      do 10 j1 = 1, k
C	 sz(j1) = zero
C	 syz(j1) = zero
C   10 continue

C      do 20 j1 = 1, index
C	 szz(j1) = zero
C   20 continue

      do 100 i = 1, nn
C         call GetRec(i, wtime, th0,
C     $              wz, woffset, wcommun, wind, wstratum, wrank, ok)
         wtime = time0(i)
	 if  (wtime .gt. zero) then
            wind = 2
            call update_null(ord1, ord2, wind, wtime,
     &           offset(i), pfixed, p, alfa, bdim, b,
     &           s, sy, syy)
         endif
         wtime = time(i)
         wind = ind(i)
         call update_null(ord1, ord2, wind, wtime,
     &        offset(i), pfixed, p, alfa, bdim, b, 
     &        s, sy, syy)
            
  100 continue

      return
      end
