      subroutine fill_boundary_flux(
     &     dir,
     &     lo0, lo1, hi0, hi1,
     &     dlo0, dlo1, dhi0, dhi1,
     &     flux, flo0, flo1, fhi0, fhi1, fnc
     &     )
      implicit none
      integer 
     &     lo0, hi0, lo1, hi1, dir,
     &     dlo0, dhi0, dlo1, dhi1,
     &     flo0, fhi0, flo1, fhi1, fnc
      double precision
     &     flux(flo0:fhi0,flo1:fhi1, fnc)

c     local variables
      integer i, j, n

      if (dir .eq. 0) then

         if (lo1 .eq. dlo1) then
c           low physical boundary
            do n = 1, fnc
               do i = lo0, hi0+1
                  flux(i,lo1-1,n) =
     &                 3.d0 * (flux(i,lo1  ,n)
     &                       - flux(i,lo1+1,n))
     &                       + flux(i,lo1+2,n)
               enddo
            enddo
         endif

         if (hi1 .eq. dhi1) then
c           high physical boundary
            do n = 1, fnc
               do i = lo0, hi0+1
                  flux(i,hi1+1,n) =
     &                 3.d0 * (flux(i,hi1  ,n)
     &                       - flux(i,hi1-1,n))
     &                       + flux(i,hi1-2,n) 
               enddo
            enddo
         endif

      else if (dir .eq. 1) then

         if (lo0 .eq. dlo0) then
c           low physical boundary
            do n = 1, fnc
               do j = lo1, hi1+1
                  flux(lo0-1,j,n) =
     &                 3.d0 * (flux(lo0  ,j,n)
     &                       - flux(lo0+1,j,n))
     &                       + flux(lo0+2,j,n)
               enddo
            enddo
         endif

         if (hi0 .eq. dhi0) then
c           high physical boundary
            do n = 1, fnc
               do j = lo1, hi1+1
                  flux(hi0+1,j,n) =
     &                 3.d0 * (flux(hi0  ,j,n)
     &                       - flux(hi0-1,j,n))
     &                       + flux(hi0-2,j,n)
            enddo
         enddo
      endif

      endif

      return
      end
