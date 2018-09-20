      subroutine compton(md, mradlopts, mthetmidpts, rad1, thet1,
     +     sb1rm2, numinterps, yintvals, xintvals, outarray, ierr,
     +     mdderiv, sb1rsup)
      
      integer md, mradlopts, mthetmidpts, numinterps, ierr, mdderiv
      real*8 rad1(mradlopts), thet1(mthetmidpts)
c      real*8 sb1rm(mthetmidpts,mradlopts)
      real*8 sb1rm2(mradlopts,mthetmidpts)
      real*8 yintvals(numinterps), xintvals(numinterps),
     +          outarray(numinterps)
c      real*8 sb1rb(1:mradlopts,1:mthetmidpts)
c      real*8 sb1rc(1:mradlopts,1:mthetmidpts)
c      real*8 sb1rd(1:mradlopts,1:mthetmidpts)
      real*8 sb1rsup(1:mradlopts,1:mthetmidpts,3)
      real*8 temp

      real*8 sbb1rb(1:mthetmidpts)
      real*8 sbb1rc(1:mthetmidpts)
      real*8 sbb1rd(1:mthetmidpts)

      real*8 sbd1rb(1:mthetmidpts)
      real*8 sbd1rc(1:mthetmidpts)
      real*8 sbd1rd(1:mthetmidpts)

      real*8 st1ryint(mthetmidpts)

      real*8 st1ryder(mthetmidpts)

      real*8 sb1ry(mradlopts)
      real*8 syvals(mradlopts), sxvals(mthetmidpts)

      real*8  dist, lodist, hidist
      integer lo, loindx, hi, hiindx, nthpts
      real*8  dr1, dr2, u
      integer nsplinepts /6/
      integer j, ix, iy, k

      real*8 seval2


      ierr = 0

      if (md.eq.1) then
      do 90 ix=1,mthetmidpts
      do 80 iy=1,mradlopts
      sb1ry(iy)=sb1rm2(iy,ix)
 80   continue
      call spline(mradlopts,rad1,sb1ry,sb1rsup(1,ix,1),sb1rsup(1,ix,2),
     +  sb1rsup(1,ix,3))
c      call spline(mradlopts,rad1,sb1ry,sb1rb(1,ix),sb1rc(1,ix),
c     +  sb1rd(1,ix))
 90   continue

      endif

c      if (numinterps.eq.1) then
      do 2000 k=1,numinterps

      lodist = 1000.d0
      hidist = 1000.d0
      loindx = 0
      hiindx = 0

      do 1018 ix=1,mthetmidpts
c      sxvals(ix) = xlo1+0.25d0*(ix-1)
      sxvals(ix) = thet1(ix)
      dist = abs(thet1(ix)-xintvals(k))
      if (dist .lt. lodist) then
        lodist = dist
        loindx = ix
      endif 
      dist = abs(thet1(ix)-xintvals(k))
      if (dist .lt. hidist) then
        hidist = dist
        hiindx = ix
      endif 
1018  continue  
 
      if (loindx-nsplinepts .lt. 1) then 
        lo = 1
      else
        lo = loindx - nsplinepts
      endif 
      if (hiindx+nsplinepts .gt. mthetmidpts) then 
        hi = mthetmidpts
      else
        hi = hiindx + nsplinepts
      endif 
 
      do 1020 iy=1,mradlopts
c      syvals(iy) = ylo1+(iy-1)
      syvals(iy) = rad1(iy)
c      print 319, iy, xvals(iy)
c 319  format (i5, f10.5)
1020  continue 

c      do 1022 ix=1,itdiv1
c      hitheta = zer
c      if (itdiv1 > 1) hitheta = (ix-one)*(othi1-otlo1)/(itdiv1-one)
c      xintvals1(ix) = otlo1+hitheta
c      print 319, ix, xintvals1(ix)
c 319  format (i5, f10.5)
c 1022 continue   

c      do 1025 iy=1,irdiv1
c      hirho = zer
c      if (irdiv1 > 1) hirho = (iy-one)*(orhi1-orlo1)/(irdiv1-one)
c      yintvals1(iy) = orlo1+hirho
c1025  continue   

      do 1090 ix=lo,hi
      do 1085 iy=1,mradlopts
      sb1ry(iy) = sb1rm2(iy,ix)
c      sb1zy(iy) = sb1zm(ix,iy)
1085  continue
c      do 1086 iy=1,1
      u = yintvals(k)
      st1ryint(ix)=seval2(mradlopts,u,syvals,sb1ry,
     + sb1rsup(1,ix,1), sb1rsup(1,ix,2), sb1rsup(1,ix,3),dr1,dr2)
c      st1ryint(ix)=seval2(mradlopts,u,syvals,sb1ry,
c     + sb1rb(1,ix), sb1rc(1,ix), sb1rd(1,ix),dr1,dr2)
       st1ryder(ix)=dr1
c      st1zyint(ix,iy)=seval2(mradlopts,u,syvals,sb1zy,
c     + sb1zb(1,ix), sb1zc(1,ix), sb1zd(1,ix),dr1,dr2)
c      st1zyder(ix,iy)=dr1
c1086  continue
1090  continue

c      ix2m1 = (ixpt2-ixpt1+2*4)*4
c       ix2m1 = mthetmidpts-1
      nthpts = hi-lo+1

c      do 1095 iy=1,irdiv1
c      do 1095 iy=1,numinterps
      call spline(nthpts,sxvals(lo),st1ryint(lo),sbb1rb(lo),
     +  sbb1rc(lo),sbb1rd(lo))
c      call spline(nthpts,sxvals(lo),st1zyint(lo,iy),sbb1zb(lo,iy),
c     +  sbb1zc(lo,iy),sbb1zd(lo,iy))
      call spline(nthpts,sxvals(lo),st1ryder(lo),sbd1rb(lo),
     +  sbd1rc(lo),sbd1rd(lo))
c      call spline(nthpts,sxvals(lo),st1zyder(lo,iy),sbd1zb(lo,iy),
c     +  sbd1zc(lo,iy),sbd1zd(lo,iy))

c      do 1092 ix=1,numinterps
c      ix = iy
      u = xintvals(k)
      
      if ((mdderiv.eq.0).or.(mdderiv.eq.2)) then
      temp=seval2(nthpts,u,sxvals(lo),st1ryint(lo),
     + sbb1rb(lo),sbb1rc(lo),sbb1rd(lo),dr1,dr2)
      if (mdderiv.eq.0) outarray(k) = temp
      if (mdderiv.eq.2) outarray(k) = dr1
      endif
c      if (mdderiv.eq.2) outarray(ix) = dr2
c       outarray(ix)=temp
c      sb1rxder(ix,iy)=dr1
c      sb1rxder2(ix,iy)=dr2
      if (mdderiv.eq.1) then 
        outarray(k) = seval2(nthpts,u,sxvals(lo),st1ryder(lo),
     +   sbd1rb(lo),sbd1rc(lo),sbd1rd(lo),dr1,dr2)
      endif
c      sb1ryxder(ix,iy)=dr1
c      sb1zyint(ix,iy)=seval2(nthpts,u,sxvals(lo),st1zyint(lo,iy),
c     + sbb1zb(lo,iy),sbb1zc(lo,iy),sbb1zd(lo,iy),dr1,dr2)
c      sb1zxder(ix,iy)=dr1
c      sb1zxder2(ix,iy)=dr2
c      sb1zyder(ix,iy)=seval2(nthpts,u,sxvals(lo),st1zyder(lo,iy),
c     + sbd1zb(lo,iy),sbd1zc(lo,iy),sbd1zd(lo,iy),dr1,dr2)
c      sb1zyxder(ix,iy)=dr1

c 1092 continue
c 1095 continue
 2000 continue
c      endif

      end
