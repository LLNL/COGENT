c
c      real x(10), y(10), b(10), c(10), d(10)
c      real s, u, seval
c      integer i, n
c
c      n=10
c      do 1 i=1,n
c      x(i)=i
c      y(i) = x(i)**3
c 1    continue
c
c      call spline (n,x,y,b,c,d)
c
c      u = 2.5
c      s = seval(n,u,x,y,b,c,d)
c      print 2, u, s
c      u = 3.9
c      s = seval(n,u,x,y,b,c,d)
c      print 2, u, s
c 2    format (2f10.5)
c      stop
c      end

      subroutine spline (n,x,y,b,c,d)
      integer n, i
      real*8 x(n), y(n), b(n), c(n), d(n)
      real*8 nux(n+2), nuy(n+2), nub(n+2), nuc(n+2), nud(n+2)
      real*8 x1, y1, x2, y2, x3, y3, q
      do 10 i=1,n
      nux(i+1) = x(i)
      nuy(i+1) = y(i)
 10   continue

      x1=x(1)
      y1=y(1)
      x2=x(2)
      y2=y(2)
      x3=x(3)
      y3=y(3)
      q = 2.0*x1-x2
      nux(1)=q
      nuy(1)=y1+(q-x1)*((y1-y2)/(x1-x2)
     +         +(q-x2)*((y2-y3)/(x2-x3)-(y1-y2)/(x1-x2))/(x3-x1))
      
      x1=x(n-2)
      y1=y(n-2)
      x2=x(n-1)
      y2=y(n-1)
      x3=x(n)
      y3=y(n)
      q = 2.0*x3-x2
      nux(n+2)=q
      nuy(n+2)=y1+(q-x1)*((y1-y2)/(x1-x2)
     +         +(q-x2)*((y2-y3)/(x2-x3)-(y1-y2)/(x1-x2))/(x3-x1))
      call spline2(n+2,nux,nuy,nub,nuc,nud)

      do 20 i=1,n
      b(i) = nub(i+1)
      c(i) = nuc(i+1)
      d(i) = nud(i+1)
 20   continue

      end
      
      subroutine spline2 (n,x,y,b,c,d)
      integer n
c      subroutine spline(CHF_INT[n],
c     & CHF_R1D[x],CHF_R1D[y],
c     & CHF_R1D[b],CHF_R1D[c],CHF_R1D[d])
      real*8 x(n), y(n), b(n), c(n), d(n)
      
      integer nm1, ib, i
      real*8 t
      nm1= n-1
      if (n.lt.2) return
      if (n.lt.3) go to 50
c      format (i5)

      d(1) = x(2)-x(1)
      c(2) = (y(2)-y(1))/d(1)
      do 10 i=2,nm1
      d(i) = x(i+1)-x(i)
      b(i) = 2.*(d(i-1)+d(i))
      c(i+1) = (y(i+1)-y(i))/d(i)
      c(i) = c(i+1) - c(i)
 10   continue

      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.
      c(n) = 0.
      if (n.eq.3) go to 15
      c(1)=c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
      c(n)=c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
      c(1)=c(1)*d(1)**2/(x(4)-x(1))
      c(n)=-c(n)*d(n-1)**2/(x(n)-x(n-3))

 15   do 20 i=2,n
      t=d(i-1)/b(i-1)
      b(i)=b(i)-t*d(i-1)
      c(i)=c(i)-t*c(i-1)
 20   continue

      c(n)=c(n)/b(n)
      do 30 ib=1,nm1
      i = n-ib
      c(i)=(c(i)-d(i)*c(i+1))/b(i)
 30   continue

      b(n)=(y(n)-y(nm1))/d(nm1)-d(nm1)*(c(nm1)+2.*c(n))
      do 40 i=1,nm1
         b(i)=(y(i+1)-y(i))/d(i)-d(i)*(c(i+1)+2.*c(i))
         d(i)=(c(i+1)-c(i))/d(i)
         c(i)=3.*c(i)
 40   continue
      c(n)=3.*c(n)
      d(n)=d(n-1)
      return

 50   b(1)=(y(2)-y(1))/(x(2)-x(1))
      c(1)=0.
      d(1)=0.
      b(2)=b(1)
      c(2)=0.
      d(2)=0.
      return
      end

c       subroutine seval(CHF_INT[n], CHF_REAL[u],
c     & CHF_R1D[x],CHF_R1D[y],CHF_R1D[b],
c     & CHF_R1D[c],CHF_R1D[d],CHF_REAL[val],CHF_REAL[der])
      real*8 function seval(n,u,x,y,b,c,d)
      integer n
      real*8 u,x(n),y(n),b(n),c(n),d(n)

      integer i,j,k
      real*8 dx
      data i/1/
      if (i.ge.n)i=1
      if (u.lt.x(i)) go to 10
      if (u.le.x(i+1)) go to 30

 10   i=1
      j=n+1
 20   k=(i+j)/2
      if (u.lt.x(k))j=k
      if(u.ge.x(k)) i=k
      if (j.gt.i+1) go to 20

 30   dx=u-x(i)
      seval = y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
c      der = b(i)+2.*dx*c(i)+3.*d(i)*dx**2
      return
      end

      real*8 function derval(n,u,x,y,b,c,d)
      integer n
      real*8 u,x(n),y(n),b(n),c(n),d(n)

      integer i,j,k
      real*8 dx
      data i/1/
      if (i.ge.n)i=1
      if (u.lt.x(i)) go to 10
      if (u.le.x(i+1)) go to 30

 10   i=1
      j=n+1
 20   k=(i+j)/2
      if (u.lt.x(k))j=k
      if(u.ge.x(k)) i=k
      if (j.gt.i+1) go to 20

 30   dx=u-x(i)
c      seval = y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
      derval = b(i)+2.*dx*c(i)+3.*d(i)*dx**2
      return
      end

      real*8 function derval2(n,u,x,y,b,c,d)
      integer n
      real*8 u,x(n),y(n),b(n),c(n),d(n)

      integer i,j,k
      real*8 dx
      data i/1/
      if (i.ge.n)i=1
      if (u.lt.x(i)) go to 10
      if (u.le.x(i+1)) go to 30

 10   i=1
      j=n+1
 20   k=(i+j)/2
      if (u.lt.x(k))j=k
      if(u.ge.x(k)) i=k
      if (j.gt.i+1) go to 20

 30   dx=u-x(i)
c      seval = y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
c      derval = b(i)+2.*dx*c(i)+3.*d(i)*dx**2
      derval2 = 2.*c(i)+6.*d(i)*dx
      return
      end

      real*8 function seval2(n,u,x,y,b,c,d,dr1,dr2)
      integer n
      real*8 u,x(n),y(n),b(n),c(n),d(n),dr1,dr2

      integer i,j,k
      real*8 dx
      data i/1/
      if (i.ge.n)i=1
      if (u.lt.x(i)) go to 10
      if (u.le.x(i+1)) go to 30

 10   i=1
      j=n+1
 20   k=(i+j)/2
      if (u.lt.x(k))j=k
      if(u.ge.x(k)) i=k
      if (j.gt.i+1) go to 20

 30   dx=u-x(i)
      dr1 = b(i)+2.*dx*c(i)+3.*d(i)*dx**2
      dr2 = 2.*c(i)+6.*d(i)*dx
      seval2 = y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
c      der = b(i)+2.*dx*c(i)+3.*d(i)*dx**2
      return
      end
