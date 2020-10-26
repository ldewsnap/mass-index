      implicit none
!gfortran gr.triangle2.f /usr/lib/libpgplot.so -o gr -lpgplot -O3 -ffast-math -march=native -ftree-vectorize!
      integer MAX,nx,ny,nz
      parameter(MAX=100000)
      parameter(nx=50,ny=50,nz=nx*ny)
      
      integer i,j,icol,irow,n,nparam,ix,iy,iz
      real*8 param(MAX,20),logz(MAX)
      real*8 cumul(MAX)
      real*4 xmin,xmax,ymin,ymax
      real*4 xx(MAX),yy(MAX),ref(20)
      real*4 tr(6),dens(nx,ny),z(nz),dx,dy,zmax,histo(nx)
      real*8 errplus,errminus,median

      nparam = 4
      
c     Set up zero point
      do i=1,nparam
         ref(i)=0.0
      end do
      ref(3)=0.0

c     Read data
      i = 0
      open(1,file='post_equal_weights.dat',status='old')
 100  continue
      read(1,*,end=101,err=101)(param(i,j),j=1,nparam),logz(i)
      do j=1,nparam
         param(i,j)=param(i,j)-ref(j)
      end do
      i=i+1
      go to 100
 101  continue
      n = i-1
      close(1)

c... Graphic part
      call pgbeg(0,'/PS',nparam,nparam)
      call pgpap(0.,1.0)
      call pgsch(3.0)
      call pgslw(1)

      do irow=1,nparam
         do icol=1,irow            

      if(irow.eq.icol) then

c... Make a histogram
         xmin = 1.e30
         xmax = -1.e30
         do i=1,n
            if(param(i,icol).gt.xmax)xmax=param(i,icol)
            if(param(i,icol).lt.xmin)xmin=param(i,icol)
         end do

c    Calculate median and errors
         do i=1,n
            cumul(i) = param(i,icol)
         end do
         call sort(n,cumul)

         median = cumul(int(0.5*n))+ref(icol)
         errminus = cumul(0.5d0*n*(1.d0-erf(1.d0/sqrt(2.d0))))+ref(icol)
         errplus = cumul(0.5d0*n*(1.d0+erf(1.d0/sqrt(2.d0))))+ref(icol) 
         write(*,*)median,errminus-median,errplus-median
  
c... Calculate density
         dx = (xmax-xmin)/float(nx)
         do ix=1,nx
            histo(ix) = 0.0
         end do
         do i=1,n
            ix = (param(i,icol)-xmin)/dx + 1
            histo(ix) = histo(ix) + (exp(logz(i)))
         end do

c... Calculate max
         ymin = 0.0
         ymax = -1.e10
         do ix=1,nx  
            if(histo(ix).gt.ymax) ymax = histo(ix)
         end do
         ymax = ymax*1.1
         do ix=1,nx  
            histo(ix) = histo(ix)/ymax
         end do
         ymax = 1.0

         call pgpanl(icol,irow)
         call pgvport(0.15,0.95,0.15,0.95)
         call pgswin(xmin,xmax,ymin,ymax)      
         
c         call pgscr(10,0.1,0.1,0.1)
c         call pgsci(10)
         call pgsfs(2)
         do ix=1,nx  
            xx(1) = xmin + (ix-1)*dx
            yy(1) = 0.0
            xx(2) = xmin + ix*dx
            yy(2) = 0.0
            xx(3) = xmin + ix*dx
            yy(3) = histo(ix)
            xx(4) = xmin + (ix-1)*dx
            yy(4) = histo(ix)
            xx(5) = xmin + (ix-1)*dx
            yy(5) = 0.0
            call pgpoly(5,xx,yy)
         end do
c         call pgsci(1)
         call pgsfs(1)

         call pgbox('bcnst', 0., 0, 'bcnst', 0., 0)
      else

c     Plot image          
         xmin = 1.e30
         xmax = -1.e30
         ymin = 1.e30
         ymax = -1.e30
         do i=1,n
            if(param(i,icol).gt.xmax)xmax=param(i,icol)
            if(param(i,icol).lt.xmin)xmin=param(i,icol)
            if(param(i,irow).gt.ymax)ymax=param(i,irow)
            if(param(i,irow).lt.ymin)ymin=param(i,irow)
         end do
         
c...  Set up plot variables
      tr(2)=(xmax-xmin)/float(nx-1)
      tr(1)=xmin-tr(2)
      tr(3)=0.0
      tr(6)=(ymax-ymin)/float(ny-1)
      tr(4)=ymin-tr(6)
      tr(5)=0.0

c... Calculate density
      dx = (xmax-xmin)/float(nx)
      dy = (ymax-ymin)/float(ny)
      do ix=1,nx
         do iy=1,ny
           dens(ix,iy) = 0.0
        end do
      end do
      do i=1,n
         ix = (param(i,icol)-xmin)/dx + 1
         if(ix.ge.1.and.ix.le.nx) then
            iy = (param(i,irow)-ymin)/dy + 1
            if(iy.ge.1.and.iy.le.ny) then
               dens(ix,iy) = dens(ix,iy) + 1.0
            end if
         end if
      end do
         
c... Move density into Z variable
      do ix=1,nx
         do iy=1,ny
            z( (iy-1)*nx + ix ) = dens(ix,iy)
         end do
      end do

c... Normalize Z variable to (0.1)
      zmax = -1.e10
      do iz=1,nz
         if(z(iz).gt.zmax) zmax = z(iz)
      end do
      do iz=1,nz
         z(iz) = z(iz)/zmax
      end do

      call pgpanl(icol,irow)
      call pgvport(0.15,0.95,0.15,0.95)
      call pgswin(xmin,xmax,ymin,ymax)      
        
c      do i=1,n
c         xx(i) = param(i,icol)
c         yy(i) = param(i,irow)
c      end do
c      call pgpt(n,xx,yy,-1)

      call pggray(z,nx,ny,1,nx,1,ny,1.0,0.0,tr)  
c      call pgimag(z,nx2,ny2,1,nx2,1,ny2,-8.0,-4.0,tr)         
   
      call pgbox('bcnst', 0., 0, 'bcnst', 0., 0)

      end if

         end do
      end do
      
      call pgend

      end

      FUNCTION erf(x)
      REAL*8 erf,x
CU    USES gammp
      REAL*8 gammp
      if(x.lt.0.d0)then
        erf=-gammp(.5d0,x**2)
      else
        erf=gammp(.5d0,x**2)
      endif
      return
      END

      FUNCTION gammp(a,x)
      REAL*8 a,gammp,x
CU    USES gcf,gser
      REAL*8 gammcf,gamser,gln
      if(x.lt.0.d0.or.a.le.0.d0)pause 'bad arguments in gammp'
      if(x.lt.a+1.d0)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.-gammcf
      endif
      return
      END

      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL*8 a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.d-7,FPMIN=1.d-30)
CU    USES gammln
      INTEGER i
      REAL*8 an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.d0-a
      c=1.d0/FPMIN
      d=1.d0/b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.d0
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.d0).lt.EPS)goto 1
 11   continue
      pause 'a too large, ITMAX too small in gcf'
 1    gammcf=exp(-x+a*log(x)-gln)*h
      return
      END

      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      REAL*8 a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.d-7)
CU    USES gammln
      INTEGER n
      REAL*8 ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.d0)then
        if(x.lt.0.d0)pause 'x < 0 in gser'
        gamser=0.d0
        return
      endif
      ap=a
      sum=1.d0/a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.d0
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
 11   continue
      pause 'a too large, ITMAX too small in gser'
 1    gamser=sum*exp(-x+a*log(x)-gln)
      return
      END

      FUNCTION gammln(xx)
      REAL*8 gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
 11   continue
      gammln=tmp+log(stp*ser/x)
      return
      END

      SUBROUTINE sort(n,arr)
      INTEGER n,M,NSTACK
      REAL*8 arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL*8 a,temp
      jstack=0
      l=1
      ir=n
 1    if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
 11      continue
         i=0
 2       arr(i+1)=a
 12   continue
      if(jstack.eq.0)return
      ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
 3      continue
          i=i+1
        if(arr(i).lt.a)goto 3
 4      continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
 5      arr(l)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END



