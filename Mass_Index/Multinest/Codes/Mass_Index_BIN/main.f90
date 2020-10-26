PROGRAM main

	use params
	use nestwrapper
      
	implicit none
	
	
        REAL(8) :: tref,dummy
        REAL(8) :: binx(5000),binw,biny(5000)

        INTEGER :: cntr, ncount,pos,j,nbin
        REAL(8) :: x, y,arr(9999999),x1,x2
      	!no parameters to wrap around
      	nest_pWrap = 0
        ! OBTAIN PRIORS
        !--------------------------------------------------------------------------
        open(1,FILE='priors.in',FORM='FORMATTED',STATUS='UNKNOWN')
        DO j=1,sdim
          read(1,*) Rmin(j),Rmax(j)
          Rflag(j)=2
          Rsol(j)=(Rmax(j)+Rmin(j))/2d0 ! the mid point
          Rdel(j)=(Rmax(j)-Rmin(j))/2d0 ! the interval size
         write(*,*) Rsol(j),Rdel(j),Rmin(j),Rmax(j),Rflag(j)
        END DO
        
        read(1,*) nest_root
        close(1)
        nest_root=trim(nest_root)//"/"
	write(*,*) nest_root
        !--------------------------------------------------------------------------


 cntr=1 ! SET COUNTER TO 1
 open(77,file="datafile",status="OLD") ! OPEN DATAFILE WITH RAW DATA
 open(78,file="plot_data") ! CREATE FILE FOR DATA PLOTTING
 open(79,file="binned_data") ! CREATE FILE WITH BINNED DATA
 
 do while (1.gt.0) ! INFINITE LOOP
   read(77,*,end=666) arr(cntr) ! READ RAW DATA INTO AN ARRAY
   if(arr(cntr).gt.0d0) cntr=cntr+1 ! ADD ONE TO THE COUNTER
 enddo
 
666 continue
! binx=0d0 
! biny=0d0 ! 
! binw=0.01d0 ! actual bin width
!   call hpsort (cntr-1,arr) ! sort the array ascending order in the amplitude
!    do ncount=1,cntr-1 ! do loop thoughout all points in the array
!      x=log10(arr(cntr-ncount)) ! log of the amplitude, we start from the highest one
!      y=log10(ncount*1d0) ! log of the counts, we start from 1 (log10 = 0)
!      pos=int(x/binw)+1 ! position of the bin
!      write(*,*) pos,x,arr(cntr-ncount),ncount,cntr
!      if(pos.gt.int(5d0/binw)) pos=int(5d0/binw) !overflow precaution when log(amplitude) > 5 (should never happen)
!      binx(pos)=binx(pos)+1 ! number of meteor in pos bin (for averaging)
!      biny(pos)=biny(pos)+y ! averaging logarithms
!      write(78,*) x,y ! write raw data output to the file
!    enddo  
!   
  
 binx=0d0 
biny=0d0 ! 
nbin=100
! binw=0.01d0 ! ACTUAL BIN WIDTH
  CALL HPSORT (cntr-1,arr) ! SORT THE ARRAY ASCENDING ORDER IN THE AMPLITUDE
  x1=log10(arr(1))
  x2=log10(arr(cntr-1))
  binw=abs(x2-x1)/(nbin*1d0)
!   write(*,*) x1,x2, "RANGE"
!   stop
  
   do ncount=1,cntr-1 ! DO LOOP THOUGHOUT ALL POINTS IN THE ARRAY
     x=log10(arr(cntr-ncount)) ! LOG OF THE AMPLITUDE, WE START FROM THE HIGHEST ONE
     y=log10(ncount*1d0) ! LOG OF THE COUNTS, WE START FROM 1 (LOG10 = 0)
     pos=int((x-x1)/binw)+1 ! POSITION OF THE BIN
!      write(*,*) pos,x,arr(cntr-ncount),ncount,cntr
     if(pos.gt.nbin) pos=nbin !OVERFLOW PRECAUTION WHEN LOG(AMPLITUDE) > 5 (SHOULD NEVER HAPPEN)
     binx(pos)=binx(pos)+1 ! NUMBER OF METEOR IN POS BIN (FOR AVERAGING)
     biny(pos)=biny(pos)+y ! AVERAGING LOGARITHMS
     write(78,*) x,y ! WRITE RAW DATA OUTPUT TO THE FILE
   enddo   
!   
  
 do ncount=1,5000 
   if(binx(ncount).gt.0d0) then  ! WRITE ONLY BINS WITH NON-ZERO VALUES FOR FITTING PURPOSES
   write(79,*) x1+ncount*binw-binw/2d0,biny(ncount)/binx(ncount),dsqrt(binx(ncount))
   endif
 enddo
 
 
 close(79)  ! CLOSE ALL FILES
 close(77)
 close(78)

	write(*,*) "CALLING NEST_SAMPLE",rflag
      	call nest_Sample
END


  ! HEAP SORT SUBROUTINE 
  
  SUBROUTINE HPSORT(N,RA)
  integer N,IR,J,L
  
  real*8 RA(N),RRA
  L=N/2+1
  IR=N
  !  (C) Copr. 1986-92 Numerical Recipes Software
  !The index L will be decremented from its initial value during the
  !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
  !will be decremented from its initial value down to 1 during the
  !"retirement-and-promotion" (heap selection) phase.
10 continue
  if(L > 1)then
    L=L-1
    RRA=RA(L)
  else
    RRA=RA(IR)
    RA(IR)=RA(1)
    IR=IR-1
    if(IR.eq.1)then
      RA(1)=RRA
      return
    end if
  end if
  I=L
  J=L+L
20 if(J.le.IR)then
  if(J < IR)then
    if(RA(J) < RA(J+1))  J=J+1
  end if
  if(RRA < RA(J))then
    RA(I)=RA(J)
    I=J; J=J+J
  else
    J=IR+1
  end if
  goto 20
  end if
  RA(I)=RRA
  goto 10
END SUBROUTINE HPSORT
