MODULE transitmod

use params
! use planmod
! use jasminemod

implicit none
      
contains
      
      
!=======================================================================
  SUBROUTINE transit(Rin,ressy)

  implicit none  


! Input variables
  REAL(8), DIMENSION(nest_nPar) :: Rin           ! input parameters

! Output
  REAL(8) :: ch2,ks                                ! merit function of fit  
  REAL(8) :: ressy                ! C-O difference    

! Internal variables 

  REAL*8 a,b,x(5000),y(5000),sigma(5000), count,xb1,xb2
! Original variables
!   REAL(8), DIMENSION(taulen) :: tauvec

	character(80) flnm,dumy,timc,flnm2	
        character(20) dummy2,cmorfile
        integer ixb1,ixb2,cntr,ncount

!---------------------------------------------------------

!  read(*,*) mean
 
 a=rin(1) ! f(x) = a * x + b
 b=rin(2)
 xb1=rin(3) ! Lower Boundary
 xb2=rin(4) ! Upper Boundary
 
 ixb1=-9999999 ! index of the lower boundary
 ixb2=-9999999 ! index of the upper boundary
 cntr=0 ! index of the lower boundary

 if (xb1.ge.xb2) goto 888 ! tell code that he's wrong when boundaries are switched
 
 open(88,file="binned_data",status="OLD") ! open binned_data
 do while (1.gt.0)
    cntr=cntr+1
   read(88,*,end=666) x(cntr),y(cntr),sigma(cntr)
   if(x(cntr).eq.0d0) goto 666  ! when x=0 leave the loop    
   if(xb1.ge.x(cntr)) ixb1=cntr ! finds lower boundary
   if(xb2.ge.x(cntr)) ixb2=cntr ! finds upper boundary
!      write(*,*) x(cntr),cntr,xb1,xb2,ixb1,ixb2
 enddo
 

 
666 close(88) ! CLOSE THE binned_data file
 if ((ixb1.le.0).or.(ixb2.le.0)) then
!   write(*,*) "going there"
 goto 888 ! tell code that he's wrong when boundaries were not found
 endif
       count=0d0 ! set number of bins to 0
       ch2=0d0 ! set ch2 to 0
           
!          write(*,*) "POCITAM CH2",rin
       do ncount=ixb1,ixb2 ! Loop from lower to upper boundary
       
   	ch2=ch2+(y(ncount)-(a*x(ncount)+b))**2d0/(sigma(ncount)**2) !! chi square (ordinary) maybe weighted least squares?
!  	ch2=ch2+(y(ncount)-(a*x(ncount)+b))**2d0/abs((a*x(ncount)+b)) !! PEARSON - doesn't take sigmas into account
!         ks=abs((y(ncount)-(a*x(ncount)+b)))*100d0
!         if(ks.gt.ch2) ch2=ks
 	count=count+1 ! Add one bin 
       enddo
       
       if(count.gt.10) then
       ch2=-ch2*(cntr/count)**2d0 ! weigthing for the number of bins
       else
       ch2=-9d19 ! high negative number to tell the code that we don't want this
       endif

888 if ((xb1.ge.xb2).or.(ixb1.le.0).or.(ixb2.le.0)) ch2=-9d19 ! high negative number to tell the code that we don't want this
!   ressy=ch2*666d0 ! To increase the precision of the fit we multiply ch2 by 333 (100 to 1000 is safe, less than 100 gives worse fits, more than 1000 could be dangerous, too restrictive)
    ressy=ch2*(10d0**y(1)) ! Works for this setup for all CMOR data we tested, binwidth = 0.01 => weight = 333d0
!   
!   write(*,*) "ENDING SUBROUTINE TRANSIT", ressy
!  read(*,*)
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 ! === 7.0 CLOSE PROGRAM ===
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

END SUBROUTINE transit
!=======================================================================

END MODULE transitmod
