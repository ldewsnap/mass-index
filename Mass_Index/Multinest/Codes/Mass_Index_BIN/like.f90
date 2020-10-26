MODULE like

use params
use transitmod
!use radialmod

implicit none
      
contains
      
!=======================================================================

SUBROUTINE slikelihood(R,slhood)
         
	implicit none
      
	REAL(8), DIMENSION(nest_nPar) :: R
        REAL(8) :: slhood, loglike
	INTEGER :: i

        ! Scaling of the parameters from hypercube
	DO i = 1, sdim
          IF( Rflag(i) .EQ. 0 .OR. Rflag(i) .EQ. 2 ) THEN ! Modes 0 and 2
            ! Uniform: Rmax = max, Rmin = min
	    R(i) = Rmin(i) + (Rmax(i)-Rmin(i))*R(i)
          ELSE IF( Rflag(i) .EQ. 3 ) THEN ! Mode 3
            ! Gaussian: Rmax = mean, Rmin = stdev
            R(i) = Rmax(i) + roottwo*Rmin(i)*inverf(-1.0D0+2.0D0*R(i))
          ELSE IF( Rflag(i) .EQ. 4 ) THEN ! Mode 4
            ! Jeffrey's: Rmax = max, Rmin = min
            R(i) = ( Rmax(i)**R(i) )*( Rmin(i)**(1.0D0-R(i)) )
          ELSE IF( Rflag(i) .EQ. 5 ) THEN ! Mode 5
            ! Modified Jefrey's: Rmax = max, Rmin = inflection point
            R(i) = -( Rmin(i)**(1.0D0-R(i)) )*( Rmin(i)**R(i) - ( Rmin(i)+Rmax(i) )**R(i) )
          END IF
	END DO

        ! Call transit to get chi^2
        call models(R,loglike,0)

	slhood = loglike
        !!write(*,*) 'R(1) = ',R(1),' yielding chi2 = ',chi2

END SUBROUTINE slikelihood
      
!=======================================================================

!! ======================================================================
SUBROUTINE models(Rvec,loglike,showpri)

!implicit none

 REAL(8), DIMENSION(nest_nPar), INTENT(IN) :: Rvec   ! Fitted-parameter vector
 INTEGER :: showpri, showrv
! REAL(8), DIMENSION(nplen) :: resP
 REAL(8) :: resP
 REAL(8) :: loglike, loglikeP, loglikeP_lc, loglikeP_sc, loglikeR
 INTEGER :: i, j, nplen_lc, nplen_sc
 REAL(8) :: jitter
 REAL(8) :: chi2RV

 ! === Call transit to primary transit ===

!call transit(Rvec,nplen,resP,&
!              showpri,tp,fp,sigfp,epochp_seq,fpwei,sigfpwei,&
!              NresamP,integ_bigP,&
!              .FALSE.)
  resp=-9d19 
!    WRITE(*,*) "CALLING TRANSIT",rvec,resP
 call transit(Rvec,resP)
!    WRITE(*,*) "AFTER TRANSIT",rvec,resP
 
 loglikeP = 0.0D0
!  if(sflag.eq.1.or.sflag.eq.2) then
   
       loglikeP = loglikeP + resP
      
!     END DO
     loglikeP = 0.5D0*loglikeP
!   write(*,*) loglikeP, resP
    
   
! if(loglikeP.gt.-1e-6) write(*,*) rvec,resP,loglikeP
!read(*,*)

! Part related to resampling was deleted

 ! === Call radial ===

! This part was deleted because there was some problem with "Vsyslen"

 loglikeR = 0.d0

 ! === Sum up all loglikes (product of likes) ===
 loglike = loglikeP + loglikeR

END SUBROUTINE models
! ======================================================================

! ======================================================================
FUNCTION inverf(x)

 implicit none

 REAL(8) :: x
 REAL(8), PARAMETER :: awil = 0.14001228868666646D0
 REAL(8), PARAMETER :: bwil = 4.546884979448289D0
 REAL(8) :: factor, xsq, inverf

 IF( x .LT. 0.0D0 ) THEN
  factor = -1.0D0
 ELSE
  factor = 1.0D0
 END IF

 xsq = 1.0D0 - x**2
 x = bwil + 0.5D0*DLOG(xsq)
 x = DSQRT( x**2 - (DLOG(xsq)/awil) ) - x
 inverf = factor*DSQRT(x)

END FUNCTION
! ======================================================================

END MODULE like

