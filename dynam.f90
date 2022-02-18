!****************************************************************
! This module includes dynamical properties: vv and entrianment 
!****************************************************************

 MODULE dynamMD
  
 use constant
 use empirical

 implicit none

CONTAINS

 SUBROUTINE binnum(aerobin)

   implicit none
    real(RLK), intent(out) :: aerobin

    aerobin=20


 END SUBROUTINE binnum

 SUBROUTINE vertvelo(ivv,v0)
 
    implicit none
    integer, intent(in) :: ivv
    real(RLK), intent(out) :: v0

    if (ivv.eq.1) then
      v0=0

    else if (ivv.eq.2) then
      v0=200.0

    endif


 END SUBROUTINE vertvelo


 SUBROUTINE entrain(ientr,temp,press,entrrt,entralpha,drytemp,dryvmr)

    implicit none
    integer, intent(in) :: ientr
    real(RLK), intent(in) :: temp, press
    real(RLK), intent(out) :: entrrt, entralpha,drytemp,dryvmr
    real*8 :: dryrh,dryvpr

    entrrt = 0.0 + 0.1 * ( real(ientr,8) - 1.0 )
    entralpha = 0.0
    drytemp = temp - 1.0
    dryrh = 0.8

    dryvpr=SVP(drytemp)*dryrh 
    dryvmr=dryvpr*eps/(press-dryvpr)
   


 END SUBROUTINE entrain


END MODULE dynamMD
