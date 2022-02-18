!****************************************************************
! This module provides dry aerosol spectra 
!****************************************************************

MODULE dropspecMD
  
 use constant
 use chemMD

 implicit none

 integer :: naero = 50

CONTAINS

!***********************************************************************************
!This subroutine provides the characteritics for soluble coomposition in
!aerosols.
!**********************************************************************************

SUBROUTINE dropspec(binnum,ast,drop_m0,drop_r0,aero_m0,aero_r0,&
                  &drop_radmin,drop_radmax)

  implicit none
  integer, intent(in) :: binnum
  integer, dimension(binnum),intent(inout) :: ast
  real(RLK), dimension(binnum),intent(out) :: drop_m0,drop_r0
  real(RLK), dimension(binnum),intent(in)  :: aero_r0,aero_m0
  real(RLK), intent(in) :: drop_radmin, drop_radmax
  integer :: iaero,ivv,ientr,ichem
  real(RLK) :: rhos,rhosap,soluba,solubb,solubc,Ms,kap,rhois,solfracm,solfracv

  real(RLK) :: drop_raddiff,term3
  real(RLK), dimension(binnum+1) :: drop_br0
  real(RLK), dimension(binnum) :: drop_r0_est,mola,rhows1
  integer :: ii
  common /caseindex/ ichem,iaero,ientr,ivv

  
  drop_raddiff = (log(drop_radmax)-log(drop_radmin))/binnum
  drop_br0(1)  = drop_radmin

  CALL chem(ichem, rhos,rhosap,soluba,solubb,solubc,Ms,kap,rhois,&
        & solfracm,solfracv)
  
  do ii=1,binnum

    drop_br0(ii+1) = exp(log(drop_br0(ii))+drop_raddiff)
    drop_r0_est(ii) = (drop_br0(ii) + drop_br0(ii+1))/2

    drop_m0(ii) = rhow*4.0/3.0*pi*drop_r0_est(ii)**3 ! assuming pure water
    mola(ii)=aero_m0(ii)*solfracm/Ms/(drop_m0(ii)-aero_m0(ii))
    rhows1(ii)=rhow+(rhosap-rhow)*mola(ii)/(rhosap/rhow/Ms+mola(ii))
    term3=0.75*(drop_m0(ii)-(1-solfracm)*aero_m0(ii))/rhows1(ii)/pi
    term3=term3+(1-solfracv)*aero_r0(ii)**3.0
    drop_r0(ii)=(abs(term3)/term3)*(abs(term3)**(1.000/3.000))

    ast(ii) = 2

  enddo

  !print *, 'Initializing Drops'
  !print *, 'dropmin', drop_radmin
  !print *, 'dropmax', drop_radmax
  !print *, 'pmass in INIDROP', drop_m0
  !print *, 'prad in INIDROP', drop_r0
  !print *, 'aero_m0',aero_m0

  !print *, ""

END SUBROUTINE dropspec

END MODULE dropspecMD