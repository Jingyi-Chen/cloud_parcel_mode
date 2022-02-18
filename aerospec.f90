!****************************************************************
! This module provides dry aerosol spectra 
!****************************************************************

MODULE aerospecMD
  
 use constant

 implicit none

 integer :: naero = 50

CONTAINS

 SUBROUTINE aerospec(iaero,aerobin,na0,rhos,solfracv,rhois,aerostate0,&
        &conc0,prad0,bprad0,pmass0)
 
    implicit none
    integer, intent(in) :: aerobin,iaero
    real(RLK), intent(in) :: rhos,solfracv,rhois
    real(RLK), intent(out) :: na0
    integer, dimension(aerobin),intent(out) :: aerostate0
    real(RLK), dimension(aerobin),intent(out) :: conc0,prad0,pmass0
    real(RLK), dimension(aerobin+1),intent(out) :: bprad0

    real(RLK) :: dev,geo,radmin,radmax,raddiff,coe1,coe2
    integer :: ii
    
    dev=1.5
    geo=0.06e-4

    if (iaero .eq. 1) then
        na0 = 200
    endif

    radmin=0.01e-4;radmax=1.0e-4
    raddiff=(log(radmax)-log(radmin))/real(aerobin,8)
    bprad0(1)=radmin

    do ii=1,aerobin
        bprad0(ii+1)=exp(log(bprad0(ii))+raddiff)
        prad0(ii)=(bprad0(ii)+bprad0(ii+1))/2.0
        
        coe1=na0/(2.0*pi)**0.5/log(dev) 
        coe2=2.0*log(dev)**2.0
        conc0(ii)=(coe1/prad0(ii))*exp(-(log(prad0(ii))-log(geo))**2.0/coe2)!#/CM3/CM
        conc0(ii)=conc0(ii)*(bprad0(ii+1)-bprad0(ii)) ! #/CM3

        pmass0(ii)=rhos*solfracv*4.0/3.0*pi*prad0(ii)**3.0
        pmass0(ii)=pmass0(ii)+rhois*(1.0-solfracv)*4.0/3.0*pi*prad0(ii)**3.0

    enddo
 
    !aerosol states
    do ii=1,aerobin
        aerostate0(ii)=0
    enddo


 END SUBROUTINE aerospec

END MODULE aerospecMD
