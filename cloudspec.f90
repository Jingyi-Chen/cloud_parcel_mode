!*********************************************************************************
! This Module provides cloud droplet spectra statistical measures
!********************************************************************************

MODULE cloudspecMD
 
 use constant

 implicit none

Contains

SUBROUTINE cloudspec(i0bin,aerobin,aerostate,concv,&
        &prad,pmass,cdconc,cdbeta,cdrmean,cdrsd,cddisp,acbin,&
        &lwc,re,rv,skew,kurt,rrainLN,ThrsFuncLN)

 implicit none

 integer,parameter :: neye=3
 integer, intent(in) :: i0bin,aerobin
 integer, dimension(aerobin), intent(in) :: aerostate
 real(RLK), dimension(aerobin), intent(in) :: prad,pmass,concv
 integer, dimension(neye),intent(out) :: acbin
 real(RLK),dimension(neye),intent(out) :: cdconc,cdbeta,cdrmean,cdrsd,cddisp,lwc
 real(RLK),dimension(neye),intent(out) :: skew,kurt,rrainLN,ThrsFuncLN,re,rv

 real(RLK),dimension(neye) :: cd_cm2, cd_cm3, cd_cm4, cd_m2, cd_m3, cd_m6
 real(RLK), dimension(neye) :: rainLN_m6, rainLN_m3, vc
 integer, dimension(neye) :: eye 
 integer :: ii, jj

!eye: 1-Reutter et al. 2009; 2-radius>1um; 3-radius>Kohler critical radius
 cdconc=0.0  !number concentration per volume
 cdrmean=0.0 !mean radius
 cdrsd=0.0  !standard deviation
 cddisp=0.0  !relative dispersion
 cdbeta=0.0

 !central moment of cloud droplet spectrum
 cd_cm2=0.0;cd_cm3=0.0;cd_cm4=0.0

 !moment of cloud droplet spectrum
 cd_m2=0.0; cd_m3=0.0; cd_m6=0.0

 acbin=0 !activated particle bin
 lwc=0.0 !liqudi water content
 vc=0.0 !total cloud water volume

 !ThrsFuncLR6=0
 ThrsFuncLN=0.0

 eye=0
 do ii=i0bin,aerobin
        
       CALL valueeye(eye,neye,aerostate(ii),prad(ii))

       do jj=1,neye
        if (eye(jj)==1) then
        
          acbin(jj)=acbin(jj)+1
          cdconc(jj)=cdconc(jj)+concv(ii)
          lwc(jj)=lwc(jj)+concv(ii)*pmass(ii)
          vc(jj)=vc(jj)+concv(ii)*4.0/3.0*pi*prad(ii)**3.0  !cm3/cm3

          cdrmean(jj)=cdrmean(jj)+prad(ii)*concv(ii)
          cd_m2(jj)=cd_m2(jj)+(prad(ii)**2.0)*concv(ii)
          cd_m3(jj)=cd_m3(jj)+(prad(ii)**3.0)*concv(ii)
          cd_m6(jj)=cd_m6(jj)+(prad(ii)**6.0)*concv(ii)
        else
          acbin(jj)=0
          cdconc(jj)=0.0
          lwc(jj)=0.0
          vc(jj)=0.0

          cdrmean(jj)=0.0
          cd_m2(jj)=0.0
          cd_m3(jj)=0.0
          cd_m6(jj)=0.0
        endif
       enddo
     
 enddo
 
 do jj=1,neye
    if(cdconc(jj)/=0.0) then
      cdrmean(jj)=cdrmean(jj)/cdconc(jj)
      rv(jj)=(vc(jj)/cdconc(jj)/4.0*3.0/pi)**(1.0/3.0)   ! volume mean radius
      re(jj)=cd_m3(jj)/cd_m2(jj)
    else 
      cdrmean(jj)=0.0
      rv(jj)=0.0
      re(jj)=0.0
    endif
 enddo

 cdbeta=re/rv

 eye=0 
 do ii=i0bin,aerobin
    
     CALL valueeye(eye,neye,aerostate(ii),prad(ii))

     do jj=1,neye
        if(eye(jj)==1) then
           cd_cm2(jj)=cd_cm2(jj)+((prad(ii)-cdrmean(jj))**2.0)*concv(ii)
           cd_cm3(jj)=cd_cm3(jj)+((prad(ii)-cdrmean(jj))**3.0)*concv(ii)
           cd_cm4(jj)=cd_cm4(jj)+((prad(ii)-cdrmean(jj))**4.0)*concv(ii)
        else
           cd_cm2(jj)=0.0
           cd_cm3(jj)=0.0
           cd_cm4(jj)=0.0
        endif
     enddo
 
 enddo
        
 do jj=1,neye
    if(cdconc(jj)/=0.0) then
        cdrsd(jj)=sqrt(cd_cm2(jj)/cdconc(jj))
        cddisp(jj)=cdrsd(jj)/cdrmean(jj)
        skew(jj)=cd_cm3(jj)/cdconc(jj)/(cdrsd(jj)**3.0)
        kurt(jj)=cd_cm4(jj)/cdconc(jj)/(cdrsd(jj)**4.0)
    else 
        cdrsd(jj)=0.0
        cddisp(jj)=0.0
        skew(jj)=0.0
        kurt(jj)=0.0
    endif
 enddo

 ! Liu et al. (2004) eq.5 and eq.7
 do jj=1,neye 
  !rrainLR6(jj)=23.72*1.0e-7/lwc(jj)**(1.0/6.0)/(rv(jj)*1.12)**(1.0/2.0)
  rrainLN(jj)=2.8522e-6*cdconc(jj)**(1.0/6.0)/lwc(jj)**(1.0/3.0)
 enddo
 
 !rainLR6_m6=0; rainLR6_m3=0
 rainLN_m6=0.0; rainLN_m3=0.0

 do ii=i0bin,aerobin

    do jj=1,neye

!       if (prad(ii)>=rrainLR6(jj)) then
!        rainLR6_m6(jj)=rainLR6_m6(jj)+(prad(ii)**6.0)*concv(ii)
!        rainLR6_m3(jj)=rainLR6_m3(jj)+(prad(ii)**3.0)*concv(ii)
!       else
!        rainLR6_m6(jj)=0
!        rainLR6_m3(jj)=0
!       endif

       if (prad(ii)>=rrainLN(jj)) then
        rainLN_m6(jj)=rainLN_m6(jj)+(prad(ii)**6.0)*concv(ii)
        rainLN_m3(jj)=rainLN_m3(jj)+(prad(ii)**3.0)*concv(ii)
       else
        rainLN_m6(jj)=0.0
        rainLN_m3(jj)=0.0
       endif

    enddo

 enddo


 do jj=1,neye

!    ThrsFuncLR6(jj)=rainLR6_m6(jj)*rainLR6_m3(jj)/cd_m6(jj)/cd_m3(jj)

    ThrsFuncLN(jj)=rainLN_m6(jj)*rainLN_m3(jj)/cd_m6(jj)/cd_m3(jj)

 enddo

 END SUBROUTINE cloudspec

 SUBROUTINE valueeye(eye,neye,aerostatei,pradi)

  implicit none

  integer, intent(in) :: neye
  integer, dimension(neye), intent(out) :: eye
  real(RLK), intent(in) :: pradi
  integer, intent(in) :: aerostatei

      if(aerostatei==2.or.aerostatei==4) then
         eye(1)=1 
      else
         eye(1)=0 
      endif
 
      if (pradi>1.0E-04) then
         eye(2)=1
      else
         eye(2)=0
      endif

      if (aerostatei==3.or.aerostatei==4) then
         eye(3)=1
      else
         eye(3)=0
      endif

 END SUBROUTINE valueeye

END MODULE cloudspecMD
