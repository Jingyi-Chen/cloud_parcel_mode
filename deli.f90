!*********************************************************************************
! The module includes deliquescences and efflorescences
!********************************************************************************

Module deliMD

 use constant  
 use empirical 
 use chemMD

 implicit none

Contains

SUBROUTINE deli(press,rhod,aerobin,prad0,pmass0,concv,temp,&
        &sat,lmr,vmr,prad,pmass,aerostate)

  implicit none
 
  real(RLK), intent(in) :: press,rhod
  integer, intent(in) :: aerobin
  real(RLK), dimension(aerobin), intent(in)  :: prad0, pmass0,concv
  real(RLK), intent(inout) :: temp,sat,lmr,vmr
  real(RLK), dimension(aerobin), intent(inout) :: prad,pmass
  integer, dimension(aerobin),intent(inout) :: aerostate

  integer :: iaero,ivv,ientr,ichem
  real(RLK) :: surtenc, kelterc, solub, aw0, surten, kelter, rhows
  real(RLK) :: wmassterm, lv, svpress, vpress
  real(RLK), dimension(aerobin) :: satdeliq, wmassat, radsatsolu
  real(RLK), dimension(aerobin) :: satefflo, mola

  integer :: ii
  real(RLK) :: rhos, rhosap,soluba,solubb,solubc,Ms,kap,rhois
  real(RLK) :: solfracm, solfracv
  common /caseindex/ ichem,iaero,ientr,ivv

  call chem(ichem,rhos,rhosap,soluba,solubb,solubc,Ms,kap,rhois,&
        & solfracm,solfracv)

  surtenc=76.1-0.155*(temp-273.15) !surface tension
  kelterc=2.0*surtenc/(rhow*gcw) !coefficient in Kelvin Term 

  !Deliquesence & Efflorescence: (Ref: CHEN 1994)
  do ii=1,aerobin
   !Deliquesence
   solub=soluba+solubb*temp+solubc*temp**2.0
   solub=1.0/(1.0+1.0/solub)
   aw0=(1.0-solub)/(1.0-solub+kap*Ms*rhow*solub/rhos/Mw)
   satdeliq(ii)=aw0*exp(kelterc/(temp*prad0(ii)))

   !Efflorescence
   !water mass of saturated solution
   wmassat(ii)=(1.0/solub-1.0)*Mw/Ms*pmass0(ii)*solfracm
   !molality of soluble solution
   mola(ii)=1.0/(1.0/solub-1.0)/Mw
   !solution surface tension: 2.17-(NH4)2SO4
   surten=76.1-0.155*(temp-273.15)+2.17*mola(ii)
   kelter=2.0*surten/(rhow*gcw)
   !solution density
   rhows=rhow+(rhosap-rhow)*mola(ii)/(rhosap/rhow/Ms+mola(ii))
   !eq. 19 in CHEN,1994
   radsatsolu(ii)=rhos/rhows*(Mw/Ms*(1.0/solub-1.0)+1.0)
   radsatsolu(ii)=(solfracv*radsatsolu(ii)+(1.0-solfracv))*prad0(ii)**3.0
   radsatsolu(ii)=radsatsolu(ii)**(1.0/3.0)
   !eq. 20 IN CHEN, 1994
   satefflo(ii)=aw0*exp(kelter/temp/radsatsolu(ii))
  enddo

  wmassterm=0.0 !total water mass change
  do ii=1,aerobin

        if (aerostate(ii)==0)   then

          if (sat>=satdeliq(ii)) then
            aerostate(ii)=1
            prad(ii)=radsatsolu(ii) 
            pmass(ii)=pmass(ii)+wmassat(ii)
            lmr=lmr+wmassat(ii)*concv(ii)/rhod
            vmr=vmr-wmassat(ii)*concv(ii)/rhod
            wmassterm=wmassterm+wmassat(ii)
          else
            continue
          endif

        elseif (aerostate(ii)==1) then

          if (sat<=satefflo(ii)) then
            aerostate(ii)=0
            prad(ii)=prad0(ii)
            pmass(ii)=pmass0(ii)
            lmr=lmr-wmassat(ii)*concv(ii)/rhod
            vmr=vmr+wmassat(ii)*concv(ii)/rhod
            wmassterm=wmassterm-wmassat(ii)
          else
            continue
          endif

        endif

   enddo

        lv=LTH(temp)/Mw
        temp=temp+lv*wmassterm/(cpa*rhod)
        svpress=SVP(temp)
        vpress=vmr*press/(eps+vmr)
        sat=vpress/svpress

 END SUBROUTINE deli

END MODULE deliMD
