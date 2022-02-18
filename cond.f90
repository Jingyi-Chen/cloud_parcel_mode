!*********************************************************************************
! This module provides the key equations
!********************************************************************************

MODULE condMD
 
 use constant 
 use empirical 
 use chemMD
 use aerospecMD
 use dynamMD


 implicit none

CONTAINS

SUBROUTINE COND(NEQ,TIME,YEQ,YDOT)

 implicit none

 !Argument
 real(RLK), intent(in) :: TIME
 integer,intent(in) :: NEQ
 real(RLK), intent(in),dimension(NEQ) :: YEQ
 real(RLK), intent(out),dimension(NEQ) :: YDOT

 integer :: wetbin,aerobin,ii,i0bin
 real(RLK) :: temp,press,lmr,vmr,nam
 real(RLK), dimension(:),allocatable :: pmass,prad,rhows,mola,concfrac,&
        &conc0,prad0,pmass0,bprad0
 integer, dimension(:),allocatable :: aerostate0
 real(RLK) :: rhos,rhosap,soluba,solubb,solubc,Ms,kap,rhois,&
        &solfracm,solfracv
 real(RLK) :: na0
 real(RLK) :: entrrt,entralpha,drytemp,dryvmr,vv
 
 real(RLK) :: mmvv,mmva,diffco,therco,mfpa,wmass,term1,term2,term3,diff,ther
 real(RLK) :: vpress,rhod,rhoa,kelvin,aw,satk,sat
 real(RLK) :: svpress,lv,surten,kelter,dvmr,dlmr
 real(RLK) :: dtemp,dpress,dnam,difflmr,cpma
 integer :: iaero,ientr,ivv,ichem
 real(RLK) :: na0_aie

 common /num_bin/ aerobin
 common /caseindex/ ichem,iaero,ientr,ivv

 wetbin=NEQ-5

 

 !READ IN 4 VARIABLES
 temp=YEQ(wetbin+1) 
 press=YEQ(wetbin+2)
 lmr=YEQ(wetbin+3)
 vmr=YEQ(wetbin+4)
 nam=YEQ(wetbin+5)
 !sat=YEQ(wetbin+6)

 allocate(pmass(aerobin),prad(aerobin),rhows(aerobin),mola(aerobin),&
        &concfrac(aerobin),aerostate0(aerobin),conc0(aerobin),prad0(aerobin),&
        &pmass0(aerobin),bprad0(aerobin+1))

 pmass=0; prad=0; rhows=0; mola=0; concfrac=0; aerostate0=0
 conc0=0; prad0=0; pmass0=0; bprad0=0

 CALL chem(ichem,rhos,rhosap,soluba,solubb,solubc,Ms,kap,rhois,&
        & solfracm,solfracv)

 na0=na0_aie
 CALL aerospec(iaero,aerobin,na0,rhos,solfracv,rhois,aerostate0,&
        &conc0,prad0,bprad0,pmass0)
 concfrac=conc0/na0
 i0bin=aerobin-wetbin+1
 
 CALL entrain(ientr,temp,press,entrrt,entralpha,drytemp,dryvmr)
 CALL vertvelo(ivv,vv)

 !MEAN MOLECULAR SPEED FOR WATER VAPOR
 mmvv=(8.0*gc*temp/pi/Mw)**0.5    
 !MEAN MOLECULAR SPEED FOR AIR
 mmva=(8.0*gc*temp/pi/Ma)**0.5 
 !GENERAL DIFFUSIVE COEFFICIENT
 diffco=0.211/(press/1013250.0)*(temp/273.15)**1.94    
 !GENERAL THERMAL CONDUCTIVITY OF AIR
 therco=(4.39+0.071*temp)*1.0E02  
 !AIR MEAN FREE PATH
 mfpa=mfpac*temp/press   
 
 
 wmass=0.0 !TOTAL WATER MASS
 lv=LTH(temp)/Mw !LATENT HEAT
 svpress=SVP(temp) !SATURATED VAPOR PRESSURE
 vpress=press*vmr/(vmr+eps) !VAPOR PRESSURE
 sat=vpress/svpress !SATURATION
 rhod=(press-vpress)/gcd/temp !DRY AIR DENSITY
 rhoa=press/gca/temp !AIR DENSITY

 !READ IN RADIUS AND MASS OF EACH PARTICLE BIN
 do ii=i0bin,aerobin
        pmass(ii)=YEQ(ii-i0bin+1)
        mola(ii)=pmass0(ii)*solfracm/Ms/(pmass(ii)-pmass0(ii))
        rhows(ii)=rhow+(rhosap-rhow)*mola(ii)/(rhosap/rhow/Ms+mola(ii))
        term3=0.75*(pmass(ii)-(1.0-solfracm)*pmass0(ii))/rhows(ii)/pi
        term3=term3+(1.0-solfracv)*prad0(ii)**3.0
        prad(ii)=(abs(term3)/term3)*(abs(term3)**(1.000/3.000))
 enddo

 dvmr=0.0
 dlmr=0.0

 !TENDENCY OF PARTICLE MASS
 do ii=i0bin,aerobin

        !MODIFIED COEFFICIENTS FOR MASS AND HEAT TRANSFER
        term1=prad(ii)/(prad(ii)+mfpa)
        diff=diffco/(term1+4.0*diffco/(alpham*mmvv*prad(ii)))
        ther=therco/(term1+4.0*therco*gc*temp/(alphat*mmva*press*cpa*Ma*prad(ii)))

        !EQUILIBRIUM VAPOR PRESSURE OVER DROPLET SURFACE
        !BASED ON KOHLER EQUATION
        aw=1.0/(1.0+kap*solfracv*prad0(ii)**3.0/(prad(ii)**3.0-prad0(ii)**3.0))
        surten=76.1-0.155*(temp-273.15)+2.17*mola(ii)
        kelter=2.0*surten/(rhow*gcw)
        kelvin=exp(kelter/temp/prad(ii))
        satk=kelvin*aw !EQUILIBRIUM SATURATION
        term2=gc*temp/(Mw*diff*svpress)+lv/(ther*temp)*(lv*Mw/gc/temp-1.0)*satk
        YDOT(ii-i0bin+1)=(sat-satk)*4.0*pi*prad(ii)/term2
        wmass=wmass+concfrac(ii)*YDOT(ii-i0bin+1)

 enddo

 !total aerosol concentration change
  dnam = -(1.0-entralpha)*nam*vv*entrrt

 !liquid water mixing ratio change due to vapro diffusion
 difflmr = wmass*nam*rhoa/rhod  !mass of condensed water per gram dry air

 !liquid water mixing ratio
 dlmr=difflmr+lmr*dnam/nam

 !water vapor mixing ratio
 dvmr = -dlmr-entrrt*vv*(vmr-dryvmr+lmr)

 !air temperature
 cpma=cpa+cpv*vmr
 dtemp = -grav*vv/cpma+&
        &(difflmr+entrrt*vv*entralpha*lmr)*lv/cpma-entrrt*vv*(temp-drytemp) 

 !air pressure
 dpress = -press*grav*vv/(gcd*temp)

 !ASSIGN THE VALUE TO NUMERICAL SOLVER VARIABLES
 YDOT(wetbin+1)=dtemp
 YDOT(wetbin+2)=dpress   
 YDOT(wetbin+3)=dlmr
 YDOT(wetbin+4)=dvmr
 YDOT(wetbin+5)=dnam
 !YDOT(wetbin+6)=dsat

 deallocate(pmass,prad,rhows,mola,concfrac,aerostate0,conc0,prad0,pmass0,bprad0)



RETURN
END SUBROUTINE COND

END MODULE condMD


