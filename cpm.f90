!*********************************************************************************
! Cloud parcel model by Jingyi Chen
! Default units are g/cm/s/K except notice.
!
! This Module calculate the updrafts of cloud parcel.
!********************************************************************************

MODULE cpmMD
 
 use constant 
 use empirical 
 use DVODE_F90_M 
 use deliMD
 use condMD
 use cloudspecMD

 implicit none

CONTAINS

SUBROUTINE cpm(dt_les,dt_cpm,temp0,press0,sat0,na0,vv,entrrt,drop_radmin,drop_radmax,&
            & rhos,rhosap,soluba,solubb,solubc,Ms,kap,rhois,solfracm,solfracv,&
            & aerostate,conc0,prad0,pmass0,prad,pmass,ioi,io_dsd) 

      ! passing from main
      integer :: aerobin
      common /num_bin/ aerobin

      real(RLK),intent(in) :: dt_les,dt_cpm,temp0,press0,sat0,na0,vv,entrrt,drop_radmin,drop_radmax
      real(RLK),intent(in) :: rhos,rhosap,soluba,solubb,solubc,Ms,kap,rhois,solfracm,solfracv
      integer, dimension(aerobin), intent(inout) :: aerostate
      real(RLK), dimension(aerobin),intent(in) :: pmass0,prad0
      real(RLK), dimension(aerobin),intent(inout) :: pmass,prad,conc0
      integer, intent(in) :: ioi,io_dsd

      ! other declarations
      integer   :: instabind,wetbin,ii,i0bin,cpmv_unit,cpmr_unit
      integer   :: totstep,step 
      real(RLK) :: temp,press,sat,height,na,nam
      real(RLK) :: svpress,vpress,vmr,lmr,totalwmr,rhod0,rhod,rhoa0,rhoa
      real(RLK), dimension(aerobin):: concv,concm,concfrac,pmass_drop,prad_drop
      integer, dimension(aerobin) :: aerostate_drop
      real(RLK) :: surten,kelter,rcri_R09
      real(RLK), dimension(aerobin) :: mola,rhows,aw,kelvin,satk,rcri_kohler

      integer, parameter :: neye=3
      integer, dimension(neye):: acbin
      real(RLK), dimension(neye) :: cdconc,cdrmean,cdrsd,cddisp,lwc,rv,re,skew,kurt
      real(RLK), dimension(neye) :: cdbeta,rrainLN,ThrsFuncLN
      CHARACTER(len=20) :: cpmv_fn, cpmr_fn

      ! for numerical solver use
      integer :: NEQ 
      real(RLK), dimension(:),allocatable :: ATOL,YEQ 
      real(RLK)  :: TOUT,TIME,RTOL,dt 
      integer :: ISTATE,ITASK 
      real(RLK), dimension(22) :: RSTATS 
      integer, dimension(31) :: ISTATS 
      TYPE(VODE_OPTS),save :: OPTIONS

      if (ioi.eq.1) then
         cpmv_unit=21
         cpmv_fn='cpmv.txt'
         open(unit=cpmv_unit,file=cpmv_fn)
         if (io_dsd.eq.1) then
             cpmr_unit=31
             cpmr_fn='cpmr.txt'
             open(unit=cpmr_unit,file=cpmr_fn)
         endif
      endif 

      totstep = dt_les/dt_cpm
      dt = dt_cpm
      concfrac=conc0/na0

      pmass_drop = pmass
      prad_drop = prad
      aerostate_drop = aerostate


      instabind=0 
      NUMCHECK: do ! for numerical instability check 

            temp=temp0; press=press0; sat=sat0
            prad=prad_drop; pmass=pmass_drop; concv=conc0
            na=na0; aerostate=aerostate_drop; height = 0

            svpress= SVP(temp)  !saturation water vapor pressure
            vpress=svpress*sat   !water vapor pressure
            vmr=vpress*eps/(press-vpress) !water vapor mixing ratio 
            lmr=0.0  !liquid water mixing ratio (g/g)
            totalwmr=vmr+lmr  !total water mixing ratio
            rhod0=(press-vpress)/gcd/temp  !dry air density
            rhod=rhod0
            rhoa0=press/gca/temp  !air density
            rhoa=rhoa0
            nam=na/rhoa   ! total particle number concentration per mass
            concm=concv/rhoa !concentration per mass in each bin

            TIME=0.0 !inticial time for numerical solver use (s)

            do ii=1,aerobin
                  lmr = lmr + concv(ii)*(pmass_drop(ii)-pmass0(ii))/rhod
             enddo


            UPLOOP: do step=int(1,8),int(totstep,8)   

                  CALL deli(press,rhod,aerobin,prad0,pmass0,concv,temp,sat,lmr,vmr,&
                        &prad,pmass,aerostate)
   
                  NEQ=5 !number of differencial equations
                  wetbin=0 !initialization of number of wet particle bins
                  do ii=1,aerobin
                        if (aerostate(ii)>0) then
                              NEQ=NEQ+1
                              wetbin=wetbin+1
                        endif
                  enddo
                  i0bin=aerobin-wetbin+1 ! the start bin for wet particles
                  

                  allocate(YEQ(NEQ),ATOL(NEQ))

                  !set tolerance for numerical solver
                  RTOL=1.0e-10
                  do ii=i0bin,aerobin
                  ATOL(ii-i0bin+1)=abs(pmass(ii))*1.0E-8
                  enddo
                  ATOL(wetbin+1)=abs(temp)*1.0E-8
                  ATOL(wetbin+2)=abs(press)*1.0E-8
                  if (lmr==0.0) then
                  ATOL(wetbin+3)=1.0E-27
                  else
                  ATOL(wetbin+3)=abs(lmr)*1.0E-8
                  endif
                  ATOL(wetbin+4)=abs(vmr)*1.0E-8
                  ATOL(wetbin+5)=abs(nam)*1.0E-8
                  !ATOL(wetbin+6)=abs(sat)*1.0E-8

                  OPTIONS=SET_NORMAL_OPTS(DENSE_J=.TRUE.,RELERR=RTOL,ABSERR_VECTOR=ATOL)

                  !put variables in to solver
                  YEQ(wetbin+1)=temp   !air temperature
                  YEQ(wetbin+2)=press  !air pressure
                  YEQ(wetbin+3)=lmr  !liquid water mixing ratio
                  YEQ(wetbin+4)=vmr  !water vapor mixing ratio
                  YEQ(wetbin+5)=nam  !total particle number concentration per mass
                  !YEQ(wetbin+6)=sat  !saturation ratio, optional

                  do ii=i0bin,aerobin
                        YEQ(ii-i0bin+1)=pmass(ii)
                  enddo

                  !numerical solver set-up
                  TOUT=TIME+dt
                  ITASK=1
                  ISTATE=1

                  CALL VODE_F90(COND,NEQ,YEQ,TIME,TOUT,ITASK,ISTATE,OPTIONS)
 
                  CALL GET_STATS (RSTATS, ISTATS)  ! check whether convergence/iteration successful
                  if (ISTATS(21)>0.or.ISTATE/=2) then
                        instabind=1
                        print *,'Instability Step',step
                        exit
                  else 
                        instabind=0
                  endif

                  !PUT VARIABLES BACK FROM NUMERICAL SOLVER:      
                  temp=YEQ(wetbin+1)
                  press=YEQ(wetbin+2)
                  lmr=YEQ(wetbin+3)
                  vmr=YEQ(wetbin+4)
                  nam=YEQ(wetbin+5)
                  sat=YEQ(wetbin+6)

                  totalwmr=vmr+lmr 
                  svpress=SVP(temp)   !saturated vapor pressure
                  vpress=press*vmr/(vmr+eps) !vapor pressure
                  sat=vpress/svpress   !saturation ratio
                  rhod=(press-vpress)/gcd/temp !dry air density 
                  rhoa=press/gca/temp !air density

                  print *,'step = ', step
                  print *,'temp = ', temp
                  print *,'sat = ', sat
                  print *,'lmr = ', lmr
                  print *,'height = ', height

                  do ii=1,aerobin
                        concm(ii)=concfrac(ii)*nam
                        concv(ii)=concm(ii)*rhoa
                  enddo
                  na=nam*rhoa

                  do ii=i0bin,aerobin

                        pmass(ii)=YEQ(ii-i0bin+1)
                        mola(ii)=pmass0(ii)*solfracv/Ms/(pmass(ii)-pmass0(ii))
                        rhows(ii)=rhow+(rhosap-rhow)*mola(ii)/(rhosap/rhow/Ms+mola(ii))
                        prad(ii)=0.75*(pmass(ii)-(1.0-solfracm)*pmass0(ii))/rhows(ii)/pi
                        prad(ii)=(prad(ii)+(1.0-solfracv)*prad0(ii)**3.0)**(1.0/3.0)

                        aw(ii)=1.0/(1.0+kap*solfracv*prad0(ii)**3.0/(prad(ii)**3.0-prad0(ii)**3.0))
                        surten=76.1-0.155*(temp-273.15)+2.17*mola(ii)
                        kelter=2.0*surten/(rhow*gcw)
                        kelvin(ii)=exp(kelter/temp/prad(ii))
                        satk(ii)=kelvin(ii)*aw(ii) !equilibrium saturation ratio

                        !critical radius as same as Retter et al. (2009) 
                        rcri_R09=2.0*kelter/temp/log(sat)/3.0
                        !critical radius for each particle based on Kohler Curve
                        rcri_kohler(ii)=(3.0*kap*solfracv*prad0(ii)**3.0/kelter*temp)**(1.0/2.0)

                        if (sat>=1.0) then
                              !cloud droplets
                              if (prad(ii)>=rcri_R09.and.prad(ii)<rcri_kohler(ii)) then
                                    aerostate(ii)=2
                              endif

                              if (prad(ii)>=rcri_kohler(ii).and.prad(ii)<rcri_R09) then
                                    aerostate(ii)=3
                              endif

                              if (prad(ii)>=rcri_kohler(ii).and.prad(ii)>=rcri_R09) then
                                    aerostate(ii)=4
                              endif

                              !interstitial particles
                              if (prad(ii)<rcri_R09.and. prad(ii)<rcri_kohler(ii)) then
                                    aerostate(ii)=1
                              endif
                              
                        else 
                              !interstitial particles if sub-saturated
                              aerostate(ii)=1  
                        endif

                  
                  enddo


                  CALL cloudspec(i0bin,aerobin,aerostate,concv,&
                        &prad,pmass,cdconc,cdbeta,cdrmean,cdrsd,cddisp,acbin,&
                        &lwc,re,rv,skew,kurt,rrainLN,ThrsFuncLN)

                  print *,'nc = ', cdconc
                  print *,'na = ', na
                  print *,''


                  if (ioi.eq.1) then
                  
                      write(cpmv_unit,'(I13.3,20E17.7E3,2I13.3)') step,TIME,height,temp,press,rhoa,&
                           &sat,vv,na,nam,entrrt,vmr,lmr,totalwmr,lwc(2),cdconc(2),cdrmean(2),lwc(3),cdconc(3),&
                           &cdrmean(3),dt,wetbin,aerobin
                           
                      if (io_dsd.eq.1) then
                          write(cpmr_unit,'(I13.3,2001E15.7E3)') step,(prad(ii),ii=1,aerobin,1)
                          write(cpmr_unit,'(I13.3,2001E15.7E3)') step,(concv(ii),ii=1,aerobin,1)
                          write(cpmr_unit,'(I13.3,2001I13.3)') step,(aerostate(ii),ii=1,aerobin,1)
                      endif
                      
                  endif

                  height = height + dt*vv

                  !OUTPUT_OUTPUT  

                  deallocate(YEQ,ATOL)
                  
                  ! stop condition
                  if(height>=60000.) then
                        exit
                  endif

            end do UPLOOP !ENDDO FOR EACH LOOP

      if (ioi.eq.1) then
         close(cpmv_unit)
         if (io_dsd.eq.1) then
             close(cpmr_unit)
         endif
      endif

      if (instabind==1) then
            deallocate(YEQ,ATOL)
            dt=dt*0.5
            exit

      else
         exit
      endif

 end do NUMCHECK  !enddo for numerical stability check
   
 
END SUBROUTINE cpm

END MODULE cpmMD


