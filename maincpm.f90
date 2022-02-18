!*********************************************************************************
! Cloud parcel model by Jingyi Chen
! Default units are g/cm/s/K except notice.
!
! This is the main program.
!********************************************************************************

PROGRAM main

    use cpmMD
    use constant
    use aerospecMD
    use dropspecMD
    use chemMD
    use dynamMD

    implicit none

    integer :: aerobin = 20
    integer :: iaero,ivv,ientr,ichem,ipredrop,ioi
    real(RLK) :: dt_les, dt_cpm
    real(RLK) :: sat_read,temp_read,press_read,temp0,press0,sat0,vv
    real(RLK) :: drop_radmin_read,drop_radmax_read,drop_radmin,drop_radmax
    real(RLK) :: na0
    real(RLK) :: rhos,rhosap,soluba,solubb,solubc,Ms,kap,rhois,solfracm,solfracv
    real(RLK) :: entrrt,entralpha,drytemp,dryvmr
    integer, dimension(:), allocatable :: aerostate
    real(RLK), dimension(:), allocatable :: conc0,prad0,pmass0,pmass,prad
    real(RLK), dimension(:), allocatable :: bprad0

    common /num_bin/ aerobin
    common /caseindex/ ichem,iaero,ientr,ivv

    allocate(aerostate(aerobin), conc0(aerobin), prad0(aerobin))
    allocate(pmass0(aerobin), pmass(aerobin), prad(aerobin))
    allocate(bprad0(aerobin))

    !###################################
    sat_read = 95 ! %
    temp_read = 295 ! K
    press_read = 1000 ! hPa

    ! the units in cpm are g,cm,s
    temp0 = temp_read
    press0 = press_read*1.0e3
    sat0 = sat_read*0.01

    !###################################
    ichem = 1 !1-sulfate ammonium; 2-sea salt; 3-organic carbon
    CALL chem(ichem,rhos,rhosap,soluba,solubb,solubc,Ms,kap,rhois,&
        & solfracm,solfracv)
   
    !###################################
    iaero=1
    CALL aerospec(iaero,aerobin,na0,rhos,solfracv,rhois,aerostate,&
        &conc0,prad0,bprad0,pmass0)

    ivv=2
    CALL vertvelo(ivv,vv)

    ientr=1
    CALL entrain(ientr,temp0,press0,entrrt,entralpha,drytemp,dryvmr)

    !###################################
    ! initialize water onto it

    ipredrop = 0

    if (ipredrop==1) then

        drop_radmin_read = 20 ! micro
        drop_radmax_read = 30 ! micro
        drop_radmin  = drop_radmin_read * 1.0e-4
        drop_radmax  = drop_radmax_read * 1.0e-4

        call dropspec(aerobin,aerostate,pmass,prad,pmass0,prad0,&
                  &drop_radmin,drop_radmax)
    else
       ! initialize with dry aerosols.
        pmass = pmass0
        prad = prad0
        
    endif
    !################################### 
    dt_les = 180 ! second
    dt_cpm = 1.0 ! second

    print *, 'before cpm'
    print *, 'prad (micro)', prad*1.0e4
    print *, 'aerostate', aerostate
    print *, ' '
   

    ! call the cpm subroutine
    ioi = 1
    CALL cpm(dt_les,dt_cpm,temp0,press0,sat0,na0,vv,entrrt,drop_radmin,drop_radmax,&
            & rhos,rhosap,soluba,solubb,solubc,Ms,kap,rhois,solfracm,solfracv,&
            & aerostate,conc0,prad0,pmass0,prad,pmass,ioi) 

    print *, 'after cpm'
    print *, 'prad (micro)', prad*1.0e4
    print *, 'aerostate', aerostate
    print *, ' '


END PROGRAM main


