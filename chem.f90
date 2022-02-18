!****************************************************************
! This module provides aerosol chemical compositions 
!****************************************************************

 MODULE chemMD

 use constant
  
 implicit none


CONTAINS

 SUBROUTINE chem(ichem,rhos,rhosap,soluba,solubb,solubc,Ms,kap,rhois,&
        & solfracm,solfracv)

  implicit none

  integer, intent(in) :: ichem
  real(RLK), intent(out) :: rhos,rhosap,soluba,solubb,solubc,&
        & Ms,kap,rhois,solfracm,solfracv 
  integer :: solid,inid

 !  soid : soluble part chemical composition: 
 !      1-sulfate ammonium; 2-sea salt; 3-organic carbon
 !  inid : insoluble part chemical composition:
 !      1-dust; 
 !  sfr : soluble fraction (mass fraction)
   solid = ichem
   inid = 1
   solfracm = 1.0

   call solcomp(solid,rhos,rhosap,soluba,solubb,solubc,Ms,kap)
   call insolcomp(inid,rhois)
      
   solfracv=solfracm*rhois/(rhos+solfracm*(rhois-rhos))

 END SUBROUTINE chem


!*************************************************************
!This subroutine provides parameters of soluble coomposition 
!**************************************************************
SUBROUTINE solcomp(idin,rhost,rhosapt,solubat,solubbt,solubct,Mst,kapt)

 implicit none
 
 integer, intent(in) :: idin
 real(RLK), intent(out) :: rhost,rhosapt,solubat,solubbt,solubct,Mst,kapt

   select case (idin)

        !SULFATE AMMONIUM
        case (1)
        rhost=1.769
        rhosapt=2.093
        solubat=0.1149
        solubbt=-4.489e-4
        solubct=1.385e-6
        Mst=132.14
        kapt=0.61

        !Sodium Chloride
        case (2) 
        rhost=2.165
        rhosapt=2.726
        solubat=0.1805
        solubbt=-5.310e-4
        solubct=9.965e-7
        Mst=58.44
        kapt=1.28

        !biomass burning aerosols
        case (3)
        rhost=1.662
        rhosapt=2.000
        solubat=0.1149
        solubbt=-4.489e-4
        solubct=1.385e-6
        Mst=111.00
        kapt=0.1
        case default
        print *,'THIS SALT DOES NOT EXIST IN THIS MODEL.'
   end select

END SUBROUTINE solcomp

!***************************************************************
!This subroutine provides parameters of insoluble composition.
!**************************************************************
SUBROUTINE insolcomp(id,rhois1)
 
 implicit none
 integer, intent(in) :: id
 real(RLK), intent(out) :: rhois1

 select case (id)
        case (1)   !dust
           rhois1=2.6
           

        case (2)   !soot BC
           rhois1=2.0
         
        
        case default
           print *,'INPUT INSOLUBLE COMPOSITION DOES NOT EXIT.'

 end select

END SUBROUTINE insolcomp

END MODULE chemMD
