!*********************************************************************************
! This module contains empirical functions
!********************************************************************************

MODULE empirical
 
 use constant  ! constants

 implicit none

CONTAINS

!***********************************************************************************************
!THIS FUNCTION CALCULATES THE SATURATED WATER VAPOR PRESSURE.
!ITS-90 FORMULATIONS ARE USED HERE. (BOB HARDY, 1988)
!*********************************************************************************************
REAL(RLK) FUNCTION SVP(temp)

 implicit none

 real(RLK),dimension(8) :: gg
 real(RLK), intent(in) :: temp

 if (temp>473.15.or.temp<173.15) then
        print *,'Temperature is out of the range of saturation vapor &
        & pressure calculation.'
       stop
 endif

 gg(1)=-2.8365744E03*temp**(-2.0)
 gg(2)=-6.028076559E03*temp**(-1.0)
 gg(3)=1.954263612E01
 gg(4)=-2.737830188E-2*temp
 gg(5)=1.6261698E-5*temp**2
 gg(6)=7.0229056E-10*temp**3
 gg(7)=-1.8680009E-13*temp**4
 gg(8)=2.7150305*log(temp)
 !UNIT: DYNE/CM2 (IN G/CM/S UNIT FORM)
 SVP=exp(sum(gg))*10.0 

 !Empirical Formula used in C Parcel Model
 !real*8 :: c0,c1,c2,c3,c4,c5,c6
 !real*8 :: tempc
 !tempc=temp-273.15
 !if (tempc<-50.or.tempc>50) them
 !        print *,'Temperature is out of the range of saturation vapor &
 !       & pressure calculation.'
 !       print *,'temp',temp
 !       stop
 !endif

 !c0=6.10780
 !c1=0.443652
 !c2=1.428946e-2
 !c3=2.650648e-4
 !c4=3.03124e-6
 !c5=2.034081e-8
 !c6=6.136821e-11

 !SVP=c0+tempc*(c1+tempc*(c2+tempc*(c3+tempc*(c4+tempc*(c5+tempc*c6)))))
 !SVP=SVP*1000  

 !Empirical Formula used in F Parcel Model
 !This does not make any difference on Nd.

 !real*8 :: ts,sr,ar,ar1,br,cr,dw,er
 !ts=373.16
 !sr=3.0057166
 !ar=ts/temp
 !ar1=ar-1.
 !br=7.90298 * ar1
 !cr=5.02808 * log10(ar)
 !dw=1.3816e-7*(10.**(11.344*(ar1/ar))-1.)
 !er=8.1328e-3*(10.**(-3.49149*ar1)-1.)
 !SVP=1000*10.**(cr-dw+er+sr-br)

END FUNCTION SVP

!***********************************************************************************************
!THIS FUNCTION CALCULATES LATENT HEAT FROM EVAPORATION/CONDENSATION
!P762 Seinfeld and Pandis 2rd
!***********************************************************************************************
REAL(RLK) FUNCTION LTH(temp)

 implicit none

 real(RLK),intent(in) :: temp
 
 LTH=18.0*2.5*(273.15/temp)**(0.167+3.67E-4*temp)  !kJ mol-1
 !UNIT: ERG/MOL (IN G/CM/S UNIT FORM) 
 LTH=LTH*1.0E10

END FUNCTION LTH

!******************************************************************
!THIS FUNCTION CALCULATES SPECIFIC HEAT OF WATER
!P762 Seinfeld and Pandis 2rd
!******************************************************************
REAL(RLK) FUNCTION WSCH(temp)

 implicit none

 real(RLK), intent(in)  :: temp

 if (temp<233.0 .or. temp>308.0) then
        print *,'Specific heat calculation error'
        print *,'temp',temp
        STOP
 elseif (temp<=273.15) then
        WSCH=4.218+3.47E-4*(temp-273.15)**2
 else
        WSCH=4.175+1.3E-5*(temp-308.0)**2+1.6E-8*(temp-308.0)**4.0
 endif

 !UNIT: ERG/(G K) (IN G/CM/S UNIT FORM)
 WSCH=WSCH*1.0E7  ! in g cm s unit form

END FUNCTION WSCH
!**********************************************************************


END MODULE empirical


