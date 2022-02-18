!*********************************************************************
!Constant definitions for cloud parcel model
!********************************************************************

 MODULE CONSTANT

 implicit none

!CONSTANTS
 integer,parameter :: RLK=KIND(1.0D0)
 real(RLK),parameter :: pi=3.1415927
 real(RLK),parameter :: eps=0.62197
 real(RLK),parameter :: rhow=1.0
 real(RLK),parameter :: Mw=18.0
 real(RLK),parameter :: Ma=28.97
 real(RLK),parameter :: gc=8.3144622E07  !UNIVERSAL GAS CONSTANT (G*CM2/S2)/(MOL*K)
 real(RLK),parameter :: gcw=gc/Mw   !GAS CONSTANT FOR LIQUID WATER (G*CM2/S2)/(G*K)
 real(RLK),parameter :: gcd=gc/28.97    !GAS CONSTANT FOR DRY AIR
 real(RLK),parameter :: gca=gc/28.97    !GAS CONSTANT FOR AIR
 real(RLK),parameter :: cpa=1004.8E04   !SPECIFIC HEAT FOR DRY AIR IN ERG G-1 K-1
 real(RLK),parameter :: cpv=1952.0E04   !SPECIFIC HEAT FOR WATER VAPOR IN ERG G-1 K-1
 real(RLK),parameter :: grav=980.0      !CM/S2
 real(RLK),parameter :: mfpac=0.02438   !AIR MEAN FREE PATH COEFFICIENT
 real(RLK),parameter :: alpham=1.0    !MASS ACCOMODATION COEFFICIENT
 real(RLK),parameter :: alphat=1.0     !THERMAL ACCOMODATION COEFFICIENT

END MODULE CONSTANT
