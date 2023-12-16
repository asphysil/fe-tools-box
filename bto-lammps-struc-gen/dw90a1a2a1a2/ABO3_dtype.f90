
 MODULE dtype
 IMPLICIT NONE
 INTEGER, PARAMETER :: dp_real=SELECTED_REAL_KIND(16), dp_int=SELECTED_INT_KIND(10)
 REAL(dp_real), PARAMETER :: charge_A=1.34730, charge_B=1.28905, charge_O=-0.87878 ! A=Ba, B=Ti, O=O BaTiO3 compound
! REAL(dp_real), PARAMETER :: charge_A=1.38177, charge_B=0.99997, charge_O=-0.79391 ! A=Pb, B=Ti, O=O PbTiO3 compound
 INTEGER(dp_int), PARAMETER :: nx=60,ny=15,nz=15 ! cell dimension
 REAL, PARAMETER :: massA=137.327, massB=47.88,massO=15.9994 !massA=Ba, massB=Ti, massO=O
 !REAL, PARAMETER :: massA=207.2, massB=47.88,massO=15.9994 ! massA=Pb, massB=Ti, massO=O
 INTEGER, PARAMETER :: angle_coeffs=200  !BaTiO3
 !INTEGER, PARAMETER :: angle_coeffs=50  !PaTiO3
 REAL(dp_real), PARAMETER :: strainx=1.5,strainy=1.5, strainz=0.0    !strain in percentage
 REAL(dp_real):: a11=4.000, a22=4.000, a33=4.109  ! a=b=c=4.005 BaTiO3 latice parameters
 !REAL(dp_real):: a11=3.86, a22=3.86, a33=4.344  ! a=b=c=3.86 PbTiO3 latice parameters 
END MODULE dtype
