MODULE constant
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp=SELECTED_REAL_KIND(15) !, int64=SELECTED_INT_KIND(15)
    !REAL(dp), PARAMETER :: zqa=4.37d0, zqb=3.49d0, zqo1=-2.61d0, zqo2=-2.61d0, zqo3=-2.61d0 ! BFO
     REAL(dp), PARAMETER :: zqpb=3.94359d0, zqw=6.64710d0, zqmg=2.83758d0, zqo1=-1.0d0, zqo2=-1.0d0, zqo3=-1.0d0
    ! qo1= parallel to T-O bond (  Born Effective value large in parallel direction due to bond elongation and compression motion)
    ! qo2 and qo3 perpendicular to T-O bond (  Born Effective value small in perpendicular direction due to bond bending motion)
    INTEGER, PARAMETER :: fileID_in=13
    INTEGER, PARAMETER :: fileID_in1=14, fileID_in2=15, fileID_in3=16, fileID_in4=17
    INTEGER, PARAMETER :: fileID_ovito=41
    INTEGER, PARAMETER :: fileID_totpol=20, fileID_locpol=21
    INTEGER, PARAMETER :: fileID_pol_ovito=22
    INTEGER, PARAMETER :: fileID_rtheta_plt=28
    INTEGER, PARAMETER :: fileID_afe=29
    INTEGER, PARAMETER :: fileID_disp_a=30, fileID_disp_b=31
    INTEGER, PARAMETER :: fileID_disp_a_sublatt=32
    INTEGER, PARAMETER :: fileID_afe_in=33
END MODULE constant
