!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE ANGLE_AMONG1(O1, O2, O3, alpha)
 USE dtype
 IMPLICIT NONE
! INTEGER, PARAMETER :: dp_real=SELECTED_REAL_KIND(16)
 REAL(dp_real), INTENT(IN) :: O1(1,3), O2(1,3), O3(1,3)
 REAL(dp_real), INTENT(OUT)::alpha
 REAl(dp_real):: l1(3), l2(3) !, degree
 !INTEGER :: i
 !degree=180.0/(4.0*ATAN(1.0))
 l1(1)=O3(1,1)-O2(1,1);l1(2)=O3(1,2)-O2(1,2);l1(3)=O3(1,3)-O2(1,3)
 l2(1)=O1(1,1)-O2(1,1);l2(2)=O1(1,2)-O2(1,2);l2(3)=O1(1,3)-O2(1,3)
! PRINT*, (O1(1,i), i=1,3)
! PRINT*, (O2(1, i),i=1,3)
! PRINT*, (O3(1, i),i=1,3)
 alpha=DOT_PRODUCT(l1,l2)/(SQRT(DOT_PRODUCT(l1,l1))* SQRT(DOT_PRODUCT(l2,l2)))
 END SUBROUTINE ANGLE_AMONG1

 SUBROUTINE ANGLE_AMONG2(O1, O2, O3, alpha)
 USE dtype
 IMPLICIT NONE
 !INTEGER, PARAMETER :: dp_real=SELECTED_REAL_KIND(16)
 REAL(dp_real), INTENT(IN) :: O1(1,3), O2(1,3), O3(1,3)
 REAL(dp_real), INTENT(OUT)::alpha
 REAl(dp_real):: l1(3), l2(3) !, degree
! INTEGER :: i
 !degree=180.0/(4.0*ATAN(1.0))
 l1(1)=O3(1,1)-O2(1,1);l1(2)=O3(1,2)-O2(1,2);l1(3)=O3(1,3)-O2(1,3)
 l2(1)=O1(1,1)-O2(1,1);l2(2)=O1(1,2)-O2(1,2);l2(3)=O1(1,3)-O2(1,3)
! PRINT*, (O1(1,i), i=1,3)
! PRINT*, (O2(1, i),i=1,3)
! PRINT*, (O3(1, i),i=1,3)
 alpha=-1*DOT_PRODUCT(l1,l2)/(SQRT(DOT_PRODUCT(l1,l1))* SQRT(DOT_PRODUCT(l2,l2)))
 END SUBROUTINE ANGLE_AMONG2
