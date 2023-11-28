!!!!!!!!!!!!!!!!!!!Domain structure aa1/aa2-aa1/aa2!!!!!!!!!!
SUBROUTINE domain_struc_A(natms, b1)
USE my_constants, ONLY : dp_real
IMPLICIT NONE

INTEGER, INTENT(IN) :: natms
REAL(dp_real), INTENT(INOUT) :: b1(4,natms)
REAL(dp_real) :: b2(3)
REAL(dp_real) :: d1, d2
INTEGER :: i
! displacement
! d1 displacement along cartesian direction
! d2 displacement along 110
d1=0.0
!d1=0.5/SQRT(3.0)
d2= d1*SQRT(2.0)



IF (d1 <0.001)THEN
        PRINT*, "****************************************"
        PRINT*, "* NO Fe atoms displace***"
        PRINT*, "***Please changed domain code if you need Fe displacement ****"
        PRINT*, "****************************************"
ENDIF

DO i =1, natms
 b2=(/b1(2,i), b1(3,i)+d2, b1(4,i)+d1/)
 b1(2:4,i)=b2(1:3)
ENDDO

END SUBROUTINE domain_struc_A

SUBROUTINE domain_struc_B(natms, b1)
   USE my_constants, ONLY : dp_real
IMPLICIT NONE
INTEGER, INTENT(IN) :: natms
REAL(dp_real), INTENT(INOUT) :: b1(4,natms)

REAL(dp_real) :: b2(3)
REAL(dp_real) :: d1, d2

INTEGER :: i

! displacement 
d1=0.0
!d1=0.2/SQRT(3.0)
d2 = d1*SQRT(2.0)

IF (d1 <0.001)THEN
        PRINT*, "****************************************"
        PRINT*, "* NO Fe atoms displace***"
        PRINT*, "***Please changed domain code if you need Fe displacement ****"
        PRINT*, "****************************************"
ENDIF
!!!!!!!!!Single domain!!!!!!!!!!!!!!!!
DO i = 1, natms
 b2=(/b1(2,i), b1(3,i)+d2, b1(4,i)+d1/)
   b1(2:4, i) = b2(1:3)
ENDDO

END SUBROUTINE domain_struc_B


SUBROUTINE domain_struc_O1(natms, b1)
   USE my_constants, ONLY : dp_real
IMPLICIT NONE

INTEGER, INTENT(IN) :: natms
REAL(dp_real), INTENT(INOUT) :: b1(4, natms)

REAL(dp_real) :: b2(3)
REAL(dp_real) :: disp

INTEGER :: i

! displacement 
disp=0.00
!!!!!!!!!Single domain!!!!!!!!!!!!!!!!
DO i = 1, natms
    b2=(/b1(2,i) + disp, b1(3,i)+disp, b1(4,i)+disp/)
   b1(2:4,i) = b2(1:3)
ENDDO
END SUBROUTINE domain_struc_O1



SUBROUTINE domain_struc_O2(natms, b1)
   USE my_constants, ONLY : dp_real
IMPLICIT NONE

INTEGER, INTENT(IN) :: natms
REAL(dp_real), INTENT(INOUT) :: b1(4, natms)

REAL(dp_real) :: b2(3)
REAL(dp_real) :: disp

INTEGER :: i

! displacement 
disp=0.00
!!!!!!!!!Single domain!!!!!!!!!!!!!!!!
DO i = 1, natms
    b2=(/b1(2,i) + disp, b1(3,i)+disp, b1(4,i)+disp/)
   b1(2:4,i) = b2(1:3)
ENDDO
END SUBROUTINE domain_struc_O2



SUBROUTINE domain_struc_O3(natms, b1)
   USE my_constants, ONLY : dp_real
IMPLICIT NONE

INTEGER, INTENT(IN) :: natms
REAL(dp_real), INTENT(INOUT) :: b1(4, natms)

REAL(dp_real) :: b2(3)
REAL(dp_real) :: disp

INTEGER :: i

! displacement 
disp=0.00
!!!!!!!!!Single domain!!!!!!!!!!!!!!!!
DO i = 1, natms
    b2=(/b1(2,i) + disp, b1(3,i)+disp, b1(4,i)+disp/)
   b1(2:4,i) = b2(1:3)
ENDDO
END SUBROUTINE domain_struc_O3


