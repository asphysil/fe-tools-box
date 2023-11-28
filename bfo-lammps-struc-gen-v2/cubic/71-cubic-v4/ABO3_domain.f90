!!!!!!!!!!!!!!!!!!!Domain structure aa1/aa2-aa1/aa2!!!!!!!!!!
SUBROUTINE domain_struc_A(natms, b1)
USE my_constants, ONLY : dp_real, dp_int
IMPLICIT NONE

INTEGER(dp_int), INTENT(IN) :: natms
REAL(dp_real), INTENT(INOUT) :: b1(4,natms)
REAL(dp_real) :: b2(3)
REAL(dp_real) :: disp
INTEGER(dp_int) :: i
! displacement
disp= 0.52/SQRT(3.0)

IF (disp <0.001)THEN
        PRINT*, "****************************************"
        PRINT*, "* NO Fe atoms displace***"
        PRINT*, "***Please changed domain code if you need Fe displacement ****"
        PRINT*, "****************************************"
ENDIF

DO i =1, natms
 b2=(/b1(2,i) + disp, b1(3,i)+disp, b1(i,4)+disp/)
 b1(i,2:4)=b2(1:3)
 ENDDO

END SUBROUTINE domain_struc_A

SUBROUTINE domain_struc_B(natms, b1)
USE my_constants, ONLY : dp_real, dp_int
IMPLICIT NONE
INTEGER(dp_int), INTENT(IN) :: natms
REAL(dp_real), INTENT(INOUT) :: b1(4,natms)

REAL(dp_real) :: b2(3)
REAL(dp_real) :: disp

INTEGER(dp_int) :: i

! displacement 
disp= 0.2/SQRT(3.0)
IF (disp <0.001)THEN
        PRINT*, "****************************************"
        PRINT*, "* NO Fe atoms displace***"
        PRINT*, "***Please changed domain code if you need Fe displacement ****"
        PRINT*, "****************************************"
ENDIF
!!!!!!!!!Single domain!!!!!!!!!!!!!!!!
DO i = 1, natms
   b2=(/b1(2,i)+disp, b1(3,i)+disp, b1(i,4)+disp/)
   b1(i,2:4) = b2(1:3)
ENDDO

END SUBROUTINE domain_struc_B


SUBROUTINE domain_struc_O1(natms, b1)
USE my_constants, ONLY : dp_real, dp_int
IMPLICIT NONE

INTEGER(dp_int), INTENT(IN) :: natms
REAL(dp_real), INTENT(INOUT) :: b1(4,natms)

REAL(dp_real) :: b2(3)
REAL(dp_real) :: disp

INTEGER(dp_int) :: i

! displacement 
disp=0.00
!!!!!!!!!Single domain!!!!!!!!!!!!!!!!
DO i = 1, natms
   b2=(/b1(2,i)+disp, b1(3,i)+disp, b1(i,4)+disp/)
   b1(i,2:4) = b2(1:3)
ENDDO
END SUBROUTINE domain_struc_O1

SUBROUTINE domain_struc_O2(natms, b1)
USE my_constants, ONLY : dp_real, dp_int
IMPLICIT NONE

INTEGER(dp_int), INTENT(IN) :: natms
REAL(dp_real), INTENT(INOUT) :: b1(4,natms)

REAL(dp_real) :: b2(3)
REAL(dp_real) :: disp

INTEGER(dp_int) :: i

! displacement
disp=0.00
!!!!!!!!!Single domain!!!!!!!!!!!!!!!!
DO i = 1, natms
   b2=(/b1(2,i)+disp, b1(3,i)+disp, b1(i,4)+disp/)
   b1(i,2:4) = b2(1:3)
ENDDO

END SUBROUTINE domain_struc_O2
