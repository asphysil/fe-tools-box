!!!!!!!!!!!!!!!!!!!Domain structure aa1/aa2-aa1/aa2!!!!!!!!!!
SUBROUTINE domain_struc_A(natms, b1)
USE my_constants, ONLY : dp_real, dp_int, nx
IMPLICIT NONE

INTEGER(dp_int), INTENT(IN) :: natms
REAL(dp_real), INTENT(INOUT) :: b1(4,natms)
REAL(dp_real) :: b2(3)
REAL(dp_real) :: d1, d2
INTEGER(dp_int) :: i, j, l

! d1 displacement along cartesian direction
! d2 displacement along 110
d1=0.55/SQRT(3.0)
d2= d1*SQRT(2.0)

IF ( MOD(nx,2) /=0) THEN
   PRINT*, "Dipole is not 50% up and 50% down"
   PRINT*,"   "
   PRINT*, " Please change nx value to even number"
ENDIF
l=nx/2
!------cell along x-axis 71DW------!
DO i = 1, natms
   j = INT(b1(1,i))

 IF (j <l) THEN
    b2=(/b1(2,i), b1(3,i)+d2, b1(4,i)+d1/)
  ELSE
   b2=(/b1(2,i), b1(3,i)+d2, b1(4,i)-d1/)
  ENDIF
 b1(2:4,i)=b2(1:3)
 ENDDO
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!!!!!!!!!Single domain!!!!!!!!!!!!!!!!
!DO i =1, natms
! b2=(/b1(i,2), b1(i,3)+(disp*SQRT(2.0)), b1(i,4)+disp/)
! b1(2:4,i)=b2(1:3)
! ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE domain_struc_A

SUBROUTINE domain_struc_B(natms, b1)
USE my_constants, ONLY : dp_real, dp_int, nx
IMPLICIT NONE
INTEGER(dp_int), INTENT(IN) :: natms
REAL(dp_real), INTENT(INOUT) :: b1(4,natms)

REAL(dp_real) :: b2(3)
REAL(dp_real) :: d1, d2

INTEGER(dp_int) :: i, j, l

d1=0.22/SQRT(3.0)
d2 = d1*SQRT(2.0)

IF ( MOD(nx,2) /=0) THEN
   PRINT*, "Dipole is not 50% up and 50% down"
   PRINT*,"   "
   PRINT*, " Please change nx value to even number"
ENDIF
l=nx/2
!------cell along x-axis 180DW------!
DO i = 1, natms
   j = INT(b1(1, i))

 IF (j <l) THEN
    b2=(/b1(2,i), b1(3,i)+d2, b1(4,i)+d1/)
  ELSE
   b2=(/b1(2,i), b1(3,i)+d2, b1(4,i)-d1/)
   ENDIF

 b1(2:4,i)=b2(1:3)
 ENDDO

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!!!!!!!!!Single domain!!!!!!!!!!!!!!!!
!DO i = 1, natms
!   b2=(/b1(i,2), b1(i,3)+(disp*SQRT(2.0)), b1(i,4)+disp/)
!   b1(2:4,i) = b2(1:3)
!ENDDO

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
     b2=(/b1(2,i), b1(3,i)+disp, b1(4,i)+disp/)
     b1(2:4,i) = b2(1:3)
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
        b2=(/b1(2,i), b1(3,i)+disp, b1(4,i)+disp/)
        b1(2:4,i) = b2(1:3)
  ENDDO

END SUBROUTINE domain_struc_O2
