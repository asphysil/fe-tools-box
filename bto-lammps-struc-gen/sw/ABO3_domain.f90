!!!!!!!!!!!!!!!!!!!Domain structure aa1/aa2-aa1/aa2!!!!!!!!!!
SUBROUTINE domain_struc_A(l, b1, b2)
USE dtype
IMPLICIT NONE
REAL(dp_real), INTENT(IN) :: b1(3)
INTEGER(dp_int), INTENT(IN) :: l
REAL(dp_real), INTENT(OUT) :: b2(3)
REAL(dp_real) :: disp

!xxxxxxxxxxxxx 90 degree domain xxxxxxxxxxxx
!!!!!!!! (Px, Py) --> (Px -Py) --> (Px Py) --> (Px -Py)!!!!!
!Pb atom displacement !!!!!
! displacement along 110
disp=0.1
!------120x10x10 cell along x-axis------!
! IF (l <=29) THEN
!    b2=(/b1(1)+disp, b1(2)+disp, b1(3)/)
!  ELSEIF(l>29 .AND.l<=59) THEN
!    b2=(/b1(1)+disp, b1(2)-disp, b1(3)/)
!  ELSEIF(l>59 .AND. l<=89) THEN
!    b2=(/b1(1)+disp, b1(2)+disp, b1(3)/)
!  ELSE
!   b2=(/b1(1)+disp, b1(2)-disp, b1(3)/)
!   ENDIF
!********60x10x10 cell*******
! IF (l <=29) THEN
!    b2=(/b1(1)+disp, b1(2)+disp, b1(3)/)
!  ELSE
!   b2=(/b1(1)+disp, b1(2)-disp, b1(3)/)
!   ENDIF
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!!!!!!!!!Single domain!!!!!!!!!!!!!!!!
 b2=(/b1(1)+disp, b1(2)+disp, b1(3)/)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!xxxxxxxxxxxxx 90 degree domain xxxxxxxxxxxx
!!!!!!!!!aa1/aa2!!!!!!!!!!!!!!!!!!!!!!!!
!!------120x10x10 cell along x-axis------!
! IF (l <=60) THEN
!    b2=(/b1(1)+disp, b1(2)+disp, b1(3)/)
!  ELSE
!   b2=(/b1(1)+disp, b1(2)-disp, b1(3)/)
!   ENDIF
!---------10x10x10-----180 degree domain------
! IF (l <=19) THEN
!    b2=(/b1(1)+disp, b1(2), b1(3)+disp/)
!  ELSE
!   b2=(/b1(1)+disp, b1(2), b1(3)-disp/)
!   ENDIF
END SUBROUTINE domain_struc_A

SUBROUTINE domain_struc_B(l, b1, b2)
USE dtype
IMPLICIT NONE
REAL(dp_real), INTENT(IN) :: b1(3)
INTEGER(dp_int), INTENT(IN) :: l
REAL(dp_real), INTENT(OUT) :: b2(3)
REAL(dp_real) :: disp

!xxxxxxxxxx 90 degree domain xxxxxxxxxxxxx
!!!!!!!! (Px, Py) --> (Px -Py) --> (Px Py) --> (Px -Py)!!!!!
!! Ti  atom displacement !!!!!
! displacement along 110
!------120x10x10 cell along x-axis------!
disp=0.2
! IF (l <=29) THEN
!    b2=(/b1(1)+disp, b1(2)+disp, b1(3)/)
!  ELSEIF(l>29 .AND. l<=59) THEN
!    b2=(/b1(1)+disp, b1(2)-disp, b1(3)/)
!  ELSEIF(l>59 .AND. l<=89) THEN
!    b2=(/b1(1)+disp, b1(2)+disp, b1(3)/)
!  ELSE
!   b2=(/b1(1)+disp, b1(2)-disp, b1(3)/)
!   ENDIF
!********60x10x10 cell*******
! IF (l <=29) THEN
!    b2=(/b1(1)+disp, b1(2)+disp, b1(3)/)
!  ELSE
!   b2=(/b1(1)+disp, b1(2)-disp, b1(3)/)
!   ENDIF
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!!!!!!!!!Single domain!!!!!!!!!!!!!!!!
b2=(/b1(1)+disp, b1(2)+disp, b1(3)/)
!xxxxxxxxxx 90 degree domain xxxxxxxxxxxxx
!!!!!!!!!!!!!!!!!!!!aa1/aa2
! IF (l <=60) THEN
!    b2=(/b1(1)+disp, b1(2)+disp, b1(3)/)
!  ELSE
!   b2=(/b1(1)+disp, b1(2)-disp, b1(3)/)
!   ENDIF
!!!!!!!180 degree domain
! IF (l <=19) THEN
!    b2=(/b1(1), b1(2), b1(3)+disp/)
!  ELSE
!   b2=(/b1(1), b1(2), b1(3)-disp/)
!   ENDIF

END SUBROUTINE domain_struc_B


SUBROUTINE domain_struc_O(l, b1, b2)
USE dtype
IMPLICIT NONE
REAL(dp_real), INTENT(IN) :: b1(3)
INTEGER(dp_int), INTENT(IN) :: l
REAL(dp_real), INTENT(OUT) :: b2(3)
REAL(dp_real) :: disp

!!!!!!!!!!!!!!!!!!!90 degree domain!!!!!!!!!!
!!!!!!!! (Px, Py) --> (Px -Py) --> (Px Py) --> (Px -Py)!!!!!
!! O atom displacement !!!!!
! displacement along 110
!------120x10x10 cell along x-axis------!
disp=0.05
! IF (l <=29) THEN
!    b2=(/b1(1)-disp, b1(2)-disp, b1(3)/)
!  ELSEIF(l>29 .AND. l<=59) THEN
!    b2=(/b1(1)-disp, b1(2)+disp, b1(3)/)
!  ELSEIF(l>59 .AND. l<=89) THEN
!    b2=(/b1(1)-disp, b1(2)-disp, b1(3)/)
!  ELSE
!   b2=(/b1(1)-disp, b1(2)+disp, b1(3)/)
!   ENDIF
!
!********60x10x10*******
! IF (l <=29) THEN
!    b2=(/b1(1)-disp, b1(2)-disp, b1(3)/)
!  ELSE
!   b2=(/b1(1)-disp, b1(2)+disp, b1(3)/)
!   ENDIF
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!!!!!!!!!Single domain!!!!!!!!!!!!!!!!
b2=(/b1(1)+disp, b1(2)+disp, b1(3)/)
!
!xxxxxxxxxx 90 degree domain xxxxxxxxxxxxx
!!!!!!!!!!!!!!!!!!!!aa1/aa2
! IF (l <=60) THEN
!    b2=(/b1(1)+disp, b1(2)-disp, b1(3)/)
!  ELSE
!   b2=(/b1(1)+disp, b1(2)+disp, b1(3)/)
!   ENDIF
!!!!!!180 degree domain
! IF (l <=19) THEN
!    b2=(/b1(1)+disp, b1(2), b1(3)-disp/)
!  ELSE
!   b2=(/b1(1)+disp, b1(2), b1(3)-disp/)
!   ENDIF

END SUBROUTINE domain_struc_O
