SUBROUTINE crys_struc_ortho(nu, nba, nti, no1, latt_vec, atms_ortho)
  USE my_constants, ONLY : dp_real
  IMPLICIT NONE 

  INTEGER, INTENT(IN) :: nu
  INTEGER,INTENT(OUT) :: nba, nti,no1
  REAL(dp_real), INTENT(OUT) :: atms_ortho(3,nu), latt_vec(3,3)
  

  INTEGER :: i,j 
  
  DO i = 1, 3
    DO j = 1, 3
    latt_vec(i,j) = 0
    ENDDO 
  ENDDO

  latt_vec(1,1) = 4.0*SQRT(2.0)
  latt_vec(2,2) = 4.0*SQRT(2.0)
  latt_vec(3,3) = 4.0

  nba =2
  nti =2
  no1 =2
   
  atms_ortho(:,1) =(/0.00000000,   0.00000000,   0.00000000/)
  atms_ortho(:,2) = (/0.50000000,   0.50000000,   0.00000000/)

  atms_ortho(:,3) = (/0.00000000,   0.50000000,   0.50000000/)
  atms_ortho(:,4) = (/0.50000000,   0.00000000,   0.50000000/)

  atms_ortho(:,5) = (/0.00000000,   0.50000000,   0.00000000/)
  atms_ortho(:,6) = (/0.50000000,   0.00000000,   0.00000000/)
        

  atms_ortho(:,7) = (/0.25000000,   0.25000000,   0.50000000/)
  atms_ortho(:,8) = (/0.75000000,   0.75000000,   0.50000000/)

  atms_ortho(:,9) = (/0.25000000,   0.75000000,   0.50000000/)
  atms_ortho(:,10) = (/0.75000000,   0.25000000,   0.50000000/)

  END SUBROUTINE crys_struc_ortho

  SUBROUTINE crys_struc_cubic(nu, nba, nti, no1, latt_vec, atms_cub)
    USE my_constants, ONLY : dp_real
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nu
    REAL(dp_real), INTENT(OUT) :: atms_cub(3, nu), latt_vec(3,3)

    INTEGER, INTENT(OUT) :: nba, nti, no1

    INTEGER :: i, j

    nba = 1  ! No. A-site atom
    nti = 1  ! No. B-site atom
    no1 = 1  ! No. O1 atoms

    DO i = 1, nu 
      DO j = 1, 3
      atms_cub(j, i) =0.0
      ENDDO 
    ENDDO
                 
    DO i=1,3
      DO j =1,3 
        latt_vec(i,j) =0.0
      ENDDO
    ENDDO

    latt_vec(1,1) = 4.0
    latt_vec(2,2) = 4.0
    latt_vec(3,3) = 4.0

    atms_cub(:,1) =(/0.00000000,   0.00000000,   0.00000000/)
  
    atms_cub(:,2) = (/0.50000000,   0.50000000,   0.50000000/)
    
    atms_cub(:,3) = (/0.50000000,   0.50000000,   0.00000000/) ! xy-plane
    atms_cub(:,4) = (/0.5000000,   0.00000000,   0.50000000/)  ! xz-plane
    atms_cub(:,5) = (/0.0000000,   0.50000000,   0.50000000/)  ! yz-plane
    

  
  END SUBROUTINE crys_struc_cubic



  SUBROUTINE crys_struc_cubic_r3c(nu, nba, nti, no1, latt_vec, atms_cub)
    USE my_constants, ONLY : dp_real
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nu
    REAL(dp_real), INTENT(OUT) :: atms_cub(3,nu), latt_vec(3,3)

    INTEGER, INTENT(OUT) :: nba, nti, no1

    INTEGER :: i, j

    nba = 8  ! No. A-site atom
    nti = 8  ! No. B-site atom
    no1 = 8  ! No. O1 atoms

    DO i = 1, nu
      DO j = 1, 3
      atms_cub(j, i) =0.0
      ENDDO
    ENDDO

    DO i=1,3
      DO j =1,3
        latt_vec(i,j) =0.0
      ENDDO
    ENDDO

OPEN(UNIT=41, FILE='struc-cubic-r3c.vasp', action='read')

READ(41,*) 
READ(41,*)

DO i=1, 3
READ(41,*) (latt_vec(i,j), j=1,3)
ENDDO

READ(41,*) 
READ(41,*)
READ(41,*) 

DO i=1, 40
READ(41,*) atms_cub(1,i), atms_cub(2,i), atms_cub(3,i)
ENDDO

CLOSE(41)

  END SUBROUTINE crys_struc_cubic_r3c


 SUBROUTINE crys_struc_cubic_r3c_opt(nu, nba, nti, no1, latt_vec, atms_cub)
  USE my_constants, ONLY : dp_real
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nu
    REAL(dp_real), INTENT(OUT) :: atms_cub(3,nu), latt_vec(3,3)

    INTEGER, INTENT(OUT) :: nba, nti, no1

    INTEGER :: i, j

    nba = 8  ! No. A-site atom
    nti = 8  ! No. B-site atom
    no1 = 8  ! No. O1 atoms

    DO i = 1, nu
      DO j = 1, 3
      atms_cub(j,i) =0.0
      ENDDO
    ENDDO

    DO i=1,3
      DO j =1,3
        latt_vec(i,j) =0.0
      ENDDO
    ENDDO

OPEN(UNIT=41, FILE='struc-cubic-r3c-opt.vasp', action='read')

READ(41,*)
READ(41,*)

DO i=1, 3
READ(41,*) (latt_vec(i,j), j=1,3)
ENDDO

READ(41,*)
READ(41,*)
READ(41,*)

DO i=1, 40
READ(41,*) atms_cub(1,i), atms_cub(2,i), atms_cub(3,i)
ENDDO

CLOSE(41)

  END SUBROUTINE crys_struc_cubic_r3c_opt
