SUBROUTINE crys_struc_ortho(nu, nba, nti, no1, latt_vec, atms_ortho)
  USE my_constants, ONLY : dp_real,a_opt_bfo
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

  
  latt_vec(1,1) = a_opt_bfo*SQRT(2.0)
  latt_vec(2,2) = a_opt_bfo*SQRT(2.0)
  latt_vec(3,3) = a_opt_bfo

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
  
  SUBROUTINE crys_struc_ortho20(nu, nba, nti, no1, latt_vec, atms_ortho)
    USE my_constants, ONLY : dp_real,a_opt_bfo
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
  
    
    latt_vec(1,1) = a_opt_bfo*SQRT(2.0)
    latt_vec(2,2) = a_opt_bfo*SQRT(2.0)
    latt_vec(3,3) = a_opt_bfo*2.0
  
    nba =4
    nti =4
    no1 =4

  atms_ortho(:,1) = (/0.000000000,         0.000000000,         0.000000000/)
  atms_ortho(:,2) = (/0.000000000,         0.000000000,         0.500000000/)
  atms_ortho(:,3) = (/0.500000000,         0.500000000,         0.000000000/)
  atms_ortho(:,4) = (/0.500000000,         0.500000000,         0.500000000/)

  atms_ortho(:,5) = (/0.000000000,         0.500000000,         0.250000000/)
  atms_ortho(:,6) = (/0.000000000,         0.500000000,         0.750000060/)
  atms_ortho(:,7) = (/0.500000000,         0.000000000,         0.250000000/)
  atms_ortho(:,8) = (/0.500000000,         0.000000000,         0.750000060/)

  atms_ortho(:,9) = (/0.000000000,         0.500000000,         0.000000000/)
  atms_ortho(:,10) = (/0.000000000,         0.500000000,         0.500000000/)
  atms_ortho(:,11) = (/0.500000000,         0.000000000,         0.000000000/)
  atms_ortho(:,12) = (/0.500000000,         0.000000000,         0.500000000/)

  atms_ortho(:,13) = (/0.250000000,         0.250000000,         0.250000000/)
  atms_ortho(:,14) = (/0.250000000,         0.250000000,         0.750000060/)
  atms_ortho(:,15) = (/0.750000000,         0.750000000,         0.250000000/)
  atms_ortho(:,16) = (/0.750000000,         0.750000000,         0.750000060/)

  atms_ortho(:,17) = (/0.250000000,         0.750000000,         0.250000000/)
  atms_ortho(:,18) = (/0.250000000,         0.750000000,         0.750000060/)
  atms_ortho(:,19) = (/0.750000000,         0.250000000,         0.250000000/)
  atms_ortho(:,20) = (/0.750000000,         0.250000000,         0.750000060/)
END SUBROUTINE crys_struc_ortho20

  SUBROUTINE crys_struc_ortho20_opt(nu, nba, nti, no1, latt_vec, atms_ortho)
    USE my_constants, ONLY : dp_real,a_opt_bfo
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
  
    
    latt_vec(1,1) = a_opt_bfo*SQRT(2.0)
    latt_vec(2,2) = a_opt_bfo*SQRT(2.0)
    latt_vec(3,3) = a_opt_bfo*2.0
  
    nba =4
    nti =4
    no1 =4

  atms_ortho(:,1)  = (/ 0.0000000755,        0.0873506407,        0.0439781769/)
  atms_ortho(:,2)  = (/-0.0000002558,        0.0873505982,        0.5439781895/)
  atms_ortho(:,3)  = (/ 0.4999998722,        0.5873506888,        0.0439781758/)
  atms_ortho(:,4)  = (/ 0.5000001690,        0.5873507650,        0.5439779726/)

  atms_ortho(:,5)  = (/ 0.0000000749,        0.5316198099,        0.2660904523/)
  atms_ortho(:,6)  = (/ 0.0000001204,        0.5316197301,        0.7660902615/)
  atms_ortho(:,7)  = (/ 0.4999999296,        0.0316199683,        0.2660903934/)
  atms_ortho(:,8)  = (/ 0.4999998301,        0.0316197202,        0.7660901824/)

  atms_ortho(:,9)  = (/ 0.0000002714,        0.4846845267,        0.0067046189/)
  atms_ortho(:,10) = (/-0.0000000694,        0.4846844949,        0.5067047042/)
  atms_ortho(:,11) = (/ 0.4999997006,       -0.0153155530,        0.0067045669/)
  atms_ortho(:,12) = (/ 0.5000001274,       -0.0153154254,        0.5067046583/)

  atms_ortho(:,13) = (/ 0.2359380772,        0.2488496857,        0.2425771494/)
  atms_ortho(:,14) = (/ 0.2359380104,        0.2488494838,        0.7425769940/)
  atms_ortho(:,15) = (/ 0.7359380064,        0.7488495949,        0.2425770137/)
  atms_ortho(:,16) = (/ 0.7359380580,        0.7488494933,        0.7425769170/)

  atms_ortho(:,17) = (/ 0.2640619518,        0.7488496538,        0.2425771696/)
  atms_ortho(:,18) = (/ 0.2640619214,        0.7488494823,        0.7425769243/)
  atms_ortho(:,19) = (/ 0.7640619952,        0.2488496419,        0.2425770288/)
  atms_ortho(:,20) = (/ 0.7640619829,        0.2488494941,        0.7425770406/)
END SUBROUTINE crys_struc_ortho20_opt

  SUBROUTINE crys_struc_cubic(nu, nba, nti, no1, latt_vec, atms_cub)
    USE my_constants, ONLY : dp_real,a_opt_bfo
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

    latt_vec(1,1) = a_opt_bfo
    latt_vec(2,2) = a_opt_bfo
    latt_vec(3,3) = a_opt_bfo

    atms_cub(:,1) =(/0.00000000,   0.00000000,   0.00000000/)

    atms_cub(:,2) = (/0.50000000,   0.50000000,   0.50000000/)

    atms_cub(:,3) = (/0.50000000,   0.50000000,   0.00000000/) ! xy-plane
    atms_cub(:,4) = (/0.5000000,   0.00000000,   0.50000000/)  ! xz-plane
    atms_cub(:,5) = (/0.0000000,   0.50000000,   0.50000000/)  ! yz-plane



  END SUBROUTINE crys_struc_cubic
