MODULE writefiles
    use constant, only : dp, PI 
    use data_structure, only : latt_vec, a_atms, b_atms, o_atms,&
                              a_disp, b_disp,&
                              a_pol, b_pol, pol_unit,&
                              sp, tndump,&
                              no_bb_neigh, bb_neigh
    IMPLICIT NONE 
 
 private    
 public ::  write_xsf, write_pol,& 
           write_dump_ovito, write_pol_ovito,&
           write_dump_ovito_LMP
contains 
SUBROUTINE write_xsf(fileID, filename,  nb, ntot)
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: fileID,  nb, ntot 
    CHARACTER(LEN=40), INTENT(IN) :: filename
    INTEGER :: no 
    REAL :: f(3)  
    INTEGER :: i, j 
    INTEGER :: atmnum(3)

    no = 3*nb 
    f = (/0.0, 0.0, 0.0/)
    atmnum(1) = 83 ! Bi
    atmnum(2) = 26 ! Fe
    atmnum(3) = 8  ! O
OPEN(UNIT=fileID, FILE=filename, ACTION="WRITE")
!!!!
!!! Polarization calculation
!!!!!!WRITEING for xcruden .xsf file FORMAT
WRITE(fileID,*) '#A B O'
WRITE(fileID,*)'CRYSTAL'
WRITE(fileID,*)'PRIMVEC'
DO i=1,3
  WRITE(fileID,300) latt_vec(i,1), latt_vec(i,2), latt_vec(i,3)
ENDDO
WRITE(fileID,*)'CONVVEC'
DO i=1,3
  WRITE(fileID,300)latt_vec(i,1), latt_vec(i,2), latt_vec(i,3)
ENDDO
!300 FORMAT(1X, 3F13.8)
WRITE(fileID,*)'PRIMCOORD'
WRITE(fileID,*) ntot, '1'

DO i =1, nb 
WRITE(fileID,301) atmnum(1), (a_atms(j,i),j=1,3), (a_pol(j,i), j=1,3)
ENDDO

DO i =1, nb 
WRITE(fileID,301) atmnum(2), (b_atms(j,i),j=1,3), (b_pol(j,i), j=1,3)
ENDDO

DO i=1, no
  WRITE(fileID,301) atmnum(3), (o_atms(j,i), j=1,3), f(1), f(2), f(3)
ENDDO
300 FORMAT(2X, 3F12.7)
301 FORMAT (2X, I5, 6F12.8)

close(fileID)
END SUBROUTINE write_xsf 

SUBROUTINE write_pol(fileID_totpol, fileID_locpol, nb)
    IMPLICIT NONE
    INTEGER, INTENT(IN) ::  fileID_totpol, fileID_locpol 
    INTEGER, INTENT(IN) ::  nb 
    REAL(dp) :: ptot(3)  
    REAL(dp) :: nb_dp
    REAL(dp) :: px, py, pz  
    INTEGER :: i, j 



nb_dp = DBLE(nb)

ptot(1)=SUM(pol_unit(1,1:nb))/nb_dp
ptot(2)=SUM(pol_unit(2,1:nb))/nb_dp
ptot(3)=SUM(pol_unit(3,1:nb))/nb_dp

! PRINT*, " ----------------------Warning-------------"
! PRINT*, " Please change total polarization calculation formula --" 
! PRINT*, 'if the following strutural infomations are  not correct'

! PRINT*, " x, y, and z axis are  along a, b and c crystalographic direction, respectively. Axis are orthogonal"
! !PRINT*, " x-axis along [1-10], y-axis along [110] and z-axis along [001]"
! PRINT*, "-----------------------------------------------"
! !px = ptot(1)/sqrt(2.0) + ptot(2)/sqrt(2.0)
! !py =  - ptot(1)/sqrt(2.0) + ptot(2)/sqrt(2.0)
! !pz = ptot(3)
! !
! PRINT*, " x-axis along [100], y-axis along [100] and z-axis along [001]"
!
px = ptot(1)
py = ptot(2)
pz = ptot(3)
WRITE(fileID_totpol, 303) px, py, pz, SQRT(px**2 + py**2 + pz**2)

DO i=1, nb 
    WRITE(fileID_locpol,302) (b_atms(j,i),j=1,3), (pol_unit(j,i),j=1,3)
ENDDO 
302 FORMAT(2X, 3F12.7, 3X, 3F12.7)
303 FORMAT(2X,  3F12.7, 3X, F12.7)

END SUBROUTINE write_pol

SUBROUTINE write_dump_ovito(fileID, nb)
IMPLICIT NONE
INTEGER, INTENT(IN) ::  fileID
INTEGER, INTENT(IN) ::  nb
INTEGER :: i, j
INTEGER :: id=1 
REAL :: x 
REAL(dp) :: ploc 
x=0.0000

write(fileID,'(A)') 'ITEM: TIMESTEP'
write(fileID,'(I7)') tndump
write(fileID,'(A)') 'ITEM: NUMBER OF ATOMS'
write(fileID, '(I7)') nb 
write(fileID, '(A)') 'ITEM: BOX BOUNDS pp pp pp'
DO i=1,3
    write(fileID,'(2F15.8)') x, latt_vec(i,i)
ENDDO
write(fileID, '(A)') 'ITEM: ATOMS x y z fx fy fz mux muy muz mu' 
DO i = 1, nb
ploc = SQRT(pol_unit(1,i)*pol_unit(1,i) + pol_unit(2,i)*pol_unit(2,i) + pol_unit(3,i)*pol_unit(3,i) )
write(fileID, 400) (b_atms(j,i), j=1,3), (sp(j,i), j=1,3), (pol_unit(j,i),j=1,3), ploc
ENDDO
400 FORMAT(1X, 3F15.10, 3F15.10, 3F15.10, F15.10)
END SUBROUTINE write_dump_ovito

SUBROUTINE write_pol_ovito(fileID, nb)

IMPLICIT NONE
INTEGER, INTENT(IN) ::  fileID
INTEGER, INTENT(IN) ::  nb
INTEGER :: i, j
INTEGER :: id=1
REAL :: x
REAL(dp) :: ploc
x=0.0000

write(fileID,'(A)') 'ITEM: TIMESTEP'
write(fileID,'(I7)') tndump
write(fileID,'(A)') 'ITEM: NUMBER OF ATOMS'
write(fileID, '(I7)') nb
write(fileID, '(A)') 'ITEM: BOX BOUNDS pp pp pp'
DO i=1,3
    write(fileID,'(2F15.8)') x, latt_vec(i,i)
ENDDO
write(fileID, '(A)') 'ITEM: ATOMS x y z mux muy muz mu'
DO i = 1, nb
ploc = SQRT(pol_unit(1,i)*pol_unit(1,i) + pol_unit(2,i)*pol_unit(2,i) + pol_unit(3,i)*pol_unit(3,i) )
write(fileID, 400) (b_atms(j,i), j=1,3), (pol_unit(j,i),j=1,3), ploc
ENDDO
400 FORMAT(1X, 3F15.10, 3F15.10, F15.10)
END SUBROUTINE write_pol_ovito

  SUBROUTINE write_rtheta(fileID, nb)
    IMPLICIT NONE
    INTEGER, INTENT(IN) ::  fileID
    INTEGER, INTENT(IN) ::  nb
    INTEGER :: i, j
    REAL :: x
    REAL(dp) :: ploc
    REAl(dp) :: u0(3), u1(3), u2(3)
    REAL(dp) :: angle1, angle2, norm, inorm 
    DO i = 1, nb
 
     norm = SQRT(a_disp(1,i)**2 + a_disp(2,i)**2 + a_disp(3,i)**2)
     inorm = 1.0/norm 
   
     u0(1) =a_disp(1,i)*inorm
     u0(2) =a_disp(2,i)*inorm
     u0(3) =a_disp(3,i)*inorm 
   
     angle1 = ACOS(DOT_PRODUCT(u0,u1))
   
     angle2 = ACOS(DOT_PRODUCT(u0,u2))     
     write(fileID,*) angle1, angle2, norm
    ENDDO
  END SUBROUTINE write_rtheta
 
 
  
SUBROUTINE write_dump_ovito_LMP(fileID, nb)
  IMPLICIT NONE
  INTEGER, INTENT(IN) ::  fileID
  INTEGER, INTENT(IN) ::  nb
  INTEGER :: i, j, k
  INTEGER ::  m 
  INTEGER :: id=1 
  
  REAL(dp) :: x, nn_bb 
  REAL(dp), DIMENSION(3) :: w1, w2, l_sum, m_sum 
  REAL(dp), DIMENSION(nb) :: ploc, m_vel, m_disp
  REAL(dp), DIMENSION(3,nb):: vel, disp     
  REAL(dp) :: bb_neigh_spin(3, no_bb_neigh, nb)
  x=0.000
  nn_bb = 6.0_dp 
! vel == M ( Ferro magnetic vector)
! disp == L (Neel vector)

  DO i =1 , nb 
    DO j = 1, no_bb_neigh
      m = bb_neigh(j,i)
      bb_neigh_spin(1:3, j, i) = sp(1:3,m)
    ENDDO
  ENDDO


  DO i =1 , nb 
    l_sum=(/0.0,0.0,0.0/)
    m_sum=(/0.0,0.0,0.0/)
    DO j = 1, no_bb_neigh
      w1(1) = sp(1,i) - bb_neigh_spin(1,j,i)
      w1(2) = sp(2,i) - bb_neigh_spin(2,j,i)
      w1(3) = sp(3,i) - bb_neigh_spin(3,j,i)

      w2(1) = sp(1,i) + bb_neigh_spin(1,j,i)
      w2(2) = sp(2,i) + bb_neigh_spin(2,j,i)
      w2(3) = sp(3,i) + bb_neigh_spin(3,j,i)

     l_sum(1) =  l_sum(1) + ( w1(2)*w2(3) - w1(3)*w2(2) )
     l_sum(2) =  l_sum(2) + ( w1(3)*w2(1) - w1(1)*w2(3) )
     l_sum(3) =  l_sum(3) + ( w1(1)*w2(2) - w1(2)*w2(1) )
     ! DO k =1, 3
    !l_sum(k) = l_sum(k) + (sp(k,i) - bb_neigh_spin(k,j,i))
    !m_sum(k) = m_sum(k) + (sp(k,i) + bb_neigh_spin(k,j,i))  
    !  ENDDO
    !  if (i==1) write(100,*) spx, spy, spz 
    ENDDO
    disp(1:3, i) = l_sum(1:3)/nn_bb
    !vel(1:3, i) = m_sum(1:3)/nn_bb
    ploc(i) = SQRT(pol_unit(1,i)**2 + pol_unit(2,i)**2+ pol_unit(3,i)**2)
    m_disp(i) = SQRT(disp(1,i)**2 + disp(2,i)**2 + disp(3,i)**2)
    !m_vel(i) = SQRT(vel(1,i)**2 + vel(2,i)**2 + vel(3,i)**2) 
  ENDDO

  write(fileID,'(A)') 'ITEM: TIMESTEP'
  write(fileID,'(I7)') tndump
  write(fileID,'(A)') 'ITEM: NUMBER OF ATOMS'
  write(fileID, '(I7)') nb 
  write(fileID, '(A)') 'ITEM: BOX BOUNDS pp pp pp'
  DO i=1,3
      write(fileID,'(2F15.8)') x, latt_vec(i,i)
  ENDDO
  write(fileID, '(A)') 'ITEM: ATOMS x y z fx fy fz f  mux muy muz mu' 

  DO i = 1, nb
  write(fileID, 400) (b_atms(j,i), j=1,3),(disp(j,i),j=1,3), m_disp(i), (pol_unit(j,i),j=1,3), ploc(i)
  ! (vel(j,i), j=1,3), m_vel(i), (disp(j,i),j=1,3), m_disp(i), (pol_unit(j,i),j=1,3), ploc(i)
  ENDDO
  400 FORMAT(1X, 3F15.10, 8F15.10)
  END SUBROUTINE write_dump_ovito_LMP
  
END MODULE writefiles
