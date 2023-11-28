MODULE readfiles
 use constant, only : dp
 use data_structure, only : latt_vec, a_atms, b_atms, o_atms, sp, tndump,&
                            angle_index 
    IMPLICIT NONE 
private
public :: read_dumpxyz, &
          read_dumpxyz_spin, read_angle_info_spin
contains

SUBROUTINE read_dumpxyz_spin(fileID,na, nb, no)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: fileID
   INTEGER, INTENT(IN) :: na, nb, no
   REAL(dp) :: x,y, start, finish, junk1, junk2
   INTEGER :: i, j, ntot
       !reading .xyz file
  
  !WRITE(15, *) '# Px   Py,   Pz,  and  P = SQRT(px^2 +py^2 + pz^2)'
  !!*********Reading files***********!!
  
     READ(fileID,*)
     READ(fileID,*) tndump
     READ(fileID,*)
  
     READ(fileID,*) ntot
     !PRINT*, ntot
  
     READ(fileID,*)
     DO i=1,3
       READ(fileID,*) x, y
       latt_vec(i,i)=ABS(y-x)
     ENDDO
  
     READ(fileID,*)
  DO i=1, na
      READ (fileID,*) junk1, a_atms(1,i), a_atms(2,i), a_atms(3,i)
   ENDDO
  
   DO i=1, nb
    READ(fileID,*) junk1,  b_atms(1,i), b_atms(2,i), b_atms(3,i), (sp(j,i), j=1,3)
  !, sp(1,i), sp(2,i), sp(3,i)
   ENDDO
  
  DO i=1, no
     READ(fileID,*) junk1,  o_atms(1,i), o_atms(2,i), o_atms(3,i)
   ENDDO
  
  END SUBROUTINE read_dumpxyz_spin 
  

  SUBROUTINE read_dumpxyz_spin_x1(fileID,na, nb, no)
   IMPLICIT NONE 
   INTEGER, INTENT(IN) :: fileID
   INTEGER, INTENT(IN) :: na, nb, no  
   REAL(dp) :: x,y, start, finish
   INTEGER :: i, j, ntot 
   INTEGER :: id !, type   
  !reading .xyz file
  !WRITE(15, *) '# Px   Py,   Pz,  and  P = SQRT(px^2 +py^2 + pz^2)'
  !!*********Reading files***********!!
     READ(fileID,*)
     READ(fileID,*) tndump
     READ(fileID,*)
     READ(fileID,*) ntot
     READ(fileID,*)
     DO i=1,3
       READ(fileID,*) x, y
       latt_vec(i,i)=ABS(y-x)
     ENDDO
  
     READ(fileID,*)
  DO i=1, na
      READ (fileID,*) id,  a_atms(1,i), a_atms(2,i), a_atms(3,i)
   ENDDO
  
   DO i=1, nb
    READ(fileID,*) id, b_atms(1,i), b_atms(2,i), b_atms(3,i), sp(1,i), sp(2,i), sp(3,i)
   ENDDO
    
  DO i=1, no
     READ(fileID,*) id,  o_atms(1,i), o_atms(2,i), o_atms(3,i)
   ENDDO
  
 
  !!*******************************!  
  END SUBROUTINE read_dumpxyz_spin_x1

  
SUBROUTINE read_dumpxyz(fileID,na, nb, no)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: fileID 
 INTEGER, INTENT(IN) :: na, nb, no
 REAL(dp) :: x,y, start, finish, junk1, junk2
 INTEGER :: i, j, ntot
     !reading .xyz file

!WRITE(15, *) '# Px   Py,   Pz,  and  P = SQRT(px^2 +py^2 + pz^2)'
!!*********Reading files***********!!

   READ(fileID,*)
   READ(fileID,*) tndump
   READ(fileID,*)

   READ(fileID,*) ntot
   !PRINT*, ntot

   READ(fileID,*)
   DO i=1,3
     READ(fileID,*) x, y
     latt_vec(i,i)=ABS(y-x)
   ENDDO

   READ(fileID,*)
DO i=1, na
    READ (fileID,*) junk1, junk2, a_atms(1,i), a_atms(2,i), a_atms(3,i)
 ENDDO

 DO i=1, nb
  READ(fileID,*) junk1, junk2, b_atms(1,i), b_atms(2,i), b_atms(3,i)
  !, sp(1,i), sp(2,i), sp(3,i)
 ENDDO

DO i=1, no
   READ(fileID,*) junk1, junk2, o_atms(1,i), o_atms(2,i), o_atms(3,i)
 ENDDO

END SUBROUTINE read_dumpxyz 


SUBROUTINE read_dumpxyz_drew(fileID,na, nb, no)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: fileID
   INTEGER, INTENT(IN) :: na, nb, no

   REAL :: x,y, z, start, finish
   REAL(dp) :: x_dp, y_dp, z_dp 
   REAL(dp) :: atm1(3,nb/2), atm2(3,nb/2) 
   INTEGER :: i, j, ntot, id, atm_type, m1, m2, m3, m4 
       !reading .xyz file
  
  !WRITE(15, *) '# Px   Py,   Pz,  and  P = SQRT(px^2 +py^2 + pz^2)'
  !!*********Reading files***********!!
  m1=0
  m2=0
  m3=0
  m4=0

     READ(fileID,*)
     READ(fileID,*) tndump
     READ(fileID,*)
  
     READ(fileID,*) ntot
     !PRINT*, ntot
  
     READ(fileID,*)
     DO i=1,3
       READ(fileID,*) x, y
       latt_vec(i,i)=ABS(y-x)
     ENDDO
  
     READ(fileID,*)
     DO i = 1, ntot 
      READ (fileID,*) id, atm_type, x,y, z
     ! PRINT*, id, atm_type, x,y, z
x_dp = x
y_dp  = y
z_dp = z
      if (atm_type==1) then 
         m1 = m1 + 1
         a_atms(1,m1) = x_dp
         a_atms(2,m1) = y_dp 
         a_atms(3,m1) = z_dp 
         
      elseif ( atm_type==2) then 
         m2 = m2 + 1 
         atm1(1, m2) = x_dp 
         atm1(2, m2) = y_dp
         atm1(3, m2) = z_dp

      elseif ( atm_type==3) then
         m3 = m3 + 1 
         atm2(1, m3) = x_dp 
         atm2(2, m3) = y_dp
         atm2(3, m3) = z_dp
      elseif ( atm_type==4) then
         m4 = m4 + 1
         o_atms(1,m4) = x_dp 
         o_atms(2,m4) = y_dp 
         o_atms(3,m4) = z_dp
      else
         print*, "Error reading files"
      endif 
   enddo  

   IF ( ( m1 .NE. na ) .OR. (m2+m3 .NE. nb) .OR.(m4 .NE. no)) then 
      print*, " Error in reading file and atom number are not match"
   endif 

   m1 =0
   DO i = 1, m2  
      m1 = m1 + 1
      b_atms(1, m1) = atm1(1,i) 
      b_atms(2, m1) = atm1(2,i) 
      b_atms(3, m1) = atm1(3,i) 
   enddo

   DO i = 1, m3  
      m1 = m1 + 1
      b_atms(1, m1) = atm2(1,i) 
      b_atms(2, m1) = atm2(2,i) 
      b_atms(3, m1) = atm2(3,i)
   enddo

  END SUBROUTINE read_dumpxyz_drew 


 
SUBROUTINE read_angle_info_spin(fileID,nbi,nangle)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: fileID
   INTEGER, INTENT(IN) :: nangle, nbi
   REAL(dp) :: x,y, start, finish, junk1, junk2
   INTEGER :: i, j, njunk
   INTEGER ::  ntot, nn 
   
   CHARACTER :: name_angle

       !reading .xyz file
  
  !WRITE(15, *) '# Px   Py,   Pz,  and  P = SQRT(px^2 +py^2 + pz^2)'
  !!*********Reading files***********!!
 nn = nbi+nbi
  ntot =nn +nangle

   njunk=20+ntot
  DO i = 1, njunk
     READ(fileID,*)
  ENDDO
  READ(fileID,*)
  READ(fileID,*) name_angle
  READ(fileID,*)

  IF (name_angle=='A' .OR. name_angle=='a') THEN 
   PRINT*, "Correctly reading  angle infomation"
  ELSE
   PRINT*, " Error reading  angle infomation"
  ENDIF
  DO i=1, nangle
      READ (fileID,*) junk1, junk2, (angle_index(j,i), j=1,3)
      angle_index(1,i) = angle_index(1,i)-nbi
      angle_index(2,i) = angle_index(2,i)- nn 
      angle_index(3,i) = angle_index(3,i)-nbi
      !print'(3I5)', (angle_index(j,i), j=1,3)
   ENDDO

  END SUBROUTINE read_angle_info_spin 

END MODULE readfiles
