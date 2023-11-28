SUBROUTINE supercell_cubic(nu, nuti, latt_vec, atmu, nsti, a_atms, b_atms, o1_atms, o2_atms, o3_atms)
    USE dtype, ONLY : dp_real, dp_int, nx, ny, nz
    IMPLICIT NONE

    INTEGER(dp_int), INTENT(IN) :: nu, nsti, nuti 
        REAL(dp_real), INTENT(IN) ::latt_vec(3,3), atmu(nu, 3)
    
        REAL (dp_real), INTENT(OUT) :: a_atms(nsti, 4), b_atms(nsti,4) !, o_atms(nsti*3, 4)

    REAL(dp_real), INTENT(OUT) :: o1_atms(nsti, 4), o2_atms(nsti, 4), o3_atms(nsti,4)

    REAL(dp_real) :: T(3), a_atmu(nuti,3), b_atmu(nuti,3),&
                      o1_atmu(nuti,3), o2_atmu(nuti,3), o3_atmu(nuti,3)

    REAL(dp_real) :: a, b, c

    INTEGER (dp_int) :: ix, iy, iz, n, i 

    n = 0
    DO i =1, nuti
       n = n + 1
    a_atmu(i,:) = atmu(n,:)
    ENDDO

    DO i =1, nuti
        n = n + 1
        b_atmu(i,:) = atmu(n,:)
    ENDDO

    DO i = 1, nuti
     n = n + 1
    o1_atmu(i,:) = atmu(n,:)
    ENDDO

    DO i = 1, nuti
      n = n + 1
       o2_atmu(i,:) = atmu(n,:)
    ENDDO

    DO i=1, nuti
    n = n + 1
    o3_atmu(i,:) = atmu(n,:)
    ENDDO

    a  = latt_vec(1,1)
    b  = latt_vec(2,2)
    c  = latt_vec(3,3)
    n = 0

print*, "*******", nx, ny ,nz, "********"

   DO i = 1, nuti


    DO iz =0, nz-1
        DO iy = 0, ny-1
            DO ix = 0, nx-1

            T=(/DBLE(ix)*a, DBLE(iy)*b, DBLE(iz)*c/)     
                n = n + 1    

            a_atms(n, 1) = DBLE(iz)
            b_atms(n, 1) = DBLE(iz)
            o1_atms(n, 1) = DBLE(iz)
            o2_atms(n, 1) = DBLE(iz)
            o3_atms(n, 1) = DBLE(iz)

            a_atms(n, 2) = a_atmu(i,1)*a + T(1)
            a_atms(n, 3) = a_atmu(i,2)*b + T(2)
            a_atms(n, 4) = a_atmu(i,3)*c + T(3)


            b_atms(n, 2) = b_atmu(i,1)*a + T(1)
            b_atms(n, 3) = b_atmu(i,2)*b + T(2)
            b_atms(n, 4) = b_atmu(i,3)*c + T(3)


            o1_atms(n, 2) = o1_atmu(i,1)*a + T(1)
            o1_atms(n, 3) = o1_atmu(i,2)*b + T(2)
            o1_atms(n, 4) = o1_atmu(i,3)*c + T(3)


            o2_atms(n, 2) = o2_atmu(i,1)*a + T(1)
            o2_atms(n, 3) = o2_atmu(i,2)*b + T(2)
            o2_atms(n, 4) = o2_atmu(i,3)*c + T(3)


            o3_atms(n, 2) = o3_atmu(i,1)*a  + T(1)
            o3_atms(n, 3) = o3_atmu(i,2)*b  + T(2)
            o3_atms(n, 4) = o3_atmu(i,3)*c  + T(3)
            
                ENDDO
            ENDDO
        ENDDO
     ENDDO
!PRINT*, 'nti total', n, 'nsti=', nsti
!n=0
!DO i =1, nsti 
!    n = n + 1
!    !o_atms(n,:) = o1_atms(n,:)
!ENDDO
!
!DO i =1, nsti 
!    n = n + 1
!    !o_atms(n,:) = o2_atms(n,:)
!ENDDO
!
!DO i =1, nsti 
!    n = n + 1
!    !o_atms(n,:) = o3_atms(n,:)
!ENDDO
!PRINT*, 'no total=', n, 'nti=', nsti*3


END SUBROUTINE supercell_cubic


