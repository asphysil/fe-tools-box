
SUBROUTINE interaction_3atms_xy1(latt_sup, nsti, b1_atms, o22_atms, angle1_ti_o_ti)

   USE my_constants, ONLY : dp_real, dp_int
   USE data_cubic40, ONLY :  neighbour_atm1_atm2,neighbour_axis,&
                  projection_car_xy1,  dm_Zint_angle, dm_int_angle

   IMPLICIT NONE

   INTEGER(dp_int), INTENT(IN) :: nsti

   REAL(dp_real), INTENT(IN) ::  b1_atms(4, nsti), o22_atms(4, nsti), latt_sup(3,3)

   INTEGER(dp_int), INTENT(OUT) :: angle1_ti_o_ti(3, nsti)

   INTEGER(dp_int) ::  neigh_ab(8, nsti), NeighAlongAxis(nsti), atmsij_tio(nsti)
   
   REAL (dp_real) :: dmin, vec_unit(3)

   INTEGER(dp_int) :: l, k,  nij_tio

   DO l =1, nsti
    angle1_ti_o_ti(1, l) = l
      DO k =2, 3
      angle1_ti_o_ti(k, l) = -1
      ENDDO 
   ENDDO

! Fe -O1
        dmin=4.0/2.0 + 1.0
        CALL neighbour_atm1_atm2(dmin, nsti, b1_atms, o22_atms, neigh_ab)
        
!          DO l =1, nsti 
!                  PRINT*, 'Ti->O', l, (neigh_ab(l, k), k=1,6)
!          ENDDO
!          PRINT*,"----------"
     
       CALL projection_car_xy1(nsti, b1_atms, o22_atms, neigh_ab, NeighAlongAxis)
       
   !PRINT*,   'Ti->O did not find O-atom if it is 0 '   
   !DO l =1, nsti 
    !  PRINT*,'2. Ti->O (+axis) ', l,  (NeighAlongAxis(l))
   !ENDDO
   !PRINT*,"----------"
      
   ! Sorting which Ti did not find O-atoms
       nij_tio = 0
        DO l =1, nsti 
            if (NeighAlongAxis(l) /=0) THEN
                angle1_ti_o_ti(2, l) = NeighAlongAxis(l)
                
            ELSE
                nij_tio = nij_tio + 1
                atmsij_tio(nij_tio) = l
            ENDIF
        ENDDO 

        ! DO l =1, nij_tio
       !      PRINT*, 'Ti did not find O -atoms', atmsij_tio(l)
        !  ENDDO

           !DO k =1, nsti
            !PRINT*, (angle_ti_o_ti(k,l), l=1,3)
           !ENDDO

        
        !PRINT*, '****Ti-----O*****'
        dmin = 4.0/2.0  + 1.0
        vec_unit = (/1.0, 0.0, 0.0/)
        CALL neighbour_axis(vec_unit, dmin, latt_sup, nsti, b1_atms, o22_atms, nij_tio, atmsij_tio, NeighAlongAxis)
        
        DO l=1, nij_tio 
           k = atmsij_tio(l)
           angle1_ti_o_ti(2, k) = NeighAlongAxis(k)
           !PRINT*, 'Ti->O (+axis)', k,  (NeighAlongAxis(k))
        ENDDO
        !PRINT*,"***----------***"
        !DO l =1, nsti 
        !    PRINT*,'2. Ti->O (+axis) ', l,  angle_ti_o_ti(l,2)
        ! ENDDO
        ! PRINT*,"***----------***"

        DO l =1, nsti
            IF (angle1_ti_o_ti(2, k) ==-1) THEN 
                PRINT*, l, 'Ti atom did find O atom, Please check program'
            ENDIF

        ENDDO
        CALL dm_int_angle(nsti, latt_sup, b1_atms, o22_atms, angle1_ti_o_ti)
    
END SUBROUTINE interaction_3atms_xy1
