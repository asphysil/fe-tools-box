MODULE constant
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp=SELECTED_REAL_KIND(15) !, int64=SELECTED_INT_KIND(15)
  !REAL(dp), PARAMETER :: zqa=4.37d0, zqb=3.49d0, zqo1=-2.61d0, zqo2=-2.61d0, zqo3=-2.61d0 ! BFO
   REAL(dp), PARAMETER :: zqa=1.0d0, zqb=1.0d0, zqo1=-1.0d0, zqo2=-1.0d0, zqo3=-1.0d0
  ! qo1= parallel to T-O bond (  Born Effective value large in parallel direction due to bond elongation and compression motion)
  ! qo2 and qo3 perpendicular to T-O bond (  Born Effective value small in perpendicular direction due to bond bending motion)
  REAL(dp),PARAMETER :: PI = 3.141592653_dp
END MODULE constant

MODULE user_fileID
  IMPLICIT NONE
  INTEGER, PARAMETER :: nmax_file=20
! file ID
  INTEGER, DIMENSION(nmax_file) :: fileID_all

  INTEGER :: fileID_read
  INTEGER :: fileID_totpol, fileID_locpol
  INTEGER :: fileID_ovito
  INTEGER :: fileID_pol_ovito
  !INTEGER ::  fileID_rtheta
  INTEGER :: file_ID_xsf1
  INTEGER :: file_ID_xsf2
  INTEGER :: fileID_angle
  INTEGER :: fileID_ovito_lmp
END MODULE user_fileID


MODULE data_structure
    use constant, only : dp 
    IMPLICIT NONE
    INTEGER, PARAMETER :: ndim=3, ntype = 3
    INTEGER, PARAMETER :: n1st_max=12, ndim3 = ndim

    !INTEGER, PARAMETER :: natms_max=1000000
    !INTEGER, PARAMETER :: na=natms_max/5
  
    INTEGER, PARAMETER :: no_aa_neigh=6, no_bb_neigh=6
    INTEGER,PARAMETER ::  no_ab_neigh=8, no_ao_neigh=12, no_bo_neigh=6

    REAL(dp) :: dab_cut =10.0
    REAL(dp) :: daa_cut =10.0
    REAL(dp) :: dao_cut =6.0
    REAL(dp) :: dbo_cut = 6.0

    REAL(dp) :: dab_min =4.7 ! 3.2
    REAL(dp) :: daa_min =4.7 ! 3.8
    REAL(dp) :: dbb_min =4.7 ! 3.8
    REAL(dp) :: dao_min =3.7 !2.81
    REAL(dp) :: dbo_min =3.7 !2.81

    REAL(dp) :: unitvol
    !INTEGER, DIMENSION(:) :: atm_info(27*natms1_max)
    !REAL(dp), DIMENSION(:) ::  d(27*natms1_max)
    !REAL(dp), DIMENSION(:,:) :: ab_neigh_coord(3, n1st_max, na)
    !
    REAL(dp),DIMENSION(ndim,ndim) :: latt_vec, ainv
    REAL(dp), DIMENSION(ndim) :: work, ipiv
    INTEGER :: info
    INTEGER :: tndump
  
    REAL(dp), DIMENSION(:,:), ALLOCATABLE ::   a_atms,b_atms,a_disp, b_disp
    REAL(dp), DIMENSION(:,:), ALLOCATABLE ::  a_pol, b_pol, pol_unit
    
    !REAL(dp), DIMENSION(:), ALLOCATABLE :: a_pol_dir

    REAL(dp), DIMENSION(:,:), ALLOCATABLE ::  o_atms
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: sp
    REAL(dp), DIMENSION(:,:,:),ALLOCATABLE ::   pol_ab 

    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: coord_fe1fe6

    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: coord_fe1ofe2
    
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: angle_index

    INTEGER, DIMENSION(:, :),ALLOCATABLE ::     aa_neigh
    INTEGER, DIMENSION(:, :),ALLOCATABLE ::     bb_neigh

    INTEGER, DIMENSION(:, :),ALLOCATABLE ::     ab_neigh
    INTEGER, DIMENSION(:, :),ALLOCATABLE ::     ao12_neigh
    INTEGER, DIMENSION(:, :),ALLOCATABLE ::     bo6_neigh

  END MODULE data_structure
