
MODULE data_structure
    use constant, only : dp 
    IMPLICIT NONE
    INTEGER, PARAMETER :: ndim=3, ntype = 3
    INTEGER, PARAMETER :: n1st_max=12, ndim3 = ndim
    
    !INTEGER, PARAMETER :: natms_max=1000000
    !INTEGER, PARAMETER :: na=natms_max/5
  
    INTEGER, PARAMETER :: no_aa_neigh=6, no_ab_neigh=8, no_ao_neigh=12, no_bo_neigh=6
    REAL(dp),PARAMETER :: PI = 3.141592653_dp
    REAL(dp) :: dab_cut =10.0
    REAL(dp) :: daa_cut =10.0
    REAL(dp) :: dao_cut =6.0
    REAL(dp) :: dbo_cut = 6.0

    REAL(dp) :: dab_min =4.7 ! 3.2
    REAL(dp) :: daa_min =4.7 ! 3.8
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
    INTEGER :: tndump, nup_spin, ndown_spin
  
    REAL(dp), DIMENSION(:,:), ALLOCATABLE ::   a_atms,b_atms,a_disp, b_disp
    REAL(dp), DIMENSION(:,:), allocatable :: a_atms_symm, b_atms_symm  
    REAL(dp), DIMENSION(:,:), ALLOCATABLE ::  a_pol, b_pol, pol_unit

    REAL(dp), DIMENSION(:,:), ALLOCATABLE ::  a_pol_up, a_pol_down
    
    REAL(dp), DIMENSION(:), ALLOCATABLE :: a_pol_dir

    REAL(dp), DIMENSION(:,:), ALLOCATABLE ::  o_atms
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: sp
    REAL(dp), DIMENSION(:,:,:),ALLOCATABLE ::   pol_ab 
    INTEGER, DIMENSION(:, :),ALLOCATABLE ::     aa_neigh
    INTEGER, DIMENSION(:, :),ALLOCATABLE ::     ab_neigh
    INTEGER, DIMENSION(:, :),ALLOCATABLE ::     ao12_neigh
    INTEGER, DIMENSION(:, :),ALLOCATABLE ::     bo6_neigh

    INTEGER, DIMENSION(:), ALLOCATABLE ::  pol_up_index, pol_down_index
  END MODULE data_structure
