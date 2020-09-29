module variables

    implicit none

    ! CONSTANTS
    double complex, parameter :: eye = dcmplx(0.d0, 1.d0)    ! Complex i
    double precision, parameter :: pi = 3.14159265359d0      ! Pie is delicious
    double precision, parameter :: sol = 137.03599908381668d0! Speed of light
    double precision, parameter :: eps0 = 0.079577471546d0   ! Vacuum permitt.

    ! SYSTEM SPECIFIC CONSTANTS AND ENERGY LEVEL ARRAY
    integer, parameter :: S = 2                              ! # of el.st.
    double precision, parameter :: mu = 1.034d0              ! Dipole moment
    double precision, parameter :: L = 236215.76557822127d0  ! Cavity length
    double precision, dimension(S), parameter :: eps = [-0.6738d0, -0.2798d0]

    ! VARIABLES
    integer :: F              ! Number of bath DoFs
    integer :: ntraj          ! Number of trajectories
    integer :: tsteps         ! Number of TS per traj.
    integer :: init_state     ! Initial electronic state
    double precision :: dt    ! TS duration
    double precision :: V0    ! State-independent potential
    character(len=14):: intgt ! Integrator type

    ! ARRAYS
    double precision :: V(S,S)                      ! Potential energy matrix
    double precision :: XE(S), PE(S)                ! Mapping variables
    double complex, allocatable :: rho(:,:,:)       ! Subsystem density matrix 
    double precision, allocatable :: G0(:)          ! State-independent force
    double precision, allocatable :: xn(:), pn(:)   ! Bath positions & momenta
    double precision, allocatable :: omega(:), c(:) ! Bath freq. & coefficients

    ! LAPACK PARAMETERS
    integer :: info, lenwork                        ! LAPACK integer parameters
    double precision, allocatable :: work(:)        ! Work array for LAPACK

end module variables
