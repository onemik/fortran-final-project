module Parameters
    use NumberKinds
    implicit none
    private
    public :: kCC, kCH, r0CC, r0CH, kAngle, theta0placeholder, Aij, Bij, V1_tors, n_tors, gamma_tors, kB

    real(KREAL), parameter :: pi = 3.141592653589793_KREAL

    real(KREAL), parameter :: r0CC = 1.54_KREAL
    real(KREAL), parameter :: r0CH = 1.09_KREAL 
    real(KREAL), parameter :: kCC = 317.0_KREAL !from the paper
    real(KREAL), parameter :: kCH = 340.0_KREAL !from the paper

    real(KREAL), parameter :: kAngle = 35.0_KREAL !from the paper 
    real(KREAL), parameter :: theta0placeholder = 109.5_KREAL * pi / 180.0_KREAL !for eq angle

    real(KREAL), parameter :: Aij = 0.1d0 !for nonbonding energy formula, like in previous exercise
    real(KREAL), parameter :: Bij = 0.2d0

    real(KREAL), parameter :: V1_tors = 0.155_KREAL !guesstimate 
    integer, parameter :: n_tors = 3 !guesstimate 
    real(KREAL), parameter :: gamma_tors = 0.0_kreal !phase=0 guesstimate 

    real(KREAL), parameter :: kB = 0.0019872041_KREAL !boltzmann constant in kcal/mol/K

end module