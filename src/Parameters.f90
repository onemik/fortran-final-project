module Parameters
    use NumberKinds
    implicit none
    private
    public :: kCC, kCH, r0CC, r0CH, kAngle, theta0placeholder

    real(KREAL), parameter :: pi = 3.141592653589793_KREAL

    real(KREAL), parameter :: r0CC = 1.54_KREAL
    real(KREAL), parameter :: r0CH = 1.09_KREAL 
    real(KREAL), parameter :: kCC = 317.0_KREAL !from the paper but double check 
    real(KREAL), parameter :: kCH = 340.0_KREAL !from the paper but double check

    real(KREAL), parameter :: kAngle = 35.0_KREAL !from the paper 
    real(KREAL), parameter :: theta0placeholder = 109.5_KREAL * pi / 180.0_KREAL !for eq angle

    real(KREAL), parameter :: Aij = 0.1d0 !for nonbonding energy formula, like in previous exercise
    real(KREAL), parameter :: Bij = 0.2d0
end module