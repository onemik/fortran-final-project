module Parameters
    use NumberKinds
    implicit none
    private
    public :: kCC, kCH, r0CC, r0CH

    real(KREAL), parameter :: r0CC = 1.54_KREAL
    real(KREAL), parameter :: r0CH = 1.09_KREAL 
    real(KREAL), parameter :: kCC = 317.0_KREAL !double check 
    real(KREAL), parameter :: kCH = 340.0_KREAL !double check

end module