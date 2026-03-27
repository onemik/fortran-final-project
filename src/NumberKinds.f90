!defining variables precision throughout the program

module NumberKinds
    implicit none
    integer, parameter :: KREAL = kind(0.d0) !double precision
end module