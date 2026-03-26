module Energy
    use NumberKinds
    use Types
    use Geometry
    use Parameters
    implicit none
    private
    public :: StretchEnergy

    contains

    !calculate stretch energy 
    function StretchEnergy(mol) result(Estretch)
        type(Molecule), intent(in) :: mol
        real(KREAL) :: Estretch

        integer :: a
        integer :: i,j 
        real(KREAL) :: r, r0, k

        Estretch = 0.0_KREAL

        do a = 1, size(mol%bonds)
            i = mol%bonds(a)%i
            j = mol%bonds(a)%j

            r = CalculateDistance(mol%atoms(i), mol%atoms(j))

            !choosing the correct constants based on bond type
            if (trim(mol%bonds(a)%kind) == 'CC') then
                r0=r0CC
                k=kCC
            else 
                r0=r0CH
                k=kCH
            end if
        
            Estretch = Estretch + k * (r-r0)**2

        end do



    end function



end module
