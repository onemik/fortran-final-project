module Builder
    use NumberKinds
    use Types
    use Parameters
    use Geometry
    implicit none
    private
    public :: BuildMolecule

    contains


    subroutine BuildMolecule(mol)
        type(Molecule), intent(inout) :: mol
        
        call BuildBonds(mol)
    

    end subroutine


    !building bonds - number of them and the types of bond 
    subroutine BuildBonds(mol)
        type(Molecule), intent(inout) :: mol

        integer :: NumberOfAtoms
        integer :: i,j
        integer :: NumberOfBonds

        type(Bond), allocatable :: tmpBonds(:) !temporary to store bonds
        real(KREAL) :: r 

        NumberOfAtoms = size(mol%atoms)

        !maximum possible number of bonds is (#atoms)*(#atoms-1)/2
        allocate(tmpBonds(NumberOfAtoms*(NumberOfAtoms-1)/2))

        NumberOfBonds = 0 !starting with 0 bonds

        !looping through the distance between each potential bond
        do i=1, NumberOfAtoms-1
            do j=i+1, NumberOfAtoms
                r = CalculateDistance(mol%atoms(i), mol%atoms(j))

                !checking 
                if (IsBonded(mol%atoms(i)%symbol, mol%atoms(j)%symbol, r)) then
                    NumberOfBonds = NumberOfBonds+1
                    tmpBonds(NumberOfBonds)%i = i
                    tmpBonds(NumberOfBonds)%j = j
                    tmpBonds(NumberOfBonds)%kind = BondKind(mol%atoms(i)%symbol, mol%atoms(j)%symbol)
                end if

            end do
        end do


        allocate(mol%bonds(NumberOfBonds)) 
        mol%bonds = tmpBonds(1:NumberOfBonds)

        deallocate(tmpBonds) !save storage space

    end subroutine


    !function to check whether two atoms are bonded
    logical function IsBonded(symbol1, symbol2, r)
        character(len=*), intent(in) :: symbol1, symbol2
        real(KREAL), intent(in) :: r 

        real(KREAL), parameter :: tol = 0.20_KREAL

        IsBonded = .false. !by default not bonded
        !unless the distance between the atoms is close enough to the bond length: 
        if (trim(symbol1) == 'C' .and. trim(symbol2) == 'C') then
            if (abs(r-r0CC)<tol) IsBonded = .true. !needs to be within tolerance 
        elseif ((trim(symbol1) == 'C' .and. trim(symbol2) == 'H') .or. &
            (trim(symbol1) == 'H' .and. trim(symbol2) == 'C')) then
            if (abs(r-r0CH)<tol) IsBonded = .true.
        end if

    end function

    !define the kind of the bond based on the atoms it connects
    function BondKind(symbol1, symbol2) result(kind)
        character(len=*), intent(in) :: symbol1, symbol2
        character(len=2) :: kind

        if (trim(symbol1) == 'C' .and. trim(symbol2) == 'C') then
            kind = 'CC'
        else 
            kind = 'CH'
        end if

    end function



end module