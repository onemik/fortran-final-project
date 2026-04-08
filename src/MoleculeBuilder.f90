module MoleculeBuilder
    use NumberKinds
    use Types
    use Parameters
    use Geometry
    implicit none
    private
    public :: BuildMolecule

    contains

    !combine all builder subroutines to build up the molecule, 
    !i.e. its bonds, angles between bonds and torsions
    subroutine BuildMolecule(mol)
        type(Molecule), intent(inout) :: mol
        
        call BuildBonds(mol)
        call BuildBondAngles(mol)
        call BuildBondTorsions(mol)
        call BuildNonBondedPairs(mol)
    
    end subroutine


    !building a list of bonds (number of them and the types of bond) in the molecule
    subroutine BuildBonds(mol)
        type(Molecule), intent(inout) :: mol

        integer :: NumberOfAtoms
        integer :: i,j
        integer :: NumberOfBonds

        type(Bond), allocatable :: tmpBonds(:) !temporary to store bonds
        real(KREAL) :: r !distance 

        NumberOfAtoms = size(mol%atoms)

        !maximum possible number of bonds is (#atoms)*(#atoms-1)/2
        allocate(tmpBonds(NumberOfAtoms*(NumberOfAtoms-1)/2))

        NumberOfBonds = 0 !starting with 0 bonds

        !looping through the distance between each potential bond
        do i=1, NumberOfAtoms-1
            do j=i+1, NumberOfAtoms
                r = CalculateDistance(mol%atoms(i), mol%atoms(j))

                !checking if a bond forms between these atoms and if yes - store
                if (IsBonded(mol%atoms(i)%symbol, mol%atoms(j)%symbol, r)) then
                    NumberOfBonds = NumberOfBonds+1
                    tmpBonds(NumberOfBonds)%i = i
                    tmpBonds(NumberOfBonds)%j = j
                    tmpBonds(NumberOfBonds)%kind = BondKind(mol%atoms(i)%symbol, mol%atoms(j)%symbol)
                end if

            end do
        end do

        !allocate the bonds from the temporary array to the bonds part of the type
        if (allocated(mol%bonds)) deallocate(mol%bonds) !in case it was allocated before
        allocate(mol%bonds(NumberOfBonds)) 
        mol%bonds = tmpBonds(1:NumberOfBonds)

        deallocate(tmpBonds) !save storage space

    end subroutine


    !function to check whether two atoms should be considered bonded
    logical function IsBonded(symbol1, symbol2, r)
        character(len=*), intent(in) :: symbol1, symbol2
        real(KREAL), intent(in) :: r 

        real(KREAL), parameter :: tol = 0.20_KREAL !tolerance for what's still considered bonded

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


    !building the list of angles in the molecule
    !angle i-j-k exists when atom j has two neighbours i and k
    subroutine BuildBondAngles(mol)
        type(Molecule), intent(inout) :: mol

        integer :: NumberOfAtoms
        integer :: NumberOfBonds
        integer :: NumberOfNeighbours
        integer :: NumberOfAngles

        integer :: j !middle atom
        integer :: b !bond index
        integer :: n1,n2 !neighbours index
        integer, allocatable :: Neighbours(:) !temporary storage for neighbours of one atom

        type(Angle), allocatable :: tmpAngles(:)

        NumberOfAtoms = size(mol%atoms)
        NumberOfBonds = size(mol%bonds)

        allocate(tmpAngles(NumberOfBonds*NumberOfBonds)) !max possible number of angles
        allocate(Neighbours(NumberOfAtoms)) !max size one atom could be bonded to
        NumberOfAngles = 0

        !make a loop over all atoms as a possible middle atom j 
        do j=1, NumberOfAtoms
            NumberOfNeighbours = 0 !re-starting with 0 per atom

            !running through each bond to determine if it's bonded with j (center atom)
            do b = 1, NumberOfBonds
                if (mol%bonds(b)%i == j) then !if the first atom of the given bond is one i'm looking for
                    NumberOfNeighbours = NumberOfNeighbours + 1
                    Neighbours(NumberOfNeighbours) = mol%bonds(b)%j !store the neighbouring atom index 
                elseif (mol%bonds(b)%j == j) then !if the second atom of the bond is one i'm looking for
                    NumberOfNeighbours = NumberOfNeighbours + 1
                    Neighbours(NumberOfNeighbours) = mol%bonds(b)%i 
                end if
            end do

            !now loop through every paid of neighbours to form the angle
            if (NumberOfNeighbours>=2) then
                do n1=1, NumberOfNeighbours-1
                    do n2=n1+1, NumberOfNeighbours

                        NumberOfAngles = NumberOfAngles + 1
                        tmpAngles(NumberOfAngles)%i = Neighbours(n1)
                        tmpAngles(NumberOfAngles)%j = j
                        tmpAngles(NumberOfAngles)%k = Neighbours(n2)

                    end do
                end do
            end if

        end do

        if (allocated(mol%angles)) deallocate(mol%angles)
        allocate(mol%angles(NumberOfAngles))
        mol%angles = tmpAngles(1:NumberOfAngles)

        !deallocate to save space
        deallocate(tmpAngles)
        deallocate(Neighbours)

    end subroutine

    !torsion is where for atoms connected in a chain i-j-k-l
    !i.e. rotation around the middle bond (between j-k), so measuring
    !the angle between planes i-j-k and j-k-l
    subroutine BuildBondTorsions(mol)
        type(Molecule), intent(inout) :: mol

        integer :: NumberOfAtoms
        integer :: NumberOfBonds
        integer :: NumberOfTorsions 

        integer :: b !bond index
        integer :: i,j,k,l
        integer :: NumberOfNeighboursJ, NumberOfNeighboursK
        integer :: n1,n2
        integer, allocatable :: NeighboursJ(:), NeighboursK(:)
        type(Torsion), allocatable :: tmpTorsions(:)

        NumberOfAtoms = size(mol%atoms)
        NumberOfBonds = size(mol%bonds)

        allocate(tmpTorsions(NumberOfBonds*NumberOfBonds)) !max possible number of torsions
        allocate(NeighboursJ(NumberOfAtoms)) !max possible number of neighbours 
        allocate(NeighboursK(NumberOfAtoms))

        NumberOfTorsions = 0 !starting value

        !loop over each bond as the central bond j-k 
        do b=1, NumberOfBonds
            !define central bond
            j = mol%bonds(b)%i
            k = mol%bonds(b)%j

            !find all neighrbours of j excluding k, ie possible "i" atoms
            call GetNeighbours(mol%bonds,j,k,NeighboursJ, NumberOfNeighboursJ)
            !same but for possible "l" atoms
            call GetNeighbours(mol%bonds,k,j,NeighboursK, NumberOfNeighboursK)

            !combine neighbours of j and k to form torsions i-j-k-l
            do n1 = 1, NumberOfNeighboursJ
                i = NeighboursJ(n1)

                do n2=1, NumberOfNeighboursK
                    l = NeighboursK(n2)

                    !store torsion defined by atoms ijkl
                    if (i/=l) then
                        NumberOfTorsions = NumberOfTorsions + 1
                        tmpTorsions(NumberOfTorsions)%i = i
                        tmpTorsions(NumberOfTorsions)%j = j
                        tmpTorsions(NumberOfTorsions)%k = k
                        tmpTorsions(NumberOfTorsions)%l = l
                    end if
                end do
            end do
        end do

        !storing values and freeing up storage space
        if (allocated(mol%torsions)) deallocate(mol%torsions)
        allocate(mol%torsions(NumberOfTorsions))
        mol%torsions = tmpTorsions(1:NumberOfTorsions)

        deallocate(tmpTorsions)
        deallocate(NeighboursJ)
        deallocate(NeighboursK)

    end subroutine

    !helper subroutine to look for all existing neighbours of e.g. atom j  
    !excluding k (because it's already connected to j from the other side)
    subroutine GetNeighbours(bonds, centerAtom, excludedAtom, neighbours, NumberOfNeighbours)
        type(Bond), intent(in) :: bonds(:)
        integer, intent(in) :: centerAtom, excludedAtom 
        integer, intent(out) :: neighbours(:)
        integer, intent(out) :: NumberOfNeighbours

        integer :: b !bond index

        !starting value
        NumberOfNeighbours = 0 

        !loop through the bonds to see whether it involves the selected center atom 
        !but that it also excludes the other atom in the bond
        do b=1, size(bonds)

            !center atom being the first atom in the bond & neighbour is the second atom
            if (bonds(b)%i == centerAtom) then 
                if (bonds(b)%j /= excludedAtom) then !only store if it's not the excluded atom
                    NumberOfNeighbours = NumberOfNeighbours + 1
                    neighbours(NumberOfNeighbours) = bonds(b)%j
                end if
            
            !center atom is the second atom in the bond & neighbour is the first atom
            elseif (bonds(b)%j == centerAtom) then
                if (bonds(b)%i /= excludedAtom) then
                    NumberOfNeighbours = NumberOfNeighbours + 1
                    neighbours(NumberOfNeighbours) = bonds(b)%i
                end if

            end if

        end do


    end subroutine

    !make a list of pairs of atoms that are not bonded (for nonbonded energy)
    subroutine BuildNonBondedPairs(mol)
        type(Molecule), intent(inout) :: mol

        integer :: NumberOfAtoms
        integer :: NumberOfPairs !nonbonded
        integer :: i,j
        type(NoBond), allocatable :: tmpPairs(:)

        NumberOfAtoms = size(mol%atoms)

        allocate(tmpPairs(NumberOfAtoms*(NumberOfAtoms-1)/2)) !max possible number of pairs
        NumberOfPairs = 0

        !loop through each possible combination of atoms and check if they're bonded
        do i=1, NumberOfAtoms-1
            do j=i+1, NumberOfAtoms
                if (.not. AreBonded(i,j,mol%bonds)) then
                    NumberOfPairs = NumberOfPairs + 1
                    tmpPairs(NumberOfPairs)%i = i
                    tmpPairs(NumberOfPairs)%j = j
                end if
            end do
        end do

        !store in the molecule type
        if (allocated(mol%nobonds)) deallocate(mol%nobonds)
        allocate(mol%nobonds(NumberOfPairs))
        mol%nobonds = tmpPairs(1:NumberOfPairs)

        deallocate(tmpPairs) !deallocate to save space

        
    end subroutine

    !function to check whether two atoms are bonded (i.e. if the bond exists)
    logical function AreBonded(i,j,bonds)
        integer, intent(in) :: i,j
        type(Bond), intent(in) :: bonds(:)

        integer :: b !bond index

        AreBonded = .false.

        do b = 1, size(bonds)
            if ((bonds(b)%i == i .and. bonds(b)%j ==j) .or. &
                (bonds(b)%i == j .and. bonds(b)%j ==i)) then
                AreBonded = .True.
                return
            end if
        end do

    end function


end module