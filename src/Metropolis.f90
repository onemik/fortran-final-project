module Metropolis
    use NumberKinds
    use Types
    use ComputeEnergy
    use Parameters
    implicit none
    private
    public :: RunMetropolis

    contains

    !metropolis monte carlo minimization algorithm 
    subroutine RunMetropolis(mol, NumberOfSteps,r, temperature, bestEnergy, acceptanceRate)

        type(Molecule), intent(inout) :: mol 
        integer, intent(in) :: NumberOfSteps
        integer :: NumberOfAtoms

        real(KREAL), intent(in) :: r !step size/radius
        real(KREAL), intent(in) :: temperature 
        real(KREAL), intent(out) :: bestEnergy !lowest energy
        real(KREAL), intent(out) :: acceptanceRate !ratio of accepted moves 

        integer :: step, acceptedMoves
        real(KREAL) :: deltaE, Ei, Ey
        real(KREAL) :: randomNumber !for acceptance 

        !to store coordinates 
        real(KREAL), allocatable :: xi(:,:) !for rejection
        real(KREAL), allocatable :: y(:,:) !for trial
        real(KREAL), allocatable :: xBest(:,:) !for best result so far


        !allocating the coordinate array sizes based on number of atoms
        NumberOfAtoms = size(mol%atoms)
        allocate(xi(NumberOfAtoms, 3))
        allocate(y(NumberOfAtoms, 3))
        allocate(xBest(NumberOfAtoms, 3))

        !start from initial configuration xi (coords taken from the datafile)
        call saveCoordinates(mol, xi)
        call saveCoordinates(mol, xbest)
        
        Ei = TotalEnergy(mol)
        bestEnergy = Ei
        acceptedMoves = 0

        !main loop for metropolis algorithm 
        do step =1, NumberOfSteps
            !store current configuration
            call saveCoordinates(mol, xi)

            !propose y = xi + r*q (q in [-1,1])
            call randomNewVector(mol, r)

            !store trial configuration from the new vector
            call saveCoordinates(mol, y)

            !compute deltaE 
            Ey = TotalEnergy(mol)
            deltaE = Ey-Ei

            !check the acceptance rule
            if (deltaE <= 0.0_KREAL) then
                Ei = Ey !always accept 
                acceptedMoves = acceptedMoves+1

            else 
                !accept with probability exp(-deltaE/kT)
                call random_number(randomNumber)

                if (randomNumber<exp(-deltaE/(kB*temperature))) then 
                    Ei = Ey
                    acceptedMoves = acceptedMoves + 1
                else 
                    !reject & store previous config
                    call restoreCoordinates(mol, xi)
                end if

            end if

            !track best config (lowest energy)
            if (Ei <bestEnergy) then 
                bestEnergy = Ei
                call saveCoordinates(mol, xBest)
            end if

        end do

        !restore the best coordinates found
        call restoreCoordinates(mol, xBest) 

        !calculate acceptance rate
        acceptanceRate = real(acceptedMoves, KREAL)/real(NumberOfSteps, KREAL)

        !saving storage 
        deallocate(xi)
        deallocate(y)
        deallocate(xBest)

    end subroutine


    !routine to save molecule coordinates into an array
    subroutine saveCoordinates(mol, coordinates)
        type(Molecule), intent(in) :: mol
        real(KREAL), intent(out) :: coordinates(:,:)

        integer :: i 

        do i = 1, size(mol%atoms)
            coordinates(i,1) = mol%atoms(i)%x
            coordinates(i,2) = mol%atoms(i)%y
            coordinates(i,3) = mol%atoms(i)%z
        end do

    end subroutine

    !routine to restore coordinates (reverse of save coordinates)
    subroutine restoreCoordinates(mol, coordinates)
        type(Molecule), intent(inout) :: mol
        real(KREAL), intent(in) :: coordinates(:,:)

        integer :: i 

        do i = 1, size(mol%atoms)
            mol%atoms(i)%x = coordinates(i,1)
            mol%atoms(i)%y = coordinates(i,2)
            mol%atoms(i)%z = coordinates(i,3)
        end do

    end subroutine

    !subroutine for calculating y with random q*r
    subroutine randomNewVector(mol, r)
        type(Molecule), intent(inout) :: mol
        real(KREAL), intent(in) :: r

        integer :: i 
        real(KREAL) :: qx,qy,qz !q vector elements

        do i=1, size(mol%atoms)
            !random numbers for vector q in [0,1) with a built in function
            call random_number(qx)
            call random_number(qy)
            call random_number(qz)

            !convert to range -1 to 1
            qx = 2.0_KREAL * qx - 1.0_KREAL
            qy = 2.0_KREAL * qy - 1.0_KREAL
            qz = 2.0_KREAL * qz - 1.0_KREAL

            !save 
            mol%atoms(i)%x = mol%atoms(i)%x + r*qx
            mol%atoms(i)%y = mol%atoms(i)%y + r*qy
            mol%atoms(i)%z = mol%atoms(i)%z + r*qz

        end do

    end subroutine

end module

