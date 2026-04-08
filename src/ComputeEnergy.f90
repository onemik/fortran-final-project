!module for calculating the energies (separate and then total)

module ComputeEnergy
    use NumberKinds
    use Types
    use Geometry
    use Parameters
    implicit none
    private
    public :: ComputeAllEnergies, TotalEnergy

    contains

    !compute total energy
    function TotalEnergy(mol) result(Etotal)
        type(Molecule), intent(in) :: mol
        real(KREAL) :: Etotal
        type(Energies) :: Energy

        Energy = ComputeAllEnergies(mol)
        Etotal = Energy%total

    end function

    !compute energy breakdown of the molecule 
    function ComputeAllEnergies(mol) result(Energy)
        type(Molecule), intent(in) :: mol
        type(Energies) :: Energy

        Energy%stretch = StretchEnergy(mol)
        Energy%bend = BendingEnergy(mol)
        Energy%nonbond = NonBondingEnergy(mol)
        Energy%torsion = TorsionalEnergy(mol)
        Energy%total = Energy%stretch + Energy%bend + Energy%nonbond + Energy%torsion

    end function

    !calculate stretch energy 
    function StretchEnergy(mol) result(Estretch)
        type(Molecule), intent(in) :: mol
        real(KREAL) :: Estretch

        integer :: a
        integer :: i,j 
        real(KREAL) :: r, r0, k

        !starting value 
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
        
            !stretch energy formula
            Estretch = Estretch + k * (r-r0)**2

        end do

    end function

    !calculate bending energy
    function BendingEnergy(mol) result(Ebending)
        type(Molecule), intent(in) :: mol
        real(KREAL) :: Ebending

        integer :: a
        integer :: i,j,k
        real(KREAL) :: theta

        !starting energy
        Ebending = 0.0_KREAL

        !looping
        do a = 1, size(mol%angles)
            i = mol%angles(a)%i
            j = mol%angles(a)%j
            k = mol%angles(a)%k

            !calculate the value of the angle for each angle that exists
            theta = CalculateAngle(mol%atoms(i), mol%atoms(j), mol%atoms(k))

            !calculate the energy with bending energy formula
            Ebending = Ebending + kAngle * (theta - theta0placeholder)**2 !FIXME 

        end do

    end function



    function TorsionalEnergy(mol) result(Etorsional)
        type(Molecule), intent(in) :: mol
        real(KREAL) :: Etorsional

        integer :: t
        integer :: i,j,k,l
        real(KREAL) :: omega_tors

        !starting energy
        Etorsional = 0.0_KREAL

        do t=1, size(mol%torsions)
            i = mol%torsions(t)%i
            j = mol%torsions(t)%j
            k = mol%torsions(t)%k
            l = mol%torsions(t)%l

            omega_tors = CalculateTorsionalAngle(mol%atoms(i), mol%atoms(j), mol%atoms(k), mol%atoms(l))

            !formula from exercise
            Etorsional = Etorsional + 0.5_KREAL * V1_tors * (1 + cos(real(n_tors,KREAL)*omega_tors - gamma_tors))

        end do

    end function



    function NonBondingEnergy(mol) result(Enonbonding)
        type(Molecule), intent(in) :: mol
        real(KREAL) :: Enonbonding

        integer :: pair
        integer :: i,j
        real(KREAL) :: r !distance
        real(KREAL) :: r6, r12 !for the formula 

        !starting energy
        Enonbonding = 0.0_KREAL

        !for every non bonded pair calculating the van der waals energy
        do pair = 1, size(mol%nobonds)
            i = mol%nobonds(pair)%i
            j = mol%nobonds(pair)%j

            r = CalculateDistance(mol%atoms(i), mol%atoms(j))

            r6 = r**6
            r12 = r6**2

            !energy formula
            Enonbonding = Enonbonding + (Aij/r12 - Bij/r6)
        end do

    end function



end module
