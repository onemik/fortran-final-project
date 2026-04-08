module Types
  use NumberKinds
  implicit none
  private
  public :: Atom, Bond, Molecule, Angle, Torsion, NoBond, Energies

  type :: Atom
    character(len=2) :: symbol !C or H
    real(KREAL) :: x,y,z !coordinates
  end type

  type :: Bond
    integer :: i,j !to describe the bond between two atoms 
    character(len=2) :: kind !either CC or CH
  end type

  type :: NoBond
    integer :: i,j !for not bonded pairs
  end type
        
  !angle term to describe which atoms form an angle
  type :: Angle
    integer :: i,j,k !three atoms
  end type

  !torsion term to describe which atoms are twisted
  type :: Torsion
    integer :: i,j,k,l !four atoms
  end type
        
  !molecule containing all sub-types
  type :: Molecule
    type(Atom), allocatable :: atoms(:) 
    type(Bond), allocatable :: bonds(:)
    type(Angle), allocatable :: angles(:)
    type(Torsion), allocatable :: torsions(:)
    type(NoBond), allocatable :: nobonds(:)
  end type

  type :: Energies
    real(KREAL) :: stretch = 0.0_KREAL
    real(KREAL) :: bend = 0.0_KREAL
    real(KREAL) :: nonbond = 0.0_KREAL
    real(KREAL) :: torsion = 0.0_KREAL
    real(KREAL) :: total = 0.0_KREAL
  end type

end module
