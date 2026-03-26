module Types
  use NumberKinds
  implicit none
  private
  public :: Atom, Bond, Molecule

  type :: Atom
    character(len=2) :: symbol !C or H
    real(KREAL) :: x,y,z !coordinates
  end type

  type :: Bond
    integer :: i,j !to describe the bond between two atoms 
    character(len=2) :: kind !either CC or CH
  end type
        
  !angle term to tell which atoms form an angle
  type :: Angle
    integer :: i,j,k !three atoms
  end type

  !torsion term to tell which atoms are twisted
  type :: Torsion
    integer :: i,j,k,l !four atoms
  end type
        
  type :: Molecule
    type(Atom), allocatable :: atoms(:) 
    type(Bond), allocatable :: bonds(:)

  end type

end module
