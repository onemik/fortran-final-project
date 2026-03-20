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


  end type
        
        
  type :: Molecule
    type(Atom), allocatable :: atoms(:) 

  end type

end module
