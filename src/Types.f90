module com
  implicit none
  public :: ComputeCenterOfMass
  private
        
        
contains

  subroutine ComputeCenterOfMass(filename)
    character(len=*), intent(in) :: filename

    print *, "from test ComputeCenterOfMass subroutine: ", filename 
  end subroutine

end module
