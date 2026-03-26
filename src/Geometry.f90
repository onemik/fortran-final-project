module Geometry
    use NumberKinds
    use Types
    implicit none
    private
    public :: CalculateDistance, CalculateAngle, CalculateTorsionalAngle

    contains

    !distance between 2 atoms (a and b)
    function CalculateDistance(a,b) result(r)
        type(Atom), intent(in) :: a,b
        real(KREAL) :: r

        !using euclidean distance formula
        r = sqrt((a%x-b%x)**2 + (a%y-b%y)**2 + (a%z-b%z)**2) 
    end function

    !angle between 3 atoms i.e. two vectors using dot product
    function CalculateAngle(a,b,c) result(theta)
        type(Atom), intent(in) :: a,b,c
        real(KREAL) :: theta

        real(KREAL) :: vector1(3), vector2(3) !for the two atoms with 3 coords each
        real(KREAL) :: cos_theta, norm1, norm2 !helpers for dot product calulation


        vector1 = [a%x - b%x, a%y - b%y, a%z - b%z] !vector from B to A
        vector2 = [c%x - b%x, c%y - b%y, c%z - b%z] !vector from B to C

        norm1 = sqrt(sum(vector1**2))
        norm2 = sqrt(sum(vector2**2))

        !getting the angle from the arccosine from dot product 
        cos_theta = sum(vector1*vector2)/(norm1*norm2)
        theta = acos(cos_theta)

    end function


    !a function of four atoms A,B,C,D and might be describe as the
    !angle between the A−B and C−D bond, looking down the B−C bond,
    !i.e. a cross product of two vectors with resulting vector perpendicular
    function CalculateTorsionalAngle(a,b,c,d) result(phi)
        type(Atom), intent(in) :: a,b,c,d
        real(KREAL) :: phi

        real(KREAL) :: vector1(3), vector2(3), vector3(3)

        vector1 = [b%x - a%x, b%y - a%y, b%z - a%z] !vector from A to B
        vector2 = [d%x - c%x, d%y - c%y, d%z - c%z] !vector from C to D
        vector3 = [c%x - b%x, c%y - b%y, c%z - b%z] !vector from B to C



    end function

    !using cross product definition
    function CrossProduct(a,b) result(c)
        real(KREAL), intent(in) :: a(3), b(3) !input vectors
        real(KREAL) :: c(3) !output vector of cross product

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    
    end function

end module