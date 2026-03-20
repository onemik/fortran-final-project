module InputOutput
    use NumberKinds
    use Types
    implicit none
    private
    public :: ReadFile

    contains

    subroutine ReadFile(mol)
        !define variables
        type(Molecule), intent(out) :: mol

        character(len=200) :: FileName 
        character(len=200) :: line !to read through the file "line by line"
        character(len=2) :: symbol !for C or H 
        integer :: unit, ios
        integer :: numberOfAtoms
        real(KREAL) :: x,y,z
        integer :: i

        !asking user for filename 
        print *, "Please enter the input filename (e.g. c4h10.xyz)"
        read (*,'(A)') FileName 

        !open file 
        unit = 9
        open(unit=unit, file=trim(FileName), status = 'old', action='read', iostat=ios)

            if (ios>0) then
                print*, "Error: could not open the file."
                stop
            end if 

            !reading first number in the file as number of atoms; first read of this file
            read(unit, *, iostat=ios) numberOfAtoms 
            
            if (ios>0) then
                print*, "Error: could not read the number of atoms. Please check the file structure (see documentation)."
                stop
            end if 


            !allocate array of atoms inside of the molecule 
            allocate(mol%atoms(numberOfAtoms))

            !read each atom line 
            i = 0
            do while (i<numberOfAtoms)
                read(unit, '(A)', iostat = ios) line !read the whole line as string

                !file errors check to not crash the code
                if (ios<0) then
                    print*, "Error: reached end of file before reading all atoms."
                    stop
                else if (ios>0) then
                    print*, "Error:problem reading file."
                    stop
                end if

                !skip if it's an empty line
                if (len_trim(line)==0) cycle 

                !read the line and split into parts for each element
                read(line, *, iostat=ios) symbol, x, y, z
                if (ios/=0) then
                    print*, "Error: invalid atom line:"
                    print*, trim(line)
                    stop
                end if

                !store data into the molecule array

                i = i+1
                mol%atoms(i)%symbol = trim(symbol) !atom 
                mol%atoms(i)%x = x
                mol%atoms(i)%y = y
                mol%atoms(i)%z = z

            end do

        close(unit)

    end subroutine



end module
