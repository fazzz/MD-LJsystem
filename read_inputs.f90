! This module contains the subroutine 
! to read the input coordinate file written in Protein Data Bank format.
module read_inputs
  use pdb
  implicit none

contains

  subroutine read_coordinate(unitnum, r, N)
    implicit none

    integer , intent(in) :: unitnum
    double precision, allocatable,dimension(:,:) , intent(inout) :: r
    integer , intent(in) :: N

    integer i

    type (PDB_atom), allocatable, dimension(:) :: PA

    allocate(PA(N))

    call read_PDB(unitnum, PA, N)

    do i=1, N, 1
       r(i, 1) = PA(i)%x
       r(i, 2) = PA(i)%y
       r(i, 3) = PA(i)%z
    end do

    deallocate(PA)
  end subroutine read_coordinate
end module read_inputs
