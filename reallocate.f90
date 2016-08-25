! This module provides the subroutine reallocate 
! which is used to re-size the array.
module reallocate
  implicit none

  contains

    subroutine reallocate_integer(A, newsizex, newsizey)
      implicit none
      integer, allocatable, dimension(:,:), intent(inout) :: A
      integer,  intent(in) :: newsizex,newsizey

      integer, allocatable, dimension(:,:) :: B

      allocate(B(lbound(A,1):ubound(A,1),lbound(A,2):ubound(A,2)))
     
      B = A

      deallocate(A)

      allocate(A(newsizex,newsizey))

      A(lbound(B,1):ubound(B,1),lbound(B,2):ubound(B,2)) = B

      deallocate(B)
    end subroutine reallocate_integer
end module reallocate
