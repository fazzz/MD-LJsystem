module pbs
  implicit none

contains

  subroutine pbs_coordinate(r, Num, L)
    implicit none

    double precision, allocatable,dimension(:,:),intent(inout) :: r
    integer, intent(in) :: Num
    double precision, intent(in) :: L

    integer i, j

    do i=1,Num,1
       do j=1,3,1
          if ( r(i,j) > L) then
             r(i,j) = r(i,j) - L
          else if (r(i,j) < 0.0d0) then
             r(i,j) = r(i,j) + L
          end if
       end do
    end do
  end subroutine pbs_coordinate

end module pbs
