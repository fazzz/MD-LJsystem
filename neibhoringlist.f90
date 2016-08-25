! This module makes neibhoring list
! with perodic boundary condition.
module neibhoringlist
  use reallocate
  implicit none

contains
  integer function make_neibhoring_list(list_pairs, list_mbox, S2, r, Num, L)
    use reallocate
    implicit none

    integer, allocatable,dimension(:,:),intent(inout) :: list_pairs,list_mbox
    integer, intent(in) :: Num
    double precision, intent(in) :: S2, L
    double precision, allocatable, dimension(:,:), intent(in) :: r

    integer  M
    double precision r2, minr2, dx, dy, dz
    integer i, j, k, x, y, z, minx, miny, minz

    M = 1
    do i=1,Num,1
       do j=i+1,Num,1

          minr2 = S2
          do x=-1,1,1
             do y=-1,1,1
                do z=-1,1,1
                   dx = r(i,1) - (r(j,1) + L*x)
                   dy = r(i,2) - (r(j,2) + L*y)
                   dz = r(i,3) - (r(j,3) + L*z)

                   r2 = dx*dx + dy*dy + dz*dz
                   if ( minr2 > r2 ) then
                      minr2 = r2
                      minx = x
                      miny = y
                      minz = z
                   end if
                end do
             end do
          end do

          if ( minr2 < S2 ) then
             list_pairs(M,1) = i
             list_pairs(M,2) = j
             list_mbox(M,1) = minx
             List_mbox(M,2) = miny
             list_mbox(M,3) = minz
             M = M + 1
             call reallocate_integer(list_pairs,M,2)
             call reallocate_integer(list_mbox,M,3)
          end if
       end do
    end do

    make_neibhoring_list = M - 1
  end function make_neibhoring_list
end module neibhoringlist
