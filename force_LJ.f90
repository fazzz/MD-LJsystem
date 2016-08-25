! This module calculate Leneord-Jones potential force and energy
! with periodic bounary conditions.
module force_LJ
  implicit none

contains
  double precision function energy_force_LJ(r, L, A, B, C2, Mum, list_p, list_b, f)
    implicit none

    integer , intent(in) :: Mum
    integer , allocatable, dimension(:,:), intent(in) :: list_p, list_b
    double precision, allocatable, dimension(:,:), intent(in) :: r
    double precision, allocatable, dimension(:,:), intent(inout) :: f
    double precision, intent(in) :: L, A, B, C2

    integer i, j, k, x, y, z
    double precision len, r2, r6, r12, r14, dx, dy, dz, ff, fx, fy, fz

    energy_force_LJ = 0.0d0
    f = 0.0d0
    do i=1,Mum,1
       j = list_p(i,1)
       k = list_p(i,2)

       x = list_b(i,1)
       y = list_b(i,2)
       z = list_b(i,3)

       dx = r(j,1) - (r(k,1) + L*x)
       dy = r(j,2) - (r(k,2) + L*y)
       dz = r(j,3) - (r(k,3) + L*z)

       r2 = dx*dx + dy*dy + dz*dz

       if ( r2 < C2 ) then
          r6 = r2 * r2 * r2
          r12 = r6 * r6
          r14 = r12 * r2

          energy_force_LJ = energy_force_LJ + (A - B * r6) / r12

          ff = ( 12.0d0 * A - 6.0d0 * B * r6 ) / r14

          fx = ff * dx
          fy = ff * dy
          fz = ff * dz

          f(j,1) = f(j,1) + fx
          f(j,2) = f(j,2) + fy
          f(j,3) = f(j,3) + fz

          f(k,1) = f(k,1) - fx
          f(k,2) = f(k,2) - fy
          f(k,3) = f(k,3) - fz
       end if
    end do
    energy_force_LJ = energy_force_LJ * 1.0d-2
  end function energy_force_LJ
end module force_LJ
