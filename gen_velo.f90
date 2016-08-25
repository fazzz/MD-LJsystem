! This module generate initial velocity
! of MD simulation.
module gen_velo
  implicit none

contains

  subroutine generate_velocity(v,Num,beta,m)
    use mtmod
    implicit none

    double precision, allocatable,dimension(:,:),intent(inout) :: v
    integer, intent(in) :: Num
    double precision, intent(in) :: beta, m

    integer i,j

    double precision beta_m, dbpi, u1, u2

    dbpi=2.0d0*acos(-1.0d0)

    beta_m = sqrt(m / beta * 1.0d-2)

    do i=1,Num,1
       do j=1,3,1
          u1=grnd()
          u2=grnd()

          v(i,j)=sqrt(-2.0d0*log(u1))*cos(dbpi*u2)*beta_m
       end do
    end do
  end subroutine generate_velocity

end module gen_velo
