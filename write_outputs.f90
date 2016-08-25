module write_outputs
  use pdb

  implicit none

! subroutine write_coordinate(r, unitnum, N)
! subroutine write_energy(e, unitnum, i)

contains

  subroutine write_coordinate(r, unitnum, N)
    implicit none

    double precision,allocatable,dimension(:,:) , intent(in) :: r
    integer , intent(in) :: unitnum
    integer , intent(in) :: N

    integer i

    type (PDB_atom), allocatable, dimension(:) :: PA

    allocate(PA(N))

    do i=1,N,1
       PA(i)%nameRecord="ATOM"
       PA(i)%serial=i
       PA(i)%atm_name="Ar"
       PA(i)%alte_l_i=""
       PA(i)%res_name=""              
       PA(i)%chain_i=""              
       PA(i)%res_seq_num=1
       PA(i)%code_for_insertion_res=""
       PA(i)%x=r(i,1)
       PA(i)%y=r(i,2)
       PA(i)%z=r(i,3)
       PA(i)%occupancy=0.0
       PA(i)%temperature_factor=0.0
       PA(i)%seq_i=""
       PA(i)%element=""   
    end do

    call write_PDB(unitnum, PA, N)

    deallocate(PA)
  end subroutine write_coordinate

  subroutine write_energy(e, unitnum, i)
    implicit none

    double precision , intent(in) :: e
    integer , intent(in) :: unitnum
    integer , intent(in) :: i

    write(unitnum,"(I4,' : 'E10.4)"),i,e
  end subroutine write_energy
end module write_outputs
