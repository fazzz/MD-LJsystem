module pdb
  implicit none

  type PDB_atom
     character nameRecord*6
     integer serial
     character atm_name*4,alte_l_i,res_name*3,chain_i

     integer res_seq_num
     character code_for_insertion_res
     real x, y, z
     real occupancy,temperature_factor
     character seq_i*2, element*2
  end type PDB_atom

! subroutine write_PDB (unitnum, PA, N) <- To write output coordinates in Protein Data Bank (PDB) format.
! subroutine read_PDB (unitnum, PA, N)  <- To read inout coordinate in PDB format.

contains

  subroutine write_PDB (unitnum, PA, N)

    implicit none
    type (PDB_atom),allocatable,dimension(:),intent(in) :: PA
    integer, intent(in) :: unitnum
    integer, intent(in) :: N

    integer i
    
    write(unitnum,'("MODEL")')
    do i = 1,N,1
       write(unitnum,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)'),&
            PA(i)%nameRecord,PA(i)%serial,PA(i)%atm_name,PA(i)%alte_l_i,PA(i)%res_name,&
            PA(i)%chain_i,PA(i)%res_seq_num,PA(i)%code_for_insertion_res,&
            PA(i)%x,PA(i)%y,PA(i)%z,&
            PA(i)%occupancy,PA(i)%temperature_factor,&
            PA(i)%seq_i,PA(i)%element
    end do
    write(unitnum,'("ENDML")')
    
  end subroutine write_PDB

  subroutine read_PDB (unitnum, PA, N)
    implicit none
    integer, intent(in) :: unitnum
    type (PDB_atom),allocatable,dimension(:),intent(inout) :: PA
    integer, intent(in) :: N
    
    character MOD,EMD
    integer i    

    read(unitnum,'(A5)'),MOD
    do i = 1,N,1
       read(unitnum,'(6X,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)'),&
            PA(i)%serial,PA(i)%atm_name,PA(i)%alte_l_i,PA(i)%res_name,PA(i)%chain_i,&
            PA(i)%res_seq_num,PA(i)%code_for_insertion_res,PA(i)%x,PA(i)%y,PA(i)%z,&
            PA(i)%occupancy,PA(i)%temperature_factor,&
            PA(i)%seq_i,PA(i)%element
    end do
    read(unitnum,'(A5)'),EMD

  end subroutine read_PDB

end module PDB
