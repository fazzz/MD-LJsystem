! This program performes Molecular Dynamics(MD) simulation
! of Lenord-Jones system with periodic boundary conditions.
program MD
  use ini_conf
  use mtmod
  use gen_velo
  use pbs
  use force_LJ
  use reallocate
  use neibhoringlist
  use read_inputs
  use write_outputs

  implicit none

  integer i,j,k

  integer Num, maxsteps, outputinterval, listinterval, seed
  ! Num: The number of atoms in the system.
  ! maxsteps: The number of steps for MD calculatin.
  ! outputinterbal: The frequency of steps to save the coordinate in the output file.
  ! listinterval: The frequenct of steps to re-calculate the neibhoring list.
  ! seed: seed of randum numbers.
  double precision T, epsilon, sigma, L_box_siz, L_cut_off, L_cut_lis, dt, m
  ! T: The temperature of the system.
  ! epsilon: The force constant of the Lenord-Jones potential.
  ! sigma: The collision radius of the atom in the system.
  ! L_box_siz: The length of the box.(A)
  ! L_cut_off: The length of the cutt off.(A)
  ! L_cut_off: The length of the list cut off.(A)
  ! dt: time step width.(ps)
  ! m:  mass of atom.(u)

  integer :: ON=1
  integer :: OFF=0

  double precision, parameter :: k_B = 2.2926d-4

  double precision, allocatable,dimension(:,:) :: r, v, f, ro, vo, fo
  ! r: coordinate, 
  ! v: velocity, 
  ! f: force, 
  ! ro: coordinate(previous step),
  ! vo: velocity(previous step),
  ! fo: force(previous step).
  double precision  dt_m_2, dt2_m_2, m_2
  ! dt_m_2 = dt/(2.0*m),
  ! dt2_m_2 = dt^2/(2.0*m),
  ! m_2 = 1.0/(2.0*m) * 0.01.
  double precision e_pote, e_kine, beta, e_total
  ! e_pote: potential energy,
  ! e_kine: kinetic energy,
  ! beta: inverse temperature,
  ! e_total = e_pote +e_kine.
  double precision A, B
  ! A = epsilon * sigma^12,
  ! B = epsilon * sigma^6.
  double precision  LC2, LS2
  ! LC2 = L_cut_off^2,
  ! LS2 = L_cut_lis^2,

  integer Mum ! M: number of interactions
  integer, allocatable,dimension(:,:) :: N_list, B_list

  integer initializeMODE

  character coordinate_filename*50
  character outputfilename_coordinate*50
  character outputfilename_energy*50

  character(10) :: argv, programname
  integer :: argc, iargc

  ! 0 reading parameters
  NAMELIST /mdcontrols/ Num, maxsteps, outputinterval, listinterval, seed, T, epsilon, sigma, L_box_siz, L_cut_off, L_cut_lis, m, dt
  open(17,file='parameters_md')
  read(17,mdcontrols)
  close(17)                                                                                   

  LC2 = L_cut_off * L_cut_off
  LS2 = L_cut_lis * L_cut_lis
  
  initializeMODE=ON

  argc=iargc()
  call getarg(0,programname)
  if (argc < 2) then
     call usage(programname)
  else if ( argc == 2 ) then
     call getarg(1,outputfilename_coordinate)
     call getarg(2,outputfilename_energy)
  else if ( argc == 3 ) then
     call getarg(1,coordinate_filename)
     call getarg(2,outputfilename_coordinate)
     call getarg(3,outputfilename_energy)
     initializeMODE = OFF
  else
     call usage(programname)
  end if

  allocate(r(Num,3))
  allocate(ro(Num,3))
  allocate(v(Num,3))
  allocate(vo(Num,3))
  allocate(f(Num,3))

  if ( initializeMODE .eq. ON ) then
     call initialize_conformation(ro, Num, L_box_siz)
  else
     open(21,file=coordinate_filename,status='old')
     call read_coordinate(21, ro, Num)
     close(21)
     call pbs_coordinate(ro, Num, L_box_siz)
  end if

  call sgrnd(seed)

  beta = 1.0d0 / (k_B * T) ! kJ/mol/K

  dt_m_2 = dt / m * 0.5d0
  dt2_m_2 = dt * dt / m * 0.5d0
  m_2 = 0.5d0 * m * 1.0d-2 ! kJ/mol/ps^2

  A = epsilon * sigma**12.0d0 / beta ! kJ/mol
  B = epsilon * sigma**6.0d0 / beta  ! kJ/mol

  ! 1 generate initial velocity
  call generate_velocity(v, Num, beta, m)

  ! 2 make initial neiboring list
  allocate(N_list(1,2))
  allocate(B_list(1,3))
  Mum = make_neibhoring_list(N_list, B_list, LS2, ro, Num, L_box_siz)

  e_pote = energy_force_LJ(ro, L_box_siz, A, B, LC2, Mum, N_list, B_list, f)
  e_kine = sum(v * v)
  e_kine = m_2 * e_kine ! kJ/mol

  r = ro + dt * v + dt2_m_2 * f
  call pbs_coordinate(r, Num, L_box_siz)

  open(19,file=outputfilename_coordinate,position='append')
  open(20,file=outputfilename_energy,position='append')

  ro = r
  vo = v
  fo = f

  do i = 1, maxsteps, 1
     ! 3 calculate force
     e_pote = energy_force_LJ(r, L_box_siz, A, B, LC2, Mum, N_list, B_list, f)

     e_kine = sum(v * v)
     e_kine = m_2 * e_kine ! kJ/mol

     ! 4 velocity-Verlet integration
     v = vo + dt_m_2 * (f + fo)
     r = ro + dt * v + dt2_m_2 * f
     call pbs_coordinate(r, Num, L_box_siz)

     ro = r
     vo = v
     fo = f

     ! 5 make neiboring list
     if (mod(i,listinterval) == 0) then
        deallocate(N_list)
        deallocate(B_list)
        allocate(N_list(1,2))
        allocate(B_list(1,3))
        Mum = make_neibhoring_list(N_list, B_list, LS2, r, Num, L_box_siz)
     end if

     ! 6 output results
     if (mod(i,outputinterval) == 0) then
        e_total = e_pote + e_kine
        call write_coordinate(r,19,Num)
        write (20,'(4F13.8)'),i*dt,e_total,e_pote,e_kine
     end if
  end do

  close(19)
  close(20)

  deallocate(r)
  deallocate(ro)
  deallocate(v)
  deallocate(vo)
  deallocate(f)
  deallocate(fo)

  deallocate(N_list)
  deallocate(B_list)

contains
  subroutine usage(programname)
    implicit none
    character(10) :: programname

    print *,"Usage: ",programname," [ inputfilename (PDB file) ] outputfilename1(coordinate) outputfilename2(energy)"
    stop
  end subroutine usage
end program MD
