!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bas_pars v.0.1 - fork of the AUGER program                          !
!                                                                     !
! MAIN Program                                                        !
! Version: 0.1                                                        !
! Purpose: Testing                                                    !
! Author:  Gilbert Grell, University of Rostock                       !
! Date:    09.08.2016                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program main
use :: constants
use :: global_data
use :: global_utils
use :: parse_inp
use :: basis_utils
IMPLICIT NONE
! initialize integer variables
! command line arguments:
! number of command line arguments, number of threads (openmp)
integer :: narg
! misc offset
integer :: off
! misc incrementors
integer :: i,j,k
! initialize character variables
!
! command line argument for the number of threads:
character(len=4) :: nthreadarg
! basis set format specifier
character(len=lfm) :: basfmt

!
! start program:
narg=0
narg=command_argument_count()
if(narg.gt.0)then
  if(narg.eq.1)then
!   read the input file name from the command line
    call get_command_argument(1,INFILE)
  elseif(narg.eq.2)then
!   read the input file name and number of OMP threads from the command line
    call get_command_argument(1,INFILE)
    call get_command_argument(2,nthreadarg)
!   write the number of threads to the globa_data variable NUM_THREADs
    read(nthreadarg,*)NUM_THREADS
  else
    write(*,"(80('*'))")
    write(*,'(a)')'Too many command line arguments!'
    write(*,'(a)')'Specify only input file and number of threads for openMP'
    write(*,'(a)')'Program stops...'
    write(*,"(80('*'))")
    error stop
  endif
else
  write(*,"(80('*'))")
  write(*,'(a)')'No input file specified.'
  write(*,'(a)')'pass the input file as a command line argument!'
  write(*,'(a)')'Program stops ...'
  write(*,"(80('*'))")
  error stop
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         Debug section starts here. Test various things           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call get_parameters(20,INFILE)
! check if the current parameters have been set correctly:
! currently works with shuffling around the keywords as well
write(*,'(a,i5)')'NATOM = ',NATOM
write(*,'(a,i5)')'NATOMIC_BASIS = ', NATOMIC_BASES
write(*,'(a,i5)')'NBAS = ', NBAS
write(*,'(a,i5)')'NBASF = ', NBASF
write(*,'(a,i5)')'NBASPRIM = ', NBASPRIM
write(*,'(a,i5)')'ENV_BASIS_DIM = ', ENV_BASIS_DIM
write(*,"(80('*'))")
! allocate global arrays
call global_par_allocate()
! get the parameter arrays, (currently ATOM,BASIS and ENV)
call get_par_arrays(20,INFILE)
!
!
! check if geometry and basis set are read correctly
! print the read geometry:
off = PTR_ENV_START
write(*,"(80('*'))")
write(*,'(a)')'GEOMETRY in ANGSTROM:'
write(*,"(80('*'))")
do i=1, NATOM,1
write(*,'(1x,i4,3(f18.12,2x))')ATOM(CHARGE_OF,i),ENV(ATOM(PTR_COORD,i)+1)*BOHR2ANGS,ENV(ATOM(PTR_COORD,i)+2)*BOHR2ANGS,&
                              &ENV(ATOM(PTR_COORD,i)+3)*BOHR2ANGS
enddo
write(*,"(80('*'))")
! check if the basis sets were read and assigned correctly to the atoms
write(*,"(80('*'))")
write(*,'(a)')'BASIS SETS as read and assigned to atoms:'
write(*,"(80('*'))")
do i=1, NBAS, 1
  if(i .gt. 1)then
    if(BASIS(ATOM_OF,i).ne.BASIS(ATOM_OF,i-1))then
      write(*,*)
     endif
  endif
!                            do not forget the 0 based atom index!
  write(*,'(i3,a,i2,a,i2,a,a)')ATOM(CHARGE_OF,BASIS(ATOM_OF,i)+1), ' (',BASIS(NPRIM_OF,i),' -> ', BASIS(NCTR_OF,i),')',&
                              &ang_label(BASIS(ANG_OF,i))
  do j=1,BASIS(NPRIM_OF,i),1
    write(basfmt,'(a,i2,a)')'(f16.8,',BASIS(NCTR_OF,i),'(f16.8,2x))'
    write(*,basfmt)ENV(BASIS(PTR_EXP,i)+j),(ENV(BASIS(PTR_COEFF,i)+j+(k-1)*BASIS(NPRIM_OF,i)),k=1,BASIS(NCTR_OF,i))
  enddo
enddo
write(*,"(80('*'))")

! 
! multiply the basis set coefficients by the normalization factors for the primitves and renormalize all contractions
call norm_contr_coeff()
! print the atomic overlap matrix if requested

! deallocate global data before finishing the run
call global_par_deallocate()


end program main