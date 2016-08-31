!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bas_pars v.0.1 - fork of the AUGER program                          !
!                                                                     !
! Module:  parseInp                                                   !
! Version: 0.1                                                        !
! Purpose: molecular geometry & basis set information                 !
! Author:  Gilbert Grell, University of Rostock                       !
! Date:    08.08.2016                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Entities to be filled with data (in the module global_data):
! BASIS:     array contiaining basis set information for every atom (2D array, 8 slots for each basis)
! ATOM:      2D array containing atom information (charge, pointer to coordinates in env, etc. 6 SLots for each atom)
! ENV:       environmental array containing geometric and basis set information
!
! NATOM:     number of atoms
! NBAS:      number of bases of different type, ie. s,p,d,etc. basis on 3 atoms make up for NBAS=3*3=9
! NBASF:     number of basis functions
!
module parse_inp
IMPLICIT NONE


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  Parsing routines to read DATA from the input- and datafiles              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! get_parameters: aquires the data necessary to allocate the global datas and fill them in a second step.
!                 for me details refer to the different parts below
subroutine get_parameters(filehandle,filename)
  use :: constants, only : lk, ll, lv
  use :: global_data, only : NATOMIC_BASES, NATOM 
  use :: parse_utils, only : count_lines, early_term, illegal_value
  use :: basis_utils, only : distinct_atoms, get_basis_parameters
  use :: utils, only : lc2uc, readnocomment, reached_EOF, read_error, set_cursor, fopen, fclose
  IMPLICIT NONE
! specify variables of the subroutine
! filehandle
  integer, intent(in) :: filehandle
! key is found in stream at position (beginning key, end key, normal key)
  integer :: bfound, efound, found
! misc incrementable variables
  integer :: count
! status of the read statement
  integer :: stat
! charges array (determined in the geometry section)
  integer, dimension(:), allocatable :: charges
! array containing the occurence of charges as value and the charges as index
  integer, dimension(:), allocatable :: chargesOccur
! keys for searching through the input file (begin, end and normal)
  character(len=lk) :: bkey
  character(len=lk) :: ekey
  character(len=lk) :: key
! stream for reading in whole lines
  character(len=ll) :: stream
! string to pass errorneous values as words to error routines:
  character(len=lv) :: errvalue
! dump certain strings
  character(len=ll) :: dummy
! name of the inputfile
  character(len=*) filename

! found at position variable needs to be at 0 for default
  bfound=0
  efound=0
  found=0
! open the input file
! open(filehandle,FILE=filename,STATUS='old',IOSTAT=ios,ERR=999)
  call  fopen(filename,filehandle,'old')
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2.) determine the number of atoms in the molecule  WORKS: yes       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
write(bkey,'(a)')'BEGIN GEOMETRY'
write(ekey,'(a)')'END GEOMETRY'
! 2.1.) determine the number of atoms
! number of atoms is equal to the number of lines between the BEGIN and END GEOMETRY statements 
! exit on error if the block is not found, since it is a mandatory block.
  NATOM = count_lines(filehandle,filename,bkey,ekey,'mandatory',stat)
!  check if the value of NATOM is >0 and throw an error if this is not the case
  if(NATOM .le. 0)then
    write(errvalue,'(I5)')NATOM
    call illegal_value(bkey, 'number of atoms' , errvalue, ' > 0')
  endif
! 2.2.) read in the charges and store them in an intermediate array
  allocate(charges(1:NATOM))
! deallocate charges array:
! set cursor to the geometry block
  call set_cursor(filehandle,filename,bkey,'mandatory',stat)
! intialize the atom count to 0
  count = 0
  efound = 0
  do
    stream = readnocomment(filehandle,stat)
    if(stat .lt. 0)then
      call reached_EOF(filename, bkey)
    elseif(stat .gt. 0)then
      call read_error(filename,bkey,stream,stat)
    endif
    efound = index(lc2uc(stream),trim(ekey))
!   make sure that the line is not the END GEOMETRY STATEMENT and that it is not commented=empty:
    if((efound .eq. 0) .and. (len(trim(stream)) .gt. 0))then
      count = count +1
!     read the charges into the array, dump the rest (3 coordinates and 1 unit specification)
      read(stream,*,IOSTAT=stat)charges(count),dummy,dummy,dummy,dummy
      if(stat .ne. 0)then
        call read_error("'streamed line from input file'",key,stream,stat)
      endif
!   if not all atoms were found, but the end block statement is there: premature termination
    elseif((efound .gt. 0) .and. (count .lt. NATOM))then
     efound = 0
      call early_term(bkey,"'atomic charges'")
!   if all atoms were found and the end block statement is there: exit the loop
    elseif((efound .gt. 0) .and. (count .eq. NATOM))then
     efound = 0
      exit
!   if too many atoms were found and the end block is there: abort
    elseif((efound .gt. 0) .and. (count .gt. NATOM))then
     efound = 0
      call early_term(bkey,"'atomic charges'")
    endif
  enddo
! 2.3.) determine : the number of distinct atoms from the charges array (equal charges: corresponding atoms)
!                   the chargesOccur array that contains the number of occurences of an atom with a certain charge in the molecule
! allocate the chargesOccur array. it must have the maximum charge in the molecule as highest dimension:
  allocate(chargesOccur(1:maxval(charges(:))))
! initialize the array to 0
  chargesOccur = 0
! the number of distint atoms is equal to the number of atomic basis to be specified in the input file
  NATOMIC_BASES  = distinct_atoms(charges,chargesOccur)
! deallocate the charges array:
  deallocate(charges)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 9.)  Get all basis set related parameters: (NBAS, NBASF, NATOMIC_BASIS, ENV_BASIS_DIM). WORKS: yes !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  call get_basis_parameters(filehandle, filename, chargesOccur)
!
! deallocate the chargesOccur array. it only served as a variable to calculate the number of angular bases (NBAS)
  deallocate(chargesOccur)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End of the parameter aquisation !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  close(filehandle)
  call fclose(filename, filehandle)
! exit with an error if the input file can not be opened correctly.
!999 continue
!  if(ios.NE.0)then
!    WRITE(*,'(A,A,A,I4)')'Error opening file: ', filename,'IOSTAT = ', ios
!    error stop 'Error upon openening the input file'
!  endif
!  return
end subroutine get_parameters
!
!
!
subroutine get_par_arrays(filehandle,filename)
  use :: constants, only : lk, ll, lu, BOHR2ANGS
  use :: global_data, only : NATOM, CHARGE_OF, PTR_COORD, ATOM, PTR_ENV_START, ENV
  use :: parse_utils, only : early_term, get_state_list
  use :: basis_utils, only : get_basis, order_basis
  use :: utils, only : lc2uc, readnocomment, reached_EOF, read_error, set_cursor, fopen, fclose
  IMPLICIT NONE
! specify variables of the subroutine
! filehandle
  integer, intent(in) :: filehandle
! key is found in stream at position (beginning key, end key, normal key)
  integer :: efound
! line counter
! offset for the global ENV array:
  integer :: off
! misc incrementable variables
  integer ::  count
! status of the read statement
  integer :: stat
! keys for searching through the input file (begin, end and normal)
  character(len=lk) :: bkey
  character(len=lk) :: ekey
  character(len=lk) :: key
! unit specifier for the geometry. must be either ANGSTROM or BOHR
  character(len=lu) :: angsorbohr
! stream for reading in whole lines
  character(len=ll) :: stream
! name of the inputfile
  character(len=*) :: filename

! inititalize the found variable to 0
  efound = 0

! open the input file
!  open(filehandle,FILE=filename,STATUS='old',IOSTAT=ios,ERR=999)
  call fopen(filename, filehandle,'old')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1.) determine the geometry of the molecule  WORKS: yes!     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  write(bkey,'(a)')'BEGIN GEOMETRY'
  write(ekey,'(a)')'END GEOMETRY'
! set cursor to the geometry block
  call set_cursor(filehandle,filename,bkey,'mandatory',stat)
! intialize the atom count to 0
  count = 0
! initialize the offset in the ENV array to the starting value
  off = PTR_ENV_START
  do
    stream = readnocomment(filehandle,stat)
    if(stat .lt. 0)then
      call reached_EOF(filename, 'GEOMETRY BLOCK')
    elseif(stat .gt. 0)then
      call read_error(filename,bkey,stream,stat)
    endif
    efound = index(lc2uc(stream),trim(ekey))
!   make sure that the line is not the END GEOMETRY STATEMENT and that it is not commented=empty:
    if((efound .eq. 0) .and. (len(trim(stream)) .gt. 0))then
      count = count + 1
!     read geometry and charges into the ENV and ATOM array, respectively.
      read(stream,*,IOSTAT=stat)ATOM(CHARGE_OF,count),ENV(off + 1),ENV(off + 2),ENV(off + 3),angsorbohr
      if(stat .ne. 0)then
        call read_error("'streamed line from input file'",key,stream,stat)
      endif
!     define pointer to the coordinates:
      ATOM(PTR_COORD,count) = off
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! possibly add here the nuclear model and zeta value if needed      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(trim(lc2uc(adjustl(angsorbohr))).eq.'ANGSTROM')then
!       convert the geometry to bohr, if present in Angstrom.
        ENV(off+1) = ENV(off+1)/BOHR2ANGS
        ENV(off+2) = ENV(off+2)/BOHR2ANGS
        ENV(off+3) = ENV(off+3)/BOHR2ANGS
      elseif(trim(lc2uc(adjustl(angsorbohr))).eq.'BOHR')then
!       do nothing, since the geometry is already specified in bohr
      else
!       if neither Angstrom nor Bohr is specified: throw an error and abort
        write(*,"(80('*'))")
        write(*,'(a)')'ERROR while reading the geometry.'
        write(*,'(a)')'the units are not specified as "Angstrom" or "Bohr"'
        write(*,"(80('*'))")
        write(*,'(a)')'program stops now ...'
        error stop 'Erroneous unit specification in the geometry section.'
      endif
!     increment the offset for the ENV array by 3 for the written coordinates
      off = off + 3
!   if not all atoms were found, but the end block statement is there: premature termination
    elseif((efound .gt. 0) .and. (count .lt. NATOM))then
     efound = 0
      call early_term(bkey,"'geometry'")
!   if all atoms were found and the end block statement is there: exit the loop
    elseif((efound .gt. 0) .and. (count .eq. NATOM))then
     efound = 0
      exit
!   if too many atoms were found and the end block is there: abort
    elseif((efound .gt. 0) .and. (count .gt. NATOM))then
     efound = 0
      call early_term(bkey,"'geometry'")
    endif
  enddo
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 9.)  Get all basis set related arrays (ENV array: Exponents and coefficients). WORKS: yes! !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call get_basis(filehandle, filename)
! 9.1.) sort the basis array such that it has the same order of the atoms in the geometry and
!       angular momentum from lowest to highest. 
  call order_basis()
! end the subroutine:
  call fclose(filename, filehandle)
end subroutine get_par_arrays
! 
subroutine get_do_array(filehandle,filename, dyorb)
   use, intrinsic :: iso_c_binding
   use :: constants, only : lk, ll, dp
   !  use :: global_data, only: ATOM, BASIS, ENV, ATOM_SLOTS, BASIS_SLOTS, NATOM, NBAS,PTR_ENV_START, ENV_BASIS_DIM
   use :: global_data, only: NBASF
   use :: parse_utils, only : early_term, get_state_list
   use :: basis_utils, only : get_basis, order_basis
   use :: utils, only : lc2uc, readnocomment, reached_EOF, read_error, set_cursor, fopen, fclose
   IMPLICIT NONE
   ! specify variables of the subroutine
   ! filehandle
   integer, intent(in) :: filehandle
   ! key is found in stream at position (beginning key, end key, normal key)
   integer :: efound, enerfound, atomfound, normfound
   ! misc incrementable variables
   integer ::  count
   ! status of the read statement
   integer :: stat
   real(kind=dp) :: energy ! TODO
   real(kind=dp) :: norm
   ! keys for searching through the input file (begin, end and normal)
   character(len=lk) :: bkey
   character(len=lk) :: ekey
   character(len=lk) :: key
   character(len=lk) :: enerkey
   character(len=lk) :: atomkey
   character(len=lk) :: normkey
   ! stream for reading in whole lines
   character(len=ll) :: stream
   ! name of the inputfile
   character(len=*) :: filename
   character(len=ll) :: dump1, dump2, dump3, dump4, dump5
   ! up-spin and down-spin dyson orbitals
   real(c_double), intent(out), dimension(0:NBASF-1) :: dyorb
   real(kind=dp), dimension(:), allocatable :: dyDO
   real(kind=dp), dimension(:), allocatable :: dyUP
   
   ! inititalize the found variable to 0
   efound = 0
   enerfound=0
   atomfound=0
   normfound=0
   
   ! initialise up-spin and down-spin orbitals:
   !  integer(c_int), intent(out), dimension(1:BASIS_SLOTS, 1:NBAS) :: dyor
   allocate(dyUP(0 :NBASF-1))
   allocate(dyDO(0 :NBASF-1))

! open the input file
!  open(filehandle,FILE=filename,STATUS='old',IOSTAT=ios,ERR=999)
   call fopen(filename, filehandle,'old')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1.) read dyson orbitals of the molecule  WORKS: no!         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   write(bkey,'(a)')'BEGIN DO'
   write(ekey,'(a)')'END DO'
   write(enerkey,'(a)')'ISTATE->FSTATE:'
   write(atomkey,'(a)')'ATOM '
   write(normkey,'(a)')'SUMNORM'
   ! set cursor to the geometry block
   call set_cursor(filehandle,filename,bkey,'mandatory',stat)
   ! intialize the atom count to 0
   count = 0
      
   do
      stream = readnocomment(filehandle,stat)
      if(stat .lt. 0)then
         call reached_EOF(filename, 'DO BLOCK')
      elseif(stat .gt. 0)then
         call read_error(filename,bkey,stream,stat)
      endif
      efound = index(lc2uc(stream),trim(ekey))
      enerfound = index(lc2uc(stream),trim(enerkey))
      atomfound = index(lc2uc(stream),trim(atomkey))
      normfound = index(lc2uc(stream),trim(normkey))
      ! make sure that the line is not the END DO statement and that it is not commented=empty:
      if((efound .eq. 0) .and. (len(trim(stream)) .gt. 0))then
         ! read geometry and charges into the ENV and ATOM array, respectively.
         if(enerfound .ne. 0) then ! energy is found
            ! read the energy
            read(stream,*,IOSTAT=stat)dump1,dump2,dump3,dump4,dump5,energy
         elseif (atomfound .ne. 0) then !the atom statement is found --> ignore it.
            ! don't read data from that line, just skip it.
            read(stream,*,IOSTAT=stat)
         elseif (normfound .ne. 0) then !the norm of DO is found 
            read(stream,*,IOSTAT=stat) dump1, dump2, norm
         else  ! this line contains the dyson orbital:
            read(stream,*,IOSTAT=stat) dump1, dump2, dyUP(count), dyDO(count)
            if(stat .ne. 0)then
               call read_error("'streamed line from input file'",key,stream,stat)
            endif
            count=count+1
            !if( dump1 .not. what I expect) then
            !  call read_error("'Fuck you!'",key, stream, stat)
            !endif
            !if( dump2 .not. what I expect) then
            !  call read_error("'Fuck you!'",key, stream, stat)
            !endif
          endif !enerfound .ne. 0
      !if not all atoms were found, but the end block statement is there: premature termination
      elseif((efound .gt. 0) .and. (count .lt. NBASF))then
         efound = 0
         call early_term(bkey,"'dyson orbitals'")
      !if all atoms were found and the end block statement is there: exit the loop
      elseif((efound .gt. 0) .and. (count .eq. NBASF))then
         efound = 0
         exit
      !if too many atoms were found and the end block is there: abort
      elseif((efound .gt. 0) .and. (count .gt. NBASF))then
         efound = 0
         call early_term(bkey,"'dyson orbitals'")
      endif
   enddo
   do count =0, NBASF-1, 1
      dyorb(count)=dyUP(count)+dyDO(count)
   enddo
   !
   ! end the subroutine:
   call fclose(filename, filehandle)
   deallocate(dyUP)
   deallocate(dyDO)
end subroutine get_do_array
! 
!
! subroutine that takes an inputfile name and returns the dimesnions of the ATOM, basis and env array:
subroutine pass_parameters (infile, length, natm, nbs, nbsf, env_bs_dim, ptr_env_st) bind (c)
  use, intrinsic :: iso_c_binding
  use :: global_data, only: NATOM, NBAS, NBASF, PTR_ENV_START, ENV_BASIS_DIM
  IMPLICIT NONE
! misc incrementors:
  integer :: i
! input file name coming from c
  character(c_char), intent(in), dimension(1:length) :: infile
! length of the c_char array infile:
  integer(c_int), intent(in) :: length
! input file name used in the fortran routines
  character(len=length) :: filename
! returned parameters:
  integer(c_int), intent(out) :: natm, nbs, nbsf, env_bs_dim, ptr_env_st
  
!  write(*,'(a,i2)')'length = ', length
  do i = 1, length, 1
    filename(i:i) = infile(i)
!    write(*,'(a)')filename(i:i) enddo
  enddo
 
! open the passed file and get the environmental parameters:
  call get_parameters(20, filename)
! pass back the now assigned global variables:
  natm       = NATOM
  nbs        = NBAS
  nbsf       = NBASF
  env_bs_dim = ENV_BASIS_DIM
  ptr_env_st = PTR_ENV_START
end subroutine pass_parameters


! function that takes an inputfile name and returns the arrays BASIS, ATOM and ENV.
subroutine pass_arrays(infile, length, atm, bas, envi) bind(C)
  use, intrinsic :: iso_c_binding
  use :: global_data, only: ATOM, BASIS, ENV, ATOM_SLOTS, BASIS_SLOTS, NATOM, NBAS, PTR_ENV_START, ENV_BASIS_DIM
  use :: global_utils, only : global_par_allocate, global_par_deallocate
  IMPLICIT NONE
! status variable
  integer(c_int) :: stat
! misc incrementors:
  integer :: i
! input file name coming from c
  character(c_char), intent(in), dimension(1:length) :: infile
! length of the c_char array infile:
  integer(c_int), intent(in) :: length
! input file name used in the fortran routines
  character(len=length) :: filename
! arrays that shall be passed back to the c routine:
  integer(c_int), intent(out), dimension(1:ATOM_SLOTS, 1:NATOM) :: atm
  integer(c_int), intent(out), dimension(1:BASIS_SLOTS, 1:NBAS) :: bas  
  real(c_double), intent(out), dimension(1:PTR_ENV_START + 3*NATOM + ENV_BASIS_DIM) :: envi
  stat = 0 
!  write(*,'(a,i2)')'length = ', length
  do i = 1, length, 1
    filename(i:i) = infile(i)
!    write(*,'(a)')filename(i:i)
  enddo
 
! open the passed file and get the environmental parameters:
!  call get_parameters(20, filename)
! allocate global arrays (ATOM, BASIS, ENV)
  call global_par_allocate()
! get the global parameter arrays:
  call get_par_arrays(20,filename)
  atm  = ATOM
  bas  = BASIS
  envi = ENV
! deallocate global arrays (ATOM, BASIS, ENV)
  call global_par_deallocate()
end subroutine pass_arrays

! function that takes an inputfile name and returns the array DYOR (HUBERT)
subroutine pass_dyor(infile, length, dyor) bind(C)
   use, intrinsic :: iso_c_binding
   ! later: check which of them are used  (HUBERT)
   !use::global_data, only: ATOM, BASIS, ENV, ATOM_SLOTS, BASIS_SLOTS, NATOM, NBAS, PTR_ENV_START, ENV_BASIS_DIM
   use :: global_data, only:  NBASF
   use :: global_utils, only : global_par_allocate, global_par_deallocate
   IMPLICIT NONE
   ! status variable
   integer(c_int) :: stat
   ! misc incrementors:
   integer :: i
   ! input file name coming from c
   character(c_char), intent(in), dimension(1:length) :: infile
   ! length of the c_char array infile:
   integer(c_int), intent(in) :: length
   ! input file name used in the fortran routines
   character(len=length) :: filename
   ! arrays that shall be passed back to the c routine:
   real(c_double), intent(out), dimension(0:NBASF-1) :: dyor
   stat = 0 
   !  write(*,'(a,i2)')'length = ', length
   do i = 1, length, 1
    filename(i:i) = infile(i)
   !    write(*,'(a)')filename(i:i)
   enddo 
   ! get the dyson orbital array:
   call get_do_array(20,filename, dyor)
end subroutine pass_dyor


end module parse_inp
