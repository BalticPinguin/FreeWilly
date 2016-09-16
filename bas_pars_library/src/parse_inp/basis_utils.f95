!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bas_pars v.0.1 - fork of the AUGER program                          !
!                                                                     !
! Module:  basis_utils                                                !
! Version: 0.1                                                        !
! Purpose: read in the basis set information into the global ENV      !
!          and BASIS arrays and store the basis set names globally    !
!          Supply all routines needed for that                        !
!                                                                     !
! Author:  Gilbert Grell, University of Rostock                       !
! Date:    08.08.2016                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module basis_utils
use, intrinsic :: iso_c_binding
IMPLICIT NONE
contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                basis set utilitites                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Returns the angular momentum quantum number that corresponds to a certain subshell label
! Throws an error, if a label other then s,p,d,f,g,h is requested.
function ang_mom(label) result(angmom) bind(C)
  use :: utils, only : lc2uc
  use, intrinsic :: iso_c_binding
  IMPLICIT NONE
  integer(c_int) :: angmom
  character(c_char), intent(in) :: label
  select case (lc2uc(label))
    case('S')
      angmom = 0
    case('P')
      angmom = 1
    case('D')
      angmom = 2
    case('F')
      angmom = 3
    case('G')
      angmom = 4
    case('H')
      angmom = 5
    case default
      write(*,"(80('*'))")
      write(*,'(a,a,a)')"Error in function 'ang(label)'. label = ",trim(label)," was illegally specified."
      write(*,'(a)')'label must correspond to either of the atomic subshells s,p,d,f,g,h, with ang. momenta l=0,1,2,3,4,5.'
      write(*,"(80('*'))")
      write(*,'(a)')'Program stops'
      error stop 'Wrong or unsupported subshell label =/= s,p,d,f,g,h was specified'
  end select
end function ang_mom
!
! Returns the angular subshell label that corresponds to a certain angular momentum number
! Throws an error, if a ang. mom. other then 0,1,2,3,4,5 is requested.
function ang_label(angmom) result(label) bind(C)
  use :: utils, only : lc2uc
  use, intrinsic :: iso_c_binding
  IMPLICIT NONE
  integer(c_int), intent(in) :: angmom
  character(c_char) :: label
  select case (angmom)
    case(0)
      label = 'S'
    case(1)
      label = 'P'
    case(2)
      label = 'D'
    case(3)
      label = 'F'
    case(4)
      label = 'G'
    case(5)
      label = 'H'
    case default
      write(*,"(80('*'))")
      write(*,'(a,i2,a)')"Error in function 'ang_label'. l = ",angmom," was illegally given as input."
      write(*,'(a)')'the given l must correspond to either of the atomic subshells s,p,d,f,g,h, with ang. momenta l=0,1,2,3,4,5.'
      write(*,"(80('*'))")
      write(*,'(a)')'Program stops'
      error stop 'Wrong or unsupported angular momentum l =/= 0,1,2,3,4,5 was specified'
  end select
end function ang_label
!

! Determine the number and corresponding charges of distinct atoms in the ATOM(CHARGE_OF,:) array.
function distinct_atoms(charges, chargesOccur) result(natomkind)
  use :: global_data, only : NATOM
  IMPLICIT NONE
  integer :: natomkind
  integer :: i
  integer, intent(in), dimension(1:NATOM) :: charges
! Array containing the charges as indexes and the occurence of that charge in the molecule as values
  integer, intent(out), dimension(1:maxval(charges(:))) :: chargesOccur
!
! Initialize the occurence array to zero. -> done in calling routine get_parameters
!  chargesOccur = 0
  natomkind = 0
  do i=1,NATOM,1
    chargesOccur(charges(i)) = chargesOccur(charges(i)) + 1
  enddo
! Count the number of nonzero elements in the chargesOccur array. It is the number of different atoms.
! Throw an error, if the number is smaller then zero.
  do i=1,maxval(charges(:)),1
    if(chargesOccur(i) .gt. 0)then
      natomkind = natomkind + 1
    elseif(chargesOccur(i) .lt. 0)then
      write(*,"(80('*'))")
      write(*,'(a,i3,a,i5,a)')'illegal value of chargesOccur(',i,') = ',chargesOccur(i), 'occured.'
      write(*,'(a)')'chargesOccur(i) contains the number of atoms with the charge i. Thus it has to be >0'
      write(*,'(a)')'Happened while estimating the number of different atoms in the molecule.'
      write(*,'(a)')'Thus counting the occurence of different atom charges in the molecule.'
      write(*,"(80('*'))")
      write(*,'(a)')'Program stops...'
      error stop 'illegal value in the chargesOccurence array'
    endif
  enddo
  return 
 end function distinct_atoms
!
! New function: within the BASIS block count the number of BEGIN & END statements that enclose each distinct basis.
! check also that the atom label of each basis is really unique (no double definitions)
function atomic_bases(filehandle,filename,startKey,stopKey) result(nAtomicBasis)
  use :: constants, only : lk, ll, lat
  use :: utils, only : set_cursor, readnocomment, lc2uc, reached_EOF, read_error
  IMPLICIT NONE
  integer :: nAtomicBasis, bfound, efound, stat, stopfound, count, bcount, ecount, i, j, idx
  integer, intent(in) :: filehandle
  character(len=*), intent(in) :: filename
  character(len=lk) :: bkey, ekey
  character(len=ll) :: dump
  character(len=lat), dimension(:), allocatable :: atomLabel
  character(len=lat) :: bAtomLabel, eAtomLabel
  character(len=*), intent(in) :: startKey, stopKey
  character(len=ll) :: stream
  write(bkey,'(a)')'BEGIN'
  write(ekey,'(a)')'END'
  bfound = 0
  efound = 0
  stopfound = 0
  bcount = 0
  ecount = 0
! 1.) determine the number of distinct atoms
! set cursor to the BEGIN BASIS statement:
  call set_cursor(filehandle,filename,startkey,'mandatory',stat)
  do
    stream = readnocomment(filehandle,stat)
    if(stat .lt. 0)then
      call reached_EOF(filename, bkey)
    elseif(stat .gt. 0)then
      call read_error(filename,bkey,stream,stat)
    endif
    bfound = index(lc2uc(stream), trim(bkey))
    efound = index(lc2uc(stream), trim(ekey))
    stopfound = index(lc2uc(stream), trim(stopKey))
    if(bfound .ne. 0)then
      bfound = 0
      bcount = bcount + 1
!     read the atom label from this line.
      read(stream,*,IOSTAT=stat)dump,bAtomLabel,dump,dump,dump,dump,dump
      if(stat .ne. 0)then
        call read_error("'streamed line from input file'",bkey,stream,stat)
      endif
    elseif(efound .ne. 0)then
      efound = 0
!     make sure that not the END BASIS statement was found.
      if(stopfound .eq. 0)then
        ecount = ecount + 1
        read(stream,*,IOSTAT=stat)dump,eAtomLabel
        if(stat .ne. 0)then
          call read_error("'streamed line from input file'",ekey,stream,stat)
        endif
      endif
!    Test if both atom labels coincide as they must
     if((bcount .eq. ecount) .and. (bcount .gt. 0))then
       if(trim(lc2uc(bAtomLabel)) .ne. trim(lc2uc(eAtomLabel)))then
         write(*,"(80('*'))")
         write(*,'(a)')'Error when reading in basis set information.'
         write(*,'(a)')"Atom labels of 'BEGIN' and 'END' statements do not coincide!!"
         write(*,'(a,a,a,a,a)')'found BEGIN ',trim(adjustl(bAtomLabel)),' and END ',trim(adjustl(eAtomLabel)),' .'
         write(*,"(80('*'))")
         write(*,'(a)')'Program stops...'
         error stop 'Wrong atomic basis set definition in the BASIS block'
       endif
     endif
! if we found the end of the basis set section, check if the obtained values are reasonable
    endif
    if(stopfound .ne. 0)then
      stopfound = 0
!     check if bcount and ecount have the same value (same number of Begin and End statements)
      if(bcount .ne. ecount)then
        write(*,"(80('*'))")
        write(*,'(a)')'Error when reading in basis set information.'
        write(*,'(a)')"Number of 'BEGIN' and 'END' statements does not coincide!!"
        write(*,'(a)')'Each atomic basis set is started and terminated with BEGIN <atom> <ang.mom>'
        write(*,'(a)')'and an END statement, respective.'
        write(*,"(80('*'))")
        write(*,'(a)')'Program stops...'
        error stop 'Wrong atomic basis set definition in the BASIS block'
      else
        nAtomicBasis = ecount
      endif
      exit
    endif
  enddo
! 2.) Check if all atomic bases are defined uniquely, i.e. that no atom is doubled.
  allocate(atomLabel(1:nAtomicBasis))
  call set_cursor(filehandle,filename,startkey,'mandatory',stat)
  count = 0
  do
    stream = readnocomment(filehandle,stat)
    if(stat .lt. 0)then
      call reached_EOF(filename, bkey)
    elseif(stat .gt. 0)then
      call read_error(filename,bkey,stream,stat)
    endif
    efound = index(lc2uc(stream), trim(ekey))
    stopfound = index(lc2uc(stream), trim(stopKey))
    if(efound .ne. 0)then
      efound = 0
      if(stopfound .eq. 0)then
        count = count + 1
        read(stream,*,IOSTAT=stat)dump,atomLabel(count)
        if(stat .ne. 0)then
          call read_error("'streamed line from input file'",ekey,stream,stat)
        endif
      else
        stopfound = 0
        exit
      endif
    endif
  enddo
! check if all atoms are defined uniquely:
! count the occurences of atom labels:
  do i = 1, nAtomicBasis - 1, 1
    idx = 0
    do j = i+1, nAtomicBasis, 1
      if(lc2uc(atomLabel(i)) .eq. lc2uc(atomLabel(j)))then
        write(*,"(80('*'))")
        write(*,'(a)')'Error when reading in basis set information.'
        write(*,'(a,a,a)')"atomic basis for atom: '", trim(atomlabel(j))," 'is defined multiple times."
        write(*,'(a)')'Only one atomic basis set per distinct atom is allowed!'
        write(*,"(80('*'))")
        write(*,'(a)')'Program stops...'
        error stop 'Wrong atomic basis set definition in the BASIS block'
      endif
    enddo
  enddo
  deallocate(atomLabel)
  return
end function atomic_bases
!
!
! subroutine that renormalizes the contraction coefficients in the ENV array to mitigate rounding errors in the basis set library.
! works only with real coefficients! no complex arithmetic
! WORKS: yes!
subroutine norm_contr_coeff() bind(c)
  use :: constants, only : dp
  use :: global_data, only : NBAS, ANG_OF, NCTR_OF, NPRIM_OF, PTR_EXP, PTR_COEFF, BASIS, CHARGE_OF, ATOM_OF, ATOM, ENV
  IMPLICIT NONE
! misc incremental indices
  integer :: i,j,k,l
! norm of a certain contraction:
  real(kind=dp) :: norm
! helper array containing the occurence of distinct charges (atoms) in the molecule.
  integer, dimension(:), allocatable :: chargesOccur
! 
! function from libCint
! calculate the inverse norm = normalization constant of a spheric GTO:
! currently only an implicit interface is used. change this to explicit asap.
!
! 1.) multiply the coefficients with the respective normalization constant. 
! allocate the chargeoccurence array:
  allocate(chargesOccur(1:maxval(Atom(CHARGE_OF,:))) )
  chargesOccur = 0
  do i = 1, NBAS, 1
!   increment the occurence number of a certain charge: (only once per atom)
    if(i .gt. 1 )then
      if(BASIS(ATOM_OF,i) .ne. BASIS(ATOM_OF,i-1))then
        chargesOccur(ATOM(CHARGE_OF, BASIS(ATOM_OF,i) + 1)) = chargesOccur(ATOM(CHARGE_OF, BASIS(ATOM_OF,i) + 1)) + 1 
      endif
    elseif(i.eq.1)then
      chargesOccur(ATOM(CHARGE_OF, BASIS(ATOM_OF,i) + 1)) = chargesOccur(ATOM(CHARGE_OF, BASIS(ATOM_OF,i) + 1)) + 1 
    endif
!   only do the normalization once for every angular subbasis of each atomic basis set
    if( chargesOccur(ATOM(CHARGE_OF, BASIS(ATOM_OF,i) + 1)) .eq. 1 )then
      do j = 1, BASIS(NCTR_OF,i), 1
!       first multiply each coefficicient that corresponds to the exp. a_k with the normalization factor N_l(a_k)
        do k = 1, BASIS(NPRIM_OF,i), 1
          ENV(BASIS(PTR_COEFF,i) + (j-1)*BASIS(NPRIM_OF,i) + k)&
       &= ENV(BASIS(PTR_COEFF,i) + (j-1)*BASIS(NPRIM_OF,i) + k)*CINTgto_norm(BASIS(ANG_OF,i), ENV(BASIS(PTR_EXP,i) + k))
        enddo
!       second: renormalize the basis set coefficients for the jth contraction to cure rounding errors in the basis set definition:
!       initialize the norm to 0.0_dp to ensure double precision
        norm = 0.0_dp
        do k = 1, BASIS(NPRIM_OF,i), 1
!         increment the norm of the jth contraction by the coefficients multiplied by 1/N_l((a_k+a_l)/2)^2
!         (result of the integral over the pair of two GTO primitives (gauss. * solid harmonic))
          do l = 1, BASIS(NPRIM_OF,i), 1
            norm = norm + ENV(BASIS(PTR_COEFF,i) + (j-1)*BASIS(NPRIM_OF,i) + k)&
                      & * ENV(BASIS(PTR_COEFF,i) + (j-1)*BASIS(NPRIM_OF,i) + l)&
                      & * CINTgto_norm( BASIS(ANG_OF,i), (ENV(BASIS(PTR_EXP,i) + k) + ENV(BASIS(PTR_EXP,i) + l))/2.0_dp )**(-2.0_dp)
          enddo
        enddo
!       loop one more time over the coefficients of the jth contraction and normalize them.
!       write(*,'(a,1x,g18.12)')'norm = ', sqrt(norm)
        do k = 1, BASIS(NPRIM_OF,i), 1
          ENV(BASIS(PTR_COEFF,i) + (j-1)*BASIS(NPRIM_OF,i) + k)&
       &= ENV(BASIS(PTR_COEFF,i) + (j-1)*BASIS(NPRIM_OF,i) + k) / sqrt(norm)
        enddo

      enddo
    endif
  enddo
deallocate(chargesOccur)


end subroutine norm_contr_coeff

! c-version of the normalization routine (c interoperable)
subroutine c_norm_contr_coeff(ATOM, BASIS, ENV) bind(c)
  use :: constants, only : dp
  use :: global_data, only : NBAS, ANG_OF, NCTR_OF, NPRIM_OF, PTR_EXP, PTR_COEFF, CHARGE_OF, ATOM_OF, BASIS_SLOTS, ATOM_SLOTS,&
                            &PTR_ENV_START, ENV_BASIS_DIM, NATOM
  IMPLICIT NONE
! misc incremental indices
  integer :: i,j,k,l
! norm of a certain contraction:
  real(kind=dp) :: norm
! helper array containing the occurence of distinct charges (atoms) in the molecule.
  integer, dimension(:), allocatable :: chargesOccur
! ATOM, BASIS, ENV arrays as passed from C:
  integer(c_int), dimension(1:BASIS_SLOTS, 1:NBAS) :: BASIS
  integer(c_int), dimension(1:ATOM_SLOTS, 1:NATOM) :: ATOM
  real(c_double), dimension(1: PTR_ENV_START + 3*NATOM + ENV_BASIS_DIM) :: ENV
! 
! function from libCint
! calculate the inverse norm = normalization constant of a spheric GTO:
! currently only an implicit interface is used. change this to explicit asap.
!
! 1.) multiply the coefficients with the respective normalization constant. 
! allocate the chargeoccurence array:
  allocate(chargesOccur(1:maxval(Atom(CHARGE_OF,:))) )
  chargesOccur = 0
  do i = 1, NBAS, 1
!   increment the occurence number of a certain charge: (only once per atom)
    if(i .gt. 1 )then
      if(BASIS(ATOM_OF,i) .ne. BASIS(ATOM_OF,i-1))then
        chargesOccur(ATOM(CHARGE_OF, BASIS(ATOM_OF,i) + 1)) = chargesOccur(ATOM(CHARGE_OF, BASIS(ATOM_OF,i) + 1)) + 1 
      endif
    elseif(i.eq.1)then
      chargesOccur(ATOM(CHARGE_OF, BASIS(ATOM_OF,i) + 1)) = chargesOccur(ATOM(CHARGE_OF, BASIS(ATOM_OF,i) + 1)) + 1 
    endif
!   only do the normalization once for every angular subbasis of each atomic basis set
    if( chargesOccur(ATOM(CHARGE_OF, BASIS(ATOM_OF,i) + 1)) .eq. 1 )then
      do j = 1, BASIS(NCTR_OF,i), 1
!       first multiply each coefficicient that corresponds to the exp. a_k with the normalization factor N_l(a_k)
        do k = 1, BASIS(NPRIM_OF,i), 1
          ENV(BASIS(PTR_COEFF,i) + (j-1)*BASIS(NPRIM_OF,i) + k)&
       &= ENV(BASIS(PTR_COEFF,i) + (j-1)*BASIS(NPRIM_OF,i) + k)*CINTgto_norm(BASIS(ANG_OF,i), ENV(BASIS(PTR_EXP,i) + k))
        enddo
!       second: renormalize the basis set coefficients for the jth contraction to cure rounding errors in the basis set definition:
!       initialize the norm to 0.0_dp to ensure double precision
        norm = 0.0_dp
        do k = 1, BASIS(NPRIM_OF,i), 1
!         increment the norm of the jth contraction by the coefficients multiplied by 1/N_l((a_k+a_l)/2)^2
!         (result of the integral over the pair of two GTO primitives (gauss. * solid harmonic))
          do l = 1, BASIS(NPRIM_OF,i), 1
            norm = norm + ENV(BASIS(PTR_COEFF,i) + (j-1)*BASIS(NPRIM_OF,i) + k)&
                      & * ENV(BASIS(PTR_COEFF,i) + (j-1)*BASIS(NPRIM_OF,i) + l)&
                      & * CINTgto_norm( BASIS(ANG_OF,i), (ENV(BASIS(PTR_EXP,i) + k) + ENV(BASIS(PTR_EXP,i) + l))/2.0_dp )**(-2.0_dp)
          enddo
        enddo
!       loop one more time over the coefficients of the jth contraction and normalize them.
!       write(*,'(a,1x,g18.12)')'norm = ', sqrt(norm)
        do k = 1, BASIS(NPRIM_OF,i), 1
          ENV(BASIS(PTR_COEFF,i) + (j-1)*BASIS(NPRIM_OF,i) + k)&
       &= ENV(BASIS(PTR_COEFF,i) + (j-1)*BASIS(NPRIM_OF,i) + k) / sqrt(norm)
        enddo

      enddo
    endif
  enddo
deallocate(chargesOccur)


end subroutine c_norm_contr_coeff


!
!
! subroutine to extract all basis set information and pass it to the environmental arrays BASIS and ENV.
! WORKS: yes!
!
subroutine get_basis(filehandle, filename)
  use :: constants, only : lbc, lbn, lat, lal, lab, ll, lk, lv, dp
  use :: global_data, only : NATOM, CHARGE_OF, ATOM, NATOMIC_BASES, BASIS_NAME, BASIS_CONTRACTION, &
                           & ATOM_OF, ANG_OF, NPRIM_OF, NCTR_OF, PTR_EXP, PTR_COEFF, BASIS, PTR_ENV_START, ENV
  use :: parse_utils, only : split_num_char, remove_chars, separator2space, early_term
  use :: utils, only : set_cursor, readnocomment, lc2uc, reached_EOF, read_error
  IMPLICIT NONE
! filehandle:    Filehandle of the input file (or any other file containing the basis set.)
  integer, intent(in) :: filehandle
! stat:          misc status 
! charge:        atomic charge
  integer ::     stat, charge
! off:           generic offset variable
! expOff:        offset for the exponents written to the ENV array
! coeffOff:      offset for the coefficients written to the ENV array
! efound:        integer indicating that an END statement 'ekey' was found
! found:         integer indicating that a generic keyword 'key' was found
  integer :: off, expOff, coeffOff, efound, found
! iterators:     i, j, k, l, n, count, bas
  integer :: i, j, k, l, n, count, bas
! nPrim:         1D integer array (1:angCount) containing the number of primitives for each ang. mom.
! nContr:        1D integer array (1:angCount) containing the number of contractions for each ang. mom.
! corrIDx:       1D integer array (1:corrAtom) containing the indices of corresponding atoms (atoms with a certain charge)
! corrAtom:      1D integer array (1:NATOMIC_BASES) containing the # of corresponding atoms, i.e. atoms with a certain charge.
! primAngCount:  1D integer array (1:NATOMIC_BASES) containing the # of ang. mom. contributions to the primitive basis
! contrAngCount: 1D integer array (1:NATOMIC_BASES) containing the # of ang. mom. contributions to the contracted basis
  integer, dimension(:), allocatable :: nPrim, nContr, corrIdx, corrAtom, primAngCount, contrAngCount
! filename:      Name of the file which is opened for reading at filhandle
  character(len=*), intent(in) :: filename
! primitive:     string to read the primitive specs from the input file
! contraction:   string to read the contraction specs from the input file
  character(len=lbc) :: primitive, contraction
! basisset:      string that contains the basisset name -> global?
  character(len=lbn) :: basisset
! atomLabel:     string containing the label of an Atom
! coeffAtomLabel:    string contianing the atom label read from the exponent/coeff. block
! termAtomLabel: string containing the atom label read from the atomic basis termination END <atomlable>
  character(len=lat) :: atomLabel, coeffAtomLabel, termAtomLabel
! coeffAngLabel:       string containing the ang. mom. for a specific exponent/coeff. block
  character(len=lal) :: coeffAngLabel
! stream:        streamed line from the inputfile
  character(len=ll) :: stream
! key:           keyword to search for within the block
! bkey:          keyword denoting the beginning of a block
! ekey:          keyword denoting the end of a block
  character(len=lk) :: key, bkey, ekey
! angBas:        1D char array (1:angCount) containing angular mom. subset specs (4s, 3p, etc.) for each ang. mom.
  character(len=lab), dimension(:), allocatable :: angBas
! primAng:       1D char array (1:angCount) containing the ang. momenta present in the primitive string
! contrAng:      1D char array (1:angCount) containing the ang. momenta present in the contraction string
  character(len=lal), dimension(:), allocatable :: primAngLabel, contrAngLabel
! 
! dummy:         long dummy string:
  character(len=ll) :: dummy
! errvalue:      string to pass errorneous values as words to error routines:
  character(len=lv) :: errvalue
!
! exponent:      1D real array of double size (1: #primitives) containing the exponents for each primitive
!                with certain ang. mom. for every distinct atom (each atomic basis)
  real(kind=dp), dimension(:), allocatable :: exponent
! contrCoeff:    2D real array of double precision (1: #primitives, 1: #contracted gaussians) containing the contraction 
!                coefficients for each contracted gaussian, for certain ang. mom of a certain atomic basis
!                (for a certain distinct atomtype)
  real(kind=dp), dimension(:,:), allocatable :: contrCoeff
! 9.2.) Read the different basis sets for different atoms. Loop over the atom kinds.
! Put the cursor to the BEGIN BASIS line and exit on error, if the line is not found, since the block is mandatory.
  write(bkey,'(a)')'BEGIN BASIS'
  write(ekey,'(a)')'END BASIS'
  call set_cursor(filehandle,filename,bkey,'mandatory',stat)
! cursor sits now at BEGIN BASIS
! define the offset for the ENV array as the starting offset + geometry xyz 
  off = PTR_ENV_START + NATOM*3
! allocate the array containing the # of atoms that correpsond to each atomicbasis:
  allocate(corrAtom(1:NATOMIC_BASES))
! allocate the arrays contianing the number of primitive and contracted angular momenta (must coincide) for each atomic basis 
  allocate(primAngCount(1:NATOMIC_BASES))
  allocate(contrAngCount(1:NATOMIC_BASES))
  do i=1,NATOMIC_BASES,1
    ! search for the next begin statment
    write(key,'(a)')'BEGIN'
    do
      stream = readnocomment(filehandle,stat)
      if(stat .lt. 0)then
        call reached_EOF(filename, key)
      elseif(stat .gt. 0)then
        call read_error(filename,key,stream,stat)
      endif
      found = index(lc2uc(stream),trim(key))
      efound = index(lc2uc(stream),trim(ekey))
!     Read in the atom type, charge, primitive basis and contraction. Needs to be specified in a format like:
!     <BEGIN> <atomlabel> <charge> <basisset> <primtive> < -> > <contracted>. 
!     Where primitive should be a comma-seperated list with no spaces in braces. Contraction should be given in square brackets.
!     like: BEGIN H 1 ano-rcc (8s,4p,3d,1f) -> [6s,4p,3d,1f]
      if(found .ne. 0)then
        found = 0
        read(stream,*,IOSTAT=stat)dummy,atomLabel,charge,basisset,primitive,dummy,contraction
        if(stat .ne. 0)then
          call read_error("'streamed line from input file'",key,stream,stat)
        endif
!       Save the name and contraction for the basis set in the global BASIS_NAME, BASIS_CONTRACTION arrays
        write(BASIS_NAME(i),'(a)')trim(basisset)
        write(BASIS_CONTRACTION(i),'(a,a,a)')trim(primitive),' --> ',trim(contraction)
!       Split the primitive and contraction arrays into the angular momentum contributions.
!       Determine the number of primitive and contracted functions for each ang. momentum basis.
!
!       9.2.1.) primitive
!
!       9.2.1.1.) Check that the primitive is given in the right format: done in get_basis_parameters()
        primitive = adjustl(primitive)
!       9.2.1.2.) Remove the braces and exchange the commas by spaces in the primitive string
        primitive = remove_chars(primitive, '()')
        primitive = adjustl(primitive)
        primitive = separator2space(primitive, ',')
!       9.2.1.3.) Determine the number of ang. momenta in the current basis set.
!       Thus count the number of spaces in the primitive array.
        primAngCount(i)=1
        do j = 1, len(trim(primitive))
          if(primitive(j:j).eq.' ')then
            primAngCount(i) = primAngCount(i) + 1
          endif
        enddo
!       allocate the number and ang. mom. label array for primitive gaussians
        allocate(nPrim(1 : primAngCount(i)))
        allocate(primAngLabel(1 : primAngCount(i)))
!       allocate the transtient angular basis array
        allocate(angBas(1: primAngCount(i)))
!       9.2.1.4.) Read in the primitive and ang. mom. specs.  for all specified ang. mom.
        read(primitive,*,IOSTAT=stat)(angBas(j),j = 1,primAngCount(i))
        if(stat .ne. 0)then
          call read_error("primitive string",key,primitive,stat)
        endif
!       9.2.1.5.) Determine the angular momentum  and # of GTOs in each ang. mom. substring.
        do j = 1, primAngCount(i), 1
          call split_num_char(angBas(j), primAngLabel(j), nPrim(j))
        enddo
!       deallocate the angular basis array (it transient)
        deallocate(angBas)
!       
!       9.2.2.) Contractions
!
!       9.2.2.1.) Check that the contraction is given in the right format: done in get_basis_paramters()
        contraction = adjustl(contraction)
!       9.2.2.2.) Remove the brackets and exchange the commas by spaces in the contraction string
        contraction = remove_chars(contraction, '[]')
        contraction = adjustl(contraction)
        contraction = separator2space(contraction, ',')
!       9.2.2.3.) Determine the number of ang. momenta in the current basis set.
!       Thus count the number of spaces in the contraction array.
        contrAngCount(i)=1
        do j = 1, len(trim(contraction))
          if(contraction(j:j).eq.' ')then
            contrAngCount(i) = contrAngCount(i) + 1
          endif
        enddo
!       If we have a different number of ang. mom. bases in the contracted basis then in the primitive: throw an error.
!       removed here, done in get_basis_parameters()
!       allocate the number and ang. mom. label arrays for contracted gaussians
        allocate(nContr(1 : contrAngCount(i)))
        allocate(contrAngLabel(1 : contrAngCount(i)))
!       allocate the transtient angular basis array
        allocate(angBas(1: contrAngCount(i)))
!       9.2.2.4.) Read in the contraction and ang. mom. specs.  for all specified ang. mom.
        read(contraction,*,IOSTAT=stat)(angBas(j),j = 1,contrAngCount(i))
        if(stat .ne. 0)then
          call read_error("contraction string",key,contraction,stat)
        endif
!       9.2.2.5.) Determine the angular momentum  and # of GTOs in each ang. mom. substring.
        do j = 1, contrAngCount(i), 1
          call split_num_char(angBas(j), contrAngLabel(j), nContr(j))
!         Check if the contraction and primitive angular mom. func. are the same
!         if not: throw an error
!         removed here, done in get_basis_paramters
        enddo
!       deallocate the angular basis array (it transient)
        deallocate(angBas)
        exit
!     Throw an error if the end block statement is reached 
      elseif(efound .ne. 0)then
        efound = 0
        call early_term(bkey,key)
      endif
!   End search for the BEGIN Keyword, if it is found.
    enddo
!
!  9.2.3.) Read in coefficients and exponents of the basis in the NWChem format (dummy.inp)
!  
!   loop over all angular momentum subsets
    do j = 1, contrAngCount(i), 1
!     allocate the transient (for every distinct atom/ atomic basis set) exponent and coefficient arrays:
      allocate(exponent(1:nPrim(j)))
      allocate(contrCoeff(1:nPrim(j), 1:nContr(j)))
!     9.2.3.1.) Search for the header of the ang. mom. subset.
      write(key,'(a)')atomLabel
      do
        stream = readnocomment(filehandle,stat)
        if(stat .lt. 0)then
          call reached_EOF(filename, key)
        elseif(stat .gt. 0)then
          call read_error(filename,key,stream,stat)
        endif
        found = index(lc2uc(stream), lc2uc(trim(key)))
        efound = index(lc2uc(stream), trim(ekey))
        
        if(found .ne. 0)then
          found = 0
!         It is problematic that atomlabels and ang. mom. labels overlap (at least for consistency checks)
!         Split the stream into two chars. First must be the atom label, second must be the angular momentum.
!         Throw an error if this is not the case.
          read(stream,*,IOSTAT=stat)coeffAtomLabel,coeffAngLabel
          if(stat .ne. 0)then
            call read_error("contraction string",key,contraction,stat)
          endif
!
!         Note the intended case sensitivity here. (no case conversion)
          if(coeffAtomLabel .ne. atomLabel)then
            write(*,"(80('*'))")
            write(*,'(a,a)')'ERROR when reading the BASIS BLOCK in the input file: ',trim(filename)
            write(*,'(a,a,a,a)')'Atom: ',trim(coeffAtomLabel), " in the coefficient header: ",trim(stream) 
            write(*,'(a,a)')'does not match the previosly defined Atom = ',trim(atomLabel)
            write(*,'(a,a,a,a,a)')'Happened in basis ',trim(stream),' of the ',trim(atomLabel),' atom'
            write(*,'(a)')'Program stops now ...'
            write(*,"(80('*'))")
            error stop 'Error while reading the basis set from input file'
          endif
!         coeffAngLabel must be exactly one time present in the contrAngLabel string-array.
          count = 0
          do k = 1, contrAngCount(i), 1
            if( trim(adjustl(lc2uc(coeffAngLabel))) .eq. trim(adjustl(lc2uc(contrAngLabel(k)))) )then
              count = count + 1
            endif
          enddo
          if(count .ne. 1)then
            write(*,"(80('*'))")
            write(*,'(a,a)')'ERROR when reading the BASIS BLOCK in the input file: ',trim(filename)
            write(*,'(a,a,a)')'Angular momentum: ',trim(adjustl(coeffAngLabel)), ' is not present in the contraction definition.'
            write(*,'(a,a,a,a,a)')'Happened in basis ',trim(stream),' of the ',trim(atomLabel),' atom'
            write(*,'(a)')'Program stops now ...'
            write(*,"(80('*'))")
            error stop 'Error while reading the basis set from input file'
          endif
!         9.2.3.2.) if everything is fine, proceed with reading in exponents and contraction coefficients:
!         loop over the number of primitives ( = number of exponents)
!          write(*,'(a,i2,a,i2)')'nPrim(',j,') = ',nPrim(j)
          do k = 1, nPrim(j), 1
            stream = readnocomment(filehandle,stat)
            if(stat .lt. 0)then
              call reached_EOF(filename, 'basis set exponents and coefficients')
            elseif(stat .gt. 0)then
              call read_error(filename,'basis set exponents and coefficients',stream,stat)
            endif
!           exclude empty = commented lines from the read:
            if(len(trim(stream)) .gt. 0)then
!             read in one exponent and nContr contraction coefficients per line:
              read(stream,*,IOSTAT=stat)exponent(k),(contrCoeff(k,l), l = 1, nContr(j))
!              write(*,'(a,i2,a,i2)')'nContr(',j,') = ',nContr(j)
              if(stat .ne. 0)then
                call read_error("'streamed line from input file'",'basis set exponents and coefficients',stream,stat)
              endif
            endif
          enddo
!         if the last block of the  ang. mom. subset was read, check for the termination of the atomic basis block:
          if(j .eq. contrAngCount(i))then
            do
              stream = readnocomment(filehandle,stat)
              if(stat .lt. 0)then
                write(errvalue,'(a,a,a)')'Termination of the ',trim(atomlabel),' - basis'
                call reached_EOF(filename, errvalue)
              elseif(stat .gt. 0)then
                call read_error(filename,errvalue,stream,stat)
              endif
!             skip all empty lines betweeen data and END <atomlabel>
              if(len(trim(stream)) .gt. 0)then
                exit
              endif
            enddo
            read(stream,*,IOSTAT=stat)dummy,termAtomLabel
            if(stat .ne. 0)then
              call read_error("'streamed line from input file'",trim(errvalue),stream,stat)
            endif
            if(trim(lc2uc(dummy)) .ne. 'END')then
              write(*,"(80('*'))")
              write(*,'(a,a,a,a)')"ERROR when reading the basis for atom '",trim(atomlabel),"' from file: ",trim(filename)
              write(*,'(a,a,a,a)')"expected termination statement: END '",trim(atomlabel),"' but got: ",trim(stream)
              write(*,'(a)')'Program stops now ...'
              write(*,"(80('*'))")
              error stop 'Error while reading the basis set from input file'
            elseif(trim(termAtomLabel) .ne. coeffAtomLabel)then
              write(*,"(80('*'))")
              write(*,'(a,a,a,a)')"ERROR when reading the basis for atom '",trim(atomlabel),"' from file: ",trim(filename)
              write(*,'(a,a,a,a)')"expected termination statement: END '",trim(atomlabel),"' but got: ",trim(stream)
              write(*,'(a,a,a)')"The atomic label from termination and begin statement: '",trim(termAtomLabel),"' and '",&
                                 &trim(coeffAtomLabel),"' do not coincide!"              
              write(*,'(a)')'Program stops now ...'
              write(*,"(80('*'))")
              error stop 'Error while reading the basis set from input file'
            endif
!           exit the loop after reading exponents and coefficients for all ang. mom subsets
            exit
          endif
!         exit the loop after reading exponents and coefficients for one ang. mom. subset
          exit
        elseif(efound .ne. 0)then
          efound = 0
          call early_term(bkey, key)
        endif
!     end search for the atom.
      enddo
!     9.2.4.) Write the data for all atoms of the current kind (all distinct atoms that share the same basis set)
!     into the BASIS and ENV arrays
!     Thereby do not repeat entries in the ENV array.
!     9.2.4.1.) Determine which atoms in the molecule correspond to the current atomic basis set
!     determine the number of atoms wich have corresponding charges:
      corrAtom(i) = 0
      do k = 1, NATOM, 1
        if(charge .eq. ATOM(CHARGE_OF,k))then
          corrAtom(i) = corrAtom(i) + 1
        endif
      enddo
!     initialize an array containing the indices of the corresponding atoms
      allocate(corrIdx(1:corrAtom(i)))
      l = 0
!     fill the array (NOTE that it has to be 0-based, since the atom indices are 0-based in libcint)
      do k = 1, NATOM, 1
        if(charge .eq. ATOM(CHARGE_OF,k))then
          l = l+1
          corrIdx(l) = k-1
        endif
      enddo
!     9.2.4.2.) Write parameters, exp., coeff. for the distinct atoms into the BASIS and ENV array.
!     define the pointers to exponent and coefficients in the ENV array:
!     exponents and coefficients are stored subsequent to the geometry information
!     exponents are stored starting from the first element after the geometry information
      expOff = off
!     write the exponents into the ENV array. (Note the 1 based iteration. First exponent sits at ENV(expOff + 1).
      do l = 1, nPrim(j), 1
        ENV(expOff + l) = exponent(l)
!        write(*,'(a,i6)')'eidx = ',expOff+l
      enddo
!     increment the total offset as:
!     + offset for the just written exponents
      off = off + nPrim(j)
!
!     coefficients are stored starting from the first element after the exponents:
      coeffOff = off
!     write the contraction coefficients into the ENV array (subseeding the exponents)
!     (Note the 1 based iteration. First coefficient sits at ENV(coeffOff + 1)

      do n = 1, nContr(j), 1
        do l = 1, nPrim(j), 1
          ENV(coeffOff + (n-1)*nPrim(j) + l) = contrCoeff(l,n)
!          write(*,'(a,i6)')'cidx = ',coeffOff+(n-1)*nPrim(j) + l
        enddo
      enddo
!     increment the total offset by the number of coefficients written to the ENV array
!     +  #contr. coefficients
      off = off +  nPrim(j)*nContr(j)
!     loop over corresponding atoms
      do k = 1, corrAtom(i), 1
!       first write the BASIS array (only the atom and basis index depend on k.)
!       calculate the 0-based basis index (each ang. mom. on each atom is a single basis in libcint.)
!       the angular momentum should be the fastest iterator followed by corr. atoms and atomic basis set.
!       bas = ang. mom. index + atom index* number of ang. mom. subsets on this atom
        if(i .eq. 1)then
          bas = j + (k-1)*contrAngCount(i)
        else
          bas = j + (k-1)*contrAngCount(i) + (i-1)*(corrAtom(i-1)*contrAngCount(i-1))
        endif
        BASIS(ATOM_OF,bas) = corrIdx(k) ! 0-based !!
        BASIS(ANG_OF,bas) = ang_mom(trim(contrAngLabel(j)))
        BASIS(NPRIM_OF,bas) = nPrim(j)
        BASIS(NCTR_OF,bas) = nContr(j)
!       set kappa to 0 always or do not define it at all
!       BASIS(KAPPA_OF,bas) = 0
!       store the pointers to exponent and coefficients in the ENV array:
        BASIS(PTR_EXP, bas) = expOff
        BASIS(PTR_COEFF,bas) = coeffOff
!        bas = bas+1
      enddo
!     deallocate the atom index array
      deallocate(corrIdx)
!     deallocate the transient (for every distinct atom & ang. mom. basis) exponent and coefficient arrays:
      deallocate(exponent)
      deallocate(contrCoeff)
!   deallocate  the number and ang. mom. label arrays
!   end loop over the ang. mom. subsets   
    enddo
    deallocate(nPrim)
    deallocate(primAngLabel)
    deallocate(nContr)
    deallocate(contrAngLabel)
!   end loop over the atomic basis sets
  enddo
! deallocate the array containing the # of atoms that correpsond to each bases:
  deallocate(corrAtom)
! deallocate the arrays containing the number of primitive and contracted angular momenta (must coincide) for each atomic basis 
  deallocate(primAngCount)
  deallocate(contrAngCount)
end subroutine get_basis

! subroutine to extract all basis set parameters, esp. the number of bases, basis functions and number of exponents/coefficients.
! ATOM array must be defined already for this routine!
! WORKS: yes
subroutine get_basis_parameters(filehandle, filename, chargesOccur)
  use :: constants, only : lbc, lbn, lat, ll, lk, lab, lal, ll ,lv
  use :: global_data, only : NATOM, NBAS, NBASF, NBASPRIM, NATOMIC_BASES, PTR_ENV_START, ENV_BASIS_DIM
  use :: parse_utils, only : split_num_char, remove_chars, separator2space, early_term, illegal_value
  use :: utils, only: lc2uc, readnocomment, reached_EOF, read_error, set_cursor
  IMPLICIT NONE
! filehandle:    Filehandle of the input file (or any other file containing the basis set.)
  integer, intent(in) :: filehandle
! nDistBases:    number of distinct atomic bases specified in the basis block
! lPrim:         length of the primitives specification
! lContr:        length of the contraction spec
! primAngCount:  number of ang. mom. contributions to the primitive basis
! contrAngCount: number of ang. mom. contributions to the contracted basis
! stat:          misc status 
! charge:        atomic charge
  integer :: nDistBases, lPrim, lContr, primAngCount, contrAngCount, stat, charge
! off:           generic offset variable
! expOff:        offset for the exponents written to the ENV array
! coeffOff:      offset for the coefficients written to the ENV array
! efound:        integer indicating that an END statement 'ekey' was found
! found:         integer indicating that a generic keyword 'key' was found
  integer :: off, efound, found
! iterators:     i, j
  integer :: i, j
! nPrim:         1D integer array (1:angCount) containing the number of primitives for each ang. mom.
! nContr:        1D integer array (1:angCount) containing the number of contractions for each ang. mom.
  integer, dimension(:), allocatable :: nPrim, nContr
! chargesOccur:  1D integer array (1: maxval(charges(:))) containing the occurences of the charges, denoted by the indices
  integer, intent(in), dimension(1:*) :: chargesOccur
! filename:      Name of the file which is opened for reading at filhandle
  character(len=*), intent(in) :: filename
! primitive:     string to read the primitive specs from the input file
! contraction:   string to read the contraction specs from the input file
  character(len=lbc) :: primitive, contraction
! basisset:      string that contains the basisset name -> global?
  character(len=lbn) :: basisset
! atomLabel:     string containing the label of an Atom
  character(len=lat) :: atomLabel
! stream:        streamed line from the inputfile
  character(len=ll) :: stream
! key:           keyword to search for within the block
! bkey:          keyword denoting the beginning of a block
! ekey:          keyword denoting the end of a block
  character(len=lk) :: key, bkey, ekey
! angBas:        1D char array (1:angCount) containing angular mom. subset specs (4s, 3p, etc.) for each ang. mom.
  character(len=lab), dimension(:), allocatable :: angBas
! primAng:       1D char array (1:angCount) containing the ang. momenta present in the primitive string
! contrAng:      1D char array (1:angCount) containing the ang. momenta present in the contraction string
  character(len=lal), dimension(:), allocatable :: primAngLabel, contrAngLabel
! 
! dummy:         long dummy string:
  character(len=ll) :: dummy
! errvalue:      string to pass errorneous values as words to error routines:
  character(len=lv) :: errvalue
!
  write(bkey,'(a)')'BEGIN BASIS'
  write(ekey,'(a)')'END BASIS'
! 9.1.) Determine the number of distinct atomic basis sets specified in the BASIS block
  nDistBases = atomic_bases(filehandle,filename,bkey,ekey) 
! check if the number does coincide with the number of distinct atoms 
! (NATOMIC_BASIS, determined in sect. 2.3 of get_parameters() in the parse_inp module)
  if(NATOMIC_BASES .ne. nDistBases)then
    write(*,"(80('*'))")
    write(*,'(a)')'Error when reading in basis set information.'
    write(*,'(a,i5,a,i5,a)')'NATOMIC_BASES = ',NATOMIC_BASES,' and nDistBases = ',nDistBases,' do not match!'
    write(*,'(a)')'This means that there is one atomic basis set too much or too few specified.'
    write(*,"(80('*'))")
    write(*,'(a)')'Program stops...'
    error stop 'Wrong contraction specification in the mandatory BASIS block'
  endif
! 9.2.) Read the different basis sets for different atoms. Loop over the atom kinds.
! Put the cursor to the BEGIN BASIS line and exit on error, if the line is not found, since the block is mandatory.
  call set_cursor(filehandle,filename,bkey,'mandatory',stat)
! cursor sits now at BEGIN BASIS
! initialize the number of angular bases to 0
  NBAS = 0
! initialize the total number of contracted gaussians to
  NBASF = 0
! initialize the total number of primitives to 0
  NBASPRIM = 0
! initialize the total number of exponents and coefficients to 0
  ENV_BASIS_DIM = 0
! define the offset for the ENV array as the starting offset + geometry xyz 
  off = PTR_ENV_START + NATOM*3
  do i=1,nDistBases,1
    ! search for the next begin statment
    write(key,'(a)')'BEGIN'
    do
      stream = readnocomment(filehandle,stat)
      if(stat .lt. 0)then
        call reached_EOF(filename, key)
      elseif(stat .gt. 0)then
        call read_error(filename,key,stream,stat)
      endif
      found = index(lc2uc(stream),trim(key))
      efound = index(lc2uc(stream),trim(ekey))
!     Read in the atom type, charge, primitive basis and contraction. Needs to be specified in a format like:
!     <BEGIN> <atomlabel> <charge> <basisset> <primtive> < -> > <contracted>. 
!     Where primitive should be a comma-seperated list with no spaces in braces. Contraction should be given in square brackets.
!     like: BEGIN H 1 ano-rcc (8s,4p,3d,1f) -> [6s,4p,3d,1f]
      if(found .ne. 0)then
        found = 0
        read(stream,*,IOSTAT=stat)dummy,atomLabel,charge,basisset,primitive,dummy,contraction
        if(stat .ne. 0)then
          call read_error("'streamed line from input file'",key,stream,stat)
        endif
!       Split the primitive and contraction arrays into the angular momentum contributions.
!       Determine the number of primitive and contracted functions for each ang. momentum basis.
!
!       9.2.1.) primitive
!
!       9.2.1.1.) Check that the primitive is given in the right format:
        primitive = adjustl(primitive)
!       Determine length of primitive string, throw an error if 0
        lPrim = len(trim(primitive))
        if (lPrim .lt. 4)then
          write(errvalue,'(i5)')lPrim
          call illegal_value(bkey,'lengths of primitive string', trim(primitive), '>4 min: (1s)')
        endif
        if((primitive(1:1) .ne. '(') .or. (primitive(lPrim:lPrim) .ne. ')'))then
          write(*,"(80('*'))")
          write(*,'(a)')'Error when reading in basis set information.'
          write(*,'(a,a,a,a)')'Error in line: ',trim(stream),' . Wrong primitive specification: ',trim(primitive)
          write(*,'(a)')"Primtive basis must be specified as a comma - separated list with enclosing parentheses '(', ')'"
          write(*,'(a)')"example: (8s,7p,5d,3f,1g)"
          write(*,"(80('*'))")
          write(*,'(a)')'Program stops...'
          error stop 'wrong primitive specification in the mandatory BASIS block'
        endif
!       9.2.1.2.) Remove the braces and exchange the commas by spaces in the primitive string
        primitive = remove_chars(primitive, '()')
        primitive = adjustl(primitive)
        primitive = separator2space(primitive, ',')
!       write(*,'(a,a)')'primitive = ',trim(primitive)
!       9.2.1.3.) Determine the number of ang. momenta in the current basis set.
!       Thus count the number of spaces in the primitive array.
        primAngCount=1
        do j = 1, len(trim(primitive))
          if(primitive(j:j) .eq. ' ')then
            primAngCount = primAngCount + 1
          endif
        enddo
!        write(*,'(a,i5)')'primangcount = ',primAngCount
!       allocate the number and ang. mom. label array for primitive gaussians
        allocate(nPrim(1 : primAngCount))
        allocate(primAngLabel(1 : primAngCount))
!       allocate the transtient angular basis array
        allocate(angBas(1: primAngCount))
!       9.2.1.4.) Read in the primitive and ang. mom. specs.  for all specified ang. mom.
        read(primitive,*,IOSTAT=stat)(angBas(j),j = 1,primAngCount)
!        write(*,*)(angbas(j),j=1,primAngCount)
        if(stat .ne. 0)then
          call read_error("primitive string",key,primitive,stat)
        endif
!       9.2.1.5.) Determine the angular momentum  and # of GTOs in each ang. mom. substring.
        do j = 1, primAngCount, 1
          call split_num_char(angBas(j), primAngLabel(j), nPrim(j))
        enddo
!       deallocate the angular basis array (it transient)
        deallocate(angBas)
!       
!       9.2.2.) Contractions
!
!       9.2.2.1.) Check that the contraction is given in the right format:
        contraction = adjustl(contraction)
!       Length of contraction string, throw an error if 0
        lContr = len(trim(contraction))
        if (lContr .lt. 4)then
          write(errvalue,'(i5)')lContr
          call illegal_value(bkey,'length of contraction string', trim(errvalue), '>4 min: (1s)')
        endif
        if((contraction(1:1) .ne. '[') .or. (contraction(lContr:lContr) .ne. ']'))then
          write(*,"(80('*'))")
          write(*,'(a)')'Error when reading in basis set information.'
          write(*,'(a,a,a,a)')'Error in line: ',trim(stream),' . Wrong contraction specification: ',trim(contraction)
          write(*,'(a)')"Contracted basis must be specified as a comma - separated list with enclosing brackets '[', ']'"
          write(*,'(a)')"example: [6s,4p,3d,2f,1g]"
          write(*,"(80('*'))")
          write(*,'(a)')'Program stops...'
          error stop 'Wrong contraction specification in the mandatory BASIS block'
        endif
!       9.2.2.2.) Remove the brackets and exchange the commas by spaces in the contraction string
        contraction = remove_chars(contraction, '[]')
        contraction = adjustl(contraction)
        contraction = separator2space(contraction, ',')
!       9.2.2.3.) Determine the number of ang. momenta in the current basis set.
!       Thus count the number of spaces in the contraction array.
        contrAngCount=1
        do j = 1, len(trim(contraction))
          if(contraction(j:j).eq.' ')then
            contrAngCount = contrAngCount + 1
          endif
        enddo
!        write(*,'(a,a)')'contraction =',trim(contraction)
!       If we have a different # of ang. mom. bases in the contracted basis then in the primitive: throw an error.
        if(contrAngCount .ne. primAngCount)then
          write(*,"(80('*'))")
          write(*,'(a)')'Error when reading in basis set information.'
          write(*,'(a,a,a)')"Contracted basis: ",trim(contraction),' does not contain the same number of ang. mom.'
          write(*,'(a,a,a)')'as the primitive: ',trim(primitive),'. This is prohibited. Only the number of functions may vary'
          write(*,"(80('*'))")
          write(*,'(a)')'Program stops...'
          error stop 'Erroneous primitve/contraction definition in the mandatory BASIS block'
        endif
!       allocate the number and angu. mom. label arrays for contracted gaussians
        allocate(nContr(1 : contrAngCount))
        allocate(contrAngLabel(1 : contrAngCount))
!       allocate the transtient angular basis array
        allocate(angBas(1: contrAngCount))
!       9.2.2.4.) Read in the contraction and ang. mom. specs.  for all specified ang. mom.
        read(contraction,*,IOSTAT=stat)(angBas(j),j = 1,contrAngCount)
        if(stat .ne. 0)then
          call read_error("contraction string",key,contraction,stat)
        endif
!       9.2.2.5.) Determine the angular momentum  and # of GTOs in each ang. mom. substring.
        do j = 1, contrAngCount, 1
          call split_num_char(angBas(j), contrAngLabel(j), nContr(j))
!         Check if the contraction and primitive angular mom. func. are the same
!         if not: throw an error
          if(contrAngLabel(j) .ne. primAngLabel(j))then
            write(*,"(80('*'))")
            write(*,'(a)')'Error when reading in basis set information.'
            write(*,'(a,a,a)')"Contracted basis: ",trim(contraction),' does not contain the same ang. mom. parts'
            write(*,'(a,a,a)')'as the primitive: ',trim(primitive),'. This is prohibited. Only the number of functions may vary'
            write(*,"(80('*'))")
            write(*,'(a)')'Program stops...'
            error stop 'Erroneous angular momentum definition in the mandatory BASIS block'
          endif
        enddo
!       increment the total number of contracted and primitive basis functions for every angular momentum
        do j = 1, contrAngCount, 1
          NBASF = NBASF + nContr(j) * ( 2*ang_mom( trim(lc2uc(contrAngLabel(j))) ) + 1 ) * chargesOccur(charge)
          NBASPRIM = NBASPRIM + nPrim(j) * ( 2*ang_mom( trim(lc2uc(primAngLabel(j))) ) + 1 ) * chargesOccur(charge)
        enddo
!       increment the dimension of the ENV array ( = total number of exponents and coefficients on all atoms)
        do j = 1, contrAngCount,1
!          ENV_BASIS_DIM = ENV_BASIS_DIM + (nPrim(j) + nContr(j)*nPrim(j)) * chargesOccur(charge)
          ENV_BASIS_DIM = ENV_BASIS_DIM + (nPrim(j) + nContr(j)*nPrim(j))
        enddo
!       deallocate some arrays
        deallocate(nPrim)
        deallocate(nContr)
        deallocate(primAngLabel)
        deallocate(contrAngLabel)
!       deallocate the angular basis array (it transient)
        deallocate(angBas)
        exit
!     Throw an error if the end block statement is reached 
      elseif(efound .ne. 0)then
        efound = 0
        call early_term(bkey,key)
      endif
!     End search for the BEGIN Keyword, if it is found.
    enddo
!   increment the number of bases by the number of ang.mom. bases in the atomic basis i 
!   NBAS is the number of entries in the BASIS array, there is one set of values for each ang. mom. basis
    NBAS = NBAS + contrAngCount*chargesOccur(charge) 
! end loop over the atomic basis sets
  enddo
!
end subroutine get_basis_parameters
!
!
! subroutine that orders the read basis such that the order of atoms in the geometry is followed.
! The slower quantum number is the angular momentum of the basis.
! Works: no!
subroutine order_basis()
  use :: global_data, only: BASIS, ATOM_OF, NBAS, BASIS_SLOTS, ANG_OF  
  ! incrementors: 
  integer :: i,n
! tmp:           temporary basis array element:
  integer, dimension(1:BASIS_SLOTS) :: tmp
! bubble sort, because NBAS is likely to be small (< 100)
!  write(*,*)'BASIS(ATOM_OF,:) ', (BASIS(ATOM_OF,i), i=1, NBAS)
!  write(*,*)'BASIS(ANG_OF,:) ', (BASIS(ANG_OF,i), i=1, NBAS)

  do n = NBAS, 2, -1
    do i = 1, n-1, 1
      if(BASIS(ATOM_OF, i) .gt. BASIS(ATOM_OF, i+1))then
        tmp(:) = BASIS(:,i+1)
        BASIS(:,i+1) = BASIS(:,i)
        BASIS(:,i) = tmp(:)
      elseif((BASIS(ATOM_OF, i) .eq. BASIS(ATOM_OF, i+1)) .and. (BASIS(ANG_OF, i) .gt. BASIS(ANG_OF, i+1)) )then
        tmp(:) = BASIS(:,i+1)
        BASIS(:,i+1) = BASIS(:,i)
        BASIS(:,i) = tmp(:)
      endif
    enddo 
!    write(*,*)'BASIS(ATOM_OF,:) ', (BASIS(ATOM_OF,i), i=1, NBAS)
!    write(*,*)'BASIS(ANG_OF,:) ', (BASIS(ANG_OF,i), i=1, NBAS)
  enddo
end subroutine

integer function factorial(n)
  implicit none
  integer, intent(in) :: n
  integer :: i,fact
  fact=1
  do i=1, n, 1
     fact=fact*i
  enddo
  factorial=fact
end function factorial

! function to compute the normalisation factor.
! Taken from misc.c of libcint.  HUBERT
real(kind=dp) function CINTgto_norm(n, a) 
  use :: constants, only : dp
  implicit none
  integer, intent(in) :: n
  real(kind=dp), intent(in) :: a
  real(kind=dp) :: nn
  nn= (2**(2*n+3)) * factorial(n+1) * ((2*a)**(n+1.5)) / (factorial(2*n+2) * 1.77245385091)
  CINTgto_norm = sqrt(nn)
end function CINTgto_norm

end module basis_utils
