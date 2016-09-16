!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bas_pars v.0.1 - fork of the AUGER program                          !
!                                                                     !
! Module:  global_data                                                !
! Version: 0.1                                                        !
! Purpose: provide  geometry, basis set and environmental variables   !
! Author:  Gilbert Grell, University of Rostock                       !
! Date:    04.08.2016                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module global_data
! include the constants module
use :: constants
use, intrinsic :: iso_c_binding
IMPLICIT NONE
!
!? maybe move all type definitions here?
!
! Everything I define here is meant to be public available, so:
public
!
! 1.) - basis set and geometry in libcint format: provide the parameters foraccessing the atom basis and ENV arrays, as well as the arrays themself
! 1.1.) initialize parameters for the ATOM, BASIS and ENV array which are interfaced to libcint
! 1.1.1.) atom array
! slots in the atom array (2 are currently unused)
integer(c_int), parameter :: ATOM_SLOTS = 6
! charge of the atom
integer(c_int), parameter   :: CHARGE_OF = 1
! pointer to coordinates (in bohr) in the ENV array
integer(c_int), parameter :: PTR_COORD = 2
! nuclear model of atom. = 2 corresponds to gaussian nuclear model (normalized GTOs)
integer(c_int), parameter :: NUC_MOD = 3 
! pointer to location of nuclear charge distribution parameter zeta in the ENV array.
integer(c_int), parameter :: PTR_ZETA = 4
!
! 1.1.2.) BASIS array
! slots in the BASIS array
integer(c_int), parameter :: BASIS_SLOTS = 8
! atom index to which the basis belongs
integer(c_int), parameter :: ATOM_OF = 1
! angular momentum of basis
integer(c_int), parameter :: ANG_OF = 2
! number of primitive GTOs in basis
integer(c_int), parameter :: NPRIM_OF = 3
! number of contracted GTOs
integer(c_int), parameter :: NCTR_OF = 4
! kappa parameter, when using spinor basis functions, not used currently.
integer(c_int), parameter :: KAPPA_OF = 5
! offset_exponent  => pointer to ENV(...) array where exponents are stored
integer(c_int), parameter :: PTR_EXP    = 6 
! offset_coefficien => pointer to ENV(...) array where coefficients are stored
integer(c_int), parameter :: PTR_COEFF  = 7 
!
! 1.1.3.) ENV array
! offset of the starting point in the ENV array.
integer(c_int), parameter :: PTR_ENV_START = 20
!
! 1.2.) initialize the actual arrays
! 1.2.1) Atom and basis
integer(c_int), dimension(:,:), allocatable :: ATOM, BASIS
! 1.2.2.) Environmental array which actually holds the data
real(c_double), dimension(:), allocatable :: ENV
!
! 1.3.) Misc basis set information
! dimension of the ENV array needed for exponents and coefficients:
integer(c_int) :: ENV_BASIS_DIM
! number of different atomic bases (number of distinct atoms in the molecule)
integer(c_int) :: NATOMIC_BASES
! basis set name array of length NATOMIC_BASES containing basis set name for each distinct atom
character(len=lbn), dimension(:), allocatable :: BASIS_NAME
! Array of length NATOMIC_BASES containing the contractions for each atomic basis (for each distinct atom)
character(len=lbc), dimension(:), allocatable :: BASIS_CONTRACTION  
! 
!
! NATOM:     number of atoms
! NBAS:      number of bases of different type, ie. s,p,d,etc. basis on 3 atoms make up for NBAS=3*3=9
! NBASF:     number of contracted basis functions
! NBASPRIM:  number of primitive basis functions
integer(c_int), bind(c) :: NATOM, NBAS, NBASF, NBASPRIM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! 3.) Define information that is read from the command line.
! INFILE:      Input file:
character(len=ld) :: INFILE
! NUM_THREADS:    Number of Threads per core to use for openmp/MPI
! NUM_NODES:      Number of cores to use for MPI
integer(c_int) :: NUM_THREADS, NUM_NODES

end module global_data