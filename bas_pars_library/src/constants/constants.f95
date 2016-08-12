!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bas_pars v.0.1 - fork of the AUGER program                          !
!                                                                     !
! Module:  constants                                                  !
! Version: 0.1                                                        !
! Purpose: globally provide physical & mathematical constants and     !
!          conversion factors, IEEE 754 standards for int and real    !
!                                                                     !
! Author:  Gilbert Grell, University of Rostock                       !
! Date:    04.08.2016                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module constants
! include environmental variables (purpose is to have int32, real32, int64, real64, real128)
use, intrinsic :: iso_fortran_env
use, intrinsic :: iso_c_binding
IMPLICIT NONE

!
! Everything I define here is meant to be public available, so:
public
!
! precision parameters
! integer precisions
! int_32
integer, parameter :: i32 = INT32
! int_64
integer, parameter :: i64 = INT64
! real precisions
! float
integer, parameter :: sp = REAL32
! double
integer, parameter :: dp = REAL64
!integer, parameter :: dp = c_double
! long
integer, parameter :: qp = REAL128
! 
! string length parameters
! reading input file lines 
integer, parameter :: ll = 32*1024
! reading long paths/filenames:
integer, parameter :: ld = 1024
! reading/writing keywords (variables with names *key)
integer, parameter :: lk = 64
! reading basis set contraction specifications like: (7s,6p,4d,2f) -> [5s,4p,2d,1f]
integer, parameter :: lbc = 64
! reading angular momentum basis specifications like 11S or 7F:
integer, parameter :: lab = 8
! reading angular momentum labels lik S and F of 11S and 7F specification, but also 5sp (hybridized basis sets:)
integer, parameter :: lal = 4
! reading basis function labels like 2px or 4d2- or similiar largest: 10h5-
integer, parameter :: lbf = 8
! reading atom labels like Fe,He, etc.
integer, parameter :: lat = 64 
! reading the basis set name like ccPVTZ or ano-rcc
integer, parameter :: lbn = 64
! format strings:
integer, parameter :: lfm = 64
! unit specifier:
integer, parameter :: lu = 64
! interval lentght in the state lists:
integer, parameter :: li = 64
! erroneous / correct values
integer, parameter :: lv = 64
! routine names:
integer, parameter :: lr = 64
! 
! conversion factors
! do not forget the precision flags after the constants! omitting them will default to single precision real!
! 1 bohr in angstroms: taken from NIST CODATA 11.05.2016
real(kind=dp), parameter :: BOHR2ANGS = 0.52917721067_dp
!                                        0.52917721067

! mathematical constants:
! pi with 30 digits:
real(kind=dp), parameter :: PI = 3.141592653589793238462643383279_dp
!

end module constants