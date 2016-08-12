!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bas_pars v.0.1 - fork of the AUGER program                          !
!                                                                     !
! Module:  libcint                                                    !
! Version: 0.1                                                        !
! Purpose: provide explicit interfaces for                            !
!          external functions from the libcint library                !
!                                                                     !
! Author:  Gilbert Grell, University of Rostock                       !
! Date:    08.08.2016                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module libCint
! intrinsic module for c interoperability:
  use, intrinsic :: iso_c_binding
  IMPLICIT NONE

!  contains
! 
! c function: FINT cintcgto_spheric_(const FINT *bas_id, const FINT *bas)
! FINT is probably the datatype of fortran integer, so no adjustment is required here
! => use standard fortran integer!
  interface
    pure integer function CINTcgto_spheric(bas_id, bas) bind(c, name='cintcgto_spheric_')
    use, intrinsic :: iso_c_binding
    use :: global_data
    integer, intent(in) :: bas_id
    integer, intent(in), dimension(1:BASIS_SLOTS,1:NBAS) :: bas
    end function CINTcgto_spheric
  end interface

! calculate the inverse norm = normalization constant of a spheric GTO:
! of ang. mom. l with the exponent a
! c function: double CINTgto_norm_(FINT *n, double *a)
! FINT is probably the datatype of fortran integer, so no adjustment is required here
! => use standard fortran integer!
  interface
    real(kind=c_double) function CINTgto_norm(l, a) bind(c, name='cintgto_norm_')
      use, intrinsic :: iso_c_binding
      real(kind=c_double), intent(in) :: a
      integer, intent(in) :: l
    end function CINTgto_norm
  end interface
!
! calculate the atomic overlap integrals between all subshells of the bases in the shells array.
! store the results in the buffer array and return 0 if the integral is 0.
! c function: FINT cint1e_ovlp_sph(double *opij, const FINT *shls, const FINT *atm, const FINT natm, const FINT *bas, const FINT nbas, const double *env)
  interface
     integer function cint1e_ovlp_sph(buf1e, shls, atm, natm, bas, nbasis, environmental) bind(c, name='cint1e_ovlp_sph_')
      use, intrinsic :: iso_c_binding
      use :: global_data
      integer, intent(in) :: natm, nbasis
      integer, intent(in), dimension(1:2) :: shls
      integer, intent(in), dimension(1:ATOM_SLOTS, 1:NATOM) :: atm
      integer, intent(in), dimension(1:BASIS_SLOTS, 1:NBAS) :: bas
      real(kind=c_double), intent(in), dimension(1:(PTR_ENV_START + 3*NATOM + ENV_BASIS_DIM)) :: environmental
      real(kind=c_double), intent(out), dimension(1:CINTcgto_spheric(shls(1), bas), 1:CINTcgto_spheric(shls(2), bas), 1:1) :: buf1e
    end function cint1e_ovlp_sph
  end interface

end module libCint
