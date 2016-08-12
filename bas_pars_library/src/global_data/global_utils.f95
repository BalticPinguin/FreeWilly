!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bas_pars v.0.1 - fork of the AUGER program                          !
!                                                                     !
! Module:  global_utils                                               !
! Version: 0.1                                                        !
! Purpose: contains different utilities.                              !
!          => routines for (de)allocation of global parameter arrays  !
!                                                                     !
! Author:  Gilbert Grell, University of Rostock                       !
! Date:    09.08.2016                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module global_utils
use :: iso_c_binding
IMPLICIT NONE

!
contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 routines for (de)allocation of global parameter arrays                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine global_par_allocate() bind(c)
  use :: global_data
  use :: utils
  IMPLICIT NONE
! status variable
  integer :: stat
! 1.) allocate the BASIS, ATOM, ENV arrays
  allocate(BASIS(1:BASIS_SLOTS, 1:NBAS), STAT=stat)
  if(stat .ne. 0)then
    call error_alloc('BASIS', '1:BASIS_SLOTS, 1:NBAS', 'global_allocate', stat)
  endif
!
  allocate(ATOM(1:ATOM_SLOTS, 1:NATOM), STAT=stat)
  if(stat .ne. 0)then
    call error_alloc('ATOM', '1:ATOM_SLOTS, 1:NATOM', 'global_allocate', stat)
  endif
!
  allocate(ENV(1:(PTR_ENV_START + 3*NATOM + ENV_BASIS_DIM)), STAT=stat)
  if(stat .ne. 0)then
    call error_alloc('ENV', '1:(PTR_ENV_START + 3*NATOM + ENV_BASIS_DIM)', 'global_allocate', stat)
  endif
!
! 1.1.) allocate arrays containing basis set names and contractions:
  allocate(BASIS_NAME(1:NATOMIC_BASES), STAT=stat)
  if(stat .ne. 0)then
    call error_alloc('BASIS_NAME', '1:NATOMIC_BASES', 'global_allocate', stat)
  endif
!
  allocate(BASIS_CONTRACTION(1:NATOMIC_BASES), STAT=stat)
  if(stat .ne. 0)then
    call error_alloc('BASIS_CONTRACTION', '1:NATOMIC_BASES, 1:NBAS', 'global_allocate', stat)
  endif
end subroutine global_par_allocate
!
! deallocate the global parameter arrays
subroutine global_par_deallocate() bind(c)
  use :: global_data
  use :: utils
  IMPLICIT NONE
  integer :: stat
! 1.) deallocate BASIS ATOM and ENV arrays
  deallocate(BASIS, STAT=stat)
  if (stat .ne. 0)then
    call error_dealloc('BASIS', 'global_deallocate', stat)
  endif
!
  deallocate(ATOM, STAT=stat)
  if (stat .ne. 0)then
    call error_dealloc('ATOM', 'global_deallocate', stat)
  endif
!
  deallocate(ENV, STAT=stat)
  if (stat .ne. 0)then
    call error_dealloc('ENV', 'global_deallocate', stat)
  endif
!
! 1.1.) allocate arrays containing basis set names and contractions:
  deallocate(BASIS_NAME, STAT=stat)
  if (stat .ne. 0)then
    call error_dealloc('BASIS_NAME', 'global_deallocate', stat)
  endif
!
  deallocate(BASIS_CONTRACTION, STAT=stat)
  if (stat .ne. 0)then
    call error_dealloc('BASIS_CONTRACTION', 'global_deallocate', stat)
  endif
!
end subroutine global_par_deallocate

end module global_utils