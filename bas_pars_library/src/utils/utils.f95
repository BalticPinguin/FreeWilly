!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bas_pars v.0.1 - fork of the AUGER program                          !
!                                                                     !
! Module:  utils                                                      !
! Version: 0.1                                                        !
! Purpose: contains different utilities that are used in many modules !
!          Currently:                                                 !
!          => contains error routines for the (de)allocation of arrays!
!          => routines for opening and closing files                  !
!                                                                     !
! Author:  Gilbert Grell, University of Rostock                       !
! Date:    05.08.2016                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module utils
IMPLICIT NONE
! List of characters for case conversion
character(len=*), parameter, private :: lower_case = 'abcdefghijklmnopqrstuvwxyz'
character(len=*), parameter, private :: upper_case = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!
contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Error section:  errors during the (de)allocation of an array              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! error during the allocation of an array
subroutine error_alloc(array, boundaries, routine, stat)
  IMPLICIT NONE
  character(len=*), intent(in) :: array
  character(len=*), intent(in) :: boundaries
  character(len=*), intent(in) :: routine
  integer, intent(in) :: stat
  write(*,"(80('*'))")
  write(*,'(a,a,a,a,a,a)')'Error on allocation of array: ', trim(array),'( ',trim(boundaries),' ) in routine: ',trim(routine)
  write(*,'(a,i8,a)')'Status value stat =  ',stat,' was received on allocation.'
  write(*,'(a)')'Program stops now ...'
  write(*,"(80('*'))")
  error stop 'Error while allocating an array'
end subroutine error_alloc
!
! error during the deallocation of an array
subroutine error_dealloc(array, routine, stat)
  IMPLICIT NONE
  character(len=*), intent(in) :: array
  character(len=*), intent(in) :: routine
  integer, intent(in) :: stat
  write(*,"(80('*'))")
  write(*,'(a,a,a,a,a,a)')'Error on deallocation of array: ', trim(array),' in routine: ',trim(routine)
  write(*,'(a,i8,a)')'Status value stat =  ',stat,' was received on deallocation.'
  write(*,'(a)')'Program stops now ...'
  write(*,"(80('*'))")
  error stop 'Error while deallocating an array'
end subroutine error_dealloc

! error upon opening a file:
subroutine fopen_error(filename,filehandle,stat)
  IMPLICIT NONE
  integer, intent(in) :: filehandle
  integer :: stat
  character(len=*), intent(in) :: filename
!
  write(*,"(80('*'))")
  write(*,'(a,1x,a,1x,a,i5)')'ERROR, when trying to OPEN file', trim(filename),'at filehandle unit = ',filehandle
  write(*,'(a,i5,1x,a)')'the open statement returned a nonzero IOSTAT = ', stat, 'value.'
  write(*,'(a,1x,a)')'Check Paths, filename and possibly the existence of the file',trim(filename)
  write(*,'(a)')'Program stops...'
  write(*,"(80('*'))")
  error stop 'Error upon opening a file. nonzero IOSTAT value'
end subroutine fopen_error
!
! error upon closing a file:
subroutine fclose_error(filename,filehandle,stat)
  IMPLICIT NONE
  integer, intent(in) :: filehandle
  integer :: stat
  character(len=*), intent(in) :: filename
!
  write(*,"(80('*'))")
  write(*,'(a,1x,a,1x,a,i5)')'ERROR, when trying to CLOSE file', trim(filename),'at filehandle unit = ',filehandle
  write(*,'(a,i5,1x,a)')'the close statement returned a nonzero IOSTAT = ', stat, 'value.'
  write(*,'(a,1x,a)')'Check Paths, filename and possibly the existence of the file',trim(filename)
  write(*,'(a)')'Program stops...'
  write(*,"(80('*'))")
  error stop 'Error upon closing a file. nonzero IOSTAT value'
end subroutine fclose_error
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   Error section 2.) errors during the read of any file                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! reached EOF statement during reading a file <filename>, looking for a certain keyword <key>
subroutine reached_EOF(filename, key)
  IMPLICIT NONE
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: key
  write(*,"(80('*'))")
  write(*,'(a,a)')'ERROR. Reached EOF during read of the file: ', trim(filename)
  write(*,'(a,a,a)')'Apparently the Keyword: ',trim(key),' is not present.'
  write(*,'(a)')'Program stops now ...'
  write(*,"(80('*'))")
  error stop 'Error while reading file'
end subroutine reached_EOF
!
! error upon a read statement (IOSTAT value =/= 0):
subroutine read_error(filename, key, line, stat)
  IMPLICIT none
  integer :: stat
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: key
  character(len=*), intent(in) :: line
  write(*,"(80('*'))")
  write(*,'(a,a,a,a)')'ERROR during read of line: ',trim(line),' in file: ',trim(filename)
  write(*,'(a,i10,a)')'IOSTAT = ',stat, ' value.'
  write(*,'(a,a)')'Happened while handling the keyword: ',trim(key)
  write(*,'(a)')'Program stops...'
  write(*,"(80('*'))")
  error stop 'Error while reading'
end subroutine read_error
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  File I/O section: routines for opening and closing files                  Work: no!   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! open a file <filename> at position (unit) <filehandle>
subroutine fopen(filename, filehandle,attribute)
  IMPLICIT NONE
  integer, intent(in) :: filehandle
  integer :: ios
  character(len=*), intent(in) :: filename
! attribute of the file, i.e. old, if it is existing, new if it should be created 
  character(len=*), intent(in) :: attribute
! check if attribute is defined correctly
  if((lc2uc(attribute) .ne. 'OLD') .and. (lc2uc(attribute) .ne. 'NEW') .and. (lc2uc(attribute) .ne. 'REPLACE'))then
!   throw an error if other values are used:
    write(*,"(80('*'))")
    write(*,'(a,1x,a,1x,a,1x,a)')'ERROR. When calling fopen to open file', trim(filename),'a wrong attribute =',trim(attribute)
    write(*,'(a)')'was specified.'
    write(*,'(a)')"Only attribute = 'old', 'new' and 'replace' are allowed"
    write(*,'(a)')'Program stops...'
    write(*,"(80('*'))")
    error stop 'Erroneous attribute specification when calling fopen to open a file'
  else
    open(filehandle, FILE=trim(filename), STATUS=trim(attribute), IOSTAT=ios, ERR=999)
  endif
  
999  if(ios .ne. 0)then
    call fopen_error(filename,filehandle,ios)
  endif
end subroutine fopen
!
! close a file <filename> at position (unit) <filehandle>
subroutine fclose(filename, filehandle)
  IMPLICIT NONE
  integer, intent(in) :: filehandle
  integer :: ios
  character(len=*), intent(in) :: filename
  
  close(UNIT=filehandle, IOSTAT=ios, ERR=999)
  
999  if(ios .ne. 0)then
    call fclose_error(filename,filehandle,ios)
  endif
end subroutine fclose

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                         general processing of texts and strings            WORKS: yes !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Allow for commentary lines in the file.
! reading routine that ignores lines with a preceeding star in the first field.
! read in lines of a file and return the first that does NOT start with a star.
function readnocomment(filehandle,stat) result(stream)
  use :: constants, only : ll
  IMPLICIT NONE
  
  integer, intent(in) :: filehandle
  integer, intent(out) :: stat
  integer :: idx, i
  character(len=ll) :: stream
  idx = 0
    read(filehandle,'(a)',IOSTAT=stat)stream
    idx = index(trim(stream),'*')
    if(idx .eq. 0)then
!     we have an uncommented line, return to the calling routine
      return
    elseif(idx .gt. 0)then
!     encountered a comment starting at the field idx of line stream.
!     exchange the * and all following characters with blanks.
      do i=idx, len(stream),1
        stream(i:i) = ''
      enddo
      return
    endif
end function readnocomment


! convert all lowercase letters in a read in string to uppercase letters
function lc2uc( Input_String ) result( Output_String )
  IMPLICIT NONE
  character(len=*), intent(in) :: Input_String
  character(len=len(Input_String)) :: Output_String
  integer :: i, n
! Copy input string
  Output_String = Input_String
! Convert case character by character
  do i = 1, len(Output_String), 1
    n = index(lower_case, Output_String(i:i))
    if ( n .ne. 0 ) Output_String(i:i) = upper_case(n:n)
  enddo
end function lc2uc

! set the cursor using standard reading of the opened file <filename> at <filehandle> to the line of a keyword <key>. 
! if attribute is set to 'mandatory', throw and error, if <key> is not found.
! if attribute is set to 'optional', throw no error, but return a status that says the keyword is not present.
! STOP on positive stat values (erroneous reading)
! STOP on negative stat values only, if the keyword has the attribute mandatory.
! otherwise pass back the stat value and
subroutine set_cursor_no_comm(filehandle, filename, key, attribute, stat)
  use :: constants, only : ll
  IMPLICIT NONE
  integer, intent(in) :: filehandle
  integer :: found
! catch the EOF mark, if possible
  integer, intent(out) :: stat
  character(len=*), intent(in) :: key
  character(len=ll) :: stream
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: attribute
! 1.) rewind the file to beginning
  rewind(filehandle)
!
! 2.) go to the requested keyword
  do 
! use comment sensitive reading:
    read(filehandle,'(a)',IOSTAT = stat)stream
    ! check the status of the read command
    if (stat .lt. 0)then
      if(trim(attribute) .eq. 'mandatory')then
        call reached_EOF(filename, key)
      elseif(trim(attribute) .eq. 'optional')then
        return
      else
        write(*,"(80('*'))")
        write(*,'(a,a,a)')'ERROR: attribute = ',trim(attribute),'was erroneusly specified'
        write(*,'(a)')"only attribute = 'mandatory' or 'optional' are allowed."
        write(*,'(a)')'Program stops now...'
        write(*,"(80('*'))")
        error stop 'erroneous attribute value in set_cursor'
      endif
    elseif(stat .eq. 0)then
!     lc2uc converts the read in string to all uppercase letters, key must be given in uppercase.
      found = index(lc2uc(stream),trim(key))
      if(found .ne. 0)then
        found = 0
        return
      endif
    else
      call read_error(filename,key,stream,stat)
    endif
  enddo
end subroutine set_cursor_no_comm

! set the cursor of the opened file <filename> at <filehandle> to the line of a keyword <key>. 
! using comment sensitive reading if attribute is set to 'mandatory', throw and error, if <key> is not found.
! if attribute is set to 'optional', throw no error, but return a status that says the keyword is not present.
! STOP on positive stat values (erroneous reading)
! STOP on negative stat values only, if the keyword has the attribute mandatory.
! otherwise pass back the stat value and
subroutine set_cursor(filehandle, filename, key, attribute, stat)
  use :: constants, only : ll
  IMPLICIT NONE
  integer, intent(in) :: filehandle
  integer :: found
! catch the EOF mark, if possible
  integer, intent(out) :: stat
  character(len=*), intent(in) :: key
  character(len=ll) :: stream
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: attribute
! 1.) rewind the file to beginning
  rewind(filehandle)
!
! 2.) go to the requested keyword
  do 
! use comment sensitive reading:
!    read(filehandle,'(a)',IOSTAT = stat)stream
     stream = readnocomment(filehandle,stat)
    ! check the status of the read command
    if (stat .lt. 0)then
      if(trim(attribute) .eq. 'mandatory')then
        call reached_EOF(filename, key)
      elseif(trim(attribute) .eq. 'optional')then
        return
      else
        write(*,"(80('*'))")
        write(*,'(a,a,a)')'ERROR: attribute = ',trim(attribute),'was erroneusly specified'
        write(*,'(a)')"only attribute = 'mandatory' or 'optional' are allowed."
        write(*,'(a)')'Program stops now...'
        write(*,"(80('*'))")
        error stop 'erroneous attribute value in set_cursor'
      endif
    elseif(stat .eq. 0)then
!     lc2uc converts the read in string to all uppercase letters, key must be given in uppercase.
      found = index(lc2uc(stream),trim(key))
      if(found .ne. 0)then
        found = 0
        return
      endif
    else
      call read_error(filename,key,stream,stat)
    endif
  enddo
end subroutine set_cursor



end module utils