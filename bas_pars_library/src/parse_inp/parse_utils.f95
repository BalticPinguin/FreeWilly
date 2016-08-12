!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bas_pars v.0.1 - fork of the AUGER program                          !
!                                                                     !
! Module:  parse_utils                                                !
! Version: 0.1                                                        !
! Purpose: contains different utilities.                              !
!          Currently:                                                 !
!          => different Error routines                                !
!          => helper routinies for input reading and processing       !
!                                                                     !
! Author:  Gilbert Grell, University of Rostock                       !
! Date:    15.05.2016                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module parse_utils
IMPLICIT NONE
! List of characters for case conversion
character(len=*), parameter, private :: lower_case = 'abcdefghijklmnopqrstuvwxyz'
character(len=*), parameter, private :: upper_case = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
! List of numbers for various purposes
character(len=*), parameter, private :: numeric_chars = '0123456789'
! List of characters for angular momentum labels
character(len=*), parameter, private :: angular_subshells = 'SPDFGH'
!
contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                Error section 1.) errors during the read of the input file                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Block termination statement reached too early, when looking for a specific keyword
subroutine early_term(bkey,key)
  use :: global_data, only : INFILE
  IMPLICIT NONE
  character(len=*), intent(in) :: bkey, key
  write(*,"(80('*'))")
  write(*,'(a,a,a,a)')'ERROR when reading the ', trim(bkey),' BLOCK in the input file: ',trim(INFILE)
  write(*,'(a,a)')'reached block termination statement before finding the keyword/data: ',trim(key)
  write(*,'(a)')'Program stops now ...'
  write(*,"(80('*'))")
  error stop 'Error while reading the input file'
end subroutine early_term
!
! keyword within block has an illegal value:
subroutine illegal_value(bkey,property,errvalue,corrvalue)
  use :: global_data, only : INFILE
  IMPLICIT NONE
  character(len=*), intent(in) :: bkey
  character(len=*), intent(in)  :: property
  character(len=*), intent(in)  :: errvalue, corrvalue
  write(*,"(80('*'))")
  write(*,'(a,a,a,a)')'ERROR when reading the ',trim(bkey),' BLOCK in the input file: ',trim(INFILE)
  write(*,'(a,a,a,a)')trim(property),' = ',trim(errvalue), ' was illegally specified.'
  write(*,'(a,a,a,a,a)')'Only ',trim(property),' = ',trim(corrvalue),' is allowed!'
  write(*,'(a)')'Program stops now ...'
  write(*,"(80('*'))")
  error stop 'Error while reading the input file'
end subroutine illegal_value
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                             Input reading & processing routines                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! convert all separators in a string to spaces.
function separator2space( inputString, separator ) result( outputString )
  character(len=*), intent(in) :: inputString, separator
  character(len=len(inputString)) :: outputString
  integer :: i,n
  ! check if the separator is a single character. If not throw an error
  if(len(separator) .ne. 1)then
    write(*,"(80('*'))")
    write(*,'(a)')'Error in function separator2space.'
    write(*,'(a)')"The separator: '",trim(separator),"' is not a single character."
    write(*,"(80('*'))")
    write(*,'(a)')'Program stops...'
    error stop 'Error in function separator2space'
  endif
  ! copy input string to output
  outputString = inputString
  ! Convert separators to spaces character by character
  do i = 1, len(inputString), 1
    n = index(outputString(i:i), trim(separator))
    if(n .ne. 0)then
      outputString(i:i) = ' '
    endif
  enddo
end function separator2space
!
!
! remove single characters as part of an passed string from a string
function remove_chars( inputString, remove) result( outputString )
  character(len=*), intent(in) :: inputString, remove
  character(len=len(inputString)) :: outputString
  integer :: i,n
  ! copy input string to output
  outputString = inputString
  ! Convert all occurences of chars in the remove array to blanks character by character
  do i = 1, len(inputString), 1
    n = index(remove, outputString(i:i))
    if(n .ne. 0)then
      outputString(i:i) = ''
    endif
  enddo
end function remove_chars
!
!
! count the lines between to specified keywords and return the number of lines as a result
! throw an error if something strange happens.
function count_lines(filehandle,filename,bkey,ekey,attribute,stat) result(line)
  use :: constants, only : ll
  use :: utils, only : set_cursor, lc2uc, readnocomment, reached_EOF, read_error
  IMPLICIT NONE
  integer, intent(in) :: filehandle
  integer :: line, efound
  integer, intent(out) :: stat
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: bkey
  character(len=*), intent(in) :: ekey
  character(len=*), intent(in) :: attribute
  character(len=ll) :: stream
!  start at zero lines
  line = 0
! set cursor to the line containing bkey (begin key) and throw an error if the key is not present and the attribute is mandatory.
! if the key is not present and the key is optional, the function shall return control to the calling routine and pass on the stat value.
! this checks and error manangements are done in the set_cursor routine.
  call set_cursor(filehandle,filename,bkey,attribute,stat)
  if(stat .lt. 0)then
    return
  elseif(stat .eq. 0)then
!   count the lines between bkey and ekey   
    do
!     use comment sensitive reading:
      stream = readnocomment(20,stat)
      if(stat .lt. 0)then
        call reached_EOF(filename,ekey)
      elseif(stat .gt. 0)then
        call read_error(filename,ekey,stream,stat)
      endif
      efound = index(lc2uc(stream),ekey)
      if(efound.ne.0)then
        efound = 0
        exit
      else
!       increment the line count if the line is not commented, i.e. empty at this point, or totally empty 
        if(len(trim(stream)) .gt. 0)then
          line = line+1
        endif
      endif
    enddo
  return
  endif
end function count_lines
!
!
! split a given string into two parts, one contianing the numerical part and another that contains the letters
! suited for splitting 4s -> 4, s or 15sp -> 15,sp
! inputString must not contain garbage like 1s5p, but always contigously numbers and characters.
subroutine split_num_char( inputString, outString, outNum)
  use :: utils, only : lc2uc, read_error
  IMPLICIT NONE
  character(len=*), intent(in) :: inputString
  character(len=*), intent(out) :: outString
! contains either c for char of n for num
  integer, intent(out) :: outNum
  integer :: i,numStart, numEnd, cl, cr, cStart, cEnd, nl,nr, stat, a
!
! First extract the angular momentum labels: search for the leftmost occurence of a character
  cStart = 0
  cEnd = 0
!  write(*,'(a,a)')'inputstring = ',lc2uc(inputString)
  do i = 1, len(upper_case), 1
!   c is the leftmost occurence of upper_case(i:i) in inputString.
    cl = index(lc2uc(inputString), upper_case(i:i))
!   cr is the rightmost occurence of uppercase(i:i) in inputstring
    cr = index(lc2uc(inputString), upper_case(i:i), .true.)
!    write(*,'(a,i5,a,i5,a,i5,a,i5)')'cr = ',cr,' cl = ',cl,'  cStart = ',cStart,' cEnd = ',cEnd
    if((cl .gt. 0) .and. (cl .lt. cStart))then
      cStart = cl
    elseif((cl .gt. 0) .and. (cStart .eq. 0))then
      cStart = cl
    endif
    if((cr .gt. 0) .and. (cr .gt. cEnd))then
      cEnd = cr
    endif
  enddo
! check now, if it is also a valid angular momentum subshell label (SPDFGH) after uppercase conversion
  do i = cStart, cEnd,1
    a = index(angular_subshells,lc2uc(inputString(i:i)))
    if(a .eq. 0)then
!   the character inputString(i:i) is not in the list of angular_subshells => ERROR
      write(*,"(80('*'))")
      write(*,'(a,a)')'ERROR in split_num_char. The definition of the #GTOs <Subshell-label>: ',trim(inputString)
      write(*,'(a)')'contains the prohibited character: ',upper_case(i:i),'. Only S,P,D,F,G,H are allowed.'
      write(*,"(80('*'))")
      write(*,'(a)')'Program stops...'
      error stop 'Error in some basis set contraction specification.'
    endif
  enddo 
! if everything is fine, write this string to the output string
  outString = lc2uc(inputString(cStart:cEnd))
! 
! second: extract the number of GTOs. Search for left and rightmost occurence of numerical chars
  numStart = 0
  numEnd = 0
  do i = 1, len(numeric_chars), 1
!   nl is the leftmost occurence of numeric_char(i:i) in inputString.
    nl = index(inputString, numeric_chars(i:i))
!   nr is the rightmost occurence of numeric_char(i:i) in inputstring
    nr = index(inputString, numeric_chars(i:i), .true.)
    if((nl .gt. 0) .and. (nl .lt. numStart))then
      numStart = nl
    elseif((nl .gt. 0) .and. (numStart .eq. 0))then
      numStart = nl
    endif
    if((nr .gt. 0) .and. (nr .gt. numEnd))then
      numEnd = nr
    endif
  enddo
! if everything is fine, write this string to the output number
  read(inputString(numStart:numEnd),*,IOSTAT=stat)outNum
  if(stat .ne. 0)then
    call read_error("inputString(numStart,numEnd)",'# of GTOS in basis',inputString(numStart:numEnd),stat)
  endif
end subroutine split_num_char
!
!
! routine to read in the lists of spin free and spin orbit states that should be taken into consideration
! WORKS: yes!
subroutine get_state_list(filename, spin, infin, stream, statelist)
  use :: constants, only : li, lv, ll, lk
  use :: utils, only : lc2uc, read_error
  IMPLICIT NONE
! specify variables of the subroutine
! key is found in stream at position (beginning key, end key, normal key)
  integer :: found
! misc incrementable variables
  integer :: comma, dash, i, j, iList
! start endend of an interval
  integer :: iStart, iEnd
! status of the read statement
  integer :: stat
! number of dummy entries in stream before the statelist starts
  integer :: nDum
! name of the inputfile
  character(len=*), intent(in) :: filename
! string containing the spin swithc (either SO or SF)
  character(len=*), intent(in) :: spin
! string containing the switch if we have initial or final states: (initial or final)
  character(len=*), intent(in) :: infin
! stream containing the line with the state list
  character(len=*), intent(inout) :: stream
! string that holds the intervals of states given in the input file (F/ISOCLIST and F/ISFLIST)
  character(len=li), dimension(:), allocatable :: interval
! string to pass errorneous values as words to error routines:
  character(len=lv) :: errvalue
! string to pass correct values as words to error routines:
  character(len=lv) :: corrvalue
! string to pass the current property to error routines:
  character(len=lv) :: propvalue
! block from which the statelist is read 
  character(len=lv) :: block
! dummy string:
  character(len=ll) :: dummy
! key string for error output:
  character(len=lk) :: key
! number of spin states that should be presenti n the list:
  integer :: nStates
! list of states that should be filled
  integer, intent(inout), dimension(1:) :: statelist

! initialize the number of states as the length of the statelist
  nStates = size(statelist,1)
! initialize the found variable to 0
  found = 0
! check if spin and infin are set accordingly
  if((lc2uc(spin) .ne. 'SO') .and. (lc2uc(spin) .ne. 'SF'))then
    write(*,"(80('*'))")
    write(*,'(a)')
    write(*,'(a)')'Program stops now'
    write(*,"(80('*'))")
    error stop 'Illegal spin value =/= SF or SO was passed to the get_state_list routine'
  endif
  if((lc2uc(infin) .ne. 'FINAL') .and. (lc2uc(infin) .ne. 'INITIAL'))then
    write(*,"(80('*'))")
    write(*,'(a)')
    write(*,'(a)')'Program stops now'
    write(*,"(80('*'))")
    error stop 'Illegal infin value =/= FINAL or INITIAL was passed to the get_state_list routine'
  endif
!
! this routine has 4 cases. initial and final for each SO and SF states
! initial and final, only affect the format of the error messages
! SO and SF affect error messages and the format of the stream.
! for SO, the stream contains two entries before the list starts :
! SOC NSOC <statelist>
! for SF, the stream contains 3 entries before the list:
! MULT, NSD, NSTATE, <STATELIST>
!
! the prelist entries are always read as dummies.
! thus the number of dummies is decided by the SO and SF switch.
  if(lc2uc(spin) .eq. 'SO')then
    nDum = 2
!   key for the error output if read errors appear
    key = 'SOC'
  elseif(spin .eq. 'SF')then
    nDum = 3
!   key for the error output, if read errors appear
    key = 'SF states'
  endif
! the infin switch marks the block from which the statelist is to be read
  if(lc2uc(infin) .eq. 'FINAL')then
    block = 'FINALWAVE'
  elseif(lc2uc(infin) .eq. 'INITIAL')then
    block = 'INITIALWAVE'
  endif  
! read in the statelist, if it is present.
! allow a space seperated list, as well as comma separated intervals: 1 2 3 4 5 = 1-5 = 1-2,3-5 = 1,2,3-5
! check if comma or dashes are present in stream:
! count the occurences of the dashes
  dash = 0
  do i = 1, len( trim(stream) ), 1
    found = 0
    found = index(stream(i:i),'-')
    if(found .ne. 0)then
      dash = dash + 1
    endif
  enddo
! count the occurences of commas
  comma = 0
  do i = 1, len( trim(stream) ), 1
    found = 0
    found = index(stream(i:i),',')
    if(found .ne. 0)then
      comma = comma + 1
    endif
  enddo
  if((dash .eq. 0) .and. (comma .eq. 0))then
!   we have either NO statelist or a full statelist seperated by spaces
!   try to read the third element from the trimmed stream array, which should correspond to the first element of the statelist
!   if the reading command returns an error, statelist is not specified.
    read(stream,*,IOSTAT = stat),(dummy, i = 1, nDum, 1), dummy
!   negative stat value means that one trys to read after the EOF. here this means after the end of the string.
    if(stat .lt. 0)then
!     no statelist was specified.
!     thus use the default convention for statelist: 1-nStates
      do i = 1, nStates, 1
        statelist(i) = i
      enddo
    elseif(stat .eq. 0)then
!     a third element is present in the stream. 
!     thus the full space seperated list of states is given
!     read it in and exit on error
      read(stream,*,IOSTAT = stat),(dummy, i = 1, nDum, 1), (statelist(i), i = 1, nStates)
      if(stat .lt. 0)then
        write(errvalue,'(a,i5)')'<', nStates
        write(corrvalue,'(I5)')nStates
        write(propvalue,'(a,a,a,a,a)')'number of listed ',trim( lc2uc(infin) ),' ',trim( lc2uc(spin) ), ' states'
        write(block,'(a)')
        call illegal_value(block, propvalue, errvalue, corrvalue)
      elseif(stat.gt.0)then
        call read_error("'streamed line from input file to read statelist'", key, stream, stat)
      endif
    elseif(stat .gt. 0)then
!     an error occured upon reading - abort
      write(*,'(a)')'too'
      call read_error("'streamed line from input file to read statelist'", key, stream, stat)
    endif
  elseif((dash .gt. 0) .or. (comma .gt. 0))then
!   statelist is specified as a comma separated list of values or intervals
    if(dash .eq. 0)then
!     if no dashes are present the list is just comma separated instead of spaces.
!     exchange commas with spaces and read it in
      stream = separator2space(stream,',')
!     read in the list and throw  an error if something happens
      read(stream,*,IOSTAT = stat),(dummy, i = 1, nDum, 1), (statelist(i), i = 1, nStates)
      if(stat .ne. 0)then
        call read_error("'streamed line from input file to read statelist'", key, stream, stat)
      endif
    elseif(comma .eq. 0)then
!     no commas are present, the whole list is specified as one interval with a dash.
!     check if the dash count is one:
      if(dash .gt. 1)then
        write(*,"(80('*'))")
        write(*,'(a,a,a,a)')'ERROR when reading the ', block, ' BLOCK in the input file: ', trim(filename)
        write(*,'(a)')'Illegal specification of the statelist.'
        write(*,'(a)')'Only space or comma separated lists of all states or comma separated'
        write(*,'(a)')'lists of intervals with dashes and individual states are allowed!'
        write(*,'(a)')'For example: 1 2 3 4 5; 1,2,3,4,5 or 1-2,3-5 or 1,2,3-5 are allowed.'
        write(*,'(a)')'Program stops now ...'
        write(*,"(80('*'))")
        error stop 'Error while reading the input file'
      endif
!     turn the dash into a space and read the interval boundaries iStart-iEnd
      stream = separator2space(stream,'-')
      read(stream,*,IOSTAT = stat), (dummy, i = 1, nDum, 1), iStart, iEnd
      if(stat .ne. 0)then
        call read_error("'streamed line from input file'",key,stream,stat)
      endif
!     check if the interval contains as many SOC states as given by nStates:
      if((iEnd-iStart+1) .ne. nStates)then
         write(errvalue,'(i6)')iEnd-iStart+1
         write(corrvalue,'(a,i6)')'nStates = ',nStates
         write(propvalue,'(a,a,a)')'# of ',lc2uc(spin),' in the interval'
        call illegal_value(block, propvalue, errvalue, corrvalue)
      endif
      do i = iStart, iEnd, 1
        statelist(i-iStart+1) = i
      enddo
    else
!     if both, dashes and commas are present, the list is specified in the most general format:
!     e.g. 1,2,3-4,5,7-9
!     turn the commas to spaces and read in the states or intervals of states. 
      stream = separator2space(stream,',')
!     allocate the char array holding the intervals or states: there are comma+1 entries.
      allocate(interval(1:comma+1))
      read(stream,*,IOSTAT=stat),(dummy, i = 1, nDum, 1), (interval(i), i = 1, comma+1)
      if(stat .ne. 0)then
        call read_error("'streamed line from input file'", key, stream, stat)
      endif
!     distinguish interval entries that are single states (one integer number and no dash in the entry) from those 
!     that are real intervals (contain exactly 2 integer values and one dash)
!     scan through the interval array:
      found = 0
      iList = 0
      do i = 1, comma+1, 1
        found = index(interval(i), '-')
        if(found .eq. 0)then
!         found a single state - read it in
!         increment the listindex by one
          iList = iList + 1
          read(interval(i),*,IOSTAT=stat)statelist(iList)
          if(stat .ne. 0)then
            write(propvalue,'(a,1x,a)')trim( lc2uc(infin) ), trim( lc2uc(spin) ),' state list'
            call read_error(propvalue, key, interval(i), stat)
          endif
        elseif(found .ne. 0)then
          found = 0
!         found an interval - count the number of dashes, check that it is only one
          dash = 0
          do j = 1, len( trim( interval(i) ) ), 1
            if( 0 .ne. index( interval(i)(j:j),  '-'))then
              dash = dash+1
            endif
          enddo
!         if more than one dash is found in the interval, throw an error
          if(dash .ne. 1)then
            write(*,"(80('*'))")
            write(*,'(a,a,a,a)')'ERROR when reading the ', block, ' BLOCK in the input file: ', trim(filename)
            write(*,'(a)')'Illegal specification of the statelist.'
            write(*,'(a)')'Only space or comma separated lists of all states or comma separated'
            write(*,'(a)')'lists of intervals with dashes and individual states are allowed!'
            write(*,'(a)')'For example: 1 2 3 4 5; 1,2,3,4,5 or 1-2,3-5 or 1,2,3-5 are allowed.'
            write(*,'(a)')'Program stops now ...'
            write(*,"(80('*'))")
            error stop 'Error while reading the input file'
          endif
!         turn the dash into a space:
          interval(i) = separator2space(interval(i), '-')
!         read in the interval borders:
          read(interval(i),*,IOSTAT=stat)iStart, iEnd
!         write the states from the interval into the list and increment the index accordingly
          do j = iStart, iEnd, 1
            iList = iList + 1
            statelist(iList) = j
          enddo
        endif
      enddo
!     check if iList is equal to nStates, if not, throw an error.
      if(nStates .ne. iList)then
        write(errvalue,'(I5)')iList
        write(corrvalue,'(I5)')nStates
        write(propvalue,'(a,a,1x,a,a)')'# of listed ', trim( lc2uc(infin) ), trim( lc2uc(spin) ),' states'
        call illegal_value(block, propvalue , errvalue, corrvalue)
      endif
!     deallocate the intervals array
      deallocate(interval)
    endif
  endif
!
end subroutine get_state_list
!

end module parse_utils