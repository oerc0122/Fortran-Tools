program whats

  use units, only: initialise_units, convert_units
  Implicit None

  Integer, Parameter :: dp = selected_real_kind(15,300)

  Character(Len=256) :: input
  Real(kind=dp) :: val
  Character(Len=256) :: unit_in, unit_out
  Real(kind=dp) :: res

  call initialise_units()

  do
     read('(A)'), input
     if (trim(input) == '') exit
     call parse_input(input, val, unit_in, unit_out)
     res = convert_units(val, unit_in, unit_out)
     print('(A,1X,"=",es12.5e2,1X,A)'), trim(input), res, trim(unit_out)
  end do

contains
  subroutine parse_input(in, val, unit_in, unit_out)
    Character(Len=256), Value :: in
    Character(Len=256) :: word
    Real(kind=dp) :: val
    Character(Len=256) :: unit_in, unit_out

    call get_word(in, word)
    read(word, *) val
    call get_word(in, word)
    read(word, '(A)') unit_in
    call get_word(in, word)
    ! "in"
    call get_word(in, word)
    read(word, '(A)') unit_out

  end subroutine parse_input

  subroutine get_word(in, word)
    Character(Len=*) :: in
    Character(Len=256) :: word
    Integer :: next_word

    next_word = scan(in, ' ')
    word = in(:next_word)
    in = adjustl(in(next_word:))

  end subroutine get_word

end program whats
