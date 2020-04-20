Module hash
  !!-----------------------------------------------------------------------
  !!
  !! Module containing hash table routines for reading control input
  !! Hash table is a fixed-size hash table with open addressing
  !!
  !! author    - j.wilkins march 2020
  !! contributions - i.j.bush april 2020
  !!-----------------------------------------------------------------------

  Use, intrinsic :: iso_fortran_env, only : dp => real64
  Implicit None

  Private

  Integer, Parameter, Public :: STR_LEN = 256
  Character(Len=*), Parameter :: BAD_VAL = "VAL_NOT_IN_KEYS"

  Type, Public :: container
     !! Generic data container
     Private
     Class( * ), Allocatable, Private :: data
   Contains
     Generic, Public :: Assignment( = ) => set, get
     Procedure,            Private :: set => set_container
     Procedure, Pass( C ), Private :: get => get_container
  End type container

  Type, Public :: hash_table
     !! Type containing hash table of parameters
     Private
     Type(container), dimension(:), allocatable :: table_data
     Character(Len=STR_LEN), dimension(:), allocatable :: table_keys
     Character(Len=STR_LEN), dimension(:), allocatable :: key_names
     Integer :: collisions = 0
     Integer :: used_keys = -1
     Integer :: size = -1
     !> Values in hash table can be overwritten: Default = False
     Logical :: can_overwrite = .false.
     Logical :: allocated = .false.
   Contains

     Private
     Procedure, Public, Pass :: init => allocate_hash_table
     Procedure, Public, Pass :: set => set_hash_value
     Generic  , Public  :: get => get_int, get_double, get_complex
     Procedure, Public, Pass :: hash => hash_value
     Procedure, Public, Pass :: keys => print_keys
     Procedure, Public, Pass(table_to) :: fill => fill_from_table
     Procedure, Public, Pass :: copy => copy_table
     Procedure, Public, Pass :: resize => resize_table
     Procedure, Public, Pass :: expand => expand_table
     Procedure, Private :: get_int, get_double, get_complex
     Procedure, Public :: get_cont => get_hash_value
     Procedure, Private, Pass :: get_loc => get_loc
     Procedure, Public, Pass :: in => contains_value
     Final :: cleanup

  End Type hash_table

  public :: error

Contains

  Subroutine set_container( C, stuff )
    !!-----------------------------------------------------------------------
    !!
    !! subroutine for setting generic data container
    !!
    !! author    - i.j.bush april 2020
    !!-----------------------------------------------------------------------

    Implicit None

    Class( container ), Intent( InOut ) :: C
    Class( *         ), Intent( In    ) :: stuff

    C%data = stuff

  End Subroutine set_container

  Subroutine get_container( stuff, C )
    !!-----------------------------------------------------------------------
    !!
    !! subroutine for getting generic data container
    !!
    !! author    - i.j.bush april 2020
    !!-----------------------------------------------------------------------

    Implicit None

    Class( container ),              Intent( In    ) :: C
    Class( *         ), Allocatable, Intent(   Out ) :: stuff

    stuff = C%data

  End Subroutine get_container

  Subroutine cleanup(table)
    !!-----------------------------------------------------------------------
    !!
    !! subroutine for deallocation of hash table
    !!
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Type(hash_table), Intent( InOut ) :: table
    Integer :: ierr

    if (allocated(table%table_data)) then
       Deallocate(table%table_data, stat=ierr)
       If (ierr /= 0) call error("Error deallocating hash%table_data in cleanup hash table")
    end if

    if (allocated(table%key_names)) then
       Deallocate(table%key_names, stat=ierr)
       If (ierr /= 0) call error("Error deallocating hash%key_names in cleanup hash table")
    end if

    if (allocated(table%table_keys)) then
       Deallocate(table%table_keys, stat=ierr)
       If (ierr /= 0) call error("Error deallocating hash%table_keys in cleanup hash table")
    end if

    table%size = -1
    table%used_keys = -1
    table%collisions = -1
    table%allocated = .false.

  End Subroutine cleanup

  Subroutine allocate_hash_table(table, nbuckets, can_overwrite)
    !!-----------------------------------------------------------------------
    !!
    !! Subroutine for allocation and initialisation of hash table
    !!
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------

    Class(hash_table), Intent( InOut ) :: table

    !> Number of buckets to allocate
    Integer, Intent( In    ) :: nbuckets
    Logical, Optional :: can_overwrite
    Integer :: ierr

    if (present(can_overwrite)) then
       table%can_overwrite = can_overwrite
    end if

    table%size = nbuckets
    table%used_keys = 0
    Allocate(table%table_data(nbuckets), stat=ierr)
    If (ierr /= 0) call error("Error allocating hash%table_data in allocate_hash_table")
    Allocate(table%key_names(nbuckets), stat=ierr)
    If (ierr /= 0) call error("Error allocating hash%key_names in allocate_hash_table")
    Allocate(table%table_keys(nbuckets), stat=ierr)
    If (ierr /= 0) call error("Error allocating hash%table_keys in allocate_hash_table")

    table%table_keys(:) = BAD_VAL
    table%collisions = 0
    table%allocated = .true.

  End Subroutine allocate_hash_table

  Function hash_value(table, input) result(output)
    !!-----------------------------------------------------------------------
    !!
    !! Function to hash string using simple sum(ord(input))%max_hash
    !!
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table) :: table
    Character(Len=*), Intent( In    ) :: input
    Integer :: output
    Integer :: i

    if (input == BAD_VAL) call error("Cannot hash value")
    output = 0
    do i = 1, len_trim(input)
       output = output + ichar(input(i:i))
    end do
    output = mod(output, table%size)

  End Function hash_value

  Function get_loc(table, input, must_find) result(location)
    !!-----------------------------------------------------------------------
    !!
    !! Find location of input or bad_val if not found
    !!
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table) :: table
    Character(Len=*), Intent(In) :: input
    Logical, Intent(In), Optional :: must_find
    Integer :: location
    Character(Len=STR_LEN) :: key


    location = table%hash(input)

    key = table%table_keys(location)
    ! Handle open addressing
    do while (trim(key) /= trim(input))
       if (key == BAD_VAL) then
          exit
       end if
       table%collisions = table%collisions + 1
       location = mod(location + 1, table%size)
       key = table%table_keys(location)
    end do

    if (present(must_find)) then
       if (must_find .and. key == BAD_VAL) call error('Key '//input//' not found in table')
    end if

  End Function get_loc

  Function contains_value(table, input) result(output)
    !!-----------------------------------------------------------------------
    !!
    !! Retrieve stored value from hash table
    !!
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table), Intent( In     ) :: table
    Character(Len=*), Intent( In    ) :: input
    Integer :: location
    Logical :: output

    if (.not. table%allocated) call error('Attempting to get from unallocated table')

    location = table%get_loc(input)
    output = table%table_keys(location) == input

  End Function contains_value

  Function get_hash_value(table, input, default) result(output)
    !!-----------------------------------------------------------------------
    !!
    !! Retrieve stored value from hash table
    !!
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table), Intent( In     ) :: table
    Character(Len=*), Intent( In    ) :: input
    Integer :: location
    Class(*), Intent( In     ), Optional :: default
    Class(*), Allocatable :: output

    if (.not. table%allocated) call error('Attempting to get from unallocated table')

    location = table%get_loc(input)

    output = table%table_data(location)
    if (table%table_keys(location) == BAD_VAL) then
       if (present(default)) then
          output = default
       else
          call error('No data pertaining to key '//input//' in table')
       end if
    end if

  End Function get_hash_value

  Subroutine set_hash_value(table, key, input)
    !!-----------------------------------------------------------------------
    !!
    !! Set table at key to input
    !!
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table), Intent( InOut ) :: table
    Character(Len=*), Intent( In    ) :: key
    Class(*), Intent( In    ) :: input
    Integer :: location

    if (.not. table%allocated) call error('Attempting to set unallocated table')

    location = table%get_loc(key)
    if (.not. table%can_overwrite .and. table%table_keys(location) /= BAD_VAL) then
       call error('Cannot overwrite key '//key)
    end if

    table%table_data(location) = input
    table%table_keys(location) = key
    table%used_keys = table%used_keys + 1
    table%key_names(table%used_keys) = key

  End Subroutine set_hash_value

  Subroutine print_keys(table)
    !!-----------------------------------------------------------------------
    !!
    !! Print all keys in table
    !!
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table), Intent( In    ) :: table
    Integer :: i

    do i = 1, table%used_keys
       print('(A)'), table%key_names(i)
    end do

  End Subroutine print_keys

  Subroutine fill_from_table(table_from, table_to)
    !!-----------------------------------------------------------------------
    !!
    !! Populate table_to with keyvals from table_from
    !!
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table), Intent( In    ) :: table_from
    Class(hash_table), Intent( InOut ) :: table_to
    Integer :: location
    Integer :: i

    do i = 1, table_from%used_keys
       location = table_from%get_loc(table_from%key_names(i))
       call table_to%set(table_from%table_keys(location), table_from%table_data(location))
    end do

  End Subroutine fill_from_table

  Subroutine resize_table(table, new_size)
    !!-----------------------------------------------------------------------
    !!
    !! Make resize table to new_size or empty table if not provided (destroys table)
    !!
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table), Intent( InOut ) :: table
    Integer, Optional, Intent( In    ) :: new_size
    Integer :: table_size

    if (present(new_size)) then
       table_size = new_size
    else
       table_size = table%size
    end if

    call cleanup(table)
    call table%init(table_size)

  End Subroutine resize_table

  Subroutine copy_table(table_from, table_to)
    !!-----------------------------------------------------------------------
    !!
    !! Make table_to a copy of table_from (destroys table_to)
    !!
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table), Intent( In    ) :: table_from
    Class(hash_table), Intent( InOut ) :: table_to

    call table_to%resize(table_from%size)
    call table_to%fill(table_from)

  End Subroutine copy_table

  Subroutine expand_table(table, new_size)
    !!-----------------------------------------------------------------------
    !!
    !! Non-destructively resize table to new_size (must be larger than old size)
    !!
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table), Intent( InOut ) :: table
    Class(hash_table), Pointer :: table_temp
    Integer, Intent( In    ) :: new_size
    Integer :: ierr

    if (new_size > table%size) then
       allocate(table_temp, stat=ierr)
       if (ierr.ne.0) call error('Error allocating table_temp in expand table')
       call table_temp%copy(table)
       call table%resize(new_size)
       call table%fill(table_temp)
       call cleanup(table_temp)
    else
       call error('New size for table less than previous size')
    end if

  End Subroutine expand_table

  Subroutine get_int( table, key, val, default )
    !!-----------------------------------------------------------------------
    !!
    !! get integer type from hash table
    !!
    !! author    - i.j.bush april 2020
    !!-----------------------------------------------------------------------

    Implicit None

    Class( hash_table ), Intent( InOut ) :: table
    Character(Len=*), Intent( In    ) :: key
    Integer               , Intent(   Out ) :: val
    Integer, Intent( In    ), Optional :: default
    Class( * ), Allocatable :: stuff

    stuff = table%get_cont(key, default)

    Select Type( stuff )
    Type is ( Integer )
       val = stuff
    Class Default
       Call error('Trying to get integer from a not integer')
    End Select

  End Subroutine get_int

  Subroutine get_double( table, key, val, default )
    !!-----------------------------------------------------------------------
    !!
    !! get double type from hash table
    !!
    !! author    - i.j.bush april 2020
    !!-----------------------------------------------------------------------

    Implicit None

    Class( hash_table ), Intent( InOut ) :: table
    Character(Len=*), Intent( In    ) :: key
    Real(kind=dp), Intent( In    ), Optional :: default
    Real(kind=dp)               , Intent(   Out ) :: val
    Class( * ), Allocatable :: stuff

    stuff = table%get_cont(key, default)

    Select Type( stuff )
    Type is ( Real( dp ) )
       val = stuff
    Class Default
       Call error('Trying to get real from a not real')
    End Select

  End Subroutine get_double

  Subroutine get_complex( table, key, val, default )
    !!-----------------------------------------------------------------------
    !!
    !! get complex type from hash table
    !!
    !! author    - i.j.bush april 2020
    !!-----------------------------------------------------------------------

    Implicit None

    Class( hash_table ), Intent( InOut ) :: table
    Character(Len=*), Intent( In    ) :: key
    Complex(kind=dp)               , Intent(   Out ) :: val
    Complex(kind=dp), Intent( In    ), Optional :: default
    Class( * ), Allocatable :: stuff

    stuff = table%get_cont(key, default)

    Select Type( stuff )
    Type is ( Complex( dp ) )
       val = stuff
    Class Default
       Call error('Trying to get complex from a not complex')
    End Select

  End Subroutine get_complex

  Subroutine error(message)
    Character(Len=*), Intent( In    ) :: message

    write(0, *) message
    stop
  end Subroutine error

end Module hash
