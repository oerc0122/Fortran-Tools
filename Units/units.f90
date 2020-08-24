Module units
  !!-----------------------------------------------------------------------
  !!
  !! Module to handle unit conversion
  !!
  !! author - j.wilkins march 2020
  !!-----------------------------------------------------------------------


  Use hash, only: hash_table, STR_LEN, error
  Use, intrinsic :: iso_fortran_env, only : dp => real64
  Implicit None

  Private

  Type, Private, Extends(hash_table) :: units_hash_table
   contains
     Generic  , Public  :: get => get_int, get_double, get_complex, get_unit
     Procedure, Private :: get_unit
  end type units_hash_table

  Type, Private :: unit_data
     !! Type containing data corresponding to units
     Character(Len=STR_LEN) :: name
     Character(Len=STR_LEN) :: abbrev
     Real(Kind=dp) :: conversion_to_si
     Integer, Dimension(7) :: dims = [0, 0, 0, 0, 0, 0, 0] ! mass length time temp mol current luminosity
   contains
     Generic, Public :: Operator(*) => unit_mult
     Generic, Public :: Operator(/) => unit_div
     Generic, Public :: Operator(**) => unit_pow
     Procedure, Pass, Private :: unit_mult, unit_div, unit_pow
     Procedure, Nopass :: init => init_unit
  end type unit_data

  Public :: initialise_units
  Public :: convert_units

  Type(unit_data), Private, Parameter :: null_unit = unit_data('', '', 1.0_dp)
  type(units_hash_table), Private, save :: units_table

contains

  Subroutine initialise_units()
    !!-----------------------------------------------------------------------
    !!
    !! Set up units table call before attempting to convert
    !!
    !! author - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Real(kind=dp), Parameter :: electron_charge = 1.602176634e-19_dp
    Real(kind=dp), Parameter :: coulomb = 1.0_dp
    Real(kind=dp), Parameter :: avogadro = 6.022140857e23_dp
    Real(kind=dp), Parameter :: boltz = 1.3806503e-23_dp ! JK-1

    Real(kind=dp), Parameter :: metre = 1.0_dp
    Real(kind=dp), Parameter :: angstrom = 1.0e-10_dp*metre
    Real(kind=dp), Parameter :: bohr = 0.52918_dp*angstrom
    Real(kind=dp), Parameter :: inch = 2.54e-2_dp*metre

    Real(kind=dp), Parameter :: parsec = 3.0856775999999956e16_dp*metre

    Real(kind=dp), Parameter :: joule = 1.0_dp
    Real(kind=dp), Parameter :: calorie = 4.1842_dp*joule
    Real(kind=dp), Parameter :: hartree = 4.359744722e-18_dp*joule

    Real(kind=dp), Parameter :: kilogram = 1.0_dp
    Real(kind=dp), Parameter :: amu = 1.660539e-21_dp*kilogram
    Real(kind=dp), Parameter :: electron_mass = 9.1093837015e-31_dp*kilogram
    Real(kind=dp), Parameter :: pound = 0.45359237_dp*kilogram

    Real(kind=dp), Parameter :: second = 1.0_dp
    Real(kind=dp), Parameter :: aut = 2.4188843265857e-17_dp*second

    Real(kind=dp), Parameter :: newton = kilogram*metre/second**2

    Real(kind=dp), Parameter :: atmosphere = 101325.0_dp

    Real(kind=dp), Parameter :: pascal = newton/metre**2

    Real(kind=dp), Parameter :: gravity = 9.81_dp*metre/second**2

    call units_table%init(100)

    ! Time

    call units_table%set("hr", init_unit(abbrev="hr", name="Hour", time=1, to_si=3600.0_dp*second))
    call units_table%set("min", init_unit(abbrev="min", name="Minute", time=1, to_si=60.0_dp*second))
    call units_table%set("s", init_unit(abbrev="s", name="Second", time=1, to_si=second))
    call units_table%set("aut", init_unit(abbrev="aut", name="Atomic Time Unit", time=1, to_si=aut))
    call units_table%set("giffy", init_unit(abbrev="giffy", name="Giffy", time=1, to_si=1e-2*second))
    call units_table%set("shake", init_unit(abbrev="shake", name="Shakes of a lamb's tail", time=1, to_si=1e-8*second))

    ! Length

    call units_table%set("ang", init_unit(abbrev="ang", name="Angstrom", length=1, to_si=angstrom))
    call units_table%set("bohr", init_unit(abbrev="bohr", name="Bohr", length=1, to_si=bohr))
    call units_table%set("m", init_unit(abbrev="m", name="Metre", length=1, to_si=metre))
    call units_table%set('in', init_unit(abbrev="in", name="Inch", length=1, to_si=inch))
    call units_table%set("ft", init_unit(abbrev="ft", name="Foot", length=1, to_si=12.0_dp*inch))
    call units_table%set("pc", init_unit(abbrev="pc", name="Parsec", length=1, to_si=parsec))
    call units_table%set("bus", init_unit(abbrev="bus", name="London Bus", length=1, to_si=2.5_dp*metre))

    ! Mass

    call units_table%set("da", init_unit(abbrev="Da", name="Atomic Mass Unit", mass=1, to_si=amu))
    call units_table%set("amu", init_unit(abbrev="amu", name="Atomic Mass Unit", mass=1, to_si=amu))
    call units_table%set("g", init_unit(abbrev="g", name="Gram", mass=1, to_si=1e-3*kilogram))
    call units_table%set("lb", init_unit(abbrev="lb", name="Pound", mass=1, to_si=pound))
    call units_table%set("oz", init_unit(abbrev="oz", name="Ounce", mass=1, to_si=pound/16.0_dp))
    call units_table%set("m_e", init_unit(abbrev="m_e", name="Electron mass", mass=1, to_si=electron_mass))
    call units_table%set("whale", init_unit(abbrev="Whale", name="Blue Whale Mass", mass=1, to_si=100e3_dp*kilogram))

    ! Charge

    call units_table%set("q_e", &
         & init_unit(abbrev="q_e", name="Elementary charge", current=1, time=1, to_si=electron_charge))
    call units_table%set("e", &
         & init_unit(abbrev="e", name="Elementary charge", current=1, time=1, to_si=electron_charge))
    call units_table%set("c", &
         & init_unit(abbrev="C", name="Coulomb", current=1, time=1, to_si=coulomb))

    ! Energy

    call units_table%set("j", init_unit(abbrev="J", name="Joule", mass=1, length=2, time=-2, to_si=joule))
    call units_table%set("cal", &
         & init_unit(abbrev="Cal", name="Calorie", mass=1, length=2, time=-2, to_si=calorie))
    call units_table%set("ha", init_unit(abbrev="Ha", name="Hartree", mass=1, length=2, time=-2, to_si=hartree))
    call units_table%set("e_h", init_unit(abbrev="Ha", name="Hartree", mass=1, length=2, time=-2, to_si=hartree))
    call units_table%set("ry", init_unit(abbrev="Ry", name="Rydberg", mass=1, length=2, time=-2, to_si=0.5_dp*hartree))
    call units_table%set("littleboy", init_unit(abbrev="littleboy", name="Little Boy", &
         & mass=1, length=2, time=-2, to_si=6.3E+13*joule))

    ! Temperature

    call units_table%set("k", init_unit(abbrev="K", name="Kelvin", temp=1, to_si=1.0_dp))

    ! Energy flow
    call units_table%set("w", init_unit(abbrev="W", name="Watt", mass=1, length=2, time=-3, to_si=joule/second))

    ! Area

    call units_table%set("wales", init_unit(abbrev="Wales", name="Are of Wales", length=2, to_si= 2.07350E10_dp*metre**2))
    call units_table%set("barn", init_unit(abbrev="barn", name="Barn", length=2, to_si=1e-28_dp*metre**2))

    ! Pressure

    call units_table%set("atm", init_unit(abbrev="atm", name="Atmosphere", mass=1, length=-1, time=-2, to_si=atmosphere))
    call units_table%set("pa", init_unit(abbrev="Pa", name="Pascal", mass=1, length=-1, time=-2, to_si=pascal))

    ! Force

    call units_table%set("n", &
         & init_unit(abbrev="N", name="Newton", mass=1, length=1, time=-2, to_si=newton))
    call units_table%set("dyn", &
         & init_unit(abbrev="dyn", name="Dyne", mass=1, length=1, time=-2, to_si=1e5_dp*newton))

    ! Velocity

    call units_table%set("auv", &
         & init_unit(abbrev="auv", name="Atomic Velocity Unit", length=1, time=-1, to_si=bohr/aut))
    call units_table%set("beard", &
         & init_unit(abbrev="beard", name="Beard growth", length=1, time=-1, to_si=5e-9*metre/second))

    ! Constants
    call units_table%set("grav", &
         & init_unit(abbrev="g_e", name="9.81 m/s^2", length=1, time=-2, to_si=gravity))
    call units_table%set("k_b", &
         & init_unit(abbrev='k_b', name="Boltzmann Constant", mass=1, length=2, time=-2, temp=-1, to_si=boltz))

    ! Voltage

    call units_table%set("v", init_unit(abbrev="V", name="Volt", mass=1, length=2, time=-3, current=-1, &
         & to_si=joule/coulomb))

    ! Mols

    call units_table%set("mol", init_unit(abbrev="mol", name="Mole", mol=1, to_si=avogadro))

    ! Unitless

    call units_table%set("%", init_unit(abbrev="%", name="%", to_si=100.0_dp))
    call units_table%set("", init_unit(abbrev="", name="", to_si=1.0_dp))

  End Subroutine initialise_units

  Function convert_units(val, from, to) result(res)
    Real(kind=dp), Intent(in) :: val
    Real(kind=dp) :: res
    Character(Len=*) :: from, to
    Type( unit_data ) :: from_unit, to_unit
    Type( unit_data ) :: output

    output = output%init("", "", 1.0_dp)
    from_unit = parse_unit_string(from)
    to_unit = parse_unit_string(to)
    output = to_unit / from_unit
    if (any(output%dims /= 0)) call error('Cannot convert between '//trim(from)//' & '//trim(to)//' different dimensions')
    res = val / output%conversion_to_si

  end Function convert_units

  Subroutine handle_decimal_prefix(string, factor)
    Character(len=*), Intent(inout) :: string
    Integer :: i
    Type(unit_data), intent(out) :: factor
    Character(Len=256) :: tmp
    Character(len=*), Parameter :: prefix_symbol = "YZEPTGMk dcmunpfazy"
    Type(unit_data), Dimension(19), Parameter :: prefix = [ &
         & unit_data(name="Yotta", abbrev="Y", conversion_to_si=1e24_dp), &
         & unit_data(name="Zetta", abbrev="Z", conversion_to_si=1e21_dp), &
         & unit_data(name="Exa",   abbrev="E", conversion_to_si=1e18_dp), &
         & unit_data(name="Peta",  abbrev="P", conversion_to_si=1e15_dp), &
         & unit_data(name="Tera",  abbrev="T", conversion_to_si=1e12_dp), &
         & unit_data(name="Giga",  abbrev="G", conversion_to_si=1e9_dp), &
         & unit_data(name="Mega",  abbrev="M", conversion_to_si=1e6_dp), &
         & unit_data(name="Kilo",  abbrev="k", conversion_to_si=1e3_dp), &
         & null_unit, &
         & unit_data(name="Deci",  abbrev="d", conversion_to_si=1e-1_dp), &
         & unit_data(name="Centi", abbrev="c", conversion_to_si=1e-2_dp), &
         & unit_data(name="Milli", abbrev="m", conversion_to_si=1e-3_dp), &
         & unit_data(name="Micro", abbrev="u", conversion_to_si=1e-6_dp), &
         & unit_data(name="Nano",  abbrev="n", conversion_to_si=1e-9_dp), &
         & unit_data(name="Pico",  abbrev="p", conversion_to_si=1e-12_dp), &
         & unit_data(name="Femto", abbrev="f", conversion_to_si=1e-15_dp), &
         & unit_data(name="Atto",  abbrev="a", conversion_to_si=1e-18_dp), &
         & unit_data(name="Zepto", abbrev="z", conversion_to_si=1e-21_dp), &
         & unit_data(name="Yocto", abbrev="y", conversion_to_si=1e-24_dp)]

    factor = null_unit
    tmp = string(2:)
    call lower_case(tmp)
    if (.not. units_table%in(tmp)) then
       i = index(prefix_symbol, string(2:2))
       tmp = string(3:)
       call lower_case(tmp)
       if (i < 1 .or. .not. units_table%in(tmp)) call error("Unit not found "//string(2:))
       factor = prefix(i)
       string = string(1:1) // string(3:)
    end if

  end Subroutine handle_decimal_prefix

  Function parse_unit_string(string) result(output)
    Type( unit_data ) :: output
    Type( unit_data ) :: tmp_unit
    Character(Len=*), Intent ( in ) :: string
    Character(Len=256) :: curr_parse
    Character(Len=256), Dimension(:), Allocatable :: parsed
    Character(Len=*), Parameter :: number = "1234567890-+"
    Type(unit_data) :: factor
    Integer :: i
    Integer :: n

    call decompose_unit_string("."//string, parsed)

    output = null_unit

    ! Handle powers first
    do i = size(parsed), 1, -1
       curr_parse = parsed(i)
       if (curr_parse(1:1) == "^") then
          if (verify(trim(curr_parse(2:)), number) /= 0) call error("Non numeric power issued")
          read(curr_parse(2:), "(i8.1)") n
          curr_parse = parsed(i-1)
          call handle_decimal_prefix(curr_parse, factor)
          call lower_case(curr_parse)
          call units_table%get(curr_parse(2:), tmp_unit)
          tmp_unit = (factor*tmp_unit) ** n
          select case (curr_parse(1:1))
          case (".")
             output = output * tmp_unit
          case ("/")
             output = output / tmp_unit
          case default
             call error("Cannot parse unit string"//string)
          end select
          parsed(i) = "."
          parsed(i-1) = "."
       end if
    end do

    do i = 1, size(parsed)
       curr_parse = parsed(i)
       call handle_decimal_prefix(curr_parse, factor)
       call lower_case(curr_parse)
       call units_table%get(curr_parse(2:), tmp_unit)

       select case (curr_parse(1:1))
       case (".")
          output = output * (factor * tmp_unit)
       case ("/")
          output = output / (factor * tmp_unit)
       case default
          call error("Cannot parse unit string "//string)
       end select
    end do

    output%abbrev = string

  end Function parse_unit_string

  Subroutine decompose_unit_string(string, output)
    Character(Len=256), Dimension(:), Allocatable :: output
    Character(Len=*), Intent( In ) :: string
    Character(Len=*), Parameter :: alphabet = "+-1234567890abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_%'"//'"'
    Character(Len=*), Parameter :: punc = "./^"
    Integer :: nParts
    Integer :: i, j, k
    Integer :: ierr

    ierr = 0
    if (allocated(output)) deallocate(output, stat=ierr)
    if (ierr /= 0) call error("Error deallocating output in decompose_unit_string")

    nParts = count([(verify(string(j:j), punc) == 0, j=1,len(string))])
    allocate(output(nParts), stat=ierr)
    if (ierr /= 0) call error("Error allocating output in decompose_unit_string")

    j = 1
    do i = 1, nParts
       k = verify(string(j+1:), alphabet)
       if (k == 0) then
          k = len(string) + 1
       else
          k = j+k
       end if
       output(i) = string(j:k-1)
       j = k
    end do
  end Subroutine decompose_unit_string


  type(unit_data) Function unit_mult(a, b)
    Class(unit_data), Intent( in ) :: a, b

    unit_mult = unit_data( &
         name = trim(a%name)//" "//b%name, &
         abbrev = trim(a%abbrev)//b%abbrev, &
         conversion_to_si = a%conversion_to_si*b%conversion_to_si, &
         dims = a%dims + b%dims &
         )
  end Function unit_mult

  type(unit_data) Function unit_div(a, b)
    Class(unit_data), Intent( in ) :: a, b

    unit_div = unit_data( &
         name = trim(a%name)//" per "//b%name, &
         abbrev = trim(a%abbrev)//"/"//b%abbrev, &
         conversion_to_si = a%conversion_to_si/b%conversion_to_si, &
         dims = a%dims - b%dims &
         )
  end Function unit_div

  type(unit_data) Function unit_pow(a, b)
    Class(unit_data), Intent( in ) :: a
    Integer, Intent( in ) :: b
    Character(Len=8) :: b_str

    write(b_str, "(i0.1)") b

    unit_pow = unit_data( &
         name = trim(a%name)//"^"//trim(b_str), &
         abbrev = trim(a%abbrev)//"^"//trim(b_str), &
         conversion_to_si = a%conversion_to_si**b, &
         dims = a%dims*b &
    )

  end Function unit_pow

  type(unit_data) Function init_unit(name, abbrev, to_si, mass, length, time, temp, mol, current, luminosity)
    Character(Len = *), Intent( In    ) :: name, abbrev
    Real(Kind = dp), Intent( In    ) :: to_si
    Integer, Optional, Intent( In    ) :: mass, length, time, temp, mol, current, luminosity

    init_unit = unit_data(name=name, abbrev=abbrev, conversion_to_si=to_si)
    if (present(mass)) then
       init_unit%dims(1) = mass
    end if

    if (present(length)) then
       init_unit%dims(2) = length
    end if

    if (present(time)) then
       init_unit%dims(3) = time
    end if

    if (present(temp)) then
       init_unit%dims(4) = temp
    end if

    if (present(mol)) then
       init_unit%dims(5) = mol
    end if

    if (present(current)) then
       init_unit%dims(6) = current
    end if

    if (present(luminosity)) then
       init_unit%dims(7) = luminosity
    end if

  end Function init_unit

  Subroutine get_unit( table, key, val, default )
    !!-----------------------------------------------------------------------
    !!
    !! get unit type from hash table
    !!
    !! author    - i.j.bush april 2020
    !!-----------------------------------------------------------------------

    Implicit None

    Class( units_hash_table ), Intent( InOut ) :: table
    Character(Len=*), Intent( In    ) :: key
    type(unit_data)            , Intent(   Out ) :: val
    type(unit_data), Intent( In    ), Optional :: default
    Class( * ), Pointer :: stuff

    call table%get_cont(key, default, stuff)

    Select Type( stuff )
    Type is ( unit_data )
       val = stuff
    Class Default
       Call error('Trying to get unit from a not unit')
    End Select

    deallocate(stuff)
    nullify(stuff)

  End Subroutine get_unit

  Subroutine lower_case(str)
    Character(len=*), Intent( InOut ) :: str
    Integer :: i

    do i = 1, len_trim(str)
       Select Case(str(i:i))
       Case ('A':'Z')
          str(i:i) = achar( ichar(str(i:i)) + 32 )
       Case default
          continue
       end Select
    end do
  end Subroutine lower_case

end Module units
