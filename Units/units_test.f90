program units_test

  use units

  Implicit None
  Integer, Parameter :: dp = selected_real_kind(15, 300)
  Real(kind=dp) :: unit_test

  call initialise_units()

  unit_test = convert_units(1.0_dp, 'ang', 'm')
  print*, "1 Ang in m = ",unit_test

  unit_test = convert_units(1.0_dp, 'e.V', 'J')
  print*, "1 eV in J = ",unit_test

  unit_test = convert_units(1.21_dp, 'GW', 'dyn.pc/aut')
  print*, "1.21GW in dyn.pc/aut = ", unit_test

end program units_test
