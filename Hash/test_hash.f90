program test_hash

  Use hash
  Use, intrinsic :: iso_fortran_env, only : wp => real64
  Implicit None

  Type( hash_table ) :: table
  Real(kind = wp) :: test_float
  Integer :: test_int
  Complex(kind = wp) :: test_comp

  open(unit=50, file="test_new_control")

  call table%init(30)

  call table%set('float', 3.1415926535897931_wp)
  call table%set('int', 42)
  call table%set('complex', (1.0_wp, 1.0_wp))

  call table%get('float', test_float)
  call table%get('int', test_int)
  call table%get('complex', test_comp)

  print*, test_float
  print*, test_int
  print*, test_comp

  if (abs(test_float - 4.0_wp*atan(1.0_wp)) > 1e-12_wp) print*, "Float retrieval failed"
  if (test_int /= 42) print*, "Int retrieval failed"
  if (abs(test_comp - (1.0_wp, 1.0_wp)) > 1e-12_wp) print*, "Complex retrieval failed"



end program test_hash
