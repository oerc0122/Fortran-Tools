program timer_test

  use timer
  type(timer_type) :: tmr

  call init_timer_system(tmr, 0)
  tmr%max_depth = 7
  call start_timer(tmr, 'Test1')
  call start_timer(tmr, 'Test2')

  call stop_timer(tmr, 'Test2')
  call stop_timer(tmr, 'Test1')
  call timer_report(tmr)

end program timer_test
