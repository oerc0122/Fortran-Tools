Module timer

#ifdef mpi
  use mpi
#endif
  Use, intrinsic :: iso_fortran_env, only : dp => real64

  Implicit None

  Private

  Integer, Parameter :: max_depth = 6, max_name = 18

  Type :: node_timer
     !!------------------------------------------------!
     !! Timer
     !!------------------------------------------------!
     Character( Len = max_name ) :: name
     Integer           :: id
     Real( Kind = dp ) :: max, min, total, last
     Real( Kind = dp ) :: start, stop
     Integer           :: calls
     Logical           :: running = .false.
  End Type node_timer

  Type :: call_stack
     !!------------------------------------------------!
     !! Call stack
     !!------------------------------------------------!
     Character ( Len = max_name ), Dimension( max_depth ) :: name
     Integer :: depth = 0
  End Type call_stack

  Type :: node
     !!------------------------------------------------!
     !! Tree node
     !!------------------------------------------------!
     Type ( node_timer ) :: time
     Type ( timer_tree ), Pointer :: tree => null()
     Type ( node ), Pointer :: child => null()
     Type ( node ), Pointer :: parent => null()
     Type ( node ), Pointer :: next_sibling => null()
  End Type node

  Type :: timer_tree
     !!------------------------------------------------!
     !! Tree structure
     !!------------------------------------------------!
     Type ( node ), Pointer :: head => null()
     Integer :: n_timers = 0
  End Type timer_tree

  Type, Public :: timer_type
     !!------------------------------------------------!
     !! Main timer system
     !!------------------------------------------------!
     Type ( timer_tree ), pointer :: tree
     Type ( call_stack ) :: stack
     Real( Kind = dp) :: elapsed
     Logical :: proc_detail = .false.
     Integer :: max_depth = max_depth
     Integer :: proc_id
     Integer :: out_unit
  End Type timer_type

  Interface timer_write
     Module Procedure timer_write_sing
     Module Procedure timer_write_mul
  End Interface timer_write

  Public :: timer_report
  Public :: timer_last_time
  Public :: start_timer
  Public :: stop_timer
  Public :: init_timer_system
  Public :: dump_call_stack
  Public :: start_timer_path
  Public :: stop_timer_path

Contains

  Subroutine dump_call_stack ( stack )
    !!------------------------------------------------!
    !!
    !! Print out a call stack to dump current location
    !!
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type ( call_stack ) :: stack
    Integer :: i

    Call timer_write('')
    Call timer_write('Process stack:')
    Do i = 1, stack%depth
       call timer_write(stack%name(i))
    End Do

  End Subroutine dump_call_stack

  Subroutine push_stack ( tmr, name )
    !!------------------------------------------------!
    !!
    !! Add a timer to the call stack
    !!
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type ( timer_type ), Intent( InOut ) :: tmr
    Character ( Len = * ), Intent ( In    ) :: name

    tmr%stack%depth = tmr%stack%depth + 1
    If ( tmr%stack%depth > max_depth ) &
         & Call timer_error(tmr, 'Call stack exceeds max depth : recursive call or unended timer?')
    tmr%stack%name(tmr%stack%depth) = name

  End Subroutine push_stack

  Subroutine pop_stack ( tmr, name )
    !!------------------------------------------------!
    !!
    !! Remove a timer from the call stack
    !!
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type ( timer_type ), Intent( InOut ) :: tmr
    Character ( Len = * ), Intent( In    ) :: name

    If ( name /= tmr%stack%name(tmr%stack%depth)) &
         & Call timer_error(tmr, 'Child timer '//name//' ended before parent')
    tmr%stack%name(tmr%stack%depth) = ''
    tmr%stack%depth = tmr%stack%depth - 1

  End Subroutine pop_stack

  Subroutine init_timer_system ( tmr, nrite )
    !!------------------------------------------------!
    !!
    !! Initialise a timer system
    !!
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type ( timer_type ), intent ( inout ) :: tmr
    Integer, intent ( in ) :: nrite
    Integer :: ierr

#ifdef mpi
    Call mpi_comm_rank(mpi_comm_world, tmr%proc_id, ierr)
    if (ierr /= 0) call timer_error(tmr, 'MPI Comm rank error')
#else
    tmr%proc_id = 0
#endif
    tmr%out_unit = nrite
    Allocate(tmr%tree)
    Allocate(tmr%tree%head)
    tmr%tree%head%tree => tmr%tree
    Call init_timer ( tmr%tree%head%time, 'Head')
    Call start_timer( tmr, 'Main' )

  End Subroutine init_timer_system

  Function find_timer ( tmr, name, stack_in ) result(current)
    !!------------------------------------------------!
    !!
    !! Locate a timer node witin a given timer system
    !!
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type ( timer_type ), Intent ( InOut ) :: tmr
    Character ( Len = * ) :: name
    Type ( call_stack ), Optional :: stack_in
    Type ( call_stack ) :: stack
    Type ( node ), Pointer :: current
    Integer :: depth

    If (Present(stack_in)) Then
       stack = stack_in
    Else
       stack = tmr%stack
    End If

    current => tmr%tree%head
    depth = 0

    Do While ( depth < stack%depth )
       If ( .not. associated(current%child)) call timer_error(tmr, 'Call stack does not match call tree (no child)')
       depth = depth + 1
       current => current%child
       Do While ( current%time%name /= stack%name(depth) )
          If ( .not. associated(current%next_sibling)) &
               & Call timer_error(tmr, 'Call stack does not match call tree (no sibling)')
          current => current%next_sibling
       End Do
    End Do

    If (current%time%name == name ) Then
       Continue
    Else If (.not. associated(current%child)) Then
       Call init_child_node(name, current)
       current => current%child
       Return
    Else
       current => current%child
       do while ( current%time%name /= name )
          If ( .not. Associated(current%next_sibling)) Then
             Call init_sibling_node(name, current)
             current => current%next_sibling
             Return
          End If
          current => current%next_sibling
       End Do
    End If

  End Function find_timer

  Subroutine init_child_node(name, parent)
    !!------------------------------------------------!
    !!
    !! Create a timer node which is a child of the parent node
    !!
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Character ( len = * ), Intent( In    ) :: name
    Type ( node ), Target :: parent

    Allocate ( parent%child )
    parent%child%parent => parent
    parent%child%tree => parent%tree
    parent%tree%n_timers = parent%tree%n_timers + 1
    Call init_timer(parent%child%time, name)

  End Subroutine init_child_node

  Subroutine init_sibling_node(name, sibling)
    !!------------------------------------------------!
    !!
    !! Create a timer node which is a child of the parent node
    !!
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Character ( len = * ) :: name
    Type ( node ) :: sibling
    Type ( node ), pointer :: child

    Allocate ( sibling%next_sibling )
    child => sibling%next_sibling
    child%parent => sibling%parent

    sibling%next_sibling%tree => sibling%tree
    sibling%tree%n_timers = sibling%tree%n_timers + 1

    Call init_timer(child%time, name)

  End Subroutine init_sibling_node

  Subroutine start_timer(tmr, name, stack)
    !!------------------------------------------------!
    !!
    !! Start a timer running on a given timer system
    !! If stack is supplied bypass standard call stack
    !!
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type ( timer_type ), Intent( InOut ) :: tmr
    Character ( Len = * ), Intent( In    )  :: name
    Type ( call_stack ), optional :: stack
    Type ( node ), Pointer :: current_timer

    current_timer => find_timer(tmr, name, stack)

    if (.not. present(stack)) Call push_stack(tmr, name)
    Call get_time(current_timer%time%start)
    current_timer%time%running = .true.

  End Subroutine start_timer

  Subroutine stop_timer(tmr, name, stack)
    !!------------------------------------------------!
    !!
    !! Stop a timer running on a given timer system
    !! If stack is supplied bypass standard call stack
    !!
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type ( timer_type ), Intent( InOut ) :: tmr
    Character ( Len = * ), Intent( In    )  :: name
    Type ( call_stack ), optional :: stack
    Type ( node ), Pointer :: current_timer

    current_timer => find_timer(tmr, name, stack)

    if (.not. present(stack)) Call pop_stack(tmr, name)
    If ( .not. current_timer%time%running ) &
         & Call timer_error(tmr, 'Timer '//Trim(current_timer%time%name)//' stopped but not running')

    Call get_time(current_timer%time%stop)

    current_timer%time%running = .false.
    current_timer%time%last = current_timer%time%stop - current_timer%time%start
    If ( current_timer%time%last > current_timer%time%max ) current_timer%time%max = current_timer%time%last
    If ( current_timer%time%last < current_timer%time%min ) current_timer%time%min = current_timer%time%last
    current_timer%time%total = current_timer%time%total + current_timer%time%last
    current_timer%time%calls = current_timer%time%calls + 1

  End Subroutine stop_timer

  Subroutine start_timer_path(tmr, name_in, start_parents)
    !!------------------------------------------------!
    !!
    !! Start a timer running on a given timer system
    !! ignoring the timer call stack
    !! If start_parents start all non-running timers
    !! which are on the path -- Default TRUE
    !! - Path should be colon separated
    !!
    !!
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type ( timer_type ), Intent( InOut ) :: tmr
    Character ( Len = * ), Intent( In    )  :: name_in
    Logical, Intent ( In    ), Optional :: start_parents
    Logical :: parents
    Character ( Len = max_name ) :: name
    Type ( call_stack ) :: stack
    Type ( node ), pointer :: is_running
    Type ( node ), Pointer :: current_timer
    Integer :: depth

    Call timer_split_stack_string(tmr, name_in, stack, name)
    current_timer => find_timer(tmr, name, stack)

    parents = .true.
    if ( present(start_parents) ) parents = start_parents

    if (parents) then
       do depth = 1, tmr%stack%depth
          tmr%stack%depth = depth-1
          is_running => find_timer(tmr, stack%name(depth), stack )
          if (.not. is_running%time%running) call start_timer(tmr, (tmr%stack%name(depth)), stack)
       end do
       tmr%stack%depth = depth-1
    end if

    call start_timer(tmr, name, stack)

  End Subroutine start_timer_path

  Subroutine stop_timer_path(tmr, name_in, stop_parents)
    !!------------------------------------------------!
    !!
    !! Stop a timer running on a given timer system
    !! ignoring the timer call stack
    !! If stop_parents stop all running timers
    !! which are on the path -- Default TRUE
    !! - Path should be colon separated
    !!
    !!
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type ( timer_type ), Intent( InOut ) :: tmr
    Character ( Len = * ), Intent( In    )  :: name_in
    Logical, Intent ( In    ), Optional :: stop_parents
    Logical :: parents
    Character ( Len = max_name ) :: name
    Type ( call_stack ) :: stack
    Type ( node ), pointer :: is_running
    Type ( node ), Pointer :: current_timer
    Integer :: depth

    Call timer_split_stack_string(tmr, name_in, stack, name)

    call stop_timer(tmr, name, stack)

    parents = .true.
    if ( present(stop_parents) ) parents = stop_parents

    if (parents) then
       do depth = tmr%stack%depth, 1, -1
          tmr%stack%depth = depth-1
          is_running => find_timer(tmr, tmr%stack%name(depth), stack )
          if (is_running%time%running) call stop_timer(tmr, (tmr%stack%name(depth)), stack)
       end do
    end if

  End Subroutine stop_timer_path

  Subroutine timer_report(tmr)
    !!------------------------------------------------!
    !!
    !! Stop the main timer system and print the full
    !! table of timers to stdout
    !!
    !!
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type( timer_type ), Intent( InOut ) :: tmr
    Type ( node ), Pointer :: current_timer
    Character( Len = 138 ), Dimension(:), Allocatable :: message
    Integer :: nprocs
    Integer :: proc
    Integer :: ierr

    Call stop_timer(tmr, 'Main')

    current_timer => tmr%tree%head%child

    Allocate(message(-2:tmr%tree%n_timers+3), stat=ierr)
    If ( ierr > 0 ) call timer_error(tmr, 'Error allocating message in timer_report')

    Call timer_print_tree(tmr, current_timer, tmr%max_depth, -1, message)
    Call timer_write(message, tmr)

#ifdef mpi
    Call mpi_comm_size(mpi_comm_world, nprocs, ierr)
    if ( ierr /= 0 ) call timer_error(tmr, 'MPI comm size error in timer_report')
    If (tmr%proc_detail) Then
       Do proc = 0, nprocs
          If (tmr%proc_id == proc) Then
             Call timer_print_tree(tmr, current_timer, tmr%max_depth, proc, message)
          End If

          If (proc /= 0) Then
             Call mpi_sendrecv(message, len(message), mpi_character, 0, timer_tag, &
                  & message, len(message), mpi_character, proc, timer_tag, mpi_comm_world, status, ierr)
             if ( ierr /= 0 ) call timer_error(tmr, 'MPI sendrecv error in timer_report')
          End If

          Call timer_write(message, tmr)
          Call mpi_barrier(mpi_comm_world, ierr)
          if ( ierr /= 0 ) call timer_error(tmr, 'MPI barrier error in timer_report')
       End Do
    End If
#endif

  End Subroutine timer_report

  Subroutine timer_print_tree(tmr, init_node, max_depth, proc_id, message)
    !!------------------------------------------------!
    !!
    !! Return a table of the given timer system to
    !! the message variable
    !!
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Implicit None
    Type ( timer_type ), Intent ( In    ) :: tmr
    Type ( node ), Target, Intent ( In     )  :: init_node
    Character( Len = 138 ), Dimension(-2:), Intent(   Out ) :: message
    Type ( node ), Pointer :: current_timer
    Integer, Intent ( In    ) :: proc_id
    Integer, Intent ( In    ) :: max_depth
    Integer :: write_node
    Real ( Kind = dp ) :: total_min, total_max, total_av
    Real ( Kind = dp ) :: call_min, call_max, call_av
    Real ( Kind = dp ) :: sum_timed, total_elapsed
    Integer :: depth, itimer
    Integer :: nprocs
    Character( Len = 8   ) :: proc_string
    Character( Len = 7   ) :: depth_symb

    message(:) = ''

    sum_timed = 0.0_dp

    If ( proc_id < 0 ) Then
       write_node = 0
       proc_string = "All"
    Else
       write_node = proc_id
       Write(proc_string,'(i8.1)') proc_id
    End If

    current_timer => init_node
    total_elapsed = current_timer%time%total

    depth = 0
    itimer = 0

    ! Write table open and header
    Write(message(-2), 100)
    Write(message(-1), 101)

    Do While (depth > -1)

       If (current_timer%time%running) &
            & Call timer_error(tmr, 'Program terminated while timer '// trim(current_timer%time%name)//' still running')

       If ( depth == 0 .and. current_timer%time%name /= "Main") Then
          depth_symb = Repeat('-', 7)
       Else If ( Associated(current_timer%child) ) Then
          depth_symb = Repeat(" ",depth)//"|v"
       Else
          depth_symb = Repeat(" ",depth)//"|-"
       End If

       total_min = current_timer%time%total
       total_max = current_timer%time%total
       total_av  = current_timer%time%total

#ifdef mpi
       Call mpi_comm_size(mpi_comm_world, nprocs, ierr)
       if ( ierr /= 0 ) call timer_error(tmr, 'MPI comm size error in timer_print_tree')
       If ( proc_id < 0 ) Then
          Call mpi_reduce(total_min, total_min, 1, mpi_double_precision, mpi_min, 0, mpi_comm_world, ierr)
          if ( ierr /= 0 ) call timer_error(tmr, 'MPI reduce error in timer_print_tree')
          Call mpi_reduce(total_max, total_max, 1, mpi_double_precision, mpi_max, 0, mpi_comm_world, ierr)
          if ( ierr /= 0 ) call timer_error(tmr, 'MPI reduce error in timer_print_tree')
          Call mpi_reduce(total_av, total_av, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
          if ( ierr /= 0 ) call timer_error(tmr, 'MPI reduce error in timer_print_tree')
          total_av = total_av / nprocs
       End If
#endif

       if (depth == 1 .and. current_timer%parent%time%name == "Main") sum_timed = sum_timed + total_av

       call_min  = current_timer%time%min
       call_max  = current_timer%time%max

#ifdef mpi
       If ( proc_id < 0 ) Then
          Call mpi_reduce(call_min, call_min, 1, mpi_double_precision, mpi_min, 0, mpi_comm_world, ierr)
          if ( ierr /= 0 ) call timer_error(tmr, 'MPI reduce error in timer_print_tree')
          Call mpi_reduce(call_max, call_max, 1, mpi_double_precision, mpi_max, 0, mpi_comm_world, ierr)
          if ( ierr /= 0 ) call timer_error(tmr, 'MPI reduce error in timer_print_tree')
       End If
#endif

       call_av   = total_av/current_timer%time%calls

       Write(message(itimer), 102 ) depth_symb,current_timer%time%name, proc_string, current_timer%time%calls, &
            & call_min, call_max, call_av,  &
            & total_min, total_max, total_av, total_av*100.0_dp/total_elapsed

       If (Associated(current_timer%child) .and. depth < max_depth ) Then
          current_timer => current_timer%child
          depth = depth + 1
       Else If (Associated(current_timer%next_sibling)) Then
          current_timer => current_timer%next_sibling
       Else If (Associated(current_timer%parent)) Then ! Recurse back up
          Do While (Associated(current_timer%parent))
             current_timer => current_timer%parent
             depth = depth - 1
             If (associated(current_timer%next_sibling)) Then
                current_timer => current_timer%next_sibling
                Exit
             End If
          End Do
       Else
          Exit
       End If
       itimer = itimer + 1

    End Do

    Write(message(itimer), 102) Repeat('-',7), "Untimed           ", proc_string, 0, &
         & 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         & total_elapsed - sum_timed , 100.0_dp - sum_timed*100.0_dp/total_elapsed
    Write(message(itimer+1), 100)
    Write(message(itimer+2), *) ''

100 Format(1X,"+",28("-"),2("+",10("-")),7("+",11("-")),"+")
101 Format(1X,"|",12X,"Name",12X,"| Process  ","|  Calls   ","| Call Min  ","| Call Max  ",&
         & "| Call Ave  ","|  Tot Min  ","|  Tot Max  ","|  Tot Ave  ","|     %     ","|")
102 Format(1X,"|",1X,A7,1X,A18,1X,"|",1X,A8,1X,"|",1X,I8,1X,"|",1X,F9.4,1X,"|",1X,F9.4,1X,&
         & "|",1X,F9.4,1X,"|",1X,F9.4,1X,"|",1X,F9.4,1X,"|",1X,F9.4,1X,"|",2X,F8.4,1X,"|")

  End Subroutine timer_print_tree

  Subroutine init_timer(current_timer, name)
    !!------------------------------------------------!
    !!
    !! Initialise a node timer
    !!
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type ( node_timer ) :: current_timer
    Character ( Len = * )  :: name

    current_timer%name  = name
    current_timer%calls = 0
    current_timer%max   = -1.0_dp
    current_timer%min   = huge(1.0_dp)
    current_timer%total = 0.0_dp
    current_timer%last  = huge(1.0_dp)

  End Subroutine init_timer

  Subroutine timer_last_time(tmr, name, screen)
    !!------------------------------------------------!
    !!
    !! Write the length of the previous call to a given
    !! node timer
    !!
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type ( timer_type ), Intent ( InOut ) :: tmr
    Character ( Len = * )  :: name
    Logical, Optional :: screen
    Logical :: to_screen
    Character ( Len = 72 ) :: message
    Type ( node ), pointer :: current_timer

    to_screen = .false.
    If (Present(screen)) to_screen = screen

    current_timer => find_timer(tmr, name)
    If (to_screen) Then
       Write(0,*) current_timer%time%name, current_timer%time%calls, current_timer%time%last
    Else
       Write(message,'(a,2(1X,i0))') current_timer%time%name, current_timer%time%calls, current_timer%time%last
       Call timer_write(message, tmr)
    End If
  End Subroutine timer_last_time

  Subroutine timer_write_mul(message, timer_in)
    !!------------------------------------------------!
    !!
    !! Write multiple lines to standard out
    !!
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Character( Len = * ), Dimension(:), Intent( In    ) :: message
    Type ( timer_type ), Optional, Intent ( In    ):: timer_in
    Integer :: i

    If (Present(timer_in)) Then

       If ( timer_in%proc_id == 0 ) Then
          Do i = 1, Size(message)
             Write(timer_in%out_unit,*) message(i)
          End Do
       End If

    Else

       Do i = 1, Size(message)
          Write(0,*) message(i)
       End Do

    End If

  End Subroutine timer_write_mul

  Subroutine timer_write_sing(message, timer_in)
    !!------------------------------------------------!
    !!
    !! Write a single line to standard out
    !!
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Character( Len = * ), Intent( In    ) :: message
    Type ( timer_type ), Optional, Intent ( In    ):: timer_in

    If (Present(timer_in)) Then

       If ( timer_in%proc_id == 0 ) Then
          Write(timer_in%out_unit,*) message
       End If

    Else

       Write(0,*) message

    End If

  End Subroutine timer_write_sing

  Subroutine timer_error(tmr, message)
    !!------------------------------------------------!
    !!
    !! Report an error with the timer system
    !!
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type ( timer_type ), Intent ( In    ):: tmr
    Character ( len = * ), Intent( In    ) :: message

    Call timer_write('')
    Call timer_write(message)
    Call timer_write('')
    Call dump_call_stack(tmr%stack)

    Stop
  End Subroutine timer_error

  Subroutine timer_split_stack_string(tmr, stack_string, newStack, name)
    !!------------------------------------------------!
    !!
    !! Split a colon-separated path string into a stack
    !!
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type ( timer_type ), Intent ( In    ):: tmr
    Character( Len = * ), Intent( In    ) :: stack_string
    Character( Len = 256 ) :: stack
    Type ( call_stack ), Intent(   Out ) :: newStack
    Character( Len = max_name ), Intent(   Out ) :: name
    Integer :: i
    Integer :: cnt

    stack = Adjustl(stack_string)
    cnt = 1
    Do i = 1, Len(stack)
       If (stack(i:i) == ":") cnt = cnt + 1
    End Do

    If (cnt > max_depth) Call timer_error(tmr, 'Stack depth greater than max depth in timer_split_stack')

    newStack%depth = 0
    Do While (Index(stack,':') > 0)
       newStack%depth = newStack%depth + 1
       newStack%name(newStack%depth) = Trim(stack(1:index(stack,':')-1))

       stack(1:index(stack,':')) = " "
       stack = adjustl(stack)

    End Do

    name = Trim(stack)

  End Subroutine timer_split_stack_string

  Subroutine get_time(out)
    Real(kind=dp), Intent(   Out ) :: out

#ifdef mpi
    call mpi_get_wtime(out)
#else
    call cpu_time(out)
#endif

  end Subroutine get_time

End Module timer
