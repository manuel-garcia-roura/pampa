program main
   
   use iso_c_binding
   
   implicit none
   
   ! Declare C function prototypes
   interface
      subroutine pampa_initialize(argc, argv, dt, ndt, error) bind(C, name="pampa_initialize")
         import :: c_int, c_double, c_ptr
         integer(c_int), value :: argc
         type(c_ptr), value :: argv
         real(c_double), pointer :: dt(:)
         integer(c_int) :: ndt
         integer(c_int) :: error
      end subroutine pampa_initialize

      subroutine pampa_solve(n, dtn, t, error) bind(C, name="pampa_solve")
         import :: c_int, c_double
         integer(c_int), value :: n
         real(c_double), value :: dtn, t
         integer(c_int) :: error
      end subroutine pampa_solve

      subroutine pampa_finalize(dt, error) bind(C, name="pampa_finalize")
         import :: c_double, c_int
         real(c_double), pointer :: dt(:)
         integer(c_int) :: error
      end subroutine pampa_finalize
   end interface

   ! Local variables
   real(c_double), pointer :: dt(:)
   integer(c_int) :: ndt
   integer :: n
   real(c_double) :: t, dtn
   integer(c_int) :: error

    integer :: i, argc = 0
    character(len=256), target :: arg
    character(len=256), target, allocatable :: arg_list(:)
    type(c_ptr), target, allocatable :: argv(:)
    type(c_ptr) :: c_argv

    ! Get the number of command line arguments
    argc = command_argument_count() + 1

    ! Allocate memory for argv array of C pointers
    allocate(argv(argc))
    allocate(arg_list(argc))

    ! First argument should be the program name
    call get_command_argument(0, arg_list(1))
    arg_list(1) = trim(adjustl(arg_list(1))) // c_null_char
    argv(1) = c_loc(arg_list(1))

    ! Populate the argv array with pointers to the Fortran strings
    do i = 2, argc
        call get_command_argument(i-1, arg_list(i))
        arg_list(i) = trim(adjustl(arg_list(i))) // c_null_char
        argv(i) = c_loc(arg_list(i))
    end do

   ! Call the C function to initialize the system
   write(*,*) "argc = ", argc
   call pampa_initialize(argc, c_loc(argv(1)), dt, ndt, error)

   ! Time-stepping loop
   t = 0.0d0
   do n = 0, ndt
      if (n == 0) then
         dtn = 0.0d0
      else
         dtn = dt(n)
      end if
      t = t + dtn
      call pampa_solve(n, dtn, t, error)
   end do

   ! Finalize the system
   call pampa_finalize(dt, error)

end program main
