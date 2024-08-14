! The pcm (Pela's conduction model) module:
module pcm
   
   implicit none
   
   ! Initialize the calculation (C function):
   interface
      subroutine pampa_initialize(argc, argv, dt, ndt, error) bind(C, name = "pampa_initialize")
         use iso_c_binding
         implicit none
         integer(c_int), value, intent(in) :: argc
         type(c_ptr), value, intent(in) :: argv
         real(c_double), pointer, intent(inout) :: dt(:)
         integer(c_int), intent(out) :: ndt
         integer(c_int), intent(out) :: error
      end subroutine pampa_initialize
   end interface
   
   ! Get the solution (C function):
   interface
      subroutine pampa_solve(n, dtn, t, error) bind(C, name = "pampa_solve")
         use iso_c_binding
         implicit none
         integer(c_int), value, intent(in) :: n
         real(c_double), value, intent(in) :: dtn
         real(c_double), value, intent(in) :: t
         integer(c_int), intent(out) :: error
      end subroutine pampa_solve
   end interface
   
   ! Finalize the calculation (C function):
   interface
      subroutine pampa_finalize(dt, error) bind(C, name = "pampa_finalize")
         use iso_c_binding
         implicit none
         real(c_double), pointer, intent(inout) :: dt(:)
         integer(c_int), intent(out) :: error
      end subroutine pampa_finalize
   end interface
   
   contains
   
   ! Initialize the calculation:
   subroutine pcm_initialize(dt, ndt, error)
      
      use iso_c_binding
      
      implicit none
      
      ! Arguments:
      real(c_double), pointer, intent(inout) :: dt(:)
      integer(c_int), intent(out) :: ndt
      integer(c_int), intent(out) :: error
      
      ! Local variables:
      integer(c_int) :: argc
      type(c_ptr), target, allocatable :: argv(:)
      character(len = 256), target, allocatable :: args(:)
      integer :: i
      
      ! Get the number of command line arguments:
      argc = command_argument_count() + 1
      
      ! Allocate memory for argv (C pointers) and args (Fortran90 strings):
      allocate(argv(argc))
      allocate(args(argc))
      
      ! Set the first argument to the program name like in C/C++:
      call get_command_argument(0, args(1))
      args(1) = trim(adjustl(args(1))) // c_null_char
      argv(1) = c_loc(args(1))
      
      ! Populate the argv array with C pointers to the Fortran90 strings:
      do i = 2, argc
         call get_command_argument(i-1, args(i))
         args(i) = trim(adjustl(args(i))) // c_null_char
         argv(i) = c_loc(args(i))
      end do
      
      ! Call the C function to initialize the calculation:
      call pampa_initialize(argc, c_loc(argv(1)), dt, ndt, error)
      if (error > 0) then; write(6, '(A)') "Error in pampa_initialize()."; return; end if
      
   end subroutine pcm_initialize
   
   ! Get the solution:
   subroutine pcm_solve(n, dtn, t, error)
      
      use iso_c_binding
      
      implicit none
      
      ! Arguments:
      integer(c_int), value, intent(in) :: n
      real(c_double), value, intent(in) :: dtn
      real(c_double), value, intent(in) :: t
      integer(c_int), intent(out) :: error
      
      ! Call the C function to get the solution:
      call pampa_solve(n, dtn, t, error)
      if (error > 0) then; write(6, '(A)') "Error in pampa_solve()."; return; end if
      
   end subroutine pcm_solve
   
   ! Finalize the calculation:
   subroutine pcm_finalize(dt, error)
      
      use iso_c_binding
      
      implicit none
      
      ! Arguments:
      real(c_double), pointer, intent(inout) :: dt(:)
      integer(c_int), intent(out) :: error
      
      ! Call the C function to finalize the calculation:
      call pampa_finalize(dt, error)
      if (error > 0) then; write(6, '(A)') "Error in pampa_finalize()."; return; end if
      
   end subroutine pcm_finalize
   
end module pcm
