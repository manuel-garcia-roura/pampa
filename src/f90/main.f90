program main
   
   implicit none
   
   ! Local variables:
   logical :: steady_state = .false.
   
   ! Run the calculation:
   if (.not. steady_state) then
      call run()
   else
      call run_steady_state()
   end if
   
   contains
   
   ! Run the calculation:
   subroutine run()
      
      use iso_c_binding
      use pcm
      
      implicit none
      
      ! Local variables:
      real(c_double), pointer :: dt(:)
      integer(c_int) :: ndt
      integer(c_int) :: n
      real(c_double) :: dtn
      real(c_double) :: t
      integer(c_int) :: error
      
      ! Initialize the calculation:
      call pcm_initialize(dt, ndt, error)
      if (error > 0) then; write(6, '(A)') "Error in pcm_initialize()."; call exit(1); end if
      
      ! Run the time-stepping loop:
      t = 0.0
      do n = 0, ndt
         if (n == 0) then
            dtn = 0.0
         else
            dtn = dt(n)
         end if
         t = t + dtn
         call pcm_solve(n, dtn, t, error)
         if (error > 0) then; write(6, '(A)') "Error in pcm_solve()."; call exit(1); end if
      end do
      
      ! Finalize the calculation:
      call pcm_finalize(dt, error)
      if (error > 0) then; write(6, '(A)') "Error in pcm_finalize()."; call exit(1); end if
      
   end subroutine run
   
   ! Run a steady-state calculation:
   subroutine run_steady_state()
      
      use iso_c_binding
      use pcm
      
      implicit none
      
      ! Local variables:
      integer(c_int) :: error
      
      ! Initialize the steady-state calculation:
      call pcm_initialize_steady_state(error)
      if (error > 0) then; write(6, '(A)') "Error in pcm_initialize_steady_state()."; call exit(1); end if
      
      ! Get the steady-state solution:
      call pcm_solve_steady_state(error)
      if (error > 0) then; write(6, '(A)') "Error in pcm_solve_steady_state()."; call exit(1); end if
      
      ! Finalize the steady-state calculation:
      call pcm_finalize_steady_state(error)
      if (error > 0) then; write(6, '(A)') "Error in pcm_finalize_steady_state()."; call exit(1); end if
      
   end subroutine run_steady_state
   
end program main
