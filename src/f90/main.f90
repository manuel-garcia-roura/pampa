program main
   
   use iso_c_binding
   use pcm
   
   implicit none
   
   ! Local variables:
   real(c_double), pointer :: dt(:)
   integer(c_int) :: ndt
   integer(c_int) :: n
   real(c_double) :: dtn
   real(c_double) :: t
   real(c_double) :: v(257040)
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
   
end program main
