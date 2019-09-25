subroutine do_react_vode(e_out,rpar,dt) &
     bind(C, name='do_react_vode')

      use amrex_constants_module, only : rt => amrex_real

      implicit none

      real(rt), intent(inout) :: e_out
      real(rt), intent(inout) :: rpar(4)
      real(rt), intent(in), value :: dt

      real(rt) :: rpar_init(4)
      real(rt) :: e_in

      rpar_init=rpar
      e_in=e_out

       call vode_wrapper(dt,rpar_init(3),rpar_init(1),rpar_init(2),e_in, &
         rpar(1),rpar(2) ,e_out)

end subroutine do_react_vode


subroutine vode_wrapper(dt, rho_in, T_in, ne_in, e_in, T_out, ne_out, e_out)

    use amrex_error_module, only : amrex_error
    use amrex_fort_module, only : rt => amrex_real
    use vode_aux_module, only: rho_vode, T_vode, ne_vode, &
                               i_vode, j_vode, k_vode

    implicit none

    real(rt), intent(in   ) :: dt
    real(rt), intent(in   ) :: rho_in, T_in, ne_in, e_in
    real(rt), intent(  out) ::         T_out,ne_out,e_out

    ! Set the number of independent variables -- this should be just "e"
    integer, parameter :: NEQ = 1
  
    ! Allocate storage for the input state
    real(rt) :: y(NEQ)

    ! Our problem is stiff, tell ODEPACK that. 21 means stiff, jacobian 
    ! function is supplied, 22 means stiff, figure out my jacobian through 
    ! differencing
    integer, parameter :: MF_ANALYTIC_JAC = 21, MF_NUMERICAL_JAC = 22

    ! Tolerance parameters:
    !
    !  itol specifies whether to use an single absolute tolerance for
    !  all variables (1), or to pass an array of absolute tolerances, one
    !  for each variable with a scalar relative tol (2), a scalar absolute
    !  and array of relative tolerances (3), or arrays for both (4)
    !  
    !  The error is determined as e(i) = rtol*abs(y(i)) + atol, and must
    !  be > 0.  
    !
    ! We will use arrays for both the absolute and relative tolerances, 
    ! since we want to be easier on the temperature than the species

    integer, parameter :: ITOL = 1
    real(rt) :: atol(NEQ), rtol(NEQ)
    
    ! We want to do a normal computation, and get the output values of y(t)
    ! after stepping though dt
    integer, PARAMETER :: ITASK = 1
  
    ! istate determines the state of the calculation.  A value of 1 meeans
    ! this is the first call to the problem -- this is what we will want.
    ! Note, istate is changed over the course of the calculation, so it
    ! cannot be a parameter
    integer :: istate

    ! we will override the maximum number of steps, so turn on the 
    ! optional arguments flag
    integer, parameter :: IOPT = 1
    
    ! declare a real work array of size 22 + 9*NEQ + 2*NEQ**2 and an
    ! integer work array of since 30 + NEQ

    integer, parameter :: LRW = 22 + 9*NEQ + 2*NEQ**2
    real(rt)   :: rwork(LRW)
    real(rt)   :: time
    ! real(rt)   :: dt4
    
    integer, parameter :: LIW = 30 + NEQ
    integer, dimension(LIW) :: iwork
    
    real(rt) :: rpar
    integer          :: ipar

    EXTERNAL jac, f_rhs
    
    logical, save :: firstCall = .true.

    T_vode   = T_in
    ne_vode  = ne_in
    rho_vode = rho_in

    ! We want VODE to re-initialize each time we call it
    istate = 1
    
    rwork(:) = 0.d0
    iwork(:) = 0
    
    ! Set the maximum number of steps allowed (the VODE default is 500)
    iwork(6) = 2000
    
    ! Initialize the integration time
    time = 0.d0
    
    ! We will integrate "e" in time. 
    y(1) = e_in

    ! Set the tolerances.  
    atol(1) = 1.d-4 * e_in
    rtol(1) = 1.d-4

    ! call the integration routine
    call dvode(f_rhs, NEQ, y, time, dt, ITOL, rtol, atol, ITASK, &
               istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_NUMERICAL_JAC, &
               rpar, ipar)

    e_out  = y(1)
    T_out  = T_vode
    ne_out = ne_vode

    if (istate < 0) then
       print *, 'istate = ', istate, 'at (i,j,k) ',i_vode,j_vode,k_vode
       call amrex_error("ERROR in vode_wrapper: integration failed")
    endif

!      print *,'Calling vode with 1/4 the time step'
!      dt4 = 0.25d0  * dt
!      y(1) = e_in

!      do n = 1,4
!         call dvode(f_rhs, NEQ, y, time, dt4, ITOL, rtol, atol, ITASK, &
!                    istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_NUMERICAL_JAC, &
!                    rpar, ipar)
!         if (istate < 0) then
!            print *, 'doing subiteration ',n
!            print *, 'istate = ', istate, 'at (i,j,k) ',i,j,k
!            call amrex_error("ERROR in vode_wrapper: sub-integration failed")
!         end if

!      end do
!   endif

end subroutine vode_wrapper
