module react_zones_module

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module

  implicit none

contains

  subroutine init_state(lo, hi, &
                        state, s_lo, s_hi, ncomp, npts) bind(C, name="init_state")

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), ncomp)
    integer, intent(in) :: npts, ncomp

    integer :: i, j, k, n

    n = 0
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             state(i, j, k, 1) = ONE - real(n, kind=rt)/real(npts**3, kind=rt)
             state(i, j, k, 2) = ZERO
             state(i, j, k, 3) = ZERO

             n = n + 1
          enddo
       enddo
    enddo

  end subroutine init_state

  AMREX_CUDA_FORT_DEVICE subroutine do_react_cuVODE(e_out,rpar,dt) &
       bind(C, name="do_react_cuVODE")

    use cuvode_parameters_module, only: MF_ANALYTIC_JAC, MF_NUMERICAL_JAC, VODE_LIW
    use cuvode_types_module, only: dvode_t, rwork_t
    use cuvode_module, only: dvode
    use amrex_constants_module, only : rt => amrex_real

    implicit none

    real(rt), intent(inout) :: e_out
    real(rt), intent(inout) :: rpar(4)
    real(rt), intent(in), value :: dt

    real(rt) :: rpar_init(4)
    real(rt) :: e_in

    ! VODE variables
    type (dvode_t) :: dvode_state
    type (rwork_t) :: rwork
    integer :: iwork(VODE_LIW)
    integer :: MF_JAC

    ! istate determines the state of the calculation.  A value of 1 meeans
    ! this is the first call to the problem -- this is what we will want.
    
    integer :: istate
    integer :: ITASK
    integer :: IOPT

#ifdef AMREX_USE_CUDA
    attributes(managed) :: dvode_state, rwork, iwork, MF_JAC, istate, ITASK, IOPT
#endif

    rpar_init=rpar
    e_in=e_out

    ITASK=1
    IOPT=0
    
    ! Use an analytic Jacobian
    MF_JAC = MF_NUMERICAL_JAC

    ! Set the absolute tolerances
    dvode_state % atol(1) = 1.d-4*e_in

    ! Set the relative tolerances
    dvode_state % rtol(1) = 1.d-4

    ! We want VODE to re-initialize each time we call it.
    dvode_state % istate = 1

    ! Initialize work arrays to zero.
    rwork % CONDOPT = ZERO
    rwork % YH   = ZERO
    rwork % WM   = ZERO
    rwork % EWT  = ZERO
    rwork % SAVF = ZERO
    rwork % ACOR = ZERO    
    iwork(:) = 0

    ! Initialize the integration time and set the final time to dt
    dvode_state % T = ZERO
    dvode_state % TOUT = dt

    ! Initialize the initial conditions
    dvode_state % y(1) = e_in

    dvode_state % rpar(1) = rpar(1)
    dvode_state % rpar(2) = rpar(2)
    dvode_state % rpar(3) = rpar(3)
    dvode_state % rpar(4) = rpar(4)

    ! Call the integration routine.
    call dvode(dvode_state, rwork, iwork, ITASK, IOPT, MF_JAC)

    ! Check if the integration failed
    if (dvode_state % istate < 0) then
#ifndef AMREX_USE_CUDA       
       print *, 'ERROR: integration failed'
       print *, 'istate = ', dvode_state % istate
       print *, 'time = ', dvode_state % T
       print *, 'Y start = ', e_in
       print *, 'Y current = ', dvode_state % y
#endif
       stop
    endif
             
    ! Store the final result
    e_out = dvode_state % y(1)

  end subroutine do_react_cuVODE

end module react_zones_module
