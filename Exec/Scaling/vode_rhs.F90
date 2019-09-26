module vode_rhs_module

contains
  
  ! The f_rhs routine provides the right-hand-side for the DVODE solver.
  AMREX_CUDA_FORT_DEVICE subroutine f_rhs(time, y, ydot, rpar)

    use, intrinsic :: iso_c_binding
    use cuvode_parameters_module, only: VODE_NEQS
!    use vode_rpar_indices, only: n_rpar_comps
    use amrex_fort_module, only: rt => amrex_real
    use f_kernel_rhs_dev

    implicit none

    real(rt), intent(INOUT) :: time, y(VODE_NEQS)
!    real(rt), intent(INOUT) :: rpar(n_rpar_comps)
    real(rt), intent(INOUT) :: rpar(4)
    real(rt), intent(INOUT) :: ydot(VODE_NEQS)

#ifdef AMREX_USE_CUDA
    attributes(managed) :: time, y, rpar, ydot
#endif

    call f_rhs_rpar(time, y, ydot, rpar)

  end subroutine f_rhs


  ! Analytical Jacobian
  AMREX_CUDA_FORT_DEVICE subroutine jac(time, y, ml, mu, pd, nrpd, rpar)

    use cuvode_parameters_module, only: VODE_NEQS
!    use vode_rpar_indices, only: n_rpar_comps
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer   , intent(IN   ) :: ml, mu, nrpd
    real(rt), intent(INOUT) :: y(VODE_NEQS), rpar(4), time
    real(rt), intent(  OUT) :: pd(VODE_NEQS,VODE_NEQS)

    !$gpu

    PD(1,1) = -.04D0
    PD(1,2) = 1.D4*Y(3)
    PD(1,3) = 1.D4*Y(2)
    PD(2,1) = .04D0
    PD(2,3) = -PD(1,3)
    PD(3,2) = 6.D7*Y(2)
    PD(2,2) = -PD(1,2) - PD(3,2)

  end subroutine jac

end module vode_rhs_module
