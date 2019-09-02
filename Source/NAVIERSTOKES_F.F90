#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <NAVIERSTOKES_F.H>
#include <AMReX_ArrayLim.H>

!=========================================================

      subroutine FORT_SET_NS_PARAMS(n_fluids,            &
                                    dyn_visc_coef,       &
                                    yield_stress,        &
                                    flow_index,          &
                                    reg_param,           &
                                    variable_vel_visc) bind(C,name="set_ns_params")
      implicit none
      integer  i
      integer  n_fluids
      integer  variable_vel_visc
      REAL_T   dyn_visc_coef(n_fluids)
      REAL_T   yield_stress(n_fluids)
      REAL_T   flow_index(n_fluids)
      REAL_T   reg_param

#include <NSCOMM_F.H>

      numflds = n_fluids
      if (.not. ((numflds .eq. 1) .or. (numflds .eq. 2))) then
        write(6,*)'FORT_SET_PARAMS : illegal value of numflds =',numflds
        call bl_abort(" ")
      end if

      do i=1,numflds
        mu_in(i) = dyn_visc_coef(i)
        tau_in(i) = yield_stress(i)
        n_in(i) = flow_index(i)

        if (mu_in(i) .lt. 0) then
          write(6,*)'FORT_SET_PARAMS : illegal value of mu = ',mu_in(i)
          call bl_abort(" ")
        end if

        if (n_in(i) .le. 0) then
          write(6,*)'FORT_SET_PARAMS : illegal value of n = ',n_in(i)
          call bl_abort(" ")
        end if

        if (tau_in(i) .lt. 0) then
          write(6,*)'FORT_SET_PARAMS : illegal value of tau = ',tau_in(i)
          call bl_abort(" ")
        end if

      end do

      eps       = reg_param
      varvisc   = variable_vel_visc

      if (eps .lt. 0) then
        write(6,*)'FORT_SET_PARAMS : illegal value of eps = ',eps
        call bl_abort(" ")
      end if

      if (.not. ((varvisc .eq. 0) .or. (varvisc .eq. 1))) then
        write(6,*)'FORT_SET_PARAMS : illegal value of varvisc =',varvisc
        call bl_abort(" ")
      end if

      end
