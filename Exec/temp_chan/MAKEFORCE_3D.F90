#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PROB_NS_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 3

module MakeForce_3D_module

  implicit none

  private

  public :: FORT_MAKEFORCE

contains


! UPSTREAM FORM:  
!    subroutine FORT_MAKEFORCE( time, &
!                              force, f_lo, f_hi,&
!                              vel, v_lo, v_hi,&
!                              scal, s_lo, s_hi,&
!                              dx, xlo, xhi, gravity, scomp, ncomp, &
!                              nscal, getForceVerbose ) &
!                              bind(C, name="FORT_MAKEFORCE") 

   !------------------------------------------------------------
   ! This routine add the forcing terms to the momentum equation
   !------------------------------------------------------------
   subroutine FORT_MAKEFORCE( time,               &
   &                          force, f_lo, f_hi,  &
   &                          vel, v_lo, v_hi,    &
   &                          scal, s_lo, s_hi,   &
   &                          rhs, r_lo, r_hi,    &
   &                          dx,xlo,xhi,gravity, &
   &                          Fx,Fy,Fz,           &
   &                          scomp,ncomp, nscal, &
   &                          getForceVerbose )   &
   &          bind(C, name="FORT_MAKEFORCE")

      implicit none

      ! In/Out
      integer :: f_lo(3), f_hi(3)
      integer :: v_lo(3), v_hi(3)
      integer :: s_lo(3), s_hi(3)
      integer :: r_lo(3), r_hi(3)
      integer :: scomp, ncomp
      integer :: nscal, getForceVerbose
      REAL_T  :: time, dx(3)
      REAL_T  :: xlo(3), xhi(3)
      REAL_T  :: gravity, Fx, Fy, Fz, rho, Tref, alpha, temp
      REAL_T, dimension(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),scomp:scomp+ncomp-1) :: force
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),0:SDIM-1)            :: vel
      REAL_T, dimension(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),0:nscal-1)           :: scal
      REAL_T, dimension(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),0:nscal-1)           :: rhs

      integer :: isioproc, n, i, j, k


      force = zero

!      return

      Tref = 323.15  ! MAKE READABLE!
      alpha = 0.0034 ! calc real values later

!      if (getForceVerbose > 0) then
!         call bl_pd_is_ioproc(isioproc)
!         if (isioproc.eq.1) then
!            do n = 1, SDIM
!               write (6,*) "No forcing applied"
!               write (6,*) "Forcing applied", n
!            enddo
!         endif
!      endif

!      print*, "SDIM:", SDIM
!      print*, "NCOMP:", ncomp
!      print*, "SCOMP:", scomp
!      print*, "NSCAL:", nscal
!      print*, "FORCING, Fi:", Fx, Fy, Fz


! From NS_getForce.cpp L45
!ncomp / scomp / forcing
!  3       0       Vel
!  X       0  "all components"
!  1       3       rho
!  X       3    scalars(all)
!  X       4      tracers
!  X       X  particular scalars

!      do k = force_l3, force_h3
!      do j = force_l2, force_h2
!      do i = force_l1, force_h1



      do k = f_lo(3), f_hi(3)
         do j = f_lo(2), f_hi(2)
            do i = f_lo(1), f_hi(1)

               ! static rho for now
               rho = scal(i,j,k,0)
               temp = scal(i,j,k,2) ! need to make these slots general
               !print*, "MAKEFORCE//rho: ", rho
               !print*, "MAKEFORCE//temp: ", temp

               ! this is terrible
               if(scomp.eq.0) then ! velocity

                 ! basic forcing term 
                 force(i,j,k,0) = rho * Fx !0.d0 ! pass a read-in generic Fi
                 force(i,j,k,1) = rho * Fy !0.d0
                 force(i,j,k,2) = rho * Fz !0.d0

                 ! Boussinesq term, hardcode to y-dir (1) fix for grav vector later
                 force(i,j,k,1) = force(i,j,k,1) + rho * gravity * alpha * (temp-Tref)

!             if(abs(rhs(i,j,k,0)).gt.0.0d0 .OR. abs(rhs(i,j,k,2)).gt.0.0d0 .OR. abs(rhs(i,j,k,2)).gt.0.0d0) then 
!                print*, " >>> POINT: ", i,j,k
                 !                print*, "     RHS: ", rhs(i,j,k,0),rhs(i,j,k,1),rhs(i,j,k,2)
                 
                 force(i,j,k,0) = force(i,j,k,0) + rho * rhs(i,j,k,0)
                 force(i,j,k,1) = force(i,j,k,1) + rho * rhs(i,j,k,1)
                 force(i,j,k,2) = force(i,j,k,2) + rho * rhs(i,j,k,2)
                 
!             endif

             ! swap back to use sum boundary?
!             rhs(i,j,k,0) = force(i,j,k,0)
!             rhs(i,j,k,1) = force(i,j,k,1)
!             rhs(i,j,k,2) = force(i,j,k,2)

!             if(rhs(i,j,k,0).gt.0.0d0) then
!                print*, " <***> RHS (vel): ", rhs(i,j,k,0),rhs(i,j,k,1),rhs(i,j,k,2)
!             endif

          elseif (scomp .EQ. 3 .AND. ncomp .EQ. 1) then ! rho only

             force(i,j,k,scomp) = 0.d0

          elseif (scomp .EQ. 3) then ! all scalars

             force(i,j,k,scomp  ) = 0.d0 ! density(ish)
             force(i,j,k,scomp+1) = 0.d0 ! tracer
             force(i,j,k,scomp+2) = 0.d0 ! temp

             force(i,j,k,scomp+2) = force(i,j,k,scomp+2) + rhs(i,j,k,2)

          elseif (scomp .EQ. 4) then ! tracer only

             force(i,j,k,scomp) = 0.d0

          elseif (scomp .EQ. 5) then ! temp only

             force(i,j,k,scomp) = 0.d0
             force(i,j,k,scomp) = force(i,j,k,scomp) + rhs(i,j,k,2)
             
!             if(rhs(i,j,k,2).gt.0.0d0) then
!                print*, " <***> RHS (scalars): ", rhs(i,j,k,0),rhs(i,j,k,1),rhs(i,j,k,2)
!             endif

          else ! who knows

             force(i,j,k,scomp) = 0.d0 

          endif

          ! get gradient at (i,j,k)
!          call get_grad(vel(:,:,:,0),i,j,k,du)
!          call get_grad(vel(:,:,:,1),i,j,k,dv)
!          call get_grad(vel(:,:,:,2),i,j,k,dw)

          ! get viscous diss
!          call get_dissipation(du,dv,dw,epsi)
!          epsi = nu_m * rho * epsi


      end do
      end do
      end do




!      call flush(6)

!      print*, " DONE WITH FORT_MAKEFORCE"

   end subroutine FORT_MAKEFORCE
  
end module MakeForce_3D_module
