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
      integer :: ixlo, ixhi, iylo, iyhi, izlo, izhi


      force = zero

      
!      if(scomp.eq.0) then
!         ixlo = max(f_lo(1),v_lo(1))
!         ixhi = min(f_hi(1),v_hi(1))
!         iylo = max(f_lo(2),v_lo(2))
!         iyhi = min(f_hi(2),v_hi(2))
!         izlo = max(f_lo(3),v_lo(3))
!         izhi = min(f_hi(3),v_hi(3))
!      else         
!         ixlo = max(f_lo(1),s_lo(1))
!         ixhi = min(f_hi(1),s_hi(1))
!         iylo = max(f_lo(2),s_lo(2))
!         iyhi = min(f_hi(2),s_hi(2))
!         izlo = max(f_lo(3),s_lo(3))
!         izhi = min(f_hi(3),s_hi(3))         
!      endif

         ixlo = f_lo(1)
         ixhi = f_hi(1)
         iylo = f_lo(2)
         iyhi = f_hi(2)
         izlo = f_lo(3)
         izhi = f_hi(3)

!         print*, " HERE I AM (1)!"

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



      do k = izlo, izhi
         do j = iylo, iyhi
            do i = ixlo, ixhi
               
!         print*, " HERE I AM (2)!"               
               ! static rho for now
               rho = scal(i,j,k,0)
               temp = scal(i,j,k,2) ! need to make these slots general
               !print*, "MAKEFORCE//rho: ", rho
               !print*, "MAKEFORCE//temp: ", temp

               !         print*, " HERE I AM (3)!"

               if(isnan(rho)) print*, " RHO NAN at:", i,j,k
               if(isnan(temp)) print*, " TEMP NAN at:", i,j,k

               

               ! this is terrible
               if(scomp.eq.0) then ! velocity

               if(isnan(force(i,j,k,0))) print*, " 1. F1 NAN at:", i,j,k
               if(isnan(force(i,j,k,1))) print*, " 1. F2 NAN at:", i,j,k
               if(isnan(force(i,j,k,2))) print*, " 1. F3 NAN at:", i,j,k
               if(isnan(rhs(i,j,k,0))) print*, " 1. RHS1 NAN at:", i,j,k
               if(isnan(rhs(i,j,k,1))) print*, " 1. RHS2 NAN at:", i,j,k
               if(isnan(rhs(i,j,k,2))) print*, " 1. RHS3 NAN at:", i,j,k                  

!         print*, " HERE I AM (4)!", i, j, k, nscal
                 ! basic forcing term 
                 force(i,j,k,0) = rho * Fx !0.d0 ! pass a read-in generic Fi
                 force(i,j,k,1) = rho * Fy !0.d0
                 force(i,j,k,2) = rho * Fz !0.d0

               if(isnan(force(i,j,k,0))) print*, " 2. F1 NAN at:", i,j,k
               if(isnan(force(i,j,k,1))) print*, " 2. F2 NAN at:", i,j,k
               if(isnan(force(i,j,k,2))) print*, " 2. F3 NAN at:", i,j,k
               if(isnan(rhs(i,j,k,0))) print*, " 2. RHS1 NAN at:", i,j,k
               if(isnan(rhs(i,j,k,1))) print*, " 2. RHS2 NAN at:", i,j,k
               if(isnan(rhs(i,j,k,2))) print*, " 2. RHS3 NAN at:", i,j,k                                   

!         print*, " HERE I AM (4a)!"                 

                 ! Boussinesq term, hardcode to y-dir (1) fix for grav vector later
               force(i,j,k,1) = force(i,j,k,1) + rho * gravity * alpha * (temp-Tref)

               if(isnan(force(i,j,k,0))) print*, " 3. F1 NAN at:", i,j,k
               if(isnan(force(i,j,k,1))) print*, " 3. F2 NAN at:", i,j,k
               if(isnan(force(i,j,k,2))) print*, " 3. F3 NAN at:", i,j,k
               if(isnan(rhs(i,j,k,0))) print*, " 3. RHS1 NAN at:", i,j,k
               if(isnan(rhs(i,j,k,1))) print*, " 3. RHS2 NAN at:", i,j,k
               if(isnan(rhs(i,j,k,2))) print*, " 3. RHS3 NAN at:", i,j,k                                 

!         print*, " HERE I AM (4b)!"          

!             if(abs(rhs(i,j,k,0)).gt.0.0d0 .OR. abs(rhs(i,j,k,2)).gt.0.0d0 .OR. abs(rhs(i,j,k,2)).gt.0.0d0) then 
!                print*, " >>> POINT: ", i,j,k
                 !                print*, "     RHS: ", rhs(i,j,k,0),rhs(i,j,k,1),rhs(i,j,k,2)

!         print*, " HERE I AM (4b0):", rhs(i,j,k,0)
!         print*, " HERE I AM (4b1):", rhs(i,j,k,1)
!         print*, " HERE I AM (4b2):", rhs(i,j,k,2)         
         
                 force(i,j,k,0) = force(i,j,k,0) + rho * rhs(i,j,k,0)
                 force(i,j,k,1) = force(i,j,k,1) + rho * rhs(i,j,k,1)
                 force(i,j,k,2) = force(i,j,k,2) + rho * rhs(i,j,k,2)

!         print*, " HERE I AM (4c)!"                  
                 
!             endif

             ! swap back to use sum boundary?
!             rhs(i,j,k,0) = force(i,j,k,0)
!             rhs(i,j,k,1) = force(i,j,k,1)
!             rhs(i,j,k,2) = force(i,j,k,2)

!             if(rhs(i,j,k,0).gt.0.0d0) then
!                print*, " <***> RHS (vel): ", rhs(i,j,k,0),rhs(i,j,k,1),rhs(i,j,k,2)
                 !             endif

               if(isnan(force(i,j,k,0))) print*, " 4. F1 NAN at:", i,j,k
               if(isnan(force(i,j,k,1))) print*, " 4. F2 NAN at:", i,j,k
               if(isnan(force(i,j,k,2))) print*, " 4. F3 NAN at:", i,j,k
               if(isnan(rhs(i,j,k,0))) print*, " 4. RHS1 NAN at:", i,j,k
               if(isnan(rhs(i,j,k,1))) print*, " 4. RHS2 NAN at:", i,j,k
               if(isnan(rhs(i,j,k,2))) print*, " 4. RHS3 NAN at:", i,j,k                                   

          elseif (scomp .EQ. 3 .AND. ncomp .EQ. 1) then ! rho only

!             print*, " HERE I AM (5)!"
             
             force(i,j,k,scomp) = 0.d0

               if(isnan(force(i,j,k,scomp))) print*, " F NAN at:", i,j,k,scomp
               if(isnan(rhs(i,j,k,0))) print*, " RHS NAN at:", i,j,k

          elseif (scomp .EQ. 3) then ! all scalars

!             print*, " HERE I AM (6)!"             

             force(i,j,k,scomp  ) = 0.d0 ! density(ish)
             force(i,j,k,scomp+1) = 0.d0 ! tracer
             force(i,j,k,scomp+2) = 0.d0 ! temp

             force(i,j,k,scomp+2) = force(i,j,k,scomp+2) + rhs(i,j,k,2)

             if(isnan(force(i,j,k,scomp))) print*, " F NAN at:", i,j,k,scomp
             if(isnan(force(i,j,k,scomp+1))) print*, " F NAN at:", i,j,k,scomp+1
             if(isnan(force(i,j,k,scomp+2))) print*, " F NAN at:", i,j,k,scomp+2             
             if(isnan(rhs(i,j,k,0))) print*, " RHS NAN at:", i,j,k
             if(isnan(rhs(i,j,k,1))) print*, " RHS NAN at:", i,j,k
             if(isnan(rhs(i,j,k,2))) print*, " RHS NAN at:", i,j,k

          elseif (scomp .EQ. 4) then ! tracer only

!             print*, " HERE I AM (7)!"             

             force(i,j,k,scomp) = 0.d0

               if(isnan(force(i,j,k,scomp))) print*, " F NAN at:", i,j,k,scomp
               if(isnan(rhs(i,j,k,0))) print*, " RHS NAN at:", i,j,k

          elseif (scomp .EQ. 5) then ! temp only

!             print*, " HERE I AM (8)!"             

             force(i,j,k,scomp) = 0.d0
             force(i,j,k,scomp) = force(i,j,k,scomp) + rhs(i,j,k,2)
             
!             if(rhs(i,j,k,2).gt.0.0d0) then
!                print*, " <***> RHS (scalars): ", rhs(i,j,k,0),rhs(i,j,k,1),rhs(i,j,k,2)
             !             endif

               if(isnan(force(i,j,k,scomp))) print*, " F NAN at:", i,j,k,scomp
               if(isnan(rhs(i,j,k,0))) print*, " RHS NAN at:", i,j,k

          else ! who knows

!             print*, " HERE I AM (9)!"             

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
