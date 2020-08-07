#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <DERIVE_F.H>
#include <PROB_NS_F.H>
#include <AMReX_ArrayLim.H>

module derive_nd_module

  use amrex_error_module, only : amrex_abort

  implicit none

  private

  public :: derdvrho, dermprho, deravgpres, &
            dermgvort, dermgdivu, &
            dergrdpx, dergrdpy, dergrdpz, & 
            derregrad, derregrad_old

contains

!===============================================================
! This file contains functions which compute derived quantities.
! All of the argument lists have the same template, shown below
! 
! INPUTS/OUTPUTS:
! 
! e         <= the quantity derived
! DIMS(e)   => index extent of e array
! nv        => number of components in e array (should be 1)
! dat       => data neded to derive e
! DIMS(dat) => index limits of dat array
! ncomp     => number of components of dat array (3)
! lo,hi     => subrange of e array where result is requested
! domlo,hi  => index extent of problem domain (cell centered)
! delta     => cell spacing
! xlo       => physical location of lower left hand
!         corner of e array
! time      => problem evolution time
! bc        => array of bndry types for component values
!              valid only if component touches bndry
!===============================================================

!=========================================================
!  Compute the amagnitude of the vorticity from the 
!  velocity field
!=========================================================

   subroutine dermgvort (e,   e_lo, e_hi, nv, &
                         dat, d_lo, d_hi, ncomp, &
                         lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                         level, grid_no) &
                         bind(C, name="dermgvort")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(inout), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      REAL_T :: uy, uz, vx, vz, wx, wy, dx, dy, dz
      REAL_T :: uycen, uzcen, uylo, uyhi, uzlo, uzhi
      REAL_T :: vxcen, vzcen, vxlo, vxhi, vzlo, vzhi
      REAL_T :: wxcen, wycen, wxlo, wxhi, wylo, wyhi
      REAL_T :: vorfun

      logical :: fixvlo_x, fixwlo_x, fixvhi_x, fixwhi_x
      logical :: fixulo_y, fixwlo_y, fixuhi_y, fixwhi_y
      logical :: fixulo_z, fixvlo_z, fixuhi_z, fixvhi_z

      integer :: i, j, k

!
!     ::::: some useful macro definitions
!
#     define U(i,j,k) dat(i,j,k,1)
#     define V(i,j,k) dat(i,j,k,2)
#     define W(i,j,k) dat(i,j,k,3)

#     define ULOY bc(2,1,1)
#     define UHIY bc(2,2,1)
#     define ULOZ bc(3,1,1)
#     define UHIZ bc(3,2,1)

#     define VLOX bc(1,1,2)
#     define VHIX bc(1,2,2)
#     define VLOZ bc(3,1,2)
#     define VHIZ bc(3,2,2)

#     define WLOX bc(1,1,3)
#     define WHIX bc(1,2,3)
#     define WLOY bc(2,1,3)
#     define WHIY bc(2,2,3)

!
!     ::::: statement functions that implement stencil
!
      uycen(i,j,k) = half*(U(i,j+1,k)-U(i,j-1,k))/dy
      uylo(i,j,k)  = (U(i,j+1,k)+three*U(i,j,k)-four*U(i,j-1,k))/(three*dy)
      uyhi(i,j,k)  =-(U(i,j-1,k)+three*U(i,j,k)-four*U(i,j+1,k))/(three*dy)

      vxcen(i,j,k) = half*(V(i+1,j,k)-V(i-1,j,k))/dx
      vxlo(i,j,k)  = (V(i+1,j,k)+three*V(i,j,k)-four*V(i-1,j,k))/(three*dx)
      vxhi(i,j,k)  =-(V(i-1,j,k)+three*V(i,j,k)-four*V(i+1,j,k))/(three*dx)

#if ( AMREX_SPACEDIM == 3 )
      uzcen(i,j,k) = half*(U(i,j,k+1)-U(i,j,k-1))/dz
      uzlo(i,j,k)  = (U(i,j,k+1)+three*U(i,j,k)-four*U(i,j,k-1))/(three*dz)
      uzhi(i,j,k)  =-(U(i,j,k-1)+three*U(i,j,k)-four*U(i,j,k+1))/(three*dz)

      vzcen(i,j,k) = half*(V(i,j,k+1)-V(i,j,k-1))/dz
      vzlo(i,j,k)  = (V(i,j,k+1)+three*V(i,j,k)-four*V(i,j,k-1))/(three*dz)
      vzhi(i,j,k)  =-(V(i,j,k-1)+three*V(i,j,k)-four*V(i,j,k+1))/(three*dz)

      wxcen(i,j,k) = half*(W(i+1,j,k)-W(i-1,j,k))/dx
      wxlo(i,j,k)  = (W(i+1,j,k)+three*W(i,j,k)-four*W(i-1,j,k))/(three*dx)
      wxhi(i,j,k)  =-(W(i-1,j,k)+three*W(i,j,k)-four*W(i+1,j,k))/(three*dx)

      wycen(i,j,k) = half*(W(i,j+1,k)-W(i,j-1,k))/dy
      wylo(i,j,k)  = (W(i,j+1,k)+three*W(i,j,k)-four*W(i,j-1,k))/(three*dy)
      wyhi(i,j,k)  =-(W(i,j-1,k)+three*W(i,j,k)-four*W(i,j+1,k))/(three*dy)
#endif

#if ( AMREX_SPACEDIM == 2 )
      vorfun(uy,uz,vx,vz,wx,wy) = vx - uy
#elif ( AMREX_SPACEDIM == 3 )
      vorfun(uy,uz,vx,vz,wx,wy) = sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)
#endif

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

      ! Init all logical tests on BC to false
      fixvlo_x = .FALSE. ; fixwlo_x = .FALSE. ; fixvhi_x = .FALSE. ; fixwhi_x = .FALSE.
      fixulo_y = .FALSE. ; fixwlo_y = .FALSE. ; fixuhi_y = .FALSE. ; fixwhi_y = .FALSE.
      fixulo_z = .FALSE. ; fixvlo_z = .FALSE. ; fixuhi_z = .FALSE. ; fixvhi_z = .FALSE.

      ! Init all vorticity comp. In 2d uz, vz, wx, wy will alway be zero
      uy = 0.0d0
      uz = 0.0d0
      vx = 0.0d0
      vz = 0.0d0
      wx = 0.0d0
      wy = 0.0d0

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               uy = uycen(i,j,k)
               vx = vxcen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wx = wxcen(i,j,k)
               wy = wycen(i,j,k)
#endif
               e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end do

      fixvlo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (VLOX .eq. EXT_DIR .or. VLOX .eq. HOEXTRAP) )
      fixvhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (VHIX .eq. EXT_DIR .or. VHIX .eq. HOEXTRAP) )
#if ( AMREX_SPACEDIM == 3 )
      fixwlo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (WLOX .eq. EXT_DIR .or. WLOX .eq. HOEXTRAP) )
      fixwhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (WHIX .eq. EXT_DIR .or. WHIX .eq. HOEXTRAP) )
#endif

      fixulo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (ULOY .eq. EXT_DIR .or. ULOY .eq. HOEXTRAP) )
      fixuhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (UHIY .eq. EXT_DIR .or. UHIY .eq. HOEXTRAP) )
#if ( AMREX_SPACEDIM == 3 )
      fixwlo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (WLOY .eq. EXT_DIR .or. WLOY .eq. HOEXTRAP) )
      fixwhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (WHIY .eq. EXT_DIR .or. WHIY .eq. HOEXTRAP) )

      fixulo_z = ( (lo(3) .eq. domlo(3)) .and. &
                  (ULOZ .eq. EXT_DIR .or. ULOZ .eq. HOEXTRAP) )
      fixuhi_z = ( (hi(3) .eq. domhi(3)) .and. &
                  (UHIZ .eq. EXT_DIR .or. UHIZ .eq. HOEXTRAP) )
      fixvlo_z = ( (lo(3) .eq. domlo(3)) .and. &
                  (VLOZ .eq. EXT_DIR .or. VLOZ .eq. HOEXTRAP) )
      fixvhi_z = ( (hi(3) .eq. domhi(3)) .and. &
                  (VHIZ .eq. EXT_DIR .or. VHIZ .eq. HOEXTRAP) )
#endif

!
!     First do all the faces
!
      if (fixvlo_x .or. fixwlo_x) then
         i = lo(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
               uy = uycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
               wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
               wy = wycen(i,j,k)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
#endif
               e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fixvhi_x .or. fixwhi_x) then
         i = hi(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
               uy = uycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
               wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
               wy = wycen(i,j,k)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
#endif
               e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fixulo_y .or. fixwlo_y) then
         j = lo(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
               wx = wxcen(i,j,k)
               wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
#endif
               e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fixuhi_y .or. fixwhi_y) then
         j = hi(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
               wx = wxcen(i,j,k)
               wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
#endif
               e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

#if ( AMREX_SPACEDIM == 3 )
      if (fixulo_z .or. fixvlo_z) then
         k = lo(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = uycen(i,j,k)
               wy = wycen(i,j,k)
               uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
               vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
               e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fixuhi_z .or. fixvhi_z) then
         k = hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = uycen(i,j,k)
               wy = wycen(i,j,k)
               uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
               vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
               e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if
#endif

!
!     Next do all the edges
!
      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y)) then
         i = lo(1)
         j = lo(2)
         do k = lo(3),hi(3)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y)) then
         i = hi(1)
         j = lo(2)
         do k = lo(3),hi(3)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y)) then
         i = lo(1)
         j = hi(2)
         do k = lo(3),hi(3)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y)) then
         i = hi(1)
         j = hi(2)
         do k = lo(3),hi(3)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         k = lo(3)
         do j = lo(2),hi(2)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            uy = uycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            wy = wycen(i,j,k)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         k = lo(3)
         do j = lo(2),hi(2)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            uy = uycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            wy = wycen(i,j,k)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         k = hi(3)
         do j = lo(2),hi(2)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            uy = uycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            wy = wycen(i,j,k)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         k = hi(3)
         do j = lo(2),hi(2)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            uy = uycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            wy = wycen(i,j,k)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixulo_y .or. fixwlo_y) .and. (fixulo_z .or. fixvlo_z)) then
         j = lo(2)
         k = lo(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = wxcen(i,j,k)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixuhi_y .or. fixwhi_y) .and. (fixulo_z .or. fixvlo_z)) then
         j = hi(2)
         k = lo(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = wxcen(i,j,k)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixulo_y .or. fixwlo_y) .and. (fixuhi_z .or. fixvhi_z)) then
         j = lo(2)
         k = hi(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = wxcen(i,j,k)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixuhi_y .or. fixwhi_y) .and. (fixuhi_z .or. fixvhi_z)) then
         j = hi(2)
         k = hi(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = wxcen(i,j,k)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

!
!     Finally do all the corners
!
      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         j = lo(2)
         k = lo(3)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         j = lo(2)
         k = lo(3)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         j = hi(2)
         k = lo(3)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         j = hi(2)
         k = lo(3)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         j = lo(2)
         k = hi(3)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         j = lo(2)
         k = hi(3)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         j = hi(2)
         k = hi(3)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         j = hi(2)
         k = hi(3)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

#     undef U
#     undef V
#     undef W
#     undef ULOY
#     undef UHIY
#     undef ULOZ
#     undef UHIZ
#     undef VLOX
#     undef VHIX
#     undef VLOZ
#     undef VHIZ
#     undef WLOX
#     undef WHIX
#     undef WLOY
#     undef WHIY

   end subroutine dermgvort

   
!=========================================================
!  Compute the gradient-based cell Reynolds number (OLD)
!=========================================================
   subroutine derregrad_old (e,   e_lo, e_hi, nv, &
                         dat, d_lo, d_hi, ncomp, &
                         lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
!                         level, grid_no, nu_m) &
                         level, grid_no) &
                         bind(C, name="derregrad_old")

      implicit none

      !  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(inout), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      REAL_T, dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3)) :: r0
      integer, intent(in) :: level, grid_no

      !  Local
      REAL_T :: ux, uy, uz, vx, vy, vz, wx, wy, wz, dx, dy, dz
      REAL_T :: uxcen, uycen, uzcen, uxlo, uxhi, uylo, uyhi, uzlo, uzhi
      REAL_T :: vxcen, vycen, vzcen, vxlo, vxhi, vylo, vyhi, vzlo, vzhi
      REAL_T :: wxcen, wycen, wzcen, wxlo, wxhi, wylo, wyhi, wzlo, wzhi
      REAL_T :: nu_m, re_g, wt1

      logical :: fixulo_x, fixvlo_x, fixwlo_x
      logical :: fixulo_y, fixvlo_y, fixwlo_y
      logical :: fixulo_z, fixvlo_z, fixwlo_z

      logical :: fixuhi_x, fixvhi_x, fixwhi_x
      logical :: fixuhi_y, fixvhi_y, fixwhi_y
      logical :: fixuhi_z, fixvhi_z, fixwhi_z

      integer :: i, j, k

!
!     ::::: some useful macro definitions
!
#     define U(i,j,k) dat(i,j,k,1)
#     define V(i,j,k) dat(i,j,k,2)
#     define W(i,j,k) dat(i,j,k,3)

#     define ULOX bc(1,1,1)
#     define UHIX bc(1,2,1)
#     define ULOY bc(2,1,1)
#     define UHIY bc(2,2,1)
#     define ULOZ bc(3,1,1)
#     define UHIZ bc(3,2,1)

#     define VLOX bc(1,1,2)
#     define VHIX bc(1,2,2)
#     define VLOY bc(2,1,2)
#     define VHIY bc(2,2,2)
#     define VLOZ bc(3,1,2)
#     define VHIZ bc(3,2,2)

#     define WLOX bc(1,1,3)
#     define WHIX bc(1,2,3)
#     define WLOY bc(2,1,3)
#     define WHIY bc(2,2,3)
#     define WLOZ bc(3,1,3)
#     define WHIZ bc(3,2,3)


!
!     ::::: statement functions that implement stencil
!
      uxcen(i,j,k) = half*(U(i+1,j,k)-U(i-1,j,k))/dx
      uxlo(i,j,k)  = (U(i+1,j,k)+three*U(i,j,k)-four*U(i-1,j,k))/(three*dx)
      uxhi(i,j,k)  =-(U(i-1,j,k)+three*U(i,j,k)-four*U(i+1,j,k))/(three*dx)

      uycen(i,j,k) = half*(U(i,j+1,k)-U(i,j-1,k))/dy
      uylo(i,j,k)  = (U(i,j+1,k)+three*U(i,j,k)-four*U(i,j-1,k))/(three*dy)
      uyhi(i,j,k)  =-(U(i,j-1,k)+three*U(i,j,k)-four*U(i,j+1,k))/(three*dy)

      vxcen(i,j,k) = half*(V(i+1,j,k)-V(i-1,j,k))/dx
      vxlo(i,j,k)  = (V(i+1,j,k)+three*V(i,j,k)-four*V(i-1,j,k))/(three*dx)
      vxhi(i,j,k)  =-(V(i-1,j,k)+three*V(i,j,k)-four*V(i+1,j,k))/(three*dx)

      vycen(i,j,k) = half*(V(i,j+1,k)-V(i,j-1,k))/dy
      vylo(i,j,k)  = (V(i,j+1,k)+three*V(i,j,k)-four*V(i,j-1,k))/(three*dy)
      vyhi(i,j,k)  =-(V(i,j-1,k)+three*V(i,j,k)-four*V(i,j+1,k))/(three*dy)

#if ( AMREX_SPACEDIM == 3 )
      uzcen(i,j,k) = half*(U(i,j,k+1)-U(i,j,k-1))/dz
      uzlo(i,j,k)  = (U(i,j,k+1)+three*U(i,j,k)-four*U(i,j,k-1))/(three*dz)
      uzhi(i,j,k)  =-(U(i,j,k-1)+three*U(i,j,k)-four*U(i,j,k+1))/(three*dz)

      vzcen(i,j,k) = half*(V(i,j,k+1)-V(i,j,k-1))/dz
      vzlo(i,j,k)  = (V(i,j,k+1)+three*V(i,j,k)-four*V(i,j,k-1))/(three*dz)
      vzhi(i,j,k)  =-(V(i,j,k-1)+three*V(i,j,k)-four*V(i,j,k+1))/(three*dz)

      wxcen(i,j,k) = half*(W(i+1,j,k)-W(i-1,j,k))/dx
      wxlo(i,j,k)  = (W(i+1,j,k)+three*W(i,j,k)-four*W(i-1,j,k))/(three*dx)
      wxhi(i,j,k)  =-(W(i-1,j,k)+three*W(i,j,k)-four*W(i+1,j,k))/(three*dx)

      wycen(i,j,k) = half*(W(i,j+1,k)-W(i,j-1,k))/dy
      wylo(i,j,k)  = (W(i,j+1,k)+three*W(i,j,k)-four*W(i,j-1,k))/(three*dy)
      wyhi(i,j,k)  =-(W(i,j-1,k)+three*W(i,j,k)-four*W(i,j+1,k))/(three*dy)

      wzcen(i,j,k) = half*(W(i,j,k+1)-W(i,j,k-1))/dz
      wzlo(i,j,k)  = (W(i,j,k+1)+three*W(i,j,k)-four*W(i,j,k-1))/(three*dz)
      wzhi(i,j,k)  =-(W(i,j,k-1)+three*W(i,j,k)-four*W(i,j,k+1))/(three*dz)
#endif

!#if ( AMREX_SPACEDIM == 2 )
!      vorfun(uy,uz,vx,vz,wx,wy) = vx - uy
!#elif ( AMREX_SPACEDIM == 3 )
!      vorfun(uy,uz,vx,vz,wx,wy) = sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)
!#endif


      ! local names
      dx = delta(1)
      dy = delta(2)
      dz = delta(3)
!      nu_m = 1.0e-4 !visc_coef[0] !ns.vel_visc_coef
      nu_m = 3.5e-4 ! FIX THIS HARD CODE
!      print*, "HERE: derregrad nu_m:", nu_m

!      wt1 = 0.9d0
      wt1 = 0.95d0

      ! Init all logical tests on BC to false
      fixulo_x = .FALSE. ; fixvlo_x = .FALSE. ; fixwlo_x = .FALSE.
      fixulo_y = .FALSE. ; fixvlo_y = .FALSE. ; fixwlo_y = .FALSE.
      fixulo_z = .FALSE. ; fixvlo_z = .FALSE. ; fixwlo_z = .FALSE.

      fixuhi_x = .FALSE. ; fixvhi_x = .FALSE. ; fixwhi_x = .FALSE.
      fixuhi_y = .FALSE. ; fixvhi_y = .FALSE. ; fixwhi_y = .FALSE.
      fixuhi_z = .FALSE. ; fixvhi_z = .FALSE. ; fixwhi_z = .FALSE.


      ! Init all grad comp. In 2d uz, vz, wx, wy will alway be zero
      ux = 0.0d0
      uy = 0.0d0
      uz = 0.0d0
      vx = 0.0d0
      vy = 0.0d0
      vz = 0.0d0
      wx = 0.0d0
      wy = 0.0d0
      wz = 0.0d0


      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ux = uxcen(i,j,k)
               uy = uycen(i,j,k)
               vx = vxcen(i,j,k)
               vy = vycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wx = wxcen(i,j,k)
               wy = wycen(i,j,k)
               wz = wzcen(i,j,k)
#endif

               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif

            end do
         end do
      end do



      fixulo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (ULOX .eq. EXT_DIR .or. ULOX .eq. HOEXTRAP) )
      fixuhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (UHIX .eq. EXT_DIR .or. UHIX .eq. HOEXTRAP) )
      fixulo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (ULOY .eq. EXT_DIR .or. ULOY .eq. HOEXTRAP) )
      fixuhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (UHIY .eq. EXT_DIR .or. UHIY .eq. HOEXTRAP) )

      fixvlo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (VLOX .eq. EXT_DIR .or. VLOX .eq. HOEXTRAP) )
      fixvhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (VHIX .eq. EXT_DIR .or. VHIX .eq. HOEXTRAP) )
      fixvlo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (VLOY .eq. EXT_DIR .or. VLOY .eq. HOEXTRAP) )
      fixvhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (VHIY .eq. EXT_DIR .or. VHIY .eq. HOEXTRAP) )

#if ( AMREX_SPACEDIM == 3 )
      fixulo_z = ( (lo(3) .eq. domlo(3)) .and. &
                  (ULOZ .eq. EXT_DIR .or. ULOZ .eq. HOEXTRAP) )
      fixuhi_z = ( (hi(3) .eq. domhi(3)) .and. &
                  (UHIZ .eq. EXT_DIR .or. UHIZ .eq. HOEXTRAP) )

      fixvlo_z = ( (lo(3) .eq. domlo(3)) .and. &
                  (VLOZ .eq. EXT_DIR .or. VLOZ .eq. HOEXTRAP) )
      fixvhi_z = ( (hi(3) .eq. domhi(3)) .and. &
                  (VHIZ .eq. EXT_DIR .or. VHIZ .eq. HOEXTRAP) )

      fixwlo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (WLOX .eq. EXT_DIR .or. WLOX .eq. HOEXTRAP) )
      fixwhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (WHIX .eq. EXT_DIR .or. WHIX .eq. HOEXTRAP) )
      fixwlo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (WLOY .eq. EXT_DIR .or. WLOY .eq. HOEXTRAP) )
      fixwhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (WHIY .eq. EXT_DIR .or. WHIY .eq. HOEXTRAP) )
      fixwlo_z = ( (lo(3) .eq. domlo(3)) .and. &
                  (WLOZ .eq. EXT_DIR .or. WLOZ .eq. HOEXTRAP) )
      fixwhi_z = ( (hi(3) .eq. domhi(3)) .and. &
                  (WHIZ .eq. EXT_DIR .or. WHIZ .eq. HOEXTRAP) )
#endif



      ! 1) First do all the faces...
      if (fixulo_x .or. fixvlo_x .or. fixwlo_x) then
         i = lo(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
               uy = uycen(i,j,k)
               vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
               vy = vycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
               wy = wycen(i,j,k)
               wz = wzcen(i,j,k)
#endif
               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
            end do
         end do
      end if

      if (fixuhi_x .or. fixvhi_x .or. fixwhi_x) then
         i = hi(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
               uy = uycen(i,j,k)
               vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
               vy = vycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
               wy = wycen(i,j,k)
               wz = uzcen(i,j,k)
#endif
               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
            end do
         end do
      end if

      if (fixulo_y .or. fixvlo_y .or. fixwlo_y) then
         j = lo(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               ux = uxcen(i,j,k)
               uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
               vx = vxcen(i,j,k)
               vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
#if ( AMREX_SPACEDIM == 3 )
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wx = wxcen(i,j,k)
               wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
               wz = wzcen(i,j,k)
#endif
               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
            end do
         end do
      end if

      if (fixuhi_y .or. fixvhi_y .or. fixwhi_y) then
         j = hi(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               ux = uxcen(i,j,k)
               uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
               vx = vxcen(i,j,k)
               vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
#if ( AMREX_SPACEDIM == 3 )
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wx = wxcen(i,j,k)
               wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
               wz = wzcen(i,j,k)
#endif
               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
            end do
         end do
      end if

#if ( AMREX_SPACEDIM == 3 )
      if (fixulo_z .or. fixvlo_z .or. fixwlo_z) then
         k = lo(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               ux = uxcen(i,j,k)
               uy = uycen(i,j,k)
               uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)

               vx = vxcen(i,j,k)
               vy = vycen(i,j,k)
               vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)

               wx = wxcen(i,j,k)
               wy = wycen(i,j,k)
               wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)

               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
            end do
         end do
      end if

      if (fixuhi_z .or. fixvhi_z .or. fixwhi_z) then
         k = hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               ux = uxcen(i,j,k)
               uy = uycen(i,j,k)
               uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)

               vx = vxcen(i,j,k)
               vy = vycen(i,j,k)
               vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)

               wx = wxcen(i,j,k)
               wy = wycen(i,j,k)
               wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)

               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
            end do
         end do
      end if
#endif

 
      ! 2) Next do all the edges...
      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixvlo_y .or. fixwlo_y)) then ! lo x-y edge
         i = lo(1)
         j = lo(2)
         do k = lo(3),hi(3)
            ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
#if ( AMREX_SPACEDIM == 3 )
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            wz = wzcen(i,j,k)
#endif
               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
         end do
      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixvlo_y .or. fixwlo_y)) then ! hi/lo x-y edge
         i = hi(1)
         j = lo(2)
         do k = lo(3),hi(3)
            ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
#if ( AMREX_SPACEDIM == 3 )
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            wz = wzcen(i,j,k)
#endif
               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
         end do
      end if

      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixvhi_y .or. fixwhi_y)) then
         i = lo(1)
         j = hi(2)
         do k = lo(3),hi(3)
            ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
#if ( AMREX_SPACEDIM == 3 )
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            wz = wzcen(i,j,k)
#endif
               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
         end do
      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixvhi_y .or. fixwhi_y)) then
         i = hi(1)
         j = hi(2)
         do k = lo(3),hi(3)
            ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
#if ( AMREX_SPACEDIM == 3 )
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            wz = wzcen(i,j,k)
#endif
               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
         end do
      end if

      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         i = lo(1)
         k = lo(3)
         do j = lo(2),hi(2)
            ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
            uy = uycen(i,j,k)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            vy = vycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            wy = wycen(i,j,k)
            wx = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
#endif
               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
         end do
      end if

      if ((fixvhi_x .or. fixuhi_x .or. fixwhi_x) .and. (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         i = hi(1)
         k = lo(3)
         do j = lo(2),hi(2)
            ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
            uy = uycen(i,j,k)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            vy = vycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            wy = wycen(i,j,k)
            wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
#endif
               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
         end do
      end if

      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         i = lo(1)
         k = hi(3)
         do j = lo(2),hi(2)
            ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
            uy = uycen(i,j,k)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            vy = uycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            wy = wycen(i,j,k)
            wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
#endif
               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
         end do
      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         i = hi(1)
         k = hi(3)
         do j = lo(2),hi(2)
            ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
            uy = uycen(i,j,k)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            vy = vycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            wy = wycen(i,j,k)
            uz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
#endif
               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
         end do
      end if

      if ((fixulo_x .or. fixvlo_y .or. fixwlo_y) .and. (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         j = lo(2)
         k = lo(3)
         do i = lo(1),hi(1)
            ux = uxcen(i,j,k)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            vx = vxcen(i,j,k)
            vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
#if ( AMREX_SPACEDIM == 3 )
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            wx = wxcen(i,j,k)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
#endif
               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
         end do
      end if

      if ((fixuhi_y .or. fixvhi_y .or. fixwhi_y) .and. (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         j = hi(2)
         k = lo(3)
         do i = lo(1),hi(1)
            ux = uxcen(i,j,k)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            vx = vxcen(i,j,k)
            vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
#if ( AMREX_SPACEDIM == 3 )
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            wx = wxcen(i,j,k)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
#endif
               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
         end do
      end if

      if ((fixulo_y .or. fixvlo_y .or. fixwlo_y) .and. (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         j = lo(2)
         k = hi(3)
         do i = lo(1),hi(1)
            ux = uxcen(i,j,k)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            vx = vxcen(i,j,k)
            vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
#if ( AMREX_SPACEDIM == 3 )
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            wx = wxcen(i,j,k)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
#endif
               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
         end do
      end if

      if ((fixuhi_y .or. fixvhi_y .or. fixwhi_y) .and. (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         j = hi(2)
         k = hi(3)
         do i = lo(1),hi(1)
            ux = uxcen(i,j,k)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            vx = vxcen(i,j,k)
            vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
#if ( AMREX_SPACEDIM == 3 )
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            wx = wxcen(i,j,k)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
#endif
               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
         end do
      end if



      ! 3) Finally do all the corners...
      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. &
          (fixulo_y .or. fixvlo_y .or. fixwlo_y) .and. &
          (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         i = lo(1)
         j = lo(2)
         k = lo(3)
         ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
#if ( AMREX_SPACEDIM == 3 )
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
#endif
               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. &
          (fixulo_y .or. fixvlo_y .or. fixwlo_y) .and. &
          (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         i = hi(1)
         j = lo(2)
         k = lo(3)
         ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
#if ( AMREX_SPACEDIM == 3 )
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
#endif
               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
      end if

      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. &
          (fixuhi_y .or. fixvhi_y .or. fixwhi_y) .and. &
          (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         i = lo(1)
         j = hi(2)
         k = lo(3)
         ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
#if ( AMREX_SPACEDIM == 3 )
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
#endif
               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. &
          (fixuhi_y .or. fixvhi_y .or. fixwhi_y) .and. &
          (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         i = hi(1)
         j = hi(2)
         k = lo(3)
         ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
#if ( AMREX_SPACEDIM == 3 )
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
#endif
               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
      end if

      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. &
          (fixulo_y .or. fixvlo_y .or. fixwlo_y) .and. &
          (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         i = lo(1)
         j = lo(2)
         k = hi(3)
         ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
#if ( AMREX_SPACEDIM == 3 )
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
#endif

               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. &
          (fixulo_y .or. fixvlo_y .or. fixwlo_y) .and. &
          (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         i = hi(1)
         j = lo(2)
         k = hi(3)
         ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
#if ( AMREX_SPACEDIM == 3 )
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
#endif

               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif
      end if

      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. &
          (fixuhi_y .or. fixvhi_y .or. fixwhi_y) .and. &
          (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         i = lo(1)
         j = hi(2)
         k = hi(3)
         ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
#if ( AMREX_SPACEDIM == 3 )
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
#endif

               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif

      end if


      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. &
          (fixuhi_y .or. fixvhi_y .or. fixwhi_y) .and. &
          (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         i = hi(1)
         j = hi(2)
         k = hi(3)
         ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
#if ( AMREX_SPACEDIM == 3 )
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
#endif

               call get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu_m,re_g)
!               if(time .LE. dt) then
                  r0(i,j,k) = re_g
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
!               endif

      end if



      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ux = uxcen(i,j,k)
               uy = uycen(i,j,k)
               vx = vxcen(i,j,k)
               vy = vycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wx = wxcen(i,j,k)
               wy = wycen(i,j,k)
               wz = wzcen(i,j,k)
#endif

               ! MUST FIX CONTAINER
!               if(time .LE. dt) then
                  e(i,j,k,1) = r0(i,j,k)
!               else
!                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*r0(i,j,k)
!               endif

            end do
         end do
      end do




#     undef U
#     undef V
#     undef W
#     undef ULOY
#     undef UHIY
#     undef ULOZ
#     undef UHIZ
#     undef VLOX
#     undef VHIX
#     undef VLOZ
#     undef VHIZ
#     undef WLOX
#     undef WHIX
#     undef WLOY
#     undef WHIY

!   print*, "   DONE WITH derregrad"

   end subroutine derregrad_old


!=========================================================
! Cleaned up version
!=========================================================
   subroutine derregrad (e,   e_lo, e_hi, nv, &
                         dat, d_lo, d_hi, ncomp, &
                         lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
!                         level, grid_no, nu_m) &
                         level, grid_no) &
                         bind(C, name="derregrad")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(inout), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      REAL_T :: nu_m, re_g, wt1, Gij(3,3)
      REAL_T :: dx, dy, dz

      integer :: i, j, k
      integer :: ii, jj, kk


      ! local names
      dx = delta(1)
      dy = delta(2)
      dz = delta(3)
!      nu_m = 1.0e-4 !visc_coef[0] !ns.vel_visc_coef
      nu_m = 3.5e-4 ! FIX THIS HARD CODE
!      print*, "HERE: derregrad nu_m:", nu_m

      wt1 = 0.9d0

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ii = i
               jj = j
               kk = k

               call get_grad (nv, dat, d_lo, d_hi, ncomp, lo, hi, domlo, domhi, delta, &
                              xlo, time, dt, bc, level, grid_no, ii, jj, kk, Gij)

               call get_reg(Gij,dx,dy,dz,nu_m,re_g)

               if(time .LE. dt) then
                  e(i,j,k,1) = re_g
               else
                  e(i,j,k,1) = wt1*e(i,j,k,1) + (1.0d0-wt1)*re_g
               endif

            end do
         end do
      end do


!   print*, "   DONE WITH derregrad"

   end subroutine derregrad



!=========================================================
! Actual Re_G comp
!=========================================================
   subroutine get_reg(duij,dx,dy,dz,nu,re_g)
      implicit none

      integer :: mm, nn, oo
      REAL_T :: dx, dy, dz
      REAL_T :: nu, re_g
      REAL_T, dimension(3,3) :: duij, Mij, Gij, r2

      ! Have full gradient tensor...
      ! calculate: Re_G = (tr[M^2*G*M^2])^1/2/nu, G=u_{k,i}u_{k,j}


      Mij(1:3,1:3) = 0.0d0
      Mij(1,1) = dx*dx
      Mij(2,2) = dy*dy
      Mij(3,3) = dz*dz

      r2(1:3,1:3) = 0.0d0      
      do mm = 1,3
         do nn = 1,3
            do oo = 1,3
               r2(mm,nn) = r2(mm,nn) + duij(mm,oo)*Mij(oo,nn)
            enddo
         enddo
      enddo

      Gij(1:3,1:3) = 0.0d0      
      do mm = 1,3
         do nn = 1,3
            do oo = 1,3
               Gij(mm,nn) = Gij(mm,nn) + r2(oo,mm)*r2(oo,nn)
            enddo
         enddo
      enddo

      ! magnitude might be better
      re_g = sqrt(Gij(1,1)+Gij(2,2)+Gij(3,3)) / nu


   end subroutine get_reg


!=========================================================
! Actual Re_G comp
!=========================================================
   subroutine get_reg_old(ux,uy,uz,vx,vy,vz,wx,wy,wz,dx,dy,dz,nu,re_g)
      implicit none
      REAL_T :: re_g
      integer :: mm, nn, oo
      REAL_T :: ux, uy, uz
      REAL_T :: vx, vy, vz
      REAL_T :: wx, wy, wz
      REAL_T :: dx, dy, dz
      REAL_T :: nu
      REAL_T, dimension(3,3) :: duij, Mij, Gij, r2

      ! Have full gradient tensor...
      ! calculate: Re_G = (tr[M^2*G*M^2])^1/2/nu, G=u_{k,i}u_{k,j}

      duij(1,1) = ux
      duij(1,2) = uy
      duij(1,3) = uz
      duij(2,1) = vx
      duij(2,2) = vy
      duij(2,3) = vz
      duij(3,1) = wx
      duij(3,2) = wy
      duij(3,3) = wz

      Mij(1:3,1:3) = 0.0d0
      Mij(1,1) = dx*dx
      Mij(2,2) = dy*dy
      Mij(3,3) = dz*dz

      r2(1:3,1:3) = 0.0d0      
      do mm = 1,3
         do nn = 1,3
            do oo = 1,3
               r2(mm,nn) = r2(mm,nn) + duij(mm,oo)*Mij(oo,nn)
            enddo
         enddo
      enddo

      Gij(1:3,1:3) = 0.0d0      
      do mm = 1,3
         do nn = 1,3
            do oo = 1,3
               Gij(mm,nn) = Gij(mm,nn) + r2(oo,mm)*r2(oo,nn)
            enddo
         enddo
      enddo

      ! magnitude might be better
      re_g = sqrt(Gij(1,1)+Gij(2,2)+Gij(3,3)) / nu


   end subroutine get_reg_old


!=========================================================
! velocity gradient
!=========================================================
   subroutine get_grad (nv, dat, d_lo, d_hi, ncomp, &
                         lo, hi, domlo, domhi, delta, &
                         xlo, time, dt, bc, &
                         level, grid_no, ii, jj, kk, Gij)

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no
      REAL_T, intent(out) :: Gij(3,3)

!  Local
      REAL_T :: ux, uy, uz, vx, vy, vz, wx, wy, wz, dx, dy, dz
      REAL_T :: uxcen, uycen, uzcen, uxlo, uxhi, uylo, uyhi, uzlo, uzhi
      REAL_T :: vxcen, vycen, vzcen, vxlo, vxhi, vylo, vyhi, vzlo, vzhi
      REAL_T :: wxcen, wycen, wzcen, wxlo, wxhi, wylo, wyhi, wzlo, wzhi
      REAL_T :: nu_m, re_g

      logical :: fixulo_x, fixvlo_x, fixwlo_x
      logical :: fixulo_y, fixvlo_y, fixwlo_y
      logical :: fixulo_z, fixvlo_z, fixwlo_z

      logical :: fixuhi_x, fixvhi_x, fixwhi_x
      logical :: fixuhi_y, fixvhi_y, fixwhi_y
      logical :: fixuhi_z, fixvhi_z, fixwhi_z

      integer :: ii, jj, kk
      integer :: i, j, k

!
!     ::::: some useful macro definitions
!
#     define U(i,j,k) dat(i,j,k,1)
#     define V(i,j,k) dat(i,j,k,2)
#     define W(i,j,k) dat(i,j,k,3)

#     define ULOX bc(1,1,1)
#     define UHIX bc(1,2,1)
#     define ULOY bc(2,1,1)
#     define UHIY bc(2,2,1)
#     define ULOZ bc(3,1,1)
#     define UHIZ bc(3,2,1)

#     define VLOX bc(1,1,2)
#     define VHIX bc(1,2,2)
#     define VLOY bc(2,1,2)
#     define VHIY bc(2,2,2)
#     define VLOZ bc(3,1,2)
#     define VHIZ bc(3,2,2)

#     define WLOX bc(1,1,3)
#     define WHIX bc(1,2,3)
#     define WLOY bc(2,1,3)
#     define WHIY bc(2,2,3)
#     define WLOZ bc(3,1,3)
#     define WHIZ bc(3,2,3)


!
!     ::::: statement functions that implement stencil
!
      uxcen(i,j,k) = half*(U(i+1,j,k)-U(i-1,j,k))/dx
      uxlo(i,j,k)  = (U(i+1,j,k)+three*U(i,j,k)-four*U(i-1,j,k))/(three*dx)
      uxhi(i,j,k)  =-(U(i-1,j,k)+three*U(i,j,k)-four*U(i+1,j,k))/(three*dx)

      uycen(i,j,k) = half*(U(i,j+1,k)-U(i,j-1,k))/dy
      uylo(i,j,k)  = (U(i,j+1,k)+three*U(i,j,k)-four*U(i,j-1,k))/(three*dy)
      uyhi(i,j,k)  =-(U(i,j-1,k)+three*U(i,j,k)-four*U(i,j+1,k))/(three*dy)

      vxcen(i,j,k) = half*(V(i+1,j,k)-V(i-1,j,k))/dx
      vxlo(i,j,k)  = (V(i+1,j,k)+three*V(i,j,k)-four*V(i-1,j,k))/(three*dx)
      vxhi(i,j,k)  =-(V(i-1,j,k)+three*V(i,j,k)-four*V(i+1,j,k))/(three*dx)

      vycen(i,j,k) = half*(V(i,j+1,k)-V(i,j-1,k))/dy
      vylo(i,j,k)  = (V(i,j+1,k)+three*V(i,j,k)-four*V(i,j-1,k))/(three*dy)
      vyhi(i,j,k)  =-(V(i,j-1,k)+three*V(i,j,k)-four*V(i,j+1,k))/(three*dy)

#if ( AMREX_SPACEDIM == 3 )
      uzcen(i,j,k) = half*(U(i,j,k+1)-U(i,j,k-1))/dz
      uzlo(i,j,k)  = (U(i,j,k+1)+three*U(i,j,k)-four*U(i,j,k-1))/(three*dz)
      uzhi(i,j,k)  =-(U(i,j,k-1)+three*U(i,j,k)-four*U(i,j,k+1))/(three*dz)

      vzcen(i,j,k) = half*(V(i,j,k+1)-V(i,j,k-1))/dz
      vzlo(i,j,k)  = (V(i,j,k+1)+three*V(i,j,k)-four*V(i,j,k-1))/(three*dz)
      vzhi(i,j,k)  =-(V(i,j,k-1)+three*V(i,j,k)-four*V(i,j,k+1))/(three*dz)

      wxcen(i,j,k) = half*(W(i+1,j,k)-W(i-1,j,k))/dx
      wxlo(i,j,k)  = (W(i+1,j,k)+three*W(i,j,k)-four*W(i-1,j,k))/(three*dx)
      wxhi(i,j,k)  =-(W(i-1,j,k)+three*W(i,j,k)-four*W(i+1,j,k))/(three*dx)

      wycen(i,j,k) = half*(W(i,j+1,k)-W(i,j-1,k))/dy
      wylo(i,j,k)  = (W(i,j+1,k)+three*W(i,j,k)-four*W(i,j-1,k))/(three*dy)
      wyhi(i,j,k)  =-(W(i,j-1,k)+three*W(i,j,k)-four*W(i,j+1,k))/(three*dy)

      wzcen(i,j,k) = half*(W(i,j,k+1)-W(i,j,k-1))/dz
      wzlo(i,j,k)  = (W(i,j,k+1)+three*W(i,j,k)-four*W(i,j,k-1))/(three*dz)
      wzhi(i,j,k)  =-(W(i,j,k-1)+three*W(i,j,k)-four*W(i,j,k+1))/(three*dz)
#endif



!      print*, "HERE: get_grad"


      ! copy over index
      i = ii
      j = jj
      k = kk

      ! local names
      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

      ! Init all logical tests on BC to false
      fixulo_x = .FALSE. ; fixvlo_x = .FALSE. ; fixwlo_x = .FALSE.
      fixulo_y = .FALSE. ; fixvlo_y = .FALSE. ; fixwlo_y = .FALSE.
      fixulo_z = .FALSE. ; fixvlo_z = .FALSE. ; fixwlo_z = .FALSE.

      fixuhi_x = .FALSE. ; fixvhi_x = .FALSE. ; fixwhi_x = .FALSE.
      fixuhi_y = .FALSE. ; fixvhi_y = .FALSE. ; fixwhi_y = .FALSE.
      fixuhi_z = .FALSE. ; fixvhi_z = .FALSE. ; fixwhi_z = .FALSE.



      ! Init all grad comp. In 2d uz, vz, wx, wy will alway be zero
      ux = 0.0d0
      uy = 0.0d0
      uz = 0.0d0
      vx = 0.0d0
      vy = 0.0d0
      vz = 0.0d0
      wx = 0.0d0
      wy = 0.0d0
      wz = 0.0d0


      ! indices passed

               ux = uxcen(i,j,k)
               uy = uycen(i,j,k)
               vx = vxcen(i,j,k)
               vy = vycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wx = wxcen(i,j,k)
               wy = wycen(i,j,k)
               wz = wzcen(i,j,k)
#endif


      fixulo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (ULOX .eq. EXT_DIR .or. ULOX .eq. HOEXTRAP) )
      fixuhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (UHIX .eq. EXT_DIR .or. UHIX .eq. HOEXTRAP) )
      fixulo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (ULOY .eq. EXT_DIR .or. ULOY .eq. HOEXTRAP) )
      fixuhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (UHIY .eq. EXT_DIR .or. UHIY .eq. HOEXTRAP) )

      fixvlo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (VLOX .eq. EXT_DIR .or. VLOX .eq. HOEXTRAP) )
      fixvhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (VHIX .eq. EXT_DIR .or. VHIX .eq. HOEXTRAP) )
      fixvlo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (VLOY .eq. EXT_DIR .or. VLOY .eq. HOEXTRAP) )
      fixvhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (VHIY .eq. EXT_DIR .or. VHIY .eq. HOEXTRAP) )

#if ( AMREX_SPACEDIM == 3 )
      fixulo_z = ( (lo(3) .eq. domlo(3)) .and. &
                  (ULOZ .eq. EXT_DIR .or. ULOZ .eq. HOEXTRAP) )
      fixuhi_z = ( (hi(3) .eq. domhi(3)) .and. &
                  (UHIZ .eq. EXT_DIR .or. UHIZ .eq. HOEXTRAP) )

      fixvlo_z = ( (lo(3) .eq. domlo(3)) .and. &
                  (VLOZ .eq. EXT_DIR .or. VLOZ .eq. HOEXTRAP) )
      fixvhi_z = ( (hi(3) .eq. domhi(3)) .and. &
                  (VHIZ .eq. EXT_DIR .or. VHIZ .eq. HOEXTRAP) )

      fixwlo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (WLOX .eq. EXT_DIR .or. WLOX .eq. HOEXTRAP) )
      fixwhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (WHIX .eq. EXT_DIR .or. WHIX .eq. HOEXTRAP) )
      fixwlo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (WLOY .eq. EXT_DIR .or. WLOY .eq. HOEXTRAP) )
      fixwhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (WHIY .eq. EXT_DIR .or. WHIY .eq. HOEXTRAP) )
      fixwlo_z = ( (lo(3) .eq. domlo(3)) .and. &
                  (WLOZ .eq. EXT_DIR .or. WLOZ .eq. HOEXTRAP) )
      fixwhi_z = ( (hi(3) .eq. domhi(3)) .and. &
                  (WHIZ .eq. EXT_DIR .or. WHIZ .eq. HOEXTRAP) )
#endif


      ! 1) First do all the faces...
      if (fixulo_x .or. fixvlo_x .or. fixwlo_x) then
         i = lo(1)
               ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
               uy = uycen(i,j,k)
               vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
               vy = vycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
               wy = wycen(i,j,k)
               wz = wzcen(i,j,k)
#endif

      end if

      if (fixuhi_x .or. fixvhi_x .or. fixwhi_x) then
         i = hi(1)

               ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
               uy = uycen(i,j,k)
               vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
               vy = vycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
               wy = wycen(i,j,k)
               wz = uzcen(i,j,k)
#endif

      end if

      if (fixulo_y .or. fixvlo_y .or. fixwlo_y) then
         j = lo(2)

               ux = uxcen(i,j,k)
               uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
               vx = vxcen(i,j,k)
               vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
#if ( AMREX_SPACEDIM == 3 )
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wx = wxcen(i,j,k)
               wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
               wz = wzcen(i,j,k)
#endif

      end if

      if (fixuhi_y .or. fixvhi_y .or. fixwhi_y) then
         j = hi(2)

               ux = uxcen(i,j,k)
               uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
               vx = vxcen(i,j,k)
               vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
#if ( AMREX_SPACEDIM == 3 )
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wx = wxcen(i,j,k)
               wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
               wz = wzcen(i,j,k)
#endif

      end if

#if ( AMREX_SPACEDIM == 3 )
      if (fixulo_z .or. fixvlo_z .or. fixwlo_z) then
         k = lo(3)

               ux = uxcen(i,j,k)
               uy = uycen(i,j,k)
               uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)

               vx = vxcen(i,j,k)
               vy = vycen(i,j,k)
               vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)

               wx = wxcen(i,j,k)
               wy = wycen(i,j,k)
               wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)

      end if

      if (fixuhi_z .or. fixvhi_z .or. fixwhi_z) then
         k = hi(3)

               ux = uxcen(i,j,k)
               uy = uycen(i,j,k)
               uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)

               vx = vxcen(i,j,k)
               vy = vycen(i,j,k)
               vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)

               wx = wxcen(i,j,k)
               wy = wycen(i,j,k)
               wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)

      end if
#endif

 
      ! 2) Next do all the edges...
      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixvlo_y .or. fixwlo_y)) then ! lo x-y edge
         i = lo(1)
         j = lo(2)

            ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
#if ( AMREX_SPACEDIM == 3 )
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            wz = wzcen(i,j,k)
#endif

      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixvlo_y .or. fixwlo_y)) then ! hi/lo x-y edge
         i = hi(1)
         j = lo(2)

            ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
#if ( AMREX_SPACEDIM == 3 )
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            wz = wzcen(i,j,k)
#endif

      end if

      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixvhi_y .or. fixwhi_y)) then
         i = lo(1)
         j = hi(2)

            ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
#if ( AMREX_SPACEDIM == 3 )
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            wz = wzcen(i,j,k)
#endif

      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixvhi_y .or. fixwhi_y)) then
         i = hi(1)
         j = hi(2)

            ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
#if ( AMREX_SPACEDIM == 3 )
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            wz = wzcen(i,j,k)
#endif

      end if

      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         i = lo(1)
         k = lo(3)

            ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
            uy = uycen(i,j,k)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            vy = vycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            wy = wycen(i,j,k)
            wx = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
#endif

      end if

      if ((fixvhi_x .or. fixuhi_x .or. fixwhi_x) .and. (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         i = hi(1)
         k = lo(3)

            ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
            uy = uycen(i,j,k)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            vy = vycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            wy = wycen(i,j,k)
            wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
#endif

      end if

      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         i = lo(1)
         k = hi(3)

            ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
            uy = uycen(i,j,k)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            vy = uycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            wy = wycen(i,j,k)
            wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
#endif

      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         i = hi(1)
         k = hi(3)

            ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
            uy = uycen(i,j,k)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            vy = vycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            wy = wycen(i,j,k)
            uz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
#endif

      end if

      if ((fixulo_x .or. fixvlo_y .or. fixwlo_y) .and. (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         j = lo(2)
         k = lo(3)

            ux = uxcen(i,j,k)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            vx = vxcen(i,j,k)
            vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
#if ( AMREX_SPACEDIM == 3 )
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            wx = wxcen(i,j,k)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
#endif

      end if

      if ((fixuhi_y .or. fixvhi_y .or. fixwhi_y) .and. (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         j = hi(2)
         k = lo(3)

            ux = uxcen(i,j,k)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            vx = vxcen(i,j,k)
            vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
#if ( AMREX_SPACEDIM == 3 )
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            wx = wxcen(i,j,k)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
#endif

      end if

      if ((fixulo_y .or. fixvlo_y .or. fixwlo_y) .and. (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         j = lo(2)
         k = hi(3)

            ux = uxcen(i,j,k)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            vx = vxcen(i,j,k)
            vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
#if ( AMREX_SPACEDIM == 3 )
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            wx = wxcen(i,j,k)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
#endif

      end if

      if ((fixuhi_y .or. fixvhi_y .or. fixwhi_y) .and. (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         j = hi(2)
         k = hi(3)

            ux = uxcen(i,j,k)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            vx = vxcen(i,j,k)
            vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
#if ( AMREX_SPACEDIM == 3 )
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            wx = wxcen(i,j,k)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
#endif

      end if



      ! 3) Finally do all the corners...
      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. &
          (fixulo_y .or. fixvlo_y .or. fixwlo_y) .and. &
          (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         i = lo(1)
         j = lo(2)
         k = lo(3)
         ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
#if ( AMREX_SPACEDIM == 3 )
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
#endif

      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. &
          (fixulo_y .or. fixvlo_y .or. fixwlo_y) .and. &
          (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         i = hi(1)
         j = lo(2)
         k = lo(3)
         ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
#if ( AMREX_SPACEDIM == 3 )
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
#endif

      end if

      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. &
          (fixuhi_y .or. fixvhi_y .or. fixwhi_y) .and. &
          (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         i = lo(1)
         j = hi(2)
         k = lo(3)
         ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
#if ( AMREX_SPACEDIM == 3 )
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
#endif

      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. &
          (fixuhi_y .or. fixvhi_y .or. fixwhi_y) .and. &
          (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         i = hi(1)
         j = hi(2)
         k = lo(3)
         ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
#if ( AMREX_SPACEDIM == 3 )
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
#endif

      end if

      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. &
          (fixulo_y .or. fixvlo_y .or. fixwlo_y) .and. &
          (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         i = lo(1)
         j = lo(2)
         k = hi(3)
         ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
#if ( AMREX_SPACEDIM == 3 )
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
#endif

      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. &
          (fixulo_y .or. fixvlo_y .or. fixwlo_y) .and. &
          (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         i = hi(1)
         j = lo(2)
         k = hi(3)
         ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
#if ( AMREX_SPACEDIM == 3 )
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
#endif

      end if

      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. &
          (fixuhi_y .or. fixvhi_y .or. fixwhi_y) .and. &
          (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         i = lo(1)
         j = hi(2)
         k = hi(3)
         ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
#if ( AMREX_SPACEDIM == 3 )
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
#endif

      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. &
          (fixuhi_y .or. fixvhi_y .or. fixwhi_y) .and. &
          (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         i = hi(1)
         j = hi(2)
         k = hi(3)
         ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
#if ( AMREX_SPACEDIM == 3 )
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
#endif

      end if

      Gij(1,1) = ux
      Gij(1,2) = uy
      Gij(1,3) = uz
      Gij(2,1) = vx
      Gij(2,2) = vy
      Gij(2,3) = vz
      Gij(3,1) = wx
      Gij(3,2) = wy
      Gij(3,3) = wz




#     undef U
#     undef V
#     undef W
#     undef ULOY
#     undef UHIY
#     undef ULOZ
#     undef UHIZ
#     undef VLOX
#     undef VHIX
#     undef VLOZ
#     undef VHIZ
#     undef WLOX
#     undef WHIX
#     undef WLOY
#     undef WHIY

!   print*, "   DONE WITH get_grad"

   end subroutine get_grad


!=========================================================
!  Compute the magnitude of the velocity divergence
!=========================================================

   subroutine dermgdivu (e,   e_lo, e_hi, nv, &
                         dat, d_lo, d_hi, ncomp, &
                         lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                         level, grid_no) &
                         bind(C, name="dermgdivu")

#if ( AMREX_SPACEDIM == 2 )                      
      use prob_2D_module, only : FORT_XVELFILL, FORT_YVELFILL
#elif ( AMREX_SPACEDIM == 3 )
      use prob_3D_module, only : FORT_XVELFILL, FORT_YVELFILL, FORT_ZVELFILL
#endif

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(inout), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      REAL_T  :: ux, vy, wz, dx, dy, dz
      REAL_T  :: uxcen, uxlo, uxhi
      REAL_T  :: vycen, vylo, vyhi
      REAL_T  :: wzcen, wzlo, wzhi
      integer :: bc_dumm(2,2,1)
      integer :: i, j, k

!
!     ::::: some useful macro definitions
!
#     define U(i,j,k) dat(i,j,k,1)
#     define V(i,j,k) dat(i,j,k,2)
#     define W(i,j,k) dat(i,j,k,3)

#     define ULOX bc(1,1,1)
#     define UHIX bc(1,2,1)
#     define VLOY bc(2,1,2)
#     define VHIY bc(2,2,2)
#     define WLOZ bc(3,1,3)
#     define WHIZ bc(3,2,3)

!
!     ::::: statement functions that implement stencil
!
      uxcen(i,j,k) = half*(U(i+1,j,k)-U(i-1,j,k))/dx
      uxlo(i,j,k) = (eight*U(i,j,k)-six*U(i+1,j,k)+U(i+2,j,k))/(three*dx)
      uxhi(i,j,k) = (eight*U(i,j,k)-six*U(i-1,j,k)+U(i-2,j,k))/(three*dx)

#if ( AMREX_SPACEDIM >= 2 )
      vycen(i,j,k) = half*(V(i,j+1,k)-V(i,j-1,k))/dy
      vylo(i,j,k) = (eight*V(i,j,k)-six*V(i,j+1,k)+V(i,j+2,k))/(three*dy)
      vyhi(i,j,k) = (eight*V(i,j,k)-six*V(i,j-1,k)+V(i,j-2,k))/(three*dy)

#if ( AMREX_SPACEDIM == 3 )
      wzcen(i,j,k) = half*(W(i,j,k+1)-W(i,j,k-1))/dz
      wzlo(i,j,k) = (eight*W(i,j,k)-six*W(i,j,k+1)+W(i,j,k+2))/(three*dz)
      wzhi(i,j,k) = (eight*W(i,j,k)-six*W(i,j,k-1)+W(i,j,k-2))/(three*dz)
#endif
#endif

#if ( AMREX_SPACEDIM == 2 )
      bc_dumm(1:2,1:2,1) = bc(1:2,1:2,1)
      call FORT_XVELFILL (dat(:,:,d_lo(3),1), d_lo(1), d_lo(2), d_hi(1), d_hi(2), & 
                          domlo, domhi, delta, xlo, time, bc_dumm(1,1,1))
      bc_dumm(1:2,1:2,1) = bc(1:2,1:2,2)
      call FORT_YVELFILL (dat(:,:,d_lo(3),2), d_lo(1), d_lo(2), d_hi(1), d_hi(2), &
                          domlo, domhi, delta, xlo, time, bc_dumm(1,1,1))
#elif ( AMREX_SPACEDIM == 3 )
      call FORT_XVELFILL (dat(:,:,:,1), d_lo(1), d_lo(2), d_lo(3), d_hi(1), d_hi(2), d_hi(3), & 
                          domlo, domhi, delta, xlo, time, bc(1,1,1))
      call FORT_YVELFILL (dat(:,:,:,2), d_lo(1), d_lo(2), d_lo(3), d_hi(1), d_hi(2), d_hi(3), &
                          domlo, domhi, delta, xlo, time, bc(1,1,2))
      call FORT_ZVELFILL (dat(:,:,:,3), d_lo(1), d_lo(2), d_lo(3), d_hi(1), d_hi(2), d_hi(3), &
                          domlo, domhi, delta, xlo, time, bc(1,1,3))
#endif

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

!
!     :: at physical bndries where an edge value is prescribed,
!     :: set the value in the outside cell so that a central
!     :: difference formula is equivalent to the higher order
!     :: one sided formula
!
      if (lo(1) == domlo(1)) then
         i = lo(1)
         if (ULOX==EXT_DIR) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  U(i-1,j,k) = two*U(i-1,j,k) - U(i,j,k)
               end do
            end do
         else if (ULOX==HOEXTRAP) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  U(i-1,j,k) = uxlo(i,j,k)
               end do
            end do
         end if
      end if
      if (hi(1) == domhi(1)) then
         i = hi(1)
         if (UHIX==EXT_DIR) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  U(i+1,j,k) = two*U(i+1,j,k) - U(i,j,k)
               end do
            end do
         else if (UHIX==HOEXTRAP) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  U(i+1,j,k) = uxhi(i,j,k)
               end do
            end do
         end if
      end if

#if ( AMREX_SPACEDIM >= 2 )
      if (lo(2) == domlo(2)) then
         j = lo(2)
         if (VLOY==EXT_DIR) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  V(i,j-1,k) = two*V(i,j-1,k) - V(i,j,k)
               end do
            end do
         else if (VLOY==HOEXTRAP) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  V(i,j-1,k) = vylo(i,j,k)
               end do
            end do
         end if
      end if
      if (hi(2) == domhi(2)) then
         j = hi(2)
         if (VHIY==EXT_DIR) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  V(i,j+1,k) = two*V(i,j+1,k) - V(i,j,k)
               end do
            end do
         else if (VHIY==HOEXTRAP) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  V(i,j+1,k) = vyhi(i,j,k)
               end do
            end do
         end if
      end if

#if ( AMREX_SPACEDIM == 3 )
      if (lo(3) == domlo(3)) then
         k = lo(3)
         if (WLOZ==EXT_DIR) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  W(i,j,k-1) = two*W(i,j,k-1) - W(i,j,k)
               end do
            end do
         else if (WLOZ==HOEXTRAP) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  W(i,j,k-1) = wzlo(i,j,k)
               end do
            end do
         end if
      end if
      if (hi(3) == domhi(3)) then
         k = hi(3)
         if (WHIZ==EXT_DIR) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  W(i,j,k+1) = two*W(i,j,k+1) - W(i,j,k)
               end do
            end do
         else if (WHIZ==HOEXTRAP) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  W(i,j,k+1) = wzhi(i,j,k)
               end do
            end do
         end if
      end if
#endif
#endif

      ux = 0.0d0
      vy = 0.0d0
      wz = 0.0d0
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ux = uxcen(i,j,k)
#if ( AMREX_SPACEDIM >= 2 )
               vy = vycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
               wz = wzcen(i,j,k)
#endif
#endif
               e(i,j,k,1) = ux + vy + wz
            end do
         end do
      end do

!
! we overwrote the ghost cells above, so set them back below
!
#if ( AMREX_SPACEDIM == 2 )
      bc_dumm(1:2,1:2,1) = bc(1:2,1:2,1)
      call FORT_XVELFILL (dat(:,:,d_lo(3),1), d_lo(1), d_lo(2), d_hi(1), d_hi(2), & 
                          domlo, domhi, delta, xlo, time, bc_dumm(1,1,1))
      bc_dumm(1:2,1:2,1) = bc(1:2,1:2,2)
      call FORT_YVELFILL (dat(:,:,d_lo(3),2), d_lo(1), d_lo(2), d_hi(1), d_hi(2), &
                          domlo, domhi, delta, xlo, time, bc_dumm(1,1,1))
#elif ( AMREX_SPACEDIM == 3 )
      call FORT_XVELFILL (dat(:,:,:,1), d_lo(1), d_lo(2), d_lo(3), d_hi(1), d_hi(2), d_hi(3), & 
                          domlo, domhi, delta, xlo, time, bc(1,1,1))
      call FORT_YVELFILL (dat(:,:,:,2), d_lo(1), d_lo(2), d_lo(3), d_hi(1), d_hi(2), d_hi(3), &
                          domlo, domhi, delta, xlo, time, bc(1,1,2))
      call FORT_ZVELFILL (dat(:,:,:,3), d_lo(1), d_lo(2), d_lo(3), d_hi(1), d_hi(2), d_hi(3), &
                          domlo, domhi, delta, xlo, time, bc(1,1,3))
#endif

#     undef U
#     undef V      
#     undef W
#     undef ULOX
#     undef UHIX
#     undef VLOY
#     undef VHIY
#     undef WLOZ
#     undef WHIZ

   end subroutine dermgdivu

!=========================================================
!  Compute C / rho
!=========================================================

   subroutine derdvrho (e,   e_lo, e_hi, nv, &
                        dat, d_lo, d_hi, ncomp, &
                        lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                        level, grid_no) &
                        bind(C, name="derdvrho")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      integer :: i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) = dat(i,j,k,2)/dat(i,j,k,1)
            end do
         end do
      end do

   end subroutine derdvrho

!=========================================================
!  Compute rho * C 
!=========================================================

   subroutine dermprho (e,   e_lo, e_hi, nv, &
                        dat, d_lo, d_hi, ncomp, &
                        lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                        level, grid_no) &
                        bind(C, name="dermprho")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      integer :: i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) = dat(i,j,k,2)*dat(i,j,k,1)
            end do
         end do
      end do

   end subroutine dermprho

!=========================================================
!  Compute cell-centered pressure as average of the 
!  surrounding nodal values 
!=========================================================

   subroutine deravgpres (e,   e_lo, e_hi, nv, &
                          dat, d_lo, d_hi, ncomp, &
                          lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                          level, grid_no) &
                          bind(C, name="deravgpres")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      REAL_T  :: factor
      integer :: i, j, k

      factor = 0.5d0
#if (AMREX_SPACEDIM >= 2 )
      factor = 0.25d0
#if (AMREX_SPACEDIM == 3 )
      factor = 0.125d0
#endif
#endif

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) =  factor * (  dat(i+1,j,k,1)     + dat(i,j,k,1)  &
#if (AMREX_SPACEDIM >= 2 )
                                       + dat(i+1,j+1,k,1)   + dat(i,j+1,k,1) &
#if (AMREX_SPACEDIM == 3 )
                                       + dat(i+1,j,k+1,1)   + dat(i,j,k+1,1)  &
                                       + dat(i+1,j+1,k+1,1) + dat(i,j+1,k+1,1) &
#endif
#endif
                                      )    
            end do
         end do
      end do

   end subroutine deravgpres

!=========================================================
!  Compute node centered pressure gradient in direction dir
!=========================================================

   subroutine gradp_dir (p, p_lo, p_hi, &
                         gp, g_lo, g_hi, &
                         lo, hi, dir, dx) &
                         bind(C, name="gradp_dir")

      implicit none

! In/Out
      integer :: lo(3),  hi(3)
      integer :: p_lo(3), p_hi(3)
      integer :: g_lo(3), g_hi(3)
      integer :: dir
      REAL_T  :: dx
      REAL_T, dimension(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3)) :: p
      REAL_T, dimension(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3)) :: gp

! Local
      integer :: i, j, k
      REAL_T  :: d

      d = one/dx
#if (AMREX_SPACEDIM >= 2 )
      d = half/dx
#if (AMREX_SPACEDIM == 3 )
      d = fourth/dx
#endif
#endif

!
!     ::::: compute gradient on interior
!
      if (dir == 0) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
#if (AMREX_SPACEDIM == 1 )
                  gp(i,j,k) = d*(p(i+1,j,k)-p(i,j,k))
#elif (AMREX_SPACEDIM == 2 )
                  gp(i,j,k) = d*(p(i+1,j,k)-p(i,j,k)+p(i+1,j+1,k)-p(i,j+1,k))
#elif (AMREX_SPACEDIM == 3 )
                  gp(i,j,k) = d*( p(i+1,j,k  )-p(i,j,k  )+p(i+1,j+1,k  )-p(i,j+1,k  ) + &
                                  p(i+1,j,k+1)-p(i,j,k+1)+p(i+1,j+1,k+1)-p(i,j+1,k+1))
#endif 

               end do
            end do
         end do
      else if (dir == 1) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
#if (AMREX_SPACEDIM == 2 )
                  gp(i,j,k) = d*(p(i,j+1,k)-p(i,j,k)+p(i+1,j+1,k)-p(i+1,j,k))
#elif (AMREX_SPACEDIM == 3 )
                  gp(i,j,k) = d*( p(i,j+1,k  )-p(i,j,k  )+p(i+1,j+1,k  )-p(i+1,j,k  ) + &
                                  p(i,j+1,k+1)-p(i,j,k+1)+p(i+1,j+1,k+1)-p(i+1,j,k+1))
#endif
               end do
            end do
         end do
      else if (dir == 2) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  gp(i,j,k) = d*( p(i,  j,k+1)-p(i,  j,k)+p(i,  j+1,k+1)-p(i,  j+1,k) + &
                                  p(i+1,j,k+1)-p(i+1,j,k)+p(i+1,j+1,k+1)-p(i+1,j+1,k))
               end do
            end do
         end do
      else
         call amrex_abort("gradp_dir: invalid dir = ")
      end if

   end subroutine gradp_dir

!=========================================================
!  Compute node centered pressure gradient in X-dir
!=========================================================

   subroutine dergrdpx (e,   e_lo, e_hi, nv, &
                        dat, d_lo, d_hi, ncomp, &
                        lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                        level, grid_no) &
                        bind(C, name="dergrdpx")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      REAL_T  :: factor
      integer :: i, j, k

      call gradp_dir ( dat, d_lo, d_hi, &
                       e,   e_lo, e_hi, &
                       lo, hi, 0, delta(1))

   end subroutine dergrdpx

!=========================================================
!  Compute node centered pressure gradient in Y-dir
!=========================================================

   subroutine dergrdpy (e,   e_lo, e_hi, nv, &
                        dat, d_lo, d_hi, ncomp, &
                        lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                        level, grid_no) &
                        bind(C, name="dergrdpy")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      REAL_T  :: factor
      integer :: i, j, k

#if (AMREX_SPACEDIM < 2 )
      call amrex_abort("dergrdpy called but AMREX_SPACEDIM<2 !")
#endif

      call gradp_dir ( dat, d_lo, d_hi, &
                       e,   e_lo, e_hi, &
                       lo, hi, 1, delta(2))

   end subroutine dergrdpy

!=========================================================
!  Compute node centered pressure gradient in Z-dir
!=========================================================

   subroutine dergrdpz (e,   e_lo, e_hi, nv, &
                        dat, d_lo, d_hi, ncomp, &
                        lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                        level, grid_no) &
                        bind(C, name="dergrdpz")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      REAL_T  :: factor
      integer :: i, j, k

#if (AMREX_SPACEDIM < 3 )
      call amrex_abort("dergrdpz called but AMREX_SPACEDIM<3 !")
#endif

      call gradp_dir ( dat, d_lo, d_hi, &
                       e,   e_lo, e_hi, &
                       lo, hi, 2, delta(3))

   end subroutine dergrdpz

end module derive_nd_module
