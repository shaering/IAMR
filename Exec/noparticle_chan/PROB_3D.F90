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

module prob_3D_module

   implicit none

   private

   public :: amrex_probinit, FORT_INITDATA, &
   !   &         FORT_AVERAGE_EDGE_STATES, &
   !   &         FORT_MAKEFORCE, FORT_DSDTFILL, &
   &         FORT_DSDTFILL, &   
   &         FORT_ADVERROR, FORT_ADV2ERROR, FORT_TEMPERROR, FORT_MVERROR, &
   &         FORT_DENFILL, FORT_ADVFILL, FORT_TEMPFILL, FORT_XVELFILL, &
   &         FORT_YVELFILL, FORT_ZVELFILL, FORT_PRESFILL, FORT_DIVUFILL, &
   &         FORT_RGERROR

   !
   ! Define some parameters
   ! These could be loaded via a probin file instead of being hardcoded
   ! like we do here
   !
   REAL_T, parameter :: U0       = 1.0
   REAL_T, parameter :: density  = 1.0
   REAL_T, parameter :: temp_bot  = 373.15
   REAL_T, parameter :: temp_top  = 273.15

   !
   ! Set problo and probhi as module variables so that they can be
   ! used everywehere in this module
   !
   REAL_T :: m_problo(SDIM), m_probhi(SDIM)

contains


   !c ::: -----------------------------------------------------------
   !c ::: This routine is called at problem initialization time
   !c ::: and when restarting from a checkpoint file.
   !c ::: The purpose is (1) to specify the initial time value
   !c ::: (not all problems start at time=0.0) and (2) to read
   !c ::: problem specific data from a namelist or other inputcdm
   !c ::: files and possibly store them or derived information
   !c ::: in FORTRAN common blocks for later use.
   !c :::
   !c ::: INPUTS/OUTPUTS:
   !c :::
   !c ::: init      => TRUE if called at start of problem run
   !c :::              FALSE if called from restart
   !c ::: name      => name of "probin" file
   !c ::: namlen    => length of name
   !c ::: strttime <=  start problem with this time variable
   !c :::
   !c ::: -----------------------------------------------------------

!   subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)
!      implicit none
!      integer init,namlen
!      integer name(namlen)
!      REAL_T  problo(SDIM), probhi(SDIM)

      !
      ! No need to read any probin file for this test
      ! We just Set module variables
      !
!      m_problo = problo
!      m_probhi = probhi

!    end subroutine amrex_probinit


      !-----------------------------------------------------------------
      subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)
      !-----------------------------------------------------------------
      implicit none
      integer init,namlen
      integer name(namlen)
      integer untin, i
      REAL_T  problo(SDIM), probhi(SDIM)

#include <probdata.H>

!c
! Dimensions of the Inflow file.
!c
      integer nCompFile
      parameter (nCompFile = 2)

      namelist /fortin/ rgerr
!c
!     Build "probin" filename -- the name of file containing fortin namelist.
!c
      integer maxlen, isioproc
      parameter (maxlen=256)

      character probin*(maxlen)

      call bl_pd_is_ioproc(isioproc)

      if (namlen .gt. maxlen) call bl_error('probin file name too long')

      do i = 1, namlen
         probin(i:i) = char(name(i))
      end do

      untin = 9
      if (namlen .eq. 0) then
         open(untin,file='probin',form='formatted',status='old')
      else
         open(untin,file=probin(1:namlen),form='formatted',status='old')
      end if

      read(untin,fortin)
      if (isioproc .eq. 1) write(6,fortin)
      close(unit=untin)

    end subroutine amrex_probinit



   !c ::: -----------------------------------------------------------
   !c ::: This routine is called at problem setup time and is used
   !c ::: to initialize data on each grid.  The velocity field you
   !c ::: provide does not have to be divergence free and the pressure
   !c ::: field need not be set.  A subsequent projection iteration
   !c ::: will define a divergence free velocity field along with a
   !c ::: consistant pressure.
   !c :::
   !c ::: NOTE:  all arrays have one cell of ghost zones surrounding
   !c :::        the grid interior.  Values in these cells need not
   !c :::        be set here.
   !c :::
   !c ::: INPUTS/OUTPUTS:
   !c :::
   !c ::: level     => amr level of grid
   !c ::: time      => time at which to init data
   !c ::: lo,hi     => index limits of grid interior (cell centered)
   !c ::: nscal     => number of scalar quantities.  You should know
   !c :::		   this already!
   !c ::: vel      <=  Velocity array
   !c ::: scal     <=  Scalar array
   !c ::: press    <=  Pressure array
   !c ::: dx     => cell size
   !c ::: xlo,xhi   => physical locations of lower left and upper
   !c :::              right hand corner of grid.  (does not include
   !c :::		   ghost region).
   !c ::: -----------------------------------------------------------
   subroutine FORT_INITDATA(level,time,lo,hi,nscal, &
   vel,scal,DIMS(state),press,DIMS(press), &
   dx,xlo,xhi) &
   bind(C, name="FORT_INITDATA")

      implicit none
      integer    level, nscal
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(state)
      integer    DIMDEC(press)
      REAL_T     time, dx(SDIM)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     vel(DIMV(state),SDIM)
      REAL_T    scal(DIMV(state),nscal)
      REAL_T   press(DIMV(press))
      REAL_T   ATG, LTG

      integer i, j, k
      REAL_T  x, y, z, yn, hx, hy, hz

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half)
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
            do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)) + half)

!               yn = y / m_probhi(2)
!               vel(i,j,k,1) = 6.0d0 * U0 * yn * (1.0 - yn)
               yn = y - 1.0d0 ! channel goes from 0 to 2
               vel(i,j,k,1) = 1.1d0 * U0 * (1.0d0 - yn**4)
               vel(i,j,k,2) = zero
               vel(i,j,k,3) = zero

               ! add perturbations
               ATG = U0/10.0d0
               LTG = 1.0d0*3.14159265359
               vel(i,j,k,1) = vel(i,j,k,1) + ATG * 1.0d0 * cos(LTG*x)*sin(LTG*yn)*sin(LTG*z)
               vel(i,j,k,2) = vel(i,j,k,2) - ATG * 3.0d0 * sin(LTG*x)*cos(LTG*yn)*sin(LTG*z)
               vel(i,j,k,3) = vel(i,j,k,3) + ATG * 2.0d0 * sin(LTG*x)*sin(LTG*yn)*cos(LTG*z)

               ! Density
               scal(i,j,k,1) = density

               ! All other scalars -- even if this case does not
               ! trace anything, we still initialize the tracers
               ! arrays to  zero
               scal(i,j,k,2:nscal) = zero

               ! Temp, not sure how this gets a 3...
               scal(i,j,k,3) = temp_bot + 0.5d0*y*(temp_top-temp_bot)

            end do
         end do
      end do

   end subroutine FORT_INITDATA

!------------------------------------------------------------

   
   !c     ::: -----------------------------------------------------------
   !c     ::: This routine will tag high error cells based on the
   !c     ::: magnitude or gradient of the density
   !c     :::
   !c     ::: INPUTS/OUTPUTS:
   !c     :::
   !c     ::: tag      <=  integer tag array
   !c     ::: DIMS(tag) => index extent of tag array
   !c     ::: set       => integer value to tag cell for refinement
   !c     ::: clear     => integer value to untag cell
   !c     ::: rho       => density array
   !c     ::: DIMS(rho) => index extent of rho array
   !c     ::: nvar      => number of components in rho array (should be 1)
   !c     ::: lo,hi     => index extent of grid
   !c     ::: domlo,hi  => index extent of problem domain
   !c     ::: dx        => cell spacing
   !c     ::: xlo       => physical location of lower left hand
   !c     :::	           corner of tag array
   !c     ::: problo    => phys loc of lower left corner of prob domain
   !c     ::: time      => problem evolution time
   !c     ::: -----------------------------------------------------------
   subroutine FORT_DENERROR (tag,DIMS(tag),set,clear, &
   rho,DIMS(rho),lo,hi,nvar, &
   domlo,domhi,dx,xlo,  &
   problo,time,level)&
   bind(C, name="FORT_DENERROR")
      implicit none
      integer   DIMDEC(tag)
      integer   DIMDEC(rho)
      integer   lo(SDIM), hi(SDIM)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    rho(DIMV(rho),nvar)

      integer   i
      integer :: isioproc

      ! For this case, refine the first half of the domain
!      do i = lo(1), hi(1)
!         if (i<domhi(1)/2) then
!            tag(i,:,:) = set
!         end if
!      end do

      call bl_pd_is_ioproc(isioproc)
      if (isioproc.eq.1) then
         stop "FORT_DENERROR not implemented for this case"
      else
         stop
      endif

   end subroutine FORT_DENERROR

   !c ::: -----------------------------------------------------------
   !c ::: This routine will tag high error cells based on the
   !c ::: magnitude of the tracer
   !c :::
   !c ::: INPUTS/OUTPUTS:
   !c :::
   !c ::: tag      <=  integer tag array
   !c ::: DIMS(tag) => index extent of tag array
   !c ::: set       => integer value to tag cell for refinement
   !c ::: clear     => integer value to untag cell
   !c ::: adv       => scalar array
   !c ::: DIMS(adv) => index extent of adv array
   !c ::: nvar      => number of components in rho array (should be 1)
   !c ::: lo,hi     => index extent of grid
   !c ::: domlo,hi  => index extent of problem domain
   !c ::: dx        => cell spacing
   !c ::: xlo       => physical location of lower left hand
   !c :::	           corner of tag array
   !c ::: problo    => phys loc of lower left corner of prob domain
   !c ::: time      => problem evolution time
   !c ::: -----------------------------------------------------------
   subroutine FORT_ADVERROR (tag,DIMS(tag),set,clear, &
   adv,DIMS(adv),lo,hi,nvar, &
   domlo,domhi,dx,xlo, &
   problo,time,level)&
   bind(C, name="FORT_ADVERROR")
      implicit none
      integer   DIMDEC(tag)
      integer   DIMDEC(adv)
      integer   lo(SDIM), hi(SDIM)
      integer   ng, nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    adv(DIMV(adv),nvar)

      integer :: isioproc

      call bl_pd_is_ioproc(isioproc)
      if (isioproc.eq.1) then
         stop "FORT_ADVERROR not implemented for this case"
      else
         stop
      endif

   end subroutine FORT_ADVERROR


   subroutine FORT_ADV2ERROR (tag,DIMS(tag),set,clear, &
   adv,DIMS(adv),lo,hi,nvar, &
   domlo,domhi,dx,xlo, &
   problo,time,level)&
   bind(C, name="FORT_ADV2ERROR")
      implicit none
      integer   DIMDEC(tag)
      integer   DIMDEC(adv)
      integer   lo(SDIM), hi(SDIM)
      integer   ng, nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    adv(DIMV(adv),nvar)

      integer :: isioproc

      call bl_pd_is_ioproc(isioproc)
      if (isioproc.eq.1) then
         stop "FORT_ADV2ERROR not implemented for this case"
      else
         stop
      endif

   end subroutine FORT_ADV2ERROR

   !c ::: -----------------------------------------------------------
   !c ::: This routine will tag high error cells based on the
   !c ::: magnitude or gradient of temperature
   !c :::
   !c ::: INPUTS/OUTPUTS:
   !c :::
   !c ::: tag      <=  integer tag array
   !c ::: DIMS(tag) => index extent of tag array
   !c ::: set       => integer value to tag cell for refinement
   !c ::: clear     => integer value to untag cell
   !c ::: temp      => density array
   !c ::: DIMS(temp)=> index extent of temp array
   !c ::: lo,hi     => index extent of grid
   !c ::: nvar      => number of components in rho array (should be 1)
   !c ::: domlo,hi  => index extent of problem domain
   !c ::: dx        => cell spacing
   !c ::: xlo       => physical location of lower left hand
   !c :::              corner of tag array
   !c ::: problo    => phys loc of lower left corner of prob domain
   !c ::: time      => problem evolution time
   !c ::: -----------------------------------------------------------
   subroutine FORT_TEMPERROR (tag,DIMS(tag),set,clear, &
   temperature,DIMS(temp),lo,hi,nvar, &
   domlo,domhi,dx,xlo, &
   problo,time,level)&
   bind(C, name="FORT_TEMPERROR")
      implicit none

      integer   DIMDEC(tag)
      integer   DIMDEC(temp)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      integer   lo(SDIM), hi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    temperature(DIMV(temp),nvar)

      integer :: isioproc

      call bl_pd_is_ioproc(isioproc)
      if (isioproc.eq.1) then
         stop "FORT_TEMPERROR not implemented for this case"
      else
         stop
      endif

   end subroutine FORT_TEMPERROR


!::: -----------------------------------------------------------
!::: This routine will tag high error cells based on the 
!::: magnitude of vorticity
!::: 
!::: INPUTS/OUTPUTS:
!::: 
!::: tag      <=  integer tag array
!::: DIMS(tag) => index extent of tag array
!::: set       => integer value to tag cell for refinement
!::: clear     => integer value to untag cell
!::: vort      => array of vorticity values
!::: DIMS(vor) => index extent of vort array
!::: nvar      => number of components in vort array (should be 1)
!::: lo,hi     => index extent of grid
!::: domlo,hi  => index extent of problem domain
!::: dx        => cell spacing
!::: xlo       => physical location of lower left hand
!:::	           corner of tag array
!::: problo    => phys lo!of lower left corner of prob domain
!::: time      => problem evolution time
!::: -----------------------------------------------------------
      subroutine FORT_MVERROR (tag,DIMS(tag),set,clear,&
                              vort,DIMS(vort),lo,hi,nvar,&
                              domlo,domhi,dx,xlo,&
     			       problo,time,level) bind(c,name="FORT_MVERROR")

      integer   i, j, k
      integer   DIMDEC(tag)
      integer   DIMDEC(vort)
      integer   nvar, set, clear, level
      integer   lo(SDIM), hi(SDIM)
      integer   domlo(SDIM), domhi(SDIM)
      integer   tag(DIMV(tag))
      REAL_T    vort(DIMV(vort),nvar)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time


#include <probdata.H>

!      do j = lo(2), hi(2)
!         do i = lo(1), hi(1)
!            tag(i,j) = merge(set,tag(i,j),abs(vort(i,j,1)).gt.vorterr*2.d0**level)
!         end do
!      end do

!      print*, " HERE!", tag(i,j,k), vorterr
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tag(i,j,k) = clear
               tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr*2.0d0**level)
!               if (abs(vort(i,j,k,1)) .GT. vorterr*2.0d0**level) print*, vort(i,j,k,1), vorterr
            end do
         end do
      end do

    end subroutine FORT_MVERROR


!c ::: -----------------------------------------------------------
      subroutine FORT_RGERROR (tag,DIMS(tag),set,clear,&
                               re_g,DIMS(re_g),lo,hi,nvar,&
                               domlo,domhi,dx,xlo,&
     			       problo,time,level) bind(c,name="FORT_RGERROR")

      integer   i, j, k
      integer   DIMDEC(tag)
      integer   DIMDEC(re_g)
      integer   nvar, set, clear, level
      integer   lo(SDIM), hi(SDIM)
      integer   domlo(SDIM), domhi(SDIM)
      integer   tag(DIMV(tag))
      REAL_T    re_g(DIMV(re_g),nvar)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time


#include <probdata.H>

!      print*, " RGERROR DUMP:"
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
!               tag(i,j,k) = clear

               tag(i,j,k) = merge(set,tag(i,j,k),re_g(i,j,k,1) .GT. rgerr )
!               tag(i,j,k) = merge(set,tag(i,j,k),re_g(i,j,k,1) .GT. rgerr / (10.0d0**level))

!               if (re_g(i,j,k,1) .GT. rgerr) print*, "TEST:", level, re_g(i,j,k,1), rgerr, tag(i,j,k), clear, set
!               if (re_g(i,j,k,1) .GT. rgerr) then 
!                  print*, "LEVEL++:", level, "Re_G:", re_g(i,j,k,1), "TAG:", tag(i,j,k)
!               else
!                  print*, "LEVEL:", level, "Re_G:", re_g(i,j,k,1), "TAG:", tag(i,j,k)
!               endif
!               if (level .eq. 2) print*, "***LEVEL 2 ACTIVE***", re_g(i,j,k,1)
            enddo
         enddo
      enddo

    end subroutine FORT_RGERROR


   !c ::: -----------------------------------------------------------
   !c ::: This routine is called during a filpatch operation when
   !c ::: the patch to be filled falls outside the interior
   !c ::: of the problem domain.  You are requested to supply the
   !c ::: data outside the problem interior in such a way that the
   !c ::: data is consistant with the types of the boundary conditions
   !c ::: you specified in the C++ code.
   !c :::
   !c ::: NOTE:  you can assume all interior cells have been filled
   !c :::        with valid data.
   !c :::
   !c ::: INPUTS/OUTPUTS:
   !c :::
   !c ::: rho      <=  density array
   !c ::: DIMS(rho) => index extent of rho array
   !c ::: domlo,hi  => index extent of problem domain
   !c ::: dx        => cell spacing
   !c ::: xlo       => physical location of lower left hand
   !c :::	           corner of rho array
   !c ::: time      => problem evolution time
   !c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
   !c ::: -----------------------------------------------------------

   subroutine FORT_DENFILL (rho,DIMS(rho),domlo,domhi,dx, &
   xlo,time,bc ) &
   bind(C, name="FORT_DENFILL")
      implicit none

      integer    DIMDEC(rho)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     rho(DIMV(rho))
      integer    bc(SDIM,2)

      integer    i, j, k, n

      ! filcc fills bc_types foextrap, hoextrap, reflect_odd and reflect_even
      call filcc(rho,DIMS(rho),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(rho).lt.domlo(1)) then
         do i = ARG_L1(rho), domlo(1)-1
            do k = ARG_L3(rho), ARG_H3(rho)
               do j = ARG_L2(rho), ARG_H2(rho)
                  rho(i,j,k) = density
               end do
            end do
         end do
      end if

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(rho).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(rho)
            do k = ARG_L3(rho), ARG_H3(rho)
               do j = ARG_L2(rho), ARG_H2(rho)
                  rho(i,j,k) = density
               end do
            end do
         end do
      end if

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(rho).lt.domlo(2)) then
         do j = ARG_L2(rho), domlo(2)-1
            do k = ARG_L3(rho), ARG_H3(rho)
               do i = ARG_L1(rho), ARG_H1(rho)
                  rho(i,j,k) = density
               end do
            end do
         end do
      end if

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(rho).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(rho)
            do k = ARG_L3(rho), ARG_H3(rho)
               do i = ARG_L1(rho), ARG_H1(rho)
                  rho(i,j,k) = density
               end do
            end do
         end do
      end if

      if (bc(3,1).eq.EXT_DIR.and.ARG_L3(rho).lt.domlo(3)) then
         do k = ARG_L3(rho), domlo(3)-1
            do j = ARG_L2(rho), ARG_H2(rho)
               do i = ARG_L1(rho), ARG_H1(rho)
                  rho(i,j,k) = density
               end do
            end do
         end do
      endif

      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(rho).gt.domhi(3)) then
         do k = domhi(3)+1, ARG_H3(rho)
            do j = ARG_L2(rho), ARG_H2(rho)
               do i = ARG_L1(rho), ARG_H1(rho)
                  rho(i,j,k) = density
               end do
            end do
         end do
      end if

   end subroutine FORT_DENFILL

   !c ::: -----------------------------------------------------------
   !c ::: This routine is called during a filpatch operation when
   !c ::: the patch to be filled falls outside the interior
   !c ::: of the problem domain.  You are requested to supply the
   !c ::: data outside the problem interior in such a way that the
   !c ::: data is consistant with the types of the boundary conditions
   !c ::: you specified in the C++ code.
   !c :::
   !c ::: NOTE:  you can assume all interior cells have been filled
   !c :::        with valid data.
   !c :::
   !c ::: INPUTS/OUTPUTS:
   !c :::
   !c ::: adv      <=  advected quantity array
   !c ::: DIMS(adv) => index extent of adv array
   !c ::: domlo,hi  => index extent of problem domain
   !c ::: dx        => cell spacing
   !c ::: xlo       => physical location of lower left hand
   !c :::	           corner of adv array
   !c ::: time      => problem evolution time
   !c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
   !c ::: -----------------------------------------------------------

   subroutine FORT_ADVFILL (adv,DIMS(adv),domlo,domhi,dx,&
   xlo,time,bc )&
   bind(C, name="FORT_ADVFILL")
      implicit none

      integer    DIMDEC(adv)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     adv(DIMV(adv))
      integer    bc(SDIM,2)

      ! Since we are not advecting any scalar besides density
      ! we do not need this subroutine to do anything

   end subroutine FORT_ADVFILL

   subroutine FORT_ADV2FILL (adv,DIMS(adv),domlo,domhi,dx, &
   xlo,time,bc )&
   bind(C, name="FORT_ADV2FILL")

      implicit none

      integer    DIMDEC(adv)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     adv(DIMV(adv))
      integer    bc(SDIM,2)
      integer    lo(SDIM), hi(SDIM)

      ! Since we are not advecting any scalar besides density
      ! we do not need this subroutine to do anything

   end subroutine FORT_ADV2FILL

   !c ::: -----------------------------------------------------------
   !c ::: This routine is called during a filpatch operation when
   !c ::: the patch to be filled falls outside the interior
   !c ::: of the problem domain.  You are requested to supply the
   !c ::: data outside the problem interior in such a way that the
   !c ::: data is consistant with the types of the boundary conditions
   !c ::: you specified in the C++ code.
   !c :::
   !c ::: NOTE:  you can assume all interior cells have been filled
   !c :::        with valid data.
   !c :::
   !c ::: INPUTS/OUTPUTS:
   !c :::
   !c ::: temp      <= temperature array
   !c ::: DIMS(temp)=> index extent of temp array
   !c ::: domlo,hi  => index extent of problem domain
   !c ::: dx        => cell spacing
   !c ::: xlo       => physical location of lower left hand
   !c :::	           corner of temp array
   !c ::: time      => problem evolution time
   !c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
   !c ::: -----------------------------------------------------------

   subroutine FORT_TEMPFILL (temp,DIMS(temp),domlo,domhi,dx,&
   xlo,time,bc )&
   bind(C, name="FORT_TEMPFILL")

      implicit none

      integer    DIMDEC(temp)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     temp(DIMV(temp))
      integer    bc(SDIM,2)
      integer    i, j, k, n


      ! filcc fills bc_types foextrap, hoextrap, reflect_odd and reflect_even
      call filcc(temp,DIMS(temp),domlo,domhi,dx,xlo,bc)

      ! really only care about bottom and top wall for this case, 
      ! others are not called

      if (bc(1,1) .eq. EXT_DIR .and. ARG_L1(temp) .lt. domlo(1)) then
         do i = ARG_L1(temp), domlo(1)-1
            do k = ARG_L3(temp), ARG_H3(temp)
               do j = ARG_L2(temp), ARG_H2(temp)
                  temp(i,j,k) = temp_bot
               end do
            end do
         end do
      end if

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(temp).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(temp)
            do k = ARG_L3(temp), ARG_H3(temp)
               do j = ARG_L2(temp), ARG_H2(temp)
                  temp(i,j,k) = temp_top
               end do
            end do
         end do
      end if

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(temp).lt.domlo(2)) then
         do j = ARG_L2(temp), domlo(2)-1
            do k = ARG_L3(temp), ARG_H3(temp)
               do i = ARG_L1(temp), ARG_H1(temp)
                  temp(i,j,k) = temp_bot
               end do
            end do
         end do
      end if

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(temp).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(temp)
            do k = ARG_L3(temp), ARG_H3(temp)
               do i = ARG_L1(temp), ARG_H1(temp)
                  temp(i,j,k) = temp_top
               end do
            end do
         end do
      end if

      if (bc(3,1).eq.EXT_DIR.and.ARG_L3(temp).lt.domlo(3)) then
         do k = ARG_L3(temp), domlo(3)-1
            do j = ARG_L2(temp), ARG_H2(temp)
               do i = ARG_L1(temp), ARG_H1(temp)
                  temp(i,j,k) = temp_bot
               end do
            end do
         end do
      endif

      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(temp).gt.domhi(3)) then
         do k = domhi(3)+1, ARG_H3(temp)
            do j = ARG_L2(temp), ARG_H2(temp)
               do i = ARG_L1(temp), ARG_H1(temp)
                  temp(i,j,k) = temp_top
               end do
            end do
         end do
      end if


   end subroutine FORT_TEMPFILL

   !c ::: -----------------------------------------------------------
   !c ::: This routine is called during a filpatch operation when
   !c ::: the patch to be filled falls outside the interior
   !c ::: of the problem domain.  You are requested to supply the
   !c ::: data outside the problem interior in such a way that the
   !c ::: data is consistant with the types of the boundary conditions
   !c ::: you specified in the C++ code.
   !c :::
   !c ::: NOTE:  you can assume all interior cells have been filled
   !c :::        with valid data.
   !c :::
   !c ::: INPUTS/OUTPUTS:
   !c :::
   !c ::: u        <=  x velocity array
   !c ::: DIMS(u)   => index extent of u array
   !c ::: domlo,hi  => index extent of problem domain
   !c ::: dx        => cell spacing
   !c ::: xlo       => physical location of lower left hand
   !c :::	           corner of rho array
   !c ::: time      => problem evolution time
   !c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
   !c ::: -----------------------------------------------------------

   subroutine FORT_XVELFILL (u,DIMS(u),domlo,domhi,dx,xlo,time,bc) &
   &          bind(C, name="FORT_XVELFILL")

      implicit none

      integer    DIMDEC(u)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     u(DIMV(u))
      integer    bc(SDIM,2)
      integer    lo(SDIM),hi(SDIM)
      integer    i, j, k
      REAL_T     y

      lo(1) = ARG_L1(u)
      lo(2) = ARG_L2(u)
      lo(3) = ARG_L3(u)
      hi(1) = ARG_H1(u)
      hi(2) = ARG_H2(u)
      hi(3) = ARG_H3(u)

      call filcc(u,DIMS(u),domlo,domhi,dx,xlo,bc)

      !
      ! At the inlet we enforce a parabolic profile
      ! At the outlet we do not need any condition on U since
      ! we have an outflow condition
      !
!      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
!         do i = lo(1), domlo(1)-1
!            do k = lo(3), hi(3)
!               do j = lo(2), hi(2)
!                  y = (xlo(2) + dx(2)*(float(j-lo(2)) + half))/m_probhi(2)
!                  u(i,j,k) = 6.0d0 * U0 * y * (1.0d0 - y)
!               end do
!            end do
!         end do
!      end if

      ! At No-slip walls (= EXT_DIR for tangential component of velocity),
      ! set tangential velocity to zero
      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  u(i,j,k) = zero
               end do
            end do
         end do
      end if

      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  u(i,j,k) = zero
               end do
            end do
         end do
      end if

      ! No need to fill along third direction because periodic conditions
      ! apply


   end subroutine FORT_XVELFILL

   !c ::: -----------------------------------------------------------
   !c ::: This routine is called during a filpatch operation when
   !c ::: the patch to be filled falls outside the interior
   !c ::: of the problem domain.  You are requested to supply the
   !c ::: data outside the problem interior in such a way that the
   !c ::: data is consistant with the types of the boundary conditions
   !c ::: you specified in the C++ code.
   !c :::
   !c ::: NOTE:  you can assume all interior cells have been filled
   !c :::        with valid data.
   !c :::
   !c ::: INPUTS/OUTPUTS:
   !c :::
   !c ::: v        <=  y velocity array
   !c ::: DIMS(v)   => index extent of v array
   !c ::: domlo,hi  => index extent of problem domain
   !c ::: dx        => cell spacing
   !c ::: xlo       => physical location of lower left hand
   !c :::	           corner of rho array
   !c ::: time      => problem evolution time
   !c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
   !c ::: -----------------------------------------------------------

   subroutine FORT_YVELFILL (v,DIMS(v),domlo,domhi,dx,xlo,time,bc) &
   &          bind(C, name="FORT_YVELFILL")

      implicit none
      integer    DIMDEC(v)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     v(DIMV(v))
      integer    bc(SDIM,2)
      integer    lo(SDIM),hi(SDIM)
      integer    i, j, k

      lo(1) = ARG_L1(v)
      lo(2) = ARG_L2(v)
      lo(3) = ARG_L3(v)
      hi(1) = ARG_H1(v)
      hi(2) = ARG_H2(v)
      hi(3) = ARG_H3(v)

      call filcc(v,DIMS(v),domlo,domhi,dx,xlo,bc)

!       if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
!         do i = lo(1), domlo(1)-1
!            do k = lo(3), hi(3)
!               do j = lo(2), hi(2)
!                  v(i,j,k) = zero
!               end do
!            end do
!         end do
!      end if

!      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
!         do i = domhi(1)+1, hi(1)
!            do k = lo(3), hi(3)
!               do j = lo(2), hi(2)
!                  v(i,j,k) = zero
!               end do
!            end do
!         end do
!      end if

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  v(i,j,k) = zero
               end do
            end do
         end do
      end if

      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  v(i,j,k) = zero
               end do
            end do
         end do
      end if

      ! No need to fill along third direction because periodic conditions
      ! apply

   end subroutine FORT_YVELFILL

   !c ::: -----------------------------------------------------------
   !c ::: This routine is called during a fillpatch operation when
   !c ::: the patch to be filled falls outside the interior
   !c ::: of the problem domain.  You are requested to supply the
   !c ::: data outside the problem interior in such a way that the
   !c ::: data is consistant with the types of the boundary conditions
   !c ::: you specified in the C++ code.
   !c :::
   !c ::: NOTE:  you can assume all interior cells have been filled
   !c :::        with valid data.
   !c :::
   !c ::: INPUTS/OUTPUTS:
   !c :::
   !c ::: w        <=  z velocity array
   !c ::: DIMS(w)   => index extent of v array
   !c ::: domlo,hi  => index extent of problem domain
   !c ::: dx        => cell spacing
   !c ::: xlo       => physical location of lower left hand
   !c :::	           corner of rho array
   !c ::: time      => problem evolution time
   !c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
   !c ::: -----------------------------------------------------------

   subroutine FORT_ZVELFILL (w,DIMS(w),domlo,domhi,dx,xlo,time,bc) &
   &          bind(C, name="FORT_ZVELFILL")

      implicit none
      integer    DIMDEC(w)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     w(DIMV(w))
      integer    bc(SDIM,2)
      integer    lo(SDIM),hi(SDIM)
      integer    i, j, k

      lo(1) = ARG_L1(w)
      lo(2) = ARG_L2(w)
      lo(3) = ARG_L3(w)
      hi(1) = ARG_H1(w)
      hi(2) = ARG_H2(w)
      hi(3) = ARG_H3(w)

      call filcc(w,DIMS(w),domlo,domhi,dx,xlo,bc)

 !     if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
 !        do i = lo(1), domlo(1)-1
 !           do k = lo(3), hi(3)
 !              do j = lo(2), hi(2)
 !                 w(i,j,k) = zero
 !              end do
 !           end do
 !        end do
 !     end if

 !     if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
 !        do i = domhi(1)+1, hi(1)
 !           do k = lo(3), hi(3)
 !              do j = lo(2), hi(2)
 !                 w(i,j,k) = zero
 !              end do
 !           end do
 !        end do
 !     end if

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  w(i,j,k) = zero
               end do
            end do
         end do
      end if

      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  w(i,j,k) = zero
               end do
            end do
         end do
      end if

      ! No need to fill along third direction because periodic conditions
      ! apply

   end subroutine FORT_ZVELFILL


   !c ::: -----------------------------------------------------------
   !c ::: This routine is called during a filpatch operation when
   !c ::: the patch to be filled falls outside the interior
   !c ::: of the problem domain.  You are requested to supply the
   !c ::: data outside the problem interior in such a way that the
   !c ::: data is consistant with the types of the boundary conditions
   !c ::: you specified in the C++ code.
   !c :::
   !c ::: NOTE:  you can assume all interior cells have been filled
   !c :::        with valid data.
   !c :::
   !c ::: INPUTS/OUTPUTS:
   !c :::
   !c ::: divu     <=  divergence of velocity array
   !c ::: DIMS(divu)=> index extent of divu array
   !c ::: domlo,hi  => index extent of problem domain
   !c ::: dx        => cell spacing
   !c ::: xlo       => physical location of lower left hand
   !c :::	           corner of rho array
   !c ::: time      => problem evolution time
   !c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
   !c ::: -----------------------------------------------------------
   subroutine FORT_DIVUFILL (divu,DIMS(divu),domlo,domhi,dx, &
   xlo,time,bc)&
   bind(C, name="FORT_DIVUFILL")
      implicit none

      integer    DIMDEC(divu)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     divu(DIMV(divu))
      integer    bc(SDIM,2)

      ! Since we are not advecting any scalar besides density
      ! we do not need this subroutine to do anything

   end subroutine FORT_DIVUFILL

   !c ::: -----------------------------------------------------------
   !c ::: This routine is called during a filpatch operation when
   !c ::: the patch to be filled falls outside the interior
   !c ::: of the problem domain.  You are requested to supply the
   !c ::: data outside the problem interior in such a way that the
   !c ::: data is consistant with the types of the boundary conditions
   !c ::: you specified in the C++ code.
   !c :::
   !c ::: NOTE:  you can assume all interior cells have been filled
   !c :::        with valid data.
   !c :::
   !c ::: INPUTS/OUTPUTS:
   !c :::
   !c ::: dsdt     <=  dsdt array
   !c ::: DIMS(dsdt)=> index extent of dsdt array
   !c ::: domlo,hi  => index extent of problem domain
   !c ::: dx        => cell spacing
   !c ::: xlo       => physical location of lower left hand
   !c :::	           corner of rho array
   !c ::: time      => problem evolution time
   !c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
   !c ::: -----------------------------------------------------------
   subroutine FORT_DSDTFILL (dsdt,DIMS(dsdt),domlo,domhi,dx, &
   xlo,time,bc) &
   bind(C, name="FORT_DSDTFILL")
      implicit none

      integer    DIMDEC(dsdt)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     dsdt(DIMV(dsdt))
      integer    bc(SDIM,2)

      ! Since we are not advecting any scalar besides density
      ! we do not need this subroutine to do anything

   end subroutine FORT_DSDTFILL

   !c ::: -----------------------------------------------------------
   !c ::: This routine is called during a filpatch operation when
   !c ::: the patch to be filled falls outside the interior
   !c ::: of the problem domain.  You are requested to supply the
   !c ::: data outside the problem interior in such a way that the
   !c ::: data is consistant with the types of the boundary conditions
   !c ::: you specified in the C++ code.
   !c :::
   !c ::: NOTE:  you can assume all interior cells have been filled
   !c :::        with valid data.
   !c :::
   !c ::: INPUTS/OUTPUTS:
   !c :::
   !c ::: p        <=  pressure array
   !c ::: lo,hi     => index extent of p array
   !c ::: domlo,hi  => index extent of problem domain
   !c ::: dx        => cell spacing
   !c ::: xlo       => physical location of lower left hand
   !c :::	           corner of rho array
   !c ::: time      => problem evolution time
   !c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
   !c ::: -----------------------------------------------------------
   subroutine FORT_PRESFILL (p,DIMS(p),domlo,domhi,dx, &
   xlo,time,bc) &
   bind(C, name="FORT_PRESFILL")
      implicit none

      integer    DIMDEC(p)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     p(DIMV(p))
      integer    bc(SDIM,2)

      ! we do not need this subroutine to do anything

   end subroutine FORT_PRESFILL


end module prob_3D_module
