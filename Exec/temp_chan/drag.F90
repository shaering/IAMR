

module drag_pic_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_real

  implicit none

contains


!===============================================================================================
! This routine computes the charge density due to the particles using cloud-in-cell 
! deposition. The Fab rho is assumed to be node-centered.
!
! Arguments:
!     particles : a pointer to the particle array-of-structs 
!     ns        : the stride length of particle struct (the size of the struct in number of reals)
!     np        : the number of particles
!     weights   : the particle weights
!     charge    : the charge of this particle species
!     drag      : a Fab that will contain the drag force on exit
!     lo        : a pointer to the lo corner of this valid box for rho, in index space
!     hi        : a pointer to the hi corner of this valid box for rho, in index space
!     plo       : the real position of the left-hand corner of the problem domain
!     dx        : the mesh spacing
!     ng        : the number of ghost cells in rho
!
!===============================================================================================
  subroutine drag_cic(ng,ns,np,lo,hi,nu_m,dx,plo,vel,rho,drg,particles) &
             bind(c,name='drag_cic')

    integer, intent(in)  :: ng 
    integer, intent(in)  :: ns
    integer, intent(in)  :: np
    integer, intent(in)  :: lo(3)
    integer, intent(in)  :: hi(3)

    real(amrex_real)  :: particles(ns,np)
    real(amrex_real)  :: nu_m
    real(amrex_real)  :: plo(3)
    real(amrex_real)  :: dx(3)
    real(amrex_real)  :: drg(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng,0:2)
    real(amrex_real)  :: vel(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng,0:2)
    real(amrex_real)  :: rho(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng)

    integer i, j, k, n, oo
    real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
    real(amrex_real) lx, ly, lz
    real(amrex_real) qp, inv_vol, Re_p, tau_p, dia_p, rho_p
    real(amrex_real) inv_dx(3)
    real(amrex_real) vdiff, CT, vp(3)
    real(amrex_real) v_ap(3), Fp(3), rho_f


    ! face-centered should go from 0:N+1
    ! cell-centered should go from 0:N

    inv_dx = 1.0d0/dx
    inv_vol = inv_dx(1) * inv_dx(2) * inv_dx(3)

    ! This will be incorrect if particle diameters are (~O) of dx
    do n = 1, np

       !qp = weights(n) * charge * inv_vol

       ! particle properties (10 is temp)
       vp(1) = particles(7,n)
       vp(2) = particles(8,n)
       vp(3) = particles(9,n)
       !T_p = particles(10,n)
       dia_p = particles(11,n)
       rho_p = particles(12,n)
       Fp(1) = particles(13,n)
       Fp(2) = particles(14,n)
       Fp(3) = particles(15,n)


       if(isnan(vp(1))) print*, " up NAN at:", n
       if(isnan(vp(2))) print*, " vp NAN at:", n
       if(isnan(vp(3))) print*, " wp NAN at:", n

       if(isnan(dia_p)) print*, " diap NAN at:", n
       if(isnan(rho_p)) print*, " rhop NAN at:", n

       if(isnan(Fp(1))) print*, " Fpx NAN at:", n
       if(isnan(Fp(2))) print*, " Fpy NAN at:", n
       if(isnan(Fp(3))) print*, " Fpz NAN at:", n

       

!       print*, "... in drag_cic (lo_x, hi_x, ng): ", lo(1),hi(1),ng

!       print*, " "
!       print*, " (DRAG) PARTICLE No. ", n, "(stride: )", ns, 2*3+3+4
!       print*, "==================="
!       print*, "    particle dump: 1 (xp)", particles(1,n)
!       print*, "    particle dump: 2 (yp)", particles(2,n)
!       print*, "    particle dump: 3 (zp)", particles(3,n)
!       print*, "    particle dump: 4 (*)", particles(4,n)
!       print*, "    particle dump: 5 (*)", particles(5,n)
!       print*, "    particle dump: 6 (*)", particles(6,n)
!       print*, "    particle dump: 7 (up)", particles(7,n)
!       print*, "    particle dump: 8 (vp)", particles(8,n)
!       print*, "    particle dump: 9 (wp)", particles(9,n)
!       print*, "    particle dump: 10 (Tp)", particles(10,n)
!       print*, "    particle dump: 11 (dp)", particles(11,n)
!       print*, "    particle dump: 12 (rhop)", particles(12,n)
!       print*, "    particle dump: 13 (Fxp)", particles(13,n)
!       print*, "    particle dump: 14 (Fyp)", particles(14,n)
!       print*, "    particle dump: 15 (Fzp)", particles(15,n)
!       print*, "    particle dump: 16 (FTp)", particles(16,n)
!!       print*, "    particle dump: 17 (?)", particles(17,n)


       ! scaled location of particle: (pos - lower corner of tile)*1/dx = x_p/dx (relative to tile)
       lx = (particles(1, n) - plo(1))*inv_dx(1)
       ly = (particles(2, n) - plo(2))*inv_dx(2)
       lz = (particles(3, n) - plo(3))*inv_dx(3)

       if(isnan(lx)) print*, " lx NAN at:", n
       if(isnan(ly)) print*, " ly NAN at:", n
       if(isnan(lz)) print*, " lz NAN at:", n
       

       ! get cell index which contains particle center
       i = floor(lx)
       j = floor(ly)
       k = floor(lz)

       ! distance(normalized by dx) from particle center to cell center (or face, not sure)
       wx_hi = lx - dble(i)
       wy_hi = ly - dble(j)
       wz_hi = lz - dble(k)

       ! takes on a weight with the 1/dx normalization, weight of force on low index
       wx_lo = 1.0d0 - wx_hi
       wy_lo = 1.0d0 - wy_hi
       wz_lo = 1.0d0 - wz_hi

       ! fluid velocity at particle
!       print*, "RANGE 1: ", lo(1)-ng,hi(1)+ng
!       print*, "RANGE 2: ", lo(2)-ng,hi(2)+ng
!       print*, "RANGE 3: ", lo(3)-ng,hi(3)+ng
!       print*, "POINT: (", i, j, k, ") drg:", drg(i,j,k,0), drg(i,j,k,1), drg(i,j,k,2)
!       print*, "POINT: (", i, j, k, ") rho_p / rho_f:", rho_p, rho(i,j,k)
!       print*, "POINT: (", i, j, k, ") Fp:", Fp(1), Fp(2),Fp(3)
!?       print*, "POINT: (", i, j, k, ") vel:", vel(i,j,k,0), vel(i,j,k,1), vel(i,j,k,2)
!       v_ap(1) = vel(i,j,k,0)
!       v_ap(2) = vel(i,j,k,1)
!       v_ap(3) = vel(i,j,k,2)

       ! for simplified drag calc (Soldati 2009), should really interp vel to be at particle...
!       vdiff = (v_ap(1)-vp(1))*(v_ap(1)-vp(1)) + &
!             & (v_ap(2)-vp(2))*(v_ap(2)-vp(2)) + &
!             & (v_ap(3)-vp(3))*(v_ap(3)-vp(3))
!       vdiff = sqrt(vdiff)




       ! sum of drags, signs flips because this is action of particle on fluid
       do oo = 0,2

!          qp = rho_f * CT/tau_p * (vp(oo) - vel(i,j,k,oo)) * inv_vol ! not sure about units here, extra rho/dV...
!          qp = rho_f * CT/tau_p * (vp(oo) - v_ap(oo)) !* inv_vol
          qp = -1.0*Fp(oo+1)

          ! drag force should be stored at face center  (0 is leftmost face)
          drg(i,   j,   k  ,oo) = drg(i,   j,   k  ,oo) + wx_lo*wy_lo*wz_lo*qp
          drg(i,   j,   k+1,oo) = drg(i,   j,   k+1,oo) + wx_lo*wy_lo*wz_hi*qp
          drg(i,   j+1, k  ,oo) = drg(i,   j+1, k  ,oo) + wx_lo*wy_hi*wz_lo*qp
          drg(i,   j+1, k+1,oo) = drg(i,   j+1, k+1,oo) + wx_lo*wy_hi*wz_hi*qp
          drg(i+1, j,   k  ,oo) = drg(i+1, j,   k  ,oo) + wx_hi*wy_lo*wz_lo*qp
          drg(i+1, j,   k+1,oo) = drg(i+1, j,   k+1,oo) + wx_hi*wy_lo*wz_hi*qp
          drg(i+1, j+1, k  ,oo) = drg(i+1, j+1, k  ,oo) + wx_hi*wy_hi*wz_lo*qp
          drg(i+1, j+1, k+1,oo) = drg(i+1, j+1, k+1,oo) + wx_hi*wy_hi*wz_hi*qp

          if(isnan(qp)) print*, " qp NAN at:", n, oo

          
!          print*, ">>> DRAG_CIC:",oo
!          print*, drg(i,   j,   k  ,oo)
!          print*, drg(i,   j,   k+1,oo)
!          print*, drg(i,   j+1, k  ,oo)
!          print*, drg(i,   j+1, k+1,oo)
!          print*, drg(i+1, j,   k  ,oo)
!          print*, drg(i+1, j,   k+1,oo)
!          print*, drg(i+1, j+1, k  ,oo)
!          print*, drg(i+1, j+1, k+1,oo)

       enddo

    end do

!    print*, "*** Dumping periodic x-face ***"
!    print*, "============================== "
!    do j=lo(2)-ng,hi(2)+ng
!       do k=lo(3)-ng,hi(3)+ng
!          print*, lo(1)-ng, hi(1), k, j, drg(lo(1)-ng, j, k,0), drg(hi(1), j , k,0)
!          print*, lo(1), hi(1)+ng, k, j, drg(lo(1), j, k,0), drg(hi(1)+ng, j , k,0)
!       enddo
!    enddo


  end subroutine drag_cic


!==================================================================================
  subroutine temp_cic(ng,ns,np,lo,hi,nu_m,dx,plo,vel,temp,hf,particles) &
             bind(c,name='temp_cic')

    integer, intent(in)  :: ng 
    integer, intent(in)  :: ns
    integer, intent(in)  :: np
    integer, intent(in)  :: lo(3)
    integer, intent(in)  :: hi(3)

    real(amrex_real)  :: particles(ns,np)
    real(amrex_real)  :: nu_m
    real(amrex_real)  :: plo(3)
    real(amrex_real)  :: dx(3)
    real(amrex_real)  :: hf(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng)
    real(amrex_real)  :: vel(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng,0:2)
    real(amrex_real)  :: temp(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng)

    integer i, j, k, n, oo
    real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
    real(amrex_real) lx, ly, lz
    real(amrex_real) qp, inv_vol, Re_p, tau_p, dia_p, temp_p, rho_p
    real(amrex_real) inv_dx(3)
    real(amrex_real) vdiff, CT, vp(3)
    real(amrex_real) v_ap(3), Fp(3), FTp, temp_f, rho_f


    inv_dx = 1.0d0/dx
    inv_vol = inv_dx(1) * inv_dx(2) * inv_dx(3)

    ! This will be incorrect if particle diameters are (~O) of dx
    do n = 1, np

       ! particle properties (10 is temp)
       vp(1) = particles(7,n)
       vp(2) = particles(8,n)
       vp(3) = particles(9,n)
       temp_p = particles(10,n)
       dia_p = particles(11,n)
       rho_p = particles(12,n)
       Fp(1) = particles(13,n)
       Fp(2) = particles(14,n)
       Fp(3) = particles(15,n)
       FTp   = particles(16,n)

!       print*, " "
!       print*, " (TEMP) PARTICLE No. ", n, "(stride: )", ns, 2*3+3+4
!       print*, "==================="
!       print*, "    particle dump: 1 (xp)", particles(1,n)
!       print*, "    particle dump: 2 (yp)", particles(2,n)
!       print*, "    particle dump: 3 (zp)", particles(3,n)
!       print*, "    particle dump: 4 (*)", particles(4,n)
!       print*, "    particle dump: 5 (*)", particles(5,n)
!       print*, "    particle dump: 6 (*)", particles(6,n)
!       print*, "    particle dump: 7 (up)", particles(7,n)
!       print*, "    particle dump: 8 (vp)", particles(8,n)
!       print*, "    particle dump: 9 (wp)", particles(9,n)
!       print*, "    particle dump: 10 (Tp)", particles(10,n)
!       print*, "    particle dump: 11 (dp)", particles(11,n)
!       print*, "    particle dump: 12 (rhop)", particles(12,n)
!       print*, "    particle dump: 13 (Fxp)", particles(13,n)
!       print*, "    particle dump: 14 (Fyp)", particles(14,n)
!       print*, "    particle dump: 15 (Fzp)", particles(15,n)
!       print*, "    particle dump: 16 (FTp)", particles(16,n)
!!       print*, "    particle dump: 17 (?)", particles(17,n)


       ! scaled location of particle: (pos - lower corner of tile)*1/dx = x_p/dx (relative to tile)
       lx = (particles(1, n) - plo(1))*inv_dx(1)
       ly = (particles(2, n) - plo(2))*inv_dx(2)
       lz = (particles(3, n) - plo(3))*inv_dx(3)

       ! get cell index which contains particle center
       i = floor(lx)
       j = floor(ly)
       k = floor(lz)

       ! distance(normalized by dx) from particle center to cell center (or face, not sure)
       wx_hi = lx - dble(i)
       wy_hi = ly - dble(j)
       wz_hi = lz - dble(k)

       ! takes on a weight with the 1/dx normalization, weight of force on low index
       wx_lo = 1.0d0 - wx_hi
       wy_lo = 1.0d0 - wy_hi
       wz_lo = 1.0d0 - wz_hi

       ! sum of drags, signs flips because this is action of particle on fluid
       qp = -1.0*FTp  !* inv_vol
       hf(i,   j,   k  ) = hf(i,   j,   k  ) + wx_lo*wy_lo*wz_lo*qp
       hf(i,   j,   k+1) = hf(i,   j,   k+1) + wx_lo*wy_lo*wz_hi*qp
       hf(i,   j+1, k  ) = hf(i,   j+1, k  ) + wx_lo*wy_hi*wz_lo*qp
       hf(i,   j+1, k+1) = hf(i,   j+1, k+1) + wx_lo*wy_hi*wz_hi*qp
       hf(i+1, j,   k  ) = hf(i+1, j,   k  ) + wx_hi*wy_lo*wz_lo*qp
       hf(i+1, j,   k+1) = hf(i+1, j,   k+1) + wx_hi*wy_lo*wz_hi*qp
       hf(i+1, j+1, k  ) = hf(i+1, j+1, k  ) + wx_hi*wy_hi*wz_lo*qp
       hf(i+1, j+1, k+1) = hf(i+1, j+1, k+1) + wx_hi*wy_hi*wz_hi*qp

    end do


  end subroutine temp_cic


!==================================================================================
  subroutine sum_fine_to_crse_nodal (lo, hi, lrat, crse, clo, chi, fine, flo, fhi) &
       bind(c, name="sum_fine_to_crse_nodal")

    integer, intent(in)             ::   lo(3),  hi(3)
    integer, intent(in)             ::  clo(3), chi(3)
    integer, intent(in)             ::  flo(3), fhi(3)
    integer, intent(in)             ::  lrat(3)
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(amrex_real), intent(in)    :: fine(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    
    integer :: i, j, k, ii, jj, kk
    
    do k        = lo(3), hi(3)
       kk       = k * lrat(3)
       do j     = lo(2), hi(2)
          jj    = j * lrat(2)
          do i  = lo(1), hi(1)
             ii = i * lrat(1)
             crse(i,j,k)  =  fine(ii,jj,kk)                              + &
! These six fine nodes are shared by two coarse nodes...
                  0.5d0   * (fine(ii-1,jj,kk)     + fine(ii+1,jj,kk)     + & 
                             fine(ii,jj-1,kk)     + fine(ii,jj+1,kk)     + &
                             fine(ii,jj,kk-1)     + fine(ii,jj,kk+1))    + &
! ... these twelve are shared by four...
                  0.25d0  * (fine(ii,jj-1,kk-1)   + fine(ii,jj+1,kk-1)   + &
                             fine(ii,jj-1,kk+1)   + fine(ii,jj+1,kk+1)   + &
                             fine(ii-1,jj,kk-1)   + fine(ii+1,jj,kk-1)   + &
                             fine(ii-1,jj,kk+1)   + fine(ii+1,jj,kk+1)   + &
                             fine(ii-1,jj-1,kk)   + fine(ii+1,jj-1,kk)   + &
                             fine(ii-1,jj+1,kk)   + fine(ii+1,jj+1,kk))  + &
! ... and these eight are shared by eight
                  0.125d0 * (fine(ii-1,jj-1,kk-1) + fine(ii-1,jj-1,kk+1) + &
                             fine(ii-1,jj+1,kk-1) + fine(ii-1,jj+1,kk+1) + &
                             fine(ii+1,jj-1,kk-1) + fine(ii+1,jj-1,kk+1) + &
                             fine(ii+1,jj+1,kk-1) + fine(ii+1,jj+1,kk+1))
! ... note that we have 27 nodes in total...
             crse(i,j,k) = crse(i,j,k) / 8.d0
          end do
       end do
    end do

  end subroutine sum_fine_to_crse_nodal


!===============================================================================================
  subroutine drag_cic_full(ng,ns,np,lo,hi,nu_m,dx,plo,vel,rho,drg,particles) &
             bind(c,name='drag_cic_full')

    integer, intent(in)  :: ng 
    integer, intent(in)  :: ns
    integer, intent(in)  :: np
    integer, intent(in)  :: lo(3)
    integer, intent(in)  :: hi(3)

    real(amrex_real)  :: particles(ns,np)
    real(amrex_real)  :: nu_m
    real(amrex_real)  :: plo(3)
    real(amrex_real)  :: dx(3)
    real(amrex_real)  :: drg(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng,3)
    real(amrex_real)  :: vel(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng,3) ! face centered (0 is left-most face, hi is right-most face)
    real(amrex_real)  :: rho(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng)   ! cell centered 


    ! face-centered should go from 0:N+1
    ! cell-centered should go from 0:N

    integer i, j, k, n, oo
    real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
    real(amrex_real) lx, ly, lz
    real(amrex_real) qp, inv_vol, Re_p, tau_p, dia_p, rho_p
    real(amrex_real) inv_dx(3)
    real(amrex_real) vdiff, CT, vp(3)
    real(amrex_real) v_ap(3), rho_f


!    print*, "... in drag_cic: ng,ns,lo,hi,nu_m,dx,plo", ng,ns,lo,hi,nu_m,dx,plo
!    print*, "... in drag_cic (lo_x, hi_x, ng): ", lo(1),hi(1),ng

    inv_dx = 1.0d0/dx
    inv_vol = inv_dx(1) * inv_dx(2) * inv_dx(3)

!     print*, "    checkpoint 1: inv_vol", inv_vol

    ! This will be incorrect if particle diameters are (~O) of dx
    do n = 1, np

       !qp = weights(n) * charge * inv_vol

       ! particle properties (10 is temp)
       vp(1) = particles(7,n)
       vp(2) = particles(8,n)
       vp(3) = particles(9,n)
       dia_p = particles(11,n)
       rho_p = particles(12,n)

!       print*, " "
!       print*, " PARTICLE No. ", n
!       print*, "==================="
!       print*, "    particle dump: 1 ", particles(1,n)
!       print*, "    particle dump: 2 ", particles(2,n)
!       print*, "    particle dump: 3 ", particles(3,n)
!       print*, "    particle dump: 4 ", particles(4,n)
!       print*, "    particle dump: 5 ", particles(5,n)
!       print*, "    particle dump: 6 ", particles(6,n)
!       print*, "    particle dump: 7 ", particles(7,n)
!       print*, "    particle dump: 8 ", particles(8,n)
!       print*, "    particle dump: 9 ", particles(9,n)
!       print*, "    particle dump: 10 ", particles(10,n)
!       print*, "    particle dump: 11 ", particles(11,n)
!       print*, "    particle dump: 12 ", particles(12,n)
!       print*, "    particle dump: 13 ", particles(13,n)
!       print*, "    particle dump: 14 ", particles(14,n)
!       print*, "    particle dump: 15 ", particles(15,n)
!       print*, "    particle dump: 16 ", particles(16,n)
!       print*, "    particle dump: 17 ", particles(17,n)
!       print*, "    checkpoint 2a: (vp)", vp
!       print*, "    checkpoint 2b: (dia_p, rho_p)", dia_p,rho_p

       ! scaled location of particle: (pos - lower corner of tile)*1/dx = x_p/dx (relative to tile)
       lx = (particles(1, n) - plo(1))*inv_dx(1)
       ly = (particles(2, n) - plo(2))*inv_dx(2)
       lz = (particles(3, n) - plo(3))*inv_dx(3)
!       print*, "    checkpoint 3: (li)", Lx,ly,lz

       ! get cell index which contains particle center
       i = floor(lx)
       j = floor(ly)
       k = floor(lz)
!       print*, "    checkpoint 3: (xi)", i,j,k

       ! distance(normalized by dx) from particle center to cell center (or face, not sure)
       wx_hi = lx - i
       wy_hi = ly - j
       wz_hi = lz - k

       ! takes on a weight with the 1/dx normalization, weight of force on low index
       wx_lo = 1.0d0 - wx_hi
       wy_lo = 1.0d0 - wy_hi
       wz_lo = 1.0d0 - wz_hi

       ! fluid velocity at particle
       v_ap(1:3) = vel(i,j,k,1:3)
!       print*, "    checkpoint 4: (v@p)", v_ap

       ! for simplified drag calc (Soldati 2009), should really interp vel to be at particle...
       vdiff = (v_ap(1)-vp(1))*(v_ap(1)-vp(1)) + &
             & (v_ap(2)-vp(2))*(v_ap(2)-vp(2)) + &
             & (v_ap(3)-vp(3))*(v_ap(3)-vp(3))
       vdiff = sqrt(vdiff)

       rho_f = rho(i,j,k)
       Re_p = vdiff * dia_p / nu_m
       tau_p = (rho_p/rho_f) * dia_p*dia_p/(18.0*nu_m)
       CT = 1.0 + 0.15*Re_p**0.687       

       ! sum of drags, signs flips because this is action of particle on fluid
       do oo = 1,3 

!          qp = rho_f * CT/tau_p * (vp(oo) - vel(i,j,k,oo)) * inv_vol ! not sure about units here, extra rho/dV...
          qp = rho_f * CT/tau_p * (vp(oo) - v_ap(oo)) !* inv_vol

          drg(i,   j,   k  ,oo) = drg(i,   j,   k  ,oo) + wx_lo*wy_lo*wz_lo*qp
          drg(i,   j,   k+1,oo) = drg(i,   j,   k+1,oo) + wx_lo*wy_lo*wz_hi*qp
          drg(i,   j+1, k  ,oo) = drg(i,   j+1, k  ,oo) + wx_lo*wy_hi*wz_lo*qp
          drg(i,   j+1, k+1,oo) = drg(i,   j+1, k+1,oo) + wx_lo*wy_hi*wz_hi*qp
          drg(i+1, j,   k  ,oo) = drg(i+1, j,   k  ,oo) + wx_hi*wy_lo*wz_lo*qp
          drg(i+1, j,   k+1,oo) = drg(i+1, j,   k+1,oo) + wx_hi*wy_lo*wz_hi*qp
          drg(i+1, j+1, k  ,oo) = drg(i+1, j+1, k  ,oo) + wx_hi*wy_hi*wz_lo*qp
          drg(i+1, j+1, k+1,oo) = drg(i+1, j+1, k+1,oo) + wx_hi*wy_hi*wz_hi*qp

       enddo

    end do

  end subroutine drag_cic_full


!========================================================================


end module




