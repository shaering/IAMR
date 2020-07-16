
#include <AMReX_FArrayBox.H>
#include <NavierStokesBase.H>
#include <AMReX_BLFort.H>
#include <PROB_NS_F.H>

using namespace amrex;

//
// Virtual access function for getting the forcing terms for the
// velocities and scalars.  The base version computes a buoyancy.
//
// As NavierStokesBase is currently implemented.  Velocities are integrated
// according to the equation
//
//     ui_t + uj ui_j = S_ui        ===> tforces = rho S_ui
//
// and scalars psi where (psi = rho q) as
//
//     psi_t + (uj psi)_j = S_psi   ===> tforces = S_psi = rho S_q
//
// q is a concentration.  This function returns a rho weighted
// source term, which requires a division by rho in the predict_velocity
// and velocity_advection routines.
//

// called for each FArrayBox in MultiFab within each level (this printed in basic advance) 

void
NavierStokesBase::getForce (FArrayBox&       force,
                            const Box&       bx,
                            int              ngrow,
                            int              scomp,
                            int              ncomp,
                            const Real       time,
                            const FArrayBox& Vel,
                            const FArrayBox& Scal,
                                  FArrayBox& rhs,
			    //                            const FArrayBox& Drag,
                            int              scalScomp,
                            int              level) //,
			    //                            int              scalDensity)
{
   if (ParallelDescriptor::IOProcessor() && getForceVerbose) {
      amrex::Print() << "NavierStokesBase::getForce(): Entered..." << std::endl 
                     << "time      = " << time << std::endl
                     << "scomp     = " << scomp << std::endl
                     << "ncomp     = " << ncomp << std::endl
                     << "ngrow     = " << ngrow << std::endl
                     << "scalScomp = " << scalScomp << std::endl;

   if (scomp==0) //                                                                      ncomp / scomp / forcing
     if  (ncomp==3) amrex::Print() << "Doing all velocity components" << std::endl;//      3       0       Vel
     else           amrex::Print() << "Doing incomplete vel components" << std::endl;//       X       0        ?
   else if (scomp==3) 
     if  (ncomp==1) amrex::Print() << "Doing density only" << std::endl;//         1       3       rho (not getting caught)
     else           amrex::Print() << "Doing all scalars" << std::endl;//          X       3    scalars(all)
   else if (scomp==4) amrex::Print() << "Doing tracer only" << std::endl;//        X       4      tracers
   else               amrex::Print() << "Doing individual scalar: "<< scomp << " " << std::endl;//  X       X    particular scalars
   
   }


   force.resize(grow(bx,ngrow),ncomp); // this needs to be adjusted for scalars (Temp)?

   int num_levels = Scal.size();
   const Real* VelDataPtr  = Vel.dataPtr(); 
   const Real* ScalDataPtr = Scal.dataPtr();
   const Real* rhoDataPtr  = Scal.dataPtr(scalScomp); // should be density, only 0 or Density gets passed
   const Real* uDataPtr    = Vel.dataPtr(0);
   const Real* vDataPtr    = Vel.dataPtr(1);
   const Real* wDataPtr    = Vel.dataPtr(2);

   const Real* dx       = geom.CellSize();
   const Real  grav     = gravity;
   const int*  f_lo     = force.loVect();
   const int*  f_hi     = force.hiVect();
   const int*  v_lo     = Vel.loVect();
   const int*  v_hi     = Vel.hiVect();
   const int*  s_lo     = Scal.loVect();
   const int*  s_hi     = Scal.hiVect();
   const int*  r_lo     = rhs.loVect();
   const int*  r_hi     = rhs.hiVect();
   const int   nscal    = NUM_SCALARS;

   if (ParallelDescriptor::IOProcessor() && getForceVerbose) {
#if (AMREX_SPACEDIM == 3)
      amrex::Print() << "NavierStokesBase::getForce()" << std::endl;
      amrex::Print() << "Force Domain:" << std::endl;
      amrex::Print() << "(" << f_lo[0] << "," << f_lo[1] << "," << f_lo[2] << ") - "
                     << "(" << f_hi[0] << "," << f_hi[1] << "," << f_hi[2] << ")" << std::endl;
      amrex::Print() << "Vel Domain:" << std::endl;
      amrex::Print() << "(" << v_lo[0] << "," << v_lo[1] << "," << v_lo[2] << ") - "
                     << "(" << v_hi[0] << "," << v_hi[1] << "," << v_hi[2] << ")" << std::endl;
      amrex::Print() << "Scal Domain:" << std::endl;
      amrex::Print() << "(" << s_lo[0] << "," << s_lo[1] << "," << s_lo[2] << ") - "
                     << "(" << s_hi[0] << "," << s_hi[1] << "," << s_hi[2] << ")" << std::endl;
      amrex::Print() << "rhs Domain:" << std::endl;
      amrex::Print() << "(" << r_lo[0] << "," << r_lo[1] << "," << r_lo[2] << ") - "
                     << "(" << r_hi[0] << "," << r_hi[1] << "," << r_hi[2] << ")" << std::endl;
#else
      amrex::Print() << "NavierStokesBase::getForce():" << std::endl;
      amrex::Print() << "Force Domain:" << std::endl;
      amrex::Print() << "(" << f_lo[0] << "," << f_lo[1] << ") - "
                     << "(" << f_hi[0] << "," << f_hi[1] << ")" << std::endl;
      amrex::Print() << "Vel Domain:" << std::endl;
      amrex::Print() << "(" << v_lo[0] << "," << v_lo[1] << ") - "
                     << "(" << v_hi[0] << "," << v_hi[1] << ")" << std::endl;
      amrex::Print() << "Scal Domain:" << std::endl;
      amrex::Print() << "(" << s_lo[0] << "," << s_lo[1] << ") - "
                     << "(" << s_hi[0] << "," << s_hi[1] << ")" << std::endl;
#endif

      Vector<Real> velmin(AMREX_SPACEDIM), velmax(AMREX_SPACEDIM);
      Vector<Real> scalmin(NUM_SCALARS), scalmax(NUM_SCALARS);

      for (int n=0; n<AMREX_SPACEDIM; n++) {
          velmin[n]= 1.e234;
          velmax[n]=-1.e234;
      }
      int ix = v_hi[0]-v_lo[0]+1;
      int jx = v_hi[1]-v_lo[1]+1;
#if (AMREX_SPACEDIM == 3)
      int kx = v_hi[2]-v_lo[2]+1;
      for (int k=0; k<kx; k++) {
#endif
         for (int j=0; j<jx; j++) {
            for (int i=0; i<ix; i++) {
               for (int n=0; n<AMREX_SPACEDIM; n++) {
#if (AMREX_SPACEDIM == 3)
                  int cell = ((n*kx+k)*jx+j)*ix+i;
#else
                  int cell = (n*jx+j)*ix+i;
#endif
                  Real v = VelDataPtr[cell];
                  if (v<velmin[n]) velmin[n] = v;
                  if (v>velmax[n]) velmax[n] = v;
               }
            }
         }
#if (AMREX_SPACEDIM == 3)
      }
#endif
      for (int n=0; n<AMREX_SPACEDIM; n++) 
         amrex::Print() << "Vel  " << n << " min/max " 
                        << velmin[n] << " / " << velmax[n] << std::endl;

      for (int n=0; n<NUM_SCALARS; n++) {
         scalmin[n]= 1.e234;
         scalmax[n]=-1.e234;
      }
      ix = s_hi[0]-s_lo[0]+1;
      jx = s_hi[1]-s_lo[1]+1;
#if (AMREX_SPACEDIM == 3)
      kx = s_hi[2]-s_lo[2]+1;
      for (int k=0; k<kx; k++) {
#endif
         for (int j=0; j<jx; j++) {
            for (int i=0; i<ix; i++) {
               for (int n=0; n<NUM_SCALARS; n++) {
#if (AMREX_SPACEDIM == 3)
                  int cell = ((n*kx+k)*jx+j)*ix+i;
#else
                  int cell = (n*jx+j)*ix+i;
#endif
                  Real s = ScalDataPtr[cell];
                  if (s<scalmin[n]) scalmin[n] = s;
                  if (s>scalmax[n]) scalmax[n] = s;
               }
            }
         }
#if (AMREX_SPACEDIM == 3)
      }
#endif
      for (int n=0; n<NUM_SCALARS; n++) 
         amrex::Print() << "Scal " << n << " min/max " << scalmin[n] 
                        << " / " << scalmax[n] << std::endl;
   }

   RealBox gridloc = RealBox(bx,geom.CellSize(),geom.ProbLo());



   //#ifdef AMREX_PARTICLES

   /*
   if (scomp==0) {
     std::cout << "Updating drag force (Sc,Nc,sSc)" << scomp << " " << ncomp << " " << scalScomp << "\n" ;
     theNSPC()->getDrag(rhs,Vel,Scal,visc_coef[0],ngrow,level);
     //     rhs.SumBoundary(0, ncomp, IntVect(1), Geom().periodicity());
     //     rhs.SumBoundary(Geom().periodicity());
   }
   else if (scomp==3) {
     //     std::cout << "Updating tempurature forcing (Sc,Nc,sSc)" << scomp << " " << ncomp << " " << scalScomp << "\n" ;
     //     theNSPC()->getTemp(rhs,Vel,Scal,visc_coef[0],ngrow,level);
   }
   else {
     std::cout << "Unknown forcing terms...\n";
   }

   std::cout << " *** get forcing okay\n";
   */

   //#endif


   // Here's the meat
   FORT_MAKEFORCE (&time,
                   BL_TO_FORTRAN_ANYD(force),
                   BL_TO_FORTRAN_ANYD(Vel),
		   //                   BL_TO_FORTRAN_N_ANYD(Scal),
                   BL_TO_FORTRAN_ANYD(Scal),
                   BL_TO_FORTRAN_ANYD(rhs),
                   dx,
                   gridloc.lo(),
                   gridloc.hi(),
                   &grav,
                   &Fx,&Fy,&Fz,
                   &scomp,&ncomp,&nscal,&getForceVerbose);

   if (ParallelDescriptor::IOProcessor() && getForceVerbose) {
      Vector<Real> forcemin(ncomp);
      Vector<Real> forcemax(ncomp);
      for (int n=0; n<ncomp; n++) {
         forcemin[n]= 1.e234;
         forcemax[n]=-1.e234;
      }

      int ix = f_hi[0]-f_lo[0]+1;
      int jx = f_hi[1]-f_lo[1]+1;
#if (AMREX_SPACEDIM == 3)
      int kx = f_hi[2]-f_lo[2]+1;
      for (int k=0; k<kx; k++) {
#endif
         for (int j=0; j<jx; j++) {
            for (int i=0; i<ix; i++) {
               for (int n=0; n<ncomp; n++) {
#if (AMREX_SPACEDIM == 3)
                  int cell = ((n*kx+k)*jx+j)*ix+i;
#else
                  int cell = (n*jx+j)*ix+i;
#endif
                  Real f = force.dataPtr()[cell];
                  if (f<forcemin[n]) forcemin[n] = f;
                  if (f>forcemax[n]) forcemax[n] = f;
               }
            }
         }
#if (AMREX_SPACEDIM == 3)
      }
#endif
      for (int n=0; n<ncomp; n++) 
         amrex::Print() << "Force " << n+scomp << " min/max " << forcemin[n] 
                        << " / " << forcemax[n] << std::endl;

      amrex::Print() << "NavierStokesBase::getForce(): Leaving..." 
                     << std::endl << "---" << std::endl;
   }
}
