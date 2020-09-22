#include <AMReX_TracerParticle_mod_K.H>
#include <AMReX_TracerParticles.H>
#include "AMReX_TracerParticles.H"
#include "AMReX_TracerParticle_mod_K.H"
#include <AMReX_Print.H>
#include <activeParticles.H>
#include "activeParticles.H"
#include "drag_F.H"
namespace amrex {

//
// Uses midpoint method to advance particles using umac.
//
void
ActiveParticleContainer::myAdvectWithUmac (MultiFab* umac, int lev, Real dt, MultiFab& rho, MultiFab& temp, Real nu_m, Real time)
{
    BL_PROFILE("AmrActiveParticleContainer::myAdvectWithUmac()");
    //    std::cout << " *** myAdvectWithUmac: lev " << lev << "\n";

    AMREX_ASSERT(OK(lev, lev, umac[0].nGrow()-1)); // segfault
    AMREX_ASSERT(lev >= 0 && lev < GetParticles().size());

    AMREX_D_TERM(AMREX_ASSERT(umac[0].nGrow() >= 1);,
                 AMREX_ASSERT(umac[1].nGrow() >= 1);,
                 AMREX_ASSERT(umac[2].nGrow() >= 1););

    AMREX_D_TERM(AMREX_ASSERT(!umac[0].contains_nan());,
                 AMREX_ASSERT(!umac[1].contains_nan());,
                 AMREX_ASSERT(!umac[2].contains_nan()););

    const Real      strttime = amrex::second();
    const Geometry& geom     = m_gdb->Geom(lev);
    const auto      plo      = geom.ProbLoArray();
    const auto      dxi      = geom.InvCellSizeArray();

    Vector<std::unique_ptr<MultiFab> > raii_umac(AMREX_SPACEDIM);
    Vector<MultiFab*> umac_pointer(AMREX_SPACEDIM);
    if (OnSameGrids(lev, umac[0]))
    {
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
	    umac_pointer[i] = &umac[i];
	}
    }
    else
    {
        for (int i = 0; i < AMREX_SPACEDIM; i++)
        {
	    int ng = umac[i].nGrow();
	    raii_umac[i].reset(new MultiFab(amrex::convert(m_gdb->ParticleBoxArray(lev),
                                                           IntVect::TheDimensionVector(i)),
					                   m_gdb->ParticleDistributionMap(lev),
					                   umac[i].nComp(), ng));


	    umac_pointer[i] = raii_umac[i].get();
	    umac_pointer[i]->copy(umac[i],0,0,umac[i].nComp(),ng,ng);
        }
    }

    //    std::cout << " *** myAdvectWithUmac: CP2 \n";

    
    
    // same thing for density
    Vector<std::unique_ptr<MultiFab> > raii_rho(1);
    Vector<MultiFab*> rho_pointer(1);
    //    FArrayBox rho = *rho_pointer;
    //    if (OnSameGrids( lev, rho[0] ))
    if (OnSameGrids( lev, rho ))
    {
      //        rho_pointer[0] = &rho[0];
        rho_pointer[0] = &rho;
    }
    else
    {
      //       int ng = rho[0].nGrow();
      int ng = rho.nGrow(); // is this it?
       raii_rho[0].reset(new MultiFab(amrex::convert(m_gdb->ParticleBoxArray(lev),
                                                     IntVect::TheDimensionVector(0)),
					             m_gdb->ParticleDistributionMap(lev),
					             rho.nComp(), ng));
				      //					             rho[0].nComp(), ng));

       rho_pointer[0] = raii_rho[0].get();
       //       rho_pointer[0]->copy(rho[0],0,0,rho[0].nComp(),ng,ng);
       rho_pointer[0]->copy(rho,0,0,rho.nComp(),ng,ng);
    }
    /**/

    
    // same thing for temp
    Vector<std::unique_ptr<MultiFab> > raii_temp(1);
    Vector<MultiFab*> temp_pointer(1);
    if (OnSameGrids( lev, temp ))
    {
        temp_pointer[0] = &temp;
    }
    else
    {
       int ng = temp.nGrow();
       raii_temp[0].reset(new MultiFab(amrex::convert(m_gdb->ParticleBoxArray(lev),
                                                     IntVect::TheDimensionVector(0)),
					             m_gdb->ParticleDistributionMap(lev),
					             temp.nComp(), ng));

       temp_pointer[0] = raii_temp[0].get();
       temp_pointer[0]->copy(temp,0,0,temp.nComp(),ng,ng);
    }
    /**/




    for (int ipass = 0; ipass < 2; ipass++)
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti)
	{
            int grid    = pti.index();
            auto& ptile = ParticlesAt(lev, pti);
            auto& aos  = ptile.GetArrayOfStructs();
            const int n = aos.numParticles();
            auto p_pbox = aos().data();
            const FArrayBox* fab[AMREX_SPACEDIM] = { AMREX_D_DECL(&((*umac_pointer[0])[grid]),
                                                                  &((*umac_pointer[1])[grid]),
                                                                  &((*umac_pointer[2])[grid])) };

            //array of these pointers to pass to the GPU
            amrex::GpuArray<amrex::Array4<const Real>, AMREX_SPACEDIM>
            const umacarr {AMREX_D_DECL((*fab[0]).array(),
                                        (*fab[1]).array(),
                                        (*fab[2]).array() )};


	    
            // same thing for density
            const FArrayBox* rfab[AMREX_SPACEDIM] = { AMREX_D_DECL(&((*rho_pointer[0])[grid]),
                                                                   &((*rho_pointer[0])[grid]),
                                                                   &((*rho_pointer[0])[grid])) };
            amrex::GpuArray<amrex::Array4<const Real>, AMREX_SPACEDIM>
	    const rhoarr {AMREX_D_DECL((*rfab[0]).array(),
                                       (*rfab[0]).array(),
                                       (*rfab[0]).array() )};
	    /**/

	    
            // same thing for temp
            const FArrayBox* tfab[AMREX_SPACEDIM] = { AMREX_D_DECL(&((*temp_pointer[0])[grid]),
                                                                   &((*temp_pointer[0])[grid]),
                                                                   &((*temp_pointer[0])[grid])) };
            amrex::GpuArray<amrex::Array4<const Real>, AMREX_SPACEDIM>
	    const temparr {AMREX_D_DECL((*tfab[0]).array(),
                                        (*tfab[0]).array(),
                                        (*tfab[0]).array() )};
	    /**/


            amrex::ParallelFor(n,
                               [=] AMREX_GPU_DEVICE (int i)
            {

                ParticleType& p = p_pbox[i];
                if (p.id() <= 0) return;
                Real v[AMREX_SPACEDIM];
                mac_interpolate(p, plo, dxi, umacarr, v);
                //std::cout << "     Mac_inerp vel okay"; //<< i << ":" << v[0] << " " << v[1] << " " << v[2] << " \n";

		
                // same thing for density
                //ParticleType& p = p_pbox[i];
                Real rho_f[AMREX_SPACEDIM];
		//std::cout << "particle pos: " << p.pos(0) << " " << p.pos(1) << " " << p.pos(2);
		//		std::cout << "box: " << p_pbox->loVect() << " " << p_pbox->hiVect();
                mac_interpolate(p, plo, dxi, rhoarr, rho_f); //segfault is apparently due to rhoarr being too large for p
                //std::cout << "     Mac_inerp rho okay"; //<< i << ":" << rho_f[0] << " " << " \n";
		/**/

		
                // same thing for temp
                Real temp_f[AMREX_SPACEDIM];
                mac_interpolate(p, plo, dxi, temparr, temp_f);
		/**/


                // soldati 2009
  	        Real dia_p = p.m_rdata.arr[AMREX_SPACEDIM+7]; // make these indices general
  	        Real rho_p = p.m_rdata.arr[AMREX_SPACEDIM+8];
                Real vdiff;
                vdiff = (v[0]-p.m_rdata.arr[2*AMREX_SPACEDIM+0])*(v[0]-p.m_rdata.arr[2*AMREX_SPACEDIM+0])
                      + (v[1]-p.m_rdata.arr[2*AMREX_SPACEDIM+1])*(v[1]-p.m_rdata.arr[2*AMREX_SPACEDIM+1])
                      + (v[2]-p.m_rdata.arr[2*AMREX_SPACEDIM+2])*(v[2]-p.m_rdata.arr[2*AMREX_SPACEDIM+2]);
                Real Re_p = sqrt(vdiff) * dia_p / nu_m;
                Real tau_p = (rho_p/rho_f[0])*dia_p*dia_p/(18.0*nu_m); // need to pass nu and rho for fluid (simple with static first)
                //Real tau_p = (rho_p/1.0)*dia_p*dia_p/(18.0*nu_m);
  	        Real CT = 1.0 + 0.15*std::pow(Re_p,0.687);
		Real wt0 = (CT/tau_p)*dt;
                wt0 = std::min(wt0,1.0);


		// add more slots in particle array and store momentum src terms for future interp? (simplifies drag_cic)



		/*		
                std::cout << " *** myAdvectWithUmac: CP5 \n";
                std::cout << "     Data (nu_m, d, rho_p, rho_f, T_f, vdiff, Re_p, tau_p, CT, wto) " 
                          << nu_m << " " 
                          << dia_p << " " 
                          << rho_p << " " 
                          << rho_f[0] << " " 
                          << temp_f[0] << " "
                          << vdiff << " " 
                          << Re_p << " " 
                          << tau_p << " " 
                          << CT << " " 
                          << wt0 << "\n";

		std::cout << " p.m_Xdata.arr dump (PRIOR) \n";
		std::cout << " ================== \n";
		std::cout << "(xp)" << p.m_rdata.pos[0] << "\n";
		std::cout << "(yp)" << p.m_rdata.pos[1] << "\n";
		std::cout << "(zp)" << p.m_rdata.pos[2] << "\n";
		std::cout << "(xp?)" << p.m_rdata.arr[0] << "\n";
		std::cout << "(yp?)" << p.m_rdata.arr[1] << "\n";
		std::cout << "(zp?)" << p.m_rdata.arr[2] << "\n";
		std::cout << "(*)" << p.m_rdata.arr[3] << "\n";
		std::cout << "(*)" << p.m_rdata.arr[4] << "\n";
		std::cout << "(*)" << p.m_rdata.arr[5] << "\n";
		std::cout << "(up)" << p.m_rdata.arr[6] << "\n";
		std::cout << "(vp)" << p.m_rdata.arr[7] << "\n";
		std::cout << "(wp)" << p.m_rdata.arr[8] << "\n";
		std::cout << "(Tp)" << p.m_rdata.arr[9] << "\n";
		std::cout << "(dp)" << p.m_rdata.arr[10] << "\n";
		std::cout << "(rhop)" << p.m_rdata.arr[11] << "\n";
		std::cout << "(Fxp)" << p.m_rdata.arr[12] << "\n";
		std::cout << "(Fyp)" << p.m_rdata.arr[13] << "\n";
		std::cout << "(Fzp)" <<  p.m_rdata.arr[14] << "\n";
		std::cout << "(FTp)" <<  p.m_rdata.arr[15] << "\n";
		*/

		
/*
Predictor/corrector update of the position. Velocity field in the IAMR case 
is face-centroid based (in PeleC, it is cell-centroid based).

1. ipass=0, the velocity is interpolated to the “old” particle position, 
and that velocity is then used to move the particles for dt/2

2. ipass=1,  the velocity is then interpolated to that intermediate location, and 
that velocity is then used to recompute the full advance, x = x + vmid*dt

Thus for each particle, you need a temporarary to hold the old position so 
that you set the position of the velocity interp on the second pass
*/

//  This makes sense but is not what is actually done, what is done makes no sense
//  is the position stored in rdata_arr used for interpolation?

/*
p.m_rdata_arr[1:SPACEDIM] = swap vector
p.m_rdata_arr[SPACEDIM+1:SPACEDIM+3] = particle velocity
p.m_rdata_arr[SPACEDIM+1:SPACEDIM+4] = particle temp
p.m_rdata_arr[SPACEDIM+5] = particle diameter
p.m_rdata_arr[SPACEDIM+6] = particle density
*/

/*
                for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                 {
                   p.m_rdata.arr[AMREX_SPACEDIM+dim] = p.m_rdata.pos[dim]; // copy old pos, not necessary any more
                   p.m_rdata.arr[2*AMREX_SPACEDIM+dim] = (1.0-wt0)*p.m_rdata.arr[2*AMREX_SPACEDIM+dim] + wt0*v[dim]; // update vel_p
		   p.m_rdata.pos[dim] += dt*p.m_rdata.arr[2*AMREX_SPACEDIM+dim];                                     // update pos
                 }
*/



                // imp-eff specific
		int hole_top_n0, hole_top_n1;
		int hole_bot_n0, hole_bot_n1;		
		const Real ground = 0.5;
	        const Real ceiling = 2.5;
	        const Real bottom = -1.49;		
        	Real r_top, r_bot, rp_top, rp_bot;
	        Real x_top, z_top, x_bot, z_bot;
		Real xp, yp, zp, norm_x,norm_z;
		Real ur_prime, ut_prime, theta_p;

		x_top = 2.5;
		z_top = 1.0;		
		r_top = 0.5;

		x_bot = 1.5;
		z_bot = 1.0;		
		r_bot = 0.5;

                // index i is particle number, gradually release particles into jet inlet
		int gradual_release_flag, p_cutoff, p_freq;
		Real t_freq;
                Real rnd_wgt;		

		t_freq = 0.1;
		gradual_release_flag = 1; // set this to zero except from gradual release in imp-eff
		p_cutoff = std::round(time/t_freq);
		//		std::cout << "*** particle: " << i << "\n";
		//		std::cout << "*** time: " << time << "\n";
		//		std::cout << "*** p_cutoff: " << p_cutoff << "\n";

		//                int p_active = p.m_idata[0];

		// is i local to the proc?  would explain freezing...
		if ( i<p_cutoff || gradual_release_flag==0 ) {
		  
		  //		  std::cout << "*** particle: " << i << " lives! (" << time << ")\n";
		
                // particles advanced here <warp>
                if (ipass == 0)
                {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {

                      // original
		      //p.m_rdata.arr[AMREX_SPACEDIM+dim] = p.m_rdata.pos[dim]; // copy old pos
		      //p.m_rdata.pos[dim] += 0.5*dt*v[dim];                    // update pos to dt/2

                      p.m_rdata.arr[AMREX_SPACEDIM+dim] = p.m_rdata.pos[dim];                                   // copy old pos
		      p.m_rdata.pos[dim] += 0.5*dt*((1.0-wt0)*p.m_rdata.arr[2*AMREX_SPACEDIM+dim]+wt0*v[dim]);  // update pos to dt/2

                    }
                }
                else
                {


		    // old position
		    xp = p.m_rdata.arr[AMREX_SPACEDIM+0]; //p.m_rdata.pos[0];
		    yp = p.m_rdata.arr[AMREX_SPACEDIM+1]; //p.m_rdata.pos[1];
		    zp = p.m_rdata.arr[AMREX_SPACEDIM+2]; //p.m_rdata.pos[2];
      	            rp_top = sqrt( (xp-x_top)*(xp-x_top) + (zp-z_top)*(zp-z_top) );
	            rp_bot = sqrt( (xp-x_bot)*(xp-x_bot) + (zp-z_bot)*(zp-z_bot) );

		    hole_top_n0 = 0;
		    hole_bot_n0 = 0;		
		    if ( rp_top<=r_top && yp>=ceiling ) hole_top_n0 = 1;
		    if ( rp_bot<=r_bot && yp<=ground )  hole_bot_n0 = 1;
		
				  
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {

                      // original
		      // p.m_rdata.pos[dim] = p.m_rdata.arr[AMREX_SPACEDIM+dim] + dt*v[dim]; // update pos to full dt
		      // p.m_rdata.arr[AMREX_SPACEDIM+dim] = v[dim];                         // copy vel to arr, why?

                      p.m_rdata.arr[2*AMREX_SPACEDIM+dim] = (1.0-wt0)*p.m_rdata.arr[2*AMREX_SPACEDIM+dim] + wt0*v[dim]; // update vel_p
		      p.m_rdata.pos[dim] = p.m_rdata.arr[AMREX_SPACEDIM+dim] + dt*p.m_rdata.arr[2*AMREX_SPACEDIM+dim];  // update pos to full dt
		      p.m_rdata.arr[AMREX_SPACEDIM+dim] = p.m_rdata.arr[2*AMREX_SPACEDIM+dim];                          // copy vel to arr, why?

                      // update force to reuse
    		      p.m_rdata.arr[2*AMREX_SPACEDIM+dim+6] = CT/tau_p * (v[dim]-p.m_rdata.arr[2*AMREX_SPACEDIM+dim]);

                    }

		    
                    // update particle temp, make alpha readable
    		    p.m_rdata.arr[9] += 0.01 * dt * (temp_f[0]-p.m_rdata.arr[9]);
    		    p.m_rdata.arr[2*AMREX_SPACEDIM+3+6] = 0.01 * (temp_f[0]-p.m_rdata.arr[9]);
		    /**/

		    // hack for reflection or sticking ==> this should be changed to some global wall distance, manual changes now 
		    //srand( (unsigned)time( NULL ) );
		    //rnd_wgt = round(rand()/RAND_MAX);
 	            rnd_wgt = 1.0;		    

		    // updated position
                    xp = p.m_rdata.pos[0];
                    yp = p.m_rdata.pos[1];
                    zp = p.m_rdata.pos[2];
  		    rp_top = sqrt( (xp-x_top)*(xp-x_top) + (zp-z_top)*(zp-z_top) );
		    rp_bot = sqrt( (xp-x_bot)*(xp-x_bot) + (zp-z_bot)*(zp-z_bot) );

		    hole_top_n1 = 0;
 		    hole_bot_n1 = 0;		
		    if ( rp_top<=r_top && yp>=ceiling ) hole_top_n1 = 1;
		    if ( rp_bot<=r_bot && yp<=ground )  hole_bot_n1 = 1;		    

		    // floor bounce
                    if (yp < ground && hole_bot_n1==0 && hole_bot_n0==0 ) {
		       p.m_rdata.pos[1] = ground;
                       p.m_rdata.arr[2*AMREX_SPACEDIM+0] =  1.0 * rnd_wgt * p.m_rdata.arr[2*AMREX_SPACEDIM+0];
                       p.m_rdata.arr[2*AMREX_SPACEDIM+1] = -1.0 * rnd_wgt * p.m_rdata.arr[2*AMREX_SPACEDIM+1];
                       p.m_rdata.arr[2*AMREX_SPACEDIM+2] =  1.0 * rnd_wgt * p.m_rdata.arr[2*AMREX_SPACEDIM+2];
	            }

		    // ceiling bounce
                    if (yp > ceiling && hole_top_n1==0 && hole_top_n0==0) {
		       p.m_rdata.pos[1] = ceiling;
                       p.m_rdata.arr[2*AMREX_SPACEDIM+0] =  1.0 * rnd_wgt * p.m_rdata.arr[2*AMREX_SPACEDIM+0];
                       p.m_rdata.arr[2*AMREX_SPACEDIM+1] = -1.0 * rnd_wgt * p.m_rdata.arr[2*AMREX_SPACEDIM+1];
                       p.m_rdata.arr[2*AMREX_SPACEDIM+2] =  1.0 * rnd_wgt * p.m_rdata.arr[2*AMREX_SPACEDIM+2];
		    }

		    
		    // upper hole
                    if (yp > ceiling && hole_top_n1==0 && hole_top_n0==1) {

		       // rotate vel
		       norm_x = (xp - x_top) / rp_top;
		       norm_z = (zp - z_top) / rp_top;
		       theta_p = atan(norm_z/norm_x);

		       ur_prime = p.m_rdata.arr[2*AMREX_SPACEDIM+0] * cos(theta_p) + p.m_rdata.arr[2*AMREX_SPACEDIM+2] * sin(theta_p);
		       ut_prime = -1.0 * p.m_rdata.arr[2*AMREX_SPACEDIM+0] * sin(theta_p) + p.m_rdata.arr[2*AMREX_SPACEDIM+2] * cos(theta_p);

		       // flip wall vel
		       ur_prime = -ur_prime; 

		       // rotate back
		       p.m_rdata.arr[2*AMREX_SPACEDIM+0] = ur_prime * cos(theta_p) - ut_prime * sin(theta_p);
                       p.m_rdata.arr[2*AMREX_SPACEDIM+2] = ur_prime * sin(theta_p) - ut_prime * cos(theta_p);

		       // back particle out to hole surface		       
                       p.m_rdata.pos[0] = r_top * norm_x;
                       p.m_rdata.pos[2] = r_top * norm_z;		       

	            }

		    // lower hole
                    if (yp < ground && hole_bot_n1==0 && hole_bot_n0==1) {

		       // rotate vel
		       norm_x = (xp - x_bot) / rp_bot;
		       norm_z = (zp - z_bot) / rp_bot;
		       theta_p = atan(norm_z/norm_x);

		       ur_prime = p.m_rdata.arr[2*AMREX_SPACEDIM+0] * cos(theta_p) + p.m_rdata.arr[2*AMREX_SPACEDIM+2] * sin(theta_p);
		       ut_prime = -1.0 * p.m_rdata.arr[2*AMREX_SPACEDIM+0] * sin(theta_p) + p.m_rdata.arr[2*AMREX_SPACEDIM+2] * cos(theta_p);

		       // flip wall vel
		       ur_prime = -ur_prime; 

		       // rotate back
		       p.m_rdata.arr[2*AMREX_SPACEDIM+0] = ur_prime * cos(theta_p) - ut_prime * sin(theta_p);
                       p.m_rdata.arr[2*AMREX_SPACEDIM+2] = ur_prime * sin(theta_p) - ut_prime * cos(theta_p);

		       // back particle out to hole surface
                       p.m_rdata.pos[0] = r_bot * norm_x;
                       p.m_rdata.pos[2] = r_bot * norm_z;		       

	            }

		    // recycle
                    if (yp < bottom) {

		       std::cout << " *** Particle Recycled: " << i << "\n";
		       p.m_rdata.pos[0] = p.m_rdata.pos[0] + 1.0; // fix these hard-codes
                       p.m_rdata.pos[1] = p.m_rdata.pos[1] + 4.48;
                       p.m_rdata.pos[2] = p.m_rdata.pos[2];

	            }		    		    

		    
                }

		}
		else { // just set forces to zero and dont update position, vel, or temp

    		      p.m_rdata.arr[2*AMREX_SPACEDIM+0+6] = 0.0;
    		      p.m_rdata.arr[2*AMREX_SPACEDIM+1+6] = 0.0;
    		      p.m_rdata.arr[2*AMREX_SPACEDIM+2+6] = 0.0;		      
		  
		}; // gradual release if

		/*
		std::cout << " p.m_Xdata.arr dump (AFTER) \n";
		std::cout << " ================== \n";
		std::cout << "(xp)" << p.m_rdata.pos[0] << "\n";
		std::cout << "(yp)" << p.m_rdata.pos[1] << "\n";
		std::cout << "(zp)" << p.m_rdata.pos[2] << "\n";
		std::cout << "(xp?)" << p.m_rdata.arr[0] << "\n";
		std::cout << "(yp?)" << p.m_rdata.arr[1] << "\n";
		std::cout << "(zp?)" << p.m_rdata.arr[2] << "\n";
		std::cout << "(*)" << p.m_rdata.arr[3] << "\n";
		std::cout << "(*)" << p.m_rdata.arr[4] << "\n";
		std::cout << "(*)" << p.m_rdata.arr[5] << "\n";
		std::cout << "(up)" << p.m_rdata.arr[6] << "\n";
		std::cout << "(vp)" << p.m_rdata.arr[7] << "\n";
		std::cout << "(wp)" << p.m_rdata.arr[8] << "\n";
		std::cout << "(Tp)" << p.m_rdata.arr[9] << "\n";
		std::cout << "(dp)" << p.m_rdata.arr[10] << "\n";
		std::cout << "(rhop)" << p.m_rdata.arr[11] << "\n";
		std::cout << "(Fxp)" << p.m_rdata.arr[12] << "\n";
		std::cout << "(Fyp)" << p.m_rdata.arr[13] << "\n";
		std::cout << "(Fzp)" <<  p.m_rdata.arr[14] << "\n";
		std::cout << "(FTp)" <<  p.m_rdata.arr[15] << "\n";
		*/
            });

        }

	//            std::cout << " *** myAdvectWithUmac: CP8 \n";
    }

    //            std::cout << " *** myAdvectWithUmac: CP9 \n";

    if (m_verbose > 1)
    {
        Real stoptime = amrex::second() - strttime;

#ifdef AMREX_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
                ParallelReduce::Max(stoptime, ParallelContext::IOProcessorNumberSub(),
                                    ParallelContext::CommunicatorSub());

                amrex::Print() << "TracerParticleContainer::AdvectWithUmac() time: " << stoptime << '\n';
#ifdef AMREX_LAZY
	});
#endif
    }
}


// Tracer version
void
ActiveParticleContainer::basicAdvectWithUmac (MultiFab* umac, int lev, Real dt)
//TracerParticleContainer::AdvectWithUmac (MultiFab* umac, int lev, Real dt, Real Re_p, Real tau_p)
{
    BL_PROFILE("TracerParticleContainer::AdvectWithUmac()");
    //    std::cout << "lev: " << lev << "\n";
    //    std::cout << "umac[0].nGrow()-1: " << umac[0].nGrow()-1 << "\n";
    AMREX_ASSERT(OK(lev, lev, umac[0].nGrow()-1)); // segfault
    AMREX_ASSERT(lev >= 0 && lev < GetParticles().size());

    AMREX_D_TERM(AMREX_ASSERT(umac[0].nGrow() >= 1);,
                 AMREX_ASSERT(umac[1].nGrow() >= 1);,
                 AMREX_ASSERT(umac[2].nGrow() >= 1););

    AMREX_D_TERM(AMREX_ASSERT(!umac[0].contains_nan());,
                 AMREX_ASSERT(!umac[1].contains_nan());,
                 AMREX_ASSERT(!umac[2].contains_nan()););

    const Real      strttime = amrex::second();
    const Geometry& geom     = m_gdb->Geom(lev);
    const auto          plo      = geom.ProbLoArray();
    const auto          dxi      = geom.InvCellSizeArray();

    Vector<std::unique_ptr<MultiFab> > raii_umac(AMREX_SPACEDIM);
    Vector<MultiFab*> umac_pointer(AMREX_SPACEDIM);
    if (OnSameGrids(lev, umac[0]))
    {
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
	    umac_pointer[i] = &umac[i];
	}
    }
    else
    {
        for (int i = 0; i < AMREX_SPACEDIM; i++)
        {
	    int ng = umac[i].nGrow();
	    raii_umac[i].reset(new MultiFab(amrex::convert(m_gdb->ParticleBoxArray(lev),
                                                           IntVect::TheDimensionVector(i)),

					                   m_gdb->ParticleDistributionMap(lev),
					                   umac[i].nComp(), ng));


	    umac_pointer[i] = raii_umac[i].get();
	    umac_pointer[i]->copy(umac[i],0,0,umac[i].nComp(),ng,ng);
        }
    }

    for (int ipass = 0; ipass < 2; ipass++)
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti)
        {
            int grid    = pti.index();
            auto& ptile = ParticlesAt(lev, pti);
            auto& aos  = ptile.GetArrayOfStructs();
            const int n = aos.numParticles();
            auto p_pbox = aos().data();
            const FArrayBox* fab[AMREX_SPACEDIM] = { AMREX_D_DECL(&((*umac_pointer[0])[grid]),
                                                                  &((*umac_pointer[1])[grid]),
                                                                  &((*umac_pointer[2])[grid])) };



            //array of these pointers to pass to the GPU
            amrex::GpuArray<amrex::Array4<const Real>, AMREX_SPACEDIM>
            const umacarr {AMREX_D_DECL((*fab[0]).array(),
                                        (*fab[1]).array(),
                                        (*fab[2]).array() )};

            amrex::ParallelFor(n,
                               [=] AMREX_GPU_DEVICE (int i)
            {

                ParticleType& p = p_pbox[i];
                if (p.id() <= 0) return;
                Real v[AMREX_SPACEDIM];
                //Real vp[AMREX_SPACEDIM];
                //Real pos_lcl[AMREX_SPACEDIM];
                mac_interpolate(p, plo, dxi, umacarr, v);

                // soldati 2009
		//  	        Real Re_p = p.m_rdata.arr[AMREX_SPACEDIM+4];
		//  	        Real tau_p = p.m_rdata.arr[AMREX_SPACEDIM+5];
  	        //Real CT = 1.0 + 0.15*pow(Re_p,0.687);
		//Real wt0 = (CT/tau_p)*dt;
                //wt0 = min(wt0,1.0);

                // particles advanced here <warp>
                if (ipass == 0) // wtf is ipass?
                {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {

                      // original
		      p.m_rdata.arr[AMREX_SPACEDIM+dim] = p.m_rdata.pos[dim]; // copy pos at n-1/2 to arr
		      p.m_rdata.pos[dim] += 0.5*dt*v[dim];                    // update pos for n

		      //                      pos_lcl[dim] = p.m_rdata.pos[dim];
		      //		      p.m_rdata.arr[AMREX_SPACEDIM+dim] = (1.0-wt0)*p.m_rdata.arr[AMREX_SPACEDIM+dim] + wt0*v[dim];
		      //                      p.m_rdata.pos[dim] += 0.5*dt*(1.0-0.5*wt0)*p.m_rdata.arr[AMREX_SPACEDIM+dim] + 0.5*wt0*v[dim]; 

                      // w/o adding fields
		      //vp[dim] = p.m_rdata.arr[AMREX_SPACEDIM+dim];  // copy vel_p
                      //vp[dim] = dt*(1.0-wt0)*vp[dim] + wt0*v[dim]; // update vel_p
		      //p.m_rdata.arr[AMREX_SPACEDIM+dim] = p.m_rdata.pos[dim];  // copy pos at n-1/2 to arr
		      //p.m_rdata.pos[dim] += 0.5*dt*vp[dim];                    // update pos for n

		      //                      p.m_rdata.arr[2*AMREX_SPACEDIM+dim] = (1.0-wt0)*p.m_rdata.arr[2*AMREX_SPACEDIM+dim] + wt0*v[dim]; // update vel_p
		      //		      p.m_rdata.pos[dim] += 0.5*dt*p.m_rdata.arr[2*AMREX_SPACEDIM+dim];                    // update pos for n

		      //		      std::cout << "ipass 0: " << dim << " " << p.m_rdata.arr[AMREX_SPACEDIM+dim] << " " << p.m_rdata.pos[dim] << " " << v[dim] <<"\n";

                    }
                }
                else
                {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {


                      // original
		      p.m_rdata.pos[dim]  = p.m_rdata.arr[AMREX_SPACEDIM+dim] + dt*v[dim]; // update pos to n+1/2
		      p.m_rdata.arr[AMREX_SPACEDIM+dim] = v[dim];                          // copy vel to arr, why?

		      //                      p.m_rdata.pos[dim] = pos_lcl[dim] + dt*p.m_rdata.arr[AMREX_SPACEDIM+dim];
		      //p.m_rdata.arr[AMREX_SPACEDIM+dim] = (1.0-wt0)*p.m_rdata.arr[AMREX_SPACEDIM+dim] + wt0*v[dim]; //already updated

                      // w/o adding fields
		      //p.m_rdata.pos[dim]  = p.m_rdata.arr[AMREX_SPACEDIM+dim] + dt*vp[dim]; // update pos to n+1/2
		      //p.m_rdata.arr[AMREX_SPACEDIM+dim] = vp[dim];                          // copy vel to arr

		      //		      p.m_rdata.pos[dim]  = p.m_rdata.arr[AMREX_SPACEDIM+dim] + dt*p.m_rdata.arr[2*AMREX_SPACEDIM+dim]; // update pos to n+1/2
		      //		      p.m_rdata.arr[AMREX_SPACEDIM+dim] = p.m_rdata.arr[2*AMREX_SPACEDIM+dim];                          // copy vel to arr, why?

		      //		      std::cout << "ipass 1: " << dim << " " << p.m_rdata.arr[AMREX_SPACEDIM+dim] << " " << p.m_rdata.pos[dim] << " " << v[dim] <<"\n";

                    }
                }
            });
        }
    }

    if (m_verbose > 1)
    {
        Real stoptime = amrex::second() - strttime;

#ifdef AMREX_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
                ParallelReduce::Max(stoptime, ParallelContext::IOProcessorNumberSub(),
                                    ParallelContext::CommunicatorSub());

                amrex::Print() << "TracerParticleContainer::AdvectWithUmac() time: " << stoptime << '\n';
#ifdef AMREX_LAZY
	});
#endif
    }
}



//
// Uses midpoint method to advance particles using cell-centered velocity
//
void
ActiveParticleContainer::myAdvectWithUcc (const MultiFab& Ucc, int lev, Real dt)
{
    BL_PROFILE("AmrActiveParticleContainer::myAdvectWithUcc()");
    AMREX_ASSERT(Ucc.nGrow() > 0);
    AMREX_ASSERT(OK(lev, lev, Ucc.nGrow()-1));
    AMREX_ASSERT(lev >= 0 && lev < GetParticles().size());
    AMREX_ASSERT(!Ucc.contains_nan());

    const Real          strttime = amrex::second();
    const Geometry&     geom     = m_gdb->Geom(lev);
    const auto          plo      = geom.ProbLoArray();
    const auto          dxi      = geom.InvCellSizeArray();

    AMREX_ASSERT(OnSameGrids(lev, Ucc));

    for (int ipass = 0; ipass < 2; ipass++)
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti)
        {
            int grid    = pti.index();
            auto& ptile = ParticlesAt(lev, pti);
            auto& aos  = ptile.GetArrayOfStructs();
            const int n          = aos.numParticles();
            const FArrayBox& fab = Ucc[grid];
            const auto uccarr = fab.array();
            auto  p_pbox = aos().data();

            amrex::ParallelFor(n,
                               [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p  = p_pbox[i];
                if (p.id() <= 0) return;
                Real v[AMREX_SPACEDIM];

                cic_interpolate(p, plo, dxi, uccarr, v);

                if (ipass == 0)
                {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.rdata(dim) = p.pos(dim);
                        p.pos(dim) += 0.5*dt*v[dim];
                    }
                }
                else
                {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.rdata(dim) = p.rdata(dim) + dt*v[dim];
                        p.rdata(dim) = v[dim];
                    }
                }
            });
        }
    }

    if (m_verbose > 1)
    {
        Real stoptime = amrex::second() - strttime;

#ifdef AMREX_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
                ParallelReduce::Max(stoptime, ParallelContext::IOProcessorNumberSub(),
                                    ParallelContext::CommunicatorSub());

                amrex::Print() << "AmrActiveParticleContainer::myAdvectWithUcc() time: " << stoptime << '\n';
#ifdef AMREX_LAZY
            });
#endif
    }
}

void
ActiveParticleContainer::Timestamp (const std::string&      basename,
				    const MultiFab&         mf,
				    int                     lev,
				    Real                    time,
				    const std::vector<int>& indices)
{
    BL_PROFILE("AmrActiveParticleContainer::Timestamp()");
    //
    // basename -> base filename for the output file
    // mf       -> the multifab
    // lev      -> level to check for particles
    // time     -> simulation time (will be recorded in Timestamp file)
    // indices  -> indices into mf that we output
    //
    AMREX_ASSERT(lev >= 0);
    AMREX_ASSERT(time >= 0);
    AMREX_ASSERT(!basename.empty());
    AMREX_ASSERT(lev <= m_gdb->finestLevel());

    const Real strttime = amrex::second();

    const int   MyProc    = ParallelDescriptor::MyProc();
    const int   NProcs    = ParallelContext::NProcsSub();
    // We'll spread the output over this many files.
    int nOutFiles(64);
    ParmParse pp("particles");
    pp.query("particles_nfiles",nOutFiles);
    if(nOutFiles == -1) {
      nOutFiles = NProcs;
    }
    nOutFiles = std::max(1, std::min(nOutFiles,NProcs));
    const int   nSets     = ((NProcs + (nOutFiles - 1)) / nOutFiles);
    const int   mySet     = (MyProc / nOutFiles);

    for (int iSet = 0; iSet < nSets; ++iSet)
      {
        if (mySet == iSet)
	  {
            //
            // Do we have any particles at this level that need writing?
            //
            bool gotwork = false;

            const auto& pmap = GetParticles(lev);
	    for (auto& kv : pmap) {
              const auto& pbox = kv.second.GetArrayOfStructs();
	      for (int k = 0; k < pbox.numParticles(); ++k)
	      {
		const ParticleType& p = pbox[k];
		if (p.id() > 0) {
		  gotwork = true;
		  break;
		}
	      }
	      if (gotwork) break;
	    }

            if (gotwork)
	      {
                std::string FileName = amrex::Concatenate(basename + '_', MyProc % nOutFiles, 2);

                VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

                std::ofstream TimeStampFile;

                TimeStampFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

                TimeStampFile.open(FileName.c_str(), std::ios::out|std::ios::app|std::ios::binary);

                TimeStampFile.setf(std::ios_base::scientific,std::ios_base::floatfield);

                TimeStampFile.precision(10);

                TimeStampFile.seekp(0, std::ios::end);

                if (!TimeStampFile.good())
                    amrex::FileOpenFailed(FileName);

                const int       M  = indices.size();
                const BoxArray& ba = mf.boxArray();

                std::vector<Real> vals(M);

		for (auto& kv : pmap) {
		  int grid = kv.first.first;
		  const auto& pbox = kv.second.GetArrayOfStructs();
		  const Box&       bx   = ba[grid];
		  const FArrayBox& fab  = mf[grid];

		  for (int k = 0; k < pbox.numParticles(); ++k)
		    {
		      const ParticleType& p = pbox[k];

		      if (p.id() <= 0) continue;

		      const IntVect& iv = Index(p,lev);

		      if (!bx.contains(iv) && !ba.contains(iv)) continue;

		      TimeStampFile << p.id()  << ' ' << p.cpu() << ' ';

		      AMREX_D_TERM(TimeStampFile << p.pos(0) << ' ';,
                                   TimeStampFile << p.pos(1) << ' ';,
                                   TimeStampFile << p.pos(2) << ' ';);

		      TimeStampFile << time;
		      //
		      // AdvectWithUmac stores the velocity in rdata ...
		      //
		      AMREX_D_TERM(TimeStampFile << ' ' << p.rdata(0);,
                                   TimeStampFile << ' ' << p.rdata(1);,
                                   TimeStampFile << ' ' << p.rdata(2););

		      if (M > 0)
                        {
			  ParticleType::Interp(p,m_gdb->Geom(lev),fab,&indices[0],&vals[0],M);

			  for (int i = 0; i < M; i++)
                            {
			      TimeStampFile << ' ' << vals[i];
                            }
                        }

		      TimeStampFile << '\n';
                    }
                }

                TimeStampFile.flush();
                TimeStampFile.close();
            }

            const int iBuff     = 0;
            const int wakeUpPID = (MyProc + nOutFiles);
            const int tag       = (MyProc % nOutFiles);

            if (wakeUpPID < NProcs)
                ParallelDescriptor::Send(&iBuff, 1, wakeUpPID, tag);
        }
        if (mySet == (iSet + 1))
        {
            //
            // Next set waits.
            //
            int       iBuff;
            const int waitForPID = (MyProc - nOutFiles);
            const int tag        = (MyProc % nOutFiles);

            ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
        }
    }

    if (m_verbose > 1)
    {
        Real stoptime = amrex::second() - strttime;

#ifdef AMREX_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
        amrex::Print() << "AmrActiveParticleContainer::Timestamp: lev: " << lev << " time: " << stoptime << '\n';
#ifdef AMREX_LAZY
        });
#endif
    }
}

  // <warp>
void
ActiveParticleContainer::getDrag(FArrayBox& rhs, const FArrayBox& vel, const FArrayBox& rho, Real nu_m, int nGrow, int level) {

    
    int num_levels = rhs.size();
    int finest_level = num_levels - 1;

    // each level deposits it's own particles => but getForce is called for each level already?
    //    const int ng; //= rhs[0]->nGrow();
    int ng;
    //    const int ng = rhs.nGrow();

    // level should be passed in
    //    for (int lev = 0; lev < num_levels; ++lev) {       

    ng = nGrow; //rhs.nGrow(); FIX

    //?        rhs.setVal(0.0, ng);

    //    std::cout<< " ...checkpoint 1, ng: " << ng << "\n";

        const auto& gm = m_gdb->Geom(level);
        const auto& ba = m_gdb->ParticleBoxArray(level);
        const Real* dx  = gm.CellSize();
        const Real* plo = gm.ProbLo();
        BoxArray nba = ba;
        nba.surroundingNodes();

	//        std::cout<< " ...checkpoint 2\n";
    
	//	for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        // iterates over particles in current level
        for (ParIterType pti(*this, level); pti.isValid(); ++pti) {

	    //const Box& box = pti.validbox();
  	    const Box& box = nba[pti];

	    //            std::cout<< " ...checkpoint 3, level:" << level << "\n";
            
	    //            auto& wp = pti.GetAttribs(PIdx::w);
            const auto& particles = pti.GetArrayOfStructs();
            const int nstride = particles.dataShape().first;
            const int np  = pti.numParticles();

	    //            std::cout<< " ...checkpoint 4, level:" << level << "\n";

	    
	    drag_cic( &ng, &nstride, &np, 
                      box.loVect(), box.hiVect(), 
                      &nu_m, dx, plo,
                      vel.dataPtr(), rho.dataPtr(), rhs.dataPtr(),
                      particles.data() );
	    

	    //            std::cout<< " ...checkpoint 5, level:" << level << "\n";

        }




	// not necessary with outer NS level loop?
    /*

    // now we average down fine to crse
    std::unique_ptr<MultiFab> crse;
    for (int lev = finest_level - 1; lev >= 0; --lev) {
      //        const BoxArray& fine_BA = rhs[lev+1]->boxArray();
        const BoxArray& fine_BA = rhs[level+1].boxArray();
	//        const DistributionMapping& fine_dm = rhs[lev+1]->DistributionMap();
        const DistributionMapping& fine_dm = rhs[level+1].DistributionMap();
      //        const BoxArray& fine_BA = rhs[lev+1].boxArray();
      //        const DistributionMapping& fine_dm = rhs[lev+1].DistributionMap();
        BoxArray coarsened_fine_BA = fine_BA;
        coarsened_fine_BA.coarsen(m_gdb->refRatio(lev));
        
        MultiFab coarsened_fine_data(coarsened_fine_BA, fine_dm, 1, 0);
        coarsened_fine_data.setVal(0.0);
        
        IntVect ratio(D_DECL(2, 2, 2));  // FIXME
        
        for (MFIter mfi(coarsened_fine_data); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            const Box& crse_box = coarsened_fine_data[mfi].box();
	    //            const Box& fine_box = (*rhs[lev+1])[mfi].box();
            const Box& fine_box = (*rhs[level+1])[mfi].box();
            sum_fine_to_crse_nodal(bx.loVect(), bx.hiVect(), ratio.getVect(),
                                   coarsened_fine_data[mfi].dataPtr(), crse_box.loVect(), crse_box.hiVect(),
                                   (*rhs[lev+1])[mfi].dataPtr(), fine_box.loVect(), fine_box.hiVect());
        }
        
	//        rhs[lev]->copy(coarsened_fine_data, m_gdb->Geom(lev).periodicity(), FabArrayBase::ADD);
        rhs[level]->copy(coarsened_fine_data, m_gdb->Geom(level).periodicity(), FabArrayBase::ADD);
	//        rhs[lev].copy(coarsened_fine_data, m_gdb->Geom(lev).periodicity(), FabArrayBase::ADD);
    }
    
    //    for (int lev = 0; lev < num_levels; ++lev) {
      //        rhs[lev]->mult(-1.0/PhysConst::ep0, ng);
      //      rhs[lev]->mult(-1.0, ng);  // not sure about this sign...
    //    }

    */

}




void
ActiveParticleContainer::getTemp(FArrayBox& rhs, const FArrayBox& vel, const FArrayBox& temp, Real nu_m, int nGrow, int level) {

    int num_levels = rhs.size();
    int finest_level = num_levels - 1;
    int ng = nGrow;

    const auto& gm = m_gdb->Geom(level);
    const auto& ba = m_gdb->ParticleBoxArray(level);
    const Real* dx  = gm.CellSize();
    const Real* plo = gm.ProbLo();
    BoxArray nba = ba;
    nba.surroundingNodes();

    for (ParIterType pti(*this, level); pti.isValid(); ++pti) {

       const Box& box = nba[pti];
       const auto& particles = pti.GetArrayOfStructs();
       const int nstride = particles.dataShape().first;
       const int np  = pti.numParticles();

       temp_cic( &ng, &nstride, &np, 
                 box.loVect(), box.hiVect(), 
                 &nu_m, dx, plo,
                 vel.dataPtr(), temp.dataPtr(), rhs.dataPtr(),
                 particles.data() );

    }
}



}
