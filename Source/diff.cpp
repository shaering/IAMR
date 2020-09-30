115a116,120
> 	// From PROB_*.F90
> 	//	subroutine FORT_INITDATA(level,time,lo,hi,nscal, &
> 	//				 vel,scal,DIMS(state),press,DIMS(press), &
> 	//				 dx,xlo,xhi) &
> 	//	  bind(C, name="FORT_INITDATA")
233a239
>     if(verbose) std::cout << " *** initParticleData okay\n";
308,310c314,318
<     MultiFab::Copy(*viscnp1_cc, *viscn_cc, 0, 0, 1, viscn_cc->nGrow());
<     MultiFab::Copy(*diffnp1_cc, *diffn_cc, 0, 0, num_diff, diffn_cc->nGrow());
<     
---
>     for (int d=0; d<AMREX_SPACEDIM; ++d) {
>       MultiFab::Copy(*viscnp1[d], *viscn[d], 0, 0, 1, viscn[d]->nGrow());
>       MultiFab::Copy(*diffnp1[d], *diffn[d], 0, 0, num_diff, diffn[d]->nGrow());
>     }
> 
322c330
<     Real dt_test = predict_velocity(dt);
---
>     Real dt_test = predict_velocity(dt,level);
421a430
> 
425,426d433
<       // original
<       //        theNSPC()->AdvectWithUmac(u_mac, level, dt);
428c435,436
<       // modified
---
>       //      theNSPC()->basicAdvectWithUmac(u_mac, level, dt);
>       
432a441
> 
434a444,448
> 
>       //FillPatchIterator S_fpi(*this,visc_terms,nGrow,prev_time,State_Type,Density,NUM_SCALARS);
>       //MultiFab& phi=S_fpi.get_mf();
>       
>       //      theNSPC()->myAdvectWithUmac(u_mac, level, dt, rho, temp, visc_coef[0]); // PARTICLES (particles) ADVANCED HERE
437c451
< 	
---
>       
439a454
> 
461c476
< NavierStokes::predict_velocity (Real  dt)
---
> NavierStokes::predict_velocity (Real  dt, int level)
537a553,596
> 
>     /*** Trying to bump drag calcs out of next MFIter loop ***/
>                                                                      
>     // Calc drag force
>     int num_levels = Umf.size();
>     int Ncomp  = 1;
> 
>     /*
>     Vector<DistributionMapping> dm(num_levels);
>     Vector< std::unique_ptr<MultiFab> > rhs(num_levels);
> 
>     BoxArray nba = grids;
>     nba.surroundingNodes();
>     for (int lev = 0; lev < num_levels; ++lev) {
>       //        BoxArray nba = grids[lev];
>       //        Box nba = grids[lev];
>       //        const Geometry& geom = m_gdb->Geom(lev);
>       //        const auto& ba = m_gdb->ParticleBoxArray(lev);
>       //        BoxArray nba = ba;
>       //        dm[lev].define(grids[lev]);
>         rhs[lev].reset(new MultiFab(nba, dm[lev], Ncomp, 1));
>         rhs[lev]->setVal(0.0); // zero
>      }
>      //   rhs.setVal(0.0,ng);
> */
> 
>    // From an already defined MultiFab
>     //   const BoxArray& ba = Umf.boxArray();
>    const BoxArray& ba = Smf.boxArray();
>    const DistributionMapping& dm = Umf.DistributionMap();
>    int ncomp = Umf.nComp();
>    int scomp = Smf.nComp();
>    int ngrow = Umf.nGrow();
>    MultiFab rhs(ba,dm,ncomp,ngrow);  // new MF with the same ncomp and ngrow
>    rhs.setVal(0.0);
>    //   MultiFab& rhs=temp.get_mf();
>    //   Vector< std::unique_ptr<MultiFab> > rhs = &temp;
> 
>     //#ifdef AMREX_PARTICLES
>     //    theNSPC()->get_Drag(rhs,Smf,Smf,Smf,Smf,visc_coef[0]);
>     //#endif
> 
> 
> 
541a601,602
> 
>       //        FArrayBox tforces(ba[0],ncomp+scomp);
545c606,607
<         for (MFIter U_mfi(Umf,true); U_mfi.isValid(); ++U_mfi)
---
> 	//AMReX provides an iterator, MFIter for looping over the FArrayBoxes in MultiFabs
>         for (MFIter U_mfi(Umf,true); U_mfi.isValid(); ++U_mfi) // this iters over FArrayBoxes in the MultiFab
548a611,615
>             FArrayBox& Sfab = Smf[U_mfi];
> 	    FArrayBox& rfab = rhs[U_mfi];
> 	    // FArrayBox& rfab = (*rhs[U_mfi]);
> 
> 	    //            FArrayBox& tforces = Umf[U_mfi];
551c618,653
<                 Print() << "---\nA - Predict velocity:\n Calling getForce...\n";
---
>                 Print() << "--------------------- \n "
>                         << "A - Predict velocity: \n" 
>                         << "--------------------- \n "
> 		        << "("<<"1"<<")...\n";
>             }
> 
> 
> #ifdef AMREX_PARTICLES
> 	    theNSPC()->getDrag(rfab,Ufab,Sfab,visc_coef[0],1,level);
> 	    std::cout << " *** get forcing okay\n";
> 
>             /* // may be an issue with uncommented sumbndry with different levels
> 	    if (level > 0) {
> 	      IntVect ghostVect(tmp_src_width*IntVect::TheUnitVector());
>               tmp_src_ptr->SumBoundary(0, tmp_src_ptr->nComp(), ghostVect, Geom(lev).periodicity());
> 	    } 
>             else {
>               tmp_src_ptr->SumBoundary(Geom(lev).periodicity());
> 	    }
> 	    */
> 
>             // seems to work
> 	    rhs.SumBoundary(0, ncomp, IntVect(1), Geom().periodicity()); 
> #endif
> 
> 
>             // Compute the total forcing (3:ngrow, 4:scomp, 5: ncomp) 4=0 / 5=3 => vel only 
>             getForce(tforces,bx,1,Xvel,BL_SPACEDIM,prev_time,Ufab,Sfab,rfab,0,level);  // EDGE VEL
> 
> 	    // (int scomp, int ncomp, IntVect const& nghost, const Periodicity& period), sensitive to cell center and ngrow of MF
> 	    // if only one arg (0, n_comp, IntVect(0), period);
>             //rhs.SumBoundary(Geom().periodicity()); 
> 	    //rhs.FillBoundary(Geom().periodicity());
> 
>             if (getForceVerbose) {
>                 Print() << "                                                    ... and done\n";
553d654
<             getForce(tforces,bx,1,Xvel,BL_SPACEDIM,prev_time,Ufab,Smf[U_mfi],0);
555,557d655
<             //
<             // Compute the total forcing.
<             //
558a657
>             //godunov->Sum_tf_gp_visc(rfab,0,visc_terms[U_mfi],0,Gp[U_mfi],0,rho_ptime[U_mfi],0);
649a749,757
>    // From an already defined MultiFab
>    const BoxArray& ba = Smf.boxArray();
>    const DistributionMapping& dm = Smf.DistributionMap();
>    int ncomp = Smf.nComp();
>    int ngrow = Smf.nGrow();
>    MultiFab rhs(ba,dm,ncomp,ngrow);  // new MF with the same ncomp and ngrow
>    rhs.setVal(0.0);
> 
> 
706a815,824
> 
> 
> #ifdef AMREX_PARTICLES
> 	    theNSPC()->getTemp(rhs[S_mfi],Umf[S_mfi],Smf[S_mfi],visc_coef[0],ngrow,level);
> 	    std::cout << " *** get forcing (temp) okay\n";
>             //rhs.SumBoundary(Geom().periodicity()); 
> 	    rhs.SumBoundary(0, ncomp, IntVect(1), Geom().periodicity()); 
> #endif
> 
> 
708,709c826,829
<                     Print() << "---" << '\n' << "C - scalar advection:" << '\n'
<                             << " Calling getForce..." << '\n';
---
>                     Print() << "---------------------\n" 
>                             << "C - Scalar advection:\n"
>                             << "---------------------\n"
>                             << " (1)..." << '\n';
711c831,839
<                 getForce(tforces,bx,nGrowF,fscalar,num_scalars,prev_time,Umf[S_mfi],Smf[S_mfi],0);
---
>                 // Compute the total forcing (3:ngrow, 4:scomp, 5: ncomp) 4=3 / 5=n(1) => scalars(density only) 
> 		//                getForce(tforces,bx,nGrowF,fscalar,num_scalars,prev_time,Umf[S_mfi],Smf[S_mfi],rhs[S_mfi],0,level); // arbitrary source terms?
>                 getForce(tforces,bx,1,fscalar,num_scalars,prev_time,Umf[S_mfi],Smf[S_mfi],rhs[S_mfi],0,level); // arbitrary source terms?
>   	        rhs.SumBoundary(Geom().periodicity()); 
>   	        //rhs.FillBoundary(Geom().periodicity());
> 
>                 if (getForceVerbose) {
>                    Print() << "                        ... and done\n";
> 	        }
735c863
< 
---
> 		// just add a hack with cfluxes zero for particle_vel scalar?
888a1017
> 	// forcing/source in delta_rhs
966d1094
< 	FluxBoxes fb_viscn, fb_viscnp1;
968a1097
> 	FluxBoxes fb_viscn, fb_viscnp1;
977c1106
< 	
---
> 
981c1110
<                                     delta_rhs,loc_viscn,viscn_cc,loc_viscnp1,viscnp1_cc);
---
>                                     delta_rhs,loc_viscn,loc_viscnp1);
1104a1234
>     Real regrad = 0.0;
1118a1249
>         regrad = std::max(regrad,ns_level.MaxVal("re_grad",time));
1132a1264
>     Print().SetPrecision(12) << "TIME= " << time << " RE_G= " << regrad << '\n';
1751c1883
<     // Delete Ucorr; we're done with it.
---
>     // fixme? clear Ucorr here? think we're done with it
1753,1755d1884
<     for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
<       delete Ucorr[idim];
< 
2320,2325c2449
< 	auto whichTime = which_time(State_Type,time);
< 	BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
< 	
< 	auto viscosityCC = (whichTime == AmrOldTime ? viscn_cc : viscnp1_cc);
< 	
<         diffusion->getTensorViscTerms(visc_terms,time,viscosity,viscosityCC,0);
---
>         diffusion->getTensorViscTerms(visc_terms,time,viscosity,0);
2384,2391c2508,2512
< // Functions calcViscosity/Diffusivity and getViscosity/Diffusivity are  
< // for calculating variable viscosity and diffusivity. Here we default to
< // constant visc/diff and set the variable viscosity and diffusivity arrays
< // to the values in visc_coef and diff_coef.
< // For variable viscosity/diffusivity, (per MSD) calcViscosity/Diffusivity
< // should compute the transport coefficients at cell centers (or cell centroids
< // for EB) and getViscosity/Diffusivity should interpolate those to faces (or
< // face-centroids for EB).
---
> // Functions for calculating the variable viscosity and diffusivity.
> // These default to setting the variable viscosity and diffusivity arrays
> // to the values in visc_coef and diff_coef.  These functions would
> // need to be replaced in any class derived from NavierStokes that
> // wants variable coefficients.
2406,2407c2527,2550
<             auto visc = (whichTime == AmrOldTime ? viscn_cc : viscnp1_cc);
< 	    visc->setVal(visc_coef[Xvel], 0, visc->nComp(), visc->nGrow());
---
>             auto visc = (whichTime == AmrOldTime ? viscn : viscnp1);
>             for (int dir=0; dir<AMREX_SPACEDIM; dir++) {
>                 visc[dir]->setVal(visc_coef[Xvel], 0, visc[dir]->nComp(), visc[dir]->nGrow());
>             }
> 
>             if (do_LES){
> 
>                FluxBoxes mu_LES(this,1,0);
>                MultiFab** mu_LES_mf = mu_LES.get();
>                for (int dir=0; dir<AMREX_SPACEDIM; dir++) {
>                 mu_LES_mf[dir]->setVal(0., 0, mu_LES_mf[dir]->nComp(), mu_LES_mf[dir]->nGrow());
>             }
> 
>                NavierStokesBase::calc_mut_LES(mu_LES_mf,time);
> 
>              for (int dir=0; dir<AMREX_SPACEDIM; dir++) {
>                 MultiFab::Add(*visc[dir], *mu_LES_mf[dir], 0, 0, 1, 0);
> 
>             }
> 
> 
>             }
> 
> 
2410c2553
< 	{
---
>         {
2432c2575
<     MultiFab* diff = (whichTime == AmrOldTime ? diffn_cc : diffnp1_cc);
---
>     MultiFab** diff = (whichTime == AmrOldTime ? diffn : diffnp1);
2441c2584,2586
< 	      diff->setVal(visc_coef[comp], diff_comp, 1, diff->nGrow());
---
>                 for (int dir=0; dir<AMREX_SPACEDIM; dir++) {
>                     diff[dir]->setVal(visc_coef[comp], diff_comp, 1, diff[dir]->nGrow());
>                 }
2455,2483c2600,2609
<     // //
<     // // Select time level to work with (N or N+1)
<     // //
<     // const TimeLevel whichTime = which_time(State_Type,time);
<     // BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
< 
<     // MultiFab *visc = (whichTime == AmrOldTime ? viscn_cc : viscnp1_cc);
<   
<     // For non-const viscosity, this where the interp from cell-center/centroid
<     // to faces would take place. 
<     // But here we simply do constant viscosity.
< 
<     for (int dir=0; dir<AMREX_SPACEDIM; dir++) {
<       viscosity[dir]->setVal(visc_coef[Xvel], 0, viscosity[dir]->nComp(), viscosity[dir]->nGrow());
<     }
< 	  
<     if (do_LES)
<     {
<       FluxBoxes mu_LES(this,1,0);
<       MultiFab** mu_LES_mf = mu_LES.get();
<       for (int dir=0; dir<AMREX_SPACEDIM; dir++) {
< 	mu_LES_mf[dir]->setVal(0., 0, mu_LES_mf[dir]->nComp(), mu_LES_mf[dir]->nGrow());
<       }
<       
<       NavierStokesBase::calc_mut_LES(mu_LES_mf,time);
<       
<       for (int dir=0; dir<AMREX_SPACEDIM; dir++) {
< 	MultiFab::Add(*viscosity[dir], *mu_LES_mf[dir], 0, 0, 1, 0);
<       }
---
>     //
>     // Select time level to work with (N or N+1)
>     //
>     const TimeLevel whichTime = which_time(State_Type,time);
>     BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
> 
>     MultiFab **visc = (whichTime == AmrOldTime ? viscn : viscnp1);
>     for (int dir = 0; dir < BL_SPACEDIM; dir++)
>     {
>         MultiFab::Copy(*viscosity[dir],*visc[dir],0,0,1,0);
2499,2509c2625,2629
<     // //
<     // // Select time level to work with (N or N+1)
<     // //
<     // const TimeLevel whichTime = which_time(State_Type,time);
<     // BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
<     
<     // MultiFab *diff = (whichTime == AmrOldTime ? diffn_cc : diffnp1_cc);
< 
<     // For non-const diffusivity, this where the interp from cell-center/centroid
<     // to faces would take place. 
<     // But here we simply do constant diffusivity.
---
>     //
>     // Select time level to work with (N or N+1)
>     //
>     const TimeLevel whichTime = which_time(State_Type,time);
>     BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
2510a2631
>     MultiFab **diff = (whichTime == AmrOldTime ? diffn : diffnp1);
2513c2634
<       diffusivity[dir]->setVal(visc_coef[diff_comp], dst_comp, ncomp, diffusivity[dir]->nGrow());
---
>         MultiFab::Copy(*diffusivity[dir],*diff[dir],diff_comp,dst_comp,ncomp,0);
