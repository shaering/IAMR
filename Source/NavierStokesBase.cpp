
#include <AMReX_ParmParse.H>
#include <AMReX_TagBox.H>
#include <AMReX_Utility.H>
#include <AMReX_PhysBCFunct.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBInterpolater.H>
#include <AMReX_EBFArrayBox.H>
#include <iamr_mol.H>
#endif

#include <iamr_godunov.H>

#include <NavierStokes.H>// okay to have this here?
#include <NavierStokesBase.H>
#include <NAVIERSTOKES_F.H>
#include <AMReX_filcc_f.H>
#include <NSB_K.H>

#include <PROB_NS_F.H>

//fixme, for writesingle level plotfile
#include<AMReX_PlotFileUtil.H>

using namespace amrex;

struct DummyFill           // Set 0.0 on EXT_DIR, nothing otherwise.
{
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
       const int* domlo = geom.Domain().loVect();
       const int* domhi = geom.Domain().hiVect();
       for (int n = 0; n < numcomp; n++ ) {
          const int* bc = bcr[n].data();
          for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
             if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {
                dest(iv, dcomp+n) = 0.0;
             }
             if ((bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and (iv[idir] > domhi[idir])) {
                dest(iv, dcomp+n) = 0.0;
             }
          }
       }
    }
};

ErrorList   NavierStokesBase::err_list;
BCRec       NavierStokesBase::phys_bc;
Projection* NavierStokesBase::projector     = 0;
MacProj*    NavierStokesBase::mac_projector = 0;

Real NavierStokesBase::init_shrink        = 1.0;
int  NavierStokesBase::init_iter          = 2;
int  NavierStokesBase::init_vel_iter      = 1;
Real NavierStokesBase::cfl                = 0.8;
Real NavierStokesBase::change_max         = 1.1;
Real NavierStokesBase::init_dt            = -1.0;
Real NavierStokesBase::fixed_dt           = -1.0;
bool NavierStokesBase::stop_when_steady   = false;
Real NavierStokesBase::steady_tol         = 1.0e-10;
int  NavierStokesBase::initial_iter       = false;
int  NavierStokesBase::initial_step       = false;
Real NavierStokesBase::dt_cutoff          = 0.0;
int  NavierStokesBase::sum_interval       = -1;
int  NavierStokesBase::turb_interval      = -1;
int  NavierStokesBase::jet_interval       = -1;
int  NavierStokesBase::jet_interval_split = 2;

int  NavierStokesBase::radius_grow = 1;
int  NavierStokesBase::verbose     = 0;
Real NavierStokesBase::gravity     = 0.0;
int  NavierStokesBase::NUM_SCALARS = 0;
int  NavierStokesBase::NUM_STATE   = 0;

Vector<AdvectionForm> NavierStokesBase::advectionType;
Vector<DiffusionForm> NavierStokesBase::diffusionType;

Vector<int>  NavierStokesBase::is_diffusive;
Vector<Real> NavierStokesBase::visc_coef;
Real        NavierStokesBase::visc_tol           = 1.0e-10;
Real        NavierStokesBase::visc_abs_tol       = 1.0e-10;
Real        NavierStokesBase::be_cn_theta        = 0.5;

int         NavierStokesBase::Tracer                    = -1;
int         NavierStokesBase::Tracer2                   = -1;
int         NavierStokesBase::Temp                      = -1;
int         NavierStokesBase::do_trac2                  = 0;
int         NavierStokesBase::do_temp                   = 0;
int         NavierStokesBase::do_cons_trac              = 0;
int         NavierStokesBase::do_cons_trac2             = 0;
int         NavierStokesBase::do_sync_proj              = 1;
int         NavierStokesBase::do_reflux                 = 1;
int         NavierStokesBase::modify_reflux_normal_vel  = 0;
int         NavierStokesBase::do_mac_proj               = 1;
int         NavierStokesBase::do_refine_outflow         = 0;
int         NavierStokesBase::do_derefine_outflow       = 1;
int         NavierStokesBase::Nbuf_outflow              = 1;
int         NavierStokesBase::do_denminmax              = 0;
int         NavierStokesBase::do_scalminmax             = 0;
int         NavierStokesBase::do_density_ref            = 0;
int         NavierStokesBase::do_tracer_ref             = 0;
int         NavierStokesBase::do_tracer2_ref            = 0;
int         NavierStokesBase::do_vorticity_ref          = 0;
int         NavierStokesBase::do_temp_ref               = 0;
int         NavierStokesBase::do_scalar_update_in_order = 0;
Vector<int>  NavierStokesBase::scalarUpdateOrder;
int         NavierStokesBase::getForceVerbose           = 0;
int         NavierStokesBase::do_LES                    = 0;
int         NavierStokesBase::getLESVerbose             = 0;
std::string NavierStokesBase::LES_model                 = "Smagorinsky";
Real        NavierStokesBase::smago_Cs_cst              = 0.18;
Real        NavierStokesBase::sigma_Cs_cst              = 1.5;

int         NavierStokesBase::particle_extra_reals      = 9;
int         NavierStokesBase::particle_extra_ints       = 0;


int  NavierStokesBase::Dpdt_Type = -1;

int  NavierStokesBase::additional_state_types_initialized = 0;
int  NavierStokesBase::Divu_Type                          = -1;
int  NavierStokesBase::Dsdt_Type                          = -1;
int  NavierStokesBase::num_state_type                     = 2;
int  NavierStokesBase::have_divu                          = 0;
int  NavierStokesBase::have_dsdt                          = 0;
Real NavierStokesBase::divu_relax_factor                  = 0.0;
int  NavierStokesBase::S_in_vel_diffusion                 = 0;
int  NavierStokesBase::do_init_vort_proj                  = 0;
int  NavierStokesBase::do_init_proj                       = 1;

int  NavierStokesBase::do_running_statistics  = 0;
Real NavierStokesBase::volWgtSum_sub_origin_x = 0;
Real NavierStokesBase::volWgtSum_sub_origin_y = 0;
Real NavierStokesBase::volWgtSum_sub_origin_z = 0;
Real NavierStokesBase::volWgtSum_sub_Rcyl     = -1;
Real NavierStokesBase::volWgtSum_sub_dx       = -1;
Real NavierStokesBase::volWgtSum_sub_dy       = -1;
Real NavierStokesBase::volWgtSum_sub_dz       = -1;

int  NavierStokesBase::do_mom_diff            = 0;
int  NavierStokesBase::predict_mom_together   = 0;
bool NavierStokesBase::def_harm_avg_cen2edge  = false;

bool NavierStokesBase::godunov_use_ppm = false;
bool NavierStokesBase::godunov_use_forces_in_trans = false;

#ifdef AMREX_USE_EB
int          NavierStokesBase::refine_cutcells     = 1;
bool         NavierStokesBase::eb_initialized      = false;
bool         NavierStokesBase::no_eb_in_domain     = true;
bool         NavierStokesBase::body_state_set      = false;
std::vector<Real> NavierStokesBase::body_state;
#endif

namespace
{
    bool initialized = false;
    int  dump_plane  = -1;
    std::string dump_plane_name("SLABS/vel-");
    bool benchmarking = false;
}

#ifdef AMREX_PARTICLES
namespace
{
    //
    // Name of subdirectory in chk???? holding checkpointed particles.
    //
    const std::string the_ns_particle_file_name("Particles");
    //
    // There's really only one of these.
    //
    AmrActiveParticleContainer* NSPC = 0;

    std::string      timestamp_dir                   ("Timestamps");
    std::vector<int> timestamp_indices;
    std::string      particle_init_file;
    std::string      particle_restart_file;
    std::string      particle_output_file;
    bool             restart_from_nonparticle_chkfile = false;
    int              pverbose                         = 2;
}

AmrActiveParticleContainer* NavierStokesBase::theNSPC () { return NSPC; }
#endif

int NavierStokesBase::DoTrac2() {return NavierStokesBase::do_trac2;}
//
BL_FORT_PROC_DECL(BL_NS_DOTRAC2,bl_ns_dotrac2)(int* dotrac2)
{
    *dotrac2 = NavierStokesBase::DoTrac2();
}

NavierStokesBase::NavierStokesBase ()
{
    rho_qtime    = 0;
    rho_tqtime   = 0;
    sync_reg     = 0;
    advflux_reg  = 0;
    viscflux_reg = 0;
    u_mac        = 0;
    aofs         = 0;
    diffusion    = 0;


    if (!additional_state_types_initialized)
        init_additional_state_types();
}

NavierStokesBase::NavierStokesBase (Amr&            papa,
				    int             lev,
				    const Geometry& level_geom,
				    const BoxArray& bl,
                                    const DistributionMapping& dm,
				    Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,dm,time)
{
    //
    // 8/2020 - Neither MLMG nor IAMR fully support rz 
    //
    if ( level_geom.IsRZ() )
      amrex::Abort("RZ geometry is not currently supported");
  
    if(!additional_state_types_initialized) {
        init_additional_state_types();
    }

    const BoxArray& P_grids = state[Press_Type].boxArray();
    //
    // Alloc old_time pressure.
    //
    state[Press_Type].allocOldData();
    //
    // Alloc space for density and temporary pressure variables.
    //
    if (level > 0)
    {
        rho_avg.define(grids,dmap,1,1,MFInfo(),Factory());
        p_avg.define(P_grids,dmap,1,0,MFInfo(),Factory());
    }

    //
    // rho_half is passed into level_project to be used as sigma in the MLMG
    // solve, but MLMG doesn't copy any ghost cells, it fills what it needs itself.
    // does rho_half still need any ghost cells?
    rho_half.define (grids,dmap,1,1,MFInfo(),Factory());
    rho_ptime.define(grids,dmap,1,1,MFInfo(),Factory());
    rho_ctime.define(grids,dmap,1,1,MFInfo(),Factory());
    rho_qtime  = 0;
    rho_tqtime = 0;
    //
    // Build metric coefficients for RZ calculations.
    // Build volume and areas.
    //
    buildMetrics();


#ifdef AMREX_USE_EB

    init_eb(level_geom, bl, dm);

    //fixme? not 100% sure this is the right place
    gradp.reset(new MultiFab(grids,dmap,BL_SPACEDIM,1, MFInfo(), Factory()));
    gradp->setVal(0.);

    //FIXME --- this fn is really similar to restart()... work on that later
#endif

    //
    // Set up reflux registers.
    //
    sync_reg = 0;
    if (level > 0 && do_sync_proj)
    {
        sync_reg = new SyncRegister(grids,dmap,crse_ratio);
    }
    advflux_reg  = 0;
    viscflux_reg = 0;
    if (level > 0 && do_reflux)
    {
        advflux_reg  = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
        viscflux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
    }
    //
    // Initialize work multifabs.
    //
    u_mac   = 0;
    aofs    = 0;

    //
    // Set up the level projector.
    //
    if (projector == 0)
    {
        projector = new Projection(parent,&phys_bc,do_sync_proj,
                                   parent->finestLevel(),radius_grow);
    }
    projector->install_level(level,this,&radius);

    //
    // Set up diffusion.
    //
    diffusion = new Diffusion(parent,this,
                              (level > 0) ? getLevel(level-1).diffusion : 0,
                              NUM_STATE,viscflux_reg,is_diffusive,visc_coef);
    //
    // Allocate the storage for variable viscosity and diffusivity
    //
    diffn_cc = new MultiFab(grids, dmap, NUM_STATE-Density-1, 1, MFInfo(), Factory());
    diffnp1_cc = new MultiFab(grids, dmap, NUM_STATE-Density-1, 1, MFInfo(), Factory());
    viscn_cc = new MultiFab(grids, dmap, 1, 1, MFInfo(), Factory());
    viscnp1_cc = new MultiFab(grids, dmap, 1, 1, MFInfo(), Factory());

    //
    // Set up the mac projector.
    //
    if (mac_projector == 0)
    {
        mac_projector = new MacProj(parent,parent->finestLevel(),
                                    &phys_bc,radius_grow);
    }
    mac_projector->install_level(level,this);
}

NavierStokesBase::~NavierStokesBase ()
{
    delete rho_qtime;
    delete rho_tqtime;
    delete sync_reg;
    delete advflux_reg;
    delete viscflux_reg;
    delete [] u_mac;

    if (mac_projector != 0)
        mac_projector->cleanup(level);
    //
    // Remove the arrays for variable viscosity and diffusivity
    // and delete the Diffusion object
    //
    delete viscn_cc;
    delete viscnp1_cc;
    delete diffn_cc;
    delete diffnp1_cc;

    delete diffusion;
}

void
NavierStokesBase::allocOldData ()
{
    bool init_pres = !(state[Press_Type].hasOldData());
    AmrLevel::allocOldData();
    if (init_pres)
        initOldPress();
}

void
NavierStokesBase::variableCleanUp ()
{
    desc_lst.clear();
    derive_lst.clear();

    err_list.clear();

    delete projector;
    projector = 0;

    delete mac_projector;
    mac_projector = 0;

#ifdef AMREX_PARTICLES
    delete NSPC;
    NSPC = 0;
#endif
}

void
NavierStokesBase::Initialize ()
{
    if (initialized) return;

    ParmParse pp("ns");

    pp.query("dump_plane",dump_plane);

    pp.query("benchmarking",benchmarking);

    pp.query("v",verbose);


    //
    // Get timestepping parameters.
    //
    pp.get("cfl",cfl);
    pp.query("init_iter",init_iter);
    pp.query("init_vel_iter",init_vel_iter);
    pp.query("init_shrink",init_shrink);
    pp.query("dt_cutoff",dt_cutoff);
    pp.query("change_max",change_max);
    pp.query("fixed_dt",fixed_dt);
    pp.query("init_dt", init_dt);
    pp.query("stop_when_steady",stop_when_steady);
    pp.query("steady_tol",steady_tol);
    pp.query("sum_interval",sum_interval);
    pp.query("turb_interval",turb_interval);
    pp.query("jet_interval",jet_interval);
    pp.query("jet_interval_split",jet_interval_split);
    pp.query("gravity",gravity);
    //
    // Get run options.
    //
    pp.query("do_temp",                  do_temp          );
    pp.query("do_trac2",                 do_trac2         );
    pp.query("do_cons_trac",             do_cons_trac     );
    pp.query("do_cons_trac2",            do_cons_trac2    );
    pp.query("do_sync_proj",             do_sync_proj     );
    pp.query("do_reflux",                do_reflux        );
    pp.query("modify_reflux_normal_vel", modify_reflux_normal_vel);
    pp.query("do_init_vort_proj",        do_init_vort_proj);
    pp.query("do_init_proj",             do_init_proj     );
    pp.query("do_mac_proj",              do_mac_proj      );
    pp.query("do_denminmax",             do_denminmax     );
    pp.query("do_scalminmax",            do_scalminmax    );
    pp.query("do_density_ref",           do_density_ref   );
    pp.query("do_tracer_ref",            do_tracer_ref    );
    pp.query("do_tracer2_ref",           do_tracer2_ref   );
    pp.query("do_vorticity_ref",         do_vorticity_ref );
    pp.query("do_temp_ref",              do_temp_ref      );

    pp.query("visc_tol",visc_tol);
    pp.query("visc_abs_tol",visc_abs_tol);

    if (modify_reflux_normal_vel)
        amrex::Abort("modify_reflux_normal_vel is no longer supported");

    pp.query("getForceVerbose",          getForceVerbose  );
    pp.query("do_LES",                   do_LES  );
    pp.query("getLESVerbose",            getLESVerbose  );
    pp.query("LES_model",                LES_model  );
    pp.query("smago_Cs_cst",             smago_Cs_cst  );
    pp.query("sigma_Cs_cst",             sigma_Cs_cst  );

#ifdef AMREX_USE_EB
    pp.query("refine_cutcells", refine_cutcells);
#endif

    pp.query("do_scalar_update_in_order",do_scalar_update_in_order );
    if (do_scalar_update_in_order) {
	    const int n_scalar_update_order_vals = pp.countval("scalar_update_order");
	    scalarUpdateOrder.resize(n_scalar_update_order_vals);
	    int got_scalar_update_order = pp.queryarr("scalar_update_order",scalarUpdateOrder,0,n_scalar_update_order_vals);
    }

    // Don't let init_shrink be greater than 1
    if (init_shrink > 1.0)
        amrex::Abort("NavierStokesBase::Initialize(): init_shrink cannot be greater than 1");

    pp.query("divu_relax_factor",divu_relax_factor);
    pp.query("S_in_vel_diffusion",S_in_vel_diffusion);
    if ( S_in_vel_diffusion ){
#ifdef AMREX_USE_EB
      // Currently, we should use the TensorOp to compute the divU terms in divtau.
      // The code is still present to use the source term S instead of a numerically
      // computed divu, however, divmusi terms isn't EB-aware.
      // Perhaps one day a comparision would be interesting.
      amrex::Abort("S_in_vel_diffusion not currently supported.\n");
#else
      amrex::Warning("WARNING: S_in_vel_diffusion is probably not what you want anymore. \nSuggested option is now to set S_in_vel_diffusion=0 to allow the tensor diffusion solver to compute divU.");
#endif
    }
    pp.query("be_cn_theta",be_cn_theta);
    if (be_cn_theta > 1.0 || be_cn_theta < .5)
        amrex::Abort("NavierStokesBase::Initialize(): Must have be_cn_theta <= 1.0 && >= .5");
    //
    // Set parameters dealing with how grids are treated at outflow boundaries.
    //
    pp.query("do_refine_outflow",do_refine_outflow);
    pp.query("do_derefine_outflow",do_derefine_outflow);
    if (do_derefine_outflow == 1 && do_refine_outflow == 1)
      amrex::Abort("NavierStokesBase::Initialize(): Cannot have both do_refine_outflow==1 and do_derefine_outflow==1");

    pp.query("Nbuf_outflow",Nbuf_outflow);
    BL_ASSERT(Nbuf_outflow >= 0);
    BL_ASSERT(!(Nbuf_outflow <= 0 && do_derefine_outflow == 1));

    //
    // Check whether we are doing running statistics.
    //
    pp.query("do_running_statistics",do_running_statistics);

    // If dx,dy,dz,Rcyl<0 (default) the volWgtSum is computed over the entire domain
    pp.query("volWgtSum_sub_origin_x",volWgtSum_sub_origin_x);
    pp.query("volWgtSum_sub_origin_y",volWgtSum_sub_origin_y);
    pp.query("volWgtSum_sub_origin_z",volWgtSum_sub_origin_z);
    pp.query("volWgtSum_sub_Rcyl",volWgtSum_sub_Rcyl);
    pp.query("volWgtSum_sub_dx",volWgtSum_sub_dx);
    pp.query("volWgtSum_sub_dy",volWgtSum_sub_dy);
    pp.query("volWgtSum_sub_dz",volWgtSum_sub_dz);

    // Are we going to do velocity or momentum update?
    pp.query("do_mom_diff",do_mom_diff);
    pp.query("predict_mom_together",predict_mom_together);

    if (do_mom_diff == 0 && predict_mom_together == 1)
    {
      amrex::Print() << "MAKES NO SENSE TO HAVE DO_MOM_DIFF=0 AND PREDICT_MOM_TOGETHER=1\n";
      exit(0);
    }

    pp.query("harm_avg_cen2edge", def_harm_avg_cen2edge);

#ifdef AMREX_PARTICLES
    read_particle_params ();
#endif

    //
    // Get godunov options
    //
    ParmParse pp2("godunov");
    pp2.query("use_ppm",             godunov_use_ppm);
    pp2.query("use_forces_in_trans", godunov_use_forces_in_trans);

    amrex::ExecOnFinalize(NavierStokesBase::Finalize);

    initialized = true;
}

// The following Initialize_specific is dedicated to read and set data
// only specific for IAMR, because it conflicts with PeleLM.
// PeleLM calls NavierStokesBase::Initialize() and its own PelelM::Initialize_specific ()
void
NavierStokesBase::Initialize_specific ()
{
    ParmParse pp("ns");

    Vector<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
    pp.getarr("lo_bc",lo_bc,0,BL_SPACEDIM);
    pp.getarr("hi_bc",hi_bc,0,BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        phys_bc.setLo(i,lo_bc[i]);
        phys_bc.setHi(i,hi_bc[i]);
    }

    read_geometry();
    //
    // Check phys_bc against possible periodic geometry
    // if periodic, must have internal BC marked.
    //
    if (DefaultGeometry().isAnyPeriodic())
    {
        //
        // Do idiot check.  Periodic means interior in those directions.
        //
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (DefaultGeometry().isPeriodic(dir))
            {
                if (lo_bc[dir] != Interior)
                {
                    std::cerr << "NavierStokesBase::variableSetUp:periodic in direction "
                              << dir
                              << " but low BC is not Interior\n";
                    amrex::Abort("NavierStokesBase::Initialize()");
                }
                if (hi_bc[dir] != Interior)
                {
                    std::cerr << "NavierStokesBase::variableSetUp:periodic in direction "
                              << dir
                              << " but high BC is not Interior\n";
                    amrex::Abort("NavierStokesBase::Initialize()");
                }
            }
        }
    }

    {
        //
        // Do idiot check.  If not periodic, should be no interior.
        //
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (!DefaultGeometry().isPeriodic(dir))
            {
              if (lo_bc[dir] == Interior)
              {
                  std::cerr << "NavierStokesBase::variableSetUp:Interior bc in direction "
                            << dir
                            << " but not defined as periodic\n";
                  amrex::Abort("NavierStokesBase::Initialize()");
              }
              if (hi_bc[dir] == Interior)
              {
                  std::cerr << "NavierStokesBase::variableSetUp:Interior bc in direction "
                            << dir
                            << " but not defined as periodic\n";
                  amrex::Abort("NavierStokesBase::Initialize()");
              }
            }
        }
    }

    //
    // Read viscous/diffusive parameters and array of viscous/diffusive coeffs.
    // NOTE: at this point, we dont know number of state variables
    //       so just read all values listed.
    //

    const int n_vel_visc_coef   = pp.countval("vel_visc_coef");
    const int n_temp_cond_coef  = pp.countval("temp_cond_coef");
    const int n_scal_diff_coefs = pp.countval("scal_diff_coefs");

    if (n_vel_visc_coef != 1)
        amrex::Abort("NavierStokesBase::Initialize(): Only one vel_visc_coef allowed");

    if (do_temp && n_temp_cond_coef != 1)
        amrex::Abort("NavierStokesBase::Initialize(): Only one temp_cond_coef allowed");

    int n_visc = BL_SPACEDIM + 1 + n_scal_diff_coefs;
    if (do_temp)
        n_visc++;
    visc_coef.resize(n_visc);
    is_diffusive.resize(n_visc);

    pp.get("vel_visc_coef",visc_coef[0]);
    for (int i = 1; i < BL_SPACEDIM; i++)
      visc_coef[i] = visc_coef[0];
    //
    // Here we set the coefficient for density, which does not diffuse.
    //
    visc_coef[Density] = -1;
    //
    // Set the coefficients for the scalars, but temperature.
    //
    Vector<Real> scal_diff_coefs(n_scal_diff_coefs);
    pp.getarr("scal_diff_coefs",scal_diff_coefs,0,n_scal_diff_coefs);

    int scalId = Density;

    // Will need to add more lines when more variables are added
    Tracer = Density+1;
    if (do_trac2)
	    Tracer2 = Density+2;

    for (int i = 0; i < n_scal_diff_coefs; i++)
    {
        visc_coef[++scalId] = scal_diff_coefs[i];
    }
    //
    // Set the coefficient for temperature.
    //
    if (do_temp)
    {
	    Temp = ++scalId;
	    pp.get("temp_cond_coef",visc_coef[Temp]);
    }
}


void
NavierStokesBase::Finalize ()
{
    initialized = false;
}

void
NavierStokesBase::read_geometry ()
{
#if (BL_SPACEDIM == 2)
    //
    // Must load coord here because Geometry hasn't read it in yet.
    //
    ParmParse pp("geometry");

    int coord;
    pp.get("coord_sys",coord);

    if ((Geometry::CoordType) coord == Geometry::RZ && phys_bc.lo(0) != Symmetry)
    {
        phys_bc.setLo(0,Symmetry);
	amrex::Print() << "\nWarning: Setting phys_bc at xlo to Symmetry\n\n";
    }
#endif
}

void
NavierStokesBase::advance_setup (Real time,
                                 Real dt,
	                         int  iteration,
                                 int  ncycle)
{
    BL_PROFILE("NavierStokesBase::advance_setup()");

    const int finest_level = parent->finestLevel();

#ifdef AMREX_USE_EB
    // incflo now uses: use_godunov ? 4 : 3;
    umac_n_grow = 4;
#else
    umac_n_grow = 1;
#endif

#ifdef AMREX_PARTICLES
    if (ncycle > umac_n_grow)
      umac_n_grow = ncycle;
#endif

    mac_projector->setup(level);
    //
    // Why are they defined here versus the constructor?
    //
    if (level < finest_level)
    {
        if (Vsync.empty())
            Vsync.define(grids,dmap,BL_SPACEDIM,1,MFInfo(),Factory());
        if (Ssync.empty())
	  Ssync.define(grids,dmap,NUM_STATE-BL_SPACEDIM,1,MFInfo(),Factory());
        Vsync.setVal(0);
        Ssync.setVal(0);
    }
    //
    // Set reflux registers to zero.
    //
    if (do_reflux && level < finest_level)
    {
        getAdvFluxReg(level+1).setVal(0);
        getViscFluxReg(level+1).setVal(0);
    }
    //
    // Alloc space for edge velocities (normal comp only).
    //
    if (u_mac == 0 || u_mac[0].nGrow() < umac_n_grow)
    {
	if (u_mac != 0) delete [] u_mac;

        u_mac = new MultiFab[BL_SPACEDIM];

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
	    const BoxArray& edgeba = getEdgeBoxArray(dir);
            u_mac[dir].define(edgeba,dmap,1,umac_n_grow,MFInfo(),Factory());
            u_mac[dir].setVal(1.e40);
        }
    }
    //
    // Alloc MultiFab to hold advective update terms.
    //
    BL_ASSERT(aofs == 0);
    aofs = new MultiFab(grids,dmap,NUM_STATE,0,MFInfo(),Factory());

    //
    // Set rho_avg.
    //
    if (!initial_step && level > 0 && iteration == 1)
        initRhoAvg(0.5/Real(ncycle));
    //
    // Set up state multifabs for the advance.
    //
    for (int k = 0; k < num_state_type; k++)
    {
	bool has_old_data = state[k].hasOldData();
        state[k].allocOldData();
	if (! has_old_data) state[k].oldData().setVal(0.0);
        state[k].swapTimeLevels(dt);
    }

    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
    {
        const Real new_press_time = .5 * (state[State_Type].prevTime() +
                                          state[State_Type].curTime());
        state[Press_Type].setNewTimeLevel(new_press_time);
    }

    make_rho_prev_time();

    // refRatio==4 is not currently supported
    //
    // If refRatio==4 to the next level coarser, and we're going to diffuse
    // scalars as SoverRho, we're going to need rho at 1/4 and 3/4 time there.
    // Make these things if need be.
    //
    // if (level > 0)
    // {
    //     bool needs_rho4 = false;

    //     if (parent->nCycle(level) == 4)
    //         for (int i = 0; i < NUM_STATE && !needs_rho4; ++i)
    //             needs_rho4 = (diffusionType[i] == Laplacian_SoverRho);

    //     if (needs_rho4)
    //     {
    //         NavierStokesBase& clevel = getLevel(level-1);
    //         const BoxArray& cgrids = clevel.boxArray();
    //         const DistributionMapping& cdmap = clevel.DistributionMap();
    //         const Real      ptime  = clevel.state[State_Type].prevTime();
    //         const Real      ctime  = clevel.state[State_Type].curTime();

    //         if (clevel.rho_qtime == 0)
    //         {
    //             const Real qtime = ptime + 0.25*(ctime-ptime);
    //             clevel.rho_qtime = new MultiFab(cgrids,cdmap,1,1);
    //             FillPatch(clevel,*(clevel.rho_qtime),1,qtime,State_Type,Density,1,0);
    //         }
    //         if (clevel.rho_tqtime == 0)
    //         {
    //             const Real tqtime = ptime + 0.75*(ctime-ptime);
    //             clevel.rho_tqtime = new MultiFab(cgrids,cdmap,1,1);
    // 		FillPatch(clevel, *(clevel.rho_tqtime), 1, tqtime, State_Type, Density, 1, 0);
    //         }
    //    }
    // }
}

//
// Clean up after the advance function.
//
void
NavierStokesBase::advance_cleanup (int iteration, int ncycle)
{
    delete aofs;
    aofs = 0;
}

void
NavierStokesBase::buildMetrics ()
{
    //
    // We "should" only need radius when we're RZ, but some 2-D code is written to
    // access it first and then "use" if if RZ.  It's easier to just always build
    // it for 2D than try to fix the underlying Fortran calls that take radius.
    //
#if (BL_SPACEDIM == 2)
    radius.resize(grids.size());

    const Real dxr = geom.CellSize()[0];

    for (int i = 0; i < grids.size(); i++)
    {
        const int ilo = grids[i].smallEnd(0)-radius_grow;
        const int ihi = grids[i].bigEnd(0)+radius_grow;
        const int len = ihi - ilo + 1;

        radius[i].resize(len);

        RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());

        const Real xlo = gridloc.lo(0) + (0.5 - radius_grow)*dxr;
        for (int j = 0; j < len; j++)
            radius[i][j] = xlo + j*dxr;
    }
#endif

    // volume and area are intentionally without EB knowledge
    volume.clear();
    volume.define(grids,dmap,1,GEOM_GROW);
    geom.GetVolume(volume);

    for (int dir = 0; dir < BL_SPACEDIM; ++dir)
    {
        area[dir].clear();
	area[dir].define(getEdgeBoxArray(dir),dmap,1,GEOM_GROW);
        geom.GetFaceArea(area[dir],dir);
    }

#ifdef AMREX_USE_EB
    // make sure dx == dy == dz
    const Real* dx = geom.CellSize();
    Print()<<"dx = "<<dx[0]<<" "<<dx[1]<<" "<<dx[2]<<" \n";
    for (int i = 1; i < BL_SPACEDIM; i++){
      if (std::abs(dx[i]-dx[i-1]) > 1.e-12*dx[0])
        amrex::Abort("EB requires dx == dy (== dz)\n");
    }

    const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());
    volfrac = &(ebfactory.getVolFrac());
    areafrac = ebfactory.getAreaFrac();


    //fixme? assume will need this part cribbed from CNS
    // level_mask.clear();
    // level_mask.define(grids,dmap,1,1);
    // level_mask.BuildMask(geom.Domain(), geom.periodicity(),
    //                      level_mask_covered,
    //                      level_mask_notcovered,
    //                      level_mask_physbnd,
    //                      level_mask_interior);

#endif
}

//
// Default dSdt is set to zero.
//
void
NavierStokesBase::calc_dsdt (Real      /*time*/,
                         Real      dt,
                         MultiFab& dsdt)
{
    if (have_divu && have_dsdt)
    {
      // Don't think we need this here, but then will have uninitialized ghost cells
      //dsdt.setVal(0);

        if (do_temp)
        {
            MultiFab& Divu_new = get_new_data(Divu_Type);
            MultiFab& Divu_old = get_old_data(Divu_Type);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
	    for (MFIter mfi(dsdt,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
	        const Box&  bx      = mfi.tilebox();
		auto const& div_new = Divu_new.array(mfi);
		auto const& div_old = Divu_old.array(mfi);
		auto const& dsdtarr = dsdt.array(mfi);

		amrex::ParallelFor(bx, [div_new, div_old, dsdtarr, dt]
	        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	        {
		  dsdtarr(i,j,k) = ( div_new(i,j,k) - div_old(i,j,k) )/ dt;
		});
            }
        }
	else
	{
	    dsdt.setVal(0);
	}
    }
}

void
NavierStokesBase::calcDpdt ()
{
    BL_ASSERT(state[Press_Type].descriptor()->timeType() == StateDescriptor::Point);

    MultiFab&  new_press   = get_new_data(Press_Type);
    MultiFab&  old_press   = get_old_data(Press_Type);
    MultiFab&  dpdt        = get_new_data(Dpdt_Type);
    const Real dt_for_dpdt = state[Press_Type].curTime()-state[Press_Type].prevTime();

    if (dt_for_dpdt != 0.0)
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(dpdt,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx        = mfi.tilebox();
            auto const& p_new    = new_press.array(mfi);
            auto const& p_old    = old_press.array(mfi);
            auto const& dpdt_arr = dpdt.array(mfi);
            Real   dt_inv = 1.0/dt_for_dpdt;
            amrex::ParallelFor(bx, [dpdt_arr, p_old, p_new, dt_inv]
            AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
               dpdt_arr(i,j,k) = dt_inv * (p_new(i,j,k) - p_old(i,j,k));
            });
        }
    }
    else
    {
        dpdt.setVal(0.0);
    }
}

void
NavierStokesBase::checkPoint (const std::string& dir,
			      std::ostream&      os,
			      VisMF::How         how,
			      bool               dump_old)
{
    AmrLevel::checkPoint(dir, os, how, dump_old);

#ifdef AMREX_PARTICLES
    if (level == 0)
    {
        if (NSPC != 0)
            NSPC->Checkpoint(dir,the_ns_particle_file_name);
    }
#endif

# ifdef AMREX_USE_EB
// Need to add gradp in the checkpoint
   std::string LevelDir, FullPath;
   LevelDirectoryNames(dir, LevelDir, FullPath);
   std::string gradp_mf_fullpath = FullPath + "/gradp";
   VisMF::Write(*gradp,gradp_mf_fullpath,how);
#endif
}

void
NavierStokesBase::computeInitialDt (int                   finest_level,
				    int                   sub_cycle,
				    Vector<int>&           n_cycle,
				    const Vector<IntVect>& ref_ratio,
				    Vector<Real>&          dt_level,
				    Real                  stop_time)
{
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;

    int i;

    Real dt_0    = 1.0e+100;
    int n_factor = 1;
    ///TODO/DEBUG: This will need to change for optimal subcycling.
    for (i = 0; i <= finest_level; i++)
    {
        dt_level[i] = getLevel(i).initialTimeStep();
        n_factor   *= n_cycle[i];
        dt_0        = std::min(dt_0,n_factor*dt_level[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    if (stop_time >= 0.0)
    {
        const Real eps      = 0.0001*dt_0;
        const Real cur_time = state[State_Type].curTime();
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        n_factor   *= n_cycle[i];
        dt_level[i] = dt_0/( (Real)n_factor );
    }
}

void
NavierStokesBase::computeNewDt (int                   finest_level,
				int                   sub_cycle,
				Vector<int>&           n_cycle,
				const Vector<IntVect>& ref_ratio,
				Vector<Real>&          dt_min,
				Vector<Real>&          dt_level,
				Real                  stop_time,
				int                   post_regrid_flag)
{
    //
    // We are at the start of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
        return;

    int i;

    Real dt_0     = 1.0e+100;
    int  n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        NavierStokesBase& adv_level = getLevel(i);
        dt_min[i] = std::min(dt_min[i],adv_level.estTimeStep());
    }

    if (fixed_dt <= 0.0)
    {
       if (post_regrid_flag == 1)
       {
          //
          // Limit dt's by pre-regrid dt
          //
          for (i = 0; i <= finest_level; i++)
          {
              dt_min[i] = std::min(dt_min[i],dt_level[i]);
          }
       }
       else
       {
          //
          // Limit dt's by change_max * old dt
          //
          for (i = 0; i <= finest_level; i++)
          {
	    if (verbose)
                 if (dt_min[i] > change_max*dt_level[i])
                 {
                     amrex::Print() << "NavierStokesBase::compute_new_dt : limiting dt at level "
				    << i << '\n';
                     amrex::Print() << " ... new dt computed: " << dt_min[i]
				    << '\n';
                     amrex::Print() << " ... but limiting to: "
				    << change_max * dt_level[i] << " = " << change_max
				    << " * " << dt_level[i] << '\n';
                 }
             dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
          }
       }
    }

    //
    // Find the minimum over all levels
    //
    for (i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_0      = std::min(dt_0,n_factor*dt_min[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps      = 0.0001*dt_0;
    const Real cur_time = state[State_Type].curTime();
    if (stop_time >= 0.0)
    {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    //
    // Set dt at each level of refinement
    //
    n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        n_factor   *= n_cycle[i];
        dt_level[i] = dt_0/( (Real)n_factor );
    }
}

void
NavierStokesBase::create_mac_rhs (MultiFab& rhs, int nGrow, Real time, Real dt)
{
    BL_PROFILE("NavierStokesBase::create_mac_rhs()");

    BL_ASSERT(rhs.nGrow()>=nGrow);
    BL_ASSERT(rhs.boxArray()==grids);

    const int sCompDivU = 0;
    const int nCompDivU = 1;
    const int sCompDsdt = 0;
    const int nCompDsdt = 1;

    if (have_divu)
    {
       FillPatch(*this,rhs,nGrow,time,Divu_Type,sCompDivU,nCompDivU,sCompDivU);
    }
    else
    {
       rhs.setVal(0.);
    }

    if (have_dsdt)
    {
       FillPatchIterator fpi(*this,rhs,nGrow,time,Dsdt_Type,sCompDsdt,nCompDsdt);
       const MultiFab& mf = fpi.get_mf();
       MultiFab::Saxpy(rhs, 0.5*dt, mf, 0, sCompDsdt, nCompDsdt, nGrow);
    }
}

void
NavierStokesBase::create_umac_grown (int nGrow)
{
    BL_PROFILE("NavierStokesBase::create_umac_grown()");

    if (level > 0)
    {
        BoxList bl = amrex::GetBndryCells(grids,nGrow);

        BoxArray f_bnd_ba(std::move(bl));

        BoxArray c_bnd_ba = f_bnd_ba; c_bnd_ba.coarsen(crse_ratio);

        c_bnd_ba.maxSize(32);

        f_bnd_ba = c_bnd_ba; f_bnd_ba.refine(crse_ratio);

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            //
            // crse_src & fine_src must have same parallel distribution.
            // We'll use the KnapSack distribution for the fine_src_ba.
            // Since fine_src_ba should contain more points, this'll lead
            // to a better distribution.
            //
            BoxArray crse_src_ba(c_bnd_ba), fine_src_ba(f_bnd_ba);

            crse_src_ba.surroundingNodes(idim);
            fine_src_ba.surroundingNodes(idim);

            const int N = fine_src_ba.size();

            std::vector<long> wgts(N);

#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < N; i++)
                wgts[i] = fine_src_ba[i].numPts();

            DistributionMapping dm;
            // This DM won't be put into the cache.
            dm.KnapSackProcessorMap(wgts,ParallelDescriptor::NProcs());

            // FIXME
            // Declaring in this way doesn't work. I think it's because the box arrays
            // have been changed and each src box is not completely contained within a
            // single box in the Factory's BA
            // For now, coarse-fine boundary doesn't intersect EB, so should be okay...
            // MultiFab crse_src(crse_src_ba, dm, 1, 0, MFInfo(), getLevel(level-1).Factory());
            // MultiFab fine_src(fine_src_ba, dm, 1, 0, MFInfo(), Factory());
            MultiFab crse_src(crse_src_ba, dm, 1, 0);
            MultiFab fine_src(fine_src_ba, dm, 1, 0);

            crse_src.setVal(1.e200);
            fine_src.setVal(1.e200);
            //
            // We want to fill crse_src from lower level u_mac including u_mac's grow cells.
            //
            const MultiFab& u_macLL = getLevel(level-1).u_mac[idim];
            crse_src.copy(u_macLL,0,0,1,1,0);

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(crse_src); mfi.isValid(); ++mfi)
            {
                const Box& box   = crse_src[mfi].box();
                auto const& crs_arr  = crse_src.array(mfi);
                auto const& fine_arr = fine_src.array(mfi);
                amrex::GpuArray<int,AMREX_SPACEDIM> c_ratio = {D_DECL(crse_ratio[0],crse_ratio[1],crse_ratio[2])};
                ParallelFor(box,[crs_arr,fine_arr,idim,c_ratio]
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                   int idx[3] = {i*c_ratio[0],j*c_ratio[1],k*c_ratio[2]};
#if ( AMREX_SPACEDIM == 2 )
                   // dim1 are the complement of idim
                   int dim1 = ( idim == 0 ) ? 1 : 0;
                   for (int n1 = 0; n1 < c_ratio[dim1]; n1++) {
                      int id[3] = {idx[0],idx[1],idx[2]};
                      id[dim1] += n1;
                      fine_arr(id[0],id[1],id[2]) = crs_arr(i,j,k);
                   }
#elif ( AMREX_SPACEDIM == 3 )
                   // dim1 and dim2 are the complements of idim
                   int dim1 = ( idim != 0 ) ? 0 : 1 ;
                   int dim2 = ( idim != 0 ) ? ( ( idim == 2 ) ? 1 : 2 ) : 2 ;
                   for (int n1 = 0; n1 < c_ratio[dim1]; n1++) {
                      for (int n2 = 0; n2 < c_ratio[dim2]; n2++) {
                         int id[3] = {idx[0],idx[1],idx[2]};
                         id[dim1] += n1;
                         id[dim2] += n2;
                         fine_arr(id[0],id[1],id[2]) = crs_arr(i,j,k);
                      }
                   }
#endif
                });
            }
            crse_src.clear();
            //
            // Replace pc-interpd fine data with preferred u_mac data at
            // this level u_mac valid only on surrounding faces of valid
            // region - this op will not fill grow region.
            //
            fine_src.copy(u_mac[idim]);
            //
            // Interpolate unfilled grow cells using best data from
            // surrounding faces of valid region, and pc-interpd data
            // on fine edges overlaying coarse edges.
            //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(fine_src); mfi.isValid(); ++mfi)
            {
                const int  nComp = 1;
                const Box& fbox  = fine_src[mfi].box();
                auto const& fine_arr = fine_src.array(mfi);
                amrex::GpuArray<int,AMREX_SPACEDIM> c_ratio = {D_DECL(crse_ratio[0],crse_ratio[1],crse_ratio[2])};
                amrex::launch(fbox, [fine_arr,nComp,c_ratio,idim]
                AMREX_GPU_DEVICE (Box const& tbx) {
                    edge_interp_k(tbx, nComp, idim, c_ratio, fine_arr);
                });
            }

            MultiFab u_mac_save(u_mac[idim].boxArray(),u_mac[idim].DistributionMap(),1,0,MFInfo(),Factory());
            u_mac_save.copy(u_mac[idim]);
            u_mac[idim].copy(fine_src,0,0,1,0,nGrow);
            u_mac[idim].copy(u_mac_save);
        }
    }
    for (int n = 0; n < BL_SPACEDIM; ++n)
    {
	u_mac[n].FillBoundary(geom.periodicity());
    }
}

void
NavierStokesBase::diffuse_scalar_setup (int sigma, int& rho_flag)
{

    rho_flag = Diffusion::set_rho_flag(diffusionType[sigma]);
}

void
NavierStokesBase::errorEst (TagBoxArray& tags,
			    int          clearval,
			    int          tagval,
			    Real         time,
			    int          n_error_buf,
			    int          ngrow)
{
    const int*  domain_lo = geom.Domain().loVect();
    const int*  domain_hi = geom.Domain().hiVect();
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();

    for (int j = 0; j < err_list.size(); j++)
    {
        auto mf = derive(err_list[j].name(), time, err_list[j].nGrow());
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
        {
	          const Box&  vbx     = mfi.tilebox();
            RealBox     gridloc = RealBox(vbx,geom.CellSize(),geom.ProbLo());
            Vector<int>  itags   = tags[mfi].tags();
            int*        tptr    = itags.dataPtr();
            const int*  tlo     = tags[mfi].box().loVect();
            const int*  thi     = tags[mfi].box().hiVect();
            const int*  lo      = vbx.loVect();
            const int*  hi      = vbx.hiVect();
            const Real* xlo     = gridloc.lo();
            FArrayBox&  fab     = (*mf)[mfi];
            Real*       dat     = fab.dataPtr();
            const int*  dlo     = fab.box().loVect();
            const int*  dhi     = fab.box().hiVect();
            const int   ncomp   = fab.nComp();

            err_list[j].errFunc()(tptr, ARLIM(tlo), ARLIM(thi), &tagval,
                                  &clearval, dat, ARLIM(dlo), ARLIM(dhi),
                                  lo,hi, &ncomp, domain_lo, domain_hi,
                                  dx, xlo, prob_lo, &time, &level);
            //
            // Don't forget to set the tags in the TagBox.
            //
            tags[mfi].tags(itags);
        }
    }

#ifdef AMREX_USE_EB
    // Enforce that the EB not cross the coarse-fine boundary
    const auto& ebfactory = dynamic_cast<amrex::EBFArrayBoxFactory const&>(Factory());
    if ( !ebfactory.isAllRegular() )
    {
      //
      // FIXME - For now, always refine cut cells
      //   Later, figure out a slick way to check if EB and CFB cross
      //   and allow !refine_cutcells
      //
      if (!refine_cutcells) amrex::Abort("For now, cutcells must always exist at finest level.");

      // Refine on cut cells
      if (refine_cutcells) // or if EB and CBF cross
      {
        const MultiFab& S_new = get_new_data(State_Type);
        amrex::TagCutCells(tags, S_new);
      }
    }
#endif
}


//
// Estimate the maximum allowable timestep at a cell center.
//
Real
NavierStokesBase::estTimeStep ()
{
    BL_PROFILE("NavierStokesBase::estTimeStep()");

    if (fixed_dt > 0.0)
    {
        Real factor = 1.0;

        if (!(level == 0))
        {
            int ratio = 1;
            for (int lev = 1; lev <= level; lev++)
            {
                ratio *= parent->nCycle(lev);
            }
            factor = 1.0/double(ratio);
        }

        return factor*fixed_dt;
    }

    const Real  small         = 1.0e-8;
    Real        estdt         = 1.0e+20;

    const Real  cur_pres_time = state[Press_Type].curTime();
    MultiFab&   S_new         = get_new_data(State_Type);

    Vector<Real> u_max(AMREX_SPACEDIM);
    Vector<Real> f_max(AMREX_SPACEDIM);

#ifdef AMREX_USE_EB
    MultiFab& Gp = getGradP();
    Gp.FillBoundary(geom.periodicity());
#else
    MultiFab Gp(grids,dmap,AMREX_SPACEDIM,1);
    getGradP(Gp, cur_pres_time);
#endif

    //
    // Find local max of velocity
    //
    u_max = S_new.norm0({AMREX_D_DECL(0,1,2)},0,true,true);

    //
    // Compute forcing terms: in this case this means external forces and grad(p)
    // Viscous terms not included since Crack-Nicholson is unconditionally stable
    // so no need to account for explicit part of viscous term
    //
    MultiFab tforces(grids,dmap,AMREX_SPACEDIM,0,MFInfo(),Factory());

      const BoxArray& ba = S_new.boxArray();
      const DistributionMapping& dm = S_new.DistributionMap();
      int ncomp = S_new.nComp();
      int ngrow = S_new.nGrow();
      MultiFab rhs(ba,dm,ncomp,ngrow);
      rhs.setVal(0.0);
    

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rho_ctime,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       const auto& bx          = mfi.tilebox();
       const auto  cur_time    = state[State_Type].curTime();
             auto& tforces_fab = tforces[mfi];


            if (getForceVerbose) {
                amrex::Print() << "------------------" << '\n'
                               << "H - est Time Step:" << '\n'
                               << "------------------" << '\n';
	    }       
       getForce(tforces_fab,bx,0,0,AMREX_SPACEDIM,cur_time,S_new[mfi],S_new[mfi],rhs[mfi],Density);

       const auto& rho   = rho_ctime.array(mfi);  
       const auto& gradp = Gp.array(mfi); 
       const auto& force = tforces.array(mfi);
       amrex::ParallelFor(bx, [rho, gradp, force] 
       AMREX_GPU_DEVICE(int i, int j, int k) noexcept
       {
          Real rho_inv = 1.0/rho(i,j,k);
          for (int n = 0; n < AMREX_SPACEDIM; n++) {
             force(i,j,k,n) -= gradp(i,j,k,n);
             force(i,j,k,n) *= rho_inv;
          } 
       });
    }

    //
    // Find local max of tforces
    //
    f_max = tforces.norm0({AMREX_D_DECL(0,1,2)},0,true,true);

    //
    // Compute local estdt
    //
    const Real* dx = geom.CellSize();

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (u_max[idim] > small)
        {
            estdt = std::min(estdt, dx[idim]/u_max[idim]);
        }

        if (f_max[idim] > small)
        {
            estdt = std::min(estdt, std::sqrt(2.0*dx[idim]/f_max[idim]));
        }
    }

    //
    // Reduce estimated dt by CFL factor and find global min
    //
    ParallelDescriptor::ReduceRealMin(estdt);

    if ( estdt < 1.0e+20) {
      //
      // timestep estimation successful
      //
      estdt = estdt * cfl;
    }
    else if (init_dt > 0 ) {
      //
      // use init_dt, scale for amr level
      //
      Real factor = 1.0;

      if (!(level == 0))
      {
         int ratio = 1;
         for (int lev = 1; lev <= level; lev++)
         {
           ratio *= parent->nCycle(lev);
         }
         factor = 1.0/double(ratio);
      }

      estdt = factor*init_dt;
    } else {
      Print()<<"\nNavierStokesBase::estTimeStep() failed to provide a good timestep "
             <<"(probably because initial velocity field is zero with no external forcing).\n"
             <<"Use ns.init_dt to provide a reasonable timestep on coarsest level.\n"
             <<"Note that ns.init_shrink will be applied to init_dt."<<std::endl;
      amrex::Abort("\n");
    }

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        ParallelDescriptor::ReduceRealMax(u_max.dataPtr(), AMREX_SPACEDIM, IOProc);

        amrex::Print() << "estTimeStep :: \n" << "LEV = " << level << " UMAX = ";
        for (int k = 0; k < AMREX_SPACEDIM; k++)
        {
            amrex::Print() << u_max[k] << "  ";
        }
        amrex::Print() << '\n';

        if (getForceVerbose) {
           ParallelDescriptor::ReduceRealMax(f_max.dataPtr(), AMREX_SPACEDIM, IOProc);
           amrex::Print() << "        FMAX = ";
           for (int k = 0; k < AMREX_SPACEDIM; k++)
           {
              amrex::Print() << f_max[k] << "  ";
           }
           amrex::Print() << '\n';
        }
        Print()<<"estimated timestep: dt = "<<estdt<<std::endl;
    }

  return estdt;
}

const MultiFab&
NavierStokesBase::get_rho (Real time)
{
    const TimeLevel whichTime = which_time(State_Type,time);

    if (whichTime == AmrOldTime)
    {
        return rho_ptime;
    }
    else if (whichTime == AmrNewTime)
    {
        return rho_ctime;
    }
    else if (whichTime == Amr1QtrTime)
    {
        BL_ASSERT(rho_qtime);
        return *rho_qtime;
    }
    else if (whichTime == Amr3QtrTime)
    {
        BL_ASSERT(rho_tqtime);
        return *rho_tqtime;
    }
    else if (whichTime == AmrHalfTime)
    {
        return get_rho_half_time();
    }
    else
    {
        amrex::Error("NavierStokesBase::get_rho(): bad time");

        return rho_ptime; // Got to return something to shut up compiler.
    }
}

MultiFab&
NavierStokesBase::get_rho_half_time ()
{
    //
    // Fill it in when needed ...
    //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rho_half,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();
        auto const& rho_h = rho_half.array(mfi);
        auto const& rho_p = rho_ptime.array(mfi);
        auto const& rho_c = rho_ctime.array(mfi);
        amrex::ParallelFor(bx, [rho_h, rho_p, rho_c] 
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
           rho_h(i,j,k) = 0.5 * (rho_p(i,j,k) + rho_c(i,j,k));
        });
    }
    return rho_half;
}

//
// Fill patch divU.
//
MultiFab*
NavierStokesBase::getDivCond (int ngrow, Real time)
{
    MultiFab* divu = 0;

    if (!have_divu)
    {
        divu = new MultiFab(grids,dmap,1,ngrow,MFInfo(),Factory());

        divu->setVal(0);
    }
    else
    {
        divu = getState(ngrow,Divu_Type,0,1,time);
    }

    return divu;
}

//
// Fill patch dSdt.
//
MultiFab*
NavierStokesBase::getDsdt (int ngrow, Real time)
{
    MultiFab* dsdt = 0;

    if (!(have_dsdt && have_divu))
    {
        dsdt = new MultiFab(grids,dmap,1,ngrow,MFInfo(),Factory());

        dsdt->setVal(0);
    }
    else
    {
        dsdt = getState(ngrow,Dsdt_Type,0,1,time);
    }

    return dsdt;
}


void
NavierStokesBase::getGradP (MultiFab& gp, Real      time)
{
    BL_PROFILE("NavierStokesBase::getGradP()");

    const int   NGrow = gp.nGrow();
    MultiFab&   P_old = get_old_data(Press_Type);
    const Real* dx    = geom.CellSize();

    if (level > 0 && state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
    {
        //
        // We want to be sure the intersection of old and new grids is
        // entirely contained within gp.boxArray()
        //
        BL_ASSERT(gp.boxArray() == grids);

        {
            const BoxArray& pBA = state[Press_Type].boxArray();
            MultiFab pMF(pBA,dmap,1,NGrow);

            if (time == getLevel(level-1).state[Press_Type].prevTime() ||
                time == getLevel(level-1).state[Press_Type].curTime())
            {
                FillCoarsePatch(pMF,0,time,Press_Type,0,1,NGrow);
            }
            else
            {
                Real crse_time;

                if (time > getLevel(level-1).state[State_Type].prevTime())
                {
                    crse_time = getLevel(level-1).state[Press_Type].curTime();
                }
                else
                {
                    crse_time = getLevel(level-1).state[Press_Type].prevTime();
                }

                FillCoarsePatch(pMF,0,crse_time,Press_Type,0,1,NGrow);

                MultiFab dpdtMF(pBA,dmap,1,NGrow);

                FillCoarsePatch(dpdtMF,0,time,Dpdt_Type,0,1,NGrow);

                Real dt_temp = time - crse_time;

                dpdtMF.mult(dt_temp,0,1,NGrow);

                pMF.plus(dpdtMF,0,1,NGrow);
            }
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(gp, true); mfi.isValid(); ++mfi)
            {
              const Box& bx=mfi.growntilebox();
              Projection::getGradP(pMF[mfi],gp[mfi],bx,dx);
            }
        }
        //
        // We've now got good coarse data everywhere in gp.
        //
        MultiFab gpTmp(gp.boxArray(),gp.DistributionMap(),1,NGrow);

        {
           FillPatchIterator P_fpi(*this,P_old,NGrow,time,Press_Type,0,1);
           MultiFab& pMF = P_fpi.get_mf();
#ifdef _OPENMP
#pragma omp parallel
#endif
           for (MFIter mfi(gpTmp, true); mfi.isValid(); ++mfi)
           {
             const Box& bx=mfi.growntilebox();
             Projection::getGradP(pMF[mfi],gpTmp[mfi],bx,dx);
           }
        }
        //
        // Now must decide which parts of gpTmp to copy to gp.
        //
        const int M = old_intersect_new.size();

        BoxArray fineBA(M);

        for (int j = 0; j < M; j++)
        {
            Box bx = old_intersect_new[j];

            for (int i = 0; i < BL_SPACEDIM; i++)
            {
                if (!geom.isPeriodic(i))
                {
                    if (bx.smallEnd(i) == geom.Domain().smallEnd(i))
                        bx.growLo(i,NGrow);
                    if (bx.bigEnd(i) == geom.Domain().bigEnd(i))
                        bx.growHi(i,NGrow);
                }
            }

            fineBA.set(j,bx);
        }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(gpTmp,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            auto isects = fineBA.intersections(mfi.growntilebox());
            auto const& gp_ar    = gp.array(mfi);
            auto const& gpTmp_ar = gpTmp.array(mfi);
            for (int ii = 0, N = isects.size(); ii < N; ii++)
            {
                const Box& ovlp = isects[ii].second;
                amrex::ParallelFor(ovlp, [gp_ar,gpTmp_ar]
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    gp_ar(i,j,k) = gpTmp_ar(i,j,k);
                });
            }
        }
        gp.EnforcePeriodicity(geom.periodicity());
    } else {
        FillPatchIterator P_fpi(*this,P_old,NGrow,time,Press_Type,0,1);
        MultiFab& pMF = P_fpi.get_mf();
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(gp, true); mfi.isValid(); ++mfi)
        {
           BL_ASSERT(amrex::grow(grids[mfi.index()],NGrow) == gp[mfi].box());
           Projection::getGradP(pMF[mfi],gp[mfi],mfi.growntilebox(),dx);
        }
    }
}

//
// Fill patch a state component.
//
MultiFab*
NavierStokesBase::getState (int  ngrow,
			    int  state_idx,
			    int  scomp,
			    int  ncomp,
			    Real time)
{
    BL_PROFILE("NavierStokesBase::getState()");

    MultiFab* mf = new MultiFab(state[state_idx].boxArray(),
                                state[state_idx].DistributionMap(),
                                ncomp,ngrow,MFInfo(),Factory());

    FillPatch(*this,*mf,ngrow,time,state_idx,scomp,ncomp,0);

    return mf;
}

void
NavierStokesBase::getOutFlowFaces (Vector<Orientation>& outFaces)
{
    outFaces.resize(0);
    for (int idir = 0; idir < BL_SPACEDIM; idir++)
    {
        if (phys_bc.lo(idir) == Outflow)
        {
            const int len = outFaces.size();
            outFaces.resize(len+1);
            outFaces[len] = Orientation(idir,Orientation::low);
        }

        if (phys_bc.hi(idir) == Outflow)
        {
            const int len = outFaces.size();
            outFaces.resize(len+1);
            outFaces[len] = Orientation(idir,Orientation::high);
        }
    }
}

void
NavierStokesBase::incrPAvg ()
{
    //
    // Increment p_avg with 1/ncycle times current pressure
    //
    MultiFab& P_new = get_new_data(Press_Type);

    Real alpha = 1.0/Real(parent->nCycle(level));

    MultiFab::Saxpy(p_avg,alpha,P_new,0,0,1,0);
}

void
NavierStokesBase::initRhoAvg (Real alpha)
{
    const MultiFab& S_new = get_new_data(State_Type);

    // Set to a ridiculous number just for debugging -- shouldn't need this otherwise
    rho_avg.setVal(1.e200);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rho_avg,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.tilebox();
       auto const& rhoavg     = rho_avg.array(mfi);
       auto const& rho_new    = S_new.array(mfi,Density);
       amrex::ParallelFor(bx, [rhoavg,rho_new,alpha]
       AMREX_GPU_DEVICE(int i, int j, int k) noexcept
       {
          rhoavg(i,j,k) = rho_new(i,j,k) * alpha;
       });
    }
}

void
NavierStokesBase::incrRhoAvg(const MultiFab& rho_incr,
                         int             sComp,
                         Real            alpha)
{
    MultiFab::Saxpy(rho_avg,alpha,rho_incr,sComp,0,1,0);
}

void
NavierStokesBase::incrRhoAvg (Real alpha)
{
    const MultiFab& S_new = get_new_data(State_Type);
    incrRhoAvg(S_new,Density,alpha);
}

//
// Fills a new level n with best level n and coarser data available.
//
void
NavierStokesBase::init (AmrLevel &old)
{
    NavierStokesBase* oldns = (NavierStokesBase*) &old;
    const Real    dt_new    = parent->dtLevel(level);
    const Real    cur_time  = oldns->state[State_Type].curTime();
    const Real    prev_time = oldns->state[State_Type].prevTime();
    const Real    dt_old    = cur_time - prev_time;
    MultiFab&     S_new     = get_new_data(State_Type);
    MultiFab&     P_new     = get_new_data(Press_Type);
    MultiFab&     P_old     = get_old_data(Press_Type);

    setTimeLevel(cur_time,dt_old,dt_new);

    const Real cur_pres_time = state[Press_Type].curTime();
    //
    // Get best state and pressure data.
    //
    FillPatch(old,S_new,0,cur_time,State_Type,0,NUM_STATE);
    //
    // Note: we don't need to worry here about using FillPatch because
    //       it will automatically use the "old dpdt" to interpolate,
    //       since we haven't yet defined a new pressure at the lower level.
    //
    {
       FillPatchIterator fpi(old,P_new,0,cur_pres_time,Press_Type,0,1);
       const MultiFab& mf_fpi = fpi.get_mf();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter mfi(mf_fpi,TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
         const Box& bx  = mfi.tilebox();
          const auto& p_arr = mf_fpi.array(mfi);
          const auto& p_o = P_old.array(mfi);
          const auto& p_n = P_new.array(mfi);
          amrex::ParallelFor(bx, [p_arr, p_o, p_n] 
          AMREX_GPU_DEVICE(int i, int j, int k) noexcept
          {
             p_o(i,j,k) = p_arr(i,j,k);
             p_n(i,j,k) = p_arr(i,j,k);
          });
       }
    }

    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
    {
        MultiFab& Dpdt_new = get_new_data(Dpdt_Type);
        FillPatch(old,Dpdt_new,0,cur_pres_time,Dpdt_Type,0,1);
    }
    //
    // Get best divu and dSdt data.
    //
    if (have_divu)
    {
        MultiFab& Divu_new = get_new_data(Divu_Type);
        FillPatch(old,Divu_new,0,cur_time,Divu_Type,0,1);

        if (have_dsdt)
        {
            MultiFab& Dsdt_new = get_new_data(Dsdt_Type);
            FillPatch(old,Dsdt_new,0,cur_time,Dsdt_Type,0,1);
        }
    }

    old_intersect_new          = amrex::intersect(grids,oldns->boxArray());
    is_first_step_after_regrid = true;
}

void
NavierStokesBase::init ()
{
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& P_new = get_new_data(Press_Type);

    BL_ASSERT(level > 0);

    const Vector<Real>& dt_amr = parent->dtLevel();
    Vector<Real>        dt_new(level+1);

    for (int lev = 0; lev < level; lev++)
        dt_new[lev] = dt_amr[lev];
    //
    // Guess new dt from new data (interpolated from coarser level).
    //
    const Real dt = dt_new[level-1]/Real(parent->MaxRefRatio(level-1));
    dt_new[level] = dt;

    parent->setDtLevel(dt_new);
    //
    // Compute dt based on old data.
    //
    NavierStokesBase& old   = getLevel(level-1);
    const Real    cur_time  = old.state[State_Type].curTime();
    const Real    prev_time = old.state[State_Type].prevTime();
    const Real    dt_old    = (cur_time-prev_time)/Real(parent->MaxRefRatio(level-1));

    setTimeLevel(cur_time,dt_old,dt);

    Real cur_pres_time = state[Press_Type].curTime();
    //
    // Get best coarse state and pressure data.
    //
    FillCoarsePatch(S_new,0,cur_time,State_Type,0,NUM_STATE);
    FillCoarsePatch(P_new,0,cur_pres_time,Press_Type,0,1);

    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
        FillCoarsePatch(get_new_data(Dpdt_Type),0,cur_time,Dpdt_Type,0,1);

    initOldPress();

    //
    // Get best coarse divU and dSdt data.
    //
    if (have_divu)
    {
        FillCoarsePatch(get_new_data(Divu_Type),0,cur_time,Divu_Type,0,1);
        if (have_dsdt)
            FillCoarsePatch(get_new_data(Dsdt_Type),0,cur_time,Dsdt_Type,0,1);
    }
    old_intersect_new = grids;
}

void
NavierStokesBase::init_additional_state_types ()
{
    additional_state_types_initialized = 1;
    //
    // Set "Temp" from user's variable setup.
    //
    int dummy_State_Type;
    int have_temp = isStateVariable("temp", dummy_State_Type, Temp);
    have_temp &= (dummy_State_Type == State_Type);
    BL_ASSERT((do_temp && have_temp)  ||  (!do_temp && !have_temp));

    int _Divu = -1;
    int dummy_Divu_Type;
    have_divu = 0;
    have_divu = isStateVariable("divu", dummy_Divu_Type, _Divu);
    have_divu = have_divu && dummy_Divu_Type == Divu_Type;
    if (verbose)
    {
        amrex::Print() << "NavierStokesBase::init_additional_state_types()::have_divu = "
                  << have_divu << '\n';
    }
    if (have_divu && _Divu!=Divu)
    {
        amrex::Print() << "divu must be 0-th Divu_Type component in the state\n";

        amrex::Abort("NavierStokesBase::init_additional_state_types()");
    }

    int _Dsdt = -1;
    int dummy_Dsdt_Type;
    have_dsdt = 0;
    have_dsdt = isStateVariable("dsdt", dummy_Dsdt_Type, _Dsdt);
    have_dsdt = have_dsdt && dummy_Dsdt_Type==Dsdt_Type;
    if (verbose)
    {
        amrex::Print() << "NavierStokesBase::init_additional_state_types()::have_dsdt = "
		       << have_dsdt << '\n';
    }
    if (have_dsdt && _Dsdt!=Dsdt)
    {
        amrex::Print() << "dsdt must be 0-th Dsdt_Type component in the state\n";

        amrex::Abort("NavierStokesBase::init_additional_state_types()");
    }
    if (have_dsdt && !have_divu)
    {
        amrex::Print() << "Must have divu in order to have dsdt\n";

        amrex::Abort("NavierStokesBase::init_additional_state_types()");
    }

    num_state_type = desc_lst.size();
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << "NavierStokesBase::init_additional_state_types: num_state_type = "
		       << num_state_type << '\n';
    }
}

Real
NavierStokesBase::initialTimeStep ()
{
    Real returnDt = init_shrink*estTimeStep();

    amrex::Print() << "Multiplying dt by init_shrink: dt = "
		   << returnDt << '\n';
    return returnDt;
}

//
// Since the pressure solver always stores its estimate of the
// pressure solver in Pnew, we need to copy it to Pold at the start.
//
void
NavierStokesBase::initOldPress ()
{
    MultiFab& P_new = get_new_data(Press_Type);
    MultiFab& P_old = get_old_data(Press_Type);

    MultiFab::Copy(P_old, P_new, 0, 0, P_old.nComp(), P_old.nGrow());
}

void
NavierStokesBase::zeroNewPress ()
{
    get_new_data(Press_Type).setVal(0);
}

void
NavierStokesBase::zeroOldPress ()
{
    get_old_data(Press_Type).setVal(0);
}

void
NavierStokesBase::level_projector (Real dt,
				   Real time,
				   int  iteration)
{
    BL_PROFILE_REGION_START("R::NavierStokesBase::level_projector()");
    BL_PROFILE("NavierStokesBase::level_projector()");

    BL_ASSERT(iteration > 0);

    MultiFab& U_old = get_old_data(State_Type);
    MultiFab& U_new = get_new_data(State_Type);
    MultiFab& P_old = get_old_data(Press_Type);
    MultiFab& P_new = get_new_data(Press_Type);

    SyncRegister* crse_ptr = 0;

    if (level < parent->finestLevel() && do_sync_proj)
    {
        crse_ptr = &(getLevel(level+1).getSyncReg());
    }

    int        crse_dt_ratio  = (level > 0) ? parent->nCycle(level) : -1;
    const Real cur_pres_time  = state[Press_Type].curTime();
    const Real prev_pres_time = state[Press_Type].prevTime();

    projector->level_project(level,time,dt,cur_pres_time,prev_pres_time,
                             geom,U_old,U_new,P_old,P_new,
                             get_rho_half_time(),crse_ptr,sync_reg,
                             crse_dt_ratio,iteration,have_divu);

    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
        calcDpdt();

    BL_PROFILE_REGION_STOP("R::NavierStokesBase::level_projector()");
}

void
NavierStokesBase::level_sync (int crse_iteration)
{
    BL_PROFILE_REGION_START("R::NavierStokesBase::level_sync()");
    BL_PROFILE("NavierStokesBase::level_sync()");

    IntVect         ratio         = parent->refRatio(level);
    const int       finest_level  = parent->finestLevel();
    int             crse_dt_ratio = parent->nCycle(level);
    Real            dt            = parent->dtLevel(level);
    MultiFab&       pres          = get_new_data(Press_Type);
    MultiFab&       vel           = get_new_data(State_Type);
    SyncRegister&   rhs_sync_reg  = getLevel(level+1).getSyncReg();
    SyncRegister*   crsr_sync_ptr = 0;
    NavierStokesBase&   fine_level    = getLevel(level+1);
    MultiFab&       pres_fine     = fine_level.get_new_data(Press_Type);
    MultiFab&       vel_fine      = fine_level.get_new_data(State_Type);
    const BoxArray& finegrids     = vel_fine.boxArray();
    const DistributionMapping& finedmap = vel_fine.DistributionMap();

    if (level > 0)
        crsr_sync_ptr = &(getLevel(level).getSyncReg());
    //
    // Get boundary conditions.
    //
    const int N = grids.size();

    Vector<int*>         sync_bc(N);
    Vector< Vector<int> > sync_bc_array(N);

    for (int i = 0; i < N; i++)
    {
        sync_bc_array[i] = getBCArray(State_Type,i,Xvel,BL_SPACEDIM);
        sync_bc[i] = sync_bc_array[i].dataPtr();
    }

    //
    // Multilevel or single-level sync projection.
    //
    MultiFab& Rh = get_rho_half_time();
    MultiFab cc_rhs_crse, cc_rhs_fine;

    cc_rhs_crse.define(    grids,    dmap,1,1,MFInfo(),           Factory());
    cc_rhs_fine.define(finegrids,finedmap,1,1,MFInfo(),fine_level.Factory());
    cc_rhs_crse.setVal(0);
    cc_rhs_fine.setVal(0);

    MultiFab&         v_fine    = fine_level.get_new_data(State_Type);
    MultiFab&       rho_fine    = fine_level.rho_avg;
    const Geometry& crse_geom   = parent->Geom(level);
    const BoxArray& P_finegrids = pres_fine.boxArray();
    const DistributionMapping& P_finedmap = pres_fine.DistributionMap();

    MultiFab phi(P_finegrids,P_finedmap,1,1,MFInfo(),fine_level.Factory());
    MultiFab V_corr(finegrids,finedmap,BL_SPACEDIM,1,MFInfo(),fine_level.Factory());

    V_corr.setVal(0);
    //
    // If periodic, enforce periodicity on Vsync.
    //
    if (crse_geom.isAnyPeriodic()) {
      Vsync.FillBoundary(0, BL_SPACEDIM, crse_geom.periodicity());
    }
    //
    // Interpolate Vsync to fine grid correction in Vcorr.
    //
    SyncInterp(Vsync, level, V_corr, level+1, ratio,
	       0, 0, BL_SPACEDIM, 0 , dt, sync_bc.dataPtr());
    //
    // The multilevel projection.  This computes the projection and
    // adds in its contribution to levels (level) and (level+1).
    //
    Real  cur_crse_pres_time = state[Press_Type].curTime();
    Real prev_crse_pres_time = state[Press_Type].prevTime();

    NavierStokesBase& fine_lev   = getLevel(level+1);
    Real  cur_fine_pres_time = fine_lev.state[Press_Type].curTime();
    Real prev_fine_pres_time = fine_lev.state[Press_Type].prevTime();

    bool first_crse_step_after_initial_iters =
      (prev_crse_pres_time > state[State_Type].prevTime());

    bool pressure_time_is_interval =
      (state[Press_Type].descriptor()->timeType() == StateDescriptor::Interval);
    projector->MLsyncProject(level,pres,vel,cc_rhs_crse,
			     pres_fine,v_fine,cc_rhs_fine,
			     Rh,rho_fine,Vsync,V_corr,
			     phi,&rhs_sync_reg,crsr_sync_ptr,
			     dt,ratio,crse_iteration,crse_dt_ratio,
			     geom,pressure_time_is_interval,
			     first_crse_step_after_initial_iters,
			     cur_crse_pres_time,prev_crse_pres_time,
			     cur_fine_pres_time,prev_fine_pres_time);
    cc_rhs_crse.clear();
    cc_rhs_fine.clear();
    //
    // Correct pressure and velocities after the projection.
    //
    const int Nf = finegrids.size();

    ratio = IntVect::TheUnitVector();

    Vector<int*>         fine_sync_bc(Nf);
    Vector< Vector<int> > fine_sync_bc_array(Nf);

    for (int i = 0; i < Nf; i++)
    {
      fine_sync_bc_array[i] = getLevel(level+1).getBCArray(State_Type,
							   i,
							   Xvel,
							   BL_SPACEDIM);
      fine_sync_bc[i] = fine_sync_bc_array[i].dataPtr();
    }

    for (int lev = level+2; lev <= finest_level; lev++)
    {
      ratio                 *= parent->refRatio(lev-1);
      NavierStokesBase& flev = getLevel(lev);
      MultiFab&     P_new    = flev.get_new_data(Press_Type);
      MultiFab&     P_old    = flev.get_old_data(Press_Type);
      MultiFab&     U_new    = flev.get_new_data(State_Type);

      SyncInterp(V_corr, level+1, U_new, lev, ratio,
		 0, 0, BL_SPACEDIM, 1 , dt, fine_sync_bc.dataPtr());
      SyncProjInterp(phi, level+1, P_new, P_old, lev, ratio,
		     first_crse_step_after_initial_iters,
		     cur_crse_pres_time, prev_crse_pres_time);
    }

    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
      calcDpdt();

    BL_PROFILE_REGION_STOP("R::NavierStokesBase::level_sync()");
}

void
NavierStokesBase::make_rho_prev_time ()
{
    const Real prev_time = state[State_Type].prevTime();

    FillPatch(*this,rho_ptime,1,prev_time,State_Type,Density,1,0);

#ifdef AMREX_USE_EB
    EB_set_covered(rho_ptime,COVERED_VAL);
#endif
}

void
NavierStokesBase::make_rho_curr_time ()
{
    const Real curr_time = state[State_Type].curTime();

    FillPatch(*this,rho_ctime,1,curr_time,State_Type,Density,1,0);

#ifdef AMREX_USE_EB
    EB_set_covered(rho_ctime,COVERED_VAL);
#endif
}

void
NavierStokesBase::mac_project (Real      time,
                               Real      dt,
                               MultiFab& Sold,
                               MultiFab* divu,
                               int       ngrow,
                               bool      increment_vel_register)
{
    BL_PROFILE_REGION_START("R::NavierStokesBase::mac_project()");
    BL_PROFILE("NavierStokesBase::mac_project()");

    if (verbose) amrex::Print() << "... mac_projection\n";

    if (verbose && benchmarking) ParallelDescriptor::Barrier();

    const Real strt_time = ParallelDescriptor::second();

    Vector<BCRec> density_math_bc = fetchBCArray(State_Type,Density,1);

    mac_projector->mac_project(level,u_mac,Sold,dt,time,*divu,have_divu,
                               density_math_bc[0], increment_vel_register);

    create_umac_grown(ngrow);

    if (verbose)
    {
        Real run_time    = ParallelDescriptor::second() - strt_time;
        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        amrex::Print() << "NavierStokesBase:mac_project(): lev: "
                       << level
                       << ", time: " << run_time << '\n';
    }
    BL_PROFILE_REGION_STOP("R::NavierStokesBase::mac_project()");
}

void
NavierStokesBase::manual_tags_placement (TagBoxArray&    tags,
                                         const Vector<IntVect>& bf_lev)
{
    Vector<Orientation> outFaces;
    getOutFlowFaces(outFaces);
    if (outFaces.size()>0)
    {
        for (int i=0; i<outFaces.size(); ++i)
        {
            const Orientation& outFace = outFaces[i];
            const int oDir = outFace.coordDir();
            const Box& crse_domain = amrex::coarsen(geom.Domain(),bf_lev[level]);
            const int mult = (outFace.isLow() ? +1 : -1);
            if (do_refine_outflow)
            {
                //
                // Refine entire outflow boundary if new boxes within grid_tol
                // from outflow
                //
                const int grid_tol = 1;

                Box outflowBox = amrex::adjCell(crse_domain,outFace,grid_tol);

                outflowBox.shift(oDir,mult*grid_tol);

                //
                // Only refine if there are already tagged cells in the outflow
                // region
                //
                bool hasTags = tags.hasTags(outflowBox);
                if (hasTags)
                    tags.setVal(BoxArray(&outflowBox,1),TagBox::SET);
	    }
            else if (do_derefine_outflow)
            {
                const int np = parent->nProper();
                //
                // Calculate the number of level 0 cells to be left uncovered
                // at the outflow.  The convoluted logic allows for the fact that
                // the number of uncovered cells must be a multiple of the level
                // blocking factor.  So, when calculating the number of coarse
                // cells below, we always round the division up.
                //
                int N_coarse_cells = Nbuf_outflow / bf_lev[0][oDir];
                if (Nbuf_outflow % bf_lev[0][oDir] != 0)
                    N_coarse_cells++;

                int N_level_cells = N_coarse_cells * bf_lev[0][oDir];

                //
                // Adjust this to get the number of cells to be left uncovered at
                // levels higher than 0
                //
                for (int j = 1; j <= level; ++j)
                {
                    /*** Calculate the minimum cells at this level ***/

                    const int rat = (parent->refRatio(j-1))[oDir];
                    N_level_cells = N_level_cells * rat + np;

                    /*** Calculate the required number of coarse cells ***/

                    N_coarse_cells = N_level_cells / bf_lev[j][oDir];
                    if (N_level_cells % bf_lev[j][oDir] != 0)
                        N_coarse_cells++;

                    /*** Calculate the corresponding number of level cells ***/

                    N_level_cells = N_coarse_cells * bf_lev[j][oDir];
                }
                //
                // Untag the cells near the outflow
                //
                if (N_coarse_cells > 0)
                {
                    //
                    // Generate box at the outflow and grow it in all directions
                    // other than the outflow.  This forces outflow cells in the
                    // ghostcells in directions other that oDir to be cleared.
                    //
                    Box outflowBox = amrex::adjCell(crse_domain, outFace, 1);
                    for (int dir = 0; dir < BL_SPACEDIM; dir++)
                        if (dir != oDir) outflowBox.grow(dir, 1);
                    //
                    // Now, grow the box into the domain (opposite direction as
                    // outFace) the number of cells we need to clear.
                    //
                    if (outFace.isLow())
                        outflowBox.growHi(oDir, N_coarse_cells);
                    else
                        outflowBox.growLo(oDir, N_coarse_cells);

                    tags.setVal(BoxArray(&outflowBox,1),TagBox::CLEAR);
                }
            }
        }
    }
}

int
NavierStokesBase::okToContinue ()
{
   //
   // Check that dt is OK across AMR levels
   //
   int okLevel = (level > 0) ? true : (parent->dtLevel(0) > dt_cutoff);

   if (stop_when_steady)
      //
      // If stop_when_steady is enabled, also check that we haven't reached
      // steady-state.
      //
      return (okLevel && !steadyState());
   else
      return okLevel;
}

int
NavierStokesBase::steadyState()
{
    if (!get_state_data(State_Type).hasOldData()) {
        return false; // If nothing to compare to, must not yet be steady :)
    }

    MultiFab&   u_old = get_old_data(State_Type);
    MultiFab&   u_new = get_new_data(State_Type);

	//
	// Estimate the maximum change in velocity magnitude since previous
	// iteration
	//
    Real max_change = 0.0;

    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    // Do not OpenMP-fy this loop for now
    // Unclear how to keep OpenMP and GPU implementation
    // from messing with each other
    for (MFIter mfi(u_old,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto& bx   = mfi.tilebox();
        const auto& uold = u_old[mfi].array();
        const auto& unew = u_new[mfi].array();

        reduce_op.eval(bx, reduce_data, [uold, unew]
        AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real uold_mag = 0.0;
            Real unew_mag = 0.0;
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
            {
                uold_mag += uold(i,j,k,d)*uold(i,j,k,d);
                unew_mag += unew(i,j,k,d)*unew(i,j,k,d);
            }

            uold_mag = std::sqrt(uold_mag);
            unew_mag = std::sqrt(unew_mag);

            return std::abs(unew_mag-uold_mag);
        });

        max_change = std::max(amrex::get<0>(reduce_data.value()),
                              max_change);
    }

    ParallelDescriptor::ReduceRealMax(max_change);

	//
	// System is classified as steady if the maximum change is smaller than
	// prescribed tolerance
	//
    bool steady = max_change < steady_tol;

    if (verbose)
    {
        amrex::Print() << "steadyState :: \n" << "LEV = " << level
                       << " MAX_CHANGE = " << max_change << std::endl;

        if (steady)
        {
            amrex::Print()
                << "System reached steady-state, stopping simulation."
                << std::endl;
        }
    }

    return steady;
}

//
// This function estimates the initial timesteping used by the model.
//
void
NavierStokesBase::post_init_estDT (Real&        dt_init,
				   Vector<int>&  nc_save,
				   Vector<Real>& dt_save,
				   Real         stop_time)
{
    const Real strt_time    = state[State_Type].curTime();
    const int  finest_level = parent->finestLevel();

    dt_init = 1.0e+100;

    int  n_factor;
    for (int k = 0; k <= finest_level; k++)
    {
        nc_save[k] = parent->nCycle(k);
        dt_save[k] = getLevel(k).initialTimeStep();

        n_factor   = 1;
        for (int m = finest_level; m > k; m--)
             n_factor *= parent->nCycle(m);
        dt_init    = std::min( dt_init, dt_save[k]/((Real) n_factor) );
    }

    Vector<Real> dt_level(finest_level+1,dt_init);
    Vector<int>  n_cycle(finest_level+1,1);

    Real dt0 = dt_save[0];
    n_factor = 1;
    for (int k = 0; k <= finest_level; k++)
    {
        n_factor *= nc_save[k];
        dt0       = std::min(dt0,n_factor*dt_save[k]);
    }

    if (stop_time >= 0.0)
    {
        const Real eps = 0.0001*dt0;
        if ((strt_time + dt0) > (stop_time - eps))
            dt0 = stop_time - strt_time;
    }

    n_factor = 1;
    for (int k = 0; k <= finest_level; k++)
    {
        n_factor  *= nc_save[k];
        dt_save[k] = dt0/( (Real) n_factor);
    }
    //
    // Hack.
    //
    parent->setDtLevel(dt_level);
    parent->setNCycle(n_cycle);
    for (int k = 0; k <= finest_level; k++)
    {
        getLevel(k).setTimeLevel(strt_time,dt_init,dt_init);
    }
}

//
// This function ensures that the state is initially consistent
// with respect to the divergence condition and fields are initially consistent
//
void
NavierStokesBase::post_init_state ()
{
    const int finest_level = parent->finestLevel();
    const Real divu_time   = have_divu ? state[Divu_Type].curTime()
                                       : state[Press_Type].curTime();

    if (do_init_vort_proj)
    {
        //
	// NOTE: this assumes have_divu == 0.
	// Only used if vorticity is used to initialize the velocity field.
        //
        BL_ASSERT(!(projector == 0));

	if (verbose) amrex::Print() << "calling initialVorticityProject" << std::endl;

	projector->initialVorticityProject(0);

	if (verbose) amrex::Print() << "done calling initialVorticityProject" << std::endl;
    }

    if (do_init_proj && projector)
    {
      //
      // Do sync project to define divergence free velocity field.
      //
      if (verbose) amrex::Print() << "calling initialVelocityProject" << std::endl;

      projector->initialVelocityProject(0,divu_time,have_divu,init_vel_iter);

      if (verbose) amrex::Print() << "done calling initialVelocityProject" << std::endl;
    }

    NavierStokesBase::initial_step = true;
    //
    // Average velocity and scalar data down from finer levels
    // so that conserved data is consistant between levels.
    //
    for (int k = finest_level-1; k>= 0; k--)
    {
      getLevel(k).avgDown();
    }
    make_rho_curr_time();

    if (do_init_proj && projector && (std::fabs(gravity)) > 0.){
      //
      // Do projection to establish initially hydrostatic pressure field.
      //
      if (verbose) amrex::Print() << "calling initialPressureProject" << std::endl;

      projector->initialPressureProject(0);

      if (verbose) amrex::Print() << "done calling initialPressureProject" << std::endl;
    }
    // make sure there's not NANs in old pressure field
    // end up with P_old = P_new as is the case when exiting initialPressureProject
    if(!do_init_proj){
      MultiFab& p_old=get_old_data(Press_Type);
      MultiFab& p_new=get_new_data(Press_Type);
      MultiFab::Copy(p_old, p_new, 0, 0, 1, p_new.nGrow());
    }

}

//
// Build any additional data structures after regrid.
//
void
NavierStokesBase::post_regrid (int lbase,
                               int new_finest)
{
#ifdef AMREX_PARTICLES
    if (NSPC && level == lbase)
    {
        NSPC->Redistribute(lbase);
    }
#endif
}

//
// Build any additional data structures after restart.
//
void
NavierStokesBase::post_restart ()
{
    make_rho_prev_time();
    make_rho_curr_time();

#ifdef AMREX_PARTICLES
    post_restart_particle ();
#endif
}

//
// Integration cycle on fine level grids is complete .
// post_timestep() is responsible for syncing levels together.
//
// The registers used for level syncing are initialized in the
// coarse level advance and incremented in the fine level advance.
// These quantities are described in comments above advance_setup.
//
void
NavierStokesBase::post_timestep (int crse_iteration)
{

  BL_PROFILE("NavierStokesBase::post_timestep()");

    const int finest_level = parent->finestLevel();

#ifdef AMREX_PARTICLES
    post_timestep_particle (crse_iteration);
#endif

    if (level == parent->finestLevel())
    {
        delete [] u_mac;
        u_mac = 0;
    }

    if (do_reflux && level < finest_level)
        reflux();

    if (level < finest_level)
        avgDown();

    if (do_mac_proj && level < finest_level)
        mac_sync();

    if (do_sync_proj && (level < finest_level))
        level_sync(crse_iteration);
    //
    // Test for conservation.
    //
    if (level==0 && sum_interval>0 && (parent->levelSteps(0)%sum_interval == 0))
    {
        sum_integrated_quantities();
    }
#if (AMREX_SPACEDIM==3)
    //
    // Derive turbulent statistics
    //
    if (level==0 && turb_interval>0 && (parent->levelSteps(0)%turb_interval == 0))
    {
        sum_turbulent_quantities();
    }
#ifdef SUMJET
    //
    // Derive turbulent statistics for the round jet
    //
    if (level==0 && jet_interval>0 && (parent->levelSteps(0)%jet_interval == 0))
    {
        sum_jet_quantities();
    }
#endif
#endif

    if (level > 0) incrPAvg();

    old_intersect_new          = grids;
    is_first_step_after_regrid = false;

    if (level == 0 && dump_plane >= 0)
    {
        Box bx = geom.Domain();

        BL_ASSERT(bx.bigEnd(AMREX_SPACEDIM-1) >= dump_plane);

        bx.setSmall(AMREX_SPACEDIM-1, dump_plane);
        bx.setBig  (AMREX_SPACEDIM-1, dump_plane);

        BoxArray ba(bx);
        DistributionMapping dm{ba};

        MultiFab mf(ba, dm, AMREX_SPACEDIM, 0, MFInfo(), Factory());

        mf.copy(get_new_data(State_Type), Xvel, 0, AMREX_SPACEDIM);

        if (ParallelDescriptor::MyProc() == mf.DistributionMap()[0])
        {
            char buf[64];
            sprintf(buf, "%14.12e", state[State_Type].curTime());

            std::string name(dump_plane_name);
            name += buf;
            name += ".fab";

            std::ofstream ofs;
            ofs.open(name.c_str(),std::ios::out|std::ios::trunc|std::ios::binary);
            if (!ofs.good())
                amrex::FileOpenFailed(name);

            mf[0].writeOn(ofs);
        }
    }
}

//
// Reset the time levels to time (time) and timestep dt.
// This is done at the start of the timestep in the pressure iteration section.
//
void
NavierStokesBase::resetState (Real time,
                              Real dt_old,
                              Real dt_new)
{
    //
    // Reset state types.
    //
    state[State_Type].reset();
    state[State_Type].setTimeLevel(time,dt_old,dt_new);

    initOldPress();
    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Interval)
    {
        state[Press_Type].setTimeLevel(time-dt_old,dt_old,dt_new);
    }
    else if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
    {
        state[Press_Type].setTimeLevel(time-.5*dt_old,dt_old,dt_old);
        state[Dpdt_Type].setTimeLevel(time-dt_old,dt_old,dt_old);
    }
    //
    // Reset state types for divu not equal to zero.
    //
    if (have_divu)
    {
        state[Divu_Type].reset();
        state[Divu_Type].setTimeLevel(time,dt_old,dt_new);
        if (have_dsdt)
        {
            //
            // Dont do this, we want to improve dsdt with press iters
            // but we do need to make sure time is set correctly..
            // state[Dsdt_Type].reset();
            state[Dsdt_Type].setTimeLevel(time,dt_old,dt_new);
        }
    }
}

void
NavierStokesBase::restart (Amr&          papa,
                           std::istream& is,
                           bool          bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

#ifdef AMREX_USE_EB
    amrex::Warning("Restart not tested with EB yet.");
#endif
    //
    // Build metric coefficients for RZ calculations.
    // Build volume and areas.
    //
    buildMetrics();

    if (projector == 0)
    {
        projector = new Projection(parent,&phys_bc,do_sync_proj,
                                   parent->finestLevel(),radius_grow);
    }
    projector->install_level(level, this, &radius);

    if (mac_projector == 0)
    {
        mac_projector = new MacProj(parent,parent->finestLevel(),
                                    &phys_bc,radius_grow);
    }
    mac_projector->install_level(level,this);

    const BoxArray& P_grids = state[Press_Type].boxArray();
#ifdef AMREX_USE_EB
    init_eb(parent->Geom(level), grids, dmap);

    //fixme? not 100% sure this is the right place
    // note --- this fn is really similar to constructor
    //  need to make sure gradp is getting properly filled for this restart case?
    //  incflo style advection does not use Gp in tracing states to edges,
    //  don't think Gp is needed until the vel update after the projection.
    // But ultimately we need gradp in the checkpoint file.
    gradp.reset(new MultiFab(grids,dmap,BL_SPACEDIM,1, MFInfo(), Factory()));

    std::string file=papa.theRestartFile();
    std::string LevelDir, FullPath;
    LevelDirectoryNames(file, LevelDir, FullPath);
    std::string gradp_mf_fullpath = FullPath;
    gradp_mf_fullpath += "/gradp";
    const char *faHeader = 0;
    VisMF::Read(*gradp, gradp_mf_fullpath, faHeader);
#endif

    //
    // Alloc space for density and temporary pressure variables.
    //
    if (level > 0)
    {
        rho_avg.define(grids,dmap,1,1,MFInfo(),Factory());
        p_avg.define(P_grids,dmap,1,0,MFInfo(),Factory());
    }
    //FIXME see similar stuff in constructor
    rho_half.define (grids,dmap,1,1,MFInfo(),Factory());
    rho_ptime.define(grids,dmap,1,1,MFInfo(),Factory());
    rho_ctime.define(grids,dmap,1,1,MFInfo(),Factory());
    rho_qtime  = 0;
    rho_tqtime = 0;

    BL_ASSERT(sync_reg == 0);
    if (level > 0 && do_sync_proj)
    {
        sync_reg = new SyncRegister(grids,dmap,crse_ratio);
    }
    BL_ASSERT(advflux_reg == 0);
    if (level > 0 && do_reflux)
    {
        advflux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
    }
    BL_ASSERT(viscflux_reg == 0);
    if (level > 0 && do_reflux)
    {
        viscflux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
    }

    if (level < parent->finestLevel())
    {
        Vsync.define(grids,dmap,AMREX_SPACEDIM,1,MFInfo(),Factory());
        Ssync.define(grids,dmap,NUM_STATE-AMREX_SPACEDIM,1,MFInfo(),Factory());
    }

    diffusion = new Diffusion(parent, this,
                              (level > 0) ? getLevel(level-1).diffusion : 0,
                              NUM_STATE, viscflux_reg,is_diffusive, visc_coef);
    //
    // Allocate the storage for variable viscosity and diffusivity
    //
    diffn_cc = new MultiFab(grids, dmap, NUM_STATE-Density-1, 1, MFInfo(), Factory());
    diffnp1_cc = new MultiFab(grids, dmap, NUM_STATE-Density-1, 1, MFInfo(), Factory());
    viscn_cc = new MultiFab(grids, dmap, 1, 1, MFInfo(), Factory());
    viscnp1_cc = new MultiFab(grids, dmap, 1, 1, MFInfo(), Factory());

    is_first_step_after_regrid = false;
    old_intersect_new          = grids;
}

void
NavierStokesBase::scalar_advection_update (Real dt,
                                           int  first_scalar,
                                           int  last_scalar)
{
    BL_PROFILE("NavierStokesBase::scalar_advection_update()");

    MultiFab&  S_old     = get_old_data(State_Type);
    MultiFab&  S_new     = get_new_data(State_Type);
    MultiFab&  Aofs      = *aofs;

    const Real prev_time = state[State_Type].prevTime();


    //
    // Compute inviscid estimate of scalars.
    // (do rho separate, as we do not have rho at new time yet)
    //
    int sComp = first_scalar;

    if (sComp == Density)
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(S_old,TilingIfNotGPU()); mfi.isValid(); ++mfi)
{
	    const Box&  bx = mfi.tilebox();
            const auto& Snew = S_new[mfi].array(Density);
            const auto& Sold = S_old[mfi].const_array(Density);
            const auto& aofs = Aofs[mfi].const_array(Density);

            amrex::ParallelFor(bx, [ Snew, Sold, aofs, dt]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
                Snew(i,j,k) = Sold(i,j,k) - dt * aofs(i,j,k);
            });
        }

        //
        // Call ScalMinMax to avoid overshoots in density.
        //
      if (do_denminmax)
      {
            //
            // Must do FillPatch here instead of MF iterator because we need the
            // boundary values in the old data (especially at inflow)
            //
            const int index_new_s   = Density;
            const int index_new_rho = Density;
            const int index_old_s   = index_new_s   - Density;
            const int index_old_rho = index_new_rho - Density;

            FillPatchIterator S_fpi(*this,S_old,1,prev_time,State_Type,Density,1);
            MultiFab& Smf=S_fpi.get_mf();

            ConservativeScalMinMax(S_new, index_new_s, index_new_rho,
                                   Smf,   index_old_s, index_old_rho);

      }
      ++sComp;
    }



    
    if (sComp <= last_scalar)
    {
        const MultiFab& rho_halftime = get_rho_half_time();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
{


	const BoxArray& ba = S_old.boxArray();
        const DistributionMapping& dm = S_old.DistributionMap();
	int ncomp = S_old.nComp();
	int ngrow = S_old.nGrow();
	MultiFab rhs(ba,dm,ncomp,ngrow);
	rhs.setVal(0.0);

       if (S_old.contains_nan(sComp,1,0))
       {
	 Print() << "Scalar (AIII) " << sComp << " contains Nans" << '\n';
	 exit(0);
       }
       if (S_new.contains_nan(sComp,1,0))
       {
	 Print() << "New scalar (AIII) " << sComp << " contains Nans" << '\n';
	 exit(0);
       }

       /*
      amrex::AllPrint() << " xxxxx proc. = " << ParallelDescriptor::MyProc() << ", tag = " <<
      ParallelDescriptor::SeqNum()
		      << " file = " << __FILE__ << " function = " << __FUNCTION__
		      << " line = " << __LINE__ << std::endl;       
       */	

#ifdef AMREX_PARTICLES
	for (MFIter Rho_mfi(rho_halftime,true); Rho_mfi.isValid(); ++Rho_mfi)
        {
            const Box& bx = Rho_mfi.tilebox();

	    // Average the mac face velocities to get cell centred velocities
            const Real halftime = 0.5*(state[State_Type].curTime()+state[State_Type].prevTime());
            FArrayBox Vel(amrex::grow(bx,0),AMREX_SPACEDIM);
            FORT_AVERAGE_EDGE_STATES(BL_TO_FORTRAN_ANYD(Vel),
                                     BL_TO_FORTRAN_ANYD(u_mac[0][Rho_mfi]),
                                     BL_TO_FORTRAN_ANYD(u_mac[1][Rho_mfi]),
#if (AMREX_SPACEDIM==3)
                                     BL_TO_FORTRAN_ANYD(u_mac[2][Rho_mfi]),
#endif
                                     &getForceVerbose);

            // Average the new and old time to get Crank-Nicholson half time approximation.
            FArrayBox Scal(amrex::grow(bx,0),NUM_SCALARS);
            Scal.copy<RunOn::Host>(S_old[Rho_mfi],bx,Density,bx,0,NUM_SCALARS);
            Scal.plus<RunOn::Host>(S_new[Rho_mfi],bx,Density,0,NUM_SCALARS);
            Scal.mult<RunOn::Host>(0.5,bx);

            if (NavierStokes::initial_iter != true) {	    
	    theNSPC()->getTemp(rhs[Rho_mfi],Vel,Scal,visc_coef[0],ngrow,level);
	    }
	}
	rhs.SumBoundary(Geom().periodicity()); 	   		
#endif



  
  
        FArrayBox  tforces;

        for (MFIter Rho_mfi(rho_halftime,TilingIfNotGPU()); Rho_mfi.isValid(); ++Rho_mfi)
        {
            const Box& bx = Rho_mfi.tilebox();

            for (int sigma = sComp; sigma <= last_scalar; sigma++)
            {
               // Need to do some funky half-time stuff
	      if (getForceVerbose) {
                  amrex::Print() << "----------------------------------------\n" 
                                 << "E - scalar advection update (half time):\n"
                                 << "----------------------------------------\n";
	      }	       

               // Average the mac face velocities to get cell centred velocities
               const Real halftime = 0.5*(state[State_Type].curTime()+state[State_Type].prevTime());
               FArrayBox Vel(amrex::grow(bx,0),AMREX_SPACEDIM);
               FORT_AVERAGE_EDGE_STATES(BL_TO_FORTRAN_ANYD(Vel),
                                        BL_TO_FORTRAN_ANYD(u_mac[0][Rho_mfi]),
                                        BL_TO_FORTRAN_ANYD(u_mac[1][Rho_mfi]),
#if (AMREX_SPACEDIM==3)
                                        BL_TO_FORTRAN_ANYD(u_mac[2][Rho_mfi]),
#endif
                                        &getForceVerbose);

               //
               // Average the new and old time to get Crank-Nicholson half time approximation.
               //
               FArrayBox Scal(amrex::grow(bx,0),NUM_SCALARS);
               Scal.copy<RunOn::Host>(S_old[Rho_mfi],bx,Density,bx,0,NUM_SCALARS);
               Scal.plus<RunOn::Host>(S_new[Rho_mfi],bx,Density,0,NUM_SCALARS);
               Scal.mult<RunOn::Host>(0.5,bx);

               if (getForceVerbose) amrex::Print() << "Calling getForce..." << '\n';
               tforces.resize(bx,1);
               getForce(tforces,bx,0,sigma,1,halftime,Vel,Scal,rhs[Rho_mfi],0);

	       const auto& Snew = S_new[Rho_mfi].array(sigma);
	       const auto& Sold = S_old[Rho_mfi].const_array(sigma);
	       const auto& aofs = Aofs[Rho_mfi].const_array(sigma);
	       const auto& tf   = tforces.const_array();
	       
	       amrex::ParallelFor(bx, [ Snew, Sold, aofs, tf, dt]
	       AMREX_GPU_DEVICE (int i, int j, int k ) noexcept
               {
		 Snew(i,j,k) = Sold(i,j,k) + dt * ( tf(i,j,k) -aofs(i,j,k) );
               });

               // Either need this synchronize here, or elixirs. Not sure if it matters which
	       amrex::Gpu::synchronize();
            }
        }
}
    }
    //
    // Call ScalMinMax to avoid overshoots in the scalars.
    //

    if ( do_scalminmax && (sComp <= last_scalar) )
    {
        const int num_scalars = last_scalar - Density + 1;
        //
        // Must do FillPatch here instead of MF iterator because we need the
        // boundary values in the old data (especially at inflow).
        //

        FillPatchIterator S_fpi(*this,S_old,1,prev_time,State_Type,Density,num_scalars);
        MultiFab& Smf=S_fpi.get_mf();

            for (int sigma = sComp; sigma <= last_scalar; sigma++)
            {
                const int index_new_s   = sigma;
                const int index_new_rho = Density;
                const int index_old_s   = index_new_s   - Density;
                const int index_old_rho = index_new_rho - Density;

                if (advectionType[sigma] == Conservative)
                {
                ConservativeScalMinMax(S_new, index_new_s, index_new_rho,
                                       Smf,   index_old_s, index_old_rho);
                }
                else if (advectionType[sigma] == NonConservative)
                {
                ConvectiveScalMinMax(S_new, index_new_s, Smf, index_old_s);
                }
            }

    }
}

//
// Set the time levels to time (time) and timestep dt.
//
void
NavierStokesBase::setTimeLevel (Real time,
                                Real dt_old,
                                Real dt_new)
{
    state[State_Type].setTimeLevel(time,dt_old,dt_new);

    if (have_divu)
    {
        state[Divu_Type].setTimeLevel(time,dt_old,dt_new);
        if (have_dsdt)
        {
            state[Dsdt_Type].setTimeLevel(time,dt_old,dt_new);
        }
    }

    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Interval)
    {
        state[Press_Type].setTimeLevel(time-dt_old,dt_old,dt_old);
    }
    else if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
    {
        state[Press_Type].setTimeLevel(time-.5*dt_old,dt_old,dt_old);
        state[Dpdt_Type].setTimeLevel(time-dt_old,dt_old,dt_old);
    }
}

void
NavierStokesBase::sync_setup (MultiFab*& DeltaSsync)
{
    BL_ASSERT(DeltaSsync == 0);

    int nconserved = 0;

    for (int comp = AMREX_SPACEDIM; comp < NUM_STATE; ++comp)
    {
        if (advectionType[comp] == Conservative)
            ++nconserved;
    }

    if (nconserved > 0 && level < parent->finestLevel())
    {
        DeltaSsync = new MultiFab(grids, dmap, nconserved, 1, MFInfo(), Factory());
        DeltaSsync->setVal(0,1);
    }
}

void
NavierStokesBase::sync_cleanup (MultiFab*& DeltaSsync)
{
    delete DeltaSsync;

    DeltaSsync = 0;
}

//
// Helper function for NavierStokesBase::SyncInterp().
//
static
void
set_bcrec_new (Vector<BCRec>  &bcrec,
               int             ncomp,
               int             src_comp,
               const Box&      box,
               const Box&      domain,
               const BoxArray& cgrids,
               int**           bc_orig_qty)

{
   for (int n = 0; n < ncomp; n++) {
      for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
      {
         int bc_index = (src_comp+n)*(2*AMREX_SPACEDIM) + dir; 
         bcrec[n].setLo(dir,INT_DIR);
         bcrec[n].setHi(dir,INT_DIR);
         if ( ( box.smallEnd(dir) < domain.smallEnd(dir) ) ||
              ( box.bigEnd(dir)   > domain.bigEnd(dir) ) ) {
            for (int crse = 0; crse < cgrids.size(); crse++) {
               const Box& crsebx = cgrids[crse];
               if ( ( box.smallEnd(dir) < domain.smallEnd(dir) ) && ( crsebx.smallEnd(dir) == domain.smallEnd(dir) ) ) {
                  bcrec[n].setLo(dir,bc_orig_qty[crse][bc_index]);
               }
               if ( ( box.bigEnd(dir) > domain.bigEnd(dir) ) && ( crsebx.bigEnd(dir) == domain.bigEnd(dir) ) ) {
                  bcrec[n].setHi(dir,bc_orig_qty[crse][bc_index+AMREX_SPACEDIM]);
               }
            }
         }
      }
   }
}

//
// Interpolate A cell centered Sync correction from a
// coarse level (c_lev) to a fine level (f_lev).
//
// This routine interpolates the num_comp components of CrseSync
// (starting at src_comp) and either increments or puts the result into
// the num_comp components of FineSync (starting at dest_comp)
// The components of bc_orig_qty corespond to the quantities of CrseSync.
//
void
NavierStokesBase::SyncInterp (MultiFab&      CrseSync,
                              int            c_lev,
                              MultiFab&      FineSync,
                              int            f_lev,
                              IntVect&       ratio,
                              int            src_comp,
                              int            dest_comp,
                              int            num_comp,
                              int            increment,
                              Real           dt_clev,
                              int**          bc_orig_qty,
                              SyncInterpType which_interp,
                              int            state_comp)
{
    BL_PROFILE("NavierStokesBase::SyncInterp()");

    BL_ASSERT(which_interp >= 0 && which_interp <= 5);

    Interpolater* interpolater = 0;

#ifdef AMREX_USE_EB
    switch (which_interp)
    {
       // As with the non-EB case, both of these point to the same interpolater
       case CellCons_T:     interpolater = &eb_cell_cons_interp;    break;
       case CellConsLin_T:  interpolater = &eb_lincc_interp;        break;
       default:
       amrex::Abort("NavierStokesBase::SyncInterp(): EB currently requires Cell Conservative interpolater. \n");
    }
#else
    switch (which_interp)
    {
       case PC_T:           interpolater = &pc_interp;           break;
       case CellCons_T:     interpolater = &cell_cons_interp;    break;
       case CellConsLin_T:  interpolater = &lincc_interp;        break;
       case CellConsProt_T: interpolater = &protected_interp;    break;
       default:
       amrex::Abort("NavierStokesBase::SyncInterp(): how did this happen \n");
    }
#endif

    NavierStokesBase& fine_level     = getLevel(f_lev);
    const BoxArray& fgrids           = fine_level.boxArray();
    const DistributionMapping& fdmap = fine_level.DistributionMap();
    const Geometry& fgeom            = parent->Geom(f_lev);
    const BoxArray& cgrids           = getLevel(c_lev).boxArray();
    const Geometry& cgeom            = parent->Geom(c_lev);
    Box             cdomain          = amrex::coarsen(fgeom.Domain(),ratio);
    const int       N                = fgrids.size();

    BoxArray cdataBA(N);

    for (int i = 0; i < N; i++) {
        cdataBA.set(i,interpolater->CoarseBox(fgrids[i],ratio));
    }
    //
    // Note: The boxes in cdataBA may NOT be disjoint !!!
    //
#ifdef AMREX_USE_EB
    // I am unsure of EBSupport and ng (set to zero here)
    auto factory = makeEBFabFactory(cgeom,cdataBA,fdmap,{0,0,0},EBSupport::basic);
    MultiFab cdataMF(cdataBA,fdmap,num_comp,0,MFInfo(),*factory);
#else
    //    ,MFInfo(),getLevel(c_lev).Factory());
    MultiFab cdataMF(cdataBA,fdmap,num_comp,0);
#endif

    cdataMF.copy(CrseSync, src_comp, 0, num_comp, cgeom.periodicity());

    //
    // Set physical boundary conditions in cdataMF.
    //
    //////////
    // Should be fine for EB for now, since EB doesn't intersect Phys BC
    // Not sure about what happens if EB intersects Phys BC
    ///////

    // tiling may not be needed here, but what the hey
    GpuBndryFuncFab<DummyFill> gpu_bndry_func(DummyFill{});
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(cdataMF,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       const Box& bx   = mfi.tilebox();
       FArrayBox& data = cdataMF[mfi];

       Vector<BCRec> bx_bcrec(num_comp);
       set_bcrec_new(bx_bcrec,num_comp,src_comp,bx,cdomain,cgrids,bc_orig_qty);
       gpu_bndry_func(bx,data,0,num_comp,cgeom,0.0,bx_bcrec,0,0);
    }

    //
    // Interpolate from cdataMF to fdata and update FineSync.
    // Note that FineSync and cdataMF will have the same distribution
    // since the length of their BoxArrays are equal.
    //
    MultiFab* fine_stateMF = 0;
    if (interpolater == &protected_interp)
    {
        fine_stateMF = &(getLevel(f_lev).get_new_data(State_Type));
    }


#ifdef AMREX_USE_EB
    //FIXME?
    // there's currently no way to reassign the EBCellFlagFab for a EBFArrayBox,
    // so the non-EB strategy of creating one FAB and resizing it for MFiters doens't work
    const FabArray<EBCellFlagFab>& flags = dynamic_cast<EBFArrayBoxFactory const&>(getLevel(f_lev).Factory()).getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
      for (MFIter mfi(FineSync,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         FArrayBox& cdata = cdataMF[mfi];
         const Box&  bx   = mfi.tilebox();
         const Box cbx    = interpolater->CoarseBox(bx,ratio);

#ifdef AMREX_USE_EB
         EBFArrayBox fdata(flags[mfi],bx,num_comp,FineSync[mfi].arena());
#else
         FArrayBox fdata(bx, num_comp);
#endif
         Elixir fdata_i = fdata.elixir();

         //
         // Set the boundary condition array for interpolation.
         //
         Vector<BCRec> bx_bcrec(num_comp);
         set_bcrec_new(bx_bcrec,num_comp,src_comp,cbx,cdomain,cgrids,bc_orig_qty);

         //ScaleCrseSyncInterp(cdata, c_lev, num_comp);

         interpolater->interp(cdata,0,fdata,0,num_comp,bx,ratio,
                              cgeom,fgeom,bx_bcrec,src_comp,State_Type,RunOn::Gpu);

         //reScaleFineSyncInterp(fdata, f_lev, num_comp);

         if (increment)
         {
            auto const& finedata    = fdata.array();
            auto const& coarsedata  = cdata.array();
            int scale_coarse = (interpolater == &protected_interp) ? 1 : 0;
            amrex::ParallelFor(bx, num_comp, [finedata,coarsedata,dt_clev, scale_coarse]
            AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
            {
               finedata(i,j,k,n) *= dt_clev;
               if ( scale_coarse ) {
                  coarsedata(i,j,k,n) *= dt_clev; 
               }
            });

            if (interpolater == &protected_interp)
            {
               FArrayBox& fine_state = (*fine_stateMF)[mfi];
               interpolater->protect(cdata,0,fdata,0,fine_state,state_comp,
                                     num_comp,bx,ratio,
                                     cgeom,fgeom,bx_bcrec,RunOn::Gpu);
            }

            auto const& fsync       = FineSync.array(mfi,dest_comp);
            amrex::ParallelFor(bx, num_comp, [finedata,fsync,coarsedata,dt_clev,scale_coarse]
            AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
            {
               if ( scale_coarse ) {
                  coarsedata(i,j,k,n) /= dt_clev;
               }
               fsync(i,j,k,n) += finedata(i,j,k,n);
            });
         }
         else
         {
            auto const& finedata    = fdata.array();
            auto const& fsync       = FineSync.array(mfi,dest_comp);
            amrex::ParallelFor(bx, num_comp, [finedata,fsync]
            AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
            {
               fsync(i,j,k,n) = finedata(i,j,k,n);
            });
         }
       }
    }
}

//
// Interpolate sync pressure correction to a finer level.
//
void
NavierStokesBase::SyncProjInterp (MultiFab& phi,
                                  int       c_lev,
                                  MultiFab& P_new,
                                  MultiFab& P_old,
                                  int       f_lev,
                                  IntVect&  ratio,
                                  bool      first_crse_step_after_initial_iters,
                                  Real      cur_crse_pres_time,
                                  Real      prev_crse_pres_time)
{
    BL_PROFILE("NavierStokesBase:::SyncProjInterp()");

    const BoxArray& P_grids = P_new.boxArray();
    const int       N       = P_grids.size();

    BoxArray crse_ba(N);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        crse_ba.set(i,node_bilinear_interp.CoarseBox(P_grids[i],ratio));

    // None  of these 3 are actually used by node_bilinear_interp()
    Vector<BCRec> bc(AMREX_SPACEDIM);
    const Geometry& fgeom   = parent->Geom(f_lev);
    const Geometry& cgeom   = parent->Geom(c_lev);

#ifdef AMREX_USE_EB
    // I am unsure of EBSupport and ng (set to 1 here)
    // need 1 ghost cell to use EB_set_covered on nodal MF
    // Factory is always CC, regardless of status of crse_ba
    auto factory = makeEBFabFactory(cgeom,crse_ba,P_new.DistributionMap(),{1,1,1},EBSupport::basic);
    MultiFab     crse_phi(crse_ba,P_new.DistributionMap(),1,0,MFInfo(),*factory);

#else
    MultiFab     crse_phi(crse_ba,P_new.DistributionMap(),1,0);
#endif

    crse_phi.setVal(1.e200);
    crse_phi.copy(phi,0,0,1);

#ifdef AMREX_USE_EB
    // FIXME -
    // For now, just zero covered fine cells. Better interpolation to come...
    ///
    EB_set_covered(crse_phi,0.);
#endif

    NavierStokesBase& fine_lev        = getLevel(f_lev);
    const Real    cur_fine_pres_time  = fine_lev.state[Press_Type].curTime();
    const Real    prev_fine_pres_time = fine_lev.state[Press_Type].prevTime();

    if (state[Press_Type].descriptor()->timeType() ==
        StateDescriptor::Point && first_crse_step_after_initial_iters)
    {
        const Real time_since_zero  = cur_crse_pres_time - prev_crse_pres_time;
        const Real dt_to_prev_time  = prev_fine_pres_time - prev_crse_pres_time;
        const Real dt_to_cur_time   = cur_fine_pres_time - prev_crse_pres_time;
        const Real cur_mult_factor  = dt_to_cur_time / time_since_zero;
        const Real prev_mult_factor = dt_to_prev_time / dt_to_cur_time;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(P_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
             const Box&  bx     = mfi.tilebox();
             FArrayBox fine_phi(bx,1);
             Elixir fine_phi_i = fine_phi.elixir();
             node_bilinear_interp.interp(crse_phi[mfi],0,fine_phi,0,1,
                                         fine_phi.box(),ratio,cgeom,fgeom,bc,
                                         0,Press_Type,RunOn::Gpu);

             auto const& f_phi    = fine_phi.array();
             auto const& p_new    = P_new.array(mfi);
             auto const& p_old    = P_old.array(mfi);
             amrex::ParallelFor(bx, [f_phi, p_old, p_new, cur_mult_factor, prev_mult_factor]
             AMREX_GPU_DEVICE(int i, int j, int k) noexcept
             {
                 p_new(i,j,k) += f_phi(i,j,k) * cur_mult_factor;
                 p_old(i,j,k) += f_phi(i,j,k) * prev_mult_factor;
             });
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(P_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
           const Box&  bx     = mfi.tilebox();
           FArrayBox fine_phi(bx,1);
           Elixir fine_phi_i = fine_phi.elixir();
           node_bilinear_interp.interp(crse_phi[mfi],0,fine_phi,0,1,
                                       fine_phi.box(),ratio,cgeom,fgeom,bc,
                                       0,Press_Type,RunOn::Gpu);

           auto const& f_phi    = fine_phi.array();
           auto const& p_new    = P_new.array(mfi);
           auto const& p_old    = P_old.array(mfi);
           amrex::ParallelFor(bx, [f_phi, p_old, p_new]
           AMREX_GPU_DEVICE(int i, int j, int k) noexcept
           {
               p_new(i,j,k) += f_phi(i,j,k);
               p_old(i,j,k) += f_phi(i,j,k);
           });
        }
    }
#ifdef AMREX_USE_EB
    // FIXME? - this can probably go after new interpolation is implemented
    EB_set_covered(P_new,0.);
    EB_set_covered(P_old,0.);
#endif

}

std::string
NavierStokesBase::thePlotFileType () const
{
    //
    // Increment this whenever the writePlotFile() format changes.
    //
    static const std::string the_plot_file_type("NavierStokes-V1.1");

    return the_plot_file_type;
}

//
// This routine advects the velocities
//
void
NavierStokesBase::velocity_advection (Real dt)
{
    BL_PROFILE("NavierStokesBase::velocity_advection()");

    if (verbose)
    {
        if (do_mom_diff == 0)
        {
            amrex::Print() << "... advect velocities\n";
        }
        else
        {
            if (predict_mom_together == 0)
            {
                amrex::Print() << "Must set predict_mom_together == 1 in NavierStokesBase." << '\n';
                exit(0);
            }
            amrex::Print() << "... advect momenta\n";
        }
    }

    const int   finest_level   = parent->finestLevel();
    const Real  prev_time      = state[State_Type].prevTime();

    //
    // Compute viscosity components.
    //
#ifdef AMREX_USE_EB
    MultiFab& Gp = getGradP();
    Gp.FillBoundary(geom.periodicity());
#else
    MultiFab Gp(grids,dmap,AMREX_SPACEDIM,1);
    getGradP(Gp, state[Press_Type].prevTime());
#endif

    MultiFab visc_terms(grids,dmap,AMREX_SPACEDIM,1,MFInfo(),Factory());

    // No need to compute this is we are using EB because we will
    // not use Godunov.
#ifndef AMREX_USE_EB
    if (be_cn_theta != 1.0)
    {
        getViscTerms(visc_terms,Xvel,AMREX_SPACEDIM,prev_time);
    }
    else
    {
        visc_terms.setVal(0.,1);
    }
#endif

    MultiFab divu_fp(grids,dmap,1,1,MFInfo(),Factory());
    create_mac_rhs(divu_fp,1,prev_time,dt);

    MultiFab fluxes[AMREX_SPACEDIM];

    if (do_reflux)
    {
        for (int i = 0; i < AMREX_SPACEDIM; i++)
        {
            const BoxArray& ba = getEdgeBoxArray(i);
            fluxes[i].define(ba, dmap, AMREX_SPACEDIM, 0, MFInfo(),Factory());
        }
    }

    //
    // Compute the advective forcing.
    //
    {
        FillPatchIterator U_fpi(*this,visc_terms,godunov_hyp_grow,prev_time,State_Type,Xvel,AMREX_SPACEDIM);
        MultiFab& Umf=U_fpi.get_mf();

#ifndef AMREX_USE_EB
        //
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>  NON-EB ALGORITHM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        //
        FillPatchIterator Rho_fpi(*this,visc_terms,godunov_hyp_grow,prev_time,State_Type,Density,1);
        FillPatchIterator S_fpi(*this,visc_terms,1,prev_time,State_Type,Density,NUM_SCALARS);
        MultiFab& Rmf=Rho_fpi.get_mf();
        MultiFab& Smf=S_fpi.get_mf();

        const Real* dx = geom.CellSize();

        MultiFab cfluxes[AMREX_SPACEDIM];
        MultiFab edgestate[AMREX_SPACEDIM];
        int ngrow = 1;

        MultiFab forcing_term( grids, dmap, AMREX_SPACEDIM, ngrow );
        MultiFab S_term( grids, dmap, AMREX_SPACEDIM,  godunov_hyp_grow);

        // Why in the original code it does this:
        //
        //
        //
        int nghost = 2;  // Do we need 2???
        for (int i = 0; i < AMREX_SPACEDIM; i++)
        {
            const BoxArray& ba = getEdgeBoxArray(i);
            cfluxes[i].define(ba, dmap, AMREX_SPACEDIM, nghost, MFInfo(), Umf.Factory());
            cfluxes[i].setVal(0.0);
            edgestate[i].define(ba, dmap, AMREX_SPACEDIM, nghost, MFInfo(), Umf.Factory());
        }


	
        //
        // Compute forcing
        //

           const BoxArray& ba = Umf.boxArray();
           const DistributionMapping& dm = Umf.DistributionMap();
	   int ncomp = Umf.nComp();
	   int ngrow_r = Umf.nGrow();
	   MultiFab rhs(ba,dm,ncomp,ngrow_r);
	   rhs.setVal(0.0);


	   /*
      amrex::AllPrint() << " xxxxx proc. = " << ParallelDescriptor::MyProc() << ", tag = " <<
      ParallelDescriptor::SeqNum()
		      << " file = " << __FILE__ << " function = " << __FUNCTION__
          	      << " line = " << __LINE__ << std::endl;
	   */
      

#ifdef AMREX_PARTICLES	   
           for (MFIter U_mfi(Umf,true); U_mfi.isValid(); ++U_mfi)
           {
 	      if (NavierStokes::initial_iter != true) {	     
	      theNSPC()->getDrag(rhs[U_mfi],Umf[U_mfi],Smf[U_mfi],visc_coef[0],1,level);
	      }
	   }
	   rhs.SumBoundary(Geom().periodicity()); 	   	   
#endif


	
	
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            for (MFIter U_mfi(Umf,TilingIfNotGPU()); U_mfi.isValid(); ++U_mfi)
            {
                const Box&  bx = U_mfi.tilebox();
                auto const gbx = U_mfi.growntilebox(ngrow);


                if (getForceVerbose)
                {
                    amrex::Print() << "-----------------------\n"
                                   << "B - velocity advection:\n"
                                   << "-----------------------\n";
                }
	        getForce(forcing_term[U_mfi],gbx,ngrow,Xvel,AMREX_SPACEDIM,prev_time,Umf[U_mfi],Smf[U_mfi],rhs[U_mfi],0);

                //
                // Compute the total forcing.
                //
                auto const& tf   = forcing_term.array(U_mfi,Xvel);
                auto const& visc = visc_terms.const_array(U_mfi,Xvel);
                auto const& gp   = Gp.const_array(U_mfi);
                auto const& rho  = rho_ptime.const_array(U_mfi);

                amrex::ParallelFor(gbx, AMREX_SPACEDIM, [tf, visc, gp, rho]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    tf(i,j,k,n) = ( tf(i,j,k,n) + visc(i,j,k,n) - gp(i,j,k,n) ) / rho(i,j,k);
                });

                //
                // Loop over the velocity components.
                //
                auto const hypbox = U_mfi.growntilebox(godunov_hyp_grow);
                S_term[U_mfi].copy<RunOn::Host>(Umf[U_mfi],0,0,AMREX_SPACEDIM);

                for (int comp = 0; comp < AMREX_SPACEDIM; ++comp )
                {
                    int use_conserv_diff = (advectionType[comp] == Conservative) ? true : false;

                    if (do_mom_diff == 1)
                    {
                        S_term[U_mfi].mult<RunOn::Host>(Rmf[U_mfi],hypbox,hypbox,0,comp,1);
                        forcing_term[U_mfi].mult<RunOn::Host>(rho_ptime[U_mfi],gbx,gbx,0,comp,1);
                    }
		}

            }
        }


        Vector<BCRec> math_bcs(AMREX_SPACEDIM);
        math_bcs = fetchBCArray(State_Type, Xvel, AMREX_SPACEDIM);

        amrex::Gpu::DeviceVector<int> iconserv;
        iconserv.resize(AMREX_SPACEDIM, 0);

        for (int comp = 0; comp < AMREX_SPACEDIM; ++comp )
        {
            iconserv[comp] = (advectionType[comp] == Conservative) ? true : false;
        }

        Godunov::ComputeAofs(*aofs, Xvel, AMREX_SPACEDIM,
                             S_term, 0,
                             AMREX_D_DECL(u_mac[0],u_mac[1],u_mac[2]),
                             AMREX_D_DECL(edgestate[0],edgestate[1],edgestate[2]), 0, false,
                             AMREX_D_DECL(cfluxes[0],cfluxes[1],cfluxes[2]), 0,
                             forcing_term, 0, divu_fp, math_bcs, geom, iconserv, dt,
                             godunov_use_ppm, godunov_use_forces_in_trans, true);

	if (do_reflux)
	{
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                MultiFab::Copy(fluxes[d], cfluxes[d], 0, 0, AMREX_SPACEDIM, 0 );
	}




#else
        //
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>  EB ALGORITHM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        //

         MultiFab cfluxes[AMREX_SPACEDIM];
         MultiFab edgstate[AMREX_SPACEDIM];
         int nghost(2);

         for (int i(0); i < AMREX_SPACEDIM; i++)
         {
             const BoxArray& ba = getEdgeBoxArray(i);
             cfluxes[i].define(ba, dmap, AMREX_SPACEDIM, nghost, MFInfo(), Umf.Factory());
             edgstate[i].define(ba, dmap, AMREX_SPACEDIM, nghost, MFInfo(), Umf.Factory());
         }

         Vector<BCRec> math_bcs(AMREX_SPACEDIM);
         math_bcs = fetchBCArray(State_Type, Xvel, AMREX_SPACEDIM);

         MOL::ComputeAofs(*aofs, Xvel, AMREX_SPACEDIM, Umf, 0,
                          D_DECL(u_mac[0],u_mac[1],u_mac[2]),
                          D_DECL(edgstate[0],edgstate[1],edgstate[2]), 0, false,
                          D_DECL(cfluxes[0],cfluxes[1],cfluxes[2]), 0,
                          math_bcs, geom  );

         // don't think this is needed here any more. Godunov sets covered vals now...
         EB_set_covered(*aofs, 0.);

         if (do_reflux)
         {
             for (int d(0); d < AMREX_SPACEDIM; d++)
                 MultiFab::Copy(fluxes[d], cfluxes[d], 0, 0, AMREX_SPACEDIM, 0 );

         }

#endif
    } //end scope of FillPatchIter

    if (do_reflux)
    {
        if (level > 0 )
        {
            for (int d = 0; d < AMREX_SPACEDIM; d++)
                advflux_reg->FineAdd(fluxes[d],d,0,0,AMREX_SPACEDIM,dt);
        }
        if(level < finest_level)
        {
            for (int i = 0; i < AMREX_SPACEDIM; i++)
                getAdvFluxReg(level+1).CrseInit(fluxes[i],i,0,0,AMREX_SPACEDIM,-dt);
        }
    }
}

//
// This subroutine updates the velocity field before the level projection.
//
// At this point in time, all we know is u^n, rho^n+1/2, and the
// general forcing terms at t^n, and after solving in this routine
// viscous forcing at t^n+1/2.  Except for a simple buoyancy term,
// b = -rho^n+1/2 g, it is usually not possible to estimate more
// general forcing terms at t^n+1/2.  Since the default getForce, handles
// this case automatically, F_new and F_old have been replaced by a single
// tforces FArrayBox.
//
// We assume that if one component of velocity is viscous that all must be.
//

void
NavierStokesBase::velocity_update (Real dt)
{
    BL_PROFILE("NavierStokesBase::velocity_update()");

    if (verbose)
    {
      if (do_mom_diff == 0)
      {
         amrex::Print() << "... update velocities \n";
      }
      else
      {
         amrex::Print() << "... update momenta \n";
      }
    }

    velocity_advection_update(dt);

    if (!initial_iter)
        velocity_diffusion_update(dt);
    else
        initial_velocity_diffusion_update(dt);

    MultiFab&  S_new     = get_new_data(State_Type);

    for (int sigma = 0; sigma < BL_SPACEDIM; sigma++)
    {
       if (S_new.contains_nan(sigma,1,0))
       {
          amrex::Print() << "New velocity " << sigma << " contains Nans" << '\n';
          exit(0);
       }
    }
}

void
NavierStokesBase::velocity_advection_update (Real dt)
{
    BL_PROFILE("NavierStokesBase::velocity_advection_update()");

    MultiFab&  U_old          = get_old_data(State_Type);
    MultiFab&  U_new          = get_new_data(State_Type);
    MultiFab&  Aofs           = *aofs;
    const Real prev_pres_time = state[Press_Type].prevTime();

#ifdef AMREX_USE_EB
    MultiFab& Gp=*gradp;
    Gp.FillBoundary(geom.periodicity());
#else
    MultiFab Gp(grids,dmap,AMREX_SPACEDIM,1);
    getGradP(Gp, prev_pres_time);
#endif

    MultiFab& Rh = get_rho_half_time();




    // particle bit
    const BoxArray& ba = U_old.boxArray();
    const DistributionMapping& dm = U_old.DistributionMap();
    int ncomp = U_old.nComp();
    int ngrow = U_old.nGrow();
    MultiFab rhs(ba,dm,ncomp,ngrow);
    rhs.setVal(0.0);

    
#ifdef AMREX_PARTICLES

  FArrayBox  tforces, S, VelFAB, ScalFAB, rhsFAB;
    for (MFIter mfi(Rh,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        VelFAB.resize(bx,AMREX_SPACEDIM);
        ScalFAB.resize(bx,NUM_SCALARS);
        //rhsFAB.resize(bx,AMREX_SPACEDIM);	

        auto const& vel  = VelFAB.array();
        auto const& scal = ScalFAB.array();
        Elixir vel_i = VelFAB.elixir();
        Elixir scal_i = ScalFAB.elixir();
        D_TERM(auto const& umac = u_mac[0].array(mfi);,
               auto const& vmac = u_mac[1].array(mfi);,
               auto const& wmac = u_mac[2].array(mfi););
        auto const& scal_o = U_old.array(mfi,Density);
        auto const& scal_n = U_new.array(mfi,Density);
        const int numscal = NUM_SCALARS;
        amrex::ParallelFor(bx, [vel,D_DECL(umac,vmac,wmac),scal, scal_o, scal_n, numscal]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
           edg2cen_average(i,j,k,D_DECL(umac,vmac,wmac),vel);
           for (int n = 0; n < numscal; n++) {
              scal(i,j,k,n) = 0.5 * ( scal_o(i,j,k,n) + scal_n(i,j,k,n) );
           }
        });

        if (getForceVerbose) amrex::Print() << "Calling getForce..." << '\n';
        const Real half_time = 0.5*(state[State_Type].prevTime()+state[State_Type].curTime());

       if (NavierStokes::initial_iter != true) {
       theNSPC()->getDrag(rhs[mfi],VelFAB,ScalFAB,visc_coef[0],1,level);
       }
    }
    rhs.SumBoundary(Geom().periodicity());    
#endif


    


    
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
{
    FArrayBox  tforces, S, VelFAB, ScalFAB;
    for (MFIter mfi(Rh,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        VelFAB.resize(bx,AMREX_SPACEDIM);
        ScalFAB.resize(bx,NUM_SCALARS);

        //
        // Need to do some funky half-time stuff.
        //
        if (getForceVerbose) {
           amrex::Print() << "------------------------------------------" << '\n'
                          << "F - velocity advection update (half time):" << '\n'
                          << "------------------------------------------" << '\n';
	}	
        //
        // Average the mac face velocities to get cell centred velocities.
        // Average the new and old time to get Crank-Nicholson half time approximation.
        //
        //FIXME - need to address this for EB
        auto const& vel  = VelFAB.array();
        auto const& scal = ScalFAB.array();
        Elixir vel_i = VelFAB.elixir();
        Elixir scal_i = ScalFAB.elixir();
        D_TERM(auto const& umac = u_mac[0].array(mfi);,
               auto const& vmac = u_mac[1].array(mfi);,
               auto const& wmac = u_mac[2].array(mfi););
        auto const& scal_o = U_old.array(mfi,Density);
        auto const& scal_n = U_new.array(mfi,Density);
        const int numscal = NUM_SCALARS;
        amrex::ParallelFor(bx, [vel,D_DECL(umac,vmac,wmac),scal, scal_o, scal_n, numscal]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
           edg2cen_average(i,j,k,D_DECL(umac,vmac,wmac),vel);
           for (int n = 0; n < numscal; n++) {
              scal(i,j,k,n) = 0.5 * ( scal_o(i,j,k,n) + scal_n(i,j,k,n) );
           }
        });

        if (getForceVerbose) amrex::Print() << "Calling getForce..." << '\n';
        const Real half_time = 0.5*(state[State_Type].prevTime()+state[State_Type].curTime());
        tforces.resize(bx,AMREX_SPACEDIM);
        Elixir tf_i = tforces.elixir();
        getForce(tforces,bx,0,Xvel,AMREX_SPACEDIM,half_time,VelFAB,ScalFAB,rhs[mfi],0);

        //
        // Do following only at initial iteration--per JBB.
        //
        if (initial_iter && is_diffusive[Xvel]) {
           auto const& force  = tforces.array();
           amrex::ParallelFor(bx, [force]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {
               force(i,j,k) = 0.0;
           });
        }

        // Update velocity
        auto const& vel_old  = U_old.array(mfi);
        auto const& vel_new  = U_new.array(mfi);
        auto const& gradp    = Gp.array(mfi);
        auto const& force    = tforces.array();
        auto const& advec    = Aofs.array(mfi);
        auto const& rho_old  = rho_ptime.array(mfi);
        auto const& rho_new  = rho_ctime.array(mfi);
        auto const& rho_Half = Rh.array(mfi);
        int mom_diff = do_mom_diff;
        amrex::ParallelFor(bx, AMREX_SPACEDIM, [vel_old,vel_new,gradp,force,advec,rho_old,rho_new,rho_Half,mom_diff,dt]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real velold = vel_old(i,j,k,n);

            if ( mom_diff ) {
               velold *= rho_old(i,j,k);
	       vel_new(i,j,k,n) = velold - dt * advec(i,j,k,n)
		                         + dt * force(i,j,k,n)
		                         - dt * gradp(i,j,k,n);

	       vel_new(i,j,k,n) /= rho_new(i,j,k);
            }
	    else
	    {
	        vel_new(i,j,k,n) = velold - dt * advec(i,j,k,n)
                                          + dt * force(i,j,k,n) / rho_Half(i,j,k)
                                          - dt * gradp(i,j,k,n) / rho_Half(i,j,k);
            }
        });
    }
}

    for (int sigma = 0; sigma < AMREX_SPACEDIM; sigma++)
    {
       if (U_old.contains_nan(sigma,1,0))
       {
          amrex::Print() << "VAU: Old velocity " << sigma << " contains Nans" << '\n';
       }
       if (U_new.contains_nan(sigma,1,0))
       {
          amrex::Print() << "VAU: New velocity " << sigma << " contains Nans" << '\n';
       }
    }
}

void
NavierStokesBase::initial_velocity_diffusion_update (Real dt)
{
    //
    // Do following only at initial iteration.
    //
    if (is_diffusive[Xvel])
    {
        MultiFab&  U_old          = get_old_data(State_Type);
        MultiFab&  U_new          = get_new_data(State_Type);
        MultiFab&  Rh             = get_rho_half_time();
        const Real prev_time      = state[State_Type].prevTime();
        const int  xvel           = Xvel;

        int   ngrow = 0;
        MultiFab visc_terms(grids,dmap,AMREX_SPACEDIM,ngrow,MFInfo(),Factory());
        MultiFab    tforces(grids,dmap,AMREX_SPACEDIM,ngrow,MFInfo(),Factory());

        //
        // Get grad(p)
        //
#ifdef AMREX_USE_EB
        MultiFab& Gp = *gradp;
#else
        const Real prev_pres_time = state[Press_Type].prevTime();
        MultiFab         Gp(grids,dmap,AMREX_SPACEDIM,ngrow,MFInfo(),Factory());
        getGradP(Gp, prev_pres_time);
#endif

        //
        // Compute additional forcing terms
        //

        // Get particle forces
        const BoxArray& ba = U_new.boxArray();
        const DistributionMapping& dm = U_new.DistributionMap();
	int ncomp = U_new.nComp();
        ngrow = U_new.nGrow();
	MultiFab rhs(ba,dm,ncomp,ngrow);
	rhs.setVal(0.0);

#ifdef AMREX_PARTICLES	  	   
        for (MFIter mfi(tforces,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
           if (NavierStokes::initial_iter != true) {	     
	   theNSPC()->getDrag(rhs[mfi],U_old[mfi],U_old[mfi],visc_coef[0],1,level);
	   }
        }
        rhs.SumBoundary(Geom().periodicity());	
#endif
	
        tforces.setVal(0.0);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(tforces,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const auto& bx = mfi.tilebox();
                  auto& tforces_fab = tforces[mfi];

            if (getForceVerbose)
            {
               amrex::Print() << "--------------------------------------\n"
                              << "G - initial velocity diffusion update:\n"
                              << "--------------------------------------\n";
            }	    
            getForce(tforces_fab,bx,0,Xvel,AMREX_SPACEDIM,prev_time,U_old[mfi],U_old[mfi],rhs[mfi],Density);
        }

        //
        // Compute viscous terms
        //
        if (be_cn_theta != 1.0)
        {
          getViscTerms(visc_terms,Xvel,AMREX_SPACEDIM,prev_time);
        }
        else
        {
          visc_terms.setVal(0.0);
        }

        //
        // Assemble RHS
        //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(tforces,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
           const Box& bx = mfi.tilebox();
           auto const& force   = tforces.array(mfi);
           auto const& viscT   = visc_terms.array(mfi);
           auto const& gradp   = Gp.array(mfi);
           auto const& rhohalf = Rh.array(mfi);
           auto const& rho_old = rho_ptime.array(mfi);
           auto const& rho_new = rho_ctime.array(mfi);
           auto const& vel_old = U_old.array(mfi,xvel);
           auto const& vel_new = U_new.array(mfi,xvel);
           auto const& advT    = aofs->array(mfi,xvel);
           int mom_diff = do_mom_diff;
           amrex::ParallelFor(bx, AMREX_SPACEDIM, [force,viscT,gradp,rhohalf,advT,rho_old,rho_new,vel_old,vel_new,mom_diff,dt]
           AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
              // Set force += (visc - Gp) / rho_half - aofs
              force(i,j,k,n) += viscT(i,j,k,n) - gradp(i,j,k,n);
              if ( !mom_diff ) {
                 force(i,j,k,n) /= rhohalf(i,j,k);
              }
              force(i,j,k,n) -= advT(i,j,k,n);
              // if mom_diff : U_new = (force* dt + U_old * rho_old) / rho_new
              // else        : U_new = U_old + force* dt
              if ( mom_diff ) {
                 vel_new(i,j,k,n) = (force(i,j,k,n) * dt + vel_old(i,j,k,n) * rho_old(i,j,k)) / rho_new(i,j,k);
              } else {
                 vel_new(i,j,k,n) = vel_old(i,j,k,n) + force(i,j,k,n) * dt;
              }
           });
        }
    }
}

Real
NavierStokesBase::volWgtSum (const std::string& name,
                             Real               time)
{
    Real  volwgtsum = 0.0;
    const Real* dx  = geom.CellSize();
    auto        mf  = derive(name,time,0);

    // First, zero covered regions
    if (level < parent->finestLevel())
    {
        BoxArray    baf;
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
{
	std::vector< std::pair<int,Box> > isects;
        for (MFIter mfi(*mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
           auto const& fabarr = mf->array(mfi);
           int          ncomp = mf->nComp(); 
           baf.intersections(grids[mfi.index()],isects);
           for (int is = 0; is < isects.size(); is++) {
              amrex::ParallelFor(isects[is].second, ncomp, [fabarr]
              AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
              {
                 fabarr(i,j,k,n) = 0.0;
              });
           }
        }
}
    }

    // Use amrex::ReduceSum
    // TODO set the cases for RZ, but it needs the radius to be somewhere in managed memory and so on
    Real vol = D_TERM(dx[0],*dx[1],*dx[2]);
#ifdef AMREX_USE_EB
    if ( (volWgtSum_sub_dz > 0 && volWgtSum_sub_Rcyl > 0) ) {
        amrex::Abort("EB volWgtSum currently only works over entire cartesian domain.");
    }
    Real sm = amrex::ReduceSum(*mf, *volfrac, 0, [vol]
    AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& mf_arr, Array4<Real const> const& vf_arr) -> Real
    {
        Real sum = 0.0;
        AMREX_LOOP_3D(bx, i, j, k,
        {
            sum += mf_arr(i,j,k) * vf_arr(i,j,k) * vol;
        });
        return sum;
    });
#else
    const Real* dom_lo = geom.ProbLo();
    const Real sub_dz = volWgtSum_sub_dz;
    const Real sub_Rcyl = volWgtSum_sub_Rcyl;
    Real sm = amrex::ReduceSum(*mf, 0, [vol, sub_dz, sub_Rcyl, dx, dom_lo]
    AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& mf_arr) -> Real
    {
        Real sum = 0.0;
        if ( sub_dz > 0 && sub_Rcyl > 0 ) {
           // TODO : test this in 2D
           const auto lo = amrex::lbound(bx);
           const auto hi = amrex::ubound(bx);
           for       (int k = lo.z; k <= hi.z; ++k) {
              Real z = dom_lo[2] + (k+0.5_rt) * dx[2];
              if ( z <= sub_dz ) {
                 for    (int j = lo.y; j <= hi.y; ++j) {
                    Real y = dom_lo[1] + (j+0.5_rt) * dx[1];
                    for (int i = lo.x; i <= hi.x; ++i) {
                       Real x = dom_lo[0] + (i+0.5_rt) * dx[0];
                       Real r = std::sqrt(x*x + y*y);
                       if ( r <= sub_Rcyl ) {
                          sum += mf_arr(i,j,k) * vol;
                       }
                    } 
                 }
              } 
           }
        } else { 
           AMREX_LOOP_3D(bx, i, j, k,
           {
               sum += mf_arr(i,j,k) * vol;
           });
        }
        return sum;
    });
#endif

    volwgtsum = sm;

    ParallelDescriptor::ReduceRealSum(volwgtsum);

    return volwgtsum;
}

#if (AMREX_SPACEDIM == 3)
void
NavierStokesBase::sum_turbulent_quantities ()
{
    Real time = state[State_Type].curTime();
    const int finestLevel = parent->finestLevel();
    const Real *dx = parent->Geom(finestLevel).CellSize();
    const int ksize(parent->Geom(finestLevel).Domain().length(2));
    const int turbVars(33);
    int refRatio(1);

    Real* turb = new Real[turbVars*ksize];

    for (int i=0; i<turbVars*ksize; i++) turb[i]=0;

    for (int lev = finestLevel; lev >= 0; lev--)
    {
	const int levKsize(parent->Geom(lev).Domain().length(2));

	Real* levTurb = new Real[turbVars*levKsize];

	for (int i=0; i<turbVars*levKsize; i++) levTurb[i]=0;

        NavierStokesBase& ns_level = getLevel(lev);
	ns_level.TurbSum(time,levTurb,levKsize,turbVars);

	if (lev<finestLevel)  refRatio *= parent->refRatio(lev)[2];
	else                  refRatio  = 1;

	for (int l=0, k=0; l<levKsize; l++)
	    for (int r=0; r<refRatio; r++, k++)
		for (int v=0; v<turbVars; v++)
		    turb[k*turbVars+v] += levTurb[l*turbVars+v];

	delete [] levTurb;
    }

    ParallelDescriptor::ReduceRealSum(&turb[0], ksize*turbVars, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
    {
        std::string DirPath = "TurbData";
        if (!amrex::UtilCreateDirectory(DirPath, 0755))
            amrex::CreateDirectoryFailed(DirPath);

        const int steps = parent->levelSteps(0);
        FILE *file;

        std::string filename = amrex::Concatenate("TurbData/TurbData_", steps, 4);
        filename += ".dat";

        file = fopen(filename.c_str(),"w");
        for (int k=0; k<ksize; k++)
        {
            fprintf(file,"%e ",dx[2]*(0.5+(double)k));
            for (int v=0; v<turbVars; v++)
                fprintf(file,"%e ",turb[k*turbVars+v]);
            fprintf(file,"\n");
        }
        fclose(file);
    }

    delete [] turb;
}

void
NavierStokesBase::TurbSum (Real time, Real *turb, int ksize, int turbVars)
{
    const Real* dx = geom.CellSize();

    const int turbGrow(0);
    const int presGrow(0);
    auto turbMF = derive("TurbVars",time,turbGrow);
    auto presMF = derive("PresVars",time,presGrow);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    std::vector< std::pair<int,Box> > isects;

    for (MFIter turbMfi(*turbMF), presMfi(*presMF);
	 turbMfi.isValid() && presMfi.isValid();
	 ++turbMfi, ++presMfi)
    {
	FArrayBox& turbFab = (*turbMF)[turbMfi];
	FArrayBox& presFab = (*presMF)[presMfi];

        if (level < parent->finestLevel())
        {
            baf.intersections(grids[turbMfi.index()],isects);

            for (int ii = 0, N = isects.size(); ii < N; ii++)
            {
                presFab.setVal<RunOn::Host>(0,isects[ii].second,0,presMF->nComp());
                turbFab.setVal<RunOn::Host>(0,isects[ii].second,0,turbMF->nComp());
            }
        }
    }

    turbMF->FillBoundary(0,turbMF->nComp(), geom.periodicity());
    presMF->FillBoundary(0,presMF->nComp(), geom.periodicity());

    for (MFIter turbMfi(*turbMF), presMfi(*presMF);
	 turbMfi.isValid() && presMfi.isValid();
	 ++turbMfi, ++presMfi)
    {
	FArrayBox& turbFab = (*turbMF)[turbMfi];
	FArrayBox& presFab = (*presMF)[presMfi];

        const Real* turbData = turbFab.dataPtr();
        const Real* presData = presFab.dataPtr();
        const int*  dlo = turbFab.loVect();
        const int*  dhi = turbFab.hiVect();
        const int*  plo = presFab.loVect();
        const int*  phi = presFab.hiVect();
	const Box& grdbx = grids[turbMfi.index()];
        const int*  lo  = grdbx.loVect();
        const int*  hi  = grdbx.hiVect();

        sumturb(turbData,presData,ARLIM(dlo),ARLIM(dhi),ARLIM(plo),ARLIM(phi),ARLIM(lo),ARLIM(hi),
		     dx,turb,&ksize,&turbVars);
   }
}

#ifdef SUMJET
void
NavierStokesBase::JetSum (Real time, Real *jetData, int levRsize,  int levKsize,  int rsize,  int ksize, int jetVars)
{
    const Real* dx = geom.CellSize();

    const int turbGrow(0);
    const int presGrow(0);

    auto turbMF = derive("JetVars",time,turbGrow);
    auto presMF = derive("JetPresVars",time,presGrow);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    std::vector< std::pair<int,Box> > isects;

    for (MFIter turbMfi(*turbMF), presMfi(*presMF);
	 turbMfi.isValid() && presMfi.isValid();
	 ++turbMfi, ++presMfi)
    {
	FArrayBox& turbFab = (*turbMF)[turbMfi];
	FArrayBox& presFab = (*presMF)[presMfi];

        if (level < parent->finestLevel())
        {
            baf.intersections(grids[turbMfi.index()],isects);

            for (int ii = 0, N = isects.size(); ii < N; ii++)
            {
                presFab.setVal(0,isects[ii].second,0,presMF->nComp());
                turbFab.setVal(0,isects[ii].second,0,turbMF->nComp());
            }
        }
    }

    turbMF->FillBoundary(0,turbMF->nComp(), geom.periodicity());
    presMF->FillBoundary(0,presMF->nComp(), geom.periodicity());

    for (MFIter turbMfi(*turbMF), presMfi(*presMF);
	 turbMfi.isValid() && presMfi.isValid();
	 ++turbMfi, ++presMfi)
    {
	FArrayBox& turbFab = (*turbMF)[turbMfi];
	FArrayBox& presFab = (*presMF)[presMfi];

        RealBox     gridloc  = RealBox(grids[turbMfi.index()],geom.CellSize(),geom.ProbLo());
        const Real* turbData = turbFab.dataPtr();
        const Real* presData = presFab.dataPtr();
        const int*  dlo = turbFab.loVect();
        const int*  dhi = turbFab.hiVect();
        const int*  plo = presFab.loVect();
        const int*  phi = presFab.hiVect();
        const int*  lo  = grids[turbMfi.index()].loVect();
        const int*  hi  = grids[turbMfi.index()].hiVect();

        sumjet(turbData,presData,ARLIM(dlo),ARLIM(dhi),ARLIM(plo),ARLIM(phi),ARLIM(lo),ARLIM(hi),
		    dx,jetData,&levRsize,&levKsize,&rsize,&ksize,&jetVars,&jet_interval_split,
		    gridloc.lo(),gridloc.hi());
    }
}

void
NavierStokesBase::sum_jet_quantities ()
{
    Real time = state[State_Type].curTime();
    const int finestLevel = parent->finestLevel();
    const Real *dx = parent->Geom(finestLevel).CellSize();
    const int isize(parent->Geom(finestLevel).Domain().length(0));
    const int ksize(parent->Geom(finestLevel).Domain().length(2));
    const int rsize=isize>>1;
    const int jetVars(104);

    amrex::Print() << "NavierStokesBase::sum_jet_quantities():" << '\n'
		   << "   jetVars: " << jetVars << '\n'
		   << "   rsize  : " << rsize << '\n'
		   << "   ksize  : " << ksize << '\n';

    Real* jetData = new Real[jetVars*ksize*rsize];

    for (int i=0; i<jetVars*ksize*rsize; i++) jetData[i]=0;

    for (int lev = finestLevel; lev >= 0; lev--)
    {
	const int levIsize(parent->Geom(lev).Domain().length(0));
	const int levKsize(parent->Geom(lev).Domain().length(2));
	const int levRsize(levIsize>>1);

        NavierStokesBase& ns_level = getLevel(lev);
	ns_level.JetSum(time,jetData,levRsize,levKsize,rsize,ksize,jetVars);
    }

    ParallelDescriptor::ReduceRealSum(&jetData[0], ksize*rsize*jetVars, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << "      Creating JetData..." << '\n';
        std::string DirPath = "JetData";
        if (!amrex::UtilCreateDirectory(DirPath, 0755))
            amrex::CreateDirectoryFailed(DirPath);

        const int steps = parent->levelSteps(0);
        FILE *file;
        std::string filename;

	Vector<Real> r(rsize);
	for (int i=0; i<rsize; i++)
	    r[i] = dx[0]*(0.5+(double)i);
	Vector<Real> z(ksize);
	for (int k=0; k<ksize; k++)
	    z[k] = dx[2]*(0.5+(double)k);

#if 0
        filename  = amrex::Concatenate("JetData/JetData_", steps, 4);
        filename += "_r.dat";

	file = fopen(filename.c_str(),"w");
	for (int i=0; i<rsize; i++)
	    fprintf(file,"%e ",r[i]);
	fclose(file);

        filename  = amrex::Concatenate("JetData/JetData_", steps, 4);
        filename += "_z.dat";

	file = fopen(filename.c_str(),"w");
	for (int k=0; k<ksize; k++)
	    fprintf(file,"%e ",dx[2]*(0.5+(double)k));
	fclose(file);

	for (int v=0; v<jetVars; v++) {

            filename  = amrex::Concatenate("JetData/JetData_", steps, 4);
            filename += amrex::Concatenate(filename + "_v", v, 4);
            filename += ".dat";

	    file = fopen(filename.c_str(),"w");
	    for (int k=0; k<ksize; k++) {
		for (int i=0; i<rsize; i++) {
		    fprintf(file,"%e ",jetData[(k*rsize+i)*jetVars+v]);
		}
		fprintf(file,"\n");
	    }
	    fclose(file);
	    amrex::Print() << "   ...done." << '\n';
	}
#else
	std::string FullPath = amrex::Concatenate("JetData/JD", steps, 4);

	if (!amrex::UtilCreateDirectory(FullPath, 0755))
	    amrex::CreateDirectoryFailed(FullPath);

        filename = FullPath;
        filename += '/';
        filename += "data.bin";

	file=fopen(filename.c_str(),"w");
	fwrite(&time,sizeof(double),1,file);
	fwrite(&rsize,sizeof(int),1,file);
	fwrite(&ksize,sizeof(int),1,file);
	fwrite(&jetVars,sizeof(int),1,file);
	fwrite(r.dataPtr(),sizeof(Real),rsize,file);
	fwrite(z.dataPtr(),sizeof(Real),ksize,file);
	fwrite(jetData,sizeof(Real),jetVars*rsize*ksize,file);
	fclose(file);
#endif
    }

    delete [] jetData;
}
#endif // SUMJET

#endif  // (BL_SPACEDIM == 3)

#ifdef AMREX_PARTICLES

void
NavierStokesBase::read_particle_params ()
{
    ParmParse ppp("particles");
    //
    // The directory in which to store timestamp files.
    //
    ppp.query("timestamp_dir", timestamp_dir);
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!amrex::UtilCreateDirectory(timestamp_dir, 0755))
            amrex::CreateDirectoryFailed(timestamp_dir);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (int nc = ppp.countval("timestamp_indices"))
    {
        timestamp_indices.resize(nc);

        ppp.getarr("timestamp_indices", timestamp_indices, 0, nc);
    }

    ppp.query("pverbose",pverbose);
    //
    // Used in initData() on startup to read in a file of particles.
    //
    ppp.query("particle_init_file", particle_init_file);
    //
    // Used in post_restart() to read in a file of particles.
    //
    ppp.query("particle_restart_file", particle_restart_file);
    //
    // This must be true the first time you try to restart from a checkpoint
    // that was written with USE_PARTICLES=FALSE; i.e. one that doesn't have
    // the particle checkpoint stuff (even if there are no active particles).
    // Otherwise the code will fail when trying to read the checkpointed particles.
    //
    ppp.query("restart_from_nonparticle_chkfile", restart_from_nonparticle_chkfile);
    //
    // Used in post_restart() to write out the file of particles.
    //
    ppp.query("particle_output_file", particle_output_file);
}

void
NavierStokesBase::initParticleData ()
{
    if (level == 0)
    {
        if (NSPC == 0)
        {
            NSPC = new AmrActiveParticleContainer(parent);
        }

        NSPC->SetVerbose(pverbose);

        if (!particle_init_file.empty())
        {
	  //NSPC->InitFromAsciiFile(particle_init_file,0);
	    NSPC->InitFromAsciiFile(particle_init_file,particle_extra_reals);	    
        }
    }
}

void
NavierStokesBase::post_restart_particle ()
{
    if (level == 0)
    {
        BL_ASSERT(NSPC == 0);

        NSPC = new AmrActiveParticleContainer(parent);

        NSPC->SetVerbose(pverbose);
        //
        // We want to be able to add new particles on a restart.
        // As well as the ability to write the particles out to an ascii file.
        //
        if (!restart_from_nonparticle_chkfile)
        {
            NSPC->Restart(parent->theRestartFile(), the_ns_particle_file_name);
        }

        if (!particle_restart_file.empty())
        {
	  //NSPC->InitFromAsciiFile(particle_restart_file,0);
	    NSPC->InitFromAsciiFile(particle_init_file,particle_extra_reals);	    
        }

        if (!particle_output_file.empty())
        {
            NSPC->WriteAsciiFile(particle_output_file);
        }
    }
}

void
NavierStokesBase::post_timestep_particle (int crse_iteration)
{
    const int ncycle = parent->nCycle(level);
    const int finest_level = parent->finestLevel();
    //
    // Don't redistribute/timestamp on the final subiteration except on the coarsest grid.
    //
    if (NSPC != 0 && (crse_iteration < ncycle || level == 0))
    {
        const Real curr_time = state[State_Type].curTime();

	int ngrow = (level == 0) ? 0 : crse_iteration;

        NSPC->Redistribute(level, finest_level, ngrow);

        if (!timestamp_dir.empty())
        {
            std::string basename = timestamp_dir;

            if (basename[basename.length()-1] != '/') basename += '/';

            basename += "Timestamp";

	    static bool first = true;
	    static int n, nextras;
	    static std::vector<int> tindices;

	    if (first)
	    {
		first = false;

		n = timestamp_indices.size();
		nextras = timestamp_num_extras();

		int sz = n + nextras;
		tindices.reserve(sz);

		for (int i = 0; i < sz; ++i) {
		    tindices.push_back(i);
		}
	    }

            for (int lev = level; lev <= finest_level; lev++)
            {
                if (NSPC->NumberOfParticlesAtLevel(lev) <= 0) continue;

		int ng = (lev == level) ? ngrow+1 : 1;

		AmrLevel& amr_level = parent->getLevel(lev);
		MultiFab& S_new = amr_level.get_new_data(State_Type);

		MultiFab tmf;

		if (tindices.size() > 0)
		{
		    tmf.define(S_new.boxArray(), S_new.DistributionMap(), tindices.size(), ng, MFInfo(), Factory());

		    if (n > 0)
		    {
		      FillPatchIterator fpi(parent->getLevel(lev), S_new,
                                            ng, curr_time, State_Type, 0, NUM_STATE);
                      const MultiFab& S = fpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel
#endif
                      for (MFIter mfi(tmf,true); mfi.isValid(); ++mfi)
                      {
                        FArrayBox& tfab = tmf[mfi];
                        const FArrayBox& sfab = S[mfi];
                        const Box& box = mfi.growntilebox();
                        for (int i = 0; i < n; ++i)
                        {
                          tfab.copy<RunOn::Host>(sfab, box, timestamp_indices[i], box, i, 1);
                        }
                      }
		    }

		    if (nextras > 0)
		    {
			timestamp_add_extras(lev, curr_time, tmf);
		    }
		}

		NSPC->Timestamp(basename, tmf, lev, curr_time, tindices);
            }
        }
    }
}

std::unique_ptr<MultiFab>
NavierStokesBase::ParticleDerive (const std::string& name,
				  Real               time,
				  int                ngrow)
{
    if (name == "particle_count" || name == "total_particle_count") {
	int ncomp = 1;
	const DeriveRec* rec = derive_lst.get(name);
	if (rec)
	{
	    ncomp = rec->numDerive();
	}

        MultiFab* ret = new MultiFab(grids, dmap, ncomp, ngrow, MFInfo(), Factory());
	ParticleDerive(name,time,*ret,0);
	return std::unique_ptr<MultiFab>{ret};
    }
    else {
	return AmrLevel::derive(name, time, ngrow);
    }
}

void
NavierStokesBase::ParticleDerive (const std::string& name,
				  Real               time,
				  MultiFab&          mf,
				  int                dcomp)
{
    if (NSPC == 0 || !(name == "particle_count" || name == "total_particle_count"))
    {
        AmrLevel::derive(name,time,mf,dcomp);
    }
    else {
	if (name == "particle_count")
	{
	    MultiFab temp_dat(grids,dmap,1,0, MFInfo(), Factory());
	    temp_dat.setVal(0);
	    NSPC->Increment(temp_dat,level);
	    MultiFab::Copy(mf,temp_dat,0,dcomp,1,0);
	}
	else if (name == "total_particle_count")
	{
	    //
	    // We want the total particle count at this level or higher.
	    //
	    ParticleDerive("particle_count",time,mf,dcomp);

	    IntVect trr(D_DECL(1,1,1));

	    for (int lev = level+1; lev <= parent->finestLevel(); lev++)
	    {
		BoxArray ba = parent->boxArray(lev);

		MultiFab temp_dat(ba,parent->DistributionMap(lev),1,0,MFInfo(),Factory());

		trr *= parent->refRatio(lev-1);

		ba.coarsen(trr);

                // FIXME? Won't work because ba has been coarsened. But don't actually need
                //  facotry here anyway...
		//MultiFab ctemp_dat(ba,parent->DistributionMap(lev),1,0,MFInfo(),Factory());
		MultiFab ctemp_dat(ba,parent->DistributionMap(lev),1,0);

		temp_dat.setVal(0);
		ctemp_dat.setVal(0);

		NSPC->Increment(temp_dat,lev);

#ifdef _OPENMP
#pragma omp parallel
#endif
		for (MFIter mfi(temp_dat,true); mfi.isValid(); ++mfi)
		{
		    const FArrayBox& ffab =  temp_dat[mfi];
		    FArrayBox&       cfab = ctemp_dat[mfi];
		    const Box&       fbx  = mfi.tilebox();

		    BL_ASSERT(cfab.box() == amrex::coarsen(fbx,trr));

		    for (IntVect p = fbx.smallEnd(); p <= fbx.bigEnd(); fbx.next(p))
		    {
		        const Real val = ffab(p);
			if (val > 0)
			    cfab(amrex::coarsen(p,trr)) += val;
		    }
		}

		temp_dat.clear();

		MultiFab dat(grids,dmap,1,0,MFInfo(),Factory());
		dat.setVal(0);
		dat.copy(ctemp_dat);

		MultiFab::Add(mf,dat,0,dcomp,1,0);
	    }
	}
	else
	{
	    amrex::Abort("NavierStokesBase::ParticleDerive: how did this happen?");
	}
    }
}

#endif  // AMREX_PARTICLES

// Boundary condition access function.
Vector<int>
NavierStokesBase::fetchBCArray (int State_Type, const Box& bx, int scomp, int ncomp)
{
    Vector<int> bc(2*BL_SPACEDIM*ncomp);
    BCRec bcr;
    const StateDescriptor* stDesc;
    const Box& domain = geom.Domain();

    for (int n = 0; n < ncomp; n++)
    {
      stDesc=state[State_Type].descriptor();
      setBC(bx,domain,stDesc->getBC(scomp+n),bcr);

      const int* b_rec = bcr.vect();
      for (int m = 0; m < 2*BL_SPACEDIM; m++)
	bc[2*BL_SPACEDIM*n + m] = b_rec[m];
    }

    return bc;
}

Vector<BCRec>
NavierStokesBase::fetchBCArray (int State_Type, int scomp, int ncomp)
{
    Vector<BCRec> bc(ncomp);
    const StateDescriptor* stDesc;
    const Box& domain = geom.Domain();

    for (int n(0); n < ncomp; ++n)
    {
      stDesc=state[State_Type].descriptor();
      setBC(domain,domain,stDesc->getBC(scomp+n), bc[n] );
    }

    return bc;
}

void
NavierStokesBase::average_down(const MultiFab& S_fine, MultiFab& S_crse,
			       int scomp, int ncomp)
{
  //
  // Choose the appropriate AMReX average_down() based on
  // whether EB or non-EB, and dimensionality
  //

#ifdef AMREX_USE_EB

  // FIXME?
  // Assume we want EB to behave the same as non-EB in regards to dimensionality
  // Not sure why we'd want 2D to be different than 3D
  // Note that 3D volume weighting doesn't exist for non-EB
#if (AMREX_SPACEDIM == 3)
    // no volume weighting
    amrex::EB_average_down(S_fine, S_crse, scomp, ncomp, fine_ratio);
#else
    // volume weighting
    amrex::EB_average_down(S_fine, S_crse, this->getLevel(level+1).Volume(),
			   *(this->getLevel(level+1).VolFrac()),
			   scomp, ncomp, fine_ratio);
#endif

#else
    // non-EB aware, uses volume weighting for 1D,2D but no volume weighting for 3D
    amrex::average_down(S_fine, S_crse,
			this->getLevel(level+1).geom, this->getLevel(level).geom,
			scomp, ncomp, fine_ratio);
#endif
}


//
//  Diagnostics functions
//
void
NavierStokesBase::printMaxVel (bool new_data)
{

    MultiFab& S = new_data? get_new_data(State_Type) : get_old_data(State_Type);

#if (AMREX_SPACEDIM==3)
    amrex::Print() << "max(abs(u/v/w))  = "
#else
        amrex::Print() << "max(abs(u/v))  = "
#endif
                   << S.norm0( Xvel,   0, false, true )
                   << "  "
		   << S.norm0( Xvel+1, 0, false, true )
#if (AMREX_SPACEDIM==3)
                   << "  "
                   << S.norm0( Xvel+2, 0, false, true )
#endif
                   << std::endl;
}


void
NavierStokesBase::printMaxGp (bool new_data)
{
#ifdef AMREX_USE_EB
    MultiFab& Gp = getGradP();
#else
    MultiFab Gp(grids,dmap,BL_SPACEDIM,1);
    const Real time = new_data ? state[Press_Type].curTime() : state[Press_Type].prevTime();
    getGradP(Gp, time);
#endif
    MultiFab& P  = new_data? get_new_data(Press_Type) : get_old_data(Press_Type);

#if (AMREX_SPACEDIM==3)
    amrex::Print() << "max(abs(gpx/gpy/gpz/p)) = "
#else
    amrex::Print() << "max(abs(gpx/gpy/p)) = "
#endif
                   << Gp.norm0( 0, 0, false, true )
                   << "  "
		   << Gp.norm0( 1, 0, false, true )
#if (AMREX_SPACEDIM==3)
                   << "  "
                   << Gp.norm0( 2, 0, false, true )
#endif
                   << "  "
                   << P.norm0(0, 0, false, true )
                   << std::endl;
}

void
NavierStokesBase::printMaxValues (bool new_data)
{
    printMaxVel(new_data);
    printMaxGp(new_data);
}


//
// Correct a conservatively-advected scalar for under-over shoots.
//
void
NavierStokesBase::ConservativeScalMinMax ( amrex::MultiFab&       Snew, const int snew_comp, const int new_density_comp,
                                           amrex::MultiFab const& Sold, const int sold_comp, const int old_density_comp )
{

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(Snew,true); mfi.isValid(); ++mfi)
    {

        const auto& bx = mfi.tilebox();

        const auto& sn   = Snew.array(mfi,snew_comp);
        const auto& so   = Sold.const_array(mfi,sold_comp);
        const auto& rhon = Snew.const_array(mfi,new_density_comp);
        const auto& rhoo = Sold.const_array(mfi,old_density_comp);

        amrex::ParallelFor(bx, [sn, so, rhon, rhoo]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real smn = 0.0;
            Real smx = 0.0;

#if (AMREX_SPACEDIM==3)
            int ks = -1;
            int ke =  1;
#else
            int ks = 0;
            int ke = 0;
#endif

            for (int kk = ks; kk <= ke; ++kk)
            {
                for (int jj = -1; jj <= 1; ++jj)
                {
                    for (int ii = -1; ii <= 1; ++ii)
                    {
                        smn =  amrex::min(smn, so(i+ii,j+jj,k+kk)/rhoo(i+ii,j+jj,k+kk));
                        smx =  amrex::max(smn, so(i+ii,j+jj,k+kk)/rhoo(i+ii,j+jj,k+kk));
                    }
                }
            }

            sn(i,j,k) = amrex::min( amrex::max(sn(i,j,k)/rhon(i,j,k), smn), smx ) * rhon(i,j,k);
        });
    }
}



//
// Correct a convectively-advected  scalar for under-over shoots.
//
void
NavierStokesBase::ConvectiveScalMinMax ( amrex::MultiFab&       Snew, const int snew_comp,
                                         amrex::MultiFab const& Sold, const int sold_comp )
{

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(Snew,true); mfi.isValid(); ++mfi)
    {

        const auto& bx = mfi.tilebox();

        const auto& sn   = Snew.array(mfi,snew_comp);
        const auto& so   = Sold.const_array(mfi,sold_comp);

        amrex::ParallelFor(bx, [sn, so]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real smn = 0.0;
            Real smx = 0.0;

#if (AMREX_SPACEDIM==3)
            int ks = -1;
            int ke =  1;
#else
            int ks = 0;
            int ke = 0;
#endif

            for (int kk = ks; kk <= ke; ++kk)
            {
                for (int jj = -1; jj <= 1; ++jj)
                {
                    for (int ii = -1; ii <= 1; ++ii)
                    {
                        smn =  amrex::min(smn, so(i+ii,j+jj,k+kk));
                        smx =  amrex::max(smn, so(i+ii,j+jj,k+kk));
                    }
                }
            }

            sn(i,j,k) = amrex::min( amrex::max(sn(i,j,k), smn), smx );
        });
    }
}
