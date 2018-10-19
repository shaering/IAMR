
#include <cstdio>

#include <AMReX_CArena.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BLProfiler.H>

#ifdef USE_AMR_MLSDC

#include <MLSDC_Amr.H>
#include <MLSDC_Context.H>
#include <MLSDC_LevelBld.H>

#endif


using namespace amrex;


#ifdef USE_AMR_MLSDC

extern "C"
{
  /* Driver function for pfasst control */
  void fmain(MLSDC_Context* mlsdc_ctx,
  	     const int*     nlevels,   
  	     const int*     niters,            
             const int      nnodes[],                                                              
             const char*    qtype_name,                            
  	     const int*     qtype_name_len, 
  	     double*        t_max,  
  	     double*        dt);
  
  /* Debug function to test the data structures */
  // void ftest(MLSDC_Context* mlsdc_ctx,
  //  	     const int*     nlevels,   
  //  	     const int*     niters,            
  //            const int      nnodes[], 
  //            const char*    qtype_name,                            
  //  	     const int*     qtype_name_len, 
  //  	     double*        t_max,  
  //  	     double*        dt);
  
}

#endif

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    BL_PROFILE_REGION_START("main()");
    BL_PROFILE_VAR("main()", pmain);

    const Real run_strt = ParallelDescriptor::second();

    int  max_step;
    int  num_steps;
    Real strt_time;
    Real stop_time;

    ParmParse pp;

    amrex::Print() << " This is my main " << std::endl;

    max_step  = -1; 
    num_steps = -1; 
    strt_time =  0.0;
    stop_time = -1.0;

    pp.query("max_step",  max_step);
    pp.query("num_steps", num_steps);
    pp.query("strt_time", strt_time);
    pp.query("stop_time", stop_time);

    stop_time = 1e-5; 

    max_step  = 64;

    if (strt_time < 0.0)
    {
        amrex::Abort("MUST SPECIFY a non-negative strt_time");
    }

    if (max_step < 0 && stop_time < 0)
    {
        amrex::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

#ifdef USE_AMR_MLSDC

    // define the libPFASST parameters
    const int n_sdc_levels  = 1;
    const int n_sdc_iters   = 6;
    const int n_sdc_nodes[n_sdc_levels] = { 3 };
    std::string nodes_type  = "SDC_GAUSS_LOBATTO";
    const int string_length = nodes_type.size();
    double t_max            = stop_time;

    // instantiate the MLSDC_Amr object
    MLSDC_Amr* amrptr        = new MLSDC_Amr;

    // instantiate the MLSDC_Context object
    MLSDC_Context* mlsdc_ctx = new MLSDC_Context(amrptr,
						 n_sdc_levels,
						 n_sdc_nodes,
						 n_sdc_iters);

    // set the pointer in the factory object
    MLSDC_LevelBld* mlsdc_level_bld = getMLSDCLevelBld();
    mlsdc_level_bld->setContextPtr(mlsdc_ctx);

#else

    // instantiate the Amr object
    Amr* amrptr = new Amr;

#endif

    // initialize the Amr object
    amrptr->init(strt_time,stop_time);   

    if (num_steps > 0)
    {
        if (max_step < 0)
        {
            max_step = num_steps + amrptr->levelSteps(0);
        }
        else
        {
            max_step = std::min(max_step, num_steps + amrptr->levelSteps(0));
        }

    	amrex::Print() << "Using effective max_step = " << max_step << '\n';
    }
    //
    // If we set the regrid_on_restart flag and if we are *not* going to take
    // a time step then we want to go ahead and regrid here.
    //
    if (amrptr->RegridOnRestart())
    {
        if (    (amrptr->levelSteps(0) >= max_step ) ||
                ( (stop_time >= 0.0) &&
                  (amrptr->cumTime() >= stop_time)  )    )
        {
            //
            // Regrid only!
            //
            amrptr->RegridOnly(amrptr->cumTime());
        }
    }

#ifdef USE_AMR_MLSDC

    num_steps = 4;
    double dt = (stop_time - strt_time)/((double)num_steps);  

    amrex::Print() << "dt = " << dt 
		   << " t_max = " << t_max 
		   << " strt_time = " << strt_time
		   << " stop_time = " << stop_time 
		   << " num_steps = " << num_steps
		   << '\n';

    // time stepping handled by libPFASST
    fmain(mlsdc_ctx,          // user defined context
	  &n_sdc_levels,      // number of sdc space-time levels
	  &n_sdc_iters,       // number of sdc iterations
	  n_sdc_nodes,        // number of sdc nodes for each sdc level
	  nodes_type.c_str(), // type of nodes
	  &string_length,     // length of nodes_type.c_str()
	  &t_max,             // simulation time
	  &dt                 // time step
	  );
    
    amrex::Print() << "###############################" << std::endl;
    amrex::Print() << "     Done with fmain           " << std::endl;
    amrex::Print() << "###############################" << std::endl;
    
#else 

    while ( amrptr->okToContinue()                            &&
           (amrptr->levelSteps(0) < max_step || max_step < 0) &&
           (amrptr->cumTime() < stop_time || stop_time < 0.0) )
    {

      amrptr->coarseTimeStep(stop_time);

    }

#endif
    //
    // Write final checkpoint and plotfile.
    //
    if (amrptr->stepOfLastCheckPoint() < amrptr->levelSteps(0))
    {
        amrptr->checkPoint();
    }

    //if (amrptr->stepOfLastPlotFile() < amrptr->levelSteps(0))
    //{
    amrptr->writePlotFile();
    //}

#ifdef USE_AMR_MLSDC
    amrex::Print() << "###############################" << std::endl;
    amrex::Print() << "     Delete context object     " << std::endl;
    amrex::Print() << "###############################" << std::endl;

    delete mlsdc_ctx;

#endif

    delete amrptr;

    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_stop = ParallelDescriptor::second() - run_strt;

    ParallelDescriptor::ReduceRealMax(run_stop,IOProc);

    amrex::Print() << "Run time = " << run_stop << std::endl;

    BL_PROFILE_VAR_STOP(pmain);
    BL_PROFILE_REGION_STOP("main()");
    BL_PROFILE_SET_RUN_TIME(run_stop);
    BL_PROFILE_FINALIZE();


    amrex::Finalize();

    return 0;
}
