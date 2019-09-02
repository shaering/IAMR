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

void
NavierStokesBase::getForce (FArrayBox&       force,
			    const Box&       bx,
			    int              ngrow,
			    int              scomp,
			    int              ncomp,
                const Real       time,
			    const FArrayBox& Vel,
			    const FArrayBox& Scal,
			    int              scalScomp)
{
    force.resize(amrex::grow(bx,ngrow),ncomp);

    const Real grav = gravity;

    for (int dc = 0; dc < ncomp; dc++)
    {
        const int sc = scomp + dc;
#if (BL_SPACEDIM == 2)
        if (sc == Yvel && std::fabs(grav) > 0.001)
#endif
#if (BL_SPACEDIM == 3)
        if (sc == Zvel && std::fabs(grav) > 0.001)
#endif
        {
            //
            // Set force to -rho*g.
            //
     	    force.copy(Scal,scalScomp,dc,1);
            force.mult(grav,dc,1);
        }
        else
        {
            force.setVal(0,dc);
        }
    }
}
