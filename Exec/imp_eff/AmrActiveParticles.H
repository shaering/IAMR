#ifndef AMREX_AmrParticles_H_
#define AMREX_AmrParticles_H_

#include <AMReX_Particles.H>
#include <AMReX_TracerParticles.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_Interpolater.H>
#include <AMReX_FillPatchUtil.H>
#include <activeParticles.H>

namespace amrex {

class AmrActiveParticleContainer
    : public ActiveParticleContainer
{
public:

    AmrActiveParticleContainer (AmrCore* amr_core)
	: ActiveParticleContainer(amr_core->GetParGDB())
    {
    }

    ~AmrActiveParticleContainer () {}
};

}

#endif
