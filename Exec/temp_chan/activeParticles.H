#ifndef AMREX_ACTIVEPARTICLES_H_
#define AMREX_ACTIVEPARTICLES_H_

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Particles.H>
#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_Interpolater.H>
#include <AMReX_FillPatchUtil.H>


namespace amrex {

// copied from electrostatic example
struct PIdx
{
  enum { // Particle Attributes stored in amrex::ParticleContainer's struct-of-arrays                                                               
    w = 0,
    vx, vy,
#if BL_SPACEDIM == 3
    vz,
#endif
    nattribs
  };
};


class MyParIter
  : public amrex::ParIter<0, 0, PIdx::nattribs>
{
public:
  using amrex::ParIter<0, 0, PIdx::nattribs>::ParIter;
  using RealVector = ParIter::ContainerType::RealVector;

  const std::array<RealVector, PIdx::nattribs>& GetAttribs () const {
    return GetStructOfArrays().GetRealData();
  }

  std::array<RealVector, PIdx::nattribs>& GetAttribs () {
    return GetStructOfArrays().GetRealData();
  }

  const RealVector& GetAttribs (int comp) const {
    return GetStructOfArrays().GetRealData(comp);
  }

  RealVector& GetAttribs (int comp) {
    return GetStructOfArrays().GetRealData(comp);
  }

};


  //using ActiveParIter = ParIter<AMREX_SPACEDIM>;
using VectorMeshData = Vector<std::array<std::unique_ptr<MultiFab>, BL_SPACEDIM> >;
using ScalarMeshData = Vector<std::unique_ptr<MultiFab> >;



class ActiveParticleContainer
  : public ParticleContainer<AMREX_SPACEDIM*2+3+4,1> //<AMREX_SPACEDIM*2+3,1>
{
public:

    ActiveParticleContainer (ParGDBBase* gdb)
      : ParticleContainer<AMREX_SPACEDIM*2+3+4,1>(gdb)
	{}

    ActiveParticleContainer (const Geometry            & geom, 
   			     const DistributionMapping & dmap,
			     const BoxArray            & ba)
      : ParticleContainer<AMREX_SPACEDIM*2+3+4,1>(geom,dmap,ba)
	{}

    ~ActiveParticleContainer () {}

  void basicAdvectWithUmac (MultiFab* umac, int level, Real dt);
  void myAdvectWithUmac (MultiFab* umac, int level, Real dt, MultiFab& rho, MultiFab& temp, Real nu_m, Real time);
  void myAdvectWithUcc (const MultiFab& ucc, int level, Real dt);
  void Timestamp (const std::string& file, const MultiFab& mf, int lev, Real time, const std::vector<int>& idx);
  void getDrag (FArrayBox& rhs, const FArrayBox& vel, const FArrayBox& rho, Real nu_m, int nGrow, int level);
  void getTemp (FArrayBox& rhs, const FArrayBox& vel, const FArrayBox& temp, Real nu_m, int nGrow, int level);

};


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


  using ActiveParIter = ParIter<AMREX_SPACEDIM>;

}

#endif

