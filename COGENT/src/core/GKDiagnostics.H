#ifndef _GKDIAGNOSTICS_H_
#define _GKDIAGNOSTICS_H_

#include "GKOps.H"

#include "NamespaceHeader.H"


class GKDiagnostics
{
  public:

    GKDiagnostics() 
    { 
      m_ops = nullptr;
      m_phase_geometry = nullptr;
      m_units = nullptr;
      m_fluid_op = nullptr;
    }

    GKDiagnostics( const GKOps* const a_ops ) 
    { 
      m_ops = a_ops;
      m_phase_geometry = m_ops->phaseGeometry();
      m_units = m_ops->units();
      m_fluid_op = m_ops->fluidOp();
    }

    ~GKDiagnostics() {}

    void getCfgVar(  CFG::LevelData<CFG::FArrayBox>&,
                     const std::string&,
                     const bool& ) const;
  
    void getCfgVar( CFG::LevelData<CFG::FArrayBox>&,
                    const std::string&,
                    const KineticSpecies& ) const;
 
    void getCfgVar( CFG::LevelData<CFG::FArrayBox>&,
                    const std::string&,
                    const KineticSpeciesPtrVect& ) const;
 
    void getPhaseVar(  LevelData<FArrayBox>&,
                       const std::string&,
                       const KineticSpecies& ) const;

    void writeDfnPlots(  const std::vector<std::string>&,
                         const std::vector<std::string>&,
                         const KineticSpecies&,
                         const Real,
                         const int,
                         const int,
                         const int,
                         const int,
                         const int ) const;

    void defineFluidCellVar( CFG::LevelData<CFG::FArrayBox>&,
                             const CFG::FluidSpecies&,
                             const string& ) const;

    void defineFluidFaceVarAtCell( CFG::LevelData<CFG::FArrayBox>&,
                                   const CFG::FluidSpecies&,
                                   const string& ) const;

    void defineFluidEdgeVarAtCell( CFG::LevelData<CFG::FArrayBox>&,
                                   const CFG::FluidSpecies&,
                                   const string& ) const;

    void defineFluidOpMember(CFG::LevelData<CFG::FArrayBox>&,
                             const CFG::FluidSpecies&,
                             const std::string& ) const;
 
    void getFluidCellVar(  CFG::LevelData<CFG::FArrayBox>&,
                           const CFG::CFGVars&,
                           const string& ) const;

    void getFluidFaceVarAtCell(  CFG::LevelData<CFG::FArrayBox>&,
                                 const CFG::CFGVars&,
                                 const string& ) const;

    void getFluidEdgeVarAtCell(  CFG::LevelData<CFG::FArrayBox>&,
                                 const CFG::CFGVars&,
                                 const string& ) const;

    void getFluidNodeVarAtCell(  CFG::LevelData<CFG::FArrayBox>&,
                                 const CFG::CFGVars&,
                                 const string& ) const;

    void getFluidOpMember( CFG::LevelData<CFG::FArrayBox>&,
                           const CFG::CFGVars&,
                           const std::string& ) const;

    template <typename T>
    void defineCfgVar( T&                  a_var,
                       const int           a_ncomp, 
                       const CFG::IntVect  a_gv = CFG::IntVect::Zero ) const
    {
      CH_assert( m_phase_geometry != NULL );
      const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
      a_var.define( mag_geom.gridsFull(), a_ncomp, a_gv );
    }
 
    template <typename T>
    void plotCfgVar( const T& a_var,
                     const std::string&            a_filename,
                     const double&                 a_time ) const
    {
      CH_assert( m_phase_geometry != NULL );
      m_phase_geometry->plotConfigurationData( a_filename.c_str(), 
                                                a_var, 
                                                a_time );
    }
 
    template <typename T>
    void plotPhaseVar( const T& a_var,
                       const PhaseGeom&    a_geom,
                       const std::string&  a_filename,
                       const double&       a_time ) const
    {
      a_geom.plotData( a_filename.c_str(), a_var, a_time );
    }
 
    template <typename T>
    void plotPhaseVarAtMu( const T& a_var,
                           const PhaseGeom&    a_geom,
                           const std::string&  a_filename,
                           const int           a_mu,
                           const double&       a_time ) const
    {
      a_geom.plotAtMuIndex( a_filename.c_str(), a_mu, a_var, a_time );
    }
 
    template <typename T>
    void plotPhaseVarVParTheta(  const T& a_var,
                                 const PhaseGeom&    a_geom,
                                 const std::string&  a_filename,
                                 const int           a_radial_index,
                                 const int           a_toroidal_index,
                                 const int           a_mu_index,
                                 const double&       a_time ) const
    {
      a_geom.plotVParPoloidalData( a_filename.c_str(), 
                                   a_radial_index,
                                   a_toroidal_index,
                                   a_mu_index, 
                                   a_var, 
                                   a_time );
    }
 
    template <typename T>
    void plotPhaseVarAtVelIndex( const T& a_var,
                                 const PhaseGeom&    a_geom,
                                 const std::string&  a_filename,
                                 const int           a_vpar_index,
                                 const int           a_mu_index,
                                 const double&       a_time ) const
    {
      VEL::IntVect vspace_index( a_vpar_index, a_mu_index );
      a_geom.plotAtVelocityIndex(  a_filename.c_str(), 
                                   vspace_index, 
                                   a_var, 
                                   a_time );
    }
 
    template <typename T>
    void plotPhaseVarAtCfgIndex( const T& a_var,
                                 const PhaseGeom&    a_geom,
                                 const std::string&  a_filename,
                                 const int           a_radial_index,
                                 const int           a_poloidal_index,
                                 const int           a_toroidal_index,
                                 const double&       a_time ) const
    {
#if CFG_DIM==3
       CFG::IntVect cspace_index(  a_radial_index, 
                                   a_poloidal_index, 
                                   a_toroidal_index);
#else
       CFG::IntVect cspace_index(  a_radial_index, 
                                   a_poloidal_index);
#endif
      a_geom.plotAtConfigurationIndex( a_filename.c_str(), 
                                       cspace_index, 
                                       a_var, 
                                       a_time );
    }

  protected:

    void getPotential( CFG::LevelData<CFG::FArrayBox>&, const bool& ) const;

    void getEField( CFG::LevelData<CFG::FArrayBox>&) const;
  
    void getExBData( CFG::LevelData<CFG::FArrayBox>& a_var ) const;
   
    void getDistributionFunction( LevelData<FArrayBox>&,
                                  const KineticSpecies& ) const;

    void getBStarParallel( LevelData<FArrayBox>&, const KineticSpecies& ) const;

    void getBStarParallelDfn( LevelData<FArrayBox>&, 
                              const KineticSpecies& ) const;

    void getDeltaF( LevelData<FArrayBox>&, const KineticSpecies& ) const;

    void getChargeDensity( CFG::LevelData<CFG::FArrayBox>&, 
                           const KineticSpecies& ) const;
 
    void getFluidVelocity( CFG::LevelData<CFG::FArrayBox>&,
                           const KineticSpecies& ) const;

    void getParallelParticleFlux( CFG::LevelData<CFG::FArrayBox>&,
                                  const KineticSpecies& ) const;

    void getPoloidalParticleFlux(CFG::LevelData<CFG::FArrayBox>&,
                                 const KineticSpecies& ) const;

    void getParallelVelocity(  CFG::LevelData<CFG::FArrayBox>&, 
                               const KineticSpecies& ) const;

    void getEnergyDensity(  CFG::LevelData<CFG::FArrayBox>&, 
                            const KineticSpecies& ) const;

    void getKineticEnergyDensity(  CFG::LevelData<CFG::FArrayBox>&, 
                                   const KineticSpecies& ) const;

    void getParallelEnergyDensity( CFG::LevelData<CFG::FArrayBox>&, 
                                   const KineticSpecies& ) const;

    void getPerpEnergyDensity( CFG::LevelData<CFG::FArrayBox>&, 
                               const KineticSpecies& ) const;

    void getPressure( CFG::LevelData<CFG::FArrayBox>&, 
                      const KineticSpecies& ) const;

    void getParallelPressure(  CFG::LevelData<CFG::FArrayBox>&, 
                               const KineticSpecies& ) const;

    void getPerpPressure(  CFG::LevelData<CFG::FArrayBox>&, 
                           const KineticSpecies& ) const;

    void getRadialForceBalance( CFG::LevelData<CFG::FArrayBox>&,
                                const KineticSpecies& ) const;

    void getTemperature( CFG::LevelData<CFG::FArrayBox>&, 
                         const KineticSpecies& ) const;

    void getParallelTemperature( CFG::LevelData<CFG::FArrayBox>&, 
                                 const KineticSpecies& ) const;

    void getPerpTemperature( CFG::LevelData<CFG::FArrayBox>&, 
                             const KineticSpecies& ) const;

    void getParallelHeatFlux( CFG::LevelData<CFG::FArrayBox>&, 
                              const KineticSpecies& ) const;

    void getTotalParallelHeatFlux( CFG::LevelData<CFG::FArrayBox>&, 
                                   const KineticSpecies& ) const;

    void getFourthMoment( CFG::LevelData<CFG::FArrayBox>&, 
                          const KineticSpecies& ) const;

    void getRadialCurrent(CFG::LevelData<CFG::FArrayBox>&,
                          const KineticSpecies& ) const;

    void getRadialParticleFlux(CFG::LevelData<CFG::FArrayBox>&,
                               const KineticSpecies& ) const;
  
    void getRadialHeatFlux(CFG::LevelData<CFG::FArrayBox>&,
                           const KineticSpecies& ) const;

    void getExBHeatFlux( CFG::LevelData<CFG::FArrayBox>&, 
                         const KineticSpecies& ) const;

    void getExBKineticEnergyFlux(  CFG::LevelData<CFG::FArrayBox>&, 
                                   const KineticSpecies& ) const;

    void getChargeDensity( CFG::LevelData<CFG::FArrayBox>&, 
                           const KineticSpeciesPtrVect& ) const;
  
    void computeRadialFlux( CFG::LevelData<CFG::FArrayBox>&,
                            const KineticSpecies&,
                            const CFG::LevelData<CFG::FArrayBox>&,
                            const CFG::EMFields&,
                            const int,
                            const Kernel<FluxBox>&) const;

    const GKOps* m_ops;
    const PhaseGeom* m_phase_geometry;
    const GKUnits* m_units;
    const CFG::GKFluidOp* m_fluid_op;

  private:

};

#include "NamespaceFooter.H"
#endif
