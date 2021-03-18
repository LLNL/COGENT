#ifndef _GYROAVGOP_H_
#define _GYROAVGOP_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <map>
#include <string>
#include <math.h>

#include "Dimensions.H"
#include "REAL.H"
#include "ParmParse.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "PhaseCoordSys.H"
#include "UnstructuredInterpUtils.H"

#define _MU_ 1

#define _point_ first
#define _weight_ second

#include "NamespaceHeader.H"

namespace CFG = CFG_NAMESPACE;
namespace VEL = VEL_NAMESPACE;

class PhaseGeom;

class GyroaverageOperator
{
  public:

    GyroaverageOperator() { m_is_defined = false; } 
    ~GyroaverageOperator() { }

    void define(const PhaseGeom&, const std::string&, Real, Real);
    void define(const GyroaverageOperator&);

    void applyOp(LevelData<FArrayBox>&, 
                 const CFG::LevelData<CFG::FArrayBox>& ) const;

    void applyOp(LevelData<FArrayBox>&, 
                 const LevelData<FArrayBox>& ) const;

    void applyOpTwice(LevelData<FArrayBox>&, 
                      const LevelData<FArrayBox>& ) const;

    void gyroaveragedEField(  LevelData<FArrayBox>&,
                              const CFG::LevelData<CFG::FArrayBox>&,
                              const int ) const;

    void gyroaveragedEField(  LevelData<FluxBox>&,
                              const CFG::LevelData<CFG::FArrayBox>&,
                              const int ) const;

    inline const LevelData<FArrayBox>& gyroradius() const
    {
      CH_assert(m_debug_mode);
      return m_gyroradius;
    }

    inline const LevelData< BaseFab<Stencil> >& op() const
    {
      return m_gyroavg_op;
    }

    inline bool isDefined() const { return m_is_defined; }

  protected:

  private:

    bool  m_is_defined;
    int   m_npts_interp;
    int   m_gpt;

    std::string m_name;
    Real        m_mass;
    Real        m_charge;
    Real        m_tolerance;
    int         m_max_muidx;
    Real        m_max_gyroradius;
    int         m_interp_order;
    bool        m_reduce_interp_order;
    bool        m_interp_mode_mapped;

    const PhaseGeom* m_phase_geom;

    DisjointBoxLayout m_grids;

    LevelData<FArrayBox>          m_gyroradius;
    LevelData< BaseFab<Stencil> > m_gyroavg_op;

    /* some debug and testing stuff while we are still verifying
     * this operator - one day, they'll disappear. */
    bool m_debug_mode;
    bool m_debug_fakegyroradius;
    Real m_debug_fakegyroradius_dr;

    void parseParams(const ParmParse&);
    void printParams() const;

    void getInterpStencil(Stencil&, 
                          const CFG::FArrayBox&, 
                          const std::vector<CFG::IntVect>&, 
                          const CFG::RealVect&,
                          const int,
                          bool&,
                          Real& );

    void applyOp( FArrayBox&, 
                  const FArrayBox&, 
                  const BaseFab<Stencil>&) const;

    inline Real computeGyroradius(Real a_La,      // Larmor number
                                  Real a_mass,    // species mass
                                  Real a_charge,  // species charge
                                  Real a_B,       // magnetic field magnitude
                                  Real a_mu       // magnetic moment
                                 ) const
    {
      if ((a_B == 0) || (a_charge == 0)) return(0.0);
      else return( (a_La/a_charge) * sqrt(a_mass * a_mu / a_B ) );
    }

    inline void evalStencil(std::vector<Real>&,
                            const FArrayBox&, 
                            const Stencil&,
                            const VEL::IntVect& ) const;

    inline Real distanceFunction( const CFG::FArrayBox& a_coords,
                                  const CFG::IntVect&   a_i,
                                  const CFG::RealVect&  a_x0 ) const
    {
      const CFG::Box& bx = a_coords.box();
      CH_assert(bx.contains(a_i));

      CFG::RealVect dx;
      for (int d=0; d<CFG_DIM; d++) {
        dx[d] = a_x0[d] - a_coords(a_i,d);
      }

      Real retval = 0.0;
      for (int d=0; d<CFG_DIM; d++) {
        retval += (dx[d]*dx[d]);
      }

      return retval;
    }

};

#include "NamespaceFooter.H"

#endif
