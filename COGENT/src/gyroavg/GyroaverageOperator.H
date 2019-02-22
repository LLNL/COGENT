#ifndef _GYROAVGOP_H_
#define _GYROAVGOP_H_

#include <iostream>
#include <vector>
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

#define _MU_ 1

#define _point_ first
#define _weight_ second

#include "NamespaceHeader.H"

namespace CFG = CFG_NAMESPACE;
namespace VEL = VEL_NAMESPACE;

typedef pair<CFG::IntVect,Real> StencilPoint;
typedef std::vector<StencilPoint> Stencil;

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

    void gyroaveragedEField(  LevelData<FArrayBox>&,
                              const CFG::LevelData<CFG::FArrayBox>&,
                              const int ) const;

    void gyroaveragedEField(  LevelData<FluxBox>&,
                              const CFG::LevelData<CFG::FArrayBox>&,
                              const int ) const;

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

    const PhaseGeom* m_phase_geom;

    DisjointBoxLayout m_grids;

    LevelData< BaseFab<Stencil> > m_gyroavg_op;

    void parseParams(const ParmParse&);
    void printParams() const;

    void getInterpStencil(Stencil&, const CFG::FArrayBox&, 
                          const CFG::IntVect&, const CFG::RealVect&);

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
      else return( a_La * sqrt(a_mass * a_mu / (a_charge*a_charge * a_B) ) );
    }

    inline void evalStencil(std::vector<Real>&    a_var_avg, 
                            const CFG::FArrayBox& a_var, 
                            const Stencil&        a_op) const
    {
      a_var_avg = std::vector<Real>(a_var.nComp(),0.0);
      for (int i=0; i<a_op.size(); i++) {
        for (int v=0; v<a_var.nComp(); v++) {
          a_var_avg[v] += a_op[i]._weight_ * a_var(a_op[i]._point_,v);
        }
      }
      return;
    }

};

#include "NamespaceFooter.H"

#endif
