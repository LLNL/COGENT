#ifndef _MBSOLVERUTIL_H_
#define _MBSOLVERUTIL_H_

#include "MultiBlockLevelGeom.H"
#include "MBSolver.H"

#include "NamespaceHeader.H"


class MBSolverUtil
{
public:

   /// Constructor with initialization.
   /**
    *  @param[in] geom the geometry
    *  @param[in] discretization_ordr discretization order
    */
   MBSolverUtil( const MBSolver& mbsolver );
      
   /// Destructor.
   /**
    */
   ~MBSolverUtil();

   void computeMappedCoefficients( const LevelData<FluxBox>&  unmapped_coefficients,
                                   LevelData<FluxBox>&        mapped_coefficients ) const;

   void getCellVolumes( LevelData<FArrayBox>& a_volume ) const;

protected:

   void fillTransverseGhosts( LevelData<FluxBox>& a_data ) const;

   void getPointwiseN( LevelData<FluxBox>& a_N ) const;

   void getPointwiseJinverse( LevelData<FluxBox>& a_NJinverse ) const;

   const MultiBlockLevelGeom& m_geometry;
   MultiBlockCoordSys* m_coord_sys_ptr;
   int m_discretization_order;

};


#include "NamespaceFooter.H"

#endif
