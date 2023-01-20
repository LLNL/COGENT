#ifndef _VELOCITYCOORDSYS_H_
#define _VELOCITYCOORDSYS_H_

#include "ParmParse.H"
#include "CartesianCS.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"

#include "NamespaceHeader.H"

/// Velocity space geometry base class
/**
 *  Base class for velocity space geometry.
 */
class VelCoordSys : public CartesianCS
{
   public:

      /// Constructors
      /**
       * Constructor with initialization.
       *
       * @param[in] parm_parse the ParmParse database.
       * @param[in] grids the box layout.
       * @param[in] domain the configuration space domain.
       * @param[in] dX the mesh spacing.
       */
      VelCoordSys( ParmParse&               parm_parse,
                   const DisjointBoxLayout& grids,
                   const ProblemDomain&     domain,
                   const RealVect&          dX );

      /// Destructor
      /**
       */
      virtual ~VelCoordSys();

      virtual void getMetricTerms( LevelData<FluxBox>& N,
                                   LevelData<FluxBox>& tanGradN ) const;

      void getCellVolumes( FArrayBox& volume ) const;

      virtual void getCellVolumes( LevelData<FArrayBox>& cellVolumes ) const;

      virtual void getJ( LevelData<FArrayBox>& J ) const;

      // Compute physical face areas
      void getFaceAreas( FluxBox& a_areas ) const;

      // Compute physical face areas
      void getFaceAreas( LevelData<FluxBox>& a_areas ) const;

      /// access function, returns mapped-space problem domain
      const ProblemDomain& domain() const {return m_domain;}

      /// access function -- returns const reference to grids
      const DisjointBoxLayout& grids() const {return m_grids;}

      void plotCellData( const string&               file_name,  
                         const LevelData<FArrayBox>& data,
                         const double&               time ) const;

      void plotFaceData( const string&             file_name,  
                         const LevelData<FluxBox>& data,
                         const double&             time ) const;

      /**
       * Object ParmParse name.
       */
      static const std::string pp_name;

   protected:

      DisjointBoxLayout m_grids;

      ProblemDomain m_domain;

      bool m_verbose;

      bool m_isDefined;
};



#include "NamespaceFooter.H"

#endif
