#ifndef _FLUIDOPUTILS_H_
#define _FLUIDOPUTILS_H_

#include "Diffusion.H"
#include "FluxSurface.H"

#include "GridFunction.H"
#include "GridFunctionLibrary.H"

#include "NamespaceHeader.H"

/// Class to support fluid operator utils.
/**
 * This class provides utilities for fluid operators, e.g., hyperviscosity.
 *
*/
class FluidOpUtils

{
public:

   /// Constructor with initialization.
   /**
    *  @param[in] pp the input database.
    *  @param[in] geom the geometry data.
    */
   FluidOpUtils(const ParmParse& pp,
                const MagGeom&   geom,
		const double     larmor,
                const int        verbosity);
      
   /// Destructor.
   /**
    */
   virtual ~FluidOpUtils();

   /// Initialize second-order (d^2/dx^2) elliptic operator for hypervoscisoty calculations
   void initializeHyperViscosityOp(const ParmParse&   pp,
                                   EllipticOpBC&      bc,
                                   bool               use_hyperviscosity_bcs);

   /// Adds 4th-order (d^4/dx^4) hyperviscosity of a_soln to a_rhs
   void addHyperViscosity(LevelData<FArrayBox>&         rhs,
                          const LevelData<FArrayBox>&   soln,
                          const bool                    is_mapped = true);
   
   
   void initializeBoundaryBuffers(const ParmParse& pp,
                                  const bool       include_full_relaxation);
   
   void suppressNonZonalCompInBoundaryBuffers(LevelData<FArrayBox>& data);
   
   void addRelaxationInBoundaryBuffers(LevelData<FArrayBox>&         rhs,
                                       const LevelData<FArrayBox>&   soln_phys,
                                       const bool                    is_mapped = true);

   void computeExBAdvection(LevelData<FArrayBox>& ExB_advection,
                            const LevelData<FArrayBox>& soln,
                            const LevelData<FluxBox>& E_field,
                            const std::string& advection_scheme,
                            const bool extrap_to_ghosts = false);
   
   void computeExBGradData(LevelData<FArrayBox>&       ExB_grad_data,
                           const LevelData<FArrayBox>& data,
                           const LevelData<FArrayBox>& E_field,
                           const std::string&          advection_scheme,
                           const bool                  extrap_to_ghosts = false);
  
   void subtractFSaverage(LevelData<FArrayBox>&  a_data);
   
   /// Computes stable dt for advection operator
   Real computeAdvectionDt(const LevelData<FluxBox>&  faceVel_mapped,
                           const std::string&         flux_type,
                           bool                       time_step_diagnostics);
   
   /// parseParameters function
   /**
    * parseParameters function
    */
   void parseParameters( const ParmParse& pp );

   /// printParameters function
   /**
    * printParameters function
    */
   void printParameters();

   /// ParmParse object name
   /**
    * ParmParse object name.  Used to locate input file parameters of interest to
    * this class.
    */
   
private:

   int m_verbosity;
   Real m_relaxation_coefficient;
   Real m_larmor;
   
   const MagGeom& m_geometry;

   FluxSurface m_flux_surface;
   
   Diffusion*   m_hyperviscosity_op;
   bool m_use_hyperviscosity_bcs;
   bool m_mapped_hyperviscosity;
   
   int m_hyperviscosity_order;
   
   LevelData<FArrayBox> m_boundary_buffer_profile;
   LevelData<FArrayBox> m_ref_soln;
   LevelData<FArrayBox> m_volume;

};



#include "NamespaceFooter.H"

#endif
