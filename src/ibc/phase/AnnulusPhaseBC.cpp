#include "AnnulusPhaseBC.H"
#include "ConstFact.H"
#include "Directions.H"
#include "FourthOrderUtil.H"
#include "KineticFunctionLibrary.H"
#include "LevelData.H"
#include "PhaseGeom.H"
#include "SPMD.H"
#include "PhaseBCUtils.H"

#include <sstream>

#include "EdgeToCell.H"

#include "FourthOrderBC.H.multidim"
#include "CodimBC.H.multidim"
#include "BCUtils.H.multidim"

//#include "inspect.H"

#include "NamespaceHeader.H"

/////// BEGIN INLINE METHODS ///////////////////////////////////////////////////

inline
KineticFunction& AnnulusPhaseBC::radialInflowFunc( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return *(m_inflow_function[RADIAL_INNER]);
   }
   return *(m_inflow_function[RADIAL_OUTER]);
}


inline
KineticFunction& AnnulusPhaseBC::vParallelInflowFunc( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return *(m_inflow_function[VPAR_LOWER]);
   }
   return *(m_inflow_function[VPAR_UPPER]);
}


inline
KineticFunction& AnnulusPhaseBC::muInflowFunc( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return *(m_inflow_function[MU_LOWER]);
   }
   return *(m_inflow_function[MU_UPPER]);
}


inline
KineticFunction& AnnulusPhaseBC::inflowFunc( const int& a_dir,
                                              const Side::LoHiSide& a_side )
{
   KineticFunction* inflow_func;
   if (a_dir==RADIAL_DIR) {
      inflow_func = &radialInflowFunc( a_side );
   }
   else if (a_dir==VPARALLEL_DIR) {
      inflow_func = &vParallelInflowFunc( a_side );
   }
   else if (a_dir==MU_DIR) {
      inflow_func = &muInflowFunc( a_side );
   }
   else {
      MayDay::Error( "AnnulusPhaseBC: Toroidal BCs not implemented!" );
   }
   return *(inflow_func);
}

/////// END INLINE METHODS /////////////////////////////////////////////////////


AnnulusPhaseBC::AnnulusPhaseBC( const std::string& a_name,
                                ParmParse& a_pp,
                                const int& a_verbosity )
   : m_name(a_name),
     m_verbosity(a_verbosity)
{
   m_inflow_function.resize( NUM_INFLOW );

   m_bdry_name.resize( NUM_INFLOW );
   m_bdry_name[RADIAL_INNER] = "radial_inner";
   m_bdry_name[RADIAL_OUTER] = "radial_outer";
   m_bdry_name[VPAR_LOWER] = "vpar_lower";
   m_bdry_name[VPAR_UPPER] = "vpar_upper";
   m_bdry_name[MU_LOWER] = "mu_lower";
   m_bdry_name[MU_UPPER] = "mu_upper";

   parseParameters( a_pp );
}


AnnulusPhaseBC::~AnnulusPhaseBC()
{
#if 0
   for (int i(0); i<m_inflow_function.size(); i++) {
      if (!(m_inflow_function[i])) {
         delete m_inflow_function[i];
      }
   }
#endif
}


inline
void AnnulusPhaseBC::fillInflowData( KineticSpeciesPtrVect& a_bdry_data,
                                     const BoundaryBoxLayoutPtrVect& a_bdry_layout,
                                     const Real& a_time )
{
   for (int i(0); i<a_bdry_layout.size(); i++) {
      const BoundaryBoxLayout& bdry_layout( *(a_bdry_layout[i]) );
      const int& dir( bdry_layout.dir() );
//      if (dir!=POLOIDAL_DIR) {
         const Side::LoHiSide& side( bdry_layout.side() );
         KineticSpecies& bdry_data( *(a_bdry_data[i]) );
         KineticFunction& inflow_func( inflowFunc( dir, side ) );
         inflow_func.assign( bdry_data, bdry_layout, a_time );
//      }
   }
}


inline
Box computeFillBox( const Box&                a_interior_box,
                    const MultiBlockCoordSys& a_coord_sys,
                    const int&                a_dir,
                    const Side::LoHiSide&     a_side,
                    const int&                a_n_ghost )
{
   const Box& domain_box( BCUtils::blockDomain( a_interior_box, a_coord_sys ) );
   return Box( BCUtils::getGhostBox( domain_box,
                                     a_interior_box,
                                     a_dir,
                                     a_side,
                                     a_n_ghost ) );
}


void AnnulusPhaseBC::apply( KineticSpecies& a_species_comp,
                            const CFG::LevelData<CFG::FArrayBox>& a_phi,
                            const LevelData<FluxBox>& a_velocity,
                            const Real& a_time )
{
   const PhaseGeom& geometry( a_species_comp.phaseSpaceGeometry() );
   const MultiBlockCoordSys& coord_sys( *(geometry.coordSysPtr()) );

   LevelData<FArrayBox>& BfJ( a_species_comp.distributionFunction() );
   const DisjointBoxLayout& grids( BfJ.disjointBoxLayout() );
   const IntVect& ghost_vect( BfJ.ghostVect() );

   BoundaryBoxLayoutPtrVect all_bdry_layouts;
   PhaseBCUtils::defineBoundaryBoxLayouts( all_bdry_layouts, grids, coord_sys, ghost_vect );

   KineticSpeciesPtrVect all_bdry_data;
   PhaseBCUtils::defineInflowDataStorage( all_bdry_data, all_bdry_layouts, a_species_comp );
   
   fillInflowData( all_bdry_data, all_bdry_layouts, a_time );

   if (m_radial_extrapolate) {

      // This option is a workaround for neoclassical simulations only.  It
      // sets the velocity on radial boundaries to zero in order to force the
      // ghost cells to be filled by extrapolation from the interior (i.e., always
      // outflow) in the fortran routine FOURTH_ORDER_OUTFLOW_BC.

      LevelData<FluxBox> tmp_velocity;
      tmp_velocity.define(a_velocity);

      for (DataIterator dit(tmp_velocity.dataIterator()); dit.ok(); ++dit) {
         tmp_velocity[dit][RADIAL_DIR].setVal(0.);
      }

      PhaseBCUtils::setInflowOutflowBC( BfJ,
                                        all_bdry_layouts,
                                        all_bdry_data,
                                        coord_sys,
                                        tmp_velocity );
   }
   else {
      PhaseBCUtils::setInflowOutflowBC( BfJ,
                                        all_bdry_layouts,
                                        all_bdry_data,
                                        coord_sys,
                                        a_velocity );
   }

   // interpolate all other codim boundaries
   CodimBC::setCodimBoundaryValues( BfJ, coord_sys );
}


void AnnulusPhaseBC::printParameters() const
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "AnnulusPhaseBC ================================" << std::endl;
      std::cout << "- variable: "  << m_name << "-------------" << std::endl;
      for (int i(0); i<m_inflow_function.size(); i++) {
         std::cout << "  " << m_bdry_name[i] << ": " << std::endl;
         m_inflow_function[i]->printParameters();
      }
      std::cout << "-----------------------------------------------" << std::endl;
      std::cout << "===============================================" << std::endl;
   }
}


inline
void AnnulusPhaseBC::parseParameters( ParmParse& a_pp )
{
   KineticFunctionLibrary* library = KineticFunctionLibrary::getInstance();
   for (int i(0); i<m_inflow_function.size(); i++) {
      std::string prefix( a_pp.prefix() );
      prefix += "." + m_bdry_name[i];
      ParmParse fpp( prefix.c_str() );
      std::string function_name;
      fpp.query( "function", function_name );
      m_inflow_function[i] = library->find( function_name );
   }

   if (a_pp.contains("radial_extrapolate")) {
      a_pp.get("radial_extrapolate", m_radial_extrapolate);
   }
   else {
      m_radial_extrapolate = false;
   }

   if (m_verbosity) {
      printParameters();
   }
}


#include "NamespaceFooter.H"

