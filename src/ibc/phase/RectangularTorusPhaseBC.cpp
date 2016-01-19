#include "RectangularTorusPhaseBC.H"
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
KineticFunction& RectangularTorusPhaseBC::radialInflowFunc( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return *(m_inflow_function[RADIAL_LOWER]);
   }
   return *(m_inflow_function[RADIAL_UPPER]);
}

inline
KineticFunction& RectangularTorusPhaseBC::poloidalInflowFunc( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return *(m_inflow_function[POLOIDAL_LOWER]);
   }
   return *(m_inflow_function[POLOIDAL_UPPER]);
}

#if CFG_DIM==3
inline
KineticFunction& RectangularTorusPhaseBC::toroidalInflowFunc( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return *(m_inflow_function[TOROIDAL_LOWER]);
   }
   return *(m_inflow_function[TOROIDAL_UPPER]);
}
#endif

inline
KineticFunction& RectangularTorusPhaseBC::vParallelInflowFunc( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return *(m_inflow_function[VPAR_LOWER]);
   }
   return *(m_inflow_function[VPAR_UPPER]);
}


inline
KineticFunction& RectangularTorusPhaseBC::muInflowFunc( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return *(m_inflow_function[MU_LOWER]);
   }
   return *(m_inflow_function[MU_UPPER]);
}


inline
KineticFunction& RectangularTorusPhaseBC::inflowFunc( const int& a_dir,
                                              const Side::LoHiSide& a_side )
{
   KineticFunction* inflow_func;
   if (a_dir==RADIAL_DIR) {
      inflow_func = &radialInflowFunc( a_side );
   }
   else if (a_dir==POLOIDAL_DIR) {
      inflow_func = &poloidalInflowFunc( a_side );
   }
#if CFG_DIM==3
   else if (a_dir==TOROIDAL_DIR) {
      inflow_func = &toroidalInflowFunc( a_side );
   }
#endif
   else if (a_dir==VPARALLEL_DIR) {
      inflow_func = &vParallelInflowFunc( a_side );
   }
   else if (a_dir==MU_DIR) {
      inflow_func = &muInflowFunc( a_side );
   }
   else {
      MayDay::Error( "RectangularTorusPhaseBC: BCs not implemented!" );
   }
   return *(inflow_func);
}

/////// END INLINE METHODS /////////////////////////////////////////////////////


RectangularTorusPhaseBC::RectangularTorusPhaseBC( const std::string& a_name,
                                ParmParse& a_pp,
                                const int& a_verbosity )
   : m_name(a_name),
     m_verbosity(a_verbosity)
{
   m_inflow_function.resize( NUM_INFLOW );

   m_bdry_name.resize( NUM_INFLOW );
   m_bdry_name[RADIAL_LOWER] = "radial_lower";
   m_bdry_name[RADIAL_UPPER] = "radial_upper";
   m_bdry_name[POLOIDAL_LOWER] = "poloidal_lower";
   m_bdry_name[POLOIDAL_UPPER] = "poloidal_upper";
#if CFG_DIM==3
   m_bdry_name[TOROIDAL_LOWER] = "toroidal_lower";
   m_bdry_name[TOROIDAL_UPPER] = "toroidal_upper";
#endif
   m_bdry_name[VPAR_LOWER] = "vpar_lower";
   m_bdry_name[VPAR_UPPER] = "vpar_upper";
   m_bdry_name[MU_LOWER] = "mu_lower";
   m_bdry_name[MU_UPPER] = "mu_upper";

   parseParameters( a_pp );
}


RectangularTorusPhaseBC::~RectangularTorusPhaseBC()
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
void RectangularTorusPhaseBC::fillInflowData( KineticSpeciesPtrVect& a_bdry_data,
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


void RectangularTorusPhaseBC::apply( KineticSpecies& a_species_comp,
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

   PhaseBCUtils::setInflowOutflowBC( BfJ,
                                     all_bdry_layouts,
                                     all_bdry_data,
                                     coord_sys,
                                     a_velocity );

   // interpolate all other codim boundaries
   CodimBC::setCodimBoundaryValues( BfJ, coord_sys );
}


void RectangularTorusPhaseBC::printParameters() const
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "RectangularTorusPhaseBC ================================" << std::endl;
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
void RectangularTorusPhaseBC::parseParameters( ParmParse& a_pp )
{
   KineticFunctionLibrary* library = KineticFunctionLibrary::getInstance();
   for (int i(0); i<m_inflow_function.size(); i++) {
      std::string prefix( a_pp.prefix() );
      prefix += "." + m_bdry_name[i];
      ParmParse fpp( prefix.c_str() );
      std::string function_name;
      fpp.query( "function", function_name );
      
      if (function_name.length()==0){ //make radial_inner or radial_lower compatible
        if( m_bdry_name[i].compare("radial_outer")==0 ){
            cout<<"m_bdry_name["<<i<<"] = \"radial_outer\" not found, trying \"radial_upper\""<<endl;
            prefix.replace(prefix.end()-12,prefix.end(),"radial_upper");
            ParmParse fpp( prefix.c_str() );
            fpp.query( "function", function_name );
        }
        else if ( m_bdry_name[i].compare("radial_upper")==0 ){
            cout<<"m_bdry_name["<<i<<"] = \"radial_upper\" not found, trying \"radial_outer\""<<endl;
            prefix.replace(prefix.end()-12,prefix.end(),"radial_outer");
            ParmParse fpp( prefix.c_str() );
            fpp.query( "function", function_name );
        }
        else if ( m_bdry_name[i].compare("radial_inner")==0 ){
            cout<<"m_bdry_name["<<i<<"] = \"radial_inner\" not found, trying \"radial_lower\""<<endl;
            prefix.replace(prefix.end()-12,prefix.end(),"radial_lower");
            ParmParse fpp( prefix.c_str() );
            fpp.query( "function", function_name );
        }
        else if ( m_bdry_name[i].compare("radial_lower")==0 ){
            cout<<"m_bdry_name["<<i<<"] = \"radial_lower\" not found, trying \"radial_inner\""<<endl;
            prefix.replace(prefix.end()-12,prefix.end(),"radial_inner");
            ParmParse fpp( prefix.c_str() );
            fpp.query( "function", function_name );
        } 
      }

      m_inflow_function[i] = library->find( function_name );
   }

   if (m_verbosity) {
      printParameters();
   }
}


#include "NamespaceFooter.H"

