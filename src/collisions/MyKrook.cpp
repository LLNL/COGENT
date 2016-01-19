#include <math.h>

#include "MyKrook.H"

#include "FourthOrderUtil.H"
#include "Directions.H"
#include "EdgeToCell.H"
#include "PhaseGeom.H"
#include "PhaseBlockCoordSys.H"
#include "CollisionsF_F.H"
#include "KineticFunctionLibrary.H"
#include "ConstFact.H"
#include "inspect.H"

#include "NamespaceHeader.H" 


MyKrook::MyKrook( ParmParse& a_ppcls, const int a_verbosity )
   : m_verbosity(a_verbosity),
     m_first_step(true),
     m_cls_freq(-1.0)
{
   parseParameters( a_ppcls );
   if (m_verbosity>0) {
      printParameters();
   }
}


MyKrook::~MyKrook()
{
}


void MyKrook::evalClsRHS( KineticSpeciesPtrVect&       a_rhs,
                          const KineticSpeciesPtrVect& a_soln,
                          const int                    a_species,
                          const Real                   a_time )
// NB: a_soln is on the computational grid
{
   // Get solution distribution function (J*Bstar_par*dfn) for the current species
   const KineticSpecies& soln_species( *(a_soln[a_species]) );
   const LevelData<FArrayBox>& soln_dfn( soln_species.distributionFunction() );

   // Create reference (J*Bstar_par*dfn_init) distribution
   KineticSpeciesPtr ref_species( soln_species.clone( IntVect::Zero, false ) );
   m_ref_func->assign( *ref_species, a_time );
   const LevelData<FArrayBox>& init_dfn( ref_species->distributionFunction() );
   
   static Real delta_F_0(0.0);
   static Real ref_val(0.0);

   if (m_first_step) {

      if (procID()==0) {
         const IntVect location( IntVect::Zero );
         DataIterator dit = soln_dfn.dataIterator();
         const FArrayBox& soln_data( soln_dfn[dit] );
         delta_F_0 = soln_data.get( location, 0 );
         const FArrayBox& ref_data( init_dfn[dit] );
         ref_val = ref_data.get( location, 0 );
         delta_F_0 -= ref_val;
      }
   
   }

   if (procID()==0) {
      const IntVect location( IntVect::Zero );
      DataIterator dit = soln_dfn.dataIterator();
      const FArrayBox& soln_data( soln_dfn[dit] );
      Real soln_val = soln_data.get( location, 0 );
      Real ext_val = ref_val + delta_F_0 * exp( -a_time * m_cls_freq );
      cout << "MyKrook: time, exact, sol, error, ref = "
           << std::scientific 
           << a_time << "\t"
           << ext_val << "\t"
           << soln_val << "\t"
           << (ext_val-soln_val)
           << "\t" << ref_val
           << endl;
//      cout << "MyKrook: time, error = " << a_time << "\t" << ext_val-soln_val << endl;
   }

   // Compute the difference from the reference (or initial) solution
//   const DisjointBoxLayout& grids( soln_dfn.getBoxes() );
//   const int n_comp( soln_dfn.nComp() );
//   LevelData<FArrayBox> delta_F( grids, n_comp, IntVect::Zero );
   KineticSpecies& rhs_species( *(a_rhs[a_species]) );
   LevelData<FArrayBox>& delta_F( rhs_species.distributionFunction() );
   
   for (DataIterator sdit(soln_dfn.dataIterator()); sdit.ok(); ++sdit) {
      delta_F[sdit].copy( soln_dfn[sdit] );
      delta_F[sdit].minus( init_dfn[sdit] );
      delta_F[sdit].mult( -m_cls_freq );
   }
   
   m_first_step = false;
}


inline
void MyKrook::parseParameters( ParmParse& a_ppcls )
{
   a_ppcls.query( "cls_freq", m_cls_freq );

   KineticFunctionLibrary* library = KineticFunctionLibrary::getInstance();
   std::string function_name;
   a_ppcls.query( "ref_function", function_name );
   m_ref_func = library->find( function_name );
}


inline
void MyKrook::printParameters()
{
   if (procID()==0) {
      std::cout << "MyKrook collisions parameters:" << std::endl;
      std::cout << "  cls_freq  =  " << m_cls_freq << std::endl;
      std::cout << "  Reference Function:" << std::endl;
      m_ref_func->printParameters();
   }
}


Real MyKrook::computeDt(const KineticSpeciesPtrVect& soln)
{
   return 1.0/m_cls_freq;
}


void MyKrook::addReferenceDfn( KineticSpecies& a_result,
                               const Real      a_time, 
                               const Real      a_scale )
{
   LevelData<FArrayBox>& result_dfn( a_result.distributionFunction() );
         
   KineticSpeciesPtr ref_species( a_result.clone( IntVect::Zero, false ) );
   m_ref_func->assign( *ref_species, a_time );
   const LevelData<FArrayBox>& ref_dfn( ref_species->distributionFunction() );
   for (DataIterator dit(result_dfn.dataIterator()); dit.ok(); ++dit) {
      result_dfn[dit].plus( ref_dfn[dit], a_scale );
   }
}

#include "NamespaceFooter.H"
