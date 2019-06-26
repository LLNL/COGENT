#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "FloorPostProcessor.H"

#include "Box.H"
#include "BoxIterator.H"
#include "BoxLayoutData.H"
#include "DataIterator.H"
#include "DisjointBoxLayout.H"
#include "Interval.H"

#include "NamespaceHeader.H"

FloorPostProcessor::FloorPostProcessor(
   const Real&    a_value,
   const bool&    a_absolute,
   const bool&    a_verbose)
{
   define( a_value, a_absolute, a_verbose);
}


void FloorPostProcessor::define(
   const Real&    a_value,
   const bool&    a_absolute,
   const bool&    a_verbose)
{
   CH_assert( a_value >= 0 );
   m_value = a_value;
   m_absolute = a_absolute;
   m_verbose = a_verbose;
}


void FloorPostProcessor::enforce( FluidSpeciesPtrVect& a_data ) const
{

   Real min_value;
   Real max_value;
   
   for (int species_id(0); species_id<a_data.size(); species_id++) {
      FluidSpecies& species( static_cast<FluidSpecies&>(*(a_data[species_id])) );
      for (int n=0; n<species.num_cell_vars(); ++n) {
         LevelData<FArrayBox>& fld( species.cell_var(n) );
         const DisjointBoxLayout& grids = fld.getBoxes();
         
         getMinMax(fld, min_value, max_value);
         
         if (m_absolute && min_value < m_value) {
            Real delta = m_value - min_value;
            for (DataIterator dit(grids); dit.ok(); ++dit) {
               fld[dit].plus(delta);
            }
         }

         if (!m_absolute && min_value/max_value < m_value) {
            Real delta = (m_value*max_value - min_value)/(1.0 - m_value);
            for (DataIterator dit(grids); dit.ok(); ++dit) {
               fld[dit].plus(delta);
            }
         }

      } //end loop over species components
   } // end loop over species

   //pout() << "  Look Ma!  No hands!" << endl;
}


void
FloorPostProcessor::getMinMax(const LevelData<FArrayBox>& a_data, double& a_min, double& a_max) const
{
   const DisjointBoxLayout & grids = a_data.disjointBoxLayout();
   
   a_min =  DBL_MAX;
   a_max = -DBL_MAX;
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      const FArrayBox& this_data = a_data[dit];
      Box box(grids[dit]);
      if (this_data.min(box) < a_min) a_min = this_data.min(box);
      if (this_data.max(box) > a_max) a_max = this_data.max(box);
   }
   
#ifdef CH_MPI
   double local_min = a_min;
   MPI_Allreduce(&local_min, &a_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
   double local_max = a_max;
   MPI_Allreduce(&local_max, &a_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
}


#include "NamespaceFooter.H"
