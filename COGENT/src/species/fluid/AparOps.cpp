#include "AparOps.H"

#include "NamespaceHeader.H"


AparOps::AparOps( const string&       a_pp_prefix,
                  const MagGeom&      a_geometry )
   : m_geometry(a_geometry)
{
}


void AparOps::inject(PS::LevelData<PS::FluxBox>&       a_injected_data,
                     const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                     const LevelData<FluxBox>&         a_data_face,
                     const LevelData<FArrayBox>&       a_data_cell)
{
  if (a_kinetic_species.size() > 0 ) {
      const PS::PhaseGeom& phaseGeom( a_kinetic_species[0]->phaseSpaceGeometry() );

#if 0
      phaseGeom.injectConfigurationToPhase(a_data_face,
                                              a_data_cell,
                                              a_injected_data);
      
#else
      //This code provides some performance optimization
      if (!phaseGeom.secondOrder()) {
         phaseGeom.injectConfigurationToPhase(a_data_face,
                                              a_data_cell,
                                              a_injected_data);
      }
      else {
         //Strip ghost since they are not needed for 2nd order
	 //Work with assertions in field calculations later 
         const DisjointBoxLayout& grids = m_geometry.gridsFull();
         LevelData<FArrayBox> cell_tmp(grids, 3, IntVect::Zero);
         LevelData<FluxBox> face_tmp(grids, 3, IntVect::Zero);
         for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
            cell_tmp[dit].copy(a_data_cell[dit]);
            face_tmp[dit].copy(a_data_face[dit]);
         }

         phaseGeom.injectConfigurationToPhase(face_tmp,
                                              cell_tmp,
                                              a_injected_data,
                                              false );

      }
#endif
   }
}


#include "NamespaceFooter.H"
