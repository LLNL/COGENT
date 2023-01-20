#ifndef _APAROPS_H_
#define _APAROPS_H_

#include "EMVorticityOp.H"


#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "GKTimeIntegration.H"
#include "KineticSpecies.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM


#include "NamespaceHeader.H"
namespace PS = PS_NAMESPACE;


/// AparOps interface.
/**
 * AparOps interface
 */
class AparOps
{
public:
   
   /// Constructor.
   /**
       */
   AparOps( const string&   a_pp_prefix,
            const MagGeom&  a_geometry );

   /// Destructor.
   /**
    */
   virtual ~AparOps() {}


   void defineAparData( EMFields&  a_EM_fields )
   {
      IntVect ghost_vect = IntVect::Unit;

      LevelData<FluxBox>&  Apar_face = a_EM_fields.getAparFace();
      Apar_face.define(m_geometry.gridsFull(), 1, ghost_vect);

      LevelData<FArrayBox>&  Apar_cell = a_EM_fields.getAparCell();
      Apar_cell.define(m_geometry.gridsFull(), 1, ghost_vect);

      LevelData<FluxBox>&  Apar_derivs_face = a_EM_fields.getAparDerivsFace();
      Apar_derivs_face.define(m_geometry.gridsFull(), 4, ghost_vect);

      LevelData<FArrayBox>&  Apar_derivs_cell = a_EM_fields.getAparDerivsCell();
      Apar_derivs_cell.define(m_geometry.gridsFull(), 4, ghost_vect);

      // N.B.: The injected data is defined during the injection
   }

   
   void computeAparDerivs( EMFields&                         a_EM_fields,
                           const PS::KineticSpeciesPtrVect&  a_kinetic_species_phys,
                           const LevelData<FArrayBox>&       a_Apar,
                           FluidOpInterface&                 a_fluid_op,
                           const bool                        a_inject )
   {
      // Set the cell-centered and face-centered Apar derivatives
      EMVorticityOp& op = (EMVorticityOp&)a_fluid_op;
      op.setAparDerivatives(a_kinetic_species_phys, a_Apar, a_EM_fields);

      if ( a_inject ) {
         inject(a_EM_fields.getAparDerivs(),
                a_kinetic_species_phys,
                a_EM_fields.getAparDerivsFace(),
                a_EM_fields.getAparDerivsCell());
      }
   }

protected:
   
   void inject( PS::LevelData<PS::FluxBox>&       injected_data,
                const PS::KineticSpeciesPtrVect&  kinetic_species,
                const LevelData<FluxBox>&         data_face,
                const LevelData<FArrayBox>&       data_cell );

   const MagGeom& m_geometry;

private:
   
};


#include "NamespaceFooter.H"

#endif
