#ifndef _EMFIELDS_H_
#define _EMFIELDS_H_

#include "LevelData.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "LevelData.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include "NamespaceHeader.H"
namespace PS = PS_NAMESPACE;

class EMFields
{
public:

   // Constructor
   EMFields()
      : m_Efield_comps(3),
        m_Apar_derivs_comps(4)
   {}


   // Copy constructor
   EMFields( const EMFields& a_EM_fields )
   {
      const LevelData<FArrayBox>& Efield_cell = a_EM_fields.m_Efield_cell;
      if ( Efield_cell.isDefined() ) {
         m_Efield_cell.define( Efield_cell.disjointBoxLayout(), Efield_cell.nComp(), Efield_cell.ghostVect());
         for (DataIterator dit(Efield_cell.dataIterator()); dit.ok(); ++dit) {
            m_Efield_cell[dit].copy(Efield_cell[dit]);
         }
      }

      const LevelData<FluxBox>& Efield_face = a_EM_fields.m_Efield_face;
      if ( Efield_face.isDefined() ) {
         m_Efield_face.define( Efield_face.disjointBoxLayout(), Efield_face.nComp(), Efield_face.ghostVect());
         for (DataIterator dit(Efield_face.dataIterator()); dit.ok(); ++dit) {
            m_Efield_face[dit].copy(Efield_face[dit]);
         }
      }

      const PS::LevelData<PS::FluxBox>& Efield_injected = a_EM_fields.m_Efield_injected;
      if ( Efield_injected.isDefined() ) {
         m_Efield_injected.define( Efield_injected.disjointBoxLayout(), Efield_injected.nComp(), Efield_injected.ghostVect());
         for (PS::DataIterator dit(Efield_injected.dataIterator()); dit.ok(); ++dit) {
            m_Efield_injected[dit].copy(Efield_injected[dit]);
         }
      }

      const LevelData<FArrayBox>& phi = a_EM_fields.m_phi;
      if ( phi.isDefined() ) {
         m_phi.define( phi.disjointBoxLayout(), phi.nComp(), phi.ghostVect());
         for (DataIterator dit(phi.dataIterator()); dit.ok(); ++dit) {
            m_phi[dit].copy(phi[dit]);
         }
      }

      const LevelData<FArrayBox>& phi_node = a_EM_fields.m_phi_node;
      if ( phi_node.isDefined() ) {
         m_phi_node.define( phi_node.disjointBoxLayout(), phi_node.nComp(), phi_node.ghostVect());
         for (DataIterator dit(phi_node.dataIterator()); dit.ok(); ++dit) {
            m_phi_node[dit].copy(m_phi_node[dit]);
         }
      }
   }


   void compare( const EMFields& a_EM_fields )
   {
      bool differences = false;
      
      const LevelData<FArrayBox>& Efield_cell = a_EM_fields.m_Efield_cell;
      if ( Efield_cell.isDefined() ) {
         LevelData<FArrayBox> diff( Efield_cell.disjointBoxLayout(), Efield_cell.nComp(), Efield_cell.ghostVect());
         for (DataIterator dit(Efield_cell.dataIterator()); dit.ok(); ++dit) {
            diff[dit].copy(Efield_cell[dit]);
            diff[dit] -= m_Efield_cell[dit];
         }
         double diff_norm = SpaceUtils::MaxNorm(diff);
         if ( diff_norm != 0. ) {
            differences = true;
            if ( procID() == 0 ) {
               cout << "EMFields::compare(): Cell-centered field is different, diff_norm = " << diff_norm << endl;
            }
         }
      }

      const LevelData<FluxBox>& Efield_face = a_EM_fields.m_Efield_face;
      if ( Efield_face.isDefined() ) {
         LevelData<FluxBox> diff( Efield_face.disjointBoxLayout(), Efield_face.nComp(), Efield_face.ghostVect());
         for (DataIterator dit(Efield_face.dataIterator()); dit.ok(); ++dit) {
            diff[dit].copy(Efield_face[dit]);
            diff[dit] -= m_Efield_face[dit];
         }
         double diff_norm = SpaceUtils::MaxNorm(diff);
         if ( diff_norm != 0. ) {
            differences = true;
            if ( procID() == 0 ) {
               cout << "EMFields::compare(): Face-centered field is different, diff_norm = " << diff_norm << endl;
            }
         }
      }

      const PS::LevelData<PS::FluxBox>& Efield_injected = a_EM_fields.m_Efield_injected;
      if ( Efield_injected.isDefined() && m_Efield_injected.isDefined() ) {
         PS::LevelData<PS::FluxBox> diff( Efield_injected.disjointBoxLayout(), Efield_injected.nComp(), Efield_injected.ghostVect());
         for (PS::DataIterator dit(Efield_injected.dataIterator()); dit.ok(); ++dit) {
            diff[dit].copy(Efield_injected[dit]);
            diff[dit] -= m_Efield_injected[dit];
         }
         double diff_norm = PS::SpaceUtils::MaxNorm(diff);
         if ( diff_norm != 0. ) {
            differences = true;
            if ( procID() == 0 ) {
               cout << "EMFields::compare(): Injected field is different, diff_norm = " << diff_norm << endl;
            }
         }
      }

      const LevelData<FArrayBox>& phi = a_EM_fields.m_phi;
      if ( phi.isDefined() ) {
         LevelData<FArrayBox> diff( phi.disjointBoxLayout(), phi.nComp(), phi.ghostVect());
         for (DataIterator dit(phi.dataIterator()); dit.ok(); ++dit) {
            diff[dit].copy(phi[dit]);
            diff[dit] -= m_phi[dit];
         }
         double diff_norm = SpaceUtils::MaxNorm(diff);
         if ( diff_norm != 0. ) {
            differences = true;
            if ( procID() == 0 ) {
               cout << "EMFields::compare(): potential is different, diff_norm = " << diff_norm << endl;
            }
         }
      }

      const LevelData<FArrayBox>& phi_node = a_EM_fields.m_phi_node;
      if ( phi_node.isDefined() ) {
         LevelData<FArrayBox> diff( phi_node.disjointBoxLayout(), phi_node.nComp(), phi_node.ghostVect());
         for (DataIterator dit(phi_node.dataIterator()); dit.ok(); ++dit) {
            diff[dit].copy(m_phi_node[dit]);
            diff[dit] -= m_phi_node[dit];
         }
         double diff_norm = SpaceUtils::MaxNorm(diff);
         if ( diff_norm != 0. ) {
            differences = true;
            if ( procID() == 0 ) {
               cout << "EMFields::compare(): Nodal potential is different, diff_norm = " << diff_norm << endl;
            }
         }
      }

      if ( differences ) {
         exit(1);
      }
   }

   void setPhi( const LevelData<FArrayBox>& a_phi ) {
      CH_assert(m_phi.isDefined());

      for (DataIterator dit(m_phi.dataIterator()); dit.ok(); ++dit) {
         m_phi[dit].copy(a_phi[dit]);
      }
   }
   
   void setSecondOrder( bool second_order ) {m_second_order = second_order;}

   
   // E field accessors
   const PS::LevelData<PS::FluxBox>& getEField() const { return m_Efield_injected;}
   const LevelData<FluxBox>& getEFieldFace() const { return m_Efield_face;}
   const LevelData<FArrayBox>& getEFieldCell() const { return m_Efield_cell;}
   PS::LevelData<PS::FluxBox>& getEField() { return m_Efield_injected;}
   LevelData<FluxBox>& getEFieldFace() { return m_Efield_face;}
   LevelData<FArrayBox>& getEFieldCell() { return m_Efield_cell;}

   // Potential accessors
   const LevelData<FArrayBox>& getPhi() const { return m_phi; }
   LevelData<FArrayBox>& getPhi() { return m_phi; }

   // Nodal potential accessors
   const LevelData<FArrayBox>& getPhiNode() const { return m_phi_node; }
   LevelData<FArrayBox>& getPhiNode() { return m_phi_node; }

   // A_parallel accessors
   const PS::LevelData<PS::FluxBox>& getApar() const { return m_Apar_injected;}
   const LevelData<FluxBox>& getAparFace() const { return m_Apar_face;}
   const LevelData<FArrayBox>& getAparCell() const { return m_Apar_cell;}
   PS::LevelData<PS::FluxBox>& getApar() { return m_Apar_injected;}
   LevelData<FluxBox>& getAparFace() { return m_Apar_face;}
   LevelData<FArrayBox>& getAparCell() { return m_Apar_cell;}

   // A_parallel derivatives accessors
   const PS::LevelData<PS::FluxBox>& getAparDerivs() const { return m_Apar_derivs_injected;}
   const LevelData<FluxBox>& getAparDerivsFace() const { return m_Apar_derivs_face;}
   const LevelData<FArrayBox>& getAparDerivsCell() const { return m_Apar_derivs_cell;}
   PS::LevelData<PS::FluxBox>& getAparDerivs() { return m_Apar_derivs_injected;}
   LevelData<FluxBox>& getAparDerivsFace() { return m_Apar_derivs_face;}
   LevelData<FArrayBox>& getAparDerivsCell() { return m_Apar_derivs_cell;}

   bool secondOrder() const {return m_second_order;}

protected:

   PS::LevelData<PS::FluxBox> m_Efield_injected;
   LevelData<FArrayBox> m_Efield_cell;
   LevelData<FluxBox> m_Efield_face;
   int m_Efield_comps;

   LevelData<FArrayBox> m_phi;
   LevelData<FArrayBox> m_phi_node;

   PS::LevelData<PS::FluxBox> m_Apar_injected;
   LevelData<FArrayBox> m_Apar_cell;
   LevelData<FluxBox> m_Apar_face;

   PS::LevelData<PS::FluxBox> m_Apar_derivs_injected;
   LevelData<FArrayBox> m_Apar_derivs_cell;
   LevelData<FluxBox> m_Apar_derivs_face;
   int m_Apar_derivs_comps;

   bool m_second_order;
};

#include "NamespaceFooter.H"
#endif
