#ifndef  _FIELDOPINTERFACE_H_
#define  _FIELDOPINTERFACE_H_

#include "REAL.H"
#include "ParmParse.H"
#include "Field.H"
#include "FluidSpecies.H"
#include "Dimensions.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "KineticSpecies.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include "NamespaceHeader.H"
namespace PS = PS_NAMESPACE;

/**
 * fieldOp interface class.
 *
 * Abstract base class the provides an interface on which algorithm classes
 * can rely to set operators on fields.
*/
class FieldOpInterface
{
   public:

      /// Destructor.
      /**
       */
      virtual ~FieldOpInterface() {;}

      /// Evaluates field RHS.
      /**
       *  Pure virtual function that will be resposible for evaluating
       *  field RHS on the given field.
       *
       *  @param[in,out] rhs           -  contribution to be added to rhs.
       *  @param[in] fields            -  current solution for fields.
       *  @param[in] fluids            -  current solution for fluids.
       *  @param[in] kinetic_specties  -  current solution for kinetic species.
       *  @param[in] fieldVecComp      - component of the field vector to which operator is applied.
       *  @param[in] time              - the time at which the field RHS is to be evaluated
       */
      virtual void evalFieldRHS( FieldPtrVect&                      rhs,
                                 const FieldPtrVect&                fields,
                                 const FluidSpeciesPtrVect&         fluids,
                                 const PS::KineticSpeciesPtrVect&   kinetic_species,
                                 const LevelData<FluxBox>&          E_field,
                                 const int                          fieldVecComp,
                                 const Real                         time) = 0;

   
      virtual Real computeDt( const FieldPtrVect&        fields,
                              const FluidSpeciesPtrVect& fluids)
      {
         return DBL_MAX;
      }

      virtual Real TimeScale( const FieldPtrVect&        fields,
                              const FluidSpeciesPtrVect& fluids)
      {
         return DBL_MAX;
      }


   
};

#include "NamespaceFooter.H"
#endif
