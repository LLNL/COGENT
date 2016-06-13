#ifndef  _FLUIDOPINTERFACE_H_
#define  _FLUIDOPINTERFACE_H_

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
 * FluidOp interface class.
 *
 * Abstract base class the provides an interface on which algorithm classes
 * can rely to set operators on fluid species.
*/
class FluidOpInterface
{
   public:

      /// Destructor.
      /**
       */
      virtual ~FluidOpInterface() {;}

      /// Evaluates fluid RHS.
      /**
       *  Pure virtual function that will be resposible for evaluating
       *  fluid RHS on the given species.
       *
       *  @param[in,out] rhs           -  contribution to be added to rhs.
       *  @param[in] fields            -  current solution for fields.
       *  @param[in] fluids            -  current solution for fluids.
       *  @param[in] kinetic_specties  -  current solution for kinetic species.
       *  @param[in] species           -  fluid species to which operator is applied.
       *  @param[in] time              -  the time at which the fluid RHS is to be evaluated.
       */
      virtual void evalFluidRHS( FluidSpeciesPtrVect&               rhs,
                                 const FieldPtrVect&                fields,
                                 const FluidSpeciesPtrVect&         fluids,
                                 const PS::KineticSpeciesPtrVect&   kinetic_species,
                                 const LevelData<FluxBox>&          E_field,
                                 const int                          species,
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
