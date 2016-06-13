#ifndef  _VORTICITY_H_
#define  _VORTICITY_H_

#include "REAL.H"
#include "FieldOpInterface.H"
#include "GKFieldOp.H"
#include "ParmParse.H"
#include <sstream>

#include "GridFunction.H"
#include "GridFunctionLibrary.H"

#include "NamespaceHeader.H"

/**
 * Krook collisions operator class.
*/
class Vorticity
   : public FieldOpInterface
{
public:

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  Vorticity( ParmParse&  a_pp,
             const int a_verbosity );

  /// Destructor.
  /**
   */
  virtual ~Vorticity();

   /// Evaluates field RHS.
   /**
    *  Pure virtual function that will be resposible for evaluating
    *  field RHS on the given field.
    *
    *  @param[in,out] rhs           -  contribution to be added to rhs.
    *  @param[in] fields            -  current solution for fields.
    *  @param[in] fluids            -  current solution for fluids.
    *  @param[in] kinetic_specties  -  current solution for kinetic species.
    *  @param[in] fieldVecComp      -  component of the field vector to which operator is applied.
    *  @param[in] time              -  the time at which the field RHS is to be evaluated
    */
   virtual void evalFieldRHS( FieldPtrVect&                     rhs,
                             const FieldPtrVect&                fields,
                             const FluidSpeciesPtrVect&         fluids,
                             const PS::KineticSpeciesPtrVect&   kinetic_species,
                             const LevelData<FluxBox>&          E_field,
                             const int                          fieldVecComp,
                             const Real                         time);
   
   /// Compute a stable time step.
   /**
    * Computes and returns an estimate of the maximum stable time step.
    */
   Real computeDt( const FieldPtrVect&         fields,
                   const FluidSpeciesPtrVect&  fluids);
   

private:

   //Computes perpendicular current density
   void computePerpCurrentDensity( LevelData<FluxBox>&                perp_current_density,
                                  const PS::KineticSpeciesPtrVect&    species,
                                  const LevelData<FluxBox>&           field) const;

   //Computes parallel current density
   void computeParallelCurrentDensity( LevelData<FluxBox>&               parallel_current_density,
                                       const PS::KineticSpeciesPtrVect&    species,
                                       const LevelData<FluxBox>&           field) const;
   
   //Computes ion charge density
   void computeIonChargeDensity( LevelData<FArrayBox>&                ion_charge_density,
                                const PS::KineticSpeciesPtrVect&      species ) const;

   //Second order version of the 4th order fourthOrderCellToFace
   void cellToFace(LevelData<FluxBox>& faceData,
                   const LevelData<FArrayBox>& cellData) const;

   //Second order extrapolation across physical boundaries
   void extrapolateAtDomainBnd(const MagGeom&  mag_geom,
                               LevelData<FArrayBox>& data) const;
   
   void setZero( LevelData<FArrayBox>& data ) const;
   
   int m_verbosity;
   bool m_first_step;
   int m_num_ghosts;
   double m_conductivity;
   double m_Te;
   
   /// Parse parameters.
   /**
    * Private method to obtain control parameters from "field_model.field_name" section
    * of ParmParse database, where "field_name" is replaced by the field name (e.g., field_model.VectPotential)
    *
    * \param[in] ppcls field_model.field input database.
    */
   void parseParameters( ParmParse& pp );
   
   /// Print parameters.
   /**
    * Private method to print parameters.
    */
   void printParameters();
   
};


#include "NamespaceFooter.H"

#endif
