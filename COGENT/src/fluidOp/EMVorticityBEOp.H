#ifndef  _EMVORTICITYBEOP_H_
#define  _EMVORTICITBEYOP_H_

#include "EMVorticityOp.H"
#include "GKEMVorticity.H"
#include "PhiOps.H"
#include "NamespaceHeader.H"

/**
 * VorticityOp operator class.
*/
class EMVorticityBEOp
   : public EMVorticityOp
{
public:
   
  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
   EMVorticityBEOp(const ParmParse&   pp,
                   const MagGeom&     geometry,
                   const double       larmor,
                   const double       ion_skin_depth,
                   const double       electron_skin_depth,
                   const int          verbosity );

  /// Destructor.
  /**
   */
   virtual ~EMVorticityBEOp();


   void initializeBEModel(const PS::KineticSpeciesPtrVect&  kinetic_species_phys,
                          const FluidSpeciesPtrVect&        fluid_species_phys,
                          const Vector<Real>&               scalar_data );

   void updatePotentialBEModel(EMFields&                         EM_fields,
                               PhiOps&                           phi_ops,
                               bool                              cbc_or_ampere,
                               const PS::KineticSpeciesPtrVect&  kinetic_species_phys,
                               const FluidSpeciesPtrVect&        fluid_species_comp,
                               const PS::ScalarPtrVect&          scalar_data,
                               EllipticOpBC&                     potential_bc,
                               const Real                        dt,
                               const Real                        time );
   
   static const char* pp_name;

private:

   
   void updateBEModelData(const PS::KineticSpeciesPtrVect&  kinetic_species_phys,
                           const FluidSpeciesPtrVect&       fluid_species_comp,
                           const Vector<Real>&              scalar_data );

   
   void advance(LevelData<FArrayBox>&   phi,
                LevelData<FArrayBox>&   Apar,
                LevelData<FArrayBox>&   rhs_phi,
                LevelData<FArrayBox>&   rhs_Apar,
                const Real              dt);
   
   /// Parse parameters.
   void parseParameters( const ParmParse& pp );
   
   /// Print parameters.
   void printParameters();
   
   GKPoisson*     m_gyropoisson_op;
   GKEMVorticity* m_Apar_soln_op;
   
   RefCountedPtr<EllipticOpBC> m_gyropoisson_op_bcs;
   RefCountedPtr<EllipticOpBC> m_Apar_soln_op_bcs;

   LevelData<FArrayBox> m_Apar;
   
   Real m_dt_sub_cycle;
   int m_update_pc_freq;

};


#include "NamespaceFooter.H"

#endif
