#ifndef  _VORTICITYBEOP_H_
#define  _VORTICITBEYOP_H_

#include "VorticityOp.H"
#include "GKVorticityBE.H"
#include "PhiOps.H"
#include "NamespaceHeader.H"

/**
 * VorticityOp operator class.
*/
class VorticityBEOp
   : public VorticityOp
{
public:
   
  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
   VorticityBEOp(const ParmParse&   pp,
                 const MagGeom&     geometry,
                 const double       larmor,
                 const int          verbosity );

  /// Destructor.
  /**
   */
   virtual ~VorticityBEOp();

   void initializeBEModel(const PS::KineticSpeciesPtrVect&  kinetic_species_phys,
                           const FluidSpeciesPtrVect&       fluid_species_phys,
                           const Vector<Real>&              scalar_data );

   void updatePotentialBEModel(EMFields&                         a_EM_fields,
                               PhiOps&                           a_phi_ops,
                               bool                              a_cbc_or_ampere,
                               const PS::KineticSpeciesPtrVect&  a_kinetic_species_phys,
                               const FluidSpeciesPtrVect&        a_fluid_species_comp,
                               const PS::ScalarPtrVect&          a_scalar_data,
                               EllipticOpBC&                     a_potential_bc,
                               const Real                        a_dt,
                               const Real                        a_time );

   static const char* pp_name;

private:

   /// Parse parameters.
   void parseParameters( const ParmParse& pp );
   
   /// Print parameters.
   void printParameters();

   GKPoisson*      m_gyropoisson_op;
   GKVorticityBE*  m_vorticity_op;

   RefCountedPtr<EllipticOpBC> m_vorticity_op_bcs;
   RefCountedPtr<EllipticOpBC> m_gyropoisson_op_bcs;
   
   bool m_split_op;
   int m_update_pc_freq;


};


#include "NamespaceFooter.H"

#endif
