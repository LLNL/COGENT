#ifndef _GKVORTICITYBE_H_
#define _GKVORTICITYBE_H_

#include "GKVorticity.H"

#include "NamespaceHeader.H"
namespace PS = PS_NAMESPACE;
namespace VEL = VEL_NAMESPACE;

/// Gyrokintic Vorticity operator class.
/**
 * This operator class solves the mapped gyrokinetic Vorticity equation:
 *
 * \f[
 *  {\nabla} \cdot \left ( N^T D N/J  {\nabla} \Phi \right ) = \rho,
 * \f]
 * where
 * \f{eqnarray*}
 * {\bf D} & \equiv & (De)^2  \mathbf{I} + \frac{(La)^2}{ B^2 }
 *  \left ( \sum_i Z_i m_i \bar{n}_i \right ) \left ( \mathbf{I} -
 *  \mathbf{b}\mathbf{b}^T \right ),\\  
 * \rho & \equiv &  n_e - \sum_i Z_i \bar{n}_i.
 * \f}
*/
class GKVorticityBE
   : public GKVorticity
{
public:

   /// Constructor with initialization.
   /**
    *  @param[in] pp the input database.
    *  @param[in] geom the geometry data.
    *  @param[in] larmor_number the Larmor number for the operator.
    *  @param[in] debye_number the Debye number for the operator.
    */
   
   GKVorticityBE(const ParmParse&   pp,
                 const ParmParse&   pp_base,
                 const MagGeom&     geom,
                 const Real         larmor_number,
                 const Real         debye_number,
                 const bool         second_order,
                 const bool         low_pollution,
                 const bool         include_pol_den_correction,
                 const bool         include_diffusion,
                 const std::string& model = "Vorticity" );
      
   /// Destructor.
   /**
    */
   virtual ~GKVorticityBE();

   
   virtual void applyOp(LevelData<FArrayBox>& out,
                        const LevelData<FArrayBox>& in,
                        bool homogeneous = false);
   
   virtual void computeBcDivergence( LevelData<FArrayBox>& out );
   
   virtual void setVorticityOperatorCoefficients( const LevelData<FArrayBox>&  ion_mass_density,
                                                 const LevelData<FArrayBox>&  ion_charge_density,
                                                 const LevelData<FArrayBox>&  electron_temeprature,
                                                 EllipticOpBC&                bc,
                                                 const bool                   update_preconditioner );

   void applyHighOrderCorrectionOp(LevelData<FArrayBox>&        out,
                                   const LevelData<FArrayBox>&  in,
                                   bool                         homogeneous );

   void  computeFluxDivergenceWithCoeff(LevelData<FArrayBox>&       out,
                                        const LevelData<FArrayBox>& in,
                                        const LevelData<FluxBox>&   coeff_mapped,
                                        const bool                  homogeneous_bcs,
                                        const bool                  subtract_fs_par_div,
                                        const bool                  extrap_to_ghosts = false);
      
   /// parseParameters function
   /**
    * parseParameters function
    */
   void parseParameters( const ParmParse& pp );

   /// printParameters function
   /**
    * printParameters function
    */
   void printParameters();

   /// ParmParse object name
   /**
    * ParmParse object name.  Used to locate input file parameters of interest to
    * this class.
    */
   static const char* pp_name;
   

private:
   
   void verifyMatrix( const MBSolver* matrix );
   
   RefCountedPtr<EllipticOpBC> m_potential_bc;
};



#include "NamespaceFooter.H"

#endif
