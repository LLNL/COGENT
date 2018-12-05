#ifndef _NEWGKPOISSONBOLTZMANN_H_
#define _NEWGKPOISSONBOLTZMANN_H_

#include "GKPoissonBoltzmann.H"
#include "SNCoreCoordSys.H"

#include "NamespaceHeader.H"

class NewGKPoissonBoltzmann
   : public GKPoisson
{
public:

   NewGKPoissonBoltzmann( ParmParse&                  pp,
                          const MagGeom&              geom,
                          const Real                  larmor_number,
                          const Real                  debye_number,
                          const LevelData<FArrayBox>& initial_ion_charge_density );

   ~NewGKPoissonBoltzmann();

   void computePotentialAndElectronDensity(
      LevelData<FArrayBox>&        phi,
      BoltzmannElectron&           ne,
      const LevelData<FArrayBox>&  ni,
      const EllipticOpBC&          bc,
      const bool                   first_step = false);

   virtual void setOperatorCoefficients( const LevelData<FArrayBox>&  ni,
                                         const EllipticOpBC&          bc ); 

   virtual void setOperatorCoefficients( const LevelData<FArrayBox>&  ion_mass_density,
                                         const EllipticOpBC&          bc,
                                         const double                 core_outer_bv,
                                         double&                      lo_value,
                                         double&                      hi_value,
                                         LevelData<FArrayBox>&        radial_gkp_divergence_average);

   void setDivertorBVs( const LevelData<FArrayBox>&  ion_charge_density,
                        const LevelData<FArrayBox>&  ion_parallel_current_density,
                        EllipticOpBC&                bc );

   ///
   /**
      Computes sheath boundary conditions; defines the corresponsding members of the BC class; 
      computes boundary data required for the single-null Boltzmann electrons
   */
   void computeBoundaryData( FArrayBox& inner_divertor_bvs,
                             FArrayBox& outer_divertor_bvs,
                             const BoltzmannElectron& ne,
                             const LevelData<FArrayBox>& Zni,
                             const LevelData<FArrayBox>& Jpar); 

   double coreOuterPotential( const MagGeom&              a_geom,
                              const LevelData<FArrayBox>& a_phi ) const;

   void getPhiTilde( const LevelData<FArrayBox>& Zni,
                     const BoltzmannElectron&    ne,
                     LevelData<FArrayBox>&       phi_tilde );

   void getFaceCenteredFieldOnCore( const LevelData<FArrayBox>& a_phi,
                                    LevelData<FluxBox>&         a_field ) const;

   static const char* pp_name;

protected:

   int getConsequtiveRadialIndex( const int a_mapped_index,
                                  const int a_block_number ) const;

   void computeSheathEllipticOpBC( const BoltzmannElectron&  ne,
                                   const double *            jpar_outer_plate, 
                                   const double *            jpar_inner_plate,
                                   const double *            bmag_outer_plate, 
                                   const double *            bmag_inner_plate,
                                   const int                 norm_dir );

   void copyToCore( const LevelData<FArrayBox>& a_in,
                    LevelData<FArrayBox>&       a_out ) const;

   void copyFromCore( const LevelData<FArrayBox>& a_in,
                      LevelData<FArrayBox>&       a_out ) const;

   double*  m_Zni_outer_plate;
   double*  m_Zni_inner_plate;
   double*  m_phi_outer_plate;
   double*  m_phi_inner_plate;

   RefCountedPtr<MagCoordSys> m_core_coord_sys;
   MagGeom* m_core_geometry;

   SNCoreEllipticOpBC* m_bc_core;

   GKPoissonBoltzmann* m_gkpb_solver;
};

#include "NamespaceFooter.H"

#endif   /* _NEWGKPOISSONBOLTZMANN_H_ */
