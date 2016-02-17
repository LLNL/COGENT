#ifndef _FIELDSOLVER_H_
#define _FIELDSOLVER_H_

#include <iostream>

#include "MagGeom.H"
#include "PotentialBC.H"
#include "MultiBlockLevelExchangeCenter.H"

#include "NamespaceHeader.H"

class FieldSolver
{
public:

   /// Constructor with initialization.
   /**
    *  @param[in] pp the input database.
    */
   FieldSolver(const ParmParse& pp, const MagGeom& geom);
      
   /// Destructor.
   /**
    */
   virtual ~FieldSolver();

   /// Computes the potential given the species.
   /**
    * Computes the pointwise, cell-centered potential field given the
    * species.  Each species is assumed to be capable of providing a
    * cell-averaged charge density.  If one of the species is a
    * BoltzmannElectron, the GKBoltzmannRelation is used to provide a
    * suitable electron charge density.
    *
    * @param[in] species vector of all SpeciesModels.
    * @param[out] phi the pointwise, cell-centered potential.
    */
   void computePotential( LevelData<FArrayBox>&       phi,
                          const PotentialBC&          bc,
                          const LevelData<FArrayBox>& charge_density );

   /// Computes the electric field given the potential field.
   /**
    * Computes the cell-centered physical-space components of the electric
    * field given the pointwise, cell-centered potential field.
    *
    * @param[in] phi the pointwise, cell-centered potential field.
    * @param[out] E_field the cell-centered electric field.
    */
   void computeElectricField( const LevelData<FArrayBox>& phi,
                              const PotentialBC&          bc,
                              LevelData<FArrayBox>& E_field ) const;

   /// Computes the electric field given the potential field.
   /**
    * Computes the face-centered physical-space components of the electric
    * field given the pointwise, cell-centered potential field.
    *
    * @param[in] phi the pointwise, cell-centered potential field.
    * @param[out] E_field the face_centered electric field.
    */
   void computeElectricField( const LevelData<FArrayBox>& phi,
                              const PotentialBC&          bc,
                              LevelData<FluxBox>&         E_field ) const;

   void computeFaceCenteredField( const LevelData<FArrayBox>& phi,
                                  const PotentialBC&          bc,
                                  LevelData<FluxBox>&         field,
                                  const bool                  homogeneousBCs = false ) const;

   void computeFaceCenteredMappedField( const LevelData<FArrayBox>& phi,
                                        const PotentialBC&          bc,
                                        LevelData<FluxBox>&         field,
                                        const bool                  homogeneousBCs = false ) const;

   void accumFieldStencil(const IntVect                           index,
                          const int                               dir,
                          const int                               side,
                          const int                               s,
                          const Tuple<BlockBoundary, 2*SpaceDim>& block_boundaries,
                          const int*                              block_bc_types,
                          const RealVect&                         dx,
                          const Box&                              domainBox,
                          const bool                              fourthOrder,
                          const bool                              dropOrder,
                          FArrayBox&                              stencil_values) const;

   void accumMappedBoundaryField( const PotentialBC& bc,
                                  LevelData<FluxBox>& field ) const;

   void setBcDivergence( const PotentialBC& bc,
                         LevelData<FArrayBox>& out ) const;

   /// Computes the operator (normalized by the shell volume) in the lHS of Eq. (5) in the "Poloidal_variations_PB_model-II" note.
   /**
    * Computes <grad(psi) * D *NJinv * a_in>/(one-cell wide shell volume) (a_in corresponds to mapped_grad(Phi) in the notes)
    * LevelData<FArrayBox> a_out does not have ghost cells, and stores the data in the radial direction  as follows
    * a_out[i] = <grad(psi) * D *NJinv* a_in[i]>/shell_volume[i]
    * a_lo_value = a_out[0]=<grad(psi) * D *NJinv* a_in[0]>/shell_volume[0]
    * a_hi_value = <grad(psi) * D *NJinv* a_in[nrad+1]>/shell_volume[nrad]
    * n_rad number of the core cells in the radial direction 
    * <...> denotes the shell-average
    */

   void computeRadialFSAverage( const LevelData<FluxBox>& a_in,
                                double&                   a_lo_value,
                                double&                   a_hi_value,
                                LevelData<FArrayBox>&     a_out ) const;
    
   void solve( const LevelData<FArrayBox>& rhs,
               const PotentialBC&          bc, 
               LevelData<FArrayBox>&       solution );

   void applyOperator( const LevelData<FArrayBox>& in,
                       const PotentialBC&          bc,
                       LevelData<FArrayBox>&       out,
                       const bool                  homogeneous_bcs = false ) const;

   int PCGSolve( const LevelData<FArrayBox>& rhs,
                 const PotentialBC&          bc,
                 LevelData<FArrayBox>&       solution );

   int OperatorBasedPCGSolve( const LevelData<FArrayBox>& rhs,
                              const PotentialBC&          bc,
                              LevelData<FArrayBox>&       solution );

   int GMRESSolve( const LevelData<FArrayBox>& rhs,
                   const PotentialBC&          bc,
                   LevelData<FArrayBox>&       solution );

   int OperatorBasedGMRESSolve( const LevelData<FArrayBox>& rhs,
                                const PotentialBC&          bc,
                                LevelData<FArrayBox>&       solution );

   virtual void solvePreconditioner( const LevelData<FArrayBox>& in, LevelData<FArrayBox>& out ) = 0;

   virtual void multiplyUnmappedCoefficients( LevelData<FluxBox>& a_data ) const = 0;

   void addBcDivergence(LevelData<FArrayBox>& rhs) const;

   void subtractBcDivergence(LevelData<FArrayBox>& rhs) const;

   void fillInternalGhosts( LevelData<FArrayBox>& phi ) const;

#if 0
   void computeElectricFieldExtrap( const LevelData<FArrayBox>& phi,
                                    const PotentialBC&          bc,
                                    LevelData<FArrayBox>& E_field ) const;

   void computeElectricFieldExtrap( const LevelData<FArrayBox>& phi,
                                    const PotentialBC&          bc,
                                    LevelData<FluxBox>&         E_field ) const;

   void computeFaceCenteredFieldExtrap( const LevelData<FArrayBox>& phi,
                                        const PotentialBC&          bc,
                                        LevelData<FluxBox>&         field,
                                        const bool                  homogeneousBCs = false ) const;

   void computeFaceCenteredMappedFieldExtrap( const LevelData<FArrayBox>& phi,
                                              const PotentialBC&          bc,
                                              LevelData<FluxBox>&         field,
                                              const bool                  homogeneousBCs = false ) const;

   void setMappedBoundaryFieldExtrap( LevelData<FluxBox>& field ) const;
#endif

   void computeElectricField( const LevelData<FArrayBox>& a_phi,
                              LevelData<FArrayBox>&       a_E_field ) const;

   void computeElectricField( const LevelData<FArrayBox>& a_phi,
                              LevelData<FluxBox>&         a_E_field ) const;

   void computeMappedElectricField( const LevelData<FArrayBox>& a_phi,
                                    LevelData<FArrayBox>&       a_E_field ) const;

   void computeMappedElectricField( const LevelData<FArrayBox>& a_phi,
                                    LevelData<FluxBox>&         a_E_field ) const;

   void getUnmappedField( const LevelData<FArrayBox>& a_mapped_E_field,
                          LevelData<FArrayBox>&       a_E_field ) const;

   void getUnmappedField(const LevelData<FluxBox>& a_mapped_E_field,
                         LevelData<FluxBox>&       a_E_field) const;

   void extrapBoundaryGhostsForCC(FArrayBox&                              data,
                                  const Box&                              interiorbox,
                                  const Box&                              domain_box,
                                  const int                               dir,
                                  const Tuple<BlockBoundary, 2*SpaceDim>& block_boundaries) const;

   void extrapBoundaryGhostsForFC(FArrayBox&                              data,
                                  const Box&                              interiorbox,
                                  const Box&                              domain_box,
                                  const int                               dir,
                                  const int                               order,
                                  const Tuple<BlockBoundary, 2*SpaceDim>& a_block_boundaries) const;

   void extractNormalComponent( const LevelData<FluxBox>& in,
                                const int                 dir,
                                LevelData<FluxBox>&       out ) const;

   int numPotentialGhostsFilled() const {return m_num_potential_ghosts_interpolated;}

   bool secondOrder() const {return m_second_order;}

protected:

    /// Set all components of a LevelData<FArrayBox> to zero
   /**
    * Set all components of a LevelData<FArrayBox> to zero
    *
    * @param[in,out]  data   LevelData to be zeroed
    */
   void setZero( LevelData<FArrayBox>& data ) const;

   /// Set all components of a LevelData<FluxBox> to zero
   /**
    * Set all components of a LevelData<FluxBox> to zero
    *
    * @param[in,out]  data   LevelData to be zeroed
    */
   void setZero( LevelData<FluxBox>& data ) const;

   /// Compute the dot product of two LevelData<FArrayBox>s
   /**
    * Compute the dot product of two LevelData<FArrayBox>s
    *
    * @param[in]  data1   First vector
    * @param[in]  data2   Second vector
    * @param[out]         Dot product of data1 and data2
    */
   double dotProduct( const LevelData<FArrayBox>& data1,
                      const LevelData<FArrayBox>& data2 ) const;

   void zeroAverage( LevelData<FArrayBox>& data ) const;

   /// Compute the L2 norm of a LevelData<FArrayBox>
   /**
    * Compute the L2 norm of a LevelData<FArrayBox>
    *
    * @param[in]  data   Input vector
    * @param[out]        L2 Norm of data
    */
   double L2Norm( const LevelData<FArrayBox>& data ) const;

   double globalMax(const double data) const;

   void Update(LevelData<FArrayBox> &x, int k, double** h, double* s, int size_s, Vector<LevelData<FArrayBox> *> &v) const;

   void GeneratePlaneRotation (double &dx, double &dy, double &cs, double &sn) const;

   void ApplyPlaneRotation (double &dx, double &dy, double &cs, double &sn) const;

   MultiBlockLevelExchangeCenter* m_mblex_potential_Ptr; 

   const MagGeom& m_geometry;
   LevelData<FArrayBox> m_volume;
   LevelData<FArrayBox> m_bc_divergence;

   int m_num_potential_ghosts_filled;

   double m_tol;
   int m_max_iter;
   double m_amg_tol;
   int m_amg_max_iter;

   bool m_verbose;
   bool m_second_order;

   int m_num_potential_ghosts_interpolated;
};



#include "NamespaceFooter.H"

#endif
