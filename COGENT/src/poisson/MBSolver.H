#ifndef _MBSOLVER_H_
#define _MBSOLVER_H_

#include "MultiBlockLevelGeom.H"
#include "MultiBlockLevelExchangeCenter.H"
#include "EllipticOpBC.H"
#include "CoDim1Stencil.H"
#include "CoDim2Stencil.H"

#include "NamespaceHeader.H"



class MBSolver
{
   friend class MBSolverUtil;

public:

   /// Constructor with initialization.
   /**
    *  @param[in] geom the geometry
    *  @param[in] discretization_ordr discretization order
    *  @param[in] mblex_ptr level exchange object pointer (leave NULL
    *                       unless you know why you don't want to)
    */
   MBSolver( const MultiBlockLevelGeom&      geom,
             const int                       discretization_order,
             MultiBlockLevelExchangeCenter*  mblex_ptr = NULL);
      
   /// Destructor.
   /**
    */
   ~MBSolver();

   virtual void multiplyMatrix( const LevelData<FArrayBox>&  in,
                                LevelData<FArrayBox>&        out ) const = 0;

   virtual void solve( const LevelData<FArrayBox>&  rhs,
                       LevelData<FArrayBox>&        solution,
                       bool                         homogeneous_bcs ) = 0;

   virtual void dumpMatrix(const string& file_name) const = 0;

   void constructMatrix( LevelData<FluxBox>&   tensor_coefficient,
                         const EllipticOpBC&   bc );

   void constructMatrix( LevelData<FArrayBox>& alpha_coefficient,
                         LevelData<FluxBox>&   tensor_coefficient,
                         const EllipticOpBC&   bc );

   void constructMatrix( LevelData<FluxBox>&   tensor_coefficient,
                         LevelData<FArrayBox>& beta_coefficient,
                         const EllipticOpBC&   bc );

   void constructMatrix( LevelData<FArrayBox>& alpha_coefficient,
                         LevelData<FluxBox>&   tensor_coefficient,
                         LevelData<FArrayBox>& beta_coefficient,
                         const EllipticOpBC&   bc );

   const LevelData<FArrayBox>& getBCRhsContrib() const {return m_rhs_from_bc;}

protected:

   virtual void constructMatrixGeneral( LevelData<FArrayBox>& alpha_coefficient,
                                        LevelData<FluxBox>&   tensor_coefficient,
                                        LevelData<FArrayBox>& beta_coefficient,
                                        const EllipticOpBC&   bc ) = 0;

   void averageAtBlockBoundaries(LevelData<FluxBox>& data) const;

   IntVectSet getInterBlockCoupledCells( const int   block_number,
                                         const int   radius,
                                         const Box&  box ) const;

   IntVectSet getBoundaryCoupledCells( const int   block_number,
                                       const int   radius,
                                       const Box&  box ) const;

   void getUnstructuredCouplings(int                                 radius,
                                 LayoutData< BaseFab<IntVectSet> >&  unstructured_couplings ) const;

   IntVectSet getUnstructuredCouplingsForCell( const Vector< BaseFab<IntVectSet>* >&  coupling,
                                               const Box&                             stencil_box,
                                               const IntVect&                         iv ) const;

   void accumStencilMatrixEntries(const IntVect    index,
                                  const int        dir,
                                  const int        side,
                                  const int        dir2,
                                  const FluxBox&   mapped_coefs,
                                  const RealVect&  dx,
                                  const bool       fourthOrder,
                                  FArrayBox&       stencil_values) const;

   void modifyStencilForBCs( const Vector<CoDim1Stencil>&  codim1_stencil,
                             const Vector<CoDim2Stencil>&  codim2_stencil,
                             const IntVect&                iv,
                             FArrayBox&                    stencil_values,
                             FArrayBox&                    rhs_from_bc,
                             const bool                    rhs_from_bc_only,
                             const bool                    force_codim2_condense ) const;

   void constructBoundaryStencils( const bool                        fourth_order,
                                   const EllipticOpBC&               bc,
                                   Vector< Vector<CoDim1Stencil> >&  codim1_stencils,
                                   Vector< Vector<CoDim2Stencil> >&  codim2_stencils ) const;

   void accumPhysicalGhosts( const Vector< Vector<CoDim1Stencil> >&  codim1_stencils,
                             const Vector< Vector<CoDim2Stencil> >&  codim2_stencils,
                             const bool                              extrapolate_from_interior,
                             const bool                              include_bvs,
                             LevelData<FArrayBox>&                   data ) const;

   void computeFluxNormalFromStencil( const LevelData<FArrayBox>&                        in,
                                      LevelData<FluxBox>&                                tensor_coefficient,
                                      const Vector< Tuple<BlockBoundary, 2*SpaceDim> >&  block_boundaries,
                                      const Vector< Vector<CoDim1Stencil> >&             codim1_stencils,
                                      const Vector< Vector<CoDim2Stencil> >&             codim2_stencils,
                                      FArrayBox&                                         stencil_values,
                                      const bool                                         fourthOrder,
                                      const bool                                         extrapolate_from_interior,
                                      const bool                                         include_bvs,
                                      LevelData<FluxBox>&                                flux_normal ) const;

   void testMatrixConstruct( LevelData<FArrayBox>&                              alpha_coefficient, 
                             LevelData<FluxBox>&                                tensor_coefficient,
                             LevelData<FArrayBox>&                              beta_coefficient,
                             const Vector< Tuple<BlockBoundary, 2*SpaceDim> >&  block_boundaries,
                             const Vector< Vector<CoDim1Stencil> >&             codim1_stencils,
                             const Vector< Vector<CoDim2Stencil> >&             codim2_stencils,
                             FArrayBox&                                         stencil_values,
                             const bool                                         fourthOrder,
                             const LevelData<FArrayBox>&                        rhs_from_bc ) const;

   void plot( const string&               file_name,  
              const LevelData<FArrayBox>& data,
              const double&               time ) const;

   void setZero( LevelData<FArrayBox>& data ) const;

   const MultiBlockLevelGeom& m_geometry;
   MultiBlockCoordSys* m_coord_sys_ptr;

   MultiBlockLevelExchangeCenter* m_mblex_potential_Ptr; 

   LevelData<FArrayBox> m_rhs_from_bc;

   int m_discretization_order;

   bool m_mblex_defined_locally;
   bool m_flux_average[SpaceDim];
};



#include "NamespaceFooter.H"

#endif
