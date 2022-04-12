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

   virtual void initializeMatrix() = 0;

   virtual void finalizeMatrix() = 0;

   virtual void setMethodParams( const ParmParse&  pp ) = 0;

   void setConvergenceParams( const ParmParse&  pp );

   void setConvergenceParams( const double  method_tol,
                              const int     method_max_iter,
                              const bool    method_verbose,
                              const double  precond_tol,
                              const int     precond_max_iter,
                              const bool    precond_verbose );

   void constructMatrix( LevelData<FluxBox>&   tensor_coefficient,
                         const EllipticOpBC&   bc );

   void constructMatrixBlock( const int             block_row,
                              const int             block_column,
                              LevelData<FluxBox>&   tensor_coefficient,
                              const EllipticOpBC&   bc );

   void constructMatrix( LevelData<FArrayBox>& alpha_coefficient,
                         LevelData<FluxBox>&   tensor_coefficient,
                         const EllipticOpBC&   bc );

   void constructMatrixBlock( const int             block_row,
                              const int             block_column,
                              LevelData<FArrayBox>& alpha_coefficient,
                              LevelData<FluxBox>&   tensor_coefficient,
                              const EllipticOpBC&   bc );

   void constructMatrix( LevelData<FluxBox>&   tensor_coefficient,
                         LevelData<FArrayBox>& beta_coefficient,
                         const EllipticOpBC&   bc );

   void constructMatrixBlock( const int             block_row,
                              const int             block_column,
                              LevelData<FluxBox>&   tensor_coefficient,
                              LevelData<FArrayBox>& beta_coefficient,
                              const EllipticOpBC&   bc );

   void constructMatrix( LevelData<FArrayBox>& alpha_coefficient,
                         LevelData<FluxBox>&   tensor_coefficient,
                         LevelData<FArrayBox>& beta_coefficient,
                         const EllipticOpBC&   bc );

   void constructMatrixBlock( const int             block_row,
                              const int             block_column,
                              LevelData<FArrayBox>& alpha_coefficient,
                              LevelData<FluxBox>&   tensor_coefficient,
                              LevelData<FArrayBox>& beta_coefficient,
                              const EllipticOpBC&   bc );

   string getMethod() const {return m_method;}
   double getTol() const {return m_method_tol;}
   int getMaxIter() const {return m_method_max_iter;}
   bool getVerbose() const {return m_method_verbose;}

   string getPrecondMethod() const {return m_precond_method;}
   double getPrecondTol() const {return m_precond_tol;}
   int getPrecondMaxIter() const {return m_precond_max_iter;}
   bool getPrecondVerbose() const {return m_precond_verbose;}

protected:

   virtual void constructMatrixGeneral( LevelData<FArrayBox>&  alpha_coefficient,
                                        LevelData<FluxBox>&    tensor_coefficient,
                                        LevelData<FArrayBox>&  beta_coefficient,
                                        const EllipticOpBC&    bc ) = 0;

   virtual void constructMatrixBlockGeneral( const int              block_row,
                                             const int              block_column,
                                             LevelData<FArrayBox>&  alpha_coefficient,
                                             LevelData<FluxBox>&    tensor_coefficient,
                                             LevelData<FArrayBox>&  beta_coefficient,
                                             const EllipticOpBC&    bc ) = 0;

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
                             const int                     rhs_from_bc_comp,
                             const bool                    rhs_from_bc_only,
                             const bool                    force_codim2_condense,
                             const bool                    update_rhs) const;

   void constructBoundaryStencils( const bool                        fourth_order,
                                   const EllipticOpBC&               bc,
                                   Vector< Vector<CoDim1Stencil> >&  codim1_stencils,
                                   Vector< Vector<CoDim2Stencil> >&  codim2_stencils ) const;

   void accumPhysicalGhosts( const Vector< Vector<CoDim1Stencil> >&  codim1_stencils,
                             const Vector< Vector<CoDim2Stencil> >&  codim2_stencils,
                             const bool                              extrapolate_from_interior,
                             const bool                              include_bvs,
                             LevelData<FArrayBox>&                   data ) const;

   double getBV( const EllipticOpBC&      bc,
                 const MagBlockCoordSys&  coord_sys,
                 const int                block_number,
                 const int                dir,
                 const Side::LoHiSide&    side,
                 const IntVect            iv_face ) const;

   void getNeumannNaturalFactor( const EllipticOpBC&  bc,
                                 LevelData<FluxBox>&  bc_factor ) const;

   void computeFluxNormalFromStencil( const LevelData<FArrayBox>&                        in,
                                      LevelData<FluxBox>&                                tensor_coefficient,
                                      const Vector< Tuple<BlockBoundary, 2*SpaceDim> >&  block_boundaries,
                                      const Vector< Vector<CoDim1Stencil> >&             codim1_stencils,
                                      const Vector< Vector<CoDim2Stencil> >&             codim2_stencils,
                                      const EllipticOpBC&                                bc,
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
                             const EllipticOpBC&                                bc,
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

   int m_discretization_order;

   string m_method;
   double m_method_tol;
   int m_method_max_iter;
   bool m_method_verbose;

   string m_precond_method;
   double m_precond_tol;
   int m_precond_max_iter;
   bool m_precond_verbose;

   bool m_method_params_set;
   bool m_convergence_params_set;

   bool m_mblex_defined_locally;
   bool m_flux_average[SpaceDim];
};



#include "NamespaceFooter.H"

#endif
