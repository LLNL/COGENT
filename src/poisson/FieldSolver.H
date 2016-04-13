#ifndef _FIELDSOLVER_H_
#define _FIELDSOLVER_H_

#include "MagGeom.H"
#include "MultiBlockLevelExchangeCenter.H"
#include "LinearSolver.H"
#include "CoDim1Stencil.H"
#include "CoDim2Stencil.H"

#include "NamespaceHeader.H"



class FieldSolver
   : public LinearOp< LevelData<FArrayBox> >
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

   /// Defines the linear solver
   /**
    * Defines the linear solver from a ParmParse object, or,
    * if no definition is provided, using defaults.
    * It is assumed that the linear solver is preconditioned,
    * and that the preconditioner may itself be a preconditioned
    * iterative method.  The linear solver, preconditioner
    * (if provided) and preconditioner preconditioner (if
    * provided) are each specified by a "method" string (e.g.,
    * "GMRES" or "AMG").  This function also sets convergence
    * tolerances, maximum number of iterations and a verbosity
    * flag for each of the three solvers.  The convergence
    * tolerance and maximum number of iterations can be
    * subsequently reset using the setConvergenceParams member.
    *
    * @param[in] a_pp  ParmParse object reference
    */
   void defineLinearSolver( const ParmParse& a_pp );

   /// Computes the potential given the species.
   /**
    * Computes the pointwise, cell-centered potential field given the
    * charge density.
    *
    * @param[in] charge_density charge density
    * @param[out] phi the pointwise, cell-centered potential.
    */
   void computePotential( LevelData<FArrayBox>&       phi,
                          const LevelData<FArrayBox>& charge_density );

   void setConvergenceParams( const double tol,
                              const int    max_iter,
                              const bool   verbose,
                              const double precond_tol,
                              const int    precond_max_iter,
                              const bool   precond_verbose,
                              const double precond_precond_tol,
                              const int    precond_precond_max_iter,
                              const bool   precond_precond_verbose )
   {
      m_tol = tol; m_max_iter = max_iter; m_verbose = verbose;
      m_precond_tol = precond_tol; m_precond_max_iter = precond_max_iter;
      m_precond_verbose = precond_verbose; m_precond_precond_tol = precond_precond_tol;
      m_precond_precond_max_iter = precond_precond_max_iter;
      m_precond_precond_verbose = precond_precond_verbose;
   }

   void solve( const LevelData<FArrayBox>& rhs,
               LevelData<FArrayBox>&       solution );

   // Computes the face-centered (i.e., pointwise) poloidal field in the physical frame
   // utilizing boundary conditions/values
   void computePoloidalFieldWithBCs( const LevelData<FArrayBox>& phi,
                                     LevelData<FluxBox>&         field,
                                     const bool                  homogeneousBCs ) const;

   // Computes the face-centered (i.e., pointwise) poloidal field in the mapped frame
   // utilizing boundary conditions/values
   void computeMappedPoloidalFieldWithBCs( const LevelData<FArrayBox>& phi,
                                           LevelData<FluxBox>&         field,
                                           const bool                  homogeneousBCs ) const;

   // Computes the cell-centered (i.e., pointwise) poloidal field in the physical frame
   // using interior potentials only (i.e., ignoring boundary conditions/values)
   void computePoloidalField( const LevelData<FArrayBox>& phi,
                              LevelData<FArrayBox>&       field ) const;

   // Computes the face-centered (i.e., pointwise) poloidal field in the physical frame
   // using interior potentials only (i.e., ignoring boundary conditions/values)
   void computePoloidalField( const LevelData<FArrayBox>& phi,
                              LevelData<FluxBox>&         field ) const;

   // Computes the cell-centered (i.e., pointwise) poloidal field in the mapped frame
   // using interior potentials only (i.e., ignoring boundary conditions/values)
   void computeMappedPoloidalField( const LevelData<FArrayBox>& phi,
                                    LevelData<FArrayBox>&       field ) const;

   // Computes the face-centered (i.e., pointwise) poloidal field in the mapped frame
   // using interior potentials only (i.e., ignoring boundary conditions/values)
   void computeMappedPoloidalField( const LevelData<FArrayBox>& phi,
                                    LevelData<FluxBox>&         field ) const;

   // Computes the cell-centered (i.e., pointwise) field in the physical frame using
   // interior potentials only (i.e., ignoring boundary conditions/values).
   void computeField( const LevelData<FArrayBox>& phi,
                      LevelData<FArrayBox>&       field ) const;

   // Computes the face-centered (i.e., pointwise) field in the physical frame using
   // interior potentials only (i.e., ignoring boundary conditions/values).
   void computeField( const LevelData<FArrayBox>& phi,
                      LevelData<FluxBox>&         field ) const;

   // Computes the cell-centered (i.e., pointwise) field in the mapped frame using
   // interior potentials only (i.e., ignoring boundary conditions/values).
   void computeMappedField( const LevelData<FArrayBox>& phi,
                            LevelData<FArrayBox>&       field ) const;

   // Computes the face-centered (i.e., pointwise) field in the mapped frame using
   // interior potentials only (i.e., ignoring boundary conditions/values).
   void computeMappedField( const LevelData<FArrayBox>& phi,
                            LevelData<FluxBox>&         field ) const;

   void computeFluxDivergence( const LevelData<FArrayBox>& in,
                               LevelData<FArrayBox>&       out,
                               const bool                  homogeneous_bcs = false );

   void computeBcDivergence( const PotentialBC&    bc,
                             LevelData<FArrayBox>& out ) const;

   void setBc( const PotentialBC& bc );

   void addBcDivergence(LevelData<FArrayBox>& rhs) const;

   void subtractBcDivergence(LevelData<FArrayBox>& rhs) const;

   /// Computes the operator (normalized by the shell volume) in the lHS of Eq. (5) in the "Poloidal_variations_PB_model-II" note.
   /**
    * Computes <grad(psi) * D *NJinv * in>/(one-cell wide shell volume) (in corresponds to mapped_grad(Phi) in the notes)
    * LevelData<FArrayBox> out does not have ghost cells, and stores the data in the radial direction  as follows
    * out[i] = <grad(psi) * D *NJinv* in[i]>/shell_volume[i]
    * lo_value = out[0]=<grad(psi) * D *NJinv* in[0]>/shell_volume[0]
    * hi_value = <grad(psi) * D *NJinv* in[nrad+1]>/shell_volume[nrad]
    * n_rad number of the core cells in the radial direction 
    * <...> denotes the shell-average
    */

   void computeRadialFSAverage( const LevelData<FluxBox>& in,
                                double&                   lo_value,
                                double&                   hi_value,
                                LevelData<FArrayBox>&     out ) const;
    
   void fillInternalGhosts( LevelData<FArrayBox>& phi ) const;

   // Convert a 3-dimensional gradient in mapped coordinates to the physical frame
   void unmapGradient( const LevelData<FArrayBox>& mapped_gradient,
                       LevelData<FArrayBox>&       gradient ) const;

   // Convert a 3-dimensional gradient in mapped coordinates to the physical frame
   void unmapGradient( const LevelData<FluxBox>& mapped_gradient,
                       LevelData<FluxBox>&       gradient ) const;

   // Convert a poloidal gradient in mapped coordinates to the physical the frame
   void unmapPoloidalGradient( const LevelData<FArrayBox>& mapped_gradient,
                               LevelData<FArrayBox>&       gradient ) const;

   // Convert a poloidal gradient in mapped coordinates to the physical the frame
   void unmapPoloidalGradient( const LevelData<FluxBox>& mapped_gradient,
                               LevelData<FluxBox>&       gradient ) const;

   void extractNormalComponent( const LevelData<FluxBox>& in,
                                const int                 dir,
                                LevelData<FluxBox>&       out ) const;

   void parseMethodAndParams( const ParmParse& pp,
                              string&          method,
                              double&          tol,
                              int&             max_iter,
                              bool&            verbose ) const;

   int numPotentialGhosts() const {return m_num_potential_ghosts;}

   bool secondOrder() const {return m_second_order;}

   // ----> Begin *this virtuals
   
   virtual void setPreconditionerConvergenceParams( const double tol,
                                                    const int    max_iter,
                                                    const double precond_tol,
                                                    const int    precond_max_iter ) = 0;

   virtual void solvePreconditioner( const LevelData<FArrayBox>& in, LevelData<FArrayBox>& out ) = 0;

   virtual void multiplyUnmappedCoefficients( LevelData<FluxBox>& data ) const = 0;

   // <---- End *this virtuals

   // ----> Begin LinOp virtuals

   virtual void preCond( LevelData<FArrayBox>&       cor,
                         const LevelData<FArrayBox>& residual );

   virtual void applyOp( LevelData<FArrayBox>&       out,
                         const LevelData<FArrayBox>& in,
                         bool                        homogeneous = false);

   virtual void residual( LevelData<FArrayBox>&       lhs,
                          const LevelData<FArrayBox>& phi,
                          const LevelData<FArrayBox>& rhs,
                          bool                        homogeneous = false );

   virtual void create( LevelData<FArrayBox>&       lhs,
                        const LevelData<FArrayBox>& rhs);

   virtual void clear(LevelData<FArrayBox>& lhs) {};

   virtual void assign( LevelData<FArrayBox>& lhs, const LevelData<FArrayBox>& rhs);

   virtual void assignLocal( LevelData<FArrayBox>& a_lhs, const LevelData<FArrayBox>& a_rhs)
   {this->assign(a_lhs, a_rhs);}

   virtual Real dotProduct(const LevelData<FArrayBox>& a_1, const LevelData<FArrayBox>& a_2);

   virtual void incr  ( LevelData<FArrayBox>&       lhs,
                        const LevelData<FArrayBox>& x,
                        Real                        scale);

   virtual void axby( LevelData<FArrayBox>&       lhs,
                      const LevelData<FArrayBox>& x,
                      const LevelData<FArrayBox>& y,
                      Real a, Real b);

   virtual void scale( LevelData<FArrayBox>& lhs, const Real& scale);

   virtual Real norm(const LevelData<FArrayBox>& rhs,
                     int ord);

   virtual Real dx() const {MayDay::Warning(" calling dx on base class\n"); return 0;}

   virtual void setToZero(LevelData<FArrayBox>& lhs);

   virtual void write(const LevelData<FArrayBox>* a, const char* filename)
   {MayDay::Warning("LinearOp::write not implemented");}

   // <---- End LinOp virtuals

protected:

   // Computes the cell-centered (i.e., pointwise) field in the mapped frame. It assumes
   // that ghost cells have already been filled by members in the public interface and
   // should not be called directly.
   void computeMappedPoloidalFieldWithGhosts( const LevelData<FArrayBox>& phi,
                                              LevelData<FArrayBox>&       field,
                                              const int                   order ) const;

   // Computes the face-centered (i.e., pointwise) field in the mapped frame. It assumes
   // that ghost cells have already been filled by members in the public interface and
   // should not be called directly.
   void computeMappedPoloidalFieldWithGhosts( const LevelData<FArrayBox>& phi,
                                              LevelData<FluxBox>&         field,
                                              const int                   order ) const;

   void constructBoundaryStencils( const bool                        fourth_order,
                                   const PotentialBC&                bc,
                                   Vector< Vector<CoDim1Stencil> >&  codim1_stencils,
                                   Vector< Vector<CoDim2Stencil> >&  codim2_stencils ) const;

   void accumPhysicalGhosts( const Vector< Vector<CoDim1Stencil> >&  codim1_stencils,
                             const Vector< Vector<CoDim2Stencil> >&  codim2_stencils,
                             const bool                              extrapolate_from_interior,
                             const bool                              include_bvs,
                             LevelData<FArrayBox>&                   data ) const;

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
                                  const Tuple<BlockBoundary, 2*SpaceDim>& block_boundaries) const;

   void injectPoloidalVector( const LevelData<FArrayBox>& poloidal_vector,
                              LevelData<FArrayBox>&       vector ) const;

   void injectPoloidalVector( const LevelData<FluxBox>& poloidal_vector,
                              LevelData<FluxBox>&       vector ) const;

   void projectPoloidalVector( const LevelData<FArrayBox>& a_vector,
                               LevelData<FArrayBox>&       a_poloidal_vector ) const;

   void projectPoloidalVector( const LevelData<FluxBox>& vector,
                               LevelData<FluxBox>&       poloidal_vector ) const;

   double L2Norm( const LevelData<FArrayBox>& data );

   double globalMax(const double data) const;

   MultiBlockLevelExchangeCenter* m_mblex_potential_Ptr; 

   const MagGeom& m_geometry;
   LevelData<FArrayBox> m_volume;
   LevelData<FArrayBox> m_volume_reciprocal;
   LevelData<FArrayBox> m_bc_divergence;

   string m_method;
   double m_tol;
   int m_max_iter;
   bool m_verbose;

   string m_precond_method;
   double m_precond_tol;
   int m_precond_max_iter;
   bool m_precond_verbose;

   string m_precond_precond_method;
   double m_precond_precond_tol;
   int m_precond_precond_max_iter;
   bool m_precond_precond_verbose;

   bool m_second_order;

   int m_num_potential_ghosts;

   LinearSolver< LevelData<FArrayBox> >* m_Chombo_solver;

   Vector< Vector<CoDim1Stencil> > m_codim1_stencils;
   Vector< Vector<CoDim2Stencil> > m_codim2_stencils;
};



#include "NamespaceFooter.H"

#endif
