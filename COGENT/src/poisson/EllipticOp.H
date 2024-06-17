#ifndef _ELLIPTICOP_H_
#define _ELLIPTICOP_H_

#include "MagGeom.H"
#include "MultiBlockLevelExchangeCenter.H"
#include "LinearSolver.H"
#include "CoDim1Stencil.H"
#include "CoDim2Stencil.H"
#include "SpaceUtils.H.multidim"


#include "NamespaceHeader.H"



class EllipticOp
   : public LinearOp< LevelData<FArrayBox> >
{
public:

   /// Constructor with initialization.
   /**
    *  @param[in] pp the input database.
    */
   EllipticOp(const ParmParse&  pp,
              const MagGeom&  geometry);
   EllipticOp( const ParmParse&  pp,
               const MagGeom&    geometry,
               const bool        second_order,
               const bool        low_pollution );
      
   /// Destructor.
   /**
    */
   virtual ~EllipticOp();

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
   void defineSolver( const ParmParse& a_pp );

   // Defines the linear solver from the argument list
   void defineSolver( const string&  a_method,
                      const double   a_tol,
                      const int      a_max_iter,
                      const bool     a_verbose );

   // Accessors to set the solver method, convergence params
   // and verbosity for subsequent solves.  These are normally
   // set upon object construction using the ParmParse version
   // above, but can be reset on the fly.
   void setMethod(const string& a_method) 
         {defineSolver(a_method, m_tol, m_max_iter, m_verbose);}
   void setTol(const double a_tol) {m_tol = a_tol;}
   void setMaxIter(const int a_max_iter) {m_max_iter = a_max_iter;}
   void setVerbose(const bool a_verbose) {m_verbose = a_verbose;}

   /// Solves the elliptic operator given a right-hand side and bcs
   /**
    * Solves the elliptic operator given a right-hand side and bcs
    *
    * @param[in] rhs   right-hand side
    * @param[out] solution the pointwise, cell-centered solution
    */
   void solveWithBCs( LevelData<FArrayBox>&        solution,
                      const LevelData<FArrayBox>&  rhs );

   /// Solves the elliptic operator given a right-hand side only
   /**
    * Solves the elliptic operator given a right-hand side only
    *
    * @param[in] rhs   right-hand side
    * @param[out] solution the pointwise, cell-centered solution
    */
   void solve( const LevelData<FArrayBox>& rhs,
               LevelData<FArrayBox>&       solution );

   // Computes the face-centered (i.e., pointwise) 3D field in the physical frame
   // utilizing boundary conditions/values
   void compute3DFieldWithBCs( const LevelData<FArrayBox>& phi,
                                     LevelData<FluxBox>&         field,
                                     const bool                  homogeneousBCs ) const;

   // Computes the cell-centered (i.e., pointwise) 3D field in the physical frame
   // utilizing boundary conditions/values
   void compute3DFieldWithBCs(const LevelData<FArrayBox>& phi,
                              LevelData<FArrayBox>&       field,
                              const bool                  homogeneousBCs ) const;


   // Computes the face-centered (i.e., pointwise) poloidal field in the physical frame
   // utilizing boundary conditions/values
   void computePoloidalFieldWithBCs( const LevelData<FArrayBox>& phi,
                                     LevelData<FluxBox>&         field,
                                     const bool                  homogeneousBCs ) const;

   // Computes the face-centered (i.e., pointwise) 3D field in the mapped frame
   // utilizing boundary conditions/values
   void computeMapped3DFieldWithBCs(const LevelData<FArrayBox>& phi,
                                    LevelData<FluxBox>&         field,
                                    const bool                  homogeneousBCs ) const;

   // Computes the face-centered (i.e., pointwise) 3D field in the mapped frame
   // utilizing boundary conditions/values
   void computeMapped3DFieldWithBCs(const LevelData<FArrayBox>& phi,
                                    LevelData<FArrayBox>&       field,
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

   // Computes the cell-centered (i.e., pointwise) 3D field in the mapped frame
   // using interior potentials only (i.e., ignoring boundary conditions/values)
   void computeMapped3DField( const LevelData<FArrayBox>& phi,
                                    LevelData<FArrayBox>&       field ) const;

   // Computes the face-centered (i.e., pointwise) poloidal field in the mapped frame
   // using interior potentials only (i.e., ignoring boundary conditions/values)
   void computeMapped3DField( const LevelData<FArrayBox>& phi,
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

   void computeFluxDivergence( const LevelData<FArrayBox>&  in,
                               LevelData<FArrayBox>&        out,
                               const bool                   homogeneous_bcs,
                               const bool                   extrap_to_ghosts = false );

   virtual void computeBcDivergence( LevelData<FArrayBox>& out );

   void updateBoundaryData(const LevelData<FluxBox>&  unmapped_coefficients,
                           EllipticOpBC&              bc,
                           const bool                 homogeneous_bcs = false);
   
   void addBcDivergence(LevelData<FArrayBox>& rhs) const;

   void subtractBcDivergence(LevelData<FArrayBox>& rhs) const;

   /// Computes the operator (normalized by the shell volume) in the lHS of Eq. (5) in the "Poloidal_variations_PB_model-II" note.
   /**
    * Computes <grad(psi) * D *NJinv * in> * S/(one-cell wide shell volume) (in corresponds to mapped_grad(Phi) in the notes)
    * LevelData<FArrayBox> out does not have ghost cells, and stores the data in the radial direction  as follows
    * out[i] = <grad(psi) * D *NJinv* in[i]> * S[i]/shell_volume[i]
    * lo_value = out[0]=<grad(psi) * D *NJinv* in[0]> * S[0]/shell_volume[0]
    * hi_value = <grad(psi) * D *NJinv* in[nrad+1]> * S[nrad]/shell_volume[nrad]
    * n_rad number of the core cells in the radial direction 
    * <...> denotes the shell-average, S is the surface area of a flux surface
    */

   void computeRadialFSAverage( const LevelData<FluxBox>& in,
                                double&                   lo_value,
                                double&                   hi_value,
                                LevelData<FArrayBox>&     out ) const;
   
   void modifyForNeumannAndNaturalBCs( const EllipticOpBC&  bc,
                                       LevelData<FluxBox>&  unmapped_coefs,
                                       LevelData<FluxBox>&  mapped_coefs ) const;

   void setNeumannNaturalFactor( const LevelData<FluxBox>&  unmapped_coefficients,
                                 EllipticOpBC&              bc ) const;

   void addNeumannAndNaturalBVContrib( const EllipticOpBC&    bc,
                                       LevelData<FArrayBox>&  data ) const;

   void fillInternalGhosts( LevelData<FArrayBox>& phi ) const;

   void extractNormalComponent( const LevelData<FluxBox>& in,
                                const int                 dir,
                                LevelData<FluxBox>&       out ) const;

   void interpToNodes( const LevelData<FArrayBox>& a_phi,
                       LevelData<FArrayBox>&       a_phi_node ) const;

   int numPotentialGhosts() const {return m_num_potential_ghosts;}

   bool secondOrder() const {return m_second_order;}
   bool lowPollution() const {return m_low_pollution;}

   // ----> Begin *this virtuals
   
   virtual void solvePreconditioner( const LevelData<FArrayBox>& in, LevelData<FArrayBox>& out ) = 0;

   virtual void multiplyCoefficients( LevelData<FluxBox>& data,
                                      const bool mapped_coeff ) const = 0;

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

   // Computes the cell-centered (i.e., pointwise) field in the mapped frame. It assumes
   // that ghost cells have already been filled by members in the public interface and
   // should not be called directly.
   void computeMapped3DFieldWithGhosts( const LevelData<FArrayBox>& phi,
                                              LevelData<FArrayBox>&       field,
                                              const int                   order ) const;

   // Computes the face-centered (i.e., pointwise) field in the mapped frame. It assumes
   // that ghost cells have already been filled by members in the public interface and
   // should not be called directly.
   void computeMapped3DFieldWithGhosts( const LevelData<FArrayBox>& phi,
                                              LevelData<FluxBox>&         field,
                                              const int                   order ) const;


   void constructBoundaryStencils( const bool                        fourth_order,
                                   const EllipticOpBC&               bc,
                                   Vector< Vector<CoDim1Stencil> >&  codim1_stencils,
                                   Vector< Vector<CoDim2Stencil> >&  codim2_stencils ) const;

   void accumPhysicalGhosts( const Vector< Vector<CoDim1Stencil> >&  codim1_stencils,
                             const Vector< Vector<CoDim2Stencil> >&  codim2_stencils,
                             const bool                              extrapolate_from_interior,
                             const bool                              include_bvs,
                             LevelData<FArrayBox>&                   data ) const;

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

   bool m_second_order;
   bool m_low_pollution;
   int m_num_potential_ghosts;
   
   LinearSolver< LevelData<FArrayBox> >* m_Chombo_solver;

   Vector< Vector<CoDim1Stencil> > m_codim1_stencils;
   Vector< Vector<CoDim2Stencil> > m_codim2_stencils;
   RefCountedPtr<EllipticOpBC> m_bc;

private:

   void init( const ParmParse& pp );
};



#include "NamespaceFooter.H"

#endif
