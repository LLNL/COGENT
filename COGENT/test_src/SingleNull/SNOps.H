#ifndef  _SNOPS_H_
#define  _SNOPS_H_

#include "Advect.H"

#include "SNSystemIC.H"
#include "SNSystemBC.H"

#include "MagGeom.H"
#include "ParmParse.H"
#include "Vector.H"
#include "CH_HDF5.H"
#include "RefCountedPtr.H"
#include "ProblemDomain.H"
#include "SNVector.H"

#include "NamespaceHeader.H"

class SNOps
{
public:

   /// Constructor with initialization.
   /**
    */
   SNOps() :
      m_operator(NULL),
      m_initial_conditions(NULL),
      m_boundary_conditions(NULL),
      m_dt(-1.0),
      m_time_scale(0.0),
      m_is_defined(false),
      m_ghost_vect(4*IntVect::Unit),
      m_verbosity(0) {;}   

   /// Destructor.
   /**
    */
   virtual ~SNOps();

   void define( const SNState& a_state, const Real a_dt );
   
   bool isDefined() const { return m_is_defined; }
   
   void resetDt( const Real a_dt )
   {
      CH_assert( a_dt>0.0 );
      CH_assert( isDefined() );
      m_dt = a_dt;
   }
   
   inline
   void initializeState( SNState& a_state, const Real a_time )
   {
      m_initial_conditions->initializeSystem( a_state, a_time );
   }

   /// Returns the time scale.
   /**
    * Returns the time scale.
    */
   Real dtScale( const SNState& state, const int step_number )
   {
      return m_time_scale;
   }
   
   /// Compute a stable time step.
   /**
    * Computes and returns a stable time step estimate.
    */
   Real stableDt( const SNState& state, const int step_number );

   void preTimeStep  (const int,const Real,const SNState&,const SNState&);
   void postTimeStep (const int,const Real,const SNState&);
   void postTimeStage(const int,const Real,const SNState&, const int);
   void postTimeStage(const int,const Real,const SNVector&, const int);

   /// Compute a stable time step for explicit time integration
   /**
    * Computes and returns a stable time step estimate.
    */
   Real stableDtExpl( const SNState& a_state, const int step_number );

   /// Compute a stable time step for IMEX time integration
   /**
    * Computes and returns a stable time step estimate.
    */
   Real stableDtImEx( const SNState& a_state, const int step_number );

   void setVelocities( FluidSpeciesPtrVect& species_vec,
                       const bool           mapped );

   /// Evaluates the RHS of the ODE.
   /**
    * Concrete implementation of pure virtual base class member function.
    * Given the current time and solution, evaluates the RHS of the ODE.
    *
    * @param[out] rhs data holder for rhs (without ghost cells).
    * @param[in] soln current solution with which to evaluate RHS.
    * @param[in] time current time.
    * @param[in] stage  stage number
    */
   void explicitOp( SNRHSData& a_rhs,
                    const Real a_time,
                    const SNState& a_state,
                    const int a_stage = 1 );

   void explicitOp( SNRHSData& a_rhs,
                    const Real a_time,
                    const SNRHSData& a_state,
                    const int a_stage = 1 );

   void explicitOp( SNVector& a_rhs,
                    const Real a_time,
                    const SNVector& a_state,
                    const int a_stage = 1 );

   void divideJ( const FluidSpeciesPtrVect& soln_mapped,
                 FluidSpeciesPtrVect&       soln_physical );

   void divideJ( const SNState& soln_mapped,
                 SNState&       soln_physical );

   
private:
   
   /// Parse parameters.
   /**
    * Private method to obtain control parameters from "gksystem" section
    * of ParmParse database.
    *
    * \param[in] ppgksys gksystem input database.
    */
   void parseParameters( ParmParse& ppgksys );

   void createTemporarySpeciesVector( FluidSpeciesPtrVect& out,
                                      const FluidSpeciesPtrVect& in );

   void createTemporaryState( SNState& out, const SNState& in );

   void fillGhostCells( FluidSpeciesPtrVect& field_phys, const Real& time );

   void applyOperator( FluidSpeciesPtrVect&       rhs,
                       const FluidSpeciesPtrVect& soln,
                       const Real&                time);

   Advect*           m_operator;
   
   const SNSystemIC* m_initial_conditions;
   const SNSystemBC* m_boundary_conditions;
   
   Real m_dt;
   Real m_time_scale;
   bool m_is_defined;

   SNState                        m_Y;
   SNRHSData                      m_rhs;

   IntVect           m_ghost_vect;

   int               m_verbosity;
};

#include "NamespaceFooter.H"

#endif
