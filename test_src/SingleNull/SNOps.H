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

#include "NamespaceHeader.H"

class SNOps
{
public:

   /// Constructor with initialization.
   /**
    */
   SNOps( ParmParse& pp,
          const MagGeom& geometry,
          const SNSystemBC& boundary_conditions,
          const SNSystemIC& initial_conditions,
          const int verbosity );

   inline
   void applyInitialConditions( LevelData<FArrayBox>& a_solution,
                                const Real a_time )
   {
      m_initial_conditions.initializeSystem( a_solution, a_time );
   }
   
   /// Destructor.
   /**
    */
   virtual ~SNOps();

   /// Compute a stable time step.
   /**
    * Computes and returns a stable time step estimate.
    */
   Real stableDt( const LevelData<FArrayBox>& soln, const int step_number );

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
   void evalRHS( LevelData<FArrayBox>&       rhs,
                 const LevelData<FArrayBox>& soln,
                 const int                   step_number,
                 const Real                  time,
                 const int                   stage );
   
   void writeCheckpointFile( HDF5Handle& handle ) const;
   
      /// Read checkpoint file.
      /**
       * Read checkpoint data from an output HDF5 file and reinitialize.
       *
       * @param[in] handle a reference to an HDF5 plot file.
       */
   void readCheckpointFile( HDF5Handle& handle, const int& cur_step  );
   
   void divideJ( const LevelData<FArrayBox>& soln_mapped,
                 LevelData<FArrayBox>&       soln_physical,
                 bool                        a_restrict_to_valid = false);

private:
   
   /// Parse parameters.
   /**
    * Private method to obtain control parameters from "gksystem" section
    * of ParmParse database.
    *
    * \param[in] ppgksys gksystem input database.
    */
   void parseParameters( ParmParse& ppgksys );

   Advect*                        m_operator;
   
   const SNSystemIC&              m_initial_conditions;
   const SNSystemBC&              m_boundary_conditions;
   const MagGeom&                 m_geometry;
   
   IntVect                        m_ghost_vect;

   int                            m_verbosity;
};

#include "NamespaceFooter.H"

#endif
