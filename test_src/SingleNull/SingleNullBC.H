#ifndef _SINGLENULLBC_
#define _SINGLENULLBC_

#include "ParmParse.H"
#include "RefCountedPtr.H"
#include "SNBC.H"
#include "GridFunctionLibrary.H"
#include "Vector.H"
#include "MagGeom.H"

#include "BoundaryBoxLayout.H.multidim"

#include "NamespaceHeader.H"

class SingleNullBC
   : public SNBC
{
public:

   SingleNullBC( const std::string& name,
                 ParmParse& pp,
                 const int& verbosity );


   ~SingleNullBC();

   /// Set the boundary conditions.
   /**
    * Sets the boundary condition on the species for the given geometry
    * at the given time.  Assumes distribution function is in computational
    * coordinates.
    *
    * @param[out] solution_comp Solution to set.
    * @param[in] velocity       Phase space velocity field.
    * @param[in] time           Time value at which to evaluate the function.
    */
   virtual void apply( LevelData<FArrayBox>& solution_comp,
                       const MagGeom& geometry,
                       const LevelData<FluxBox>& velocity,
                       const Real& time );
   
   virtual void printParameters() const;

private:

   // prohibit copying
   SingleNullBC( const SNBC& );
   SingleNullBC& operator=( const SNBC& );

   void fillInflowData( Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_bdry_data,
                        const LevelData<FArrayBox>& solution_comp,
                        const BoundaryBoxLayoutPtrVect& a_bdry_layout,
                        const MagGeom& a_geometry,
                        const Real& a_time );

   inline void parseParameters( ParmParse& a_pp );

   inline GridFunction& inflowFunc( const int& a_dir,
                                    const Side::LoHiSide& a_side,
                                    const int& a_block_type );

   inline GridFunction& radialInflowFunc( const Side::LoHiSide& a_side,
                                          const int& a_block_type );

   inline GridFunction& poloidalInflowFunc( const Side::LoHiSide& a_side );

   const std::string m_name;
   int m_verbosity;

   enum {INVALID=-1,
         RADIAL_CORE, RADIAL_SOL, RADIAL_PF,
         POLOIDAL_INNER_DIV, POLOIDAL_OUTER_DIV,
         NUM_INFLOW};

   Vector<RefCountedPtr<GridFunction> > m_inflow_function;
   Vector<std::string> m_bdry_name;
};

#include "NamespaceFooter.H"

#endif
