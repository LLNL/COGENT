#ifndef _POTENTIALBC_
#define _POTENTIALBC_

#include "MultiBlockLevelGeom.H"
#include "GridFunction.H"

#include "NamespaceHeader.H"

class PotentialBC
{
public:

   PotentialBC(const int num_boundaries);

   virtual ~PotentialBC();

   virtual void setBCType( const int block_number,
                           const int dir,
                           const int side,
                           const int type ) = 0;
   
   virtual int getBCType( const int block_number,
                          const int dir,
                          const int side ) const = 0;
   
   virtual void setBCValue( const int block_number,
                            const int dir,
                            const int side,
                            const double value) = 0;

   virtual double getBCValue( const int block_number,
                              const int dir,
                              const int side ) const = 0;
   
   virtual void setBCFunction( const int block_number,
                               const int dir,
                               const int side,
                               const RefCountedPtr<GridFunction>& function ) = 0;

   virtual RefCountedPtr<GridFunction> getBCFunction( const int block_number,
                                                      const int dir,
                                                      const int side ) const = 0;

   virtual void apply( const MultiBlockLevelGeom& a_geom,
                       const Box&                 a_coord_sys_box,
                       const double&              a_time,
                       const int                  a_dir,
                       const int                  a_side,
                       FArrayBox&                 a_phi ) const = 0;

   bool hasNeumannCondition() const;

   static const int DIRICHLET = 0;
   static const int NEUMANN = 1;
   static const int UNDEFINED = -1;

protected:

   Vector<RefCountedPtr<GridFunction> > m_bc_function;
   Vector<std::string> m_bdry_name;
   Vector<int> m_bc_type;
   Vector<double> m_bc_value;
};

#include "NamespaceFooter.H"

#endif
