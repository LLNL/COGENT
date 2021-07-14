#ifndef _ELLIPTICOPBC_
#define _ELLIPTICOPBC_

#include "MultiBlockLevelGeom.H"
#include "GridFunction.H"

#include "NamespaceHeader.H"

class EllipticOpBC
{
public:

   EllipticOpBC(const int num_boundaries,
                const int num_blocks);

   virtual ~EllipticOpBC();

   virtual void setBCType( const int block_number,
                           const int dir,
                           const int side,
                           const int type ) = 0;
   
   virtual int getBCType( const int block_number,
                          const int dir,
                          const int side ) const = 0;
 
   virtual std::string getBCSubType(const int block_number,
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

   virtual void apply( const MultiBlockLevelGeom&  geom,
                       const Box&                  coord_sys_box,
                       const double&               time,
                       const int                   dir,
                       const int                   side,
                       FArrayBox&                  phi ) const = 0;

   virtual RefCountedPtr<EllipticOpBC> clone( const bool extrapolated = false ) const = 0;

   bool hasNeumannCondition() const;
   
   bool hasCoupledBoundary() const;

   enum EllipticBCType {UNDEFINED = -1, DIRICHLET, NEUMANN, EXTRAPOLATED};

protected:
   
   RefCountedPtr<GridFunction> getBlockBCData( const int  block_number,
                                               const int  dir,
                                               const int  side ) const
   {
      return m_bc_block_data[block_number*2*SpaceDim + dir + side*SpaceDim];
   }

   Vector<RefCountedPtr<GridFunction> > m_bc_function;
   Vector<RefCountedPtr<GridFunction> > m_bc_block_data;
   Vector<std::string> m_bdry_name;
   Vector<int> m_bc_type;
   Vector<std::string> m_bc_subtype;
   Vector<double> m_bc_value;

   int m_num_boundaries;
   int m_num_blocks;
};

#include "NamespaceFooter.H"

#endif
