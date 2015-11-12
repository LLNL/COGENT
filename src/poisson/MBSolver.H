#ifndef _MBSOLVER_H_
#define _MBSOLVER_H_

#include "MagGeom.H"
#include "MultiBlockLevelExchangeCenter.H"
#include "SparseCoupling.H"

#include "NamespaceHeader.H"

class MBSolver
{
public:

   /// Constructor with initialization.
   /**
    *  @param[in] pp the input database.
    *  @param[in] geom the geometry
    */
   MBSolver( const MagGeom& geom,
             const int      discretization_order );
      
   /// Destructor.
   /**
    */
   ~MBSolver() {;}

   virtual void multiplyMatrix( const LevelData<FArrayBox>& in,
                                LevelData<FArrayBox>&       out ) = 0;

   virtual void updateMatrixDiagonals( const LevelData<FArrayBox>& diagonal_increment ) = 0;

   const LevelData<FArrayBox>& getARadial() const {return m_A_radial;}

protected:

   IntVectSet getInterBlockCoupledCells( const int  block_number,
                                         const int  radius,
                                         const Box& box ) const;

   void getUnstructuredCouplings(int                                radius,
                                 LayoutData< BaseFab<IntVectSet> >& unstructured_couplings ) const;

   IntVectSet getUnstructuredCouplingsForCell( const Vector< BaseFab<IntVectSet>* >& coupling,
                                               const Box&                            stencil_box,
                                               const IntVect&                        iv ) const;

   void accumStencilMatrixEntries(const IntVect                           index,
                                  const int                               dir,
                                  const int                               side,
                                  const int                               s,
                                  const Tuple<BlockBoundary, 2*SpaceDim>& block_boundaries,
                                  const int*                              block_bc_type,
                                  const FluxBox&                          mapped_coefs,
                                  const RealVect&                         dx,
                                  const Box&                              domainBox,
                                  const bool                              fourthOrder,
                                  const bool                              dropOrder,
                                  FArrayBox&                              stencil_values) const;

   void setZero( LevelData<FArrayBox>& data ) const;

   const MagGeom& m_geometry;

   MultiBlockLevelExchangeCenter* m_mblex_potential_Ptr; 

   LevelData<FArrayBox> m_volume;
   LevelData<FArrayBox> m_A_diagonal;
   LevelData<FArrayBox> m_A_radial;
   LevelData<FArrayBox> m_diagonal_increment;

   int m_num_potential_ghosts_filled;
   int m_discretization_order;
   bool m_dropOrder;
};



#include "NamespaceFooter.H"

#endif
