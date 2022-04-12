#ifndef MAGCOORDSYS_H
#define MAGCOORDSYS_H

#include "MagBlockCoordSys.H"
#include "MultiBlockCoordSys.H"
#include "EllipticOpBC.H"

#include "NamespaceHeader.H"

class MagCoordSys : public MultiBlockCoordSys
{
public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor.
  */
   MagCoordSys()
      : m_sheared_geometry(false) {;}

  ///
  /** Destructor.
   */
   virtual ~MagCoordSys() {;}

  /*@}*/

   /// Return the mapping which corresponds to the block with the given index.
   virtual const MagBlockCoordSys* getCoordSys(int a_n) const {return (const MagBlockCoordSys*)m_coordSysVect[a_n];}

   /// Return the mapping which corresponds to the given box.
   /// \param a_box The Box for which the coordinate system is retrieved. This
   ///              Box must be contained completely within a block or this
   ///              method will fail.
   virtual const MagBlockCoordSys* getCoordSys(const Box& a_box) const
     { return (const MagBlockCoordSys*)m_coordSysVect[whichBlock(a_box)]; }

   void defineCoordSystemsAndBoundaries( Vector<MagBlockCoordSys *>& geom_vec);

   virtual void defineBoundaries() = 0;

   /// general block mapping conversion function
   /** Given a_xiSrc in mapped space of block a_nSrc,
       and block index a_nDst,
       returns coordinates in mapped space of block a_nDst.
   */
   virtual RealVect blockRemappingGeneral(int a_nDst,
                                          const RealVect& a_xiSrc,
                                          int a_nSrc) const;

   virtual IntVect getDecomposition(int block) const = 0;

   virtual void postProcessMetricData( LevelData<FArrayBox>& data ) const {;}

   virtual void postProcessMetricData( LevelData<FluxBox>& data ) const {;}
   
   virtual void postProcessMetricData( LevelData<EdgeDataBox>& data ) const {;}

   virtual EllipticOpBC* getDivergenceCleaningBC() {return NULL;}

   virtual void defineEllipticOpBC( EllipticOpBC& bc ) const {};

   virtual std::string type() const = 0;

   virtual bool containsPhysicalBoundary( int                    block_number,
                                          int                    dir,
                                          const Side::LoHiSide&  side ) const;

   virtual void setBoundaryBoxes();

   const Vector<Box>& getBoundaryBoxes( const int              block_number,
                                        const int              dir,
                                        const Side::LoHiSide&  side )
   {
      return m_boundary_boxes[block_number*2*SpaceDim + dir + side*SpaceDim];
   }

#if CFG_DIM ==3
   /// toroidal block mapping conversion function
   /** Given a_xiSrc in mapped space of block a_nSrc,
    and a_side, returns coordinate in mapped space of the
    destination block (determined from a_side argument)
    and the  second-order interpolation stencil.
    */

   virtual void toroidalBlockRemapping( IntVect&               ivDst,
                                        int&                   nDst,
                                        Vector<Real>&          interpStecil,
                                        Vector<int>&           interpStecilOffsets,
                                        const RealVect&        xiSrc,
                                        const int              nSrc,
                                        const Side::LoHiSide&  side ) const = 0;

#endif


protected:

   bool m_sheared_geometry;

   Vector<Vector<Box> > m_boundary_boxes;
};



#include "NamespaceFooter.H"

#endif
