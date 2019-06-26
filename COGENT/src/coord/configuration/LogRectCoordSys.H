#ifndef LOGRECTCOORDSYS_H
#define LOGRECTCOORDSYS_H

#include "MagBlockCoordSys.H"
#include "MagCoordSys.H"
#include "LogRectEllipticOpBC.H"

#include "NamespaceHeader.H"

#define BLOCK_SEPARATION 8

class LogRectCoordSys : public MagCoordSys
{
public:

   /**
      \name Constructors, destructor and defines
   */
   /*@{*/
   LogRectCoordSys(ParmParse&               a_pp,
                   const string&            a_mag_geom_type,
                   const std::vector<int>&  a_numcells,
                   const std::vector<bool>& a_is_periodic,
                   const std::vector<int>&  a_decomposition );

   ///
   /** Destructor.
    */
   virtual ~LogRectCoordSys();

   /*@}*/

   /// block mapping conversion function
   /** given xi_k in block k, returns which block (k_valid)
       contains that location in its valid region and its
       location in that block's computational coordinate space.
   */
   virtual void blockRemapping(RealVect&       xi_valid,
                               int&            n_valid,
                               const RealVect& xiSrc,
                               int             nSrc) const;

   virtual void defineBoundaries();

   virtual IntVect getDecomposition(int block) const {return m_decomposition;}

   virtual const EllipticOpBC& getDivergenceCleaningBC() const {return m_divergence_cleaning_bc;}

   bool isPeriodic( const int dir ) const;
   
   // Object ParmParse name.
   static const std::string pp_name;
   
#if CFG_DIM ==3
   /// toroidal block mapping conversion function
   /** Given a_xiSrc in mapped space of block a_nSrc,
    and a_side, returns coordinate in mapped space of the
    destination block (determined from a_side argument)
    and the  second-order interpolation stencil.
    */
   
   virtual void toroidalBlockRemapping(IntVect& a_ivDst,
                                       Vector<Real>& a_interpStecil,
                                       const RealVect& a_xiSrc,
                                       const int a_nSrc,
                                       const Side::LoHiSide& a_side) const;
#endif

   virtual std::string type() const { return "LogicallyRectangular"; }

protected:
   
   virtual void postProcessMetricData( LevelData<FArrayBox>& data ) const;
   
   virtual void postProcessMetricData( LevelData<FluxBox>& data ) const;
   
   void spreadRadially( LevelData<FArrayBox>& data ) const;
   
   void spreadRadially( LevelData<FluxBox>& data ) const;

   
private:

   bool m_spread_radially;
   IntVect m_decomposition;
   string m_mag_geom_type;
   LogRectEllipticOpBC m_divergence_cleaning_bc;

   void operator = (const LogRectCoordSys&){;}
   int m_mb_dir;
   int m_num_blocks;

   Vector<MagBlockCoordSys *> m_coord_vec;
   
#if CFG_DIM ==3 
   //applies toroidal periodicity
   void applyPeriodicity(RealVect& a_x) const;

   //enforces cut conditions in the poloidal direction
   //for the case of the toroidal geometry
   void enforcePoloidalCut(const RealVect& xiSrc,
                           double& xiNW_dst,
                           double& xiNE_dst,
                           double& xiSW_dst,
                           double& xiSE_dst,
			   double& xi0_dst) const;

   void getInterpolationCoefficients(Vector<Real>&    coeff,
				     const RealVect&  xi0_dst,
				     const IntVect&   iv0_dst,
				     const RealVect&  dx_dst,
				     const int        order) const;


   void getInterpolationCoefficients(Vector<Real>&   coeff,
				     const Real&     xiNW_dst,
				     const Real&     xiNE_dst,
				     const Real&     xiSW_dst,
				     const Real&     xiSE_dst,
				     const IntVect&  iv0_dst,
				     const RealVect& dx_dst) const;


#endif
};



#include "NamespaceFooter.H"

#endif
