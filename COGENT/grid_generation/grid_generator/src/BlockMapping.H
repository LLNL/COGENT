#ifndef _BLOCKMAPPING_H_
#define _BLOCKMAPPING_H_

#define USE_DCT_FIELD

#include "ParmParse.H"
#include "RealVect.H"
#include "FArrayBox.H"
#include "LevelData.H"

#include "NamespaceHeader.H"

class BlockMapping
{
  public:

   BlockMapping(const string& field_coefficients_file_name,
             const int     block_number,
             const double  psiNorm_bigEnd,
             const int     n_extrp,
             const double  trans_length,
             const double  blending_factor);

   ~BlockMapping();

   void init(const double a_trans_rad);

   void getPhysicalCoordinates(FArrayBox&                physical_coordinates,
                               const Vector<RealVect*>&  mapping) const;
   
   void getArcLengthMapping(Vector<RealVect*>&  mapping,
                            double&             block_length,
                            const int           block_id) const;

   Vector<Real> getMagFieldData(const RealVect& physical_coordinate ) const;
   
   double getMagneticFluxFromDCT(const RealVect& physical_coordinate ) const;

   RealVect findCriticalPoint(const RealVect& init_guess,
                              const double    init_step) const;
   
   void getMinimalCellSize(Real&             size,
                           RealVect&         location,
                           const FArrayBox&  coords,
                           const int         dir) const;
   
   int m_block_number;
   RealVect m_Xpoint, m_magAxis;

   double m_lengthCORE; //2*Pi poloidal arc length of the separatrix
   double m_nthetaCORE; //number of poloidal cell in the 2Pi core region

  private:

   void assembleDecomposedData(FArrayBox&                    data,
                               const LevelData<FArrayBox>&   data_decomp,
                               const Vector<Box>&            grid_boxes,
                               const int&                    decomp_num ) const;

   void getPoloidalDisjointBoxLayout(DisjointBoxLayout&   grids,
                                     Vector<Box>&         boxes,
                                     int&                 decomp_num,
                                     const Box&           box) const;
   
   void getMagneticFluxFromDCT(const FArrayBox& a_physical_coordinates,
                               FArrayBox&       a_magnetic_flux ) const;
   
   void getNormMagneticFlux(const FArrayBox& physical_coordinates,
                            FArrayBox&       magnetic_flux ) const;
   
   double getNormMagneticFlux(const RealVect& physical_coordinate ) const;
   
   
   double getLinearMagneticFlux( const RealVect& physical_coordinate ) const;

   
   void readFieldFromDCTFile(const string& file_name);
  
   void getNormFluxCoordinates(FArrayBox&       flux_coords) const;
   
   double getPoloidalArcLength(const int       pol_index,
                               const string&   option) const;
   
   void getRadialGridLine(FArrayBox&               radialGridLine,
                          const Vector<RealVect*>& mapping,
                          const double             length) const;
   
   void extrapolateGrid(FArrayBox& physical_coordinates, const Box& box_extrp) const;

   void radialRayTracing( RealVect&       X,
                          const double    psi_norm,
                          const double    arc_length,
                          const RealVect& cut_dir) const;
   
   void pushToFluxSurface(RealVect& X, const double psi_norm) const;
   
   void pushToFluxSurface(RealVect& X, const double psi_norm, const RealVect& dir) const;
   
   void initializeBlockCut(FArrayBox&  block_cut_coords,
                           FArrayBox&  block_cut_dir) const;
   
   void getCutDir(const Vector<Real>&  sep_dir,
                  RealVect&            core_cut_dir,
                  RealVect&            sol_cut_dir ) const;

   
   RealVect getStraightCutCoords( const double a_psi ) const;
   
   Vector <Real> getSeparatrixDirAnalytic() const;

   void getPsiLinearCoeff(RealVect& vec_1,
                          RealVect& vec_2 ) const;
   
   RealVect interpolateFromArcLengthMapping(const Vector<RealVect*>& mapping,
                                            const double             length) const;
   
   RealVect getRadialDir(const RealVect&  X,
                         const double     length,
                         const RealVect&  cut_dir) const;
   
   RealVect getPsiNorm(const RealVect& X,
                       const bool      norm_dir) const;
   
   RealVect getPsiTang(const RealVect& X,
                       const int       block_number) const;
   
   double vectProd(const RealVect& vect_1,
                   const RealVect& vect_2) const;
   
   double vectMag(const RealVect& vect) const;
   
   double projection(const RealVect& vec,
                     const RealVect& targ) const;
   
   RealVect rotation(const RealVect&   vec,
                     const int         dir) const;
   
   FArrayBox m_psi_coefs;
   
   mutable FArrayBox m_block_cut_dir;
   mutable FArrayBox m_block_cut_coords;
   
   double m_psi_smallEnd, m_psi_bigEnd;
   double m_physFluxOnAxis, m_physFluxOnSep;
   int m_nExtrp;

   double m_tol, m_dl;
   int m_max_iterations, m_psi_expansion_order;
   
   double m_Rmin, m_Rmax, m_Zmin, m_Zmax;
   RealVect m_core_cut_dir, m_sol_cut_dir, m_lin_1, m_lin_2;
   
   Vector<Real> m_sep_dir;

   bool m_dealigned_grid;
   bool m_straight_block_bndry;
   double m_trans_rad; //transition radius between the rectalinear and original flux functions
   double m_trans_length; // sets the region beyond which the generated grid is locally orthogonal
   double m_blending_factor; // scaling factor for blending the real and block-aligned flux 

#if 0
   Vector <Real> getSeparatrixDirNumeric() const;

   RealVect poloidaCoordRefPoint(const double a_theta,
                                 const double a_psi,
                                 const int    a_side ) const;
#endif
   
};

#include "NamespaceFooter.H"

#endif
