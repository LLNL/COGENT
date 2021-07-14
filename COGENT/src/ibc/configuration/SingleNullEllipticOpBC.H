#ifndef _SINGLENULLELLIPTICOPBC_
#define _SINGLENULLELLIPTICOPBC_

#include "ParmParse.H"
#include "RefCountedPtr.H"
#include "EllipticOpBC.H"
#include "GridFunctionLibrary.H"

#include "NamespaceHeader.H"


class SingleNullEllipticOpBC
   : public EllipticOpBC
{
public:

   SingleNullEllipticOpBC(const int& nblocks);

   SingleNullEllipticOpBC( const std::string&  name,
                           ParmParse&          pp,
                           const int&          nblocks,
                           const int&          poloidal_blocks_per_sector,
                           const int&          verbosity );

   virtual ~SingleNullEllipticOpBC() {};

   virtual void setBCType( const int block_number,
                           const int dir,
                           const int side,
                           const int type );

   virtual int getBCType( const int block_number,
                          const int dir,
                          const int side ) const;
   
   virtual std::string getBCSubType(const int block_number,
                                    const int dir,
                                    const int side ) const;

   virtual void setBCValue( const int block_number,
                            const int dir,
                            const int side,
                            const double value );

   virtual double getBCValue( const int block_number,
                              const int dir,
                              const int side ) const;

   virtual void setBCFunction( const int block_number,
                               const int dir,
                               const int side,
                               const RefCountedPtr<GridFunction>& function );

   virtual RefCountedPtr<GridFunction> getBCFunction( const int block_number,
                                                      const int dir,
                                                      const int side ) const;

   virtual void apply( const MultiBlockLevelGeom& a_geom,
                       const Box&                 a_coord_sys_box,
                       const double&              a_time,
                       const int                  a_dir,
                       const int                  a_side,
                       FArrayBox&                 a_phi ) const;

   virtual RefCountedPtr<EllipticOpBC> clone( const bool extrapolated = false ) const;

   virtual void printParameters() const;

   RefCountedPtr<GridFunction> getCoreFunction() {return m_bc_function[RADIAL_CORE];}
   const RefCountedPtr<GridFunction> getCoreFunction() const {return m_bc_function[RADIAL_CORE];}

   RefCountedPtr<GridFunction> getSOLFunction() {return m_bc_function[RADIAL_SOL];}
   const RefCountedPtr<GridFunction> getSOLFunction() const {return m_bc_function[RADIAL_SOL];}
   
   RefCountedPtr<GridFunction> getPFFunction() {return m_bc_function[RADIAL_PF];}
   const RefCountedPtr<GridFunction> getPFFunction() const {return m_bc_function[RADIAL_PF];}

   RefCountedPtr<GridFunction> getInnerDivertorFunction() {return m_bc_function[POLOIDAL_INNER_DIV];}
   const RefCountedPtr<GridFunction> getInnerDivertorFunction() const {return m_bc_function[POLOIDAL_INNER_DIV];}

   RefCountedPtr<GridFunction> getOuterDivertorFunction() {return m_bc_function[POLOIDAL_OUTER_DIV];}
   const RefCountedPtr<GridFunction> getOuterDivertorFunction() const {return m_bc_function[POLOIDAL_OUTER_DIV];}

   void setPoloidalBlocksPerSector(const int& poloidal_blocks_per_sector) {m_poloidal_blocks_per_sector = poloidal_blocks_per_sector;}

private:

   // prohibit copying
   SingleNullEllipticOpBC( const SingleNullEllipticOpBC& );
   SingleNullEllipticOpBC& operator=( const SingleNullEllipticOpBC& );

   void setNames();

   inline void parseParameters( ParmParse& a_pp );

   std::string m_name;
   int m_poloidal_blocks_per_sector;
   int m_verbosity;

   enum {INVALID=-1,
         RADIAL_CORE, RADIAL_SOL, RADIAL_PF,
         POLOIDAL_INNER_DIV, POLOIDAL_OUTER_DIV,
#if CFG_DIM==3
         TOROIDAL_CORE, TOROIDAL_SOL, TOROIDAL_PF,
         TOROIDAL_INNER_DIV, TOROIDAL_OUTER_DIV,
#endif
         NUM_BOUNDARIES};
};

#include "NamespaceFooter.H"

#endif
