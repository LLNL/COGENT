#ifndef _SLABPOTENTIALBC_
#define _SLABPOTENTIALBC_

#include "ParmParse.H"
#include "RefCountedPtr.H"
#include "PotentialBC.H"

#include "NamespaceHeader.H"


class SlabPotentialBC
   : public PotentialBC
{
public:

   SlabPotentialBC();

   SlabPotentialBC( const std::string& name,
                    ParmParse& pp,
                    const int& verbosity );

   virtual ~SlabPotentialBC() {};

   virtual void setBCType( const int block_number,
                           const int dir,
                           const int side,
                           const int type);

   virtual int getBCType( const int block_number,
                          const int dir,
                          const int side ) const;

   virtual void setBCValue( const int block_number,
                            const int dir,
                            const int side,
                            const double value );

   virtual double getBCValue( const int block_number,
                              const int dir,
                              const int side ) const;

   virtual void setBCFunction( const int                          block_number,
                               const int                          dir,
                               const int                          side,
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

   virtual void printParameters() const;

private:

   // prohibit copying
   SlabPotentialBC( const SlabPotentialBC& );
   SlabPotentialBC& operator=( const SlabPotentialBC& );

   void setNames();

   inline void parseParameters( ParmParse& a_pp );

   std::string m_name;
   int m_verbosity;

   enum {INVALID=-1,
         RADIAL_LOWER, RADIAL_UPPER,
         POLOIDAL_LOWER, POLOIDAL_UPPER,
#if CFG_DIM==3
         TOROIDAL_LOWER, TOROIDAL_UPPER,
#endif
         NUM_BOUNDARIES};
};

#include "NamespaceFooter.H"

#endif
