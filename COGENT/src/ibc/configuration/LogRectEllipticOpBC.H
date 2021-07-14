#ifndef _LOGRECTELLIPTICOPBC_
#define _LOGRECTELLIPTICOPBC_

#include "ParmParse.H"
#include "RefCountedPtr.H"
#include "EllipticOpBC.H"

#include "NamespaceHeader.H"


class LogRectEllipticOpBC
   : public EllipticOpBC
{
public:

   LogRectEllipticOpBC(const int& nblocks);

   LogRectEllipticOpBC( const std::string& name,
                        ParmParse& pp,
                        const int& nblocks,
                        const int& verbosity );

   virtual ~LogRectEllipticOpBC() {};

   virtual void setBCType( const int block_number,
                           const int dir,
                           const int side,
                           const int type);

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

   virtual RefCountedPtr<EllipticOpBC> clone( const bool extrapolated = 0 ) const;

   virtual void printParameters() const;

private:

   // prohibit copying
   LogRectEllipticOpBC( const LogRectEllipticOpBC& );
   LogRectEllipticOpBC& operator=( const LogRectEllipticOpBC& );

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
