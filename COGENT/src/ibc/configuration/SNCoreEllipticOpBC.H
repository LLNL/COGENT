#ifndef _SNCOREELLIPTICOPBC_
#define _SNCOREELLIPTICOPBC_

#include "ParmParse.H"
#include "RefCountedPtr.H"
#include "EllipticOpBC.H"
#include "GridFunctionLibrary.H"
#include "SingleNullEllipticOpBC.H"

#include "NamespaceHeader.H"


class SNCoreEllipticOpBC
   : public EllipticOpBC
{
public:

   SNCoreEllipticOpBC();

   SNCoreEllipticOpBC( const std::string&  name,
                       ParmParse&          pp,
                       const int&          verbosity );

   SNCoreEllipticOpBC( const SingleNullEllipticOpBC&  bc,
                       const int                      outer_radial_bc_type );

   virtual ~SNCoreEllipticOpBC() {};

   virtual void setBCType( const int block_number,
                           const int dir,
                           const int side,
                           const int type );

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

   virtual void printParameters() const;

   RefCountedPtr<GridFunction> getRadialInnerFunction() {return m_bc_function[RADIAL_INNER];}

   RefCountedPtr<GridFunction> getRadialOuterFunction() {return m_bc_function[RADIAL_OUTER];}

private:

   // prohibit copying
   SNCoreEllipticOpBC( const SNCoreEllipticOpBC& );
   SNCoreEllipticOpBC& operator=( const SNCoreEllipticOpBC& );

   void setNames();

   inline void parseParameters( ParmParse& a_pp );

   std::string m_name;
   int m_verbosity;

   enum {INVALID=-1,
         RADIAL_INNER, RADIAL_OUTER,
         NUM_BOUNDARIES};
};

#include "NamespaceFooter.H"

#endif
