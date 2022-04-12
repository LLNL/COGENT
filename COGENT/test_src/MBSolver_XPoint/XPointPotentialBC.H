#ifndef _XPOINTPOTENTIALBC_
#define _XPOINTPOTENTIALBC_

#include "ParmParse.H"
#include "RefCountedPtr.H"
//#include "PotentialBC.H"
#include "EllipticOpBC.H"

#include "NamespaceHeader.H"


class XPointPotentialBC
   : public EllipticOpBC
{
public:

   XPointPotentialBC();

   XPointPotentialBC( const std::string& name,
                          ParmParse& pp,
                          const int& verbosity );

   virtual ~XPointPotentialBC() {};

   virtual void setBCType( const int block_number,
                           const int dir,
                           const int side,
                           const int type );

   virtual int getBCType( const int block_number,
                          const int dir,
                          const int side ) const;

   // STUB
   virtual std::string getBCSubType(const int block_number,
                                    const int dir,
                                    const int side ) const {return "";}

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

   // STUB
   virtual RefCountedPtr<EllipticOpBC> clone( const bool extrapolated = false ) const {return RefCountedPtr<EllipticOpBC>(NULL);}

   virtual void printParameters() const;

private:

   // prohibit copying
   XPointPotentialBC( const XPointPotentialBC& );
   XPointPotentialBC& operator=( const XPointPotentialBC& );

   void setNames();

   inline void parseParameters( ParmParse& a_pp );

   std::string m_name;
   int m_verbosity;

   enum {INVALID=-1,
         BLOCK_0, BLOCK_1, BLOCK_2, BLOCK_3, BLOCK_4, BLOCK_5, BLOCK_6, BLOCK_7,
         NUM_BOUNDARIES};
};

#include "NamespaceFooter.H"

#endif
