#ifndef _SINGLENULLPOTENTIALBC_
#define _SINGLENULLPOTENTIALBC_

#include "ParmParse.H"
#include "RefCountedPtr.H"
#include "PotentialBC.H"
#include "GridFunctionLibrary.H"

#include "NamespaceHeader.H"


class SingleNullPotentialBC
   : public PotentialBC
{
public:

   SingleNullPotentialBC();

   SingleNullPotentialBC( const std::string& name,
                          ParmParse& pp,
                          const int& verbosity );

   virtual ~SingleNullPotentialBC() {};

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

private:

   // prohibit copying
   SingleNullPotentialBC( const SingleNullPotentialBC& );
   SingleNullPotentialBC& operator=( const SingleNullPotentialBC& );

   void setNames();

   inline void parseParameters( ParmParse& a_pp );

   std::string m_name;
   int m_verbosity;

   enum {INVALID=-1,
         RADIAL_CORE, RADIAL_SOL, RADIAL_PF,
         POLOIDAL_INNER_DIV, POLOIDAL_OUTER_DIV,
         NUM_BOUNDARIES};
};

#include "NamespaceFooter.H"

#endif
