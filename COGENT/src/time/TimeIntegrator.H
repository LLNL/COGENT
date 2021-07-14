#ifndef _TIME_INTEGRATOR_H_
#define _TIME_INTEGRATOR_H_

#include <string>

#include "REAL.H"
#include "ParmParse.H"
#include "parstream.H"
#include "TiDefinitions.H"

#include "NamespaceHeader.H"

template <class Vec, class Ops>
class TimeIntegrator
{
   
public:

   virtual ~TimeIntegrator<Vec,Ops>() {}
   
   virtual void define( ParmParse&,
                        const std::string&,
                        const Vec&,
                        Ops* const ) = 0;

   virtual void initialize( const Real, const Vec& ) = 0;

   virtual void advance( const Real& a_time, Vec& Y ) = 0;
   
   virtual Ops& getOperators() = 0;
   
   virtual void setTimeStep( const int& a_n ) = 0;

   virtual void getTimeStep ( int& a_n ) const = 0;

   virtual void setTimeStepSize( const Real& a_dt) = 0;

   virtual void getTimeStepSize( Real& a_dt) const = 0;

   virtual void setCurrentTime( const Real& a_time ) = 0;

   virtual void getCurrentTime( Real& a_time ) const = 0;

   virtual TimeIntegratorType type() const = 0;
   
   virtual void printCounts() const = 0;

};

#include "NamespaceFooter.H"

#endif
