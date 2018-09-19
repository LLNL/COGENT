#include "SNSystemIC.H"

#include "Box.H"
#include "FourthOrderUtil.H"
#include "GridFunctionLibrary.H"

#include "NamespaceHeader.H"

SNSystemIC::SNSystemIC( ParmParse& a_pp,
                        const MagGeom& a_geometry,
                        const LevelData<FArrayBox>& a_solution )
   : m_mag_geometry( a_geometry ),
     m_nghosts(2),
     m_verbosity(0)
{
   a_pp.query( "snsystem.verbosity", m_verbosity );
   parseInput( a_pp, a_solution );
}


SNSystemIC::~SNSystemIC()
{
}


void SNSystemIC::initializeSystem( LevelData<FArrayBox>& a_solution,
                                   const Real& a_time ) const
{
   const GridFunction& ic( initialCondition( "solution" ) );
   ic.assign( a_solution, m_mag_geometry, a_time );
//   const MagGeom& mag_geometry( static_cast<const MagGeom&>(a_geometry) );
   m_mag_geometry.multJonValid( a_solution );
   a_solution.exchange();
}


const GridFunction& SNSystemIC::initialCondition( const std::string& a_name ) const
{
   ICMap::const_iterator it = m_ics.find( a_name );
   if ( it == m_ics.end() ) {
      const std::string msg( "SNSystemIC: Initial condition " + a_name + " not found!" );
      MayDay::Error( msg.c_str() );
   }
   return *(it->second);
}


void SNSystemIC::parseInput( ParmParse& a_pp,
                             const LevelData<FArrayBox>& a_solution )
{
   GridFunctionLibrary* library = GridFunctionLibrary::getInstance();
   ParmParse ppsp( "IC" );
   std::string function_name;
   ppsp.query( "function", function_name );
   RefCountedPtr<GridFunction> ic( library->find( function_name ) );
   m_ics.insert( ICMap::value_type( "solution", ic ) );
}


#include "NamespaceFooter.H"
