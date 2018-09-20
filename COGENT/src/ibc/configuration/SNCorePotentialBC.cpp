#include "SNCorePotentialBC.H"
#include "SNCoreBlockCoordSys.H"
#include "SingleNullBlockCoordSys.H"
#include "GridFunctionLibrary.H"
#include "MagGeom.H"
#include "Directions.H"
#include "DataArray.H"

#include "NamespaceHeader.H"



SNCorePotentialBC::SNCorePotentialBC()
   : PotentialBC(NUM_BOUNDARIES)
{
   setNames();
}



SNCorePotentialBC::SNCorePotentialBC( const std::string& a_name,
                                      ParmParse&         a_pp,
                                      const int&         a_verbosity )
   : PotentialBC(NUM_BOUNDARIES),
     m_name(a_name),
     m_verbosity(a_verbosity)
{
   setNames();
   parseParameters( a_pp );
}



SNCorePotentialBC::SNCorePotentialBC( const SingleNullPotentialBC& a_bc,
                                      const int                    a_outer_radial_bc_type )
   : PotentialBC(NUM_BOUNDARIES)
{
   CH_assert(a_outer_radial_bc_type == DIRICHLET || a_outer_radial_bc_type == NEUMANN);
   setNames();

   // Use the same lower radial conditions as the initializing full single null
   m_bc_type[RADIAL_INNER] = a_bc.getBCType(LCORE, RADIAL_DIR, 0);
   m_bc_value[RADIAL_INNER] = a_bc.getBCValue(LCORE, RADIAL_DIR, 0);
   m_bc_function[RADIAL_INNER] = a_bc.getBCFunction(LCORE, RADIAL_DIR, 0);

   // Assume outer radial conditions provided by a grid function
   m_bc_type[RADIAL_OUTER] = a_outer_radial_bc_type;
   m_bc_function[RADIAL_OUTER] = RefCountedPtr<GridFunction>( new DataArray( false ) ); 
}


void
SNCorePotentialBC::setNames()
{
   m_bdry_name[RADIAL_INNER] = "radial_inner";
   m_bdry_name[RADIAL_OUTER] = "radial_outer";
}



void
SNCorePotentialBC::setBCType( const int a_block_number,
                              const int a_dir,
                              const int a_side,
                              const int a_type )
{
   CH_assert(a_side == 0 || a_side == 1);

   switch( a_block_number )
      {
      case M_CORE:
      case L_CORE:
      case R_CORE:
         if ( a_dir == RADIAL_DIR ) {
            if ( a_side == 0 ) {
               m_bc_type[RADIAL_INNER] = a_type;
            }
            else if ( a_side == 1 ) {
               m_bc_type[RADIAL_OUTER] = a_type;
            }
            else {
               MayDay::Error("SNCorePotentialBC::setBCType(): Invalid side argument");
            }
         }
         else {
            MayDay::Error("SNCorePotentialBC::setBCType(): Invalid direction argument");
         }
         break;
      default:
         MayDay::Error("SNCorePotentialBC::setBCType(): Unrecognized block number");
      }
}


int
SNCorePotentialBC::getBCType( const int a_block_number,
                              const int a_dir,
                              const int a_side ) const
{
   CH_assert(a_side == 0 || a_side == 1);
   int bc_type = UNDEFINED;

   switch( a_block_number )
      {
      case M_CORE:
      case L_CORE:
      case R_CORE:
         if ( a_dir == RADIAL_DIR ) {
            if ( a_side == 0 ) {
               bc_type = m_bc_type[RADIAL_INNER];
            }
            else if ( a_side == 1 ) {
               bc_type = m_bc_type[RADIAL_OUTER];
            }
            else {
               MayDay::Error("SNCorePotentialBC::getBCType(): Invalid side argument");
            }
         }
         else {
            MayDay::Error("SNCorePotentialBC::getBCType(): Invalid direction argument");
         }
         break;
      default:
         MayDay::Error("SNCorePotentialBC::getBCType(): Unrecognized block number");
      }

   return bc_type;
}



void
SNCorePotentialBC::setBCValue( const int    a_block_number,
                               const int    a_dir,
                               const int    a_side,
                               const double a_value )
{
   CH_assert(a_side == 0 || a_side == 1);

   switch( a_block_number )
      {
      case M_CORE:
      case L_CORE:
      case R_CORE:
         if ( a_dir == RADIAL_DIR ) {
            if ( a_side == 0 ) {
               m_bc_value[RADIAL_INNER] = a_value;
            }
            else if ( a_side == 1 ) {
               m_bc_value[RADIAL_OUTER] = a_value;
            }
            else {
               MayDay::Error("SNCorePotentialBC::setBCValue(): Invalid side argument");
            }
         }
         else {
            MayDay::Error("SNCorePotentialBC::setBCValue(): Invalid direction argument");
         }
         break;
      default:
         MayDay::Error("SNCorePotentialBC::setBCValue(): Unrecognized block number");
      }
}



double
SNCorePotentialBC::getBCValue( const int a_block_number,
                               const int a_dir,
                               const int a_side ) const
{
   CH_assert(a_side == 0 || a_side == 1);
   double bc_value = BASEFAB_REAL_SETVAL;

   switch( a_block_number )
      {
      case M_CORE:
      case L_CORE:
      case R_CORE:
         if ( a_dir == RADIAL_DIR ) {
            if ( a_side == 0 ) {
               bc_value = m_bc_value[RADIAL_INNER];
            }
            else if ( a_side == 1 ) {
               bc_value = m_bc_value[RADIAL_OUTER];
            }
            else {
               MayDay::Error("SNCorePotentialBC::getBCValue(): Invalid side argument");
            }
         }
         else {
            MayDay::Error("SNCorePotentialBC::getBCValue(): Invalid direction argument");
         }
         break;
      default:
         MayDay::Error("SNCorePotentialBC::getBCValue(): Unrecognized block number");
      }

   return bc_value;
}



void
SNCorePotentialBC::setBCFunction( const int                          a_block_number,
                                  const int                          a_dir,
                                  const int                          a_side,
                                  const RefCountedPtr<GridFunction>& a_function )
{
   CH_assert(a_side == 0 || a_side == 1);

   switch( a_block_number )
      {
      case M_CORE:
      case L_CORE:
      case R_CORE:
         if ( a_dir == RADIAL_DIR ) {
            if ( a_side == 0 ) {
               m_bc_function[RADIAL_INNER] = a_function;
            }
            else if ( a_side == 1 ) {
               m_bc_function[RADIAL_OUTER] = a_function;
            }
            else {
               MayDay::Error("SNCorePotentialBC::setBCFunction(): Invalid side argument");
            }
         }
         else {
            MayDay::Error("SNCorePotentialBC::setBCFunction(): Invalid direction argument");
         }
         break;
      default:
         MayDay::Error("SNCorePotentialBC::setBCFunction(): Unrecognized block number");
      }
}



RefCountedPtr<GridFunction>
SNCorePotentialBC::getBCFunction( const int a_block_number,
                                  const int a_dir,
                                  const int a_side ) const
{
   CH_assert(a_side == 0 || a_side == 1);
   RefCountedPtr<GridFunction> function;

   switch( a_block_number )
      {
      case M_CORE:
      case L_CORE:
      case R_CORE:
         if ( a_dir == RADIAL_DIR ) {
            if ( a_side == 0 ) {
               function = m_bc_function[RADIAL_INNER];
            }
            else if ( a_side == 1 ) {
               function = m_bc_function[RADIAL_OUTER];
            }
            else {
               MayDay::Error("SNCorePotentialBC::getBCFunction(): Invalid side argument");
            }
         }
         else {
            MayDay::Error("SNCorePotentialBC::getBCFunction(): Invalid direction argument");
         }
         break;
      default:
         MayDay::Error("SNCorePotentialBC::getBCFunction(): Unrecognized block number");
      }

   return function;
}



void
SNCorePotentialBC::apply( const MultiBlockLevelGeom& a_geom,
                          const Box&                 a_coord_sys_box,
                          const double&              a_time,
                          const int                  a_dir,
                          const int                  a_side,
                          FArrayBox&                 a_phi ) const
{
   RefCountedPtr<GridFunction> function;
   double value;

   const MagCoordSys* mag_coord_sys = ((MagGeom&)a_geom).getCoordSys();
   int block_number = mag_coord_sys->whichBlock(a_coord_sys_box);

   switch( block_number )
      {
      case M_CORE:
      case L_CORE:
      case R_CORE:
         if ( a_dir == RADIAL_DIR ) {
            if ( a_side == 0 ) {
               function = m_bc_function[RADIAL_INNER];
               value = m_bc_value[RADIAL_INNER];
            }
            else if ( a_side == 1 ) {
               function = m_bc_function[RADIAL_OUTER];
               value = m_bc_value[RADIAL_OUTER];
            }
            else {
               MayDay::Error("SNCorePotentialBC::apply(): Invalid side argument");
            }
         }
         else {
            MayDay::Error("SNCorePotentialBC::apply(): Invalid direction argument");
         }
         break;
      default:
         MayDay::Error("SNCorePotentialBC::apply(): Unrecognized block number");
      }

   if ( !function.isNull() ) {
      function->assign(a_phi, a_geom, a_coord_sys_box, a_time, false);
   }
   else {
      a_phi.setVal(value);
   }
}



void SNCorePotentialBC::printParameters() const
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "SNCorePotentialBC ================================" << std::endl;
      std::cout << "- variable: "  << m_name << "-------------" << std::endl;
      for (int i(0); i<m_bc_function.size(); i++) {
         std::cout << "  " << m_bdry_name[i] << ": " << std::endl;
         if (m_bc_function[i]) m_bc_function[i]->printParameters();
         std::cout << "     bc_type  = " << m_bc_type[i] << std::endl;
         std::cout << "     bc_value = " << m_bc_value[i] << std::endl;
      }
      std::cout << "-----------------------------------------------" << std::endl;
      std::cout << "===============================================" << std::endl;
   }
}



inline
void SNCorePotentialBC::parseParameters( ParmParse& a_pp )
{
   GridFunctionLibrary* library = GridFunctionLibrary::getInstance();

   for (int i(0); i<m_bc_type.size(); i++) {
      std::string prefix( a_pp.prefix() );
      prefix += "." + m_bdry_name[i];
      ParmParse fpp( prefix.c_str() );
      std::string bc_type;
      fpp.query( "type", bc_type );

      if (bc_type == "dirichlet") {
         m_bc_type[i] = DIRICHLET;
      }
      else if (bc_type == "neumann") {
         m_bc_type[i] = NEUMANN;
      }
      else {
         MayDay::Error("SNCorePotentialBC::parseParameter(): Unrecognized potential bc type");
      }

      bool value_specified = fpp.contains("value");

      if (value_specified) {
         fpp.query( "value", m_bc_value[i] );
      }

      bool function_specified = fpp.contains("function");

      if (function_specified) {
         std::string function_name;
         fpp.query( "function", function_name );
         m_bc_function[i] = library->find( function_name );
      }

      if (value_specified && function_specified) {
         MayDay::Error("SNCorePotentialBC::parseParameters(): Please specify either a value or a function, but not both");
      }
   }

   if (m_verbosity) {
      printParameters();
   }
}



#include "NamespaceFooter.H"
