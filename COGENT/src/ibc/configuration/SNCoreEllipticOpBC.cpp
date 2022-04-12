#include "SNCoreEllipticOpBC.H"
#include "SNCoreBlockCoordSys.H"
#include "MagGeom.H"
#include "GridFunctionLibrary.H"
#include "Directions.H"
#include "DataArray.H"

#include "NamespaceHeader.H"



SNCoreEllipticOpBC::SNCoreEllipticOpBC(const int&  a_nblocks)
   : EllipticOpBC(NUM_BOUNDARIES, a_nblocks)
{
   setNames();
}



SNCoreEllipticOpBC::SNCoreEllipticOpBC( const std::string&  a_name,
                                        ParmParse&          a_pp,
                                        const int&          a_nblocks,
                                        const int&          a_verbosity )
   : EllipticOpBC(NUM_BOUNDARIES, a_nblocks),
     m_name(a_name),
     m_verbosity(a_verbosity)
{
   setNames();
   parseParameters( a_pp );
   
   if (hasCoupledBoundary()) {
      for (int i=0; i<m_bc_block_data.size(); ++i) {
         int verbosity = 0;
         m_bc_block_data[i] = RefCountedPtr<GridFunction>( new DataArray( verbosity ) );
      }
   }
}


SNCoreEllipticOpBC::SNCoreEllipticOpBC( const SingleNullEllipticOpBC&  a_bc,
                                        const int                      a_outer_radial_bc_type,
                                        const int                      a_num_blocks )
  : EllipticOpBC(NUM_BOUNDARIES, a_num_blocks)
{
   CH_assert(a_outer_radial_bc_type == DIRICHLET || a_outer_radial_bc_type == MAPPED_NEUMANN);
   setNames();

   // Use the same lower radial conditions as the initializing full single null
   m_bc_type[RADIAL_INNER] = a_bc.getBCType(SNCoreBlockCoordSys::LCORE, RADIAL_DIR, 0);
   m_bc_value[RADIAL_INNER] = a_bc.getBCValue(SNCoreBlockCoordSys::LCORE, RADIAL_DIR, 0);
   m_bc_function[RADIAL_INNER] = a_bc.getBCFunction(SNCoreBlockCoordSys::LCORE, RADIAL_DIR, 0);

   // Assume outer radial conditions provided by a grid function
   m_bc_type[RADIAL_OUTER] = a_outer_radial_bc_type;
   m_bc_function[RADIAL_OUTER] = RefCountedPtr<GridFunction>( new DataArray( false ) ); 
}


void
SNCoreEllipticOpBC::setNames()
{
   m_bdry_name[RADIAL_INNER] = "radial_inner";
   m_bdry_name[RADIAL_OUTER] = "radial_outer";
}



void
SNCoreEllipticOpBC::setBCType( const int  a_block_number,
                               const int  a_dir,
                               const int  a_side,
                               const int  a_type )
{
   CH_assert(a_side == 0 || a_side == 1);

   switch( a_block_number )
      {
      case SNCoreBlockCoordSys::MCORE:
      case SNCoreBlockCoordSys::LCORE:
      case SNCoreBlockCoordSys::RCORE:
         if ( a_dir == RADIAL_DIR ) {
            if ( a_side == 0 ) {
               m_bc_type[RADIAL_INNER] = a_type;
            }
            else if ( a_side == 1 ) {
               m_bc_type[RADIAL_OUTER] = a_type;
            }
            else {
               MayDay::Error("SNCoreEllipticOpBC::setBCType(): Invalid side argument");
            }
         }
         else {
            MayDay::Error("SNCoreEllipticOpBC::setBCType(): Invalid direction argument");
         }
         break;
      default:
         MayDay::Error("SNCoreEllipticOpBC::setBCType(): Unrecognized block number");
      }
}


int
SNCoreEllipticOpBC::getBCType( const int  a_block_number,
                               const int  a_dir,
                               const int  a_side ) const
{
   CH_assert(a_side == 0 || a_side == 1);
   int bc_type = UNDEFINED;

   switch( a_block_number )
      {
      case SNCoreBlockCoordSys::MCORE:
      case SNCoreBlockCoordSys::LCORE:
      case SNCoreBlockCoordSys::RCORE:
         if ( a_dir == RADIAL_DIR ) {
            if ( a_side == 0 ) {
               bc_type = m_bc_type[RADIAL_INNER];
            }
            else if ( a_side == 1 ) {
               bc_type = m_bc_type[RADIAL_OUTER];
            }
            else {
               MayDay::Error("SNCoreEllipticOpBC::getBCType(): Invalid side argument");
            }
         }
         break;
      default:
         MayDay::Error("SNCoreEllipticOpBC::getBCType(): Unrecognized block number");
      }

   return bc_type;
}

std::string
SNCoreEllipticOpBC::getBCSubType(const int  a_block_number,
                                 const int  a_dir,
                                 const int  a_side ) const
{
   CH_assert(a_side == 0 || a_side == 1);
   std::string bc_subtype;

   switch( a_block_number )
      {
      case SNCoreBlockCoordSys::MCORE:
      case SNCoreBlockCoordSys::LCORE:
      case SNCoreBlockCoordSys::RCORE:
         if ( a_dir == RADIAL_DIR ) {
            if ( a_side == 0 ) {
               bc_subtype = m_bc_subtype[RADIAL_INNER];
            }
            else if ( a_side == 1 ) {
               bc_subtype = m_bc_subtype[RADIAL_OUTER];
            }
            else {
               MayDay::Error("SNCoreEllipticOpBC::getBCSubType(): Invalid side argument");
            }
         }
         else {
            MayDay::Error("SNCoreEllipticOpBC::getBCSubType(): Invalid direction argument");
         }
         break;
      default:
         MayDay::Error("SNCoreEllipticOpBC::getBCSubType(): Unrecognized block number");
      }

   return bc_subtype;
}

void
SNCoreEllipticOpBC::setBCValue( const int     a_block_number,
                                const int     a_dir,
                                const int     a_side,
                                const double  a_value )
{
   CH_assert(a_side == 0 || a_side == 1);

   switch( a_block_number )
      {
      case SNCoreBlockCoordSys::MCORE:
      case SNCoreBlockCoordSys::LCORE:
      case SNCoreBlockCoordSys::RCORE:
         if ( a_dir == RADIAL_DIR ) {
            if ( a_side == 0 ) {
               m_bc_value[RADIAL_INNER] = a_value;
            }
            else if ( a_side == 1 ) {
               m_bc_value[RADIAL_OUTER] = a_value;
            }
            else {
               MayDay::Error("SNCoreEllipticOpBC::setBCValue(): Invalid side argument");
            }
         }
         else {
            MayDay::Error("SNCoreEllipticOpBC::setBCValue(): Invalid direction argument");
         }
         break;
      default:
         MayDay::Error("SNCoreEllipticOpBC::setBCValue(): Unrecognized block number");
      }
}



double
SNCoreEllipticOpBC::getBCValue( const int  a_block_number,
                                const int  a_dir,
                                const int  a_side ) const
{
   CH_assert(a_side == 0 || a_side == 1);
   double bc_value = BASEFAB_REAL_SETVAL;

   switch( a_block_number )
      {
      case SNCoreBlockCoordSys::MCORE:
      case SNCoreBlockCoordSys::LCORE:
      case SNCoreBlockCoordSys::RCORE:
         if ( a_dir == RADIAL_DIR ) {
            if ( a_side == 0 ) {
               bc_value = m_bc_value[RADIAL_INNER];
            }
            else if ( a_side == 1 ) {
               bc_value = m_bc_value[RADIAL_OUTER];
            }
            else {
               MayDay::Error("SNCoreEllipticOpBC::getBCValue(): Invalid side argument");
            }
         }
         else {
            MayDay::Error("SNCoreEllipticOpBC::getBCValue(): Invalid direction argument");
         }
         break;
      default:
         MayDay::Error("SNCoreEllipticOpBC::getBCValue(): Unrecognized block number");
      }

   return bc_value;
}



void
SNCoreEllipticOpBC::setBCFunction( const int                           a_block_number,
                                   const int                           a_dir,
                                   const int                           a_side,
                                   const RefCountedPtr<GridFunction>&  a_function )
{
   CH_assert(a_side == 0 || a_side == 1);

   switch( a_block_number )
      {
      case SNCoreBlockCoordSys::MCORE:
      case SNCoreBlockCoordSys::LCORE:
      case SNCoreBlockCoordSys::RCORE:
         if ( a_dir == RADIAL_DIR ) {
            if ( a_side == 0 ) {
               m_bc_function[RADIAL_INNER] = a_function;
            }
            else if ( a_side == 1 ) {
               m_bc_function[RADIAL_OUTER] = a_function;
            }
            else {
               MayDay::Error("SNCoreEllipticOpBC::setBCFunction(): Invalid side argument");
            }
         }
         else {
            MayDay::Error("SNCoreEllipticOpBC::setBCFunction(): Invalid direction argument");
         }
         break;
      default:
         MayDay::Error("SNCoreEllipticOpBC::setBCFunction(): Unrecognized block number");
      }
}



RefCountedPtr<GridFunction>
SNCoreEllipticOpBC::getBCFunction( const int  a_block_number,
                                   const int  a_dir,
                                   const int  a_side ) const
{
   CH_assert(a_side == 0 || a_side == 1);
   RefCountedPtr<GridFunction> function;

   if (getBCSubType(a_block_number, a_dir, a_side) == "coupled" ) {
      function = getBlockBCData(a_block_number, a_dir, a_side);
   }
   else {
      switch( a_block_number )
      {
         case SNCoreBlockCoordSys::MCORE:
         case SNCoreBlockCoordSys::LCORE:
         case SNCoreBlockCoordSys::RCORE:
            if ( a_dir == RADIAL_DIR ) {
               if ( a_side == 0 ) {
                  function = m_bc_function[RADIAL_INNER];
               }
               else if ( a_side == 1 ) {
                  function = m_bc_function[RADIAL_OUTER];
               }
               else {
                  MayDay::Error("SNCoreEllipticOpBC::getBCFunction(): Invalid side argument");
               }
            }
            else {
               MayDay::Error("SNCoreEllipticOpBC::getBCFunction(): Invalid direction argument");
            }
            break;
         default:
            MayDay::Error("SNCoreEllipticOpBC::getBCFunction(): Unrecognized block number");
      }
   }
   return function;
}



void
SNCoreEllipticOpBC::apply( const MultiBlockLevelGeom&  a_geom,
                           const Box&                  a_coord_sys_box,
                           const double&               a_time,
                           const int                   a_dir,
                           const int                   a_side,
                           FArrayBox&                  a_phi ) const
{
   RefCountedPtr<GridFunction> function;
   double value;

   const MagCoordSys* mag_coord_sys = ((MagGeom&)a_geom).getCoordSys();
   int block_number = mag_coord_sys->whichBlock(a_coord_sys_box);

   switch( block_number )
      {
      case SNCoreBlockCoordSys::MCORE:
      case SNCoreBlockCoordSys::LCORE:
      case SNCoreBlockCoordSys::RCORE:
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
               MayDay::Error("SNCoreEllipticOpBC::apply(): Invalid side argument");
            }
         }
         else {
            MayDay::Error("SNCoreEllipticOpBC::apply(): Invalid direction argument");
         }
         break;
      default:
         MayDay::Error("SNCoreEllipticOpBC::apply(): Unrecognized block number");
      }

   if ( !function.isNull() ) {
      // Get the block id
      const MultiBlockCoordSys& coord_sys( *(a_geom.coordSysPtr()) );
      const int block_number( coord_sys.whichBlock( a_coord_sys_box ) );
      FArrayBox dummy;
      function->assign(a_phi, a_geom, dummy, dummy, block_number, a_time, false);
   }
   else {
      a_phi.setVal(value);
   }
}



void SNCoreEllipticOpBC::printParameters() const
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "SNCoreEllipticOpBC ================================" << std::endl;
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
void SNCoreEllipticOpBC::parseParameters( ParmParse& a_pp )
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
      else if (bc_type == "mapped_neumann") {
         m_bc_type[i] = MAPPED_NEUMANN;
      }
      else if (bc_type == "natural") {
         m_bc_type[i] = NATURAL;
      }
      else {
         MayDay::Error("SNCoreEllipticOpBC::parseParameter(): Unrecognized potential bc type");
      }

      fpp.query( "subtype", m_bc_subtype[i] );
      
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
         MayDay::Error("SNCoreEllipticOpBC::parseParameters(): Please specify either a value or a function, but not both");
      }
   }

   if (m_verbosity) {
      printParameters();
   }
}


RefCountedPtr<EllipticOpBC>
SNCoreEllipticOpBC::clone( const bool a_extrapolated ) const
{
   // Also sets m_bdry_name
   RefCountedPtr<SNCoreEllipticOpBC> result
      = RefCountedPtr<SNCoreEllipticOpBC>(new SNCoreEllipticOpBC(m_num_blocks));

   // Copy the data members set by the full constructor (i.e., the one with the ParmParse
   // argument).  
   result->m_name = m_name;
   result->m_verbosity = m_verbosity;

   // Copy the rest of the data owned by the EllipticOpBC base class.
   result->copyBaseData(*this);

   if ( a_extrapolated ) {
      result->setExtrapolatedType();
   }

   return result;
}


#include "NamespaceFooter.H"
