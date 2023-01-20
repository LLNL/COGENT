#include "LogRectEllipticOpBC.H"
#include "Directions.H"
#include "GridFunctionLibrary.H"
#include "DataArray.H"

#include "NamespaceHeader.H"


LogRectEllipticOpBC::LogRectEllipticOpBC(const int& a_nblocks)
   : EllipticOpBC(NUM_BOUNDARIES, a_nblocks)
{
   setNames();
}


LogRectEllipticOpBC::LogRectEllipticOpBC( const std::string&  a_name,
                                          ParmParse&          a_pp,
                                          const int           a_nblocks,
                                          const bool*         a_is_periodic,
                                          const int&          a_verbosity )
   : EllipticOpBC(NUM_BOUNDARIES, a_nblocks),
     m_name(a_name),
     m_verbosity(a_verbosity)
{
   setNames();
   parseParameters( a_pp, a_is_periodic );
   
   if (hasCoupledBoundary()) {
      for (int i=0; i<m_bc_block_data.size(); ++i) {
         int verbosity = 0;
         m_bc_block_data[i] = RefCountedPtr<GridFunction>( new DataArray( verbosity ) );
      }
   }
}


void
LogRectEllipticOpBC::setNames()
{
   m_bdry_name[RADIAL_LOWER] = "radial_lower";
   m_bdry_name[RADIAL_UPPER] = "radial_upper";
   m_bdry_name[POLOIDAL_LOWER] = "poloidal_lower";
   m_bdry_name[POLOIDAL_UPPER] = "poloidal_upper";
#if CFG_DIM==3
   m_bdry_name[TOROIDAL_LOWER] = "toroidal_lower";
   m_bdry_name[TOROIDAL_UPPER] = "toroidal_upper";
#endif
}



void
LogRectEllipticOpBC::setBCType( const int  a_block_number,
                                const int  a_dir,
                                const int  a_side,
                                const int  a_type )
{
   if ( a_dir == RADIAL_DIR ) {
      if ( a_side == 0 ) {
         m_bc_type[RADIAL_LOWER] = a_type;
      }
      else if ( a_side == 1 ) {
         m_bc_type[RADIAL_UPPER] = a_type;
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::setBCType(): Invalid side argument");
      }
   }
   else if ( a_dir == POLOIDAL_DIR ) {
      if ( a_side == 0 ) {
         m_bc_type[POLOIDAL_LOWER] = a_type;
      }
      else if ( a_side == 1 ) {
         m_bc_type[POLOIDAL_UPPER] = a_type;
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::setBCType(): Invalid side argument");
      }
   }
#if CFG_DIM==3
   else if ( a_dir == TOROIDAL_DIR ) {
      if ( a_side == 0 ) {
         m_bc_type[TOROIDAL_LOWER] = a_type;
      }
      else if ( a_side == 1 ) {
         m_bc_type[TOROIDAL_UPPER] = a_type;
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::setBCType(): Invalid side argument");
      }
   }
#endif
   else {
      MayDay::Error("LogRectEllipticOpBC::setBCType(): Invalid direction argument");
   }
}



int
LogRectEllipticOpBC::getBCType( const int  a_block_number,
                                const int  a_dir,
                                const int  a_side ) const
{
   int bc_type = UNDEFINED;

   if ( a_dir == RADIAL_DIR ) {
      if ( a_side == 0 ) {
         bc_type = m_bc_type[RADIAL_LOWER];
      }
      else if ( a_side == 1 ) {
         bc_type = m_bc_type[RADIAL_UPPER];
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::getBCType(): Invalid side argument");
      }
   }
   else if ( a_dir == POLOIDAL_DIR ) {
      if ( a_side == 0 ) {
         bc_type = m_bc_type[POLOIDAL_LOWER];
      }
      else if ( a_side == 1 ) {
         bc_type = m_bc_type[POLOIDAL_UPPER];
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::getBCType(): Invalid side argument");
      }
   }
#if CFG_DIM==3
   else if ( a_dir == TOROIDAL_DIR ) {
      if ( a_side == 0 ) {
         bc_type = m_bc_type[TOROIDAL_LOWER];
      }
      else if ( a_side == 1 ) {
         bc_type = m_bc_type[TOROIDAL_UPPER];
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::getBCType(): Invalid side argument");
      }
   }
#endif

   return bc_type;
}

std::string
LogRectEllipticOpBC::getBCSubType(const int  a_block_number,
                                  const int  a_dir,
                                  const int  a_side ) const
{
   std::string bc_subtype;

   if ( a_dir == RADIAL_DIR ) {
      if ( a_side == 0 ) {
         bc_subtype = m_bc_subtype[RADIAL_LOWER];
      }
      else if ( a_side == 1 ) {
         bc_subtype = m_bc_subtype[RADIAL_UPPER];
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::getBCSubType(): Invalid side argument");
      }
   }
   else if ( a_dir == POLOIDAL_DIR ) {
      if ( a_side == 0 ) {
         bc_subtype = m_bc_subtype[POLOIDAL_LOWER];
      }
      else if ( a_side == 1 ) {
         bc_subtype = m_bc_subtype[POLOIDAL_UPPER];
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::getBCSubType(): Invalid side argument");
      }
   }
#if CFG_DIM==3
   else if ( a_dir == TOROIDAL_DIR ) {
      if ( a_side == 0 ) {
         bc_subtype = m_bc_subtype[TOROIDAL_LOWER];
      }
      else if ( a_side == 1 ) {
         bc_subtype = m_bc_subtype[TOROIDAL_UPPER];
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::getBCSubType(): Invalid side argument");
      }
   }
#endif
   else {
      MayDay::Error("LogRectEllipticOpBC::getBCSubType(): Invalid direction argument");
   }

   return bc_subtype;
}


void
LogRectEllipticOpBC::setBCValue( const int     a_block_number,
                                 const int     a_dir,
                                 const int     a_side,
                                 const double  a_value )
{
   if ( a_dir == RADIAL_DIR ) {
      if ( a_side == 0 ) {
         m_bc_value[RADIAL_LOWER] = a_value;
      }
      else if ( a_side == 1 ) {
         m_bc_value[RADIAL_UPPER] = a_value;
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::setBCValue(): Invalid side argument");
      }
   }
   else if ( a_dir == POLOIDAL_DIR ) {
      if ( a_side == 0 ) {
         m_bc_value[POLOIDAL_LOWER] = a_value;
      }
      else if ( a_side == 1 ) {
         m_bc_value[POLOIDAL_UPPER] = a_value;
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::setBCValue(): Invalid side argument");
      }
   }
#if CFG_DIM==3
   else if ( a_dir == TOROIDAL_DIR ) {
      if ( a_side == 0 ) {
         m_bc_value[TOROIDAL_LOWER] = a_value;
      }
      else if ( a_side == 1 ) {
         m_bc_value[TOROIDAL_UPPER] = a_value;
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::setBCValue(): Invalid side argument");
      }
   }
#endif
   else {
      MayDay::Error("LogRectEllipticOpBC::setBCValue(): Invalid direction argument");
   }
}



double
LogRectEllipticOpBC::getBCValue( const int  a_block_number,
                                 const int  a_dir,
                                 const int  a_side ) const
{
   double value = BASEFAB_REAL_SETVAL;

   if ( a_dir == RADIAL_DIR ) {
      if ( a_side == 0 ) {
         value = m_bc_value[RADIAL_LOWER];
      }
      else if ( a_side == 1 ) {
         value = m_bc_value[RADIAL_UPPER];
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::getBCValue(): Invalid side argument");
      }
   }
   else if ( a_dir == POLOIDAL_DIR ) {
      if ( a_side == 0 ) {
         value = m_bc_value[POLOIDAL_LOWER];
      }
      else if ( a_side == 1 ) {
         value = m_bc_value[POLOIDAL_UPPER];
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::getBCValue(): Invalid side argument");
      }
   }
#if CFG_DIM==3
   else if ( a_dir == TOROIDAL_DIR ) {
      if ( a_side == 0 ) {
         value = m_bc_value[TOROIDAL_LOWER];
      }
      else if ( a_side == 1 ) {
         value = m_bc_value[TOROIDAL_UPPER];
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::getBCValue(): Invalid side argument");
      }
   }
#endif
   else {
      MayDay::Error("LogRectEllipticOpBC::getBCValue(): Invalid direction argument");
   }

   return value;
}



void
LogRectEllipticOpBC::setBCFunction( const int                           a_block_number,
                                    const int                           a_dir,
                                    const int                           a_side,
                                    const RefCountedPtr<GridFunction>&  a_function )
{
   if ( a_dir == RADIAL_DIR ) {
      if ( a_side == 0 ) {
         m_bc_function[RADIAL_LOWER] = a_function;
      }
      else if ( a_side == 1 ) {
         m_bc_function[RADIAL_UPPER] = a_function;
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::setBCFunction(): Invalid side argument");
      }
   }
   else if ( a_dir == POLOIDAL_DIR ) {
      if ( a_side == 0 ) {
         m_bc_function[POLOIDAL_LOWER] = a_function;
      }
      else if ( a_side == 1 ) {
         m_bc_function[POLOIDAL_UPPER] = a_function;
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::setBCFunction(): Invalid side argument");
      }
   }
#if CFG_DIM==3
   else if ( a_dir == TOROIDAL_DIR ) {
      if ( a_side == 0 ) {
         m_bc_function[TOROIDAL_LOWER] = a_function;
      }
      else if ( a_side == 1 ) {
         m_bc_function[TOROIDAL_UPPER] = a_function;
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::setBCFunction(): Invalid side argument");
      }
   }
#endif
   else {
      MayDay::Error("LogRectEllipticOpBC::setBCFunction(): Invalid direction argument");
   }
}



RefCountedPtr<GridFunction>
LogRectEllipticOpBC::getBCFunction( const int  a_block_number,
                                    const int  a_dir,
                                    const int  a_side ) const
{
   RefCountedPtr<GridFunction> function;

   if (getBCSubType(a_block_number, a_dir, a_side) == "coupled" ) {
      function = getBlockBCData(a_block_number, a_dir, a_side);
   }
   
   else {
      if ( a_dir == RADIAL_DIR ) {
         if ( a_side == 0 ) {
            function = m_bc_function[RADIAL_LOWER];
         }
         else if ( a_side == 1 ) {
            function = m_bc_function[RADIAL_UPPER];
         }
         else {
            MayDay::Error("LogRectEllipticOpBC::getBCFunction(): Invalid side argument");
         }
      }
      else if ( a_dir == POLOIDAL_DIR ) {
         if ( a_side == 0 ) {
            function = m_bc_function[POLOIDAL_LOWER];
         }
         else if ( a_side == 1 ) {
            function = m_bc_function[POLOIDAL_UPPER];
         }
         else {
            MayDay::Error("LogRectEllipticOpBC::getBCFunction(): Invalid side argument");
         }
      }
#if CFG_DIM==3
      else if ( a_dir == TOROIDAL_DIR ) {
         if ( a_side == 0 ) {
            function = m_bc_function[TOROIDAL_LOWER];
         }
         else if ( a_side == 1 ) {
            function = m_bc_function[TOROIDAL_UPPER];
         }
         else {
            MayDay::Error("LogRectEllipticOpBC::getBCFunction(): Invalid side argument");
         }
      }
#endif
      else {
         MayDay::Error("LogRectEllipticOpBC::getBCFunction(): Invalid direction argument");
      }
   }
   return function;
}



void
LogRectEllipticOpBC::apply( const MultiBlockLevelGeom&  a_geom,
                            const Box&                  a_coord_sys_box,
                            const double&               a_time,
                            const int                   a_dir,
                            const int                   a_side,
                            FArrayBox&                  a_phi ) const
{
   RefCountedPtr<GridFunction> func;
   double value;

   if ( a_dir == RADIAL_DIR ) {
      if ( a_side == 0 ) {
         func = m_bc_function[RADIAL_LOWER];
         value = m_bc_value[RADIAL_LOWER];
      }
      else if ( a_side == 1 ) {
         func = m_bc_function[RADIAL_UPPER];
         value = m_bc_value[RADIAL_UPPER];
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::apply(): Invalid side argument");
      }
   }
   else if ( a_dir == POLOIDAL_DIR ) {
      if ( a_side == 0 ) {
         func = m_bc_function[POLOIDAL_LOWER];
         value = m_bc_value[POLOIDAL_LOWER];
      }
      else if ( a_side == 1 ) {
         func = m_bc_function[POLOIDAL_UPPER];
         value = m_bc_value[POLOIDAL_UPPER];
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::apply(): Invalid side argument");
      }
   }
#if CFG_DIM==3
   else if ( a_dir == TOROIDAL_DIR ) {
      if ( a_side == 0 ) {
         func = m_bc_function[TOROIDAL_LOWER];
         value = m_bc_value[TOROIDAL_LOWER];
      }
      else if ( a_side == 1 ) {
         func = m_bc_function[TOROIDAL_UPPER];
         value = m_bc_value[TOROIDAL_UPPER];
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::apply(): Invalid side argument");
      }
   }
#endif
   else {
      MayDay::Error("LogRectEllipticOpBC::apply(): Invalid direction argument");
   }

   if ( !func.isNull() ) {
      // Get the block id
      const MultiBlockCoordSys& coord_sys( *(a_geom.coordSysPtr()) );
      const int block_number( coord_sys.whichBlock( a_coord_sys_box ) );
      FArrayBox dummy;
      func->assign(a_phi, a_geom, dummy, dummy, block_number, a_time, false);
   }
   else {
      a_phi.setVal(value);
   }
}



void LogRectEllipticOpBC::printParameters() const
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "LogRectEllipticOpBC ================================" << std::endl;
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
void LogRectEllipticOpBC::parseParameters( ParmParse&   a_pp,
                                           const bool*  a_is_periodic )
{
   GridFunctionLibrary* library = GridFunctionLibrary::getInstance();
   
   for (int i=0; i<m_bc_type.size(); i++) {
      std::string prefix( a_pp.prefix() );
      prefix += "." + m_bdry_name[i];
      ParmParse fpp( prefix.c_str() );
      std::string bc_type;
      fpp.query( "type", bc_type );

      m_bc_type[i] = UNDEFINED;

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
      else if (bc_type.length() == 0) {
         switch(i) 
            {
            case RADIAL_LOWER:
               if ( !a_is_periodic[RADIAL_DIR] ) {
                  MayDay::Error("LogRectEllipticOpBC::parseParameters(): No boundary condition specified on lower radial boundary!");
               }
               break;
            case RADIAL_UPPER:
               if ( !a_is_periodic[RADIAL_DIR] ) {
                  MayDay::Error("LogRectEllipticOpBC::parseParameters(): No boundary condition specified on upper radial boundary!");
               }
               break;
            case POLOIDAL_LOWER:
               if ( !a_is_periodic[POLOIDAL_DIR] ) {
                  MayDay::Error("LogRectEllipticOpBC::parseParameters(): No boundary condition specified on lower poloidal boundary!");
               }
               break;
            case POLOIDAL_UPPER:
               if ( !a_is_periodic[POLOIDAL_DIR] ) {
                  MayDay::Error("LogRectEllipticOpBC::parseParameters(): No boundary condition specified on upper poloidal boundary!");
               }
               break;
#if CFG_DIM==3
            case TOROIDAL_LOWER:
               if ( !a_is_periodic[TOROIDAL_DIR] ) {
                  MayDay::Error("LogRectEllipticOpBC::parseParameters(): No boundary condition specified on lower toroidal boundary!");
               }
               break;
            case TOROIDAL_UPPER:
               if ( !a_is_periodic[TOROIDAL_DIR] ) {
                  MayDay::Error("LogRectEllipticOpBC::parseParameters(): No boundary condition specified on upper toroidal boundary!");
               }
               break;
#endif
            }
      }
      else {
         MayDay::Error("LogRectEllipticOpBC::parseParameter(): Unknown potential bc type");
      }

      // If a boundary condition was set, then make sure it wasn't at a periodic boundary
      if ( m_bc_type[i] != UNDEFINED ) {

         switch(i) 
            {
            case RADIAL_LOWER:
               if ( a_is_periodic[RADIAL_DIR] ) {
                  MayDay::Error("LogRectEllipticOpBC::parseParameters(): Boundary condition specified on periodic lower radial boundary!");
               }
               break;
            case RADIAL_UPPER:
               if ( a_is_periodic[RADIAL_DIR] ) {
                  MayDay::Error("LogRectEllipticOpBC::parseParameters(): Boundary condition specified on periodic upper radial boundary!");
               }
               break;
            case POLOIDAL_LOWER:
               if ( a_is_periodic[POLOIDAL_DIR] ) {
                  MayDay::Error("LogRectEllipticOpBC::parseParameters(): Boundary condition specified on periodic lower poloidal boundary!");
               }
               break;
            case POLOIDAL_UPPER:
               if ( a_is_periodic[POLOIDAL_DIR] ) {
                  MayDay::Error("LogRectEllipticOpBC::parseParameters(): Boundary condition specified on periodic upper poloidal boundary!");
               }
               break;
#if CFG_DIM==3
            case TOROIDAL_LOWER:
               if ( a_is_periodic[TOROIDAL_DIR] ) {
                  MayDay::Error("LogRectEllipticOpBC::parseParameters(): Boundary condition specified on periodic lower toroidal boundary!");
               }
               break;
            case TOROIDAL_UPPER:
               if ( a_is_periodic[TOROIDAL_DIR] ) {
                  MayDay::Error("LogRectEllipticOpBC::parseParameters(): Boundary condition specified on periodic upper toroidal boundary!");
               }
               break;
#endif
            }
      }


      fpp.query( "subtype", m_bc_subtype[i] );
      
      bool value_specified = fpp.contains("value");
      
      if (fpp.contains("value")) {
         fpp.query("value", m_bc_value[i]);
      }
      else {
         m_bc_value[i] = 0.;
      }
      
      bool function_specified = fpp.contains("function");
      
      if (fpp.contains("function")) {
         std::string function_name;
         fpp.query( "function", function_name );
         m_bc_function[i] = library->find( function_name );
      }
      
      if (value_specified && function_specified) {
         MayDay::Error("LogRectEllipticOpBC::parseParameters(): Please specify either a value or a function, but not both");
      }
   }

   if (m_verbosity) {
      printParameters();
   }
}


RefCountedPtr<EllipticOpBC>
LogRectEllipticOpBC::clone( const bool a_extrapolated ) const
{
   // This sets m_bdry_name too
   RefCountedPtr<LogRectEllipticOpBC> result
      = RefCountedPtr<LogRectEllipticOpBC>(new LogRectEllipticOpBC(m_num_blocks));

   // Copy the data members set by the full constructor (i.e., the one with the ParmParse
   // argument).  
   result->m_name = m_name;
   result->m_verbosity = m_verbosity;

   // Copy the rest of the data owned by the EllipticOpBC base class
   result->copyBaseData(*this);

   if ( a_extrapolated ) {
      result->setExtrapolatedType();
   }

   return result;
}

#include "NamespaceFooter.H"
