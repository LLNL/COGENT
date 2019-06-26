#include "LogRectEllipticOpBC.H"
#include "Directions.H"

#include "NamespaceHeader.H"


LogRectEllipticOpBC::LogRectEllipticOpBC()
   : EllipticOpBC(NUM_BOUNDARIES)
{
   setNames();
}


LogRectEllipticOpBC::LogRectEllipticOpBC( const std::string&  a_name,
                                          ParmParse&          a_pp,
                                          const int&          a_verbosity )
   : EllipticOpBC(NUM_BOUNDARIES),
     m_name(a_name),
     m_verbosity(a_verbosity)
{
   setNames();
   parseParameters( a_pp );
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
   else {
      MayDay::Error("LogRectEllipticOpBC::getBCType(): Invalid direction argument");
   }

   return bc_type;
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
void LogRectEllipticOpBC::parseParameters( ParmParse& a_pp )
{
   for (int i=0; i<m_bc_type.size(); i++) {
      std::string prefix( a_pp.prefix() );
      prefix += "." + m_bdry_name[i];
      ParmParse fpp( prefix.c_str() );
      std::string bc_type;
      fpp.query( "type", bc_type );

      if (bc_type.length()==0){ //make radial_inner or radial_lower compatible
        if( m_bdry_name[i].compare("radial_outer")==0 ){
            if (procID()==0) cout<<"LogRectEllipticOpBC: m_bdry_name["<<i<<"] = \"radial_outer\" not found, trying \"radial_upper\""<<endl;
            prefix.replace(prefix.end()-12,prefix.end(),"radial_upper");
            ParmParse fpp( prefix.c_str() );
            fpp.query( "type", bc_type );
        }
        else if ( m_bdry_name[i].compare("radial_upper")==0 ){
            if (procID()==0) cout<<"LogRectEllipticOpBC: m_bdry_name["<<i<<"] = \"radial_upper\" not found, trying \"radial_outer\""<<endl;
            prefix.replace(prefix.end()-12,prefix.end(),"radial_outer");
            ParmParse fpp( prefix.c_str() );
            fpp.query( "type", bc_type);
        }
        else if ( m_bdry_name[i].compare("radial_inner")==0 ){
            if (procID()==0) cout<<"LogRectEllipticOpBC: m_bdry_name["<<i<<"] = \"radial_inner\" not found, trying \"radial_lower\""<<endl;
            prefix.replace(prefix.end()-12,prefix.end(),"radial_lower");
            ParmParse fpp( prefix.c_str() );
            fpp.query( "type", bc_type);
        }
        else if ( m_bdry_name[i].compare("radial_lower")==0 ){
            if (procID()==0) cout<<"LogRectEllipticOpBC: m_bdry_name["<<i<<"] = \"radial_lower\" not found, trying \"radial_inner\""<<endl;
            prefix.replace(prefix.end()-12,prefix.end(),"radial_inner");
            ParmParse fpp( prefix.c_str() );
            fpp.query( "type", bc_type);
        } 
      }


      if (bc_type == "dirichlet") {
         m_bc_type[i] = DIRICHLET;
      }

      else if (bc_type == "neumann") {
         m_bc_type[i] = NEUMANN;
      }
      
      //Default option
      else if (bc_type.length() == 0) {
         if ( procID()==0 ) {
            std::cout << "LogRectEllipticOpBC: No potential BC is specified at the " << m_bdry_name[i]
            <<" boundary: a zero Dirichlet BC will be used if the direction is not periodic" << std::endl;
         }
         m_bc_type[i] = DIRICHLET;
      }

      else {
         MayDay::Error("LogRectEllipticOpBC::parseParameter(): Illegal potential bc type");
      }

      if (fpp.contains("value")) {
         fpp.query("value", m_bc_value[i]);
      }
      else {
         m_bc_value[i] = 0.;
      }

   }

   if (m_verbosity) {
      printParameters();
   }
}



#include "NamespaceFooter.H"
