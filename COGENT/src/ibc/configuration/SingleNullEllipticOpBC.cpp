#include "SingleNullEllipticOpBC.H"
#include "SingleNullBlockCoordSys.H"
#include "GridFunctionLibrary.H"
#include "MagGeom.H"
#include "Directions.H"

#include "NamespaceHeader.H"



SingleNullEllipticOpBC::SingleNullEllipticOpBC()
   : EllipticOpBC(NUM_BOUNDARIES)
{
   setNames();
}



SingleNullEllipticOpBC::SingleNullEllipticOpBC( const std::string&  a_name,
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
SingleNullEllipticOpBC::setNames()
{
   m_bdry_name[RADIAL_CORE] = "radial_core";
   m_bdry_name[RADIAL_SOL] = "radial_sol";
   m_bdry_name[RADIAL_PF] = "radial_pf";
   m_bdry_name[POLOIDAL_INNER_DIV] = "poloidal_inner_div";
   m_bdry_name[POLOIDAL_OUTER_DIV] = "poloidal_outer_div";
}



void
SingleNullEllipticOpBC::setBCType( const int  a_block_number,
                                   const int  a_dir,
                                   const int  a_side,
                                   const int  a_type )
{
   CH_assert(a_side == 0 || a_side == 1);

   switch( a_block_number )
      {
      case MCORE:
      case LCORE:
      case RCORE:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            m_bc_type[RADIAL_CORE] = a_type;
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::setBCType(): Invalid argument");
         }
         break;
      case MCSOL:
      case LCSOL:
      case RCSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            m_bc_type[RADIAL_SOL] = a_type;
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::setBCType(): Invalid argument");
         }
         break;
      case LSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            m_bc_type[RADIAL_SOL] = a_type;
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 1) {
            m_bc_type[POLOIDAL_INNER_DIV] = a_type;
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::setBCType(): Invalid argument");
         }
         break;
      case RSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            m_bc_type[RADIAL_SOL] = a_type;
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 0) {
            m_bc_type[POLOIDAL_OUTER_DIV] = a_type;
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::setBCType(): Invalid argument");
         }
         break;
      case LPF:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            m_bc_type[RADIAL_PF] = a_type;
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 1) {
            m_bc_type[POLOIDAL_INNER_DIV] = a_type;
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::setBCType(): Invalid argument");
         }
         break;
      case RPF:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            m_bc_type[RADIAL_PF] = a_type;
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 0) {
            m_bc_type[POLOIDAL_OUTER_DIV] = a_type;
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::setBCType(): Invalid argument");
         }
         break;
      default:
         MayDay::Error("SingleNullEllipticOpBC::setBCType(): Unrecognized block number");
      }
}


int
SingleNullEllipticOpBC::getBCType( const int  a_block_number,
                                   const int  a_dir,
                                   const int  a_side ) const
{
   CH_assert(a_side == 0 || a_side == 1);
   int bc_type = UNDEFINED;

   switch( a_block_number )
      {
      case MCORE:
      case LCORE:
      case RCORE:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            bc_type = m_bc_type[RADIAL_CORE];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::getBCType(): Invalid argument");
         }
         break;
      case MCSOL:
      case LCSOL:
      case RCSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            bc_type = m_bc_type[RADIAL_SOL];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::getBCType(): Invalid argument");
         }
         break;
      case LSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            bc_type = m_bc_type[RADIAL_SOL];
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 1) {
            bc_type = m_bc_type[POLOIDAL_INNER_DIV];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::getBCType(): Invalid argument");
         }
         break;
      case RSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            bc_type = m_bc_type[RADIAL_SOL];
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 0) {
            bc_type = m_bc_type[POLOIDAL_OUTER_DIV];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::getBCType(): Invalid argument");
         }
         break;
      case LPF:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            bc_type = m_bc_type[RADIAL_PF];
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 1) {
            bc_type = m_bc_type[POLOIDAL_INNER_DIV];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::getBCType(): Invalid argument");
         }
         break;
      case RPF:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            bc_type = m_bc_type[RADIAL_PF];
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 0) {
            bc_type = m_bc_type[POLOIDAL_OUTER_DIV];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::getBCType(): Invalid argument");
         }
         break;
      default:
         MayDay::Error("SingleNullEllipticOpBC::getBCType(): Unrecognized block number");
      }

   return bc_type;
}



void
SingleNullEllipticOpBC::setBCValue( const int     a_block_number,
                                    const int     a_dir,
                                    const int     a_side,
                                    const double  a_value )
{
   CH_assert(a_side == 0 || a_side == 1);

   switch( a_block_number )
      {
      case MCORE:
      case LCORE:
      case RCORE:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            m_bc_value[RADIAL_CORE] = a_value;
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::setBCValue(): Invalid argument");
         }
         break;
      case MCSOL:
      case LCSOL:
      case RCSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            m_bc_value[RADIAL_SOL] = a_value;
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::setBCValue(): Invalid argument");
         }
         break;
      case LSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            m_bc_value[RADIAL_SOL] = a_value;
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 1) {
            m_bc_value[POLOIDAL_INNER_DIV] = a_value;
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::setBCValue(): Invalid argument");
         }
         break;
      case RSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            m_bc_value[RADIAL_SOL] = a_value;
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 0) {
            m_bc_value[POLOIDAL_OUTER_DIV] = a_value;
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::setBCValue(): Invalid argument");
         }
         break;
      case LPF:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            m_bc_value[RADIAL_PF] = a_value;
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 1) {
            m_bc_value[POLOIDAL_INNER_DIV] = a_value;
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::setBCValue(): Invalid argument");
         }
         break;
      case RPF:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            m_bc_value[RADIAL_PF] = a_value;
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 0) {
            m_bc_value[POLOIDAL_OUTER_DIV] = a_value;
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::setBCValue(): Invalid argument");
         }
         break;
      default:
         MayDay::Error("SingleNullEllipticOpBC::setBCValue(): Unrecognized block number");
      }
}



double
SingleNullEllipticOpBC::getBCValue( const int  a_block_number,
                                    const int  a_dir,
                                    const int  a_side ) const
{
   CH_assert(a_side == 0 || a_side == 1);
   double bc_value = BASEFAB_REAL_SETVAL;

   switch( a_block_number )
      {
      case MCORE:
      case LCORE:
      case RCORE:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            bc_value = m_bc_value[RADIAL_CORE];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::getBCValue(): Invalid argument");
         }
         break;
      case MCSOL:
      case LCSOL:
      case RCSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            bc_value = m_bc_value[RADIAL_SOL];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::getBCValue(): Invalid argument");
         }
         break;
      case LSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            bc_value = m_bc_value[RADIAL_SOL];
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 1) {
            bc_value = m_bc_value[POLOIDAL_INNER_DIV];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::getBCValue(): Invalid argument");
         }
         break;
      case RSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            bc_value = m_bc_value[RADIAL_SOL];
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 0) {
            bc_value = m_bc_value[POLOIDAL_OUTER_DIV];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::getBCValue(): Invalid argument");
         }
         break;
      case LPF:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            bc_value = m_bc_value[RADIAL_PF];
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 1) {
            bc_value = m_bc_value[POLOIDAL_INNER_DIV];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::getBCValue(): Invalid argument");
         }
         break;
      case RPF:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            bc_value = m_bc_value[RADIAL_PF];
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 0) {
            bc_value = m_bc_value[POLOIDAL_OUTER_DIV];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::getBCValue(): Invalid argument");
         }
         break;
      default:
         MayDay::Error("SingleNullEllipticOpBC::getBCValue(): Unrecognized block number");
      }

   return bc_value;
}



void
SingleNullEllipticOpBC::setBCFunction( const int                           a_block_number,
                                       const int                           a_dir,
                                       const int                           a_side,
                                       const RefCountedPtr<GridFunction>&  a_function )
{
   CH_assert(a_side == 0 || a_side == 1);

   switch( a_block_number )
      {
      case MCORE:
      case LCORE:
      case RCORE:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            m_bc_function[RADIAL_CORE] = a_function;
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::setBCFunction(): Invalid argument");
         }
         break;
      case MCSOL:
      case LCSOL:
      case RCSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            m_bc_function[RADIAL_SOL] = a_function;
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::setBCFunction(): Invalid argument");
         }
         break;
      case LSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            m_bc_function[RADIAL_SOL] = a_function;
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 1) {
            m_bc_function[POLOIDAL_INNER_DIV] = a_function;
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::setBCFunction(): Invalid argument");
         }
         break;
      case RSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            m_bc_function[RADIAL_SOL] = a_function;
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 0) {
            m_bc_function[POLOIDAL_OUTER_DIV] = a_function;
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::setBCFunction(): Invalid argument");
         }
         break;
      case LPF:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            m_bc_function[RADIAL_PF] = a_function;
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 1) {
            m_bc_function[POLOIDAL_INNER_DIV] = a_function;
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::setBCFunction(): Invalid argument");
         }
         break;
      case RPF:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            m_bc_function[RADIAL_PF] = a_function;
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 0) {
            m_bc_function[POLOIDAL_OUTER_DIV] = a_function;
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::setBCFunction(): Invalid argument");
         }
         break;
      default:
         MayDay::Error("SingleNullEllipticOpBC::setBCFunction(): Unrecognized block number");
      }
}



RefCountedPtr<GridFunction>
SingleNullEllipticOpBC::getBCFunction( const int  a_block_number,
                                       const int  a_dir,
                                       const int  a_side ) const
{
   CH_assert(a_side == 0 || a_side == 1);
   RefCountedPtr<GridFunction> function;

   switch( a_block_number )
      {
      case MCORE:
      case LCORE:
      case RCORE:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            function = m_bc_function[RADIAL_CORE];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::getBCFunction(): Invalid argument");
         }
         break;
      case MCSOL:
      case LCSOL:
      case RCSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            function = m_bc_function[RADIAL_SOL];
         }
         break;
      case LSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            function = m_bc_function[RADIAL_SOL];
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 1) {
            function = m_bc_function[POLOIDAL_INNER_DIV];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::getBCFunction(): Invalid argument");
         }
         break;
      case RSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            function = m_bc_function[RADIAL_SOL];
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 0) {
            function = m_bc_function[POLOIDAL_OUTER_DIV];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::getBCFunction(): Invalid argument");
         }
         break;
      case LPF:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            function = m_bc_function[RADIAL_PF];
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 1) {
            function = m_bc_function[POLOIDAL_INNER_DIV];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::getBCFunction(): Invalid argument");
         }
         break;
      case RPF:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            function = m_bc_function[RADIAL_PF];
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 0) {
            function = m_bc_function[POLOIDAL_OUTER_DIV];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::getBCFunction(): Invalid argument");
         }
         break;
      default:
         MayDay::Error("SingleNullEllipticOpBC::getBCFunction(): Unrecognized block number");
      }

   return function;
}



void
SingleNullEllipticOpBC::apply( const MultiBlockLevelGeom&  a_geom,
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
      case MCORE:
      case LCORE:
      case RCORE:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            function = m_bc_function[RADIAL_CORE];
            value = m_bc_value[RADIAL_CORE];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::apply(): Invalid argument");
         }
         break;
      case MCSOL:
      case LCSOL:
      case RCSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            function = m_bc_function[RADIAL_SOL];
            value = m_bc_value[RADIAL_SOL];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::apply(): Invalid argument");
         }
         break;
      case LSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            function = m_bc_function[RADIAL_SOL];
            value = m_bc_value[RADIAL_SOL];
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 1) {
            function = m_bc_function[POLOIDAL_INNER_DIV];
            value = m_bc_value[POLOIDAL_INNER_DIV];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::apply(): Invalid argument");
         }
         break;
      case RSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            function = m_bc_function[RADIAL_SOL];
            value = m_bc_value[RADIAL_SOL];
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 0) {
            function = m_bc_function[POLOIDAL_OUTER_DIV];
            value = m_bc_value[POLOIDAL_OUTER_DIV];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::apply(): Invalid argument");
         }
         break;
      case LPF:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            function = m_bc_function[RADIAL_PF];
            value = m_bc_value[RADIAL_PF];
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 1) {
            function = m_bc_function[POLOIDAL_INNER_DIV];
            value = m_bc_value[POLOIDAL_INNER_DIV];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::apply(): Invalid argument");
         }
         break;
      case RPF:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            function = m_bc_function[RADIAL_PF];
            value = m_bc_value[RADIAL_PF];
         }
         else if (a_dir == POLOIDAL_DIR && a_side == 0) {
            function = m_bc_function[POLOIDAL_OUTER_DIV];
            value = m_bc_value[POLOIDAL_OUTER_DIV];
         }
         else {
            MayDay::Error("SingleNullEllipticOpBC::apply(): Invalid argument");
         }
         break;
      default:
         MayDay::Error("SingleNullEllipticOpBC::apply(): Unrecognized block number");
      }

   if ( !function.isNull() ) {
      function->assign(a_phi, a_geom, a_coord_sys_box, a_time, false);
   }
   else {
      a_phi.setVal(value);
   }
}



void SingleNullEllipticOpBC::printParameters() const
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "SingleNullEllipticOpBC ================================" << std::endl;
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
void SingleNullEllipticOpBC::parseParameters( ParmParse& a_pp )
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
         MayDay::Error("SingleNullEllipticOpBC::parseParameter(): Unrecognized potential bc type");
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
         MayDay::Error("SingleNullEllipticOpBC::parseParameters(): Please specify either a value or a function, but not both");
      }
   }

   if (m_verbosity) {
      printParameters();
   }
}



#include "NamespaceFooter.H"
