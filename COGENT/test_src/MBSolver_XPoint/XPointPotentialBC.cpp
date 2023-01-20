#include "XPointPotentialBC.H"

#include "NamespaceHeader.H"



XPointPotentialBC::XPointPotentialBC()
   : EllipticOpBC(NUM_BOUNDARIES, 1)
{
   setNames();
}



XPointPotentialBC::XPointPotentialBC( const std::string& a_name,
                                      ParmParse&         a_pp,
                                      const int&         a_verbosity )
   : EllipticOpBC(NUM_BOUNDARIES, 1),
     m_name(a_name),
     m_verbosity(a_verbosity)
{
   setNames();
   parseParameters( a_pp );
}



void
XPointPotentialBC::setNames()
{
   m_bdry_name[BLOCK_0] = "block_0";
   m_bdry_name[BLOCK_1] = "block_1";
   m_bdry_name[BLOCK_2] = "block_2";
   m_bdry_name[BLOCK_3] = "block_3";
   m_bdry_name[BLOCK_4] = "block_4";
   m_bdry_name[BLOCK_5] = "block_5";
   m_bdry_name[BLOCK_6] = "block_6";
   m_bdry_name[BLOCK_7] = "block_7";
}



void
XPointPotentialBC::setBCType( const int a_block_number,
                                       const int a_dir,
                                       const int a_side,
                                       const int a_type )
{
   CH_assert(a_side == 0 || a_side == 1);
   CH_assert(a_block_number>=0 && a_block_number<NUM_BOUNDARIES);

   m_bc_type[a_block_number] = a_type;
}


int
XPointPotentialBC::getBCType( const int a_block_number,
                                       const int a_dir,
                                       const int a_side ) const
{
   CH_assert(a_side == 0 || a_side == 1);
   int bc_type = UNDEFINED;

   bc_type = m_bc_type[a_block_number];

   return bc_type;
}



void
XPointPotentialBC::setBCValue( const int    a_block_number,
                                        const int    a_dir,
                                        const int    a_side,
                                        const double a_value )
{
   CH_assert(a_side == 0 || a_side == 1);

   m_bc_value[a_block_number] = a_value;

}



double
XPointPotentialBC::getBCValue( const int a_block_number,
                                        const int a_dir,
                                        const int a_side ) const
{
   CH_assert(a_side == 0 || a_side == 1);
   double bc_value = BASEFAB_REAL_SETVAL;

   bc_value = m_bc_value[a_block_number];

   return bc_value;
}



void
XPointPotentialBC::setBCFunction( const int                          a_block_number,
                                           const int                          a_dir,
                                           const int                          a_side,
                                           const RefCountedPtr<GridFunction>& a_function )
{
   CH_assert(a_side == 0 || a_side == 1);

   m_bc_function[a_block_number] = a_function;
}



RefCountedPtr<GridFunction>
XPointPotentialBC::getBCFunction( const int a_block_number,
                                           const int a_dir,
                                           const int a_side ) const
{
   CH_assert(a_side == 0 || a_side == 1);
   RefCountedPtr<GridFunction> function;

   function = m_bc_function[a_block_number];

   return function;
}



void
XPointPotentialBC::apply( const MultiBlockLevelGeom& a_geom,
                          const Box&                 a_coord_sys_box,
                          const double&              a_time,
                          const int                  a_dir,
                          const int                  a_side,
                          FArrayBox&                 a_phi ) const
{
   MayDay::Error("Fix XPointPotentialBC::apply");
#if 0
   RefCountedPtr<GridFunction> function;
   double value;

   const MultiBlockCoordSys* coord_sys = a_geom.coordSysPtr();
   int block_number = coord_sys->whichBlock(a_coord_sys_box);

   value = m_bc_value[block_number];

   if ( !function.isNull() ) {
      function->assign(a_phi, a_geom, a_coord_sys_box, a_time, false);
   }
   else {
      a_phi.setVal(value);
   }
#endif
}



void XPointPotentialBC::printParameters() const
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "XPointPotentialBC ================================" << std::endl;
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
void XPointPotentialBC::parseParameters( ParmParse& a_pp )
{
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
         MayDay::Error("XPointPotentialBC::parseParameter(): Unrecognized potential bc type");
      }

      fpp.query( "value", m_bc_value[i] );
   }

   if (m_verbosity) {
      printParameters();
   }
}



#include "NamespaceFooter.H"
