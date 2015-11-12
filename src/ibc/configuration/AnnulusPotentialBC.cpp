#include "AnnulusPotentialBC.H"
#include "Directions.H"

#include "NamespaceHeader.H"


AnnulusPotentialBC::AnnulusPotentialBC()
   : PotentialBC(NUM_BOUNDARIES)
{
   setNames();
}


AnnulusPotentialBC::AnnulusPotentialBC( const std::string& a_name,
                                        ParmParse& a_pp,
                                        const int& a_verbosity )
   : PotentialBC(NUM_BOUNDARIES),
     m_name(a_name),
     m_verbosity(a_verbosity)
{
   setNames();
   parseParameters( a_pp );
}



void
AnnulusPotentialBC::setNames()
{
   m_bdry_name[RADIAL_INNER] = "radial_inner";
   m_bdry_name[RADIAL_OUTER] = "radial_outer";
}



void
AnnulusPotentialBC::setBCType(const int a_block_number,
                              const int a_dir,
                              const int a_side,
                              const int a_type)
{
   if ( a_dir == RADIAL_DIR ) {
      if ( a_side == 0 ) {
         m_bc_type[RADIAL_INNER] = a_type;
      }
      else if ( a_side == 1 ) {
         m_bc_type[RADIAL_OUTER] = a_type;
      }
      else {
         MayDay::Error("AnnulusPotentialBC::setBCType(): Invalid side argument");
      }
   }
   else {
      MayDay::Error("AnnulusPotentialBC::setBCType(): Invalid direction argument");
   }
}



int
AnnulusPotentialBC::getBCType(const int a_block_number,
                              const int a_dir,
                              const int a_side) const
{
   int bc_type = UNDEFINED;

   if ( a_dir == RADIAL_DIR ) {
      if ( a_side == 0 ) {
         bc_type = m_bc_type[RADIAL_INNER];
      }
      else if ( a_side == 1 ) {
         bc_type = m_bc_type[RADIAL_OUTER];
      }
      else {
         MayDay::Error("AnnulusPotentialBC::getBCType(): Invalid side argument");
      }
   }
   else {
      MayDay::Error("AnnulusPotentialBC::getBCType(): Invalid direction argument");
   }

   return bc_type;
}



void
AnnulusPotentialBC::setBCValue(const int    a_block_number,
                               const int    a_dir,
                               const int    a_side,
                               const double a_value)
{
   if ( a_dir == RADIAL_DIR ) {
      if ( a_side == 0 ) {
         m_bc_value[RADIAL_INNER] = a_value;
      }
      else if ( a_side == 1 ) {
         m_bc_value[RADIAL_OUTER] = a_value;
      }
      else {
         MayDay::Error("AnnulusPotentialBC::setBCValue(): Invalid side argument");
      }
   }
   else {
      MayDay::Error("AnnulusPotentialBC::setBCValue(): Invalid direction argument");
   }
}



double
AnnulusPotentialBC::getBCValue(const int a_block_number,
                               const int a_dir,
                               const int a_side) const
{
   double value = BASEFAB_REAL_SETVAL;

   if ( a_dir == RADIAL_DIR ) {
      if ( a_side == 0 ) {
         value = m_bc_value[RADIAL_INNER];
      }
      else if ( a_side == 1 ) {
         value = m_bc_value[RADIAL_OUTER];
      }
      else {
         MayDay::Error("AnnulusPotentialBC::getBCValue(): Invalid side argument");
      }
   }
   else {
      MayDay::Error("AnnulusPotentialBC::getBCValue(): Invalid direction argument");
   }

   return value;
}



void
AnnulusPotentialBC::setBCFunction(const int                          a_block_number,
                                  const int                          a_dir,
                                  const int                          a_side,
                                  const RefCountedPtr<GridFunction>& a_function)
{
   if ( a_dir == RADIAL_DIR ) {
      if ( a_side == 0 ) {
         m_bc_function[RADIAL_INNER] = a_function;
      }
      else if ( a_side == 1 ) {
         m_bc_function[RADIAL_OUTER] = a_function;
      }
      else {
         MayDay::Error("AnnulusPotentialBC::setBCFunction(): Invalid side argument");
      }
   }
   else {
      MayDay::Error("AnnulusPotentialBC::setBCFunction(): Invalid direction argument");
   }
}



RefCountedPtr<GridFunction>
AnnulusPotentialBC::getBCFunction( const int                          a_block_number,
                                   const int                          a_dir,
                                   const int                          a_side ) const
{
   RefCountedPtr<GridFunction> function;

   if ( a_dir == RADIAL_DIR ) {
      if ( a_side == 0 ) {
         function = m_bc_function[RADIAL_INNER];
      }
      else if ( a_side == 1 ) {
         function = m_bc_function[RADIAL_OUTER];
      }
      else {
         MayDay::Error("AnnulusPotentialBC::getBCFunction(): Invalid side argument");
      }
   }
   else {
      MayDay::Error("AnnulusPotentialBC::getBCFunction(): Invalid direction argument");
   }

   return function;
}



void
AnnulusPotentialBC::apply( const MultiBlockLevelGeom& a_geom,
                           const Box&                 a_coord_sys_box,
                           const double&              a_time,
                           const int                  a_dir,
                           const int                  a_side,
                           FArrayBox&                 a_phi ) const
{
   RefCountedPtr<GridFunction> func;
   double value;

   if ( a_dir == RADIAL_DIR ) {
      if ( a_side == 0 ) {
         func = m_bc_function[RADIAL_INNER];
         value = m_bc_value[RADIAL_INNER];
      }
      else if ( a_side == 1 ) {
         func = m_bc_function[RADIAL_OUTER];
         value = m_bc_value[RADIAL_OUTER];
      }
      else {
         MayDay::Error("AnnulusPotentialBC::apply(): Invalid side argument");
      }
   }
   else {
      MayDay::Error("AnnulusPotentialBC::apply(): Invalid direction argument");
   }

   if ( !func.isNull() ) {
      func->assign(a_phi, a_geom, a_coord_sys_box, a_time, false);
   }
   else {
      a_phi.setVal(value);
   }
}



void AnnulusPotentialBC::printParameters() const
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "AnnulusPotentialBC ================================" << std::endl;
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
void AnnulusPotentialBC::parseParameters( ParmParse& a_pp )
{
   for (int i=0; i<m_bc_type.size(); i++) {
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
         MayDay::Error("AnnulusPotentialBC::parseParameter(): Illegal potential bc type");
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
