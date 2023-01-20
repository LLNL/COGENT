#include "Potentials.H"
#include "SpaceUtils.H.multidim"

#include "NamespaceHeader.H"


Potentials::Potentials(const string&     a_pp_prefix,
                       const string&     a_name,
                       const MagGeom&    a_geometry,
                       const IntVect&    a_ghost_vect,
                       const bool        a_evolve_A_parallel)
   : CFGVars(a_pp_prefix, a_name, a_geometry),
     m_evolve_A_parallel(a_evolve_A_parallel)
{
   ParmParse pp(a_pp_prefix.c_str());

   addCellVar("potential", 1, a_ghost_vect);
   
   if(m_evolve_A_parallel) {
     addCellVar("A_parallel", 1, a_ghost_vect);
   }
}


RefCountedPtr<CFGVars>
Potentials::clone( const IntVect&  a_ghost_vect,
                   const bool      a_copy_data ) const
{
   CFGVars* field_ptr = new Potentials( m_pp_prefix, m_name, m_geometry, a_ghost_vect, m_evolve_A_parallel );

   if ( a_copy_data ) {
      field_ptr->copy(*this);
   }

   return RefCountedPtr<CFGVars>(field_ptr);
}


#include "NamespaceFooter.H"
