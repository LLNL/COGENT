#ifndef _CELLVAR_H_
#define _CELLVAR_H_

#include "CFGVars.H"

#include "NamespaceHeader.H"


/// Cell-averaged or cell-centered variable
/**
 * Cell-averaged or cell-centered variable
 */
class CellVar
   : public CFGVars
{
   public:
   
      /// Constructor.
      /**
       */
      CellVar( const string&       pp_prefix,
               const std::string&  name,
               const MagGeom&      geometry,
               const int           depth,
               const IntVect       ghosts )
         : CFGVars(pp_prefix, name, geometry),
           m_depth(depth),
           m_ghosts(ghosts)
      {
         addCellVar(name, depth, ghosts);
      }

      virtual RefCountedPtr<CFGVars> clone( const IntVect&  a_ghost_vect = IntVect::Zero,
                                            const bool      a_copy_data = true ) const
      {
         CFGVars* field_ptr = new CellVar( m_pp_prefix, m_name, m_geometry, m_depth, m_ghosts );

         if (a_copy_data) {
            field_ptr->copy(*this);
         }

         return RefCountedPtr<CFGVars>(field_ptr);
      }

      protected:

      int m_depth;
      IntVect m_ghosts;
};


#include "NamespaceFooter.H"

#endif
