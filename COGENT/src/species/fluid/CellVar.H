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
         : CFGVars(pp_prefix, name, geometry)
      {
         addCellVar(name, depth, ghosts);
      }
};


#include "NamespaceFooter.H"

#endif
