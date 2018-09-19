#ifndef _CELLVAR_H_
#define _CELLVAR_H_

#include "CFGVar.H"

#include "NamespaceHeader.H"


/// Cell-averaged or cell-centered variable
/**
 * Cell-averaged or cell-centered variable
 */
class CellVar
   : public CFGVar
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
         : CFGVar(pp_prefix, name, geometry)
      {
         cell_data().define(geometry.gridsFull(), depth, ghosts);
      }
};


#include "NamespaceFooter.H"

#endif
