#ifndef _ADVECTIVE_COORD_SYS_
#define _ADVECTIVE_COORD_SYS_

#include "NamespaceHeader.H"

/// Coordinate system interface for a slice of phase space.
/**
 * A reification of the \verb+NewFourthOrderCoordSys+ for
 * a planar subset of phase space at a particular magnetic moment $\mu$.
 */
class AdvectiveCoordSys : public NewFourthOrderCoordSys
{
      /// Constructor with definition.
      /**
       * Constructor with definition.
       *
       * @param[in] phase_space_coord_sys phase space coordinate system.
       * @param[in] dx phase space mesh size.
       */
      AdvectiveCoordSys(
         const PhaseSpaceCoordSys& phase_space_coord_sys,
         const Real               dx );
};

#include "NamespaceFooter.H"

#endif
