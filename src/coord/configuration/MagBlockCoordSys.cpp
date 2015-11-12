#include "MagBlockCoordSys.H"
#include "MagBlockCoordSysF_F.H"

#include "Directions.H"
#include "FourthOrderUtil.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"

#include "NamespaceHeader.H"

inline Real computePointwiseJ( const Vector<Real>& Jmtrx )
{
   CH_assert(Jmtrx.size()==SpaceDim*SpaceDim);
   Real detJ = 0.0;
   if ( SpaceDim == 1 )
   {
      detJ = Jmtrx[0];
   }
   else if ( SpaceDim == 2 )
   {
      detJ = Jmtrx[0] * Jmtrx[3] - Jmtrx[1] * Jmtrx[2];
   }
   else if ( SpaceDim == 3 )
   {
      detJ = Jmtrx[0] * ( Jmtrx[4] * Jmtrx[8] - Jmtrx[5] * Jmtrx[7] )
           + Jmtrx[1] * ( Jmtrx[5] * Jmtrx[6] - Jmtrx[3] * Jmtrx[8] )
           + Jmtrx[2] * ( Jmtrx[3] * Jmtrx[7] - Jmtrx[4] * Jmtrx[6] );
   }
   else
   {
      MayDay::Error("NewFourthOrderCoordSys::pointwiseJ not defined for SpaceDim > 3");
   }
   return detJ;
}



MagBlockCoordSys::MagBlockCoordSys( ParmParse& a_parm_parse )
   : m_axisymmetric(true),
     m_verbose(false)
{
   if (a_parm_parse.contains("verbose")) {
      a_parm_parse.get("verbose", m_verbose);
   }

   if (a_parm_parse.contains("axisymmetric")) {
      a_parm_parse.get("axisymmetric", m_axisymmetric);
   }

   if (a_parm_parse.contains("inner_radial_bdry")) {
      a_parm_parse.get("inner_radial_bdry", m_inner_radial_boundary);
   }

   if (a_parm_parse.contains("outer_radial_bdry")) {
      a_parm_parse.get("outer_radial_bdry", m_outer_radial_boundary);
   }
}



void MagBlockCoordSys::define( const ProblemDomain& a_domain,
                               const RealVect& a_cellSpacing )
{
   m_domain = a_domain;
   m_dx = a_cellSpacing;

   m_mapped_cell_volume = 1.;
   for (int dir=0; dir<SpaceDim; ++dir) {
      m_mapped_cell_volume *= m_dx[dir];
      m_mapped_face_area[dir] = 1.;
      for (int tdir=0; tdir<SpaceDim; ++tdir) {
         if (tdir != dir) m_mapped_face_area[dir] *= m_dx[tdir];
      }
   }
}



void
MagBlockCoordSys::cellVol( FArrayBox&     a_vol,
                           const FluxBox& a_N,
                           const Box&     a_box ) const
{
  /*
    This code routine was borrowed from NewFourthOrderCoordSys::cellVol(),
    modified to incorporate the axisymmetric option.
  */

  // The flux function F needs to have two more ghost cells than cellVolumes
  // because we need to take derivatives of <F>, and because we
  // need to take derivatives in order to compute 4th-order avg of F
  Box fbox(a_box);
  fbox.grow(2);
  FluxBox F(fbox, SpaceDim);

  F.setVal(0.0);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      RealVect offset = 0.5*m_dx;
      offset[dir] = 0.0;

      FArrayBox& thisFdir = F[dir];
      BoxIterator bit(thisFdir.box());
      // this is going to be slow, but we can
      // eventually move this into fortran
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect mappedLoc = m_dx*iv;
          mappedLoc += offset;
          RealVect realLoc = realCoord(mappedLoc);

          double R = realLoc[0];
          thisFdir(iv,0) = R;

          if ( m_axisymmetric ) {
            thisFdir(iv,0) *= 0.5 * R;
          }

          thisFdir(iv,1) = realLoc[1];
        }
    } // end loop over directions


  // convert point values to 4th-order face averages
  fourthOrderAverageFace(F);

  Interval divInterval = a_vol.interval();
  computeDivergence(a_vol, F, a_N, a_box, divInterval);

  a_vol *= 0.5;
}



void MagBlockCoordSys::getPointwiseMajorRadius( FluxBox& a_R ) const
{
  for (int dir=0; dir<SpaceDim; ++dir) {
    FArrayBox& this_R = a_R[dir];
    const Box& box = this_R.box();

    FArrayBox xi(box,SpaceDim);
    getFaceCenteredMappedCoords(dir, xi);

    BoxIterator bit(box);
    for (bit.begin();bit.ok();++bit) {
      const IntVect& iv = bit();

      RealVect this_xi;
      for (int comp=0; comp<SpaceDim; ++comp) {
        this_xi[comp] = xi(iv,comp);
      }

      this_R(iv,0) = majorRadiusConstMinorRadAdj(this_xi);
    }
  }
}



void MagBlockCoordSys::getCellCentereddXdxi( FArrayBox& a_dXdxi ) const
{
   CH_assert(a_dXdxi.nComp() == SpaceDim*SpaceDim);
   const Box& box = a_dXdxi.box();

   FArrayBox xi(box,SpaceDim);
   getCellCenteredMappedCoords(xi);

   dXdXi(a_dXdxi, xi, 0, 0, 0, box);
   dXdXi(a_dXdxi, xi, 1, 0, 1, box);
   dXdXi(a_dXdxi, xi, 2, 1, 0, box);
   dXdXi(a_dXdxi, xi, 3, 1, 1, box);
}



void MagBlockCoordSys::getFaceCentereddXdxi( FluxBox& a_dXdxi ) const
{
   CH_assert(a_dXdxi.nComp() == SpaceDim*SpaceDim);

   for (int dir=0; dir<SpaceDim; ++dir) {
      FArrayBox& this_dXdxi = a_dXdxi[dir];
      const Box& box = this_dXdxi.box();

      FArrayBox xi(box,SpaceDim);
      getFaceCenteredMappedCoords(dir, xi);

      dXdXi(this_dXdxi, xi, 0, 0, 0, box);
      dXdXi(this_dXdxi, xi, 1, 0, 1, box);
      dXdXi(this_dXdxi, xi, 2, 1, 0, box);
      dXdXi(this_dXdxi, xi, 3, 1, 1, box);
   }
}



void MagBlockCoordSys::getPointwiseN( FArrayBox& a_N ) const
{
   const Box& box = a_N.box();

   FArrayBox xi(box,SpaceDim);
   getCellCenteredMappedCoords(xi);

   dXdXi(a_N, xi, 0, 1, 1, box);
   dXdXi(a_N, xi, 1, 1, 0, box);
   dXdXi(a_N, xi, 2, 0, 1, box);
   dXdXi(a_N, xi, 3, 0, 0, box);
   a_N.negate(1,2);

   if ( m_axisymmetric ) {
      BoxIterator bit(box);
      for (bit.begin();bit.ok();++bit) {
         const IntVect& iv = bit();

         RealVect this_xi;
         this_xi[RADIAL_DIR] = xi(iv,RADIAL_DIR);
         this_xi[POLOIDAL_DIR] = xi(iv,POLOIDAL_DIR);

         double TwoPiRmaj = 2. * Pi * majorRadius(this_xi);

         for (int comp=0; comp<4; ++comp) {
            a_N(iv,comp) *= TwoPiRmaj;
         }
      }
   }
}



void MagBlockCoordSys::getPointwiseN( FluxBox& a_N ) const
{
   for (int dir=0; dir<SpaceDim; ++dir) {
      FArrayBox& this_N = a_N[dir];
      const Box& box = this_N.box();

      FArrayBox xi(box,SpaceDim);
      getFaceCenteredMappedCoords(dir, xi);

      dXdXi(this_N, xi, 0, 1, 1, box);
      dXdXi(this_N, xi, 1, 1, 0, box);
      dXdXi(this_N, xi, 2, 0, 1, box);
      dXdXi(this_N, xi, 3, 0, 0, box);
      this_N.negate(1,2);

      if ( m_axisymmetric ) {
         BoxIterator bit(box);
         for (bit.begin();bit.ok();++bit) {
            const IntVect& iv = bit();

            RealVect this_xi;
            this_xi[RADIAL_DIR] = xi(iv,RADIAL_DIR);
            this_xi[POLOIDAL_DIR] = xi(iv,POLOIDAL_DIR);
            
            double TwoPiRmaj = 2. * Pi * majorRadius(this_xi);

            for (int comp=0; comp<4; ++comp) {
               this_N(iv,comp) *= TwoPiRmaj;
            }
         }
      }
   }
}



void MagBlockCoordSys::getPointwiseNJInverse( FArrayBox& a_NJInverse ) const
{
   const Box& box = a_NJInverse.box();

   FArrayBox xi(box,SpaceDim);
   getCellCenteredMappedCoords(xi);

   dXdXi(a_NJInverse, xi, 0, 1, 1, box);
   dXdXi(a_NJInverse, xi, 1, 1, 0, box);
   dXdXi(a_NJInverse, xi, 2, 0, 1, box);
   dXdXi(a_NJInverse, xi, 3, 0, 0, box);
   a_NJInverse.negate(1,2);

   BoxIterator bit(box);
   for (bit.begin();bit.ok();++bit) {
      const IntVect& iv = bit();

      double J = a_NJInverse(iv,0)*a_NJInverse(iv,3) - a_NJInverse(iv,1)*a_NJInverse(iv,2);

      for (int comp=0; comp<4; ++comp) {
         a_NJInverse(iv,comp) /= J;
      }
   }
}



void MagBlockCoordSys::getPointwiseNJInverse( FluxBox& a_NJInverse ) const
{
   for (int dir=0; dir<SpaceDim; ++dir) {
      FArrayBox& NJInverse_dir = a_NJInverse[dir];
      const Box& box = NJInverse_dir.box();

      FArrayBox xi(box,SpaceDim);
      getFaceCenteredMappedCoords(dir, xi);

      dXdXi(NJInverse_dir, xi, 0, 1, 1, box);
      dXdXi(NJInverse_dir, xi, 1, 1, 0, box);
      dXdXi(NJInverse_dir, xi, 2, 0, 1, box);
      dXdXi(NJInverse_dir, xi, 3, 0, 0, box);
      NJInverse_dir.negate(1,2);

      BoxIterator bit(box);
      for (bit.begin();bit.ok();++bit) {
         const IntVect& iv = bit();

         double J = NJInverse_dir(iv,0)*NJInverse_dir(iv,3) - NJInverse_dir(iv,1)*NJInverse_dir(iv,2);

         for (int comp=0; comp<4; ++comp) {
            NJInverse_dir(iv,comp) /= J;
         }
      }
   }
}



void MagBlockCoordSys::pointwiseJ(FArrayBox& a_J,
                                  const FArrayBox& a_Xi,
                                  const Box& a_box) const
{
  NewFourthOrderCoordSys::pointwiseJ(a_J, a_Xi, a_box);

  if ( m_axisymmetric ) {

    BoxIterator bit(a_box);
    for (bit.begin();bit.ok();++bit) {
      const IntVect& iv = bit();

      RealVect this_xi;
      this_xi[RADIAL_DIR] = a_Xi(iv,RADIAL_DIR);
      this_xi[POLOIDAL_DIR] = a_Xi(iv,POLOIDAL_DIR);

      a_J(iv) *= 2. * Pi * majorRadius(this_xi);
    }
  }
}



void MagBlockCoordSys::getTwoPirRmaj( FArrayBox& a_TwoPirRmaj ) const
{
   const Box& box = a_TwoPirRmaj.box();

   FArrayBox xi(box,SpaceDim);
   getCellCenteredMappedCoords(xi);

   FArrayBox dXdXi_fab(box,SpaceDim*SpaceDim);
   dXdXi(dXdXi_fab, xi, 0, 1, 1, box);
   dXdXi(dXdXi_fab, xi, 1, 1, 0, box);
   dXdXi(dXdXi_fab, xi, 2, 0, 1, box);
   dXdXi(dXdXi_fab, xi, 3, 0, 0, box);
   dXdXi_fab.negate(1,2);

   BoxIterator bit(box);
   for (bit.begin();bit.ok();++bit) {
      const IntVect& iv = bit();

      // N is the inverse of the 2x2 matrix dXdXi times its determinant J
      double r = dXdXi_fab(iv,0)*dXdXi_fab(iv,3) - dXdXi_fab(iv,1)*dXdXi_fab(iv,2);

      a_TwoPirRmaj(iv,0) = r;

      if ( m_axisymmetric ) {
        RealVect this_xi;
        this_xi[RADIAL_DIR] = xi(iv,RADIAL_DIR);
        this_xi[POLOIDAL_DIR] = xi(iv,POLOIDAL_DIR);

        a_TwoPirRmaj(iv,0) *= 2. * Pi * majorRadius(this_xi);
      }
   }
}



void
MagBlockCoordSys::getCellCenteredMappedCoords( FArrayBox& a_xi ) const
{
  FORT_GET_CC_MAPPED_COORDS(CHF_BOX(a_xi.box()),
                                   CHF_CONST_REALVECT(m_dx),
                                   CHF_FRA(a_xi));
}



void
MagBlockCoordSys::getFaceCenteredMappedCoords( const int a_dir,
                                               FArrayBox&  a_xi ) const
{
  FORT_GET_FC_MAPPED_COORDS(CHF_BOX(a_xi.box()),
                            CHF_CONST_INT(a_dir),
                            CHF_CONST_REALVECT(m_dx),
                            CHF_FRA(a_xi));
}


void
MagBlockCoordSys::getNodeCenteredRealCoords( FArrayBox&  a_X ) const
{
   for (BoxIterator bit(a_X.box()); bit.ok(); ++bit) {
      IntVect iv = bit();

      RealVect mapped_loc = m_dx*iv;
      RealVect real_loc = realCoord(mapped_loc);
      for (int n=0; n<SpaceDim; ++n) {
         a_X(iv,n) = real_loc[n];
      }
   }
}


void
MagBlockCoordSys::getCellCenteredRealCoords( FArrayBox& a_x ) const
{
   const Box& box( a_x.box() );
   FArrayBox xi_array( box, a_x.nComp() );
   getCellCenteredMappedCoords( xi_array );

   for (BoxIterator bit( box ); bit.ok(); ++bit) {
      IntVect iv( bit() );
      RealVect xi;
      for (int dir(0); dir<SpaceDim; ++dir) {
         xi[dir] = xi_array( iv, dir );
      }
      RealVect real_loc( realCoord( xi ) );
      for (int dir(0); dir<SpaceDim; ++dir) {
         a_x( iv, dir ) = real_loc[dir];
      }
   }
}


void
MagBlockCoordSys::getFaceCenteredRealCoords( const int a_dir, FArrayBox& a_x ) const
{
   const Box& box( a_x.box() );
   FArrayBox xi_array( box, SpaceDim );
   getFaceCenteredMappedCoords( a_dir, xi_array );

   for (BoxIterator bit( box ); bit.ok(); ++bit) {
      IntVect iv( bit() );
      RealVect xi;
      for (int dir(0); dir<SpaceDim; ++dir) {
         xi[dir] = xi_array( iv, dir );
      }
      RealVect real_loc( realCoord( xi ) );
      for (int dir(0); dir<SpaceDim; ++dir) {
         a_x( iv, dir ) = real_loc[dir];
      }
   }
}



void
MagBlockCoordSys::getCellCenteredFluxCoords( FArrayBox& a_xi ) const
{
   const Box& box( a_xi.box() );
   FArrayBox xi_array( box, a_xi.nComp() );
   MagBlockCoordSys::getCellCenteredMappedCoords( xi_array );

   RealVect offset( RealVect::Zero );
   offset[RADIAL_DIR] = getInnerFluxLabel();

   for (BoxIterator bit( box ); bit.ok(); ++bit) {
      IntVect iv( bit() );
      for (int dir(0); dir<SpaceDim; ++dir) {
         a_xi( iv, dir ) = xi_array( iv, dir ) + offset[dir];
      }
   }
}



void
MagBlockCoordSys::computeFluxSurfaceUnitTangent( FArrayBox& a_data ) const
{
   CH_assert(a_data.nComp() == SpaceDim);

   const Box& box = a_data.box();
   FArrayBox dXdxi(box, SpaceDim*SpaceDim);
   getCellCentereddXdxi(dXdxi);

   FORT_UNIT_FS_TANGENT(CHF_BOX(box),
                        CHF_CONST_FRA(dXdxi),
                        CHF_FRA(a_data));
}



void
MagBlockCoordSys::computeFluxSurfaceUnitTangent( FluxBox& a_data ) const
{
   CH_assert(a_data.nComp() == SpaceDim);

   FluxBox dXdxi(a_data.box(), SpaceDim*SpaceDim);
   getFaceCentereddXdxi(dXdxi);

   for (int dir=0; dir<SpaceDim; ++dir) {
      FArrayBox& this_data = a_data[dir];
      FORT_UNIT_FS_TANGENT(CHF_BOX(this_data.box()),
                           CHF_CONST_FRA(dXdxi[dir]),
                           CHF_FRA(this_data));
   }
}



void
MagBlockCoordSys::projectOntoFluxSurface( FArrayBox& a_data ) const
{
   FArrayBox fs_unit_tangent(a_data.box(), a_data.nComp());
   computeFluxSurfaceUnitTangent(fs_unit_tangent);

   FORT_MAG_BLOCK_PROJECT(CHF_BOX(a_data.box()),
                          CHF_CONST_FRA(fs_unit_tangent),
                          CHF_FRA(a_data));
}



void
MagBlockCoordSys::projectOntoFluxSurface( FluxBox& a_data ) const
{
   for (int dir=0; dir<SpaceDim; ++dir) {
      projectOntoFluxSurface(a_data[dir]);
   }
}



void
MagBlockCoordSys::computeFluxSurfaceUnitNormal( FArrayBox& a_data ) const
{
   CH_assert(a_data.nComp() == SpaceDim);

   const Box& box = a_data.box();
   FArrayBox dXdxi(box, SpaceDim*SpaceDim);
   getCellCentereddXdxi(dXdxi);

   FORT_UNIT_FS_NORMAL(CHF_BOX(box),
                       CHF_CONST_FRA(dXdxi),
                       CHF_FRA(a_data));
}



void
MagBlockCoordSys::computeFluxSurfaceUnitNormal( FluxBox& a_data ) const
{
   CH_assert(a_data.nComp() == SpaceDim);

   FluxBox dXdxi(a_data.box(), SpaceDim*SpaceDim);
   getFaceCentereddXdxi(dXdxi);

   for (int dir=0; dir<SpaceDim; ++dir) {
      FArrayBox& this_data = a_data[dir];
      FORT_UNIT_FS_NORMAL(CHF_BOX(this_data.box()),
                          CHF_CONST_FRA(dXdxi[dir]),
                          CHF_FRA(this_data));
   }
}



void
MagBlockCoordSys::projectOntoFluxSurfaceNormal( FArrayBox& a_data ) const
{
   FArrayBox fs_unit_normal(a_data.box(), a_data.nComp());
   computeFluxSurfaceUnitNormal(fs_unit_normal);

   FORT_MAG_BLOCK_PROJECT(CHF_BOX(a_data.box()),
                          CHF_CONST_FRA(fs_unit_normal),
                          CHF_FRA(a_data));
}



void
MagBlockCoordSys::projectOntoFluxSurfaceNormal( FluxBox& a_data ) const
{
   for (int dir=0; dir<SpaceDim; ++dir) {
      projectOntoFluxSurfaceNormal(a_data[dir]);
   }
}

void
MagBlockCoordSys::computePsiThetaProjections( FArrayBox& a_data ) const
{
   FArrayBox fs_unit_normal(a_data.box(), a_data.nComp());
   computeFluxSurfaceUnitNormal(fs_unit_normal);

   FArrayBox fs_unit_tangent(a_data.box(), a_data.nComp());
   computeFluxSurfaceUnitTangent(fs_unit_tangent);

   FORT_MAG_BLOCK_PSITHETA_PROJECTIONS(CHF_BOX(a_data.box()),
                                       CHF_CONST_FRA(fs_unit_normal),
                                       CHF_CONST_FRA(fs_unit_tangent),
                                       CHF_FRA(a_data));

}

void
MagBlockCoordSys::computePsiThetaProjections( FluxBox& a_data ) const
{
   for (int dir=0; dir<SpaceDim; ++dir) {
      computePsiThetaProjections(a_data[dir]);
   }
}


void
MagBlockCoordSys::computeUnitRadial( FArrayBox& a_data ) const
{
   CH_assert(a_data.nComp() == SpaceDim);

   const Box& box = a_data.box();
   FArrayBox dXdxi(box, SpaceDim*SpaceDim);
   getCellCentereddXdxi(dXdxi);

   FORT_UNIT_RADIAL(CHF_BOX(box),
                    CHF_CONST_FRA(dXdxi),
                    CHF_FRA(a_data));
}



void
MagBlockCoordSys::computeUnitRadial( FluxBox& a_data ) const
{
   CH_assert(a_data.nComp() == SpaceDim);

   FluxBox dXdxi(a_data.box(), SpaceDim*SpaceDim);
   getFaceCentereddXdxi(dXdxi);

   for (int dir=0; dir<SpaceDim; ++dir) {
      FArrayBox& this_data = a_data[dir];
      FORT_UNIT_RADIAL(CHF_BOX(this_data.box()),
                       CHF_CONST_FRA(dXdxi[dir]),
                       CHF_FRA(this_data));
   }
}



void
MagBlockCoordSys::projectOntoRadial( FArrayBox& a_data ) const
{
   FArrayBox unit_radial(a_data.box(), a_data.nComp());
   computeUnitRadial(unit_radial);

   FORT_MAG_BLOCK_PROJECT(CHF_BOX(a_data.box()),
                          CHF_CONST_FRA(unit_radial),
                          CHF_FRA(a_data));
}



void
MagBlockCoordSys::projectOntoRadial( FluxBox& a_data ) const
{
   for (int dir=0; dir<SpaceDim; ++dir) {
      projectOntoRadial(a_data[dir]);
   }
}



void
MagBlockCoordSys::computeGradFFactor( FArrayBox& a_data ) const
{
   CH_assert(a_data.nComp() == 1);

   const Box& box = a_data.box();
   FArrayBox dXdxi(box, SpaceDim*SpaceDim);
   getCellCentereddXdxi(dXdxi);

   FORT_GRADF_FACTOR(CHF_BOX(box),
                     CHF_CONST_FRA(dXdxi),
                     CHF_FRA1(a_data,0));
}



void
MagBlockCoordSys::computeGradFFactor( FluxBox& a_data ) const
{
   CH_assert(a_data.nComp() == 1);

   FluxBox dXdxi(a_data.box(), SpaceDim*SpaceDim);
   getFaceCentereddXdxi(dXdxi);

   for (int dir=0; dir<SpaceDim; ++dir) {
      FArrayBox& this_data = a_data[dir];
      FORT_GRADF_FACTOR(CHF_BOX(this_data.box()),
                        CHF_CONST_FRA(dXdxi[dir]),
                        CHF_FRA1(this_data,0));
   }
}



#if 1
void
MagBlockCoordSys::incrementFaceMetricWithEdgeTerm(FArrayBox& a_N,
                                                  int a_faceDir,
                                                  int a_edgeDir,
                                                  const Box& a_box,
                                                  const Side::LoHiSide& a_side) const
{
  // make sure we have enough components
  CH_assert(a_N.nComp() >= SpaceDim);

  // this is where we do the integration along the edges
  Real mult = 1.0;
  if (a_side == Side::Lo) mult = -1.0;

  // this is the offset dx*(i,j,k) for the faces in a_faceMetrics
  RealVect faceOffset = m_dx;
  faceOffset *= 0.5;
  faceOffset[a_faceDir] = 0.0;

  RealVect edgeOffset = faceOffset;
  edgeOffset[a_edgeDir] += mult*0.5*m_dx[a_edgeDir];

  const Vector<QuadratureElement > quadPts = m_quadraturePtr->coefficients(a_faceDir,
                                                                           a_edgeDir);
  Real weightMult = m_quadraturePtr->weightMult(m_dx, a_faceDir, a_edgeDir);

  // this is gonna be slow, but we can hopefully make it faster at
  // some point down the road
  BoxIterator bit(a_box);

  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect faceCenter = m_dx*iv + faceOffset;

      RealVect edgeCenter = m_dx*iv + edgeOffset;

      for (int sDir = 0; sDir<SpaceDim; sDir++)
        {
          Real edgeVal = 0.0;

          // since 2D is just nodes, do this separately
          if (SpaceDim == 2)
            {
              if ( m_axisymmetric ) {
                if (sDir == 0) {
                  double Z = realCoord(edgeCenter)[1];
                  edgeVal = (1 - 2*a_faceDir) * Z * 2. * Pi;
                }
                else {
                  double R = majorRadius(edgeCenter);
                  edgeVal = (2*a_faceDir - 1) * Pi * R * R;
                }
              }
              else {
                edgeVal = NewFourthOrderCoordSys::getN(edgeCenter, sDir, a_faceDir, a_edgeDir);
              }
            }
          else if (SpaceDim == 3)
            {
              // integrate along tangentDir using quadrature
              for (int i=0; i<quadPts.size(); i++)
                {
                  RealVect mappedLoc(edgeCenter);
                  mappedLoc += (quadPts[i].location)*m_dx/2.0;
                  // convert from mapped->real space
                  //RealVect realLoc = realCoord(mappedLoc);
                  Real Nvalue = NewFourthOrderCoordSys::getN(mappedLoc, sDir, a_faceDir, a_edgeDir);
                  edgeVal += Nvalue*quadPts[i].weight*weightMult;
                }

            }
          else
            {
              MayDay::Error("NewFourthOrderCoordSys::faceMetrics not defined for SpaceDim > 3");
            }

          // since this is for normal components, faceDir = d
          int component = getNcomponent(sDir,a_faceDir);
          a_N(iv, component) += mult*edgeVal;

        } // end loop over s directions
    } // end loop over faces in faceMetrics

}
#endif


#if 1
void
MagBlockCoordSys::computeTransverseFaceMetric(FArrayBox& a_faceMetrics,
                                              const Box& a_box,
                                              int a_faceDir,
                                              int a_dDir) const
{
  // first, compute face area
  Real faceArea = m_dx[0];
  for (int dir=1; dir<SpaceDim; dir++)
    {
      faceArea *= m_dx[dir];
    }
  faceArea /= m_dx[a_faceDir];


  // temp storage for nodal values
  Box nodeBox(a_box);
  Box LapBox(nodeBox);
  // switch box to node centering and grow by one in the
  // transverse directions (Laplacian doesn't need the extra cells)
  for (int dir=0; dir<SpaceDim; dir++)
    {
      if (dir != a_faceDir)
        {
          nodeBox.surroundingNodes(dir);
          nodeBox.grow(dir, 1);
          LapBox.surroundingNodes(dir);
        }
    }

  FArrayBox nodeN(nodeBox, SpaceDim);


  // compute N at nodes
  BoxIterator bit(nodeBox);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect mappedLoc = m_dx*iv;
      RealVect realLoc = realCoord(mappedLoc);

      // now, compute derivatives for Jacobian matrix.  Note that
      // we will not need anything in the dDir column. so don't
      // bother to compute them. All other entries wind up being used 2x
      // (which is why we're precomputing them)
      // do this with an eye toward using the pointwiseJ function
      // NOTE -- since eqn(6) in Phil's notes is actually for N^T instead
      // of N, use transpose of Jacobian matrix (dXdXi(row,col) instead
      // of dxdXi(col,row) as one might expect.  This results in computing
      // (N^T)^T, which is, of course, N.
      Vector<Real> Jmatrix(SpaceDim*SpaceDim, 0.0);

      for (int col=0; col<SpaceDim; col++)
        {
          if (col != a_dDir)
            {
              for (int row=0; row<SpaceDim; row++)
                {
                  Jmatrix[row + SpaceDim*col] = dXdXi(mappedLoc, row, col);
                }
            }
        } // end loop over columns for jacobian matrix computation


      for (int sDir=0; sDir<SpaceDim; sDir++)
        {
          // make use of the computePointwiseJ function to take the
          // determinant
          // first, copy the derivative values
          Vector<Real> jmtrx = Jmatrix;
          // now, replace the sDir'th row with the dth unit vector
          for (int j=0; j<SpaceDim; j++)
            {
              if (j == a_dDir)
                {
                  jmtrx[sDir + SpaceDim*j] = 1.0;
                }
              else
                {
                  jmtrx[sDir + SpaceDim*j] = 0.0;
                }
            }
          Real determinant = computePointwiseJ(jmtrx);
          nodeN(iv,sDir) = determinant;
          if (m_axisymmetric) nodeN(iv,sDir) *= 2. * Pi * majorRadius(mappedLoc);
        }
    } // end loop over cells


  // compute node-centered tangential Laplacian in fortran
  FArrayBox lapN(LapBox,SpaceDim);
  lapN.setVal(0.0);

  for (int dir=0; dir<SpaceDim; dir++)
    {
      if (dir != a_faceDir)
        {
          // simple undivided difference approximation to Laplacian
          Real factor = 1.0;
          FORT_INCREMENTLAPLACIAN2(CHF_FRA(lapN),
                                  CHF_CONST_FRA(nodeN),
                                  CHF_BOX(LapBox),
                                  CHF_CONST_INT(dir),
                                  CHF_CONST_REAL(factor));
        }
    }

  // now average to cell centers from nodes.

  // this box represents which nodes will be averaged to
  // the face center
  Box averageBox(IntVect::Zero, IntVect::Unit);
  averageBox.setBig(a_faceDir, 0);
  BoxIterator averageIterator(averageBox);

  Real averageFactor = D_TERM6(1.0, *0.5, *0.5, *0.5, *0.5, *0.5);

  // this is another candidate for fortran...
  BoxIterator faceBit(a_faceMetrics.box());
  for (faceBit.begin(); faceBit.ok(); ++faceBit)
    {
      IntVect faceIV = faceBit();
      // start by just computing simple average to face
      // do this by looping through the averageBox
      for (averageIterator.begin(); averageIterator.ok(); ++averageIterator)
        {
          IntVect nodeLoc = faceIV + averageIterator();
          for (int sDir=0; sDir<SpaceDim; sDir++)
            {
              int metricComp = getNcomponent(sDir, a_dDir);
              a_faceMetrics(faceIV,metricComp) += averageFactor*nodeN(nodeLoc,sDir);
              // increment with Laplacian term as well
              // this is a fourth-order point value
              //a_faceMetrics(faceIV,metricComp) -= 0.125*averageFactor*lapN(nodeLoc,sDir);
              // this is a fourth-order face average
              a_faceMetrics(faceIV,metricComp) -= averageFactor*lapN(nodeLoc,sDir)/12.0;

            }
        }

    }
  // multiply by face area to get the integral over the face
  for (int sDir=0; sDir<SpaceDim; sDir++)
    {
      int metricComp = getNcomponent(sDir, a_dDir);
      //
      a_faceMetrics.mult(faceArea, metricComp, 1);
    }
}
#endif


#include "NamespaceFooter.H"
