#include "MagBlockCoordSys.H"
#include "MagBlockCoordSysF_F.H"

#include "Directions.H"
#include "FourthOrderUtil.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "inspect.H"

#include "NamespaceHeader.H"

inline Real computeMatrixDeterminant( const Vector<Real>& Jmtrx )
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
     m_field_aligned(false),
     m_pointwise_metrics(false),
     m_verbose(false),
     m_subgrid_geometry(false)
{
   if (a_parm_parse.contains("verbose")) {
      a_parm_parse.get("verbose", m_verbose);
   }

   if (a_parm_parse.contains("axisymmetric")) {
      a_parm_parse.get("axisymmetric", m_axisymmetric);

      if (m_axisymmetric && CFG_DIM==3) {
         MayDay::Error("MagBlockCoordSys::MagBlockCoordSys(): Cannot specify axisymmetry in 3D");
      }
   }
   
   if (a_parm_parse.contains("field_aligned")) {
      a_parm_parse.get("field_aligned", m_field_aligned);
   }

   if (a_parm_parse.contains("pointwise_metrics")) {
      a_parm_parse.get("pointwise_metrics", m_pointwise_metrics);
   }

   if (a_parm_parse.contains("subgrid_geometry")) {
      a_parm_parse.get("subgrid_geometry", m_subgrid_geometry);
   }

  
   Vector<int> ref_level(SpaceDim,1);
   if (a_parm_parse.contains("subgrid_refinement")) {
      a_parm_parse.getarr("subgrid_refinement",ref_level,0,SpaceDim);
   }
   for (int dir=0; dir<SpaceDim; ++dir) {
      m_N_mr[dir] = ref_level[dir];
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

inline void MagBlockCoordSys::TRM1x2less3x4( FArrayBox&        a_out,
                                             const FArrayBox&  a_in,
                                             const Box&        a_box,
                                             int               a_destComp,
                                             int i1, int i2, int i3, int i4 ) const
{
   // Sets destcomp of a_out to a_in[i1]*a_in[i2] - a_in[i3]*a_in[i4] where the index denotes
   // the component of a_in (which has 9 components for 3D calculations)
   FArrayBox tmp(a_box,1);
   a_out.copy(a_in,i1,a_destComp);
   a_out.mult(a_in,i2,a_destComp);
   tmp.copy(a_in,i3,0);
   tmp.mult(a_in,i4,0);
   a_out.minus(tmp,0,a_destComp);
}


inline Real MagBlockCoordSys::TRM1x2less3x4( const Vector<Real>&  a_in,
                                             int i1, int i2, int i3, int i4 ) const
{
   // Retruns a_in[i1]*a_in[i2] - a_in[i3]*a_in[i4] where the index denotes
   // the component of a_in (which has 9 components for 3D calculations)

   return a_in[i1] * a_in[i2] - a_in[i3] * a_in[i4];
}


void
MagBlockCoordSys::compute_dXdXi (FArrayBox&        a_dXdXiArr,
                                 const FArrayBox&  a_xi,
                                 const Box&        a_box  ) const
{
       dXdXi(a_dXdXiArr, a_xi, 0, 0, 0, a_box);
       dXdXi(a_dXdXiArr, a_xi, 1, 0, 1, a_box);
   if ( SpaceDim == 2 ){
       dXdXi(a_dXdXiArr, a_xi, 2, 1, 0, a_box);
       dXdXi(a_dXdXiArr, a_xi, 3, 1, 1, a_box);
   }
   else if ( SpaceDim == 3 ){
       dXdXi(a_dXdXiArr, a_xi, 2, 0, 2, a_box);
       dXdXi(a_dXdXiArr, a_xi, 3, 1, 0, a_box);
       dXdXi(a_dXdXiArr, a_xi, 4, 1, 1, a_box);
       dXdXi(a_dXdXiArr, a_xi, 5, 1, 2, a_box);
       dXdXi(a_dXdXiArr, a_xi, 6, 2, 0, a_box);
       dXdXi(a_dXdXiArr, a_xi, 7, 2, 1, a_box);
       dXdXi(a_dXdXiArr, a_xi, 8, 2, 2, a_box);
   }
   // NOTE ABOVE WILL PRODUCE ERROR IF spacedim = 1
   else{
       MayDay::Error("NewFourthOrderCoordSys::dXdXi  defined only for SpaceDim = 2 or 3");
   }
}

Vector<Real>
MagBlockCoordSys::compute_dXdXi (const RealVect&  a_Xi) const
{
  Vector<Real> result(SpaceDim*SpaceDim);
  result[0] = dXdXi(a_Xi, 0, 0);
  result[1] = dXdXi(a_Xi, 0, 1);

  if ( SpaceDim == 2 ){
    result[2] = dXdXi(a_Xi, 1, 0);
    result[3] = dXdXi(a_Xi, 1, 1);
  }
  
  else if ( SpaceDim == 3 ){

    result[2] = dXdXi(a_Xi, 0, 2);
    result[3] = dXdXi(a_Xi, 1, 0);
    result[4] = dXdXi(a_Xi, 1, 1);
    result[5] = dXdXi(a_Xi, 1, 2);
    result[6] = dXdXi(a_Xi, 2, 0);
    result[7] = dXdXi(a_Xi, 2, 1);
    result[8] = dXdXi(a_Xi, 2, 2);
    
  }
  // NOTE ABOVE WILL PRODUCE ERROR IF spacedim = 1
  else{
    MayDay::Error("MagBlockCoordSys::dXdXi  defined only for SpaceDim = 2 or 3");
  }

  return result;
  
}


void 
MagBlockCoordSys::compute_dXdXiCofac ( FArrayBox& a_cofarr,
                                       const FArrayBox& a_xi,
                                       const Box& a_box  ) const
{
   CH_TIMERS("MagBlockCoordSys::compute_dXdXiCofac");
   CH_TIMER("compute_dXdXi", t_compute_dxdxi);

   // Compute cofactors of dXdXi with negation
   if ( SpaceDim == 2 ){
       dXdXi(a_cofarr, a_xi, 0, 1, 1, a_box);
       dXdXi(a_cofarr, a_xi, 1, 1, 0, a_box);
       dXdXi(a_cofarr, a_xi, 2, 0, 1, a_box);
       dXdXi(a_cofarr, a_xi, 3, 0, 0, a_box);
       a_cofarr.negate(1,2);
   }
   else if ( SpaceDim == 3 ){
       CH_START(t_compute_dxdxi);
       FArrayBox tmp(a_box,9);
       for (int icomp = 0; icomp < 9; ++icomp){
          int ind1 = icomp/3;
          int ind2 = icomp % 3;
          dXdXi(tmp, a_xi, icomp, ind1,ind2,a_box);
       }
       CH_STOP(t_compute_dxdxi);

       // Get components of a_cofarr
       TRM1x2less3x4(a_cofarr, tmp, a_box, 0, 4, 8, 7, 5);
       TRM1x2less3x4(a_cofarr, tmp, a_box, 1, 5, 6, 3, 8);
       TRM1x2less3x4(a_cofarr, tmp, a_box, 2, 3, 7, 4, 6);
       TRM1x2less3x4(a_cofarr, tmp, a_box, 3, 2, 7, 1, 8);
       TRM1x2less3x4(a_cofarr, tmp, a_box, 4, 0, 8, 2, 6);
       TRM1x2less3x4(a_cofarr, tmp, a_box, 5, 1, 6, 0, 7);
       TRM1x2less3x4(a_cofarr, tmp, a_box, 6, 1, 5, 2, 4);
       TRM1x2less3x4(a_cofarr, tmp, a_box, 7, 2, 3, 0, 5);
       TRM1x2less3x4(a_cofarr, tmp, a_box, 8, 0, 4, 1, 3);
   }
   // NOTE ABOVE WILL PRODUCE ERROR IF spacedim = 1
   else{
       MayDay::Error("NewFourthOrderCoordSys::dXdXi  defined only for SpaceDim = 2 or 3");
   }
}

Vector<Real>
MagBlockCoordSys::compute_dXdXiCofac ( const RealVect& a_Xi) const
{
   CH_TIMERS("MagBlockCoordSys::compute_dXdXiCofac_point");
   CH_TIMER("compute_dXdXi", t_compute_dxdxi_point);

   Vector<Real> cofarr(SpaceDim*SpaceDim,0);
   
   // Compute cofactors of dXdXi with negation
   if ( SpaceDim == 2 ){
   
      cofarr[0] = dXdXi(a_Xi, 1, 1);
      cofarr[1] = -dXdXi(a_Xi, 1, 0);
      cofarr[2] = -dXdXi(a_Xi, 0, 1);
      cofarr[3] = dXdXi(a_Xi, 0, 0);
   }
   else if ( SpaceDim == 3 ){
      Vector<Real> tmp(SpaceDim*SpaceDim,0);
      CH_START(t_compute_dxdxi_point);
      for (int icomp = 0; icomp < 9; ++icomp){
         int ind1 = icomp/3;
         int ind2 = icomp % 3;
         tmp[icomp] = dXdXi(a_Xi, ind1, ind2);
      }
      CH_STOP(t_compute_dxdxi_point);
      
      // Get components of a_cofarr
      cofarr[0] = TRM1x2less3x4(tmp, 4, 8, 7, 5);
      cofarr[1] = TRM1x2less3x4(tmp, 5, 6, 3, 8);
      cofarr[2] = TRM1x2less3x4(tmp, 3, 7, 4, 6);
      cofarr[3] = TRM1x2less3x4(tmp, 2, 7, 1, 8);
      cofarr[4] = TRM1x2less3x4(tmp, 0, 8, 2, 6);
      cofarr[5] = TRM1x2less3x4(tmp, 1, 6, 0, 7);
      cofarr[6] = TRM1x2less3x4(tmp, 1, 5, 2, 4);
      cofarr[7] = TRM1x2less3x4(tmp, 2, 3, 0, 5);
      cofarr[8] = TRM1x2less3x4(tmp, 0, 4, 1, 3);
      
   }
   // NOTE ABOVE WILL PRODUCE ERROR IF spacedim = 1
   else{
      MayDay::Error("NewFourthOrderCoordSys::dXdXi  defined only for SpaceDim = 2 or 3");
   }
   
   return cofarr;
}

void
MagBlockCoordSys::cellVol( FArrayBox&     a_vol,
                           const FluxBox& a_N,
                           const Box&     a_box ) const
{
   
   if (!m_subgrid_geometry && !m_pointwise_metrics) {
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
      for (int dir=0; dir<SpaceDim; dir++) {
         RealVect offset = 0.5*m_dx;
         offset[dir] = 0.0;

         FArrayBox& thisFdir = F[dir];
         BoxIterator bit(thisFdir.box());
         // this is going to be slow, but we can
         // eventually move this into fortran
         for (bit.begin(); bit.ok(); ++bit) {
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
#if CFG_DIM==3
            thisFdir(iv,2) = realLoc[2];
#endif
         }
      } // end loop over directions

      // convert point values to 4th-order face averages
      fourthOrderAverageFace(F);

      Interval divInterval = a_vol.interval();
      computeDivergence(a_vol, F, a_N, a_box, divInterval);

      a_vol /= CFG_DIM;
   }
   
   if (m_pointwise_metrics) {

      FArrayBox J(a_box,1);
      FArrayBox Xi(a_box,SpaceDim);
      getCellCenteredMappedCoords(Xi);
      pointwiseJ(J, Xi, a_box);
      J *= getMappedCellVolume();
      a_vol.copy(J);
   }
   
   if (m_subgrid_geometry) {

      /*this will overwrite the pointwise calculation
       if both flags (pointwise_metric and subgrid are ON)
       because the latter is more accurate and does
       not access data outside of the given cell
       */
      FArrayBox Xi(a_box,SpaceDim);
      getCellCenteredMappedCoords(Xi);
      cellVolSubGrid(a_vol, Xi);
   }
}

void
MagBlockCoordSys::cellVolSubGrid(FArrayBox&        a_vol,
                                 const FArrayBox&  a_Xi) const
{
   //check the validity of the subgrid parameters
   for (int dir=0; dir<SpaceDim; dir++) {
      if (m_N_mr[dir] < 2) {
         MayDay::Error("MagBlockCoordSys::cellVolSubGrid refinement must be not less than 2 in all dirs");
      }

      if (m_N_mr[dir] % 2 != 0) {
         MayDay::Error("MagBlockCoordSys::cellVolSubGrid refinement must be even in all dirs");
      }
   }
   
   const Box& box = a_Xi.box();
   BoxIterator bit(box);
   for (bit.begin();bit.ok();++bit) {
      const IntVect& iv = bit();
      RealVect this_xi;
      double integrated_volume = 0.0;

#if CH_SPACEDIM==2
      for (int i = -m_N_mr[0]/2; i < m_N_mr[0]/2; i++) {
         for (int j = -m_N_mr[1]/2; j < m_N_mr[1]/2; j++) {
            
            //get subgrid cell centers
            this_xi[RADIAL_DIR] = a_Xi(iv,RADIAL_DIR) + (i + 0.5) * m_dx[0]/m_N_mr[0];
            this_xi[POLOIDAL_DIR] = a_Xi(iv,POLOIDAL_DIR) + (j + 0.5) * m_dx[1]/m_N_mr[1];

            RealVect Xi00(this_xi[0]-0.5*m_dx[0]/m_N_mr[0], this_xi[1]-0.5*m_dx[1]/m_N_mr[1]);
            RealVect Xi10(this_xi[0]+0.5*m_dx[0]/m_N_mr[0], this_xi[1]-0.5*m_dx[1]/m_N_mr[1]);
            RealVect Xi01(this_xi[0]-0.5*m_dx[0]/m_N_mr[0], this_xi[1]+0.5*m_dx[1]/m_N_mr[1]);
            RealVect Xi11(this_xi[0]+0.5*m_dx[0]/m_N_mr[0], this_xi[1]+0.5*m_dx[1]/m_N_mr[1]);

            RealVect X00(realCoord(Xi00));
            RealVect X10(realCoord(Xi10));
            RealVect X01(realCoord(Xi01));
            RealVect X11(realCoord(Xi11));
            
            //compute cell face area (2D area) as the 2D area of two triangles (parallelogram shape is not assume)
            double area_lo = 0.5 * fabs(X00[0]*(X01[1]-X10[1])+X01[0]*(X10[1]-X00[1])+X10[0]*(X00[1]-X01[1]));
            double area_hi = 0.5 * fabs(X11[0]*(X10[1]-X01[1])+X10[0]*(X01[1]-X11[1])+X01[0]*(X11[1]-X10[1]));

            double fac = 1.0;
            if (m_axisymmetric) {
              fac = 2. * Pi * majorRadius(this_xi);
            }

            integrated_volume = integrated_volume + fac*(area_lo + area_hi);
         }
      }
#endif

    
#if CH_SPACEDIM==3
      for (int i = -m_N_mr[0]/2; i < m_N_mr[0]/2; i++) {
         for (int j = -m_N_mr[1]/2; j < m_N_mr[1]/2; j++) {
            for (int k = -m_N_mr[2]/2; k < m_N_mr[2]/2; k++) {
            
               //get subgrid cell centers
               this_xi[RADIAL_DIR] = a_Xi(iv,RADIAL_DIR) + (i + 0.5) * m_dx[0]/m_N_mr[0];
               this_xi[TOROIDAL_DIR] = a_Xi(iv,TOROIDAL_DIR) + (j + 0.5) * m_dx[1]/m_N_mr[1];
               this_xi[POLOIDAL_DIR] = a_Xi(iv,POLOIDAL_DIR) + (k + 0.5) * m_dx[2]/m_N_mr[2];
            
               RealVect Xi00(this_xi[0]-0.5*m_dx[0]/m_N_mr[0], this_xi[1]-0.5*m_dx[1]/m_N_mr[1], this_xi[2]-0.5*m_dx[2]/m_N_mr[2]);
               RealVect Xi10(this_xi[0]+0.5*m_dx[0]/m_N_mr[0], this_xi[1]-0.5*m_dx[1]/m_N_mr[1], this_xi[2]-0.5*m_dx[2]/m_N_mr[2]);
               RealVect Xi01(this_xi[0]-0.5*m_dx[0]/m_N_mr[0], this_xi[1]-0.5*m_dx[1]/m_N_mr[1], this_xi[2]+0.5*m_dx[2]/m_N_mr[2]);
               RealVect Xi11(this_xi[0]-0.5*m_dx[0]/m_N_mr[0], this_xi[1]+0.5*m_dx[1]/m_N_mr[1], this_xi[2]-0.5*m_dx[2]/m_N_mr[2]);

               RealVect X00(realCoord(Xi00));
               RealVect X10(realCoord(Xi10));
               RealVect X01(realCoord(Xi01));
               RealVect X11(realCoord(Xi11));
               
               //making X00 the origin of a Cartesian coordinate system
               X10 -= X00;
               X01 -= X00;
               X11 -= X00;
               
               //compute cell volume (3D volume) as the volume of a parallelepiped
               double this_vol(0.0);
               this_vol += X10[0] * (X01[1]*X11[2] - X11[1]*X01[2]);
               this_vol -= X10[1] * (X01[0]*X11[2] - X11[0]*X01[2]);
               this_vol += X10[2] * (X01[0]*X11[1] - X11[0]*X01[1]);
               
               integrated_volume = integrated_volume + abs(this_vol);
            }
         }
      }

#endif

      a_vol(iv,0) = integrated_volume;
   }
}

void MagBlockCoordSys::getPointwiseMajorRadius( const FArrayBox& a_xi,
                                                FArrayBox&       a_R ) const
{
   const Box& box = a_R.box();
   CH_assert(box == a_xi.box());

   BoxIterator bit(box);
   for (bit.begin();bit.ok();++bit) {
      const IntVect& iv = bit();

      RealVect this_xi;
      for (int comp=0; comp<SpaceDim; ++comp) {
         this_xi[comp] = a_xi(iv,comp);
      }

      a_R(iv,0) = majorRadiusConstMinorRadAdj(this_xi);
   }
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
   CH_TIME("MagBlockCoordSys::getCellCentereddXdxi()");
   CH_assert(a_dXdxi.nComp() == SpaceDim*SpaceDim);
   const Box& box = a_dXdxi.box();

   FArrayBox xi(box,SpaceDim);
   getCellCenteredMappedCoords(xi);

   compute_dXdXi(a_dXdxi, xi, box);
}


void MagBlockCoordSys::getFaceCentereddXdxi( FluxBox& a_dXdxi ) const
{
   CH_TIME("MagBlockCoordSys::getFaceCentereddXdxi()");
   CH_assert(a_dXdxi.nComp() == SpaceDim*SpaceDim);

   for (int dir=0; dir<SpaceDim; ++dir) {
      FArrayBox& this_dXdxi = a_dXdxi[dir];
      const Box& box = this_dXdxi.box();

      FArrayBox xi(box,SpaceDim);
      getFaceCenteredMappedCoords(dir, xi);

      compute_dXdXi(this_dXdxi, xi, box);
   }
}


void MagBlockCoordSys::getEdgeCentereddXdxi( EdgeDataBox& a_dXdxi ) const
{
   CH_TIME("MagBlockCoordSys::getEdgeCentereddXdxi()");
   CH_assert(a_dXdxi.nComp() == SpaceDim*SpaceDim);

   for (int dir=0; dir<SpaceDim; ++dir) {
      FArrayBox& this_dXdxi = a_dXdxi[dir];
      const Box& box = this_dXdxi.box();

      FArrayBox xi(box,SpaceDim);
      getEdgeCenteredMappedCoords(dir, xi);

      compute_dXdXi(this_dXdxi, xi, box);
   }
}



void MagBlockCoordSys::getPointwiseN( FArrayBox& a_N ) const
{
   CH_TIME("MagBlockCoordSys::getPointwiseN()");
   const Box& box = a_N.box();
 
   FArrayBox xi(box,SpaceDim);
   getCellCenteredMappedCoords(xi);

   compute_dXdXiCofac( a_N, xi, box);
   // replaces dXdXi(a_N, xi, 0, 1, 1, box), etc
   // and a_N.negate(1,2);

   if ( m_axisymmetric ) {
      BoxIterator bit(box);
      for (bit.begin();bit.ok();++bit) {
         const IntVect& iv = bit();

         RealVect this_xi;
         this_xi[RADIAL_DIR] = xi(iv,RADIAL_DIR);
         this_xi[POLOIDAL_DIR] = xi(iv,POLOIDAL_DIR);

         double TwoPiRmaj = 2. * Pi * majorRadius(this_xi);

         for (int comp=0; comp<(SpaceDim*SpaceDim); ++comp) {
            a_N(iv,comp) *= TwoPiRmaj;
         }
      }
   }
}



void MagBlockCoordSys::getPointwiseN( FluxBox& a_N ) const
{
   CH_TIMERS("MagBlockCoordSys::getPointwiseN");
   CH_TIMER("compute_dXdXiCoFac", t_compute_dxdxi_cofac);

   for (int dir=0; dir<SpaceDim; ++dir) {
      FArrayBox& this_N = a_N[dir];
      const Box& box = this_N.box();

      FArrayBox xi(box,SpaceDim);
      getFaceCenteredMappedCoords(dir, xi);

      CH_START(t_compute_dxdxi_cofac);
      compute_dXdXiCofac( this_N, xi, box );
      CH_STOP(t_compute_dxdxi_cofac);

      if ( m_axisymmetric ) {
         BoxIterator bit(box);
         for (bit.begin();bit.ok();++bit) {
            const IntVect& iv = bit();

            RealVect this_xi;
            this_xi[RADIAL_DIR] = xi(iv,RADIAL_DIR);
            this_xi[POLOIDAL_DIR] = xi(iv,POLOIDAL_DIR);
            
            double TwoPiRmaj = 2. * Pi * majorRadius(this_xi);

            for (int comp=0; comp<(SpaceDim*SpaceDim); ++comp) {
               this_N(iv,comp) *= TwoPiRmaj;
            }
         }
      }
   }
}


void MagBlockCoordSys::getPointwiseMetrics( FluxBox& a_N ) const
{
   /* This function is supposed to be used instead of Chombo's NewFourthOrderCoordSys::getN()
    for the case where pointwise metrics is set to true. Becasue the indexing
    of the entries in matrices computed by Chombo, i.e., the result
    of Chombo's NewFourthOrderCoordSys::getN(), is different from how the pointwise metrics
    data in indexed in the COGENT classes, we need to transpose the result.
    In other words, the result of this function corresponds to 
    NT*faceArea in the COGENT indexing convention, but to N*faceArea
    in Chombo indexing convetion. Recall: div_X(Flux) = div_Xi(faceArea*NT*F)
    NB: If axisymetry is ON the result of NewFourthOrderCoordSys::getN() will contain
    extra 1/R factor in some of compenets (since computeTransverseFaceMetric 
    and incrementFaceMetricWithEdgeTerm are overwritten in this class). This function does
    not have those extra multipliers, thus applyAxisymmetryCorrections should do nothing 
    if pointwise metrics is set to true and this function is used.
    */
   
   for (int dir=0; dir<SpaceDim; ++dir) {
      FArrayBox& this_N = a_N[dir];
      const Box& box = this_N.box();
      
      FArrayBox xi(box,SpaceDim);
      getFaceCenteredMappedCoords(dir, xi);
      
      FArrayBox tmp(box,SpaceDim*SpaceDim);
      compute_dXdXiCofac( tmp, xi, box );
    
      //perform transpose
      if (SpaceDim == 2) {
      
         this_N.copy(tmp,0,0,1);
         this_N.copy(tmp,2,1,1);
         this_N.copy(tmp,1,2,1);
         this_N.copy(tmp,3,3,1);
         
         if ( m_axisymmetric ) {
            BoxIterator bit(box);
            for (bit.begin();bit.ok();++bit) {
               const IntVect& iv = bit();
            
               RealVect this_xi;
               this_xi[RADIAL_DIR] = xi(iv,RADIAL_DIR);
               this_xi[POLOIDAL_DIR] = xi(iv,POLOIDAL_DIR);
            
               double TwoPiRmaj = 2. * Pi * majorRadius(this_xi);
               for (int comp=0; comp<(SpaceDim*SpaceDim); ++comp) {
                  this_N(iv,comp) *= TwoPiRmaj;
               }
            }
         }
      }
      
      //perform transpose
      if (SpaceDim == 3) {
         this_N.copy(tmp,0,0,1);
         this_N.copy(tmp,3,1,1);
         this_N.copy(tmp,6,2,1);
         this_N.copy(tmp,1,3,1);
         this_N.copy(tmp,4,4,1);
         this_N.copy(tmp,7,5,1);
         this_N.copy(tmp,2,6,1);
         this_N.copy(tmp,5,7,1);
         this_N.copy(tmp,8,8,1);
      }
   }

   
   RealVect faceArea = getMappedFaceArea();
   for (int dir=0; dir<SpaceDim; dir++)
   {
      FArrayBox& this_N = a_N[dir];
      this_N.mult(faceArea[dir]);
   }

}

Vector<Real>
MagBlockCoordSys::getPointwiseNT(const RealVect& a_Xi) const
{
   
   Vector<Real> tmp = compute_dXdXiCofac(a_Xi);
   Vector<Real> result(tmp);

   //perform transpose
   if (SpaceDim == 2) {
      result[1] = tmp[2];
      result[2] = tmp[1];

      if ( m_axisymmetric ) {
         for (int nComp=0; nComp<SpaceDim*SpaceDim; ++nComp) {
            result[nComp] = result[nComp] * 2. * Pi * majorRadius(a_Xi);
         }
      }
   }

   //perform transpose
   if (SpaceDim == 3) {
      result[1] = tmp[3];
      result[2] = tmp[6];
      result[3] = tmp[1];
      result[5] = tmp[7];
      result[6] = tmp[2];
      result[7] = tmp[5];
   }

   return result;
}


void MagBlockCoordSys::getPointwiseNJInverse( FArrayBox& a_NJInverse ) const
{
   const Box& box = a_NJInverse.box();

   FArrayBox xi(box,SpaceDim);
   getCellCenteredMappedCoords(xi);

   compute_dXdXiCofac( a_NJInverse, xi, box );

   //get dXdXi matrix for Jacobian calc
   FArrayBox temp_dXdxi(box, SpaceDim*SpaceDim);
   getCellCentereddXdxi(temp_dXdxi);

   BoxIterator bit(box);
   for (bit.begin();bit.ok();++bit) {
      const IntVect& iv = bit();
      double J;
      if ( SpaceDim == 2 ){
        J = a_NJInverse(iv,0)*a_NJInverse(iv,3) - a_NJInverse(iv,1)*a_NJInverse(iv,2);
      }
      else if ( SpaceDim == 3 ){
        J =   temp_dXdxi(iv,0)*(temp_dXdxi(iv,4)*temp_dXdxi(iv,8) - temp_dXdxi(iv,5)*temp_dXdxi(iv,7));
        J = J-temp_dXdxi(iv,1)*(temp_dXdxi(iv,3)*temp_dXdxi(iv,8) - temp_dXdxi(iv,5)*temp_dXdxi(iv,6));
        J = J+temp_dXdxi(iv,2)*(temp_dXdxi(iv,3)*temp_dXdxi(iv,7) - temp_dXdxi(iv,4)*temp_dXdxi(iv,6));
      }
      else{
        MayDay::Error("MagBlockCoordSys::getPointwiseNJInverse  defined only for SpaceDim = 2 or 3");
      }
      for (int comp=0; comp<(SpaceDim*SpaceDim); ++comp) {
         a_NJInverse(iv,comp) /= J;
      }
   }
}



void MagBlockCoordSys::getPointwiseNJInverse( FluxBox& a_NJInverse ) const
{
   FluxBox temp_dXdxi(a_NJInverse.box(), SpaceDim*SpaceDim);
   getFaceCentereddXdxi(temp_dXdxi);

   for (int dir=0; dir<SpaceDim; ++dir) {
      FArrayBox& NJInverse_dir = a_NJInverse[dir];
      FArrayBox& temp_dXdxi_dir = temp_dXdxi[dir];

      const Box& box = NJInverse_dir.box();

      FArrayBox xi(box,SpaceDim);
      getFaceCenteredMappedCoords(dir, xi);

      compute_dXdXiCofac( NJInverse_dir, xi, box );

      BoxIterator bit(box);
      for (bit.begin();bit.ok();++bit) {
         const IntVect& iv = bit();
         double J; 
         if ( SpaceDim == 2 ){
           J = NJInverse_dir(iv,0)*NJInverse_dir(iv,3) - NJInverse_dir(iv,1)*NJInverse_dir(iv,2);
         }
         else if ( SpaceDim == 3 ){
           J =   temp_dXdxi_dir(iv,0)*(temp_dXdxi_dir(iv,4)*temp_dXdxi_dir(iv,8) - temp_dXdxi_dir(iv,5)*temp_dXdxi_dir(iv,7));
           J = J-temp_dXdxi_dir(iv,1)*(temp_dXdxi_dir(iv,3)*temp_dXdxi_dir(iv,8) - temp_dXdxi_dir(iv,5)*temp_dXdxi_dir(iv,6));
           J = J+temp_dXdxi_dir(iv,2)*(temp_dXdxi_dir(iv,3)*temp_dXdxi_dir(iv,7) - temp_dXdxi_dir(iv,4)*temp_dXdxi_dir(iv,6));
         }
         else{
           MayDay::Error("MagBlockCoordSys::getPointwiseNJInverse  defined only for SpaceDim = 2 or 3");
         }

         for (int comp=0; comp<(SpaceDim*SpaceDim); ++comp) {
            NJInverse_dir(iv,comp) /= J;
         }
      }
   }
}


Vector<Real> 
MagBlockCoordSys::getPointwiseNJInverse(const RealVect& a_Xi) const
{

  Vector<Real> NJInverse = compute_dXdXiCofac(a_Xi);

  Vector<Real> dXdXi = compute_dXdXi(a_Xi);

  double J = computeMatrixDeterminant(dXdXi);
  
   for (int nComp=0; nComp<SpaceDim*SpaceDim; ++nComp) {
      NJInverse[nComp] = NJInverse[nComp] / J;
   }

   return NJInverse;
   
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

void MagBlockCoordSys::pointwiseJ( FluxBox& a_J ) const
{
   CH_assert(a_J.nComp() == 1);
   for (int dir=0; dir<SpaceDim; ++dir) {
      FArrayBox& J_on_dir( a_J[dir] );
      FArrayBox Xi( a_J.box(),SpaceDim );
      getFaceCenteredMappedCoords( dir,Xi );
      pointwiseJ(J_on_dir, Xi, a_J.box());
   }
}

void MagBlockCoordSys::pointwiseJ( EdgeDataBox& a_J ) const
{
   CH_assert(a_J.nComp() == 1);
   for (int dir=0; dir<SpaceDim; ++dir) {
      FArrayBox& J_on_dir( a_J[dir] );
      FArrayBox Xi( a_J.box(),SpaceDim );
      getEdgeCenteredMappedCoords( dir,Xi );
      pointwiseJ(J_on_dir, Xi, a_J.box());
   }
}



void MagBlockCoordSys::getPointwisePoloidalJ( FArrayBox& a_poloidal_J ) const
{
   const Box& box = a_poloidal_J.box();

   FArrayBox xi(box,SpaceDim);
   getCellCenteredMappedCoords(xi);

   FArrayBox dXdXi_fab(box,SpaceDim*SpaceDim);
   compute_dXdXiCofac( dXdXi_fab, xi, box ); //this may be a bug for cfg_dim>2
   
   FArrayBox temp_dXdxi(box, SpaceDim*SpaceDim);
   getCellCentereddXdxi(temp_dXdxi);

   if ( SpaceDim == 3 ) {
      for (BoxIterator bit(box); bit.ok(); ++bit) {
         const IntVect& iv = bit();

         // N is the inverse of the 2x2 matrix dXdXi times its determinant J
         double r;
         r =   temp_dXdxi(iv,0)*(temp_dXdxi(iv,4)*temp_dXdxi(iv,8) - temp_dXdxi(iv,5)*temp_dXdxi(iv,7));
         r = r-temp_dXdxi(iv,1)*(temp_dXdxi(iv,3)*temp_dXdxi(iv,8) - temp_dXdxi(iv,5)*temp_dXdxi(iv,6));
         r = r+temp_dXdxi(iv,2)*(temp_dXdxi(iv,3)*temp_dXdxi(iv,7) - temp_dXdxi(iv,4)*temp_dXdxi(iv,6));

         a_poloidal_J(iv,0) = r;
      }
   }
   else {
      for (BoxIterator bit(box); bit.ok(); ++bit) {
         const IntVect& iv = bit();

         // N is the inverse of the 2x2 matrix dXdXi times its determinant J
         double r = dXdXi_fab(iv,0)*dXdXi_fab(iv,3) - dXdXi_fab(iv,1)*dXdXi_fab(iv,2);

         a_poloidal_J(iv,0) = r;

         if ( m_axisymmetric ) {
            RealVect this_xi;
            this_xi[RADIAL_DIR] = xi(iv,RADIAL_DIR);
            this_xi[POLOIDAL_DIR] = xi(iv,POLOIDAL_DIR);
            
            a_poloidal_J(iv,0) *= 2. * Pi * majorRadius(this_xi);
         }
      }
   }
}


void
MagBlockCoordSys::getNodeCenteredMappedCoords( FArrayBox& a_xi ) const
{
  FORT_GET_NC_MAPPED_COORDS(CHF_BOX(a_xi.box()),
                                   CHF_CONST_REALVECT(m_dx),
                                   CHF_FRA(a_xi));
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
MagBlockCoordSys::getEdgeCenteredMappedCoords( const int a_dir,
                                               FArrayBox&  a_xi ) const
{
  FORT_GET_EC_MAPPED_COORDS(CHF_BOX(a_xi.box()),
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
   CH_TIME("MagBlockCoordSys::getCellCenteredRealCoords()");
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
MagBlockCoordSys::getEdgeCenteredRealCoords( const int a_dir, FArrayBox& a_x ) const
{
   const Box& box( a_x.box() );
   FArrayBox xi_array( box, SpaceDim );
   getEdgeCenteredMappedCoords( a_dir, xi_array );

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


double
MagBlockCoordSys::mappingError( const Box& a_box ) const
{
   CH_assert(a_box.cellCentered());

   double max_diff = 0.;

   RealVect offset = 0.5*m_dx;

   for (BoxIterator bit(a_box); bit.ok(); ++bit) {
      IntVect iv = bit();
      RealVect xi = m_dx*iv;
      xi += offset;

      RealVect new_xi = mappedCoord(realCoord(xi));
      RealVect delta_xi =  new_xi - xi;

      double diff = delta_xi.vectorLength();
      double normalization = 0.5 * (xi.vectorLength() + new_xi.vectorLength());
      if ( normalization > 0. ) diff /= normalization;
      
      if (diff > max_diff) max_diff = diff;
   }

   return max_diff;
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


/*In order to preserve axisymmetry, the Chombo functions
 computeTransverseFaceMetric() and incrementFaceMetricWithEdgeTerm()
 (used in the calculation of Chombo getN) are overriden here with
 some metric terms being multiplied by 1/R. To compensate for that
 we need to call applyAxisymmetricCorrection() [which multiples the
 radial component of an input vector by R] before calling mappedGridDivergence().
 */
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
          Real determinant = computeMatrixDeterminant(jmtrx);
          nodeN(iv,sDir) = determinant;
          if (m_axisymmetric) {
             nodeN(iv,sDir) *= 2. * Pi;
             if (sDir != 0) nodeN(iv,sDir) *= majorRadius(mappedLoc);
          }
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
