#include <float.h>
#include <fstream>
#include "SingleNullBlockCoordSysModel.H"
#include "MagBlockCoordSysF_F.H"
#include "SingleNullBlockCoordSysF_F.H"
#include "HermiteInterp.H"
#include "SplineInterp.H"
#include "BSplineInterp.H"
#include "Directions.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"

#undef REPORT_NEWTON_FAILURE

enum SingleNullBlockType {LCORE,RCORE,LCSOL,RCSOL,LSOL,RSOL,LPF,RPF,MCORE,MCSOL,NUM_SINGLE_NULL_BLOCKS};
enum interpolation_package {HERMITE, SPLINE, BSPLINE};
enum interpolation_mode {FUNCTION, RADIAL_DERIVATIVE, POLOIDAL_DERIVATIVE, SECOND_RADIAL_DERIVATIVE, RADIAL_POLOIDAL_DERIVATIVE, SECOND_POLOIDAL_DERIVATIVE};
enum field_vars {RBR_VAR, RBZ_VAR, RB_VAR, UNIT_BR_VAR, UNIT_BPHI_VAR, UNIT_BZ_VAR, NUM_FIELD_VARS};
enum rz_vars {R_VAR, Z_VAR};

#include "NamespaceHeader.H"

const std::string SingleNullBlockCoordSysModel::pp_name = "singlenullModelGeom";

SingleNullBlockCoordSysModel::SingleNullBlockCoordSysModel( ParmParse&     a_parm_parse,
                                                  const ProblemDomain&     a_domain,
                                                  const RealVect&          a_dx,
                                                  const int                a_block_type )
   : MagBlockCoordSys(a_parm_parse),
     m_block_type(a_block_type),
     m_poloidally_truncated(false)
{
   if (SpaceDim != 2) {
      MayDay::Error("SingleNullBlockCoordSysModel is only two-dimensional");
   }

   if (m_verbose && procID()==0) {
      cout << "Constructing single null ";
      switch (a_block_type)
         {
         case LCORE:
           cout << "left core";
           break;
         case MCORE:
           cout << "middle core";
           break;
         case RCORE:
           cout << "right core";
           break;
         case LCSOL:
           cout << "left center scrape-off layer";
           break;
         case MCSOL:
           cout << "middle center scrape-off layer";
           break;
         case RCSOL:
           cout << "right center scrape-off layer";
           break;
         case LSOL:
           cout << "left scrape_off layer";
           break;
         case RSOL:
           cout << "right scrape-off layer";
           break;
         case LPF:
           cout << "left private flux";
           break;
         case RPF:
           cout << "right private flux";
           break;
         default:
           MayDay::Error("SingleNullBlockCoordSysModelModel::SingleNullBlockCoordSysModel(): Invalid block_type encountered");
         }

      cout << " block with global index space domain box = " << a_domain.domainBox() << endl;
   }

   define( a_domain, a_dx );

   if (m_verbose && procID()==0) {
      switch (a_block_type)
         {
         case LCORE:
           cout << "Left core mapped domain: " << lowerMappedCoordinate(0) << " < xi_0 < " << upperMappedCoordinate(0) << ", "
                << lowerMappedCoordinate(1) << " < xi_1 < " << upperMappedCoordinate(1);
            break;
         case MCORE:
           cout << "Middle core mapped domain: " << lowerMappedCoordinate(0) << " < xi_0 < " << upperMappedCoordinate(0) << ", "
                << lowerMappedCoordinate(1) << " < xi_1 < " << upperMappedCoordinate(1);
            break;
         case RCORE:
            cout << "Right core mapped domain: " << lowerMappedCoordinate(0) << " < xi_0 < " << upperMappedCoordinate(0) << ", "
                 << lowerMappedCoordinate(1) << " < xi_1 < " << upperMappedCoordinate(1);
            break;
         case LCSOL:
            cout << "Left center scrape-off layer mapped domain: " << lowerMappedCoordinate(0) << " < xi_0 < " << upperMappedCoordinate(0)
                 << ", " << lowerMappedCoordinate(1) << " < xi_1 < " << upperMappedCoordinate(1);
            break;
         case MCSOL:
            cout << "Middle center scrape-off layer mapped domain: " << lowerMappedCoordinate(0) << " < xi_0 < " << upperMappedCoordinate(0)
                 << ", " << lowerMappedCoordinate(1) << " < xi_1 < " << upperMappedCoordinate(1);
            break;
         case RCSOL:
            cout << "Right center scrape-off layer mapped domain: " << lowerMappedCoordinate(0) << " < xi_0 < " << upperMappedCoordinate(0)
                 << ", " << lowerMappedCoordinate(1) << " < xi_1 < " << upperMappedCoordinate(1);
            break;
         case LSOL:
            cout << "Left scrape_off layer mapped domain: " << lowerMappedCoordinate(0) << " < xi_0 < " << upperMappedCoordinate(0)
                 << ", " << lowerMappedCoordinate(1) << " < xi_1 < " << upperMappedCoordinate(1);
            break;
         case RSOL:
            cout << "Right scrape-off layer mapped domain: " << lowerMappedCoordinate(0) << " < xi_0 < " << upperMappedCoordinate(0)
                 << ", " << lowerMappedCoordinate(1) << " < xi_1 < " << upperMappedCoordinate(1);
            break;
         case LPF:
            cout << "Left private flux mapped domain: " << lowerMappedCoordinate(0) << " < xi_0 < " << upperMappedCoordinate(0)
                 << ", " << lowerMappedCoordinate(1) << " < xi_1 < " << upperMappedCoordinate(1);
            break;
         case RPF:
            cout << "Right private flux mapped domain " << lowerMappedCoordinate(0) << " < xi_0 < " << upperMappedCoordinate(0)
                 << ", " << lowerMappedCoordinate(1) << " < xi_1 < " << upperMappedCoordinate(1);
            break;
         default:
           MayDay::Error("SingleNullBlockCoordSysModel::SingleNullBlockCoordSysModel(): Invalid block_type encountered");
         }

      cout << endl;
   }

   init( a_parm_parse );
}


SingleNullBlockCoordSysModel::~SingleNullBlockCoordSysModel()
{
}

void
SingleNullBlockCoordSysModel::init( ParmParse& a_pp )
{
   
   //Default values
   m_a =  1.2;
   m_b =  0.9;
   m_c =  0.7;
   m_R0 = 1.6;
   m_Zx = -acos(m_c/m_b);
   m_Zc  = acos(m_c/m_b);
   m_Z0  = 0.4-m_Zx;
   m_core_width = 0.066;
   m_sol_width  = 0.066;
   m_pf_width   = 0.066;
   m_pol_ref_length =  PI;
   m_pol_ref_length_mp = PI/2.0;
   m_div_leg_length = 1.0/2.0 * PI;
   m_pol_ref_z_lo = m_Zx+0.4;
   m_magAxis[0] = m_R0;
   m_magAxis[1] = m_Zc + m_Z0;
   m_Xpoint[0] = m_R0;
   m_Xpoint[1] = m_Zx + m_Z0;
   
   
   //Get the geometry parameters
   a_pp.query( "shape_factor_a",                   m_a );
   a_pp.query( "shape_factor_b",                   m_b );
   a_pp.query( "shape_factor_c",                   m_c );
   a_pp.query( "core_width",                       m_core_width );
   a_pp.query( "sol_width",                        m_sol_width );
   a_pp.query( "pf_width",                         m_pf_width );
   a_pp.query( "divertor_leg_poloidal_length",     m_div_leg_length );
   a_pp.query( "reference_poloidal_plane",         m_pol_ref_z_lo );
   
   //Important to define this solver parameters before doing any calculations!
   m_max_iterations = 1000;
   m_tol = 1.1e-7;
   
   // Get the toroidal field component scale factor
   if (a_pp.contains("Btor_scale")) {
      a_pp.get("Btor_scale", m_RB_toroidal);
   }
   else {
      m_RB_toroidal = 0.;  // default
   }
   
   
   // Get the poloidal field component scale factor
   if (a_pp.contains("Bpol_scale")) {
      a_pp.get("Bpol_scale", m_BpScale);
   }
   else {
      m_BpScale = m_RB_toroidal / m_R0;
   }
   
   normalizeBpScale();
   
   //Derived parameters
   double psiVal = psiAtR0(m_Zx+m_Z0);
   m_Zsep_hi = getZhi(psiVal) - m_Z0;
   m_pol_ref_psi = psiAtR0(m_pol_ref_z_lo + m_Z0);
   
   //Create an interpolator for the cylindrical coordinates
   //Here, we just use this interpolation utility as an effective storage for cell vertices and cell ceneters values (computed analytically)
   //(otherwise we need to compute coordinates (by solving algebraic equations) each time step, and it takes time)
   string interp_method = "bspline";
   if (a_pp.contains("interpolation_method")) {
      a_pp.get("interpolation_method", interp_method);
   }
   
   int m_ghost = 4;
   IntVect interp_box_smallEnd;
   interp_box_smallEnd[RADIAL_DIR] = m_N_mr[RADIAL_DIR] * (m_domain.domainBox().smallEnd()[RADIAL_DIR] - m_ghost);
   interp_box_smallEnd[POLOIDAL_DIR] = m_N_mr[POLOIDAL_DIR] * (m_domain.domainBox().smallEnd()[POLOIDAL_DIR] - m_ghost);
   
   IntVect interp_box_bigEnd;
   interp_box_bigEnd[RADIAL_DIR] = m_N_mr[RADIAL_DIR] * ( m_domain.domainBox().bigEnd()[RADIAL_DIR] + (m_ghost + 1) );
   interp_box_bigEnd[POLOIDAL_DIR] = m_N_mr[POLOIDAL_DIR] * ( m_domain.domainBox().bigEnd()[POLOIDAL_DIR] + (m_ghost + 1) );
   
   
   Box interp_box(interp_box_smallEnd,interp_box_bigEnd);
   
   FArrayBox interp_node_coords(interp_box, SpaceDim);
   
   FArrayBox RZ_data(interp_box, 2);
   
   BoxIterator bit(interp_box);
   for (bit.begin();bit.ok();++bit) {
      IntVect iv = bit();
      RealVect mapped_coord;
      for (int dir=0; dir<SpaceDim; ++dir) {
         mapped_coord[dir] = (m_dx[dir] * iv[dir]) / m_N_mr[dir];
         RZ_data(iv,dir) = realCoordPointwise(mapped_coord)[dir];
         interp_node_coords(iv,dir) = mapped_coord[dir];
      }
   }
   
   if (interp_method == "bspline") {
      m_RZ_interp = new BSplineInterp(a_pp, interp_node_coords, RZ_data);
   }
   else if (interp_method == "hermite") {
      m_RZ_interp = new HermiteInterp(a_pp, interp_node_coords, RZ_data);
   }
   else if (interp_method == "spline") {
      m_RZ_interp = new SplineInterp(a_pp, interp_node_coords, RZ_data);
   }
 
}

RealVect
SingleNullBlockCoordSysModel::realCoord( const RealVect& a_xi ) const
{
   RealVect x;
   D_TERM(x[0] = m_RZ_interp->interpolate(R_VAR, FUNCTION, a_xi);,
          x[1] = m_RZ_interp->interpolate(Z_VAR, FUNCTION, a_xi);,
          x[2] = 0.;)
   
   return x;
}

//Keep it so we not break the code structure, but make trivial
RealVect
SingleNullBlockCoordSysModel::mappedCoord( const RealVect& a_X ) const
{
   return a_X;
}

RealVect
SingleNullBlockCoordSysModel::realCoordPointwise( const RealVect& a_xi ) const
{
   
   double epsilon = 0.1*m_tol;
   RealVect x;
   double Lmin =0.8;
   
   //In computations of dXdXi we may slighly (with the machine presicion accuracy) step over the block boundaries (e.g., lower core boundary) in the mapped coordinate.
   //Becasue realCoord() function have difficulty near all block boundaries (the theta function is undefined at the block boundaries), we want to make sure that the
   //small overstepping is anihilated and special treatment of the block boundaries (see below) is used
   RealVect a_xi_tmp(a_xi);
   RealVect vicinity = epsilon * m_dx;
   for (int dir = 0; dir < SpaceDim; dir++) {
      if (a_xi[dir] > (lowerMappedCoordinate(dir) - vicinity[dir]) && a_xi[dir] < (lowerMappedCoordinate(dir) + vicinity[dir])) a_xi_tmp[dir]=lowerMappedCoordinate(dir);
      if (a_xi[dir] > (upperMappedCoordinate(dir) - vicinity[dir]) && a_xi[dir] < (upperMappedCoordinate(dir) + vicinity[dir])) a_xi_tmp[dir]=upperMappedCoordinate(dir);
   }
   //Same argumets apply for the midplane where theta function is undefined
   double midPlaneMappedPolCoord = 0.5*(lowerMappedCoordinate(POLOIDAL_DIR) + upperMappedCoordinate(POLOIDAL_DIR));
   if (a_xi[1] > (midPlaneMappedPolCoord - vicinity[1]) && a_xi[1] < (midPlaneMappedPolCoord + vicinity[1])) a_xi_tmp[1]=midPlaneMappedPolCoord;
   
   RealVect mapped_block_width;
   mapped_block_width[RADIAL_DIR] = upperMappedCoordinate(RADIAL_DIR) - lowerMappedCoordinate(RADIAL_DIR);
   mapped_block_width[POLOIDAL_DIR] = upperMappedCoordinate(POLOIDAL_DIR) - lowerMappedCoordinate(POLOIDAL_DIR);
   
   if (m_block_type == RCORE || m_block_type == LCORE) {
      
      double z_hi = (m_Zsep_hi - m_core_width + m_Z0) + m_core_width * (a_xi_tmp[RADIAL_DIR]-lowerMappedCoordinate(RADIAL_DIR))/mapped_block_width[RADIAL_DIR];
      double psiVal = psiAtR0(z_hi);
      
      double length = m_pol_ref_length * (a_xi_tmp[POLOIDAL_DIR]-lowerMappedCoordinate(POLOIDAL_DIR)) / mapped_block_width[POLOIDAL_DIR];
      length = length*(1.0 - Lmin*abs(cos(length/2.0))); //mesh packing
      
      if (m_block_type == LCORE)  {
         length = m_pol_ref_length * (upperMappedCoordinate(POLOIDAL_DIR) - a_xi_tmp[POLOIDAL_DIR]) / mapped_block_width[POLOIDAL_DIR];
         length = length*(1.0 - Lmin*abs(cos(length/2.0)));
      }
      
      if (length < 0) length = -length;
      if (length > m_pol_ref_length) length = m_pol_ref_length - (length - m_pol_ref_length);
      
      if (length < epsilon) {
         if (a_xi_tmp[RADIAL_DIR] <= upperMappedCoordinate(RADIAL_DIR) + epsilon) {
            x[0] = m_R0 + epsilon;
            x[1] = getZlo(psiVal);
         }
         else {
            x[0] = 2.0*m_R0 - getXlo(psiVal);
            x[1] = getZlo(psiVal);
         }
      }
      
      else if ((2.0 * length > (m_pol_ref_length - epsilon)) && (2.0 * length < (m_pol_ref_length + epsilon))) {
         x[0] = m_R0 + acos(psiVal - m_b * sin(m_Zc) + m_c * m_Zc)/m_a;
         x[1] = m_Zc + m_Z0;
      }
      
      else if (length > m_pol_ref_length - epsilon) {
         x[0] = m_R0 + epsilon;
         x[1] = z_hi;
      }
      
      else {
         
         RealVect refPolPlaneInters = arcLengthInverse(length);
         double thetaVal = theta(refPolPlaneInters);
         int side = (length < m_pol_ref_length_mp) ? -1 : 1;
         x = gridLinesIntersection(psiVal, thetaVal, side);
         
         if ( (a_xi_tmp[POLOIDAL_DIR]>upperMappedCoordinate(POLOIDAL_DIR)) || (a_xi_tmp[POLOIDAL_DIR]<lowerMappedCoordinate(POLOIDAL_DIR)) ) {
            x[0] = m_R0 - (x[0]-m_R0);
         }
      }
      
      if ( m_block_type == LCORE ) x[0] = m_R0 - (x[0]-m_R0);
      
   }
   
   if (m_block_type == RCSOL || m_block_type == LCSOL) {
      
      double z_hi = (m_Zsep_hi + m_Z0) + m_sol_width * (a_xi_tmp[RADIAL_DIR]-lowerMappedCoordinate(RADIAL_DIR))/mapped_block_width[RADIAL_DIR];
      double psiVal = psiAtR0(z_hi);
      
      double length = m_pol_ref_length * (a_xi_tmp[POLOIDAL_DIR]-lowerMappedCoordinate(POLOIDAL_DIR)) / mapped_block_width[POLOIDAL_DIR];
      length = length*(1.0 - Lmin*abs(cos(length/2.0)));
      
      if (m_block_type == LCSOL) {
         length = m_pol_ref_length * (upperMappedCoordinate(POLOIDAL_DIR) - a_xi_tmp[POLOIDAL_DIR]) / mapped_block_width[POLOIDAL_DIR];
         length = length*(1.0 - Lmin*abs(cos(length/2.0)));
      }
      if (length < 0) length = - length;
      if (length > m_pol_ref_length) length = m_pol_ref_length - (length - m_pol_ref_length);
      
      if (length < epsilon) {
         if (a_xi_tmp[RADIAL_DIR] >= lowerMappedCoordinate(RADIAL_DIR) - epsilon) {
            x[1] = getZlo(psiVal);
            x[0] = getXlo(psiVal);
         }
         else {
            x[0] = m_R0 + epsilon;
            x[1] = getZlo(psiVal);
         }
         
         if ( m_block_type == LCSOL ) x[0] = m_R0 - (x[0]-m_R0);
      }
      
      else if ((2.0 * length > (m_pol_ref_length - epsilon)) && (2.0 * length < (m_pol_ref_length + epsilon))) {
         
         x[0] = m_R0 + acos(psiVal - m_b * sin(m_Zc) + m_c * m_Zc)/m_a;
         x[1] = m_Zc + m_Z0;
         
         if ( m_block_type == LCSOL ) x[0] = m_R0 - (x[0]-m_R0);
      }
      
      else if (length > m_pol_ref_length - epsilon) {
         x[0] = m_R0 + epsilon;
         x[1] = z_hi;
         
         if ( m_block_type == LCSOL ) x[0] = m_R0 - (x[0]-m_R0);
      }
      
      else {
         
         RealVect refPolPlaneInters = arcLengthInverse(length);
         double thetaVal = theta(refPolPlaneInters);
         int side = (length < m_pol_ref_length_mp) ? -1 : 1;
         x = gridLinesIntersection(psiVal, thetaVal, side);
         
         if ( ((a_xi_tmp[POLOIDAL_DIR] <= upperMappedCoordinate(POLOIDAL_DIR)) && (a_xi_tmp[POLOIDAL_DIR]>=lowerMappedCoordinate(POLOIDAL_DIR))) &&  m_block_type == LCSOL) {
            x[0] = m_R0 - (x[0]-m_R0);
         }
         
         if ( (a_xi_tmp[POLOIDAL_DIR] > upperMappedCoordinate(POLOIDAL_DIR)) &&  m_block_type == RCSOL) x[0] = m_R0 - (x[0]-m_R0);
         if ( (a_xi_tmp[POLOIDAL_DIR] < lowerMappedCoordinate(POLOIDAL_DIR)) &&  m_block_type == RCSOL) x[1] = (m_Z0 +m_Zx) - (x[1]-m_Z0-m_Zx);
         if ( (a_xi_tmp[POLOIDAL_DIR] > upperMappedCoordinate(POLOIDAL_DIR)) &&  m_block_type == LCSOL) {
            x[1] = (m_Z0 +m_Zx) - (x[1]-m_Z0-m_Zx);
            x[0] = m_R0 - (x[0]-m_R0);
         }
      }
   }
   
   
   if (m_block_type == RSOL || m_block_type == LSOL) {
      
      double z_hi = (m_Zsep_hi + m_Z0) + m_sol_width * (a_xi_tmp[RADIAL_DIR]-lowerMappedCoordinate(RADIAL_DIR))/mapped_block_width[RADIAL_DIR];
      double psiVal = psiAtR0(z_hi);
      
      double length = m_div_leg_length * (upperMappedCoordinate(POLOIDAL_DIR) - a_xi_tmp[POLOIDAL_DIR]) / mapped_block_width[POLOIDAL_DIR];
      length = length*(1.0 - Lmin*abs(cos(length/2.0)));
      
      if (m_block_type == LSOL) {
         length = m_div_leg_length * (a_xi_tmp[POLOIDAL_DIR] - lowerMappedCoordinate(POLOIDAL_DIR)) / mapped_block_width[POLOIDAL_DIR];
         length = length*(1.0 - Lmin*abs(cos(length/2.0)));
      }
      if (length < 0) length = - length;
      
      RealVect refPolPlaneInters = arcLengthInverse(length);
      //      double thetaVal = theta(refPolPlaneInters);
      
      if (length < epsilon) {
         if (a_xi_tmp[RADIAL_DIR] >= lowerMappedCoordinate(RADIAL_DIR) - epsilon) {
            x[1] = getZlo(psiVal);
            x[0] = getXlo(psiVal);
         }
         else {
            x[1] = getZlo(psiVal);
            x[0] = m_R0 + epsilon;
         }
         
         if ( m_block_type == LSOL ) x[0] = m_R0 - (x[0]-m_R0);
         
         x[1] = (m_Z0 +m_Zx) - (x[1]-m_Z0-m_Zx);
      }
      
      else if ((2.0 * length > (m_pol_ref_length - epsilon)) && (2.0 * length < (m_pol_ref_length + epsilon))) {
         x[0] = m_R0 + acos(psiVal - m_b * sin(m_Zc) + m_c * m_Zc)/m_a;
         x[1] = m_Zc + m_Z0;
         if ( m_block_type == LSOL ) x[0] = m_R0 - (x[0]-m_R0);
         x[1] = (m_Z0 +m_Zx) - (x[1]-m_Z0-m_Zx);
      }
      
      else {
         
         RealVect refPolPlaneInters = arcLengthInverse(length);
         double thetaVal = theta(refPolPlaneInters);
         int side = (length < m_pol_ref_length_mp) ? -1 : 1;
         x = gridLinesIntersection(psiVal, thetaVal, side);
         
         x[1] = (m_Z0 +m_Zx) - (x[1]-m_Z0-m_Zx);
         
         if ( (a_xi_tmp[POLOIDAL_DIR]>=lowerMappedCoordinate(POLOIDAL_DIR)) &&  m_block_type == LSOL) {
            x[0] = m_R0 - (x[0]-m_R0);
         }
         
         else if ((a_xi_tmp[POLOIDAL_DIR]<lowerMappedCoordinate(POLOIDAL_DIR)) &&  m_block_type == LSOL) {
            x[0] = m_R0 - (x[0]-m_R0);
            x[1] = (m_Z0 +m_Zx) - (x[1]-m_Z0-m_Zx); //another reflection to bring the vertical coordinate back to upper blocks
         }
         else if ((a_xi_tmp[POLOIDAL_DIR]>upperMappedCoordinate(POLOIDAL_DIR)) &&  m_block_type == RSOL) {
            x[1] = (m_Z0 +m_Zx) - (x[1]-m_Z0-m_Zx); //another reflection to bring the vertical coordinate back to upper blocks
         }
         
      }
      
   }
   
   if (m_block_type == RPF || m_block_type == LPF) {
      
      double z_hi = (m_Zsep_hi - m_pf_width + m_Z0) + m_pf_width * (a_xi_tmp[RADIAL_DIR]-lowerMappedCoordinate(RADIAL_DIR))/mapped_block_width[RADIAL_DIR];
      double psiVal = psiAtR0(z_hi);
      
      double length = m_div_leg_length * (upperMappedCoordinate(POLOIDAL_DIR) - a_xi_tmp[POLOIDAL_DIR]) / mapped_block_width[POLOIDAL_DIR];
      length = length*(1.0 - Lmin*abs(cos(length/2.0)));
      
      if (m_block_type == LPF) {
         length = m_div_leg_length * (a_xi_tmp[POLOIDAL_DIR] - lowerMappedCoordinate(POLOIDAL_DIR)) / mapped_block_width[POLOIDAL_DIR];
         length = length*(1.0 - Lmin*abs(cos(length/2.0)));
      }
      if (length < 0) length = - length;
      
      if (length < epsilon) {
         if (a_xi_tmp[RADIAL_DIR] <= upperMappedCoordinate(RADIAL_DIR) + epsilon) {
            x[0] = m_R0 + epsilon;
            x[1] = getZlo(psiVal);
         }
         else {
            x[0] = 2.0*m_R0 - getXlo(psiVal);
            x[1] = getZlo(psiVal);
         }
      }
      
      else if ((2.0 * length > (m_pol_ref_length - epsilon)) && (2.0 * length < (m_pol_ref_length + epsilon))) {
         x[0] = m_R0 + acos(psiVal - m_b * sin(m_Zc) + m_c * m_Zc)/m_a;
         x[1] = m_Zc + m_Z0;
      }
      
      else {
         
         RealVect refPolPlaneInters = arcLengthInverse(length);
         double thetaVal = theta(refPolPlaneInters);
         int side = (length < m_pol_ref_length_mp) ? -1 : 1;
         x = gridLinesIntersection(psiVal, thetaVal, side);
         
         if ( ((a_xi_tmp[POLOIDAL_DIR]>upperMappedCoordinate(POLOIDAL_DIR)) && m_block_type == RPF)
             || ((a_xi_tmp[POLOIDAL_DIR]<lowerMappedCoordinate(POLOIDAL_DIR)) && m_block_type == LPF ) ) {
            
            x[0] = m_R0 - (x[0]-m_R0);
         }
         
      }
      
      x[1] = (m_Z0 +m_Zx) - (x[1]-m_Z0-m_Zx);
      
      if ( m_block_type == LPF ) x[0] = m_R0 - (x[0]-m_R0);
      
   }
   
   return x;
}

Real
SingleNullBlockCoordSysModel::dXdXi( const RealVect& a_Xi,
                                     int             a_dirX,
                                     int             a_dirXi ) const
{
   
   RealVect a_Xi_tmp(a_Xi);
   RealVect x = realCoord(a_Xi_tmp);
   
   RealVect Xi_hi(a_Xi_tmp);
   RealVect Xi_lo(a_Xi_tmp);
   
   Xi_hi[a_dirXi] = a_Xi[a_dirXi] +  m_dx[a_dirXi]/m_N_mr[a_dirXi];
   Xi_lo[a_dirXi] = a_Xi[a_dirXi] -  m_dx[a_dirXi]/m_N_mr[a_dirXi];
   
   RealVect x_hi = realCoord(Xi_hi);
   RealVect x_lo = realCoord(Xi_lo);
   
   RealVect metrics;
   for (int dir(0); dir<SpaceDim; ++dir) {
      metrics[dir] = (x_hi[dir] - x_lo[dir])/m_dx[a_dirXi] * m_N_mr[a_dirXi]/2.0;
   }
   
   return metrics[a_dirX];
   
}


void
SingleNullBlockCoordSysModel::dXdXi(FArrayBox&       a_dXdXi,
                                    const FArrayBox& a_Xi,
                                    int              a_destComp,
                                    int              a_dirX,
                                    int              a_dirXi,
                                    const Box&       a_box ) const
{
   BoxIterator bit(a_box);
   for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      RealVect this_Xi;
      for (int dir=0; dir<SpaceDim; ++dir) {
         this_Xi[dir] = a_Xi(iv,dir);
      }
      a_dXdXi(iv,a_destComp) = dXdXi(this_Xi,a_dirX,a_dirXi);
   }
}


bool
SingleNullBlockCoordSysModel::isValid(const RealVect& xi) const
{
   return xi[RADIAL_DIR]   >= lowerMappedCoordinate(RADIAL_DIR)   &&
   xi[RADIAL_DIR]   <= upperMappedCoordinate(RADIAL_DIR)   &&
   xi[POLOIDAL_DIR] >= lowerMappedCoordinate(POLOIDAL_DIR) &&
   xi[POLOIDAL_DIR] <= upperMappedCoordinate(POLOIDAL_DIR);
}


void
SingleNullBlockCoordSysModel::getMagneticFlux( const FArrayBox& a_physical_coordinates,
                                               FArrayBox&       a_magnetic_flux ) const
{
   const Box& box(a_physical_coordinates.box());
   CH_assert(a_magnetic_flux.box().contains(box));

   for (BoxIterator bit(box); bit.ok(); ++bit) {
      IntVect iv = bit();
      
      RealVect phys_coord(iv);
      for (int dir=0; dir<SpaceDim; dir++) {
         phys_coord[dir] = a_physical_coordinates(iv,dir);
      }
      
      a_magnetic_flux(iv,0) = psi(phys_coord);
   }
}


double
SingleNullBlockCoordSysModel::getMagneticFlux( const RealVect& a_physical_coordinate ) const
{
   return psi(a_physical_coordinate);
}


void SingleNullBlockCoordSysModel::getNormMagneticFlux( const FArrayBox& a_physical_coordinates,
                                                        FArrayBox&       a_magnetic_flux ) const

{
   double physFluxOnAxis = getMagneticFlux(m_magAxis);
   double physFluxOnSep = getMagneticFlux(m_Xpoint);
   
   getMagneticFlux(a_physical_coordinates, a_magnetic_flux);
   a_magnetic_flux.plus(-physFluxOnAxis);
   a_magnetic_flux.divide(physFluxOnSep - physFluxOnAxis);
   
}

double SingleNullBlockCoordSysModel::getNormMagneticFlux( const RealVect& a_physical_coordinate ) const
{
   
   double physFluxOnAxis = getMagneticFlux(m_magAxis);
   double physFluxOnSep = getMagneticFlux(m_Xpoint);
   double physFlux = getMagneticFlux(a_physical_coordinate);
   
   return (physFlux - physFluxOnAxis)/(physFluxOnSep - physFluxOnAxis);
}

   
double SingleNullBlockCoordSysModel::getOuterMidplaneCoord(const double psiNorm) const
{
   bool residual_tolerance_satisfied = false;
   int max_iterations = 10000;
   double tol = 1.0e-5;
   
   double tmp_lo_R = m_magAxis[0];

   double tmp_hi_R, half_R, f_lo_R, f_hi_R, f_half_R;

   tmp_hi_R = m_R0 + PI/m_a;
   
      
   int iter = 0;
   while ( iter++ < max_iterations && !residual_tolerance_satisfied ) {
         
      half_R = 0.5 * (tmp_lo_R + tmp_hi_R);
      
#if CH_SPACEDIM == 2
      RealVect lo_R_vect(tmp_lo_R,m_magAxis[1]);
      RealVect hi_R_vect(tmp_hi_R,m_magAxis[1]);
      RealVect half_R_vect(half_R,m_magAxis[1]);
#endif
      
#if CH_SPACEDIM == 3
      
      //Assume that the third direction is a toroidal angle
      //and that the geoemtry is axisymmetric
      double phi = 0.0;
      
      RealVect lo_R_vect(tmp_lo_R,m_magAxis[1],phi);
      RealVect hi_R_vect(tmp_hi_R,m_magAxis[1],phi);
      RealVect half_R_vect(half_R,m_magAxis[1],phi);
#endif
      
      f_lo_R = getNormMagneticFlux(lo_R_vect) - psiNorm;
      f_hi_R = getNormMagneticFlux(hi_R_vect) - psiNorm;
      f_half_R = getNormMagneticFlux(half_R_vect) - psiNorm;
      
      if (f_lo_R * f_half_R < 0) {
         tmp_hi_R = half_R;
      }
      
      else {
         tmp_lo_R = half_R;
      }
      
#if CH_SPACEDIM == 2
      RealVect R_hi_test(half_R + tol,m_magAxis[1]);
      RealVect R_lo_test(half_R - tol,m_magAxis[1]);
#endif
      
#if CH_SPACEDIM == 3
      RealVect R_hi_test(half_R + tol,m_magAxis[1],phi);
      RealVect R_lo_test(half_R - tol,m_magAxis[1],phi);
#endif
      
      if ( (getNormMagneticFlux(R_hi_test) - psiNorm) * (getNormMagneticFlux(R_lo_test) - psiNorm) < 0 ) {
         residual_tolerance_satisfied = true;
      }
      
   }

   return half_R;
}
   
double SingleNullBlockCoordSysModel::getOuterRsep() const
{
   double psiNorm_atSep = getNormMagneticFlux(m_Xpoint);

   return getOuterMidplaneCoord(psiNorm_atSep);
}
   

void
SingleNullBlockCoordSysModel::computeFieldData( const int  a_dir,
                                          FArrayBox& a_BField,
                                          FArrayBox& a_BFieldMag,
                                          FArrayBox& a_BFieldDir,
                                          FArrayBox& a_gradBFieldMag,
                                          FArrayBox& a_curlBFieldDir,
                                          FArrayBox& a_BFieldDirdotcurlBFieldDir,
                                          const bool a_derived_data_only ) const
{
   
   // Get box intersection
   Box box = a_BField.box();
   box    &= a_BFieldMag.box();
   box    &= a_BFieldDir.box();
   box    &= a_gradBFieldMag.box();
   box    &= a_curlBFieldDir.box();
   box    &= a_BFieldDirdotcurlBFieldDir.box();
   
   FArrayBox RZ(box,SpaceDim);
   if (a_dir>=0 && a_dir<SpaceDim) {
      getFaceCenteredRealCoords(a_dir, RZ);
   }
   else {
      getCellCenteredRealCoords(RZ);
   }
   
   
   for (BoxIterator bit( box ); bit.ok(); ++bit) {
      IntVect iv( bit() );
      
      //Compute unshifted coordinate
      double r = RZ(iv,0)-m_R0;
      double z = ((RZ(iv,1)-m_Z0)>m_Zx) ? (RZ(iv,1)-m_Z0) : (m_Zx + (m_Zx+m_Z0 - RZ(iv,1)));
      
      //We only have psi and thus B defined in the upper blocks (z>m_Zx), and need to use symmetry for the lower blocks.
      double sign = ((RZ(iv,1)-m_Z0)>m_Zx) ? 1.0 : -1.0;
      
      //Compute field data
      a_BField(iv,0) = sign * m_BpScale * (m_c - m_b * cos(z))/(r + m_R0);
      a_BField(iv,1) = m_RB_toroidal/(r + m_R0);
      a_BField(iv,2) = - m_a * m_BpScale * sin(m_a * r) /(r + m_R0);
      
      a_BFieldMag(iv,0) = sqrt(pow(a_BField(iv,0),2) + pow(a_BField(iv,1),2) + pow(a_BField(iv,2),2));
      
      a_BFieldDir(iv,0) = a_BField(iv,0)/a_BFieldMag(iv,0);
      a_BFieldDir(iv,1) = a_BField(iv,1)/a_BFieldMag(iv,0);
      a_BFieldDir(iv,2) = a_BField(iv,2)/a_BFieldMag(iv,0);
      
      double coeff = pow(m_BpScale*m_c,2) + pow(m_RB_toroidal,2) - 2.0 * m_b * pow(m_BpScale,2) * m_c * cos(z) + pow(m_b * m_BpScale * cos(z), 2) + pow(m_a * m_BpScale * sin(m_a * r),2);
      
      a_gradBFieldMag(iv,0) = (-coeff + 1.0/2.0 * pow(m_a,3) * pow(m_BpScale,2) * (r + m_R0) * sin(2.0 * m_a * r)) / (pow(r+m_R0,2)*sqrt(coeff));
      a_gradBFieldMag(iv,1) = 0.0;
      a_gradBFieldMag(iv,2) = sign * m_b * pow(m_BpScale,2) * (m_c - m_b * cos(z)) * sin(z)/ ((r+m_R0)*sqrt(coeff));
      
      a_curlBFieldDir(iv, 0) = sign * m_b * pow(m_BpScale,2) * m_RB_toroidal * (m_c - m_b * cos(z)) * sin(z) / pow(coeff,3.0/2.0);
      
      a_curlBFieldDir(iv, 1) = m_BpScale * (m_a * m_a * cos(m_a * r) * (pow(m_BpScale * m_c,2) + pow(m_RB_toroidal,2) - 2.0 * m_b * pow(m_BpScale,2) * m_c * cos(z) + pow(m_b * m_BpScale * cos(z),2)) + m_b * (pow(m_RB_toroidal,2) + pow(m_a * m_BpScale * sin(m_a * r),2)) * sin(z)) / pow(coeff,3.0/2.0);
      
      a_curlBFieldDir(iv, 2) = m_RB_toroidal * (coeff - 1.0/2.0 * pow(m_a,3) * pow(m_BpScale,2) * (r+m_R0) * sin(2.0 * m_a * r)) / ((r+m_R0) * pow(coeff,3.0/2.0));
      
      a_BFieldDirdotcurlBFieldDir(iv, 0) =  m_BpScale * m_RB_toroidal * (m_a * m_a * (r+m_R0) * cos(m_a *r) - m_a * sin(m_a * r) + m_b * (r+m_R0) * sin(z)) / ((r+m_R0) * coeff);
      
   }
}

void
SingleNullBlockCoordSysModel::getNodalFieldData(FArrayBox& a_points,
                                                FArrayBox& a_A,
                                                FArrayBox& a_b,
                                                FArrayBox& a_Bmag) const
{

  Box box = a_points.box();

  FArrayBox RZ(box,SpaceDim);
  for (BoxIterator bit(box); bit.ok(); ++bit) {
    IntVect iv = bit();

    RealVect xi;
    for (int n=0; n<SpaceDim; ++n) {
      xi[n] = a_points(iv,n);
    }

    RealVect X = realCoord(xi);
    for (int n=0; n<SpaceDim; ++n) {
      RZ(iv,n) = X[n];
    }
  }


  for (BoxIterator bit( box ); bit.ok(); ++bit) {
    IntVect iv( bit() );

    //Compute unshifted coordinate                                                                                                                                                            
    double r = RZ(iv,0)-m_R0;
    double z = ((RZ(iv,1)-m_Z0)>m_Zx) ? (RZ(iv,1)-m_Z0) : (m_Zx + (m_Zx+m_Z0 - RZ(iv,1)));

    //We only have psi and thus B defined in the upper blocks (z>m_Zx), and need to use symmetry for the lower blocks.                                                                        
    double sign = ((RZ(iv,1)-m_Z0)>m_Zx) ? 1.0 : -1.0;

    //Compute field data                                                                                                                                                                      
    double BR = sign * m_BpScale * (m_c - m_b * cos(z))/(r + m_R0);
    double Bphi = m_RB_toroidal/(r + m_R0);
    double BZ = - m_a * m_BpScale * sin(m_a * r) /(r + m_R0);

    a_Bmag(iv,0) = sqrt(pow(BR,2) + pow(Bphi,2) + pow(BZ,2));

    a_b(iv,0) = BR/a_Bmag(iv,0);
    a_b(iv,1) = Bphi/a_Bmag(iv,0);
    a_b(iv,2) = BZ/a_Bmag(iv,0);

    a_A(iv,0)=RZ(iv,1)/RZ(iv,0)*m_RB_toroidal;
    a_A(iv,1)= m_BpScale * (cos(m_a * r) + m_b * sin(z) - m_c * z) / (r + m_R0);
    a_A(iv,2)=0.0;

  }
}

Vector<Real>
SingleNullBlockCoordSysModel::computeBField(const RealVect& a_X) const
{
   
   Vector<Real> result(3,0);
   
   //Compute unshifted coordinate
   double r = a_X[0]-m_R0;
   double z = ((a_X[1]-m_Z0)>m_Zx) ? (a_X[1]-m_Z0) : (m_Zx + (m_Zx+m_Z0 - a_X[1]));
   
   //We only have psi and thus B defined in the upper blocks (z>m_Zx), and need to use symmetry for the lower blocks.
   double sign = ((a_X[1]-m_Z0)>m_Zx) ? 1.0 : -1.0;
   
   //Compute field data
   double BField0 = sign * m_BpScale * (m_c - m_b * cos(z))/(r + m_R0);
   double BField1 = m_RB_toroidal/(r + m_R0);
   double BField2 = - m_a * m_BpScale * sin(m_a * r) /(r + m_R0);
   
   result[0] = BField0;
   result[1] = BField1;
   result[2] = BField2;
   
   return result;
}

RealVect
SingleNullBlockCoordSysModel::gridLinesIntersection( const double psiVal, const double thetaVal, const int side ) const
{
   
   // We do not solve here for the intersections on the "critical edges"
   //(i.e., vertical, horizontal -X , and horizonta - midPlane)
   // The critical edges are directly handled in the realCoord function
   double epsilon = m_tol;
   
   bool residual_tolerance_satisfied = false;
   
   RealVect half_Z_sol, tmp_lo_Z_sol, tmp_hi_Z_sol;
   double f_lo_Z, f_hi_Z, f_half_Z, half_Z;
   
   double tmp_lo_Z = getZlo(psiVal);
   double tmp_hi_Z = getZhi(psiVal);
   
   if (side==-1) tmp_hi_Z = m_Zc + m_Z0 - epsilon;
   if (side==1) tmp_lo_Z = m_Zc + m_Z0 + epsilon;
   
   int iter = 0;
   while ( iter++ < m_max_iterations && !residual_tolerance_satisfied ) {
      
      half_Z = 0.5 * (tmp_lo_Z + tmp_hi_Z);
      half_Z_sol[0] = m_R0 + acos(psiVal - m_b * sin(half_Z - m_Z0) + m_c * (half_Z - m_Z0)) / m_a;
      half_Z_sol[1] = half_Z;
      
      tmp_lo_Z_sol[0] = m_R0 + acos(psiVal - m_b * sin(tmp_lo_Z - m_Z0) + m_c * (tmp_lo_Z - m_Z0)) / m_a;
      tmp_lo_Z_sol[1] = tmp_lo_Z;
      
      tmp_hi_Z_sol[0] = m_R0 + acos(psiVal - m_b * sin(tmp_hi_Z - m_Z0) + m_c * (tmp_hi_Z - m_Z0)) / m_a;
      tmp_hi_Z_sol[1] = tmp_hi_Z;
      
      
      f_lo_Z = theta(tmp_lo_Z_sol) - thetaVal;
      f_hi_Z = theta(tmp_hi_Z_sol) - thetaVal;
      f_half_Z = theta(half_Z_sol) - thetaVal;
      
      if (f_lo_Z * f_half_Z < 0) {
         tmp_hi_Z = half_Z;
      }
      
      else {
         tmp_lo_Z = half_Z;
      }
      
      if ( (tmp_hi_Z - tmp_lo_Z) < m_tol ) residual_tolerance_satisfied = true;
      
   }
   
   
   RealVect x;
   
   D_TERM(x[0] = acos(psiVal - m_b * sin(half_Z - m_Z0) + m_c * (half_Z - m_Z0)) / m_a + m_R0;,
          x[1] = half_Z;,
          x[2] = 0.;)
   
   return x;
}

double SingleNullBlockCoordSysModel::getXlo(double psiVal ) const
{
   CH_assert(psiVal <= psiAtR0(m_Zsep_hi + m_Z0));
   double epsilon = m_tol;
   double x;
   
   if (psiVal == psiAtR0(m_Zsep_hi + m_Z0)) {
      x = m_R0 + epsilon;
   }
   
   else {
      x= m_R0 + acos(psiVal - m_b * sin(m_Zx) + m_c * m_Zx)/m_a;
   }
   
   return x;
}

double SingleNullBlockCoordSysModel::getZlo(double psiVal ) const
{
   
   bool residual_tolerance_satisfied = false;
   
   double epsilon = m_tol;
   double half_Z;
   
   if (psiVal < psiAtR0(m_Zsep_hi + m_Z0)) {
      half_Z = m_Zx + m_Z0;
      residual_tolerance_satisfied = true;
   }
   
   double tmp_lo_Z = m_Zx + m_Z0;
   double tmp_hi_Z = m_Zc + m_Z0;
   double f_lo_Z, f_hi_Z, f_half_Z;
   
   int iter = 0;
   while ( iter++ < m_max_iterations && !residual_tolerance_satisfied ) {
      
      half_Z = 0.5 * (tmp_lo_Z + tmp_hi_Z);
      f_lo_Z = psiAtR0(tmp_lo_Z) - psiVal;
      f_hi_Z = psiAtR0(tmp_hi_Z) - psiVal;
      f_half_Z = psiAtR0(half_Z) - psiVal;
      
      if (f_lo_Z * f_half_Z < 0) {
         tmp_hi_Z = half_Z;
      }
      
      else {
         tmp_lo_Z = half_Z;
      }
      
      if ( (psiAtR0(half_Z + m_tol) - psiVal) * (psiAtR0(half_Z - m_tol) - psiVal) < 0 ) {
         residual_tolerance_satisfied = true;
      }
      
   }
   
   // add epsilon to move away from z=Zx (where theta is undefined) and also to make sure that for x_unshifted~0, there will be a solution for x.
   return (half_Z + epsilon);
}

double SingleNullBlockCoordSysModel::getZhi(double psiVal ) const
{
   bool residual_tolerance_satisfied = false;
   
   double tmp_lo_Z = m_Zc + m_Z0;
   double tmp_hi_Z = (3.0*m_Zc + m_Z0); //assume the simulation domain with the ghost extensions be smaller than tmp_hi_Z
   double half_Z, f_lo_Z, f_hi_Z, f_half_Z;
   
   int iter = 0;
   while ( iter++ < m_max_iterations && !residual_tolerance_satisfied ) {
      
      half_Z = 0.5 * (tmp_lo_Z + tmp_hi_Z);
      f_lo_Z = psiAtR0(tmp_lo_Z) - psiVal;
      f_hi_Z = psiAtR0(tmp_hi_Z) - psiVal;
      f_half_Z = psiAtR0(half_Z) - psiVal;
      
      if (f_lo_Z * f_half_Z < 0) {
         tmp_hi_Z = half_Z;
      }
      
      else {
         tmp_lo_Z = half_Z;
      }
      
      if ( (psiAtR0(half_Z + m_tol) - psiVal) * (psiAtR0(half_Z - m_tol) - psiVal) < 0 ) {
         residual_tolerance_satisfied = true;
      }
      
   }
   
   //Check if our solution is outside the analytical limit (i.e., the subsequent evalution of theta fails)
   if ( (psiAtR0(half_Z) - psiVal) < 0.0 ) half_Z += -m_tol;
   
   return half_Z;
}

RealVect
SingleNullBlockCoordSysModel::arcLengthInverse(const double length) const
{
   //Here, length is the poloidal angle that goes counter clockwise with zero being at the bottom core boundary
   
   RealVect x;
   
   double epsilon = m_tol;
   
   double psiVal = m_pol_ref_psi;
   double z_polar_axis = m_Zc + m_Z0; //shifted from the toakamak axis where theta is not defined
   
   bool residual_tolerance_satisfied = false;
   
   double half_r, tmp_lo_r, tmp_hi_r, f_lo_r, f_hi_r, f_half_r, x_lo, x_hi, half_x, z_lo, z_hi, half_z;
   
#if CH_SPACEDIM == 3
   double phi;
#endif

   if ( length > (PI/2.0 - epsilon) && length < (PI/2.0 + epsilon) ) {
      x[0] = acos(m_pol_ref_psi - m_b * sin(z_polar_axis) + m_c * z_polar_axis)/m_a + m_R0;
      x[1] = z_polar_axis;

#if CH_SPACEDIM == 3
      //Assume that the third direction coordinate is a toroidal angle
      //and that the geoemtry is axisymmetric
      phi = 0.0;
      x[2] = phi;
#endif

      residual_tolerance_satisfied = true;
   }
   
   tmp_lo_r = 0.0;
   if ( length > PI/2.0 ) {
      tmp_hi_r = 3.0*m_Zc; //assume that this radius encompasses the top part of the domain
   }
   
   if ( length < PI/2.0 ) {
      tmp_hi_r = fabs(m_Zx + m_Z0 - z_polar_axis); //assume that this radius encompasses the lower half of the domain
   }
   
   int iter = 0;
   while ( iter++ < m_max_iterations && !residual_tolerance_satisfied ) {
      
      x_lo = tmp_lo_r * sin(length) + m_R0;
      z_lo = z_polar_axis - tmp_lo_r * cos(length);
      
      x_hi = tmp_hi_r * sin(length) + m_R0;
      z_hi = z_polar_axis - tmp_hi_r * cos(length);
      
      half_r = 0.5 * (tmp_lo_r + tmp_hi_r);
      half_x = half_r * sin(length) + m_R0;
      half_z = z_polar_axis - half_r * cos(length);

#if CH_SPACEDIM == 2
      f_lo_r = psi(RealVect(x_lo,z_lo)) - psiVal;
      f_hi_r = psi(RealVect(x_hi,z_hi)) - psiVal;
      f_half_r = psi(RealVect(half_x,half_z)) - psiVal;
#endif

#if CH_SPACEDIM == 3
      f_lo_r = psi(RealVect(x_lo,z_lo,phi)) - psiVal;
      f_hi_r = psi(RealVect(x_hi,z_hi,phi)) - psiVal;
      f_half_r = psi(RealVect(half_x,half_z,phi)) - psiVal;
#endif

      if (f_lo_r * f_half_r < 0) {
         tmp_hi_r = half_r;
      }
      
      else {
         tmp_lo_r = half_r;
      }
      
      double x_test_hi = (half_r+m_tol) * sin(length) + m_R0;
      double x_test_lo = (half_r-m_tol) * sin(length) + m_R0;
      double z_test_hi = z_polar_axis - (half_r + m_tol) * cos(length);
      double z_test_lo = z_polar_axis - (half_r - m_tol) * cos(length);

#if CH_SPACEDIM == 2
      if ( (psi(RealVect(x_test_hi,z_test_hi)) - psiVal) * (psi(RealVect(x_test_lo,z_test_lo)) - psiVal) < 0 ) {
         residual_tolerance_satisfied = true;
      }
#endif

#if CH_SPACEDIM == 3
      if ( (psi(RealVect(x_test_hi,z_test_hi,phi)) - psiVal) * (psi(RealVect(x_test_lo,z_test_lo,phi)) - psiVal) < 0 ) {
         residual_tolerance_satisfied = true;
      }
#endif

      x[0] = half_r * sin(length) + m_R0;
      x[1] = z_polar_axis - half_r * cos(length);
      
   }
   
   return x;
   
}


//NB: psi, theta and derivatives experssions are used only for z>m_Zx and symmetry is used to obtain those for z<m_Zx (we assume that psi is symmetric)
double SingleNullBlockCoordSysModel::psi(const RealVect& a_x) const
{
   double x = a_x[0] - m_R0;
   double z = a_x[1] - m_Z0;
   
   if (z<m_Zx) z = 2.0 * m_Zx - z;
   
   //   double psiMin = m_b * sin(m_Zc) - m_c * m_Zc;
   double psi = cos(m_a * x) + m_b * sin(z) - m_c * z;
   
   //Bacause psi is periodic (though periodicity is seen only far outside the domain)
   //we flatten it out here, baising to a constant after the "first period boundary".
   //This is the extra cautionary measure for various Newton solvers.
   //   double psiTruncated = (psi < psiMin) ? psiMin : psi;
   return psi; //psiTruncated;
}

double SingleNullBlockCoordSysModel::theta(const RealVect& a_x) const
{
   double x = a_x[0] - m_R0;
   double z = a_x[1] - m_Z0;
   
   //CH_assert(z >= m_Zx);
   
   double Fx = -log(tan((m_a * abs(x))/2.0))/pow(m_a,2);
   
   double sqr = sqrt(m_b * m_b - m_c * m_c);
   double FyMath = (2.0/sqr) * atanh((m_b + m_c) * tan(z/2.0) / sqr);
   double FyGR = (1.0/sqr) * log((sqr*tan(z/2.0) + m_b - m_c)/(sqr*tan(z/2.0) - m_b + m_c));
   double Fy = ( abs(z) < abs(m_Zx) ) ? FyMath : FyGR;
   double theta_orig = Fx - Fy;
   
   double z_lin = m_Zx + m_a*x/pow((m_b*m_b - m_c*m_c),(1.0/4.0));
   double FyMath_lin = (2.0/sqr) * atanh((m_b + m_c) * tan(z_lin/2.0) / sqr);
   double FyGR_lin = (1.0/sqr) * log((sqr*tan(z_lin/2.0) + m_b - m_c)/(sqr*tan(z_lin/2.0) - m_b + m_c));
   double Fy_lin = ( abs(z_lin) < abs(m_Zx) ) ? FyMath_lin : FyGR_lin;
   double theta_lin_core = Fx - Fy_lin;
   
   
   double psiSep = 1.0 + m_b * sin(m_Zx) - m_c * m_Zx;
   double psi = cos(m_a * x) + m_b * sin(z) - m_c * z;
   
   double Rw = 0.01;
   double Xw = 0.2;
   double r = sqrt(x*x + pow(z-m_Zx,2));
   //double fac_core = exp(-pow(x/Xw,2))*(1-exp(-pow(r/Rw,2)));
   double fac_core = exp(-pow(x/Xw,2)) * (1-exp(-pow(r/Rw,2))) * 0.5*(1.0+tanh((psi-psiSep)/0.01)) * (1.0-tanh((z-m_Zx)*(z-m_Zx)/1.0));
   double theta_new_core =  (1-fac_core) * theta_orig + fac_core * theta_lin_core;
   
#if 0   
   double Zw = 0.2;
   double x_lin = (z - m_Zx)*pow((m_b*m_b - m_c*m_c),(1.0/4.0))/m_a;
   double Fx_lin = -log(tan((m_a * abs(x_lin))/2.0))/pow(m_a,2);
   double fac_sol = exp(-pow((z-m_Zx)/Zw,2))*(1-exp(-pow(r/Rw,2)));
   double theta_lin_sol = Fx_lin - Fy;
   double theta_new_sol =  (1-fac_sol) * theta_orig + fac_sol * theta_lin_sol;
#endif
   
   double theta_new;
   if (psi > psiSep) {
      theta_new = theta_new_core;
   }
   else {
      theta_new = theta_orig;
   }
   
   return theta_new_core;
}


double SingleNullBlockCoordSysModel::psiAtR0(const double a_Z) const
{
   double z = a_Z - m_Z0;
   //CH_assert(z >= m_Zx);
   return 1.0 + m_b*sin(z) - m_c*(z);
}

double SingleNullBlockCoordSysModel::psiZ(const RealVect& a_x) const
{
   double z = a_x[1] - m_Z0;
   //CH_assert(z >= m_Zx);
   return -m_c + m_b * cos(z);
}

double SingleNullBlockCoordSysModel::psiX(const RealVect& a_x) const
{
   double x = a_x[0] - m_R0;
   return -m_a * sin(m_a * x);
}

double SingleNullBlockCoordSysModel::thetaZ(const RealVect& a_x) const
{
   double z = a_x[1] - m_Z0;
   //CH_assert(z >= m_Zx);
   return 1.0/(m_c - m_b * cos(z));
}

double SingleNullBlockCoordSysModel::thetaX(const RealVect& a_x) const
{
   double x = a_x[0] - m_R0;
   return -1.0/(m_a*sin(m_a*x));
}

void SingleNullBlockCoordSysModel::normalizeBpScale() {
   
   double psiSep = psiAtR0(m_Zx + m_Z0);
   double sepTopZ = getZhi(psiSep);

#if CH_SPACEDIM == 2
   RealVect sepTop(m_R0, sepTopZ);
#endif
   
#if CH_SPACEDIM == 3
   //Assume that the third direction is a toroidal angle
   //and that the geoemtry is axisymmetric
   RealVect sepTop(m_R0, sepTopZ, 0.0);
#endif
   
   m_BpScale *= m_R0/sqrt(pow(psiX(sepTop),2) + pow(psiZ(sepTop),2));
}



#include "NamespaceFooter.H"
