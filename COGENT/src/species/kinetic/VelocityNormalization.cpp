#include "VelocityNormalization.H"
#include "FourthOrderUtil.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "FourthOrderUtil.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

VelocityNormalization::VelocityNormalization(const RefCountedPtr<CFG::MagGeom>&     a_mag_geom,
                                             const RefCountedPtr<VEL::VelCoordSys>& a_vel_coords,
                                             const ParmParse&                       a_ppspecies)
  : m_mag_geom( a_mag_geom ),
    m_vel_coords( a_vel_coords),
    m_include_Bfield_factor(false),
    m_norm_func(NULL)
{
  m_normalization_type = "None";
  parseParameters(a_ppspecies);
    
  if (m_normalization_type == "local") {
    setNormalizationData();
  }
}

void
VelocityNormalization::setNormalizationData()
{
  /*
    Set normalization data. The normalization function plays
    a role of the effective normalization temperature (Teff).
    That is mu_norm = Teff and vpar_norm = sqrt(Teff).
   */
  
  int order = (m_mag_geom->secondOrder()) ? 2 : 4;
  
  const CFG::DisjointBoxLayout& mag_grids = m_mag_geom->grids();

  // We use fourthOrderCellToFace below thus need at least two ghost
  // Decide later if four cells are still needed for the 4th order
  CFG::IntVect ghostVect = order*CFG::IntVect::Unit;
  CFG::LevelData<CFG::FArrayBox> norm_profile( mag_grids, 1, ghostVect );

  if (m_norm_func != NULL) {
    m_norm_func->assign(norm_profile, *m_mag_geom, 0., false);
    m_mag_geom->fillInternalGhosts(norm_profile);
  }
  
  // Set cell-centered data
  m_vel_norm.define( mag_grids, 2, ghostVect );

  const CFG::LevelData<CFG::FArrayBox>& B_field = m_mag_geom->getCCBFieldMag();
  for (CFG::DataIterator dit(mag_grids); dit.ok(); ++dit) {
    CFG::FArrayBox& this_norm_profile = norm_profile[dit];
    CFG::FArrayBox& this_vel_norm = m_vel_norm[dit];
    this_vel_norm.copy(this_norm_profile,0,1,1);
    for (CFG::BoxIterator bit(this_norm_profile.box()); bit.ok(); ++bit) {
      CFG::IntVect iv = bit();
      this_vel_norm(iv,0) = sqrt(this_norm_profile(iv));
       if (m_include_Bfield_factor) {
          this_vel_norm(iv,1) *= 1.0/B_field[dit](iv,0);
       }
    }
  }
  
  // Set face-centered data
  ghostVect = (order == 4) ? CFG::IntVect::Unit : CFG::IntVect::Zero;
  m_vel_norm_face.define( mag_grids, 2, ghostVect );
  fourthOrderCellToFaceCenters(m_vel_norm_face, m_vel_norm );
  //m_mag_geom->convertCellToFace(m_vel_norm_face, m_vel_norm);
}

void
VelocityNormalization::parseParameters(const ParmParse& a_pp)
{
  a_pp.query("vel_normalization_type", m_normalization_type);
   
  a_pp.query("include_Bfield_factor", m_include_Bfield_factor);
  
  CFG::GridFunctionLibrary* grid_library = CFG::GridFunctionLibrary::getInstance();
  std::string grid_function_name;
  
  if (a_pp.contains("vel_normalization_function")) {
     a_pp.get("vel_normalization_function", grid_function_name );
     m_norm_func = grid_library->find( grid_function_name );
  }
  
}

#include "NamespaceFooter.H"
