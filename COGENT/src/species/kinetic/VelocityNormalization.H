#ifndef _VELOCITY_NORMALIZATION_H_
#define _VELOCITY_NORMALIZATION_H_

#include "Dimensions.H"
#include "FArrayBox.H"
#include "LevelData.H"
//#include "MultiBlockCoordSys.H"
#include "PhaseCoordSys.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "GridFunction.H"
#include "GridFunctionLibrary.H"
#include "MagGeom.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"
//namespace CFG = CFG_NAMESPACE;
//namespace VEL = VEL_NAMESPACE;

/// Velocity normalization
/**
 * This class stores information about velocity spatial normalziation
 * and performs distribution function remapping
 */
class VelocityNormalization
{
   public:

      /// Constructor.
      /**
       */
      VelocityNormalization(const RefCountedPtr<CFG::MagGeom>&     mag_geom,
                            const RefCountedPtr<VEL::VelCoordSys>& vel_coords,
                            const ParmParse&                       ppspecies);

      /// Destructor.
      /**
       */
      virtual ~VelocityNormalization() {;}
   
      void setNormalizationData();
   
      const string get_velocity_norm_type() const {return m_normalization_type;}

      const CFG::LevelData<CFG::FArrayBox>& get_norm_vels() {return m_vel_norm;}

      const CFG::LevelData<CFG::FluxBox>& get_norm_vels_face() {return m_vel_norm_face;}

   private:
   
      void parseParameters(const ParmParse& a_pp);
      
      const RefCountedPtr<CFG::MagGeom> m_mag_geom;
      const RefCountedPtr<VEL::VelCoordSys> m_vel_coords;
   
      CFG::LevelData<CFG::FArrayBox> m_vel_norm;
      CFG::LevelData<CFG::FluxBox> m_vel_norm_face;
      RefCountedPtr<CFG::GridFunction> m_norm_func;
   
      std::string m_normalization_type;
      bool m_include_Bfield_factor;
};

#include "NamespaceFooter.H"

#endif
