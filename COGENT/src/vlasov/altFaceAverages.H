#ifndef _ALTFACEAVERAGES_H_
#define _ALTFACEAVERAGES_H_

#include "PhaseGeom.H"
#include "NamespaceHeader.H"

class AltFaceAverages
{
public:
   
   /// Constructor with initialization.
   /**
    *  @param[in] pp the input database.
    *  @param[in] larmor_number the Larmor number for the operator.
    */
   AltFaceAverages( const PhaseGeom&  phase_geometry);
   
   /// Destructor.
   /**
    */
   //virtual ~altFaceAverages() {};
   virtual ~AltFaceAverages();

   
   void computeFaceAverages(LevelData<FluxBox>&         face_phi,
                            const LevelData<FArrayBox>& cell_phi,
                            const LevelData<FluxBox>&   face_vel,
                            const PhaseGeom&            geom,
                            const std::string&          method,
                            const bool                  mapped_velocity );
   
private:
   LevelData<FluxBox> m_normal_vel;

};

#include "NamespaceFooter.H"

#endif
