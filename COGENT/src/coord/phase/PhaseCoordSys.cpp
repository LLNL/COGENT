#include "PhaseCoordSys.H"

#include "NamespaceHeader.H"


PhaseCoordSys::PhaseCoordSys( const RefCountedPtr<CFG::MagCoordSys>&  a_mag_coords,
                              const RefCountedPtr<VEL::VelCoordSys>&  a_vel_coords,
                              const Vector<ProblemDomain>&            a_domains )
   : m_mag_coords(a_mag_coords),
     m_vel_coords(a_vel_coords)
{
  int num_blocks = a_domains.size();

  m_coordSysVect.resize(num_blocks, NULL);
  m_mappingBlocks.resize(num_blocks);

  for (int block=0; block<num_blocks; ++block) {
    const CFG::MagBlockCoordSys* mag_block_coords
      = (CFG::MagBlockCoordSys *)a_mag_coords->getCoordSys(block);
    PhaseBlockCoordSys * phase_block_coords
      = new PhaseBlockCoordSys( *mag_block_coords, *m_vel_coords, a_domains[block]);

    m_coordSysVect[block] = (NewCoordSys *)phase_block_coords;
    m_mappingBlocks[block] = a_domains[block].domainBox();
  }

  m_gotCoordSysVect = m_gotMappingBlocks = true;
}



PhaseCoordSys::~PhaseCoordSys()
{
  for (int block=0; block<m_coordSysVect.size(); ++block) {
    delete m_coordSysVect[block];
  }
}



#include "NamespaceFooter.H"
