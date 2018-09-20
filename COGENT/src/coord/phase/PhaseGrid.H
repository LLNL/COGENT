#ifndef PHASEGRID_H
#define PHASEGRID_H

#include "Dimensions.H"
#include "List.H"
#include "DisjointBoxLayout.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "Box.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#undef CH_SPACEDIM
#define CH_SPACEDIM VEL_DIM
#include "Box.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"
namespace CFG = CFG_NAMESPACE;
namespace VEL = VEL_NAMESPACE;

using namespace std;

class PhaseGrid
{
public:

   /**
      Constructors, destructor and defines
   */
   PhaseGrid(const Vector<ProblemDomain>&  domains,
             const Vector<IntVect>&        decomps,
             const string&                 mag_geom_type);

   /**
      Destructor.
    */
   virtual ~PhaseGrid();

   const DisjointBoxLayout& disjointBoxLayout() const {return m_dbl;}

   const ProblemDomain& domain() const {return m_dbl.physDomain();}

   void print(const IntVect& ghostVect) const;

   unsigned numConfigBoxes() const {return m_local_config_boxes.size();}

   const CFG::Box& configBox(int box_number) const {CH_assert(box_number < numConfigBoxes());
      return m_local_config_boxes[box_number];}

#ifdef CH_MPI
   const MPI_Comm& configBoxComm(int box_number) const {CH_assert(box_number < numConfigBoxes());
      return m_local_comm[box_number];}
#endif

   const List<VEL::Box>& velocitySlice(int box_number) const {CH_assert(box_number < numConfigBoxes());
      return m_local_velocity_slices[box_number];}

private:

   void getConfigBoxes(const Vector<Box>& boxes,
                       List<CFG::Box>&    config_boxes) const;

#ifdef CH_MPI
   void assignPhaseDecompositionToProcs(const Vector<Box>& boxes,
                                        List<CFG::Box>&    config_boxes,
                                        Vector<int>&       procMap,
                                        Vector<CFG::Box>&  local_config_boxes,
                                        Vector<MPI_Comm>&  comm) const;
#else
   void assignPhaseDecompositionToProcs(const Vector<Box>& boxes,
                                        List<CFG::Box>&    config_boxes,
                                        Vector<int>&       procMap,
                                        Vector<CFG::Box>&  local_config_boxes) const;
#endif

   void getVelocitySlices(const DisjointBoxLayout&  dbl,
                          const Vector<CFG::Box>&   config_boxes,
                          Vector< List<VEL::Box> >& velocity_slices) const;

   void projectPhaseToConfiguration( const Box& a_src,
                                     CFG::Box&  a_dst ) const;

   void projectPhaseToVelocity( const Box& a_src,
                                VEL::Box&  a_dst ) const;

   CFG::IntVect config_restrict( const IntVect& a_x ) const;

   VEL::IntVect vel_restrict( const IntVect& a_x ) const;

   DisjointBoxLayout m_dbl;

   Vector<CFG::Box> m_local_config_boxes;
#ifdef CH_MPI
   Vector<MPI_Comm> m_local_comm;
#endif
   Vector< List<VEL::Box> > m_local_velocity_slices;
};

#include "NamespaceFooter.H"

#endif
