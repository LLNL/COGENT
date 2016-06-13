#include "MomentOp.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "CartesianCoordSys.H"
#undef CH_SPACEDIM
#include "Injection.H.transdim"
#include "Slicing.H.transdim"
#include "ReductionOps.H.multidim"
#include "ReductionCopier.H.multidim"
#ifdef CH_SPACEDIM
#undef CH_SPACEDIM
#endif
#define CH_SPACEDIM PDIM

#include "REAL.H"
#include "IntVect.H"
#include "DataIterator.H"
#include "DisjointBoxLayout.H"
#include "ProblemDomain.H"

#include "KineticSpecies.H"
#include "PhaseBlockCoordSys.H"
#include "Directions.H"

#include "NamespaceHeader.H"
namespace VEL = VEL_NAMESPACE;
using namespace CH_MultiDim;



MomentOp& MomentOp::instance() {
   static MomentOp inst;
   return inst;
}


inline
void MomentOp::computeIntegrand( LevelData<FArrayBox>& a_integrand,
                                 const KineticSpecies& a_kinetic_species,
                                 const Kernel&         a_kernel ) const
{
   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   int dfn_ncomp = dfn.nComp();
   int kernel_ncomp = a_kernel.nComponents();
   a_integrand.define( dfn.getBoxes(), dfn_ncomp*kernel_ncomp, dfn.ghostVect());

   // Initialize the integrand with the distribution function.
   DataIterator dit = a_integrand.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
     for (int comp=0; comp<kernel_ncomp; ++comp) {
       a_integrand[dit].copy(dfn[dit],0,comp,1);
     }
   }

   a_kernel.eval( a_integrand, a_kinetic_species );
}

// modifed computeIntegrand to work for arbitrary LevelData<FArrayBox>
inline
void MomentOp::computeIntegrand( LevelData<FArrayBox>& a_integrand,
                                 const KineticSpecies& a_kinetic_species,
                                 const LevelData<FArrayBox>& a_function,
                                 const Kernel& a_kernel ) const
{
   int func_ncomp = a_function.nComp();
   int kernel_ncomp = a_kernel.nComponents();
   a_integrand.define( a_function.getBoxes(), func_ncomp*kernel_ncomp, a_function.ghostVect());

   // Initialize the integrand with the a_function
   DataIterator dit = a_integrand.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      for (int ncomp=0; ncomp<func_ncomp; ++ncomp) {
         a_integrand[dit].copy(a_function[dit],ncomp,ncomp*kernel_ncomp,kernel_ncomp);
      }
   }

   a_kernel.eval( a_integrand, a_kinetic_species );
}

inline
void MomentOp::partialIntegralMu( CP1::LevelData<CP1::FArrayBox>& a_result,
                                  const LevelData<FArrayBox>&     a_integrand,
                                  const ProblemDomain&            a_domain,
                                  const SliceSpec&                a_slice_mu) const
{
   if (a_result.isDefined()) {
      MayDay::Error("MomentOp::partialIntegralMu(): a_result is already defined");
   }

   DisjointBoxLayout degenerate_grids;
   adjCellLo(degenerate_grids, a_integrand.getBoxes(), MU_DIR, -1);

   IntVect degenerate_ghosts(a_integrand.ghostVect());
   degenerate_ghosts[MU_DIR] = 0;

   LevelData<FArrayBox> degenerate_integrand(degenerate_grids, a_integrand.nComp(), degenerate_ghosts);

   // The SumOp below should initialize the destination, but doesn't
   DataIterator dit = degenerate_integrand.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      degenerate_integrand[dit].setVal(0.);
   }

   // sum reduce onto thin LevelData
   const ReductionCopier reduceCopier( a_integrand.getBoxes(),
                                       degenerate_integrand.getBoxes(),
                                       a_domain,
                                       MU_DIR );

   const SumOp op_mu( MU_DIR );
   a_integrand.copyTo( a_integrand.interval(),
                       degenerate_integrand,
                       degenerate_integrand.interval(),
                       reduceCopier,
                       op_mu );

   sliceLevelData( a_result, degenerate_integrand, a_slice_mu );
}


inline
void MomentOp::partialIntegralVp( CFG::LevelData<CFG::FArrayBox>&       a_result,
                                  const CP1::LevelData<CP1::FArrayBox>& a_integrand,
                                  const CP1::ProblemDomain&             a_domain,
                                  const CP1::SliceSpec&                 a_slice_vp)       const
{
   if (a_result.isDefined()) {
      MayDay::Error("MomentOp::partialIntegralVp(): a_result is already defined");
   }

   CP1::DisjointBoxLayout degenerate_grids;
   adjCellLo(degenerate_grids, a_integrand.getBoxes(), VPARALLEL_DIR, -1);

   CP1::IntVect degenerate_ghosts(a_integrand.ghostVect());
   degenerate_ghosts[VPARALLEL_DIR] = 0;

   CP1::LevelData<CP1::FArrayBox> degenerate_integrand(degenerate_grids, a_integrand.nComp(), degenerate_ghosts);

   // The SumOp below should initialize the destination, but doesn't
   CP1::DataIterator dit = degenerate_integrand.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      degenerate_integrand[dit].setVal(0.);
   }

   const CP1::ReductionCopier reduceCopier( a_integrand.getBoxes(),
                                            degenerate_integrand.getBoxes(),
                                            a_domain,
                                            VPARALLEL_DIR );
   const CP1::SumOp op_vp( VPARALLEL_DIR );
   a_integrand.copyTo( a_integrand.interval(),
                       degenerate_integrand,
                       degenerate_integrand.interval(),
                       reduceCopier,
                       op_vp );

   sliceLevelData( a_result, degenerate_integrand, a_slice_vp );
}


inline Real velocitySpaceArea( const VEL::CartesianCS& a_vel_coords )
{
   const Real detJ = a_vel_coords.pointwiseJ( VEL::RealVect::Zero );
   const VEL::RealVect& dx( a_vel_coords.dx() );
   return detJ * dx.product();
}


inline
void scaleResult( CFG::LevelData<CFG::FArrayBox>& a_result,
                  const Real&                     a_factor,
                  const PhaseGeom&                a_coords )
{
   const CFG::DisjointBoxLayout& boxes = a_result.getBoxes();
   CFG::DataIterator dit = boxes.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      a_result[dit].mult( a_factor );
   }
}



void MomentOp::compute( CFG::LevelData<CFG::FArrayBox>& a_result,
                        const KineticSpecies&           a_kinetic_species,
                        const Kernel&                   a_kernel ) const
{
   CH_assert(a_result.nComp()==a_kernel.nComponents()*a_kinetic_species.distributionFunction().nComp());
   LevelData<FArrayBox> integrand;
   computeIntegrand( integrand, a_kinetic_species, a_kernel );

   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();
   const VEL::VelCoordSys& vel_coords = geometry.velSpaceCoordSys();
   const ProblemDomain& domain = geometry.domain();

   SliceSpec slice_mu(MU_DIR, domain.domainBox().smallEnd(MU_DIR));
   CP1::SliceSpec slice_vp(VPARALLEL_DIR, domain.domainBox().smallEnd(VPARALLEL_DIR));

   CP1::LevelData<CP1::FArrayBox> partial_integrand;
   partialIntegralMu( partial_integrand, integrand, domain, slice_mu );

   CP1::ProblemDomain partial_domain( sliceDomain( domain, slice_mu ) );

   CFG::LevelData<CFG::FArrayBox> moment;
   partialIntegralVp( moment, partial_integrand, partial_domain, slice_vp );

   const Real scale = a_kernel.scale( a_kinetic_species );
   const Real area = velocitySpaceArea( vel_coords );
   const Real mass = a_kinetic_species.mass();
   const Real factor = scale * area / mass;
   scaleResult( moment, factor, geometry );

   moment.copyTo( a_result );
}

// modified compute to take moments of arbitrary LevelData<FArrayBox>
void MomentOp::compute( CFG::LevelData<CFG::FArrayBox>& a_result,
                        const KineticSpecies& a_kinetic_species,
                        const LevelData<FArrayBox>& a_function,
                        const Kernel& a_kernel ) const
{
   CH_assert(a_result.nComp()==a_kernel.nComponents()*a_function.nComp());
   LevelData<FArrayBox> integrand;
   computeIntegrand( integrand, a_kinetic_species, a_function, a_kernel );

   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();
   const VEL::VelCoordSys& vel_coords = geometry.velSpaceCoordSys();
   const ProblemDomain& domain = geometry.domain();

   SliceSpec slice_mu(MU_DIR, domain.domainBox().smallEnd(MU_DIR));
   CP1::SliceSpec slice_vp(VPARALLEL_DIR, domain.domainBox().smallEnd(VPARALLEL_DIR));

   CP1::LevelData<CP1::FArrayBox> partial_integrand;
   partialIntegralMu( partial_integrand, integrand, domain, slice_mu );

   CP1::ProblemDomain partial_domain( sliceDomain( domain, slice_mu ) );

   CFG::LevelData<CFG::FArrayBox> moment;
   partialIntegralVp( moment, partial_integrand, partial_domain, slice_vp );

   const Real scale = a_kernel.scale( a_kinetic_species );
   const Real area = velocitySpaceArea( vel_coords );
   const Real mass = a_kinetic_species.mass();
   const Real factor = scale * area / mass;
   scaleResult( moment, factor, geometry );

   moment.copyTo( a_result );
}

#include "NamespaceFooter.H"
