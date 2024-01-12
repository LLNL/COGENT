#include "Advection.H"
#include "FourthOrderUtil.H"
#include "ToroidalBlockLevelExchangeCenter.H"

#include "NamespaceHeader.H"

const char* Advection::pp_name = {"advection"};


Advection::Advection( const ParmParse&  a_pp,
                      const MagGeom&    a_geom )
   : EllipticOp(a_pp, a_geom)
{

   // N.B. Presently the operator only works for homogenous BCs
   // this needs to be eventually fixed
   
   m_verbosity = true;
   parseParameters( a_pp );
   if (m_verbosity>0) {
      printParameters();
   }
   
   const DisjointBoxLayout& grids = a_geom.grids();
   m_velocity.define(grids, SpaceDim, IntVect::Unit);

   m_volume_reciprocal.define(grids, 1, IntVect::Zero);
   m_geometry.getCellVolumes(m_volume_reciprocal);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_volume_reciprocal[dit].invert(1.0);
   }
   
   int preconditioner_order = 2;

   if ( m_geometry.shearedMBGeom() ) {
     m_mblx_ptr = new ToroidalBlockLevelExchangeCenter(m_geometry, preconditioner_order, preconditioner_order);
   }
   else {
      if ( m_geometry.extrablockExchange() ) {
         m_mblx_ptr = new MagFluxAlignedMBLEXCenter;
         m_mblx_ptr->define(&m_geometry, preconditioner_order, preconditioner_order);
      }
      else {
         m_mblx_ptr = NULL;
      }
   }


   m_preconditioner      = new MBHypreSolver(a_geom, 1, preconditioner_order, m_mblx_ptr);

   ParmParse pp_precond( ((string)a_pp.prefix() + ".linear_solver.precond").c_str());
   m_preconditioner->setMethodParams(pp_precond);
   m_preconditioner->setConvergenceParams(pp_precond);
   
   if ( pp_precond.query("verbose", m_precond_verbose) == 0 ) m_precond_verbose = false;

   ParmParse pp_precond_precond( ((string)pp_precond.prefix() + ".precond").c_str());
   if ( pp_precond_precond.query("verbose", m_precond_precond_verbose) == 0 ) m_precond_precond_verbose = false;
}
      


Advection::~Advection()
{
   if (m_mblx_ptr) delete m_mblx_ptr;
   if (m_preconditioner) delete m_preconditioner;
}


void
Advection::applyOp(LevelData<FArrayBox>&       a_out,
                   const LevelData<FArrayBox>& a_in,
                   bool                        a_homogeneous )
{
   
   // N.B. Presently we assume homogeneous BCs!!!
   // and use PC matrix to define operator
   
   m_preconditioner->multiplyMatrix(a_in, a_out);
}


void
Advection::applyPCOp(const LevelData<FArrayBox>& a_in_vec,
                     LevelData<FArrayBox>&       a_matvec )
{
   m_preconditioner->multiplyMatrix(a_in_vec, a_matvec);
}

void
Advection::computeBcDivergence( LevelData<FArrayBox>& a_out )
{
   // Assume homogeneous BCs
   
   const DisjointBoxLayout & grids = m_geometry.grids();
   
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      a_out[dit].setVal(0.);
   }
}

void
Advection::setOperatorCoefficients( const LevelData<FluxBox>&  a_velocity_phys)
{
   CH_TIME("Advection::setOperatorCoefficients");
   
   const DisjointBoxLayout & grids = m_geometry.grids();

   CH_assert(a_velocity_phys.ghostVect() >= IntVect::Unit);
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      m_velocity[dit].copy(a_velocity_phys[dit]);
   }
   
   computeBcDivergence(m_bc_divergence);
}

void
Advection::updatePreconditioner(const EllipticOpBC&   a_bc )
{
   
   const DisjointBoxLayout & grids = m_geometry.grids();

   m_preconditioner->initializeMatrix();
   
   LevelData<FluxBox> NTV_normal(grids, 1, IntVect::Unit);
   m_geometry.computeMetricTermProductAverage(NTV_normal, m_velocity, false);

   const LevelData<FluxBox>& par_coeff_mapped = m_geometry.getEllipticOpParCoeffMapped();
   LevelData<FluxBox> tensor_coefficient(grids, SpaceDim*SpaceDim, IntVect::Unit);
   
   LevelData<FArrayBox> beta_coefficient(grids, 1, IntVect::Unit);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      tensor_coefficient[dit].copy(par_coeff_mapped[dit]);
      tensor_coefficient[dit] *= m_linear_shift;
      beta_coefficient[dit].setVal(m_diffusion_coeff);
   }
   m_preconditioner->constructMatrixBlock(0, 0,  m_volume_reciprocal, tensor_coefficient, beta_coefficient, a_bc);
   m_preconditioner->addAdvectionMatrixBlock(0, 0, m_volume_reciprocal, NTV_normal, a_bc);

   m_preconditioner->finalizeMatrix();

}


void
Advection::solvePreconditioner( const LevelData<FArrayBox>& a_in,
                                LevelData<FArrayBox>&       a_out )
{
   m_preconditioner->solve(a_in, a_out, true);
}

void
Advection::parseParameters( const ParmParse&   a_pp )
{
   if ( a_pp.query("linear_shift", m_linear_shift) == 0 ) m_linear_shift = 0.0;
   if ( a_pp.query("diffusion_coeff", m_diffusion_coeff) == 0 ) m_diffusion_coeff = 0.0;
}


void
Advection::printParameters()
{
   if (procID()==0) {
      std::cout << "Advection operator ad-hoc linear shift: " << m_linear_shift << std::endl;
      std::cout << "Advection operator ad-hoc parallel diffusion coefficient: " << m_diffusion_coeff << std::endl;
   }
}

#include "NamespaceFooter.H"
