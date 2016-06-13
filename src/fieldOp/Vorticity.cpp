#include <math.h>
#include "Vorticity.H"

#include "FourthOrderUtil.H"
#include "Directions.H"
#include "EdgeToCell.H"
#include "ConstFact.H"
#include "FieldOpF_F.H"

#include "inspect.H"
#include "NamespaceHeader.H" 


Vorticity::Vorticity( ParmParse& a_pp, const int a_verbosity )
   : m_verbosity(a_verbosity),
     m_first_step(true),
     m_num_ghosts(4),
     m_conductivity(-1.0),
     m_Te(1.0)
{
   parseParameters( a_pp );
   if (m_verbosity>0) {
      printParameters();
   }
}


Vorticity::~Vorticity()
{
}


void Vorticity::evalFieldRHS( FieldPtrVect&                      a_rhs,
                              const FieldPtrVect&                a_fields,
                              const FluidSpeciesPtrVect&         a_fluids,
                              const PS::KineticSpeciesPtrVect&   a_kinetic_species,
                              const LevelData<FluxBox>&          a_E_field,
                              const int                          a_fieldVecComp,
                              const Real                         a_time)

// NB: kinetic species input corresponds fo the physical species
{
   // Get solution distribution function (Bstar_par*dfn) for the current species
   // and check the number of ghosts
   const PS::KineticSpecies& soln_species( *(a_kinetic_species[0]) );
   const PS::LevelData<PS::FArrayBox>& soln_dfn( soln_species.distributionFunction() );
   CH_assert(soln_dfn.ghostVect() == m_num_ghosts * PS::IntVect::Unit);

   const DisjointBoxLayout& grids( a_E_field.getBoxes() );
   
   //Compute perpendicular current density
   LevelData<FluxBox> perp_current_density(grids, SpaceDim, IntVect::Unit);
   computePerpCurrentDensity(perp_current_density, a_kinetic_species, a_E_field);
   
   //Compute parallel current density
   LevelData<FluxBox> par_current_density(grids, SpaceDim, IntVect::Unit);
   computeParallelCurrentDensity(par_current_density, a_kinetic_species, a_E_field);

   //Compute total current density
   LevelData<FluxBox> current_density(grids, SpaceDim, IntVect::Unit);
   for (DataIterator dit(current_density.dataIterator()); dit.ok(); ++dit) {
      current_density[dit].copy( perp_current_density[dit] );
      current_density[dit] += par_current_density[dit] ;
   }
   
   //Compute field rhs
   Field& rhs_field( *(a_rhs[a_fieldVecComp]) );
   LevelData<FArrayBox>& rhs_data( rhs_field.data() );
   
   const MagGeom& mag_geom = rhs_field.configurationSpaceGeometry();

   mag_geom.applyAxisymmetricCorrection( current_density );
   mag_geom.computeMappedGridDivergence( current_density, rhs_data, false);
   
   for (DataIterator dit( rhs_data.dataIterator() ); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = mag_geom.getBlockCoordSys(grids[dit]);
      double fac = 1. / block_coord_sys.getMappedCellVolume();
      rhs_data[dit].mult(fac);
   }
   
   m_first_step = false;
}


void
Vorticity::computePerpCurrentDensity( LevelData<FluxBox>&               a_perp_current_density,
                                      const PS::KineticSpeciesPtrVect&  a_species,
                                      const LevelData<FluxBox>&         a_field) const
{
   
   LevelData<FArrayBox> perp_current_density_cell(a_perp_current_density.disjointBoxLayout(),
                                                  SpaceDim,
                                                  m_num_ghosts * IntVect::Unit);
   setZero( perp_current_density_cell );
   
   // Container for individual species charge density
   LevelData<FArrayBox> species_perp_current_density;
   species_perp_current_density.define(perp_current_density_cell);
   
   for (int species(0); species<a_species.size(); species++) {
      
      const PS::KineticSpecies& this_species( *(a_species[species]) );
      const PS::PhaseGeom& phase_geom = this_species.phaseSpaceGeometry();
   
      // Compute the perpendicular current density for this species
      PS::LevelData<PS::FluxBox> field_inj;
      phase_geom.injectConfigurationToPhase( a_field, field_inj);
      this_species.perpCurrentDensity( species_perp_current_density, field_inj );
      
      DataIterator dit( a_perp_current_density.dataIterator() );
      for (dit.begin(); dit.ok(); ++dit) {
         perp_current_density_cell[dit].plus( species_perp_current_density[dit] );
      }
   }
   
   //Get geometry parameters
   const PS::KineticSpecies& species( *(a_species[0]) );
   const PS::PhaseGeom& phase_geom = species.phaseSpaceGeometry();
   const MagGeom& mag_geom = phase_geom.magGeom();
   
   extrapolateAtDomainBnd(mag_geom, perp_current_density_cell);
   
   //NB: we loose two layers of ghost cells by doing this
   //fourthOrderCellToFace(a_perp_current_density, perp_current_density_cell);
   cellToFace(a_perp_current_density, perp_current_density_cell);

   DataIterator dit( a_perp_current_density.dataIterator() );
   for (dit.begin(); dit.ok(); ++dit) {
      for (int dir = 0; dir < SpaceDim; dir++) {
      }
   }
}

void
Vorticity::computeParallelCurrentDensity(LevelData<FluxBox>&               a_parallel_current,
                                         const PS::KineticSpeciesPtrVect&  a_species,
                                         const LevelData<FluxBox>&         a_field) const
{

   //Get geometry parameters
   const PS::KineticSpecies& species( *(a_species[0]) );
   const PS::PhaseGeom& phase_geom = species.phaseSpaceGeometry();
   const MagGeom& mag_geom = phase_geom.magGeom();
   const DisjointBoxLayout& grids( a_field.getBoxes() );
   

   //Compute parallel E-field
   LevelData<FluxBox> E_parallel(grids, SpaceDim, IntVect::Unit);

#if CFG_DIM == 2
   LevelData<FluxBox> E_poloidal_vector(grids, SpaceDim, IntVect::Unit);
   mag_geom.projectPoloidalVector(a_field, E_poloidal_vector);
   for (DataIterator dit(E_parallel.dataIterator()); dit.ok(); ++dit) {
     E_parallel[dit].copy( E_poloidal_vector[dit] );
   }
#else
   for (DataIterator dit(E_parallel.dataIterator()); dit.ok(); ++dit) {
     E_parallel[dit].copy( a_field[dit] );
   }
#endif

   mag_geom.projectOntoParallel(E_parallel);

   //Compute electron density
   LevelData<FArrayBox> electron_density(grids, 1, m_num_ghosts * IntVect::Unit);
   computeIonChargeDensity(electron_density, a_species);
   extrapolateAtDomainBnd(mag_geom, electron_density);
   
   //Compute electron density on faces
   LevelData<FluxBox> electron_density_face(grids, 1,IntVect::Unit);
   //NB: we would loose two layers of ghost cells if did forth-order
   //fourthOrderCellToFace(electron_density_face, electron_density);
   cellToFace(electron_density_face, electron_density);
   
   //Increase the number of components in the density object (copy density into all components)
   LevelData<FluxBox> electron_density_face_tmp(grids, SpaceDim, IntVect::Unit);
   for (DataIterator dit(electron_density.dataIterator()); dit.ok(); ++dit) {
       for (int comp = 0; comp < SpaceDim; comp++) {
          electron_density_face_tmp[dit].copy(electron_density_face[dit], 0, comp, 1);
       }
   }
   
   //Compute electron pressure
   LevelData<FArrayBox> electron_pressure(grids, 1, m_num_ghosts * IntVect::Unit);
   for (DataIterator dit(electron_pressure.dataIterator()); dit.ok(); ++dit) {
      electron_pressure[dit].copy(electron_density[dit]);
      electron_pressure[dit].mult(m_Te);
   }

   //Compute electron pressure parallel gradient
   LevelData<FArrayBox> pressure_grad_mapped(grids, 2, IntVect::Unit);
   //Do second order, loos information in one layer of ghosts
   mag_geom.computeMappedPoloidalGradientWithGhosts(electron_pressure, pressure_grad_mapped,2);

   LevelData<FArrayBox> pressure_grad_phys(grids, 2, IntVect::Unit);
   mag_geom.unmapPoloidalGradient(pressure_grad_mapped, pressure_grad_phys);
   mag_geom.projectOntoParallel(pressure_grad_phys);
   
   //Computing pressure gradient on faces
   LevelData<FluxBox> pressure_grad_phys_face(grids, 2, IntVect::Unit);
   extrapolateAtDomainBnd(mag_geom, pressure_grad_phys);
   cellToFace(pressure_grad_phys_face, pressure_grad_phys);
   
   //Compute parallel current
   for (DataIterator dit(a_parallel_current.dataIterator()); dit.ok(); ++dit) {
      for (int dir = 0; dir < SpaceDim; dir++) {
         a_parallel_current[dit][dir].copy( pressure_grad_phys_face[dit][dir] );
         a_parallel_current[dit][dir].divide( electron_density_face_tmp[dit][dir], 0, 0, SpaceDim );
         a_parallel_current[dit][dir].plus( E_parallel[dit][dir] );
         a_parallel_current[dit][dir].mult( m_conductivity );
      }
   }
}

void
Vorticity::computeIonChargeDensity( LevelData<FArrayBox>&               a_ion_charge_density,
                                    const PS::KineticSpeciesPtrVect&    a_species ) const
{
   // Container for individual species charge density
   LevelData<FArrayBox> ion_species_charge;
   ion_species_charge.define(a_ion_charge_density);
   
   setZero( a_ion_charge_density );
   
   for (int species(0); species<a_species.size(); species++) {
      
      const PS::KineticSpecies& this_species( *(a_species[species]) );
      if ( this_species.charge() < 0.0 ) continue;
      
      // Compute the charge density for this species
      this_species.chargeDensity( ion_species_charge );
      
      DataIterator dit( a_ion_charge_density.dataIterator() );
      for (dit.begin(); dit.ok(); ++dit) {
         a_ion_charge_density[dit].plus( ion_species_charge[dit] );
      }
   }
}


//Second-order version of the 4th order fourthOrderCellToFace
void Vorticity::cellToFace(LevelData<FluxBox>& a_faceData,
                           const LevelData<FArrayBox>& a_cellData) const
{
   
   const DisjointBoxLayout& grids = a_faceData.getBoxes();
   
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      for (int dir=0; dir<SpaceDim; dir++)
      {
         FArrayBox& thisFaceDataDir = a_faceData[dit][dir];
         // compute box over which we can do the averaging
         Box testBox(a_cellData[dit].box());
         // need one cell in each direction for averaging
         testBox.grow(dir, -1);
         testBox.surroundingNodes(dir);
         testBox &= thisFaceDataDir.box();
         FORT_CELLTOFACE2NDORDER(CHF_FRA(thisFaceDataDir),
                                 CHF_CONST_FRA(a_cellData[dit]),
                                 CHF_BOX(testBox),
                                 CHF_INT(dir));
      }
   }
}

void
Vorticity::extrapolateAtDomainBnd(const MagGeom&  a_mag_geom,
                                  LevelData<FArrayBox>& a_data) const
                                  
{

   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
   
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& coord_sys = a_mag_geom.getBlockCoordSys(grids[dit]);
      const ProblemDomain& domain = coord_sys.domain();
      
      secondOrderCellExtrapAtDomainBdry(a_data[dit], grids[dit], domain);

   }
   
   a_mag_geom.fillInternalGhosts(a_data);
}


inline
void Vorticity::setZero( LevelData<FArrayBox>& a_data ) const
{
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      a_data[dit].setVal(0.);
   }
}

inline
void Vorticity::parseParameters( ParmParse& a_pp )
{
   a_pp.query( "conductivity", m_conductivity );
}


inline
void Vorticity::printParameters()
{
   if (procID()==0) {
      std::cout << "Vorticity collisions parameters:" << std::endl;
      std::cout << "  conductivity  =  " << m_conductivity << std::endl;

   }
}

Real Vorticity::computeDt( const FieldPtrVect&         fields,
                           const FluidSpeciesPtrVect&  fluids)
{
   return 1.0/m_conductivity;
}


#include "NamespaceFooter.H"
