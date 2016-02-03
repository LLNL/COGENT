#include "FieldSolver.H"
#include "FieldSolverF_F.H"
#include "BlockRegister.H"
#include "EdgeToCell.H"
#include "FourthOrderUtil.H"
#include "Directions.H"
#include "FluxSurface.H"


#include "NamespaceHeader.H"


FieldSolver::FieldSolver( const ParmParse& a_pp,
                          const MagGeom&   a_geom )
   : m_geometry(a_geom),
     m_num_potential_ghosts_filled(3),
     m_tol(1.e-6),
     m_max_iter(20),
     m_amg_tol(0.),
     m_amg_max_iter(1),
     m_verbose(false),
     m_num_potential_ghosts_interpolated(3)
{
   const DisjointBoxLayout& grids = m_geometry.grids();

   // Check minimum box size.  Needed for boundary extrapolation.
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         if(grids[dit].size(dir) < 3) {
            MayDay::Error( "Poisson: Minimum box dimension is 3!" );
         }
      }
   }

   if (a_pp.contains("verbose")) {
      a_pp.get("verbose", m_verbose);
   }

   if (a_pp.contains("max_iter")) {
      a_pp.get("max_iter", m_max_iter);
   }

   if (a_pp.contains("tol")) {
      a_pp.get("tol", m_tol);
   }

   if (a_pp.contains("amg_max_iter")) {
      a_pp.get("amg_max_iter", m_amg_max_iter);
   }

   if (a_pp.contains("amg_tol")) {
      a_pp.get("amg_tol", m_amg_tol);
   }

   if (a_pp.contains("second_order")) {
     a_pp.get("second_order", m_second_order);
   }
   else {
     m_second_order = false;
   }

   // If there is more than one block, construct the multiblock exchange object
   if ( m_geometry.coordSysPtr()->numBlocks() > 1 ) {
     m_mblex_potential_Ptr = new MultiBlockLevelExchangeCenter();
     int spaceOrder = 4;
     m_mblex_potential_Ptr->define(&a_geom, m_num_potential_ghosts_filled, spaceOrder);
   }
   else {
     m_mblex_potential_Ptr = NULL;
   }

   m_volume.define(grids, 1, IntVect::Zero);
   m_geometry.getCellVolumes(m_volume);
   m_bc_divergence.define(grids, 1, IntVect::Zero);
}
      


FieldSolver::~FieldSolver()
{
   if (m_mblex_potential_Ptr) delete m_mblex_potential_Ptr;
}



void
FieldSolver::computePotential( LevelData<FArrayBox>&       a_phi,
                               const PotentialBC&          a_bc,
                               const LevelData<FArrayBox>& a_charge_density )
{
   LevelData<FArrayBox> rhs;
   rhs.define(a_charge_density);

   subtractBcDivergence(rhs);

   setZero( a_phi );
   solve( rhs, a_bc, a_phi );
}



void
FieldSolver::computeElectricField( const LevelData<FArrayBox>& a_phi,
                                   const PotentialBC&          a_bc,
                                   LevelData<FArrayBox>&       a_E_field ) const
{
   CH_assert(a_E_field.nComp() == 3);

   const DisjointBoxLayout & grids = m_geometry.grids();

   LevelData<FluxBox> flux_face(grids, SpaceDim, a_E_field.ghostVect());
   computeFaceCenteredField(a_phi, a_bc, flux_face);

   LevelData<FArrayBox> tmp(grids, SpaceDim*SpaceDim, a_E_field.ghostVect());

   // Need to make this fourth-order
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      EdgeToCell(flux_face[dit], tmp[dit]);
   }

   // Average the directions
   LevelData<FArrayBox> flux_cell(grids, SpaceDim, a_E_field.ghostVect());

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      FArrayBox& this_tmp = tmp[dit];
      FArrayBox& this_flux = flux_cell[dit];

      this_flux.setVal(0.);
      for (int comp=0; comp<SpaceDim; comp++) {
        for (int dir=0; dir<SpaceDim; dir++) {
          this_flux.plus(this_tmp, comp*SpaceDim + dir, comp, 1);
        }
      }
      this_flux.mult(1./SpaceDim);
   }

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      FArrayBox& this_E_field = a_E_field[dit];
      FArrayBox& this_flux = flux_cell[dit];

      this_E_field.copy(this_flux, 0, 0, 1);
      this_E_field.setVal(0., 1);
      this_E_field.copy(this_flux, 1, 2, 1);
   }
}



void
FieldSolver::computeElectricField( const LevelData<FArrayBox>& a_phi,
                                   const PotentialBC&          a_bc,
                                   LevelData<FluxBox>&         a_E_field ) const
{
   CH_assert(a_E_field.nComp() == 3);

   LevelData<FluxBox> flux(a_E_field.disjointBoxLayout(), SpaceDim, a_E_field.ghostVect());
   computeFaceCenteredField(a_phi, a_bc, flux);

   DataIterator dit = a_E_field.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FluxBox& this_E_field = a_E_field[dit];
      FluxBox& this_flux = flux[dit];
      for (int dir=0; dir<SpaceDim; ++dir) {
         FArrayBox& this_E_field_dir = this_E_field[dir];
         FArrayBox& this_flux_dir = this_flux[dir];
         this_E_field_dir.copy(this_flux_dir, 0, 0, 1);
         this_E_field_dir.setVal(0., 1);
         this_E_field_dir.copy(this_flux_dir, 1, 2, 1);
      }
   }
}



void
FieldSolver::computeFaceCenteredField( const LevelData<FArrayBox>& a_phi,
                                       const PotentialBC&          a_bc,
                                       LevelData<FluxBox>&         a_field,
                                       const bool                  a_homogeneousBCs ) const
{
   /*
     Computes the face-centered (i.e., pointwise) field:

          a_field = - NJinverse grad a_phi

     It is assumed that the ghost cells of a_phi have already been filled.        

   */

   CH_assert(a_field.ghostVect() == IntVect::Unit);
   LevelData<FluxBox> mapped_field(m_geometry.grids(), SpaceDim, a_field.ghostVect());
   computeFaceCenteredMappedField(a_phi, a_bc, mapped_field, a_homogeneousBCs);

   // Multiply mapped field by NJinverse to put it in the physical frame
   LevelData<FluxBox> NJinverse(m_geometry.grids(), SpaceDim*SpaceDim, a_field.ghostVect());
   m_geometry.getPointwiseNJinverse(NJinverse);

   for (DataIterator dit(a_field.dataIterator()); dit.ok(); ++dit) {

      for (int dir=0; dir<SpaceDim; ++dir) {
         const FArrayBox& this_NJInverse = NJinverse[dit][dir];
         const FArrayBox& this_mapped_field = mapped_field[dit][dir];
         FArrayBox& this_field = a_field[dit][dir];
         FArrayBox tmp(this_field.box(), 1);

         this_field.setVal(0.);
         for (int row=0; row<SpaceDim; ++row) {
            for (int col=0; col<SpaceDim; ++col) {
               tmp.copy(this_NJInverse,SpaceDim*row+col,0,1);
               tmp.mult(this_mapped_field,col,0,1);
               this_field.plus(tmp,0,row,1);
            }
         }
      }
   }
}



void
FieldSolver::computeFaceCenteredMappedField( const LevelData<FArrayBox>& a_phi,
                                             const PotentialBC&          a_bc,
                                             LevelData<FluxBox>&         a_field,
                                             const bool                  a_homogeneousBCs ) const
{
   /*
     Computes the face-centered (i.e., pointwise) field in mapped coordinates:

          a_field = - grad a_phi

     It is assumed that the ghost cells of a_phi have already been filled.        

   */

   CH_assert(a_field.ghostVect() == IntVect::Unit);
   const MagCoordSys* mag_coord_sys = m_geometry.getCoordSys();
   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = mag_coord_sys->boundaries();
   const DisjointBoxLayout& grids = m_geometry.grids();
   const IntVect& field_ghosts = a_field.ghostVect();

   // Compute second-order for now
   //   bool fourthOrder = true;
   //   bool dropOrder = true;
   bool fourthOrder = false;
   bool dropOrder = false;

   IntVect stencil_box_lo(IntVect::Zero);
   IntVect stencil_box_hi;
   if( fourthOrder ) {
      stencil_box_hi = 4*IntVect::Unit;
      CH_assert(a_phi.ghostVect() >= field_ghosts + 2*IntVect::Unit);
   }
   else {
      stencil_box_hi = 2*IntVect::Unit;
      CH_assert(a_phi.ghostVect() >= field_ghosts + IntVect::Unit);
   }
   Box stencil_box(stencil_box_lo, stencil_box_hi);
          
   FArrayBox stencil_values(stencil_box, 1);

   int radius = (stencil_box.size(0)-1)/2;

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
     const Box & box = grids[dit];
     int block_number = mag_coord_sys->whichBlock(box);
     const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(block_number);
     const Box& domainBox = block_coord_sys.domain().domainBox();
     RealVect dx = block_coord_sys.dx();

     const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];

     Box fill_domain_box(domainBox);
     int block_bc_types[2*SpaceDim];
     for (int dir=0; dir<SpaceDim; ++dir) {
        for (SideIterator sit; sit.ok(); ++sit) {
           Side::LoHiSide side = sit();

           if ( !this_block_boundaries[dir+side*SpaceDim].isDomainBoundary() ) {
              fill_domain_box.growDir(dir, side, field_ghosts[dir]);
           }
           else {
              block_bc_types[dir+2*side] = a_bc.getBCType(block_number, dir, side);
           }
        }
     }

     Box fill_box = grow(box,field_ghosts) & fill_domain_box;

     const FArrayBox& this_phi = a_phi[dit];

     for (int dir=0; dir<SpaceDim; ++dir) {
        FArrayBox& this_flux = a_field[dit][dir];

        Box template_box(IntVect(D_DECL6(-radius,-radius,-radius,-radius,-radius,-radius)),
                         IntVect(D_DECL6(radius,radius,radius,radius,radius,radius)));
        template_box.growHi(dir, -1);

        BoxIterator bit(surroundingNodes(fill_box,dir));
        for (bit.begin();bit.ok();++bit) {
           IntVect iv_face = bit();

           IntVect iv(iv_face);
           int side;
           if(fill_box.contains(iv)) {
              side = 0;
           }
           else {
              side = 1;
              iv.shift(dir,-1);
              CH_assert(fill_box.contains(iv));
           }

           Box phi_box = template_box + iv_face;

           for( int s = 0; s < SpaceDim; ++s ) {
              stencil_values.setVal(0.);
              accumFieldStencil(iv, dir, side, s, this_block_boundaries, block_bc_types,
                                dx, domainBox, fourthOrder, dropOrder, stencil_values);

              this_flux(iv_face,s) = 0.;

              BoxIterator bit2(phi_box);
              for (bit2.begin();bit2.ok();++bit2) {
                 IntVect iv2 = bit2();
                 this_flux(iv_face,s) += stencil_values(iv2 - iv + radius*IntVect::Unit) * this_phi(iv2);
              }
           }
        }
     }
   }

   if ( !a_homogeneousBCs) {
      accumMappedBoundaryField(a_bc, a_field);
   }

   a_field.exchange();
}



void
FieldSolver::accumMappedBoundaryField( const PotentialBC&  a_bc,
                                       LevelData<FluxBox>& a_field ) const
{
   // N.B.: This is only second-order for now

   const MagCoordSys* mag_coord_sys = m_geometry.getCoordSys();
   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = mag_coord_sys->boundaries();
   const DisjointBoxLayout& grids = m_geometry.grids();

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
     int block_number = mag_coord_sys->whichBlock(grids[dit]);
     const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(block_number);
     const Box& domainBox = block_coord_sys.domain().domainBox();
     const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];
     const RealVect& dx = block_coord_sys.dx();

     FluxBox& this_field = a_field[dit];

     for (int dir=0; dir<SpaceDim; dir++) {
        FArrayBox& this_field_dir = this_field[dir];

        IntVect trans_dirs = IntVect::Unit;
        trans_dirs[dir] = 0;

        for (int side=0; side<2; ++side) {

           if ( this_block_boundaries[dir+side*SpaceDim].isDomainBoundary() ) {
              Box boundary_box;
              if (side == 0) {
                 boundary_box = bdryLo(domainBox,dir,1);
              }
              else {
                 boundary_box = bdryHi(domainBox,dir,1);
              }
                      
              ProblemDomain box_domain(boundary_box);

              Box overlap = boundary_box & this_field_dir.box();
              if (overlap.ok()) {

                 int bc_type = a_bc.getBCType(block_number, dir, side);
                 CH_assert(bc_type == PotentialBC::DIRICHLET || bc_type == PotentialBC::NEUMANN);

                 Box grown_box = grow(overlap,trans_dirs);
                 FArrayBox bv(grown_box,1);
                 a_bc.apply(m_geometry, grids[dit], 0., dir, side, bv);

                 secondOrderTransExtrapAtDomainBdry(bv, dir, overlap, box_domain);

                 for (int comp=0; comp<SpaceDim; ++comp) {

                    FORT_ACCUM_FIELD_BV_CONTRIB2( CHF_BOX(overlap),
                                                  CHF_CONST_INT(dir),
                                                  CHF_CONST_INT(comp),
                                                  CHF_CONST_INT(side),
                                                  CHF_CONST_INT(bc_type),
                                                  CHF_CONST_REALVECT(dx),
                                                  CHF_CONST_FRA1(bv,0),    
                                                  CHF_FRA1(this_field_dir,comp) );
                 }
              }
           }
        }
     }
   }

   a_field.exchange();
}



void
FieldSolver::addBcDivergence( LevelData<FArrayBox>& data ) const
{
   for (DataIterator dit(data.dataIterator()); dit.ok(); ++dit) {
      data[dit] += m_bc_divergence[dit];
   }
}



void
FieldSolver::subtractBcDivergence( LevelData<FArrayBox>& data ) const
{
   for (DataIterator dit(data.dataIterator()); dit.ok(); ++dit) {
      data[dit] -= m_bc_divergence[dit];
   }
}



void
FieldSolver::accumFieldStencil( const IntVect                           a_index,
                                const int                               a_dir,
                                const int                               a_side,
                                const int                               a_sprime,
                                const Tuple<BlockBoundary, 2*SpaceDim>& a_block_boundaries,
                                const int*                              a_block_bc_types,
                                const RealVect&                         a_dx,
                                const Box&                              a_domainBox,
                                const bool                              a_fourthOrder,
                                const bool                              a_dropOrder,
                                FArrayBox&                              a_stencil_values ) const
{
   int dist2physbdry[] = {2, 2, 2, 2};

   for (int dir=0; dir<SpaceDim; dir++) {
      if ( a_block_boundaries[dir].isDomainBoundary() ) {
         if (a_index[dir]==a_domainBox.smallEnd(dir)) {
            dist2physbdry[dir] = 0;
            dist2physbdry[dir+SpaceDim] = -1;
         }
         if (a_index[dir]==a_domainBox.smallEnd(dir)+1) dist2physbdry[dir] = -1;
      }
      if ( a_block_boundaries[dir+SpaceDim].isDomainBoundary() ) {
         if (a_index[dir]==a_domainBox.bigEnd(dir)) {
            dist2physbdry[dir+SpaceDim] = 0;
            dist2physbdry[dir] = 1;
         }
         if (a_index[dir]==a_domainBox.bigEnd(dir)-1) dist2physbdry[dir+SpaceDim] = 1;
      }
   }

   // Add stencil contribution from flux on side in this direction 

   if (a_fourthOrder) {
        
      int dropOrder = a_dropOrder? 1: 0;
      
      FORT_ACCUM_FIELD_STENCIL4(CHF_CONST_INT(a_dir),
                                CHF_CONST_INT(a_sprime),
                                CHF_CONST_INT(a_side),
                                CHF_CONST_I1D(dist2physbdry,2*SpaceDim),
                                CHF_CONST_REALVECT(a_dx),
                                CHF_CONST_INTVECT(a_index),
                                CHF_CONST_INT(dropOrder),
                                CHF_FRA1(a_stencil_values,0));
   }
   else {

      FORT_ACCUM_FIELD_STENCIL2(CHF_CONST_INT(a_dir),
                                CHF_CONST_INT(a_sprime),
                                CHF_CONST_INT(a_side),
                                CHF_CONST_I1D(dist2physbdry,2*SpaceDim),
                                CHF_CONST_I1D(a_block_bc_types,2*SpaceDim),
                                CHF_CONST_REALVECT(a_dx),
                                CHF_CONST_INTVECT(a_index),
                                CHF_FRA1(a_stencil_values,0));
   }
}



void
FieldSolver::applyOperator( const LevelData<FArrayBox>& a_in,
                            const PotentialBC&          a_bc,
                            LevelData<FArrayBox>&       a_out,
                            const bool                  a_homogeneous_bcs ) const
{
   const DisjointBoxLayout& grids = a_in.disjointBoxLayout();

   LevelData<FArrayBox> phi(grids, 1, 3*IntVect::Unit);

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      phi[dit].copy(a_in[dit]);
   }
   fillInternalGhosts(phi);

   LevelData<FluxBox> flux(grids, SpaceDim, IntVect::Unit);
   computeFaceCenteredField(phi, a_bc, flux, a_homogeneous_bcs);

   // Multiply the field by the unmapped, face-centered GKP coefficients
   multiplyUnmappedCoefficients(flux);

   m_geometry.applyAxisymmetricCorrection(flux);

   m_geometry.fillTransversePhysicalGhosts(flux);

   // Convert to face-averaged
   if (!m_second_order) fourthOrderAverage(flux);

   m_geometry.computeMappedGridDivergence(flux, a_out, !m_second_order);

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      a_out[dit] /= m_volume[dit];
   }
}



void
FieldSolver::setBcDivergence( const PotentialBC&    a_bc,
                              LevelData<FArrayBox>& a_out ) const
{
   const DisjointBoxLayout& grids = a_out.disjointBoxLayout();

   LevelData<FluxBox> mapped_flux(grids, SpaceDim, IntVect::Unit);
   setZero(mapped_flux);
   accumMappedBoundaryField(a_bc, mapped_flux);

   // Multiply mapped flux by NJinverse to put it in the physical frame

   LevelData<FluxBox> NJinverse(m_geometry.grids(), SpaceDim*SpaceDim, IntVect::Unit);
   m_geometry.getPointwiseNJinverse(NJinverse);

   LevelData<FluxBox> flux(grids, SpaceDim, IntVect::Unit);

   for (DataIterator dit(flux.dataIterator()); dit.ok(); ++dit) {

      for (int dir=0; dir<SpaceDim; ++dir) {
         const FArrayBox& this_NJInverse = NJinverse[dit][dir];
         const FArrayBox& this_mapped_flux = mapped_flux[dit][dir];
         FArrayBox& this_flux = flux[dit][dir];
         FArrayBox tmp(this_flux.box(), 1);

         this_flux.setVal(0.);
         for (int row=0; row<SpaceDim; ++row) {
            for (int col=0; col<SpaceDim; ++col) {
               tmp.copy(this_NJInverse,SpaceDim*row+col,0,1);
               tmp.mult(this_mapped_flux,col,0,1);
               this_flux.plus(tmp,0,row,1);
            }
         }
      }
   }

   // Multiply the flux by the unmapped, face-centered coefficients
   multiplyUnmappedCoefficients(flux);

   m_geometry.applyAxisymmetricCorrection(flux);

   // Convert to face-averaged
   if (!m_second_order) fourthOrderAverage(flux);

   // We only compute this to second-order now since the flux has only
   // been computed to second-order in accumMappedBoundaryField().
   //   m_geometry.computeMappedGridDivergence(flux, a_out, true);
   m_geometry.computeMappedGridDivergence(flux, a_out, false);

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      a_out[dit] /= m_volume[dit];
   }
}



void
FieldSolver::computeRadialFSAverage( const LevelData<FluxBox>& a_in,
                                     double&                   a_lo_value,
                                     double&                   a_hi_value,
                                     LevelData<FArrayBox>&     a_out ) const
{
   CH_assert(a_in.nComp() == 1);
   CH_assert(a_in.ghostVect() == IntVect::Zero);
   CH_assert(a_out.ghostVect() == IntVect::Zero);

   const DisjointBoxLayout& grids = m_geometry.grids();
   const MagCoordSys& coords = *m_geometry.getCoordSys();
    
   LevelData<FluxBox> mapped_flux_even(grids, SpaceDim, IntVect::Unit);
   LevelData<FluxBox> mapped_flux_odd(grids, SpaceDim, IntVect::Unit);
    
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      for (int dir = 0; dir < SpaceDim; dir++) {
            
         const FArrayBox& this_flux_input_dir = a_in[dit][dir];
         FArrayBox& this_flux_even_dir = mapped_flux_even[dit][dir];
         FArrayBox& this_flux_odd_dir = mapped_flux_odd[dit][dir];
            
         Box box( this_flux_input_dir.box() );
            
         this_flux_even_dir.setVal(0.);
         this_flux_odd_dir.setVal(0.);
            
         int block_number = coords.whichBlock(grids[dit]);
            
         if ( dir == RADIAL_DIR ) {
                
            if ( block_number < 2 ) {
                    
               BoxIterator bit(box);
               for (bit.begin(); bit.ok(); ++bit) {
                  IntVect iv = bit();
                  if (iv[0]%2 == 0 ) {
                     this_flux_even_dir(iv,0) = this_flux_input_dir(iv,0);
                  }
                  else {
                     this_flux_odd_dir(iv,0) = this_flux_input_dir(iv,0);
                  }
               }
            }
         }
      }
   }
   mapped_flux_even.exchange();
   mapped_flux_odd.exchange();
    
   // Multiply mapped flux by NJinverse to put it in the physical frame
    
   LevelData<FluxBox> NJinverse(m_geometry.grids(), SpaceDim*SpaceDim, IntVect::Unit);
   m_geometry.getPointwiseNJinverse(NJinverse);
    
   LevelData<FluxBox> flux_even(grids, SpaceDim, IntVect::Unit);
   LevelData<FluxBox> flux_odd(grids, SpaceDim, IntVect::Unit);
    
   for (DataIterator dit(flux_even.dataIterator()); dit.ok(); ++dit) {
        
      for (int dir=0; dir<SpaceDim; ++dir) {
         const FArrayBox& this_NJInverse = NJinverse[dit][dir];
         const FArrayBox& this_mapped_flux_even = mapped_flux_even[dit][dir];
         const FArrayBox& this_mapped_flux_odd  = mapped_flux_odd[dit][dir];
         FArrayBox& this_flux_even = flux_even[dit][dir];
         FArrayBox& this_flux_odd  = flux_odd[dit][dir];
         FArrayBox tmp(this_flux_even.box(), 1);
            
         this_flux_even.setVal(0.);
         for (int row=0; row<SpaceDim; ++row) {
            for (int col=0; col<SpaceDim; ++col) {
               tmp.copy(this_NJInverse,SpaceDim*row+col,0,1);
               tmp.mult(this_mapped_flux_even,col,0,1);
               this_flux_even.plus(tmp,0,row,1);
            }
         }
            
         this_flux_odd.setVal(0.);
         for (int row=0; row<SpaceDim; ++row) {
            for (int col=0; col<SpaceDim; ++col) {
               tmp.copy(this_NJInverse,SpaceDim*row+col,0,1);
               tmp.mult(this_mapped_flux_odd,col,0,1);
               this_flux_odd.plus(tmp,0,row,1);
            }
         }
      }
   }
    
   // Multiply the flux by the unmapped, face-centered coefficients
   multiplyUnmappedCoefficients(flux_even);
   multiplyUnmappedCoefficients(flux_odd);
   
   m_geometry.applyAxisymmetricCorrection(flux_even);
   m_geometry.applyAxisymmetricCorrection(flux_odd);
 
   // At this point, we have face-centered fluxes including one
   // ghost cell layer.

   if (!m_second_order) {
     // Convert to face-averaged
     fourthOrderAverage(flux_even);
     fourthOrderAverage(flux_odd);
   }
    
   // Now we have face-averaged fluxes on valid faces.  In one
   // layer of ghost cells, we still have the face-centered values,
   // which approximate face-averaged values to second-order.

   LevelData<FArrayBox> divergence_even(grids, 1, IntVect::Zero);
   m_geometry.computeMappedGridDivergence(flux_even, divergence_even, !m_second_order);
    
   LevelData<FArrayBox> divergence_odd(grids, 1, IntVect::Zero);
   m_geometry.computeMappedGridDivergence(flux_odd, divergence_odd, !m_second_order);
    
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      divergence_even[dit] /= m_volume[dit];
      divergence_odd[dit] /= m_volume[dit];
   }

   FluxSurface fs(m_geometry);
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      a_out[dit].copy(divergence_even[dit]);  // using a_out as a temp here
   }
   fs.averageAndSpread(a_out, divergence_even);

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      a_out[dit].copy(divergence_odd[dit]);  // using a_out as a temp here
   }
   fs.averageAndSpread(a_out, divergence_odd);

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      for (BoxIterator bit(grids[dit]); bit.ok(); ++bit) {
         IntVect iv = bit();
         if (iv[0]%2 == 0 ) {
            a_out[dit](iv,0) = -divergence_even[dit](iv,0);
         }
         else { 
            a_out[dit](iv,0) = -divergence_odd[dit](iv,0);
         }
      }
   }

   const ProblemDomain& domain0 = m_geometry.getBlockCoordSys(0).domain();
   const Box& domain_box = domain0.domainBox();

   a_lo_value = a_hi_value = -DBL_MAX;
   const MagCoordSys& coord_sys( *(m_geometry.getCoordSys()) );
   
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      const Box& box = grids[dit];
      int block_number = coord_sys.whichBlock(box);
      if (block_number < 2) {
         IntVect iv = box.smallEnd();
         if ( iv[RADIAL_DIR] == domain_box.smallEnd(RADIAL_DIR) ) {
            if (iv[0]%2 == 0 ) {
               a_lo_value = -divergence_even[dit](iv,0);
            }
            else { 
               a_lo_value = -divergence_odd[dit](iv,0);
            }
         }
         iv = box.bigEnd();
         if ( iv[RADIAL_DIR] == domain_box.bigEnd(RADIAL_DIR) ) {
            if (iv[0]%2 == 0 ) {
               a_hi_value = divergence_odd[dit](iv,0);
            }
            else { 
               a_hi_value = divergence_even[dit](iv,0);
            }
         }
      }
   }
    
   a_lo_value = globalMax(a_lo_value);
   a_hi_value = globalMax(a_hi_value);
}


void
FieldSolver::solve( const LevelData<FArrayBox>& a_rhs,
                    const PotentialBC&          a_bc,
                    LevelData<FArrayBox>&       a_solution )
{
   //   if ( a_bc.hasNeumannCondition() ) {
      GMRESSolve(a_rhs, a_bc, a_solution);
      //   }
      //   else {
      //      PCGSolve(a_rhs, a_bc, a_solution);
      //   }
}



int
FieldSolver::PCGSolve( const LevelData<FArrayBox>& a_rhs,
                       const PotentialBC&          a_bc,
                       LevelData<FArrayBox>&       a_solution )
{
   return OperatorBasedPCGSolve(a_rhs, a_bc, a_solution);
}



int
FieldSolver::OperatorBasedPCGSolve( const LevelData<FArrayBox>& a_rhs,
                                    const PotentialBC&          a_bc,
                                    LevelData<FArrayBox>&       a_solution )
{
   double eps = m_tol;
   int itmax = m_max_iter;

   const DisjointBoxLayout& grids = m_geometry.grids();

   LevelData<FArrayBox> z(grids, 1, IntVect::Zero);
   LevelData<FArrayBox> p(grids, 1, IntVect::Zero);
   LevelData<FArrayBox> r(grids, 1, IntVect::Zero);

   double bnorm;

   // Calculate norm of right hand size.  If 0, return null solution
   if( (bnorm = L2Norm(a_rhs)) == 0.0 ) {
      setZero(a_solution);
      return 1;
   }

   // Calculate r = b - Ax and initial error
   applyOperator(a_solution, a_bc, r, true);

   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      r[dit].negate();
      r[dit].plus(a_rhs[dit]);
   }

   double perr = L2Norm(r) / bnorm;

   if(m_verbose && procID()==0) {
      cout << "   PCG iteration = 0, relative residual norm = " << perr << endl;
   }

   if( perr < eps ) return 0;

   // Iterate
   for( int i=0; i<itmax; i++ ) {

      // Solve Pz = r
      setZero(z);
      solvePreconditioner(r, z);

      // Calculate bknum = (z, Mz) and p
      double bknum = dotProduct( z, r );
      double bkden;
      if( i == 0 ) {                         // First iteration
         // p = z */
         for (dit.begin(); dit.ok(); ++dit) {
            p[dit].copy(z[dit]);
         }
      }
      else {
         // p = z + bk*p
         for (dit.begin(); dit.ok(); ++dit) {
            p[dit] *= bknum/bkden;
            p[dit].plus(z[dit]);
         }
      }
      bkden = bknum;

      // Calculate z = Ap
      applyOperator(p, a_bc, z, true);

      // Compute akden = (p, Ap) and ak
      double akden = dotProduct( p, z );
      double ak    = bknum/akden;
          
      // Update x and r. Calculate error
      for (dit.begin(); dit.ok(); ++dit) {
         a_solution[dit].plus(p[dit], ak);
         r[dit].plus(z[dit], -ak);
      }
      perr = L2Norm(r) / bnorm;

      if (m_verbose && procID()==0) {
         cout << "   PCG iteration = " << i+1 << ", relative residual norm = " << perr << endl;
      }

      if( perr < eps ) return i+1;
  }

  return itmax;
}



int
FieldSolver::GMRESSolve( const LevelData<FArrayBox>& a_rhs,
                         const PotentialBC&          a_bc,
                         LevelData<FArrayBox>&       a_solution )
{
   return OperatorBasedGMRESSolve(a_rhs, a_bc, a_solution);
}



int
FieldSolver::OperatorBasedGMRESSolve( const LevelData<FArrayBox>& a_rhs,
                                      const PotentialBC&          a_bc,
                                      LevelData<FArrayBox>&       a_solution )
{
   int max_krylov_dim = 10;
   int max_num_restarts = 5;
   bool printit = m_verbose > 0 && procID() == 0;

   LevelData<FArrayBox> w, av, r;
   w.define(a_solution);
   av.define(a_solution);
   r.define(a_solution);

   int m = max_krylov_dim;

   double** H = new double*[m+1];
   for (int i=0; i<=m; ++i) {
      H[i] = new double[m];
   }
   double* s = new double[m+1];
   double* cs = new double[m+1];
   double* sn = new double[m+1];

   setZero(w);
   solvePreconditioner(a_rhs, w);
   double normb = L2Norm(w); // normb = ||M b||
   if (normb == 0.0) {
      normb = 1;
   }

   double tol = m_tol * normb;

   applyOperator(a_solution, a_bc, r, true);
   for (DataIterator dit(a_solution.dataIterator()); dit.ok(); ++dit) {
      w[dit].copy(a_rhs[dit]);   // w = b - r
      w[dit] -= r[dit];
   }
   setZero(r);
   solvePreconditioner(w, r);
   double beta = L2Norm(r);  // beta = ||r||

   double resid = beta;
   int itmax = 0;

   if (resid <= tol) {
      if (printit) {
         cout << "   GMRES: initial residual " << resid << " less than tolerance " << tol << endl;
         }
      for (int i=0; i<=m; i++) {
         if (H[i]) delete [] H[i];
      }
      delete [] sn;
      delete [] cs;
      delete [] s;
      delete [] H;
      return itmax;
   }

   if (printit)
      cout << "   GMRES pass " << setw(2) << 1
           << ", iteration " << setw(3) << 0
           << ":  relative residual = " << resid/normb << endl;

   Vector<LevelData<FArrayBox> *> v(m+1);
   for (int i=0; i<=m; ++i) {
      v[i] = NULL;
   }

   int j = 1;
   while (j <= max_num_restarts) {

      itmax++;

      if (v[0] == NULL)  {
         v[0] = new LevelData<FArrayBox>;
         v[0]->define(a_solution);
      }
      for (DataIterator dit(a_solution.dataIterator()); dit.ok(); ++dit) {
         v[0]->operator[](dit).copy(r[dit]);
         v[0]->operator[](dit) /= beta;
      }

      for (int k=0; k<m+1; ++k) s[k] = 0.;
      s[0] = beta;

      int i;
      for (i=0; i<m && j<=max_num_restarts; ++i) {

         applyOperator(*v[i], a_bc, av, true);
         setZero(w);
         solvePreconditioner(av, w);    // w = M A v[i]

         for (int k=0; k<=i; ++k) {
            H[k][i] = dotProduct(w, *v[k]);    // H(k,i) = w * v[k]
            for (DataIterator dit(a_solution.dataIterator()); dit.ok(); ++dit) {
               w[dit].plus((*(v[k]))[dit], -H[k][i]);  // w -= H(k,i) * v[k]
            }
         }

         H[i+1][i] = L2Norm(w);     // H(i+1,i) = ||w||
         if (v[i+1] == NULL)  {
            v[i+1] = new LevelData<FArrayBox>;
            v[i+1]->define(a_solution);
         }

         for (DataIterator dit(a_solution.dataIterator()); dit.ok(); ++dit) {
            v[i+1]->operator[](dit).copy(w[dit]);
            v[i+1]->operator[](dit) /= H[i+1][i];
         }

         for (int k=0; k<i; ++k)
            ApplyPlaneRotation(H[k][i], H[k+1][i], cs[k], sn[k]);

         GeneratePlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
         ApplyPlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
         ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);

         resid = fabs(s[i+1]);
         if (printit)
            cout << "   GMRES pass " << setw(2) << j
                 << ", iteration " << setw(3) << i+1
                 << ": relative residual = " << resid/normb << endl;

         if ( resid < tol) {

            Update(a_solution, i, H, s, m+1, v);

            for (int i=0; i<=m; ++i) {
               if (v[i])  delete v[i];
               if (H[i]) delete [] H[i];
            }
            delete [] sn;
            delete [] cs;
            delete [] s;
            delete [] H;
            return itmax;
         }
      }

      //      if (printit)
         //         cout << "Restarting..." << endl;

      Update(a_solution, i-1, H, s, m+1, v);

      applyOperator(a_solution, a_bc, r, true);
      for (DataIterator dit(a_solution.dataIterator()); dit.ok(); ++dit) {
         w[dit].copy(a_rhs[dit]);     // w = b - r
         w[dit] -= r[dit];
      }
      setZero(r);
      solvePreconditioner(w, r);
      beta = L2Norm(r);         // beta = ||r||
      if ( resid < tol) {
         for (int i=0; i<=m; ++i) {
            if (v[i])  delete v[i];
            if (H[i]) delete [] H[i];
         }
         delete [] sn;
         delete [] cs;
         delete [] s;
         delete [] H;
         return itmax;
      }

      j++;
   }

   for (int i=0; i<=m; ++i) {
      if (v[i])  delete v[i];
      if (H[i]) delete [] H[i];
   }
   delete [] sn;
   delete [] cs;
   delete [] s;
   delete [] H;

   return itmax;
}


void
FieldSolver::Update(LevelData<FArrayBox> &x, int k, double** h, double* s, int size_s, Vector<LevelData<FArrayBox> *> &v) const
{
   double* y = new double[size_s];
   for (int i=0; i<size_s; ++i) y[i] = s[i];

   // Backsolve:
   for (int i = k; i >= 0; i--) {
      y[i] /= h[i][i];
      for (int j = i - 1; j >= 0; j--)
         y[j] -= h[j][i] * y[i];
   }

   for (int j = 0; j <= k; j++) {
      for(DataIterator dit(x.dataIterator()); dit.ok(); ++dit) {
         x[dit].plus((*v[j])[dit], y[j]);
      }
   }

   delete [] y;
}


void
FieldSolver::GeneratePlaneRotation (double &dx, double &dy, double &cs, double &sn) const
{
   if (dy == 0.0) {
      cs = 1.0;
      sn = 0.0;
   } else if (fabs(dy) > fabs(dx)) {
      double temp = dx / dy;
      sn = 1.0 / sqrt( 1.0 + temp*temp );
      cs = temp * sn;
   } else {
      double temp = dy / dx;
      cs = 1.0 / sqrt( 1.0 + temp*temp );
      sn = temp * cs;
   }
}



void
FieldSolver::ApplyPlaneRotation (double &dx, double &dy, double &cs, double &sn) const
{
   double temp  =  cs * dx + sn * dy;
   dy = -sn * dx + cs * dy;
   dx = temp;
}



void
FieldSolver::fillInternalGhosts( LevelData<FArrayBox>& a_phi ) const
{
   a_phi.exchange();

   if (m_mblex_potential_Ptr && (!m_geometry.extrablockExchange())) {
     m_mblex_potential_Ptr->interpGhosts(a_phi);
   }
   else if (m_geometry.extrablockExchange()) {
     m_geometry.exchangeExtraBlockGhosts(a_phi);
   }
}



void
FieldSolver::setZero( LevelData<FArrayBox>& a_data ) const
{
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      a_data[dit].setVal(0.);
   }
}



void
FieldSolver::setZero( LevelData<FluxBox>& a_data ) const
{
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      a_data[dit].setVal(0.);
   }
}



double
FieldSolver::dotProduct( const LevelData<FArrayBox>& a_data1,
                         const LevelData<FArrayBox>& a_data2 ) const
{
   const DisjointBoxLayout& grids = a_data1.getBoxes();

   double local_sum = 0.;
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      local_sum += a_data1[dit].dotProduct(a_data2[dit],grids[dit]);
   }

   double global_sum;
#ifdef CH_MPI
   MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
   global_sum = local_sum;
#endif

   return global_sum;
}



void
FieldSolver::zeroAverage( LevelData<FArrayBox>& a_data ) const
{
   const DisjointBoxLayout & grids = a_data.disjointBoxLayout();

   double local_sum = 0.;
   double local_volume = 0.;

   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      const Box & box = grids[dit];
      FArrayBox temp(box,1);

      temp.copy(a_data[dit]);
      temp *= m_volume[dit];

      local_sum += temp.sum(box,0);
      local_volume += m_volume[dit].sum(box,0);
   }
   double sum;
   double volume;
#ifdef CH_MPI
   MPI_Allreduce(&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&local_volume, &volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
   sum = local_sum;
   volume = local_volume;
#endif
   double average = sum / volume;

   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      a_data[dit] -= average;
   }
}



double
FieldSolver::L2Norm( const LevelData<FArrayBox>& a_data ) const
{
   return sqrt(dotProduct(a_data, a_data));
}

#if 0
void
FieldSolver::computeElectricFieldExtrap( const LevelData<FArrayBox>& a_phi,
                                         const PotentialBC&          a_bc,
                                         LevelData<FArrayBox>&       a_E_field ) const
{
   CH_assert(a_E_field.nComp() == 3);

   const DisjointBoxLayout & grids = m_geometry.grids();

   LevelData<FluxBox> flux_face(grids, SpaceDim, a_E_field.ghostVect());
   computeFaceCenteredFieldExtrap(a_phi, a_bc, flux_face);

   LevelData<FArrayBox> tmp(grids, SpaceDim*SpaceDim, a_E_field.ghostVect());

   // Need to make this fourth-order
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      EdgeToCell(flux_face[dit], tmp[dit]);
   }

   // Average the directions
   LevelData<FArrayBox> flux_cell(grids, SpaceDim, a_E_field.ghostVect());

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      FArrayBox& this_tmp = tmp[dit];
      FArrayBox& this_flux = flux_cell[dit];

      this_flux.setVal(0.);
      for (int comp=0; comp<SpaceDim; comp++) {
        for (int dir=0; dir<SpaceDim; dir++) {
          this_flux.plus(this_tmp, comp*SpaceDim + dir, comp, 1);
        }
      }
      this_flux.mult(1./SpaceDim);
   }

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      FArrayBox& this_E_field = a_E_field[dit];
      FArrayBox& this_flux = flux_cell[dit];

      this_E_field.copy(this_flux, 0, 0, 1);
      this_E_field.setVal(0., 1);
      this_E_field.copy(this_flux, 1, 2, 1);
   }
}



void
FieldSolver::computeElectricFieldExtrap( const LevelData<FArrayBox>& a_phi,
                                         const PotentialBC&          a_bc,
                                         LevelData<FluxBox>&         a_E_field ) const
{
   CH_assert(a_E_field.nComp() == 3);

   LevelData<FluxBox> flux(a_E_field.disjointBoxLayout(), SpaceDim, a_E_field.ghostVect());
   computeFaceCenteredFieldExtrap(a_phi, a_bc, flux);

   DataIterator dit = a_E_field.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FluxBox& this_E_field = a_E_field[dit];
      FluxBox& this_flux = flux[dit];
      for (int dir=0; dir<SpaceDim; ++dir) {
         FArrayBox& this_E_field_dir = this_E_field[dir];
         FArrayBox& this_flux_dir = this_flux[dir];
         this_E_field_dir.copy(this_flux_dir, 0, 0, 1);
         this_E_field_dir.setVal(0., 1);
         this_E_field_dir.copy(this_flux_dir, 1, 2, 1);
      }
   }
}



void
FieldSolver::computeFaceCenteredFieldExtrap( const LevelData<FArrayBox>& a_phi,
                                             const PotentialBC&          a_bc,
                                             LevelData<FluxBox>&         a_field,
                                             const bool                  a_homogeneousBCs ) const
{
   /*
     Computes the face-centered (i.e., pointwise) field:

          a_field = - NJinverse grad a_phi

     It is assumed that the ghost cells of a_phi have already been filled.        

   */

   LevelData<FluxBox> mapped_field(m_geometry.grids(), SpaceDim, IntVect::Unit);
   computeFaceCenteredMappedFieldExtrap(a_phi, a_bc, mapped_field, a_homogeneousBCs);

   // Multiply mapped field by NJinverse to put it in the physical frame
   LevelData<FluxBox> NJinverse(m_geometry.grids(), SpaceDim*SpaceDim, IntVect::Unit);
   m_geometry.getPointwiseNJinverse(NJinverse);

   for (DataIterator dit(a_field.dataIterator()); dit.ok(); ++dit) {

      for (int dir=0; dir<SpaceDim; ++dir) {
         const FArrayBox& this_NJInverse = NJinverse[dit][dir];
         const FArrayBox& this_mapped_field = mapped_field[dit][dir];
         FArrayBox& this_field = a_field[dit][dir];
         FArrayBox tmp(this_field.box(), 1);

         this_field.setVal(0.);
         for (int row=0; row<SpaceDim; ++row) {
            for (int col=0; col<SpaceDim; ++col) {
               tmp.copy(this_NJInverse,SpaceDim*row+col,0,1);
               tmp.mult(this_mapped_field,col,0,1);
               this_field.plus(tmp,0,row,1);
            }
         }
      }
   }
}



void
FieldSolver::computeFaceCenteredMappedFieldExtrap( const LevelData<FArrayBox>& a_phi,
                                                   const PotentialBC&          a_bc,
                                                   LevelData<FluxBox>&         a_field,
                                                   const bool                  a_homogeneousBCs ) const
{
   /*
     Computes the face-centered (i.e., pointwise) field in mapped coordinates:

          a_field = - grad a_phi

     It is assumed that the ghost cells of a_phi have already been filled.        

   */

   const MagCoordSys* mag_coord_sys = m_geometry.getCoordSys();
   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = mag_coord_sys->boundaries();
   const DisjointBoxLayout& grids = m_geometry.grids();

   // Compute second-order for now
   //   bool fourthOrder = true;
   //   bool dropOrder = true;
   bool fourthOrder = false;
   bool dropOrder = false;

   IntVect stencil_box_lo(IntVect::Zero);
   IntVect stencil_box_hi;
   if( fourthOrder ) {
      stencil_box_hi = 4*IntVect::Unit;
   }
   else {
      stencil_box_hi = 2*IntVect::Unit;
   }
   Box stencil_box(stencil_box_lo, stencil_box_hi);
          
   FArrayBox stencil_values(stencil_box, 1);

   int radius = (stencil_box.size(0)-1)/2;

   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
     const Box & box = grids[dit];
     int block_number = mag_coord_sys->whichBlock(box);
     const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(block_number);
     const Box& domainBox = block_coord_sys.domain().domainBox();
     RealVect dx = block_coord_sys.dx();

     const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];
     bool lo_block_interface[SpaceDim];
     bool hi_block_interface[SpaceDim];
     for (int dir=0; dir<SpaceDim; ++dir) {
        lo_block_interface[dir] = this_block_boundaries[dir].isInterface();
        hi_block_interface[dir] = this_block_boundaries[dir + SpaceDim].isInterface();
     }

     int block_bc_types[2*SpaceDim];
     for (int dir=0; dir<SpaceDim; ++dir) {
        for (int side=0; side<2; ++side) {
           block_bc_types[dir+2*side] = a_bc.getBCType(block_number, dir, side);
        }
     }

     const FArrayBox& this_phi = a_phi[dit];

     for (int dir=0; dir<SpaceDim; ++dir) {
        FArrayBox& this_flux = a_field[dit][dir];

        Box template_box(IntVect(D_DECL6(-radius,-radius,-radius,-radius,-radius,-radius)),
                         IntVect(D_DECL6(radius,radius,radius,radius,radius,radius)));
        template_box.growHi(dir, -1);

        BoxIterator bit(surroundingNodes(box,dir));
        for (bit.begin();bit.ok();++bit) {
           IntVect iv_face = bit();

           IntVect iv(iv_face);
           int side;
           if(box.contains(iv)) {
              side = 0;
           }
           else {
              side = 1;
              iv.shift(dir,-1);
              CH_assert(box.contains(iv));
           }

           Box phi_box = template_box + iv_face;

           for( int s = 0; s < SpaceDim; ++s ) {
              stencil_values.setVal(0.);
              accumFieldStencil(iv, dir, side, s, this_block_boundaries, block_bc_types,
                                dx, domainBox, fourthOrder, dropOrder, stencil_values);

              this_flux(iv_face,s) = 0.;

              BoxIterator bit2(phi_box);
              for (bit2.begin();bit2.ok();++bit2) {
                 IntVect iv2 = bit2();
                 this_flux(iv_face,s) += stencil_values(iv2 - iv + radius*IntVect::Unit) * this_phi(iv2);
              }
           }
        }
     }
   }

   if ( !a_homogeneousBCs) {
      setMappedBoundaryFieldExtrap(a_field);
   }

   a_field.exchange();
}


void
FieldSolver::setMappedBoundaryFieldExtrap( LevelData<FluxBox>& a_field ) const
{
   // N.B.: This is only second-order for now

   const MagCoordSys* mag_coord_sys = m_geometry.getCoordSys();
   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = mag_coord_sys->boundaries();
   const DisjointBoxLayout& grids = m_geometry.grids();

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
     int block_number = mag_coord_sys->whichBlock(grids[dit]);
     const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(block_number);
     const Box& domainBox = block_coord_sys.domain().domainBox();
     const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];
     const RealVect& dx = block_coord_sys.dx();

     FluxBox& this_field = a_field[dit];

     Box domain_box_lo_radial = adjCellLo(domainBox,RADIAL_DIR,-1);
     Box domain_box_hi_radial = adjCellHi(domainBox,RADIAL_DIR,-1);

     for (int dir=0; dir<SpaceDim; dir++) {
        FArrayBox& this_field_dir = this_field[dir];

        for (int radial_side=0; radial_side<2; ++radial_side) {

           Box boundary_box;
           if (radial_side == 0) {
              if (dir == RADIAL_DIR) {
                 boundary_box = domain_box_lo_radial;
                 boundary_box.shiftHalf(dir,-1);
              }
              else {
                 boundary_box = surroundingNodes(domain_box_lo_radial,dir);
              }
           }
           else {
              if (dir == RADIAL_DIR) {
                 boundary_box = domain_box_hi_radial;
                 boundary_box.shiftHalf(dir,1);
              }
              else {
                 boundary_box = surroundingNodes(domain_box_hi_radial,dir);
              }
           }

           Box overlap = boundary_box & this_field_dir.box();
           if (overlap.ok()) {

              for (int comp=0; comp<SpaceDim; ++comp) {

                 FORT_SET_EXTRAP_RADIAL_FIELD_BV2( CHF_BOX(overlap),
                                                   CHF_CONST_INT(dir),
                                                   CHF_CONST_INT(comp),
                                                   CHF_CONST_INT(radial_side),
                                                   CHF_CONST_REALVECT(dx),
                                                   CHF_FRA1(this_field_dir,comp) );
              }
           }
        }
     }
   }

   a_field.exchange();
}

#endif


void
FieldSolver::computeElectricField( const LevelData<FArrayBox>& a_phi,
                                   LevelData<FArrayBox>&       a_E_field ) const
{
   CH_assert(a_E_field.nComp() == 3);
   LevelData<FArrayBox> mapped_E_field(a_E_field.disjointBoxLayout(), 3, a_E_field.ghostVect());

   computeMappedElectricField( a_phi, mapped_E_field );
   getUnmappedField(mapped_E_field, a_E_field);
}



void
FieldSolver::computeElectricField( const LevelData<FArrayBox>& a_phi,
                                   LevelData<FluxBox>&         a_E_field ) const
{
   CH_assert(a_E_field.nComp() == 3);
   LevelData<FluxBox> mapped_E_field(a_E_field.disjointBoxLayout(), 3, a_E_field.ghostVect());

   computeMappedElectricField( a_phi, mapped_E_field );
   getUnmappedField(mapped_E_field, a_E_field);
}


void
FieldSolver::getUnmappedField( const LevelData<FArrayBox>& a_mapped_E_field,
                               LevelData<FArrayBox>&       a_E_field ) const
{
   const DisjointBoxLayout & grids = m_geometry.grids();

   // Multiply by NJInverse to get the unmapped field
   DataIterator dit = a_mapped_E_field.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);

      Box box = a_mapped_E_field[dit].box();
      FArrayBox NJInverse(box, SpaceDim*SpaceDim);
      block_coord_sys.getPointwiseNJInverse(NJInverse);

      FORT_MULT_NJINVERSE(CHF_BOX(box),
                          CHF_CONST_FRA(a_mapped_E_field[dit]),
                          CHF_CONST_FRA(NJInverse),
                          CHF_FRA(a_E_field[dit]));
   }
}



void
FieldSolver::getUnmappedField(const LevelData<FluxBox>& a_mapped_E_field,
                              LevelData<FluxBox>&       a_E_field) const
{
   const DisjointBoxLayout & grids = m_geometry.grids();

   // Multiply by NJInverse to get the unmapped field
   DataIterator dit = a_mapped_E_field.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);

      FluxBox NJInverse(a_mapped_E_field[dit].box(), SpaceDim*SpaceDim);
      block_coord_sys.getPointwiseNJInverse(NJInverse);

      for (int dir=0; dir<SpaceDim; dir++) {
         FORT_MULT_NJINVERSE(CHF_BOX(a_mapped_E_field[dit][dir].box()),
                             CHF_CONST_FRA(a_mapped_E_field[dit][dir]),
                             CHF_CONST_FRA(NJInverse[dir]),
                             CHF_FRA(a_E_field[dit][dir]));
      }
   }
}



void
FieldSolver::computeMappedElectricField( const LevelData<FArrayBox>& a_phi,
                                         LevelData<FArrayBox>&       a_E_field ) const
{
   /*
     This function computes the field

          E = - grad phi

     to fourth-order accuracy in valid cells.  It also computes second-order accurate
     values in one layer of interior/extrablock ghost cells and one layer of ghost cells
     at the domain boundary.

     Field values not set to valid values are deliberately set to NaNs to help
     detect their inadvertent use.

   */

   CH_assert(a_E_field.ghostVect() == IntVect::Unit);
   CH_assert(a_E_field.nComp() == 3);
   CH_assert(m_num_potential_ghosts_interpolated >= 2);

   const DisjointBoxLayout & grids = m_geometry.grids();
   LevelData<FArrayBox> phi_cell(grids, 1, m_num_potential_ghosts_interpolated*IntVect::Unit);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      phi_cell[dit].copy(a_phi[dit]);
   }

   fillInternalGhosts(phi_cell);

   // Extrapolate the potential to two ghost cells at physical boundaries

   const MagCoordSys* mag_coord_sys = m_geometry.getCoordSys();
   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries = mag_coord_sys->boundaries();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = mag_coord_sys->whichBlock(grids[dit]);
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(block_number);
      const ProblemDomain& domain = block_coord_sys.domain();
      for (int dir=0; dir<SpaceDim; ++dir) {
         if ( !domain.isPeriodic(dir) ) {
            IntVect grow_vec = 2*IntVect::Unit;
            grow_vec[dir] = 0;
            Box interior_box = grow(grids[dit],grow_vec);
            Box domain_box = grow(domain.domainBox(),grow_vec);

            extrapBoundaryGhostsForCC(phi_cell[dit], interior_box, domain_box, dir, boundaries[block_number]);
         }
      }
   }

   // We now have potential values in two ghost cell layers.  We next compute the field to
   // second-order at valid cell centers and one ghost cell layer.

   int order = 2;

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      RealVect dx = block_coord_sys.dx();
      Box box = grow(grids[dit],1);

      // Initializing to NaN to ensure that the field is actually computed everywhere it's used
      a_E_field[dit].setVal(1./0.);

      for (int dir=0; dir<SpaceDim; dir++) {
         FORT_CELL_CENTERED_FIELD_COMPONENT(CHF_BOX(box),
                                            CHF_CONST_INT(dir),
                                            CHF_CONST_FRA1(phi_cell[dit],0),
                                            CHF_CONST_REALVECT(dx),
                                            CHF_CONST_INT(order),
                                            CHF_FRA1(a_E_field[dit],(2*dir)%3));
      }
   }

   a_E_field.exchange();

   // Recompute the field at valid cell centers to fourth-order

   order = 4;

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      RealVect dx = block_coord_sys.dx();

      for (int dir=0; dir<SpaceDim; ++dir) {
         FORT_CELL_CENTERED_FIELD_COMPONENT(CHF_BOX(grids[dit]),
                                            CHF_CONST_INT(dir),
                                            CHF_CONST_FRA1(phi_cell[dit],0),
                                            CHF_CONST_REALVECT(dx),
                                            CHF_CONST_INT(order),
                                            CHF_FRA1(a_E_field[dit],(2*dir)%3));
      }
   }
}



void
FieldSolver::computeMappedElectricField( const LevelData<FArrayBox>& a_phi,
                                         LevelData<FluxBox>&         a_E_field ) const
{
   /*
     This function computes the field

          E = - grad phi

     to fourth-order accuracy at valid cell faces.  It also computes second-order
     accurate values on transverse faces at all block boundaries.

     Since not all entries of a_E_field are set, those not set are returned
     as NaN to help detect any inadvertent attempt to access them.

   */

   CH_assert(a_E_field.ghostVect() == IntVect::Unit);
   CH_assert(a_E_field.nComp() == 3);
   CH_assert(m_num_potential_ghosts_interpolated >= 2);

   const DisjointBoxLayout & grids = m_geometry.grids();

   LevelData<FArrayBox> phi_cell(grids, 1, m_num_potential_ghosts_interpolated*IntVect::Unit);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      phi_cell[dit].copy(a_phi[dit]);
   }

   fillInternalGhosts(phi_cell);

   // Fourth-order extrapolate the potential to two cells at physical boundaries

   int order = 4;

   const MagCoordSys* mag_coord_sys = m_geometry.getCoordSys();
   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries = mag_coord_sys->boundaries();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = mag_coord_sys->whichBlock(grids[dit]);
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(block_number);
      const ProblemDomain& domain = block_coord_sys.domain();
      for (int dir=0; dir<SpaceDim; ++dir) {
         if ( !domain.isPeriodic(dir) ) {
            IntVect grow_vec = phi_cell.ghostVect();
            grow_vec[dir] = 0;
            Box interior_box = grow(grids[dit],grow_vec);
            Box domain_box = grow(domain.domainBox(),grow_vec);

            extrapBoundaryGhostsForFC(phi_cell[dit], interior_box, domain_box, dir, order, boundaries[block_number]);
         }
      }
   }

   // We now have potential values in two extrablock ghost cell layers and two physical ghost cell layers.
   // We next compute the potential to fourth-order at valid cell faces and two layers of transverse cell faces.

   LevelData<FluxBox> phi_face(grids, 1, 2*IntVect::Unit);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         IntVect grow_vect = 2*IntVect::Unit;
         grow_vect[dir] = 0;
         Box box = grow(grids[dit],grow_vect);

         FORT_FACE_INTERPOLATE(CHF_CONST_INT(dir),
                               CHF_BOX(surroundingNodes(box,dir)),
                               CHF_CONST_INT(order),
                               CHF_CONST_FRA1(phi_cell[dit],0),
                               CHF_FRA1(phi_face[dit][dir],0));
      }
   }
   phi_face.exchange();

   // Compute the field to second-order including one layer of transverse faces

   order = 2;

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      RealVect dx = block_coord_sys.dx();
      Box box = grow(grids[dit],1);

      // Initializing to NaN to ensure that the field is actually computed everywhere it's used
      a_E_field[dit].setVal(1./0.);

      for (int dir=0; dir<SpaceDim; dir++) {
         Box box_dir = surroundingNodes(box,dir);
         FORT_FACE_CENTERED_FIELD_COMPONENT(CHF_BOX(box_dir),
                                            CHF_CONST_INT(dir),
                                            CHF_CONST_FRA1(phi_cell[dit],0),
                                            CHF_CONST_REALVECT(dx),
                                            CHF_CONST_INT(order),
                                            CHF_FRA1(a_E_field[dit][dir],(2*dir)%3));

         for (int tdir=0; tdir<SpaceDim; ++tdir) {
            if (tdir != dir) {
               FORT_CELL_CENTERED_FIELD_COMPONENT(CHF_BOX(box_dir),
                                                  CHF_CONST_INT(tdir),
                                                  CHF_CONST_FRA1(phi_face[dit][dir],0),
                                                  CHF_CONST_REALVECT(dx),
                                                  CHF_CONST_INT(order),
                                                  CHF_FRA1(a_E_field[dit][dir],(2*tdir)%3));
            }
         }
      }
   }
   a_E_field.exchange();

   // Second-order extrapolate the field to the transverse physical boundary ghosts
   m_geometry.fillTransversePhysicalGhosts(a_E_field);

   // Recompute the field to fourth-order on valid cell faces

   order = 4;

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      RealVect dx = block_coord_sys.dx();

      for (int dir=0; dir<SpaceDim; dir++) {
         Box box_dir = surroundingNodes(grids[dit],dir);
         FORT_FACE_CENTERED_FIELD_COMPONENT(CHF_BOX(box_dir),
                                            CHF_CONST_INT(dir),
                                            CHF_CONST_FRA1(phi_cell[dit],0),
                                            CHF_CONST_REALVECT(dx),
                                            CHF_CONST_INT(order),
                                            CHF_FRA1(a_E_field[dit][dir],(2*dir)%3));

         for (int tdir=0; tdir<SpaceDim; ++tdir) {
            if (tdir != dir) {
               FORT_CELL_CENTERED_FIELD_COMPONENT(CHF_BOX(box_dir),
                                                  CHF_CONST_INT(tdir),
                                                  CHF_CONST_FRA1(phi_face[dit][dir],0),
                                                  CHF_CONST_REALVECT(dx),
                                                  CHF_CONST_INT(order),
                                                  CHF_FRA1(a_E_field[dit][dir],(2*tdir)%3));
            }
         }
      }
   }

   a_E_field.exchange();
}



void
FieldSolver::extrapBoundaryGhostsForCC(FArrayBox&                              a_data,
                                       const Box&                              a_interiorbox,
                                       const Box&                              a_domain_box,
                                       const int                               a_dir,
                                       const Tuple<BlockBoundary, 2*SpaceDim>& a_block_boundaries) const
{
   // This function fourth-order extrapolates to fill two layers of a_data ghost cells
   // at physical boundaries in the direction a_dir.

   const Box& bx= a_data.box();
   CH_assert(bx.contains(a_interiorbox));

   for (int side=-1; side<2; side+=2) {

      if ( a_block_boundaries[a_dir + (side+1)*SpaceDim/2].isDomainBoundary() ) {

         Box dstbox;
         bool isbc = false;
      
         switch(side) 
            {
            case -1:
               // Handle the low side
               isbc = (a_interiorbox.smallEnd(a_dir) == a_domain_box.smallEnd(a_dir));
               if (isbc) {
                  dstbox = adjCellLo(a_interiorbox, a_dir, 2);
               }
               break;
            case 1:
               // Handle high side
               isbc = (a_interiorbox.bigEnd(a_dir) == a_domain_box.bigEnd(a_dir));
               if (isbc) {
                  dstbox = adjCellHi(a_interiorbox, a_dir, 2);
               }
               break;
            }

         if (isbc) {
            CH_assert(a_interiorbox.size(a_dir)>=5);
            CH_assert(bx.contains(dstbox));
            FORT_EXTRAP_FOR_CC_OPS(CHF_CONST_INT(a_dir),
                                   CHF_CONST_INT(side),
                                   CHF_BOX(dstbox),
                                   CHF_BOX(a_interiorbox),
                                   CHF_FRA(a_data));
         }
      }
   }
} 


void
FieldSolver::extrapBoundaryGhostsForFC(FArrayBox&                              a_data,
                                       const Box&                              a_interiorbox,
                                       const Box&                              a_domain_box,
                                       const int                               a_dir,
                                       const int                               a_order,
                                       const Tuple<BlockBoundary, 2*SpaceDim>& a_block_boundaries) const
{
   // This function fourth-order extrapolates to fill two layers of a_data ghost cells
   // at physical boundaries in the direction a_dir.

   const Box& bx= a_data.box();
   CH_assert(bx.contains(a_interiorbox));

   for (int side=-1; side<2; side+=2) {

      if ( a_block_boundaries[a_dir + (side+1)*SpaceDim/2].isDomainBoundary() ) {

         Box dstbox;
         bool isbc = false;
      
         switch(side) 
            {
            case -1:
               // Handle the low side
               isbc = (a_interiorbox.smallEnd(a_dir) == a_domain_box.smallEnd(a_dir));
               if (isbc) {
                  dstbox = adjCellLo(a_interiorbox, a_dir, 2);
               }
               break;
            case 1:
               // Handle high side
               isbc = (a_interiorbox.bigEnd(a_dir) == a_domain_box.bigEnd(a_dir));
               if (isbc) {
                  dstbox = adjCellHi(a_interiorbox, a_dir, 2);
               }
               break;
            }

         if (isbc) {
            CH_assert( (a_order==4 && a_interiorbox.size(a_dir)>=5)
                    || (a_order==2 && a_interiorbox.size(a_dir)>=3));
            CH_assert(bx.contains(dstbox));
            FORT_EXTRAP_FOR_FC_OPS(CHF_CONST_INT(a_dir),
                                   CHF_CONST_INT(side),
                                   CHF_CONST_INT(a_order),
                                   CHF_BOX(dstbox),
                                   CHF_BOX(a_interiorbox),
                                   CHF_FRA(a_data));
         }
      }
   }
} 


void
FieldSolver::extractNormalComponent( const LevelData<FluxBox>& a_in,
                                     const int                 a_dir,
                                     LevelData<FluxBox>&       a_out ) const
{
   CH_assert(a_out.nComp() == 1);
   CH_assert(a_in.nComp() >= a_dir);

   for (DataIterator dit(a_out.dataIterator()); dit.ok(); ++dit) {
      FluxBox& this_fb = a_out[dit];
      this_fb.setVal(0.);
      this_fb[a_dir].copy(a_in[dit][a_dir],a_dir,0,1);
   }
}



double
FieldSolver::globalMax(const double a_data) const
{
   double global_max;

#ifdef CH_MPI
   double local_data = a_data;
   MPI_Allreduce(&local_data, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
   global_max = a_data;
#endif

   return global_max;
}


#include "NamespaceFooter.H"
