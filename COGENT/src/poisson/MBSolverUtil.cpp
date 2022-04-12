#include "MBSolverUtil.H"
#include "MBSolverUtilF_F.H"
#include "FourthOrderUtil.H"
#include "NewCoordSys.H"

#include "NamespaceHeader.H"



MBSolverUtil::MBSolverUtil( const MBSolver& a_mb_solver )
   : m_geometry(a_mb_solver.m_geometry),
     m_coord_sys_ptr(m_geometry.coordSysPtr()),
     m_discretization_order(a_mb_solver.m_discretization_order)
{
   CH_assert(m_discretization_order == 2 || m_discretization_order == 4);
}
      


MBSolverUtil::~MBSolverUtil()
{
}



void
MBSolverUtil::computeMappedCoefficients( const LevelData<FluxBox>&  a_unmapped_coefficients,
                                         LevelData<FluxBox>&        a_mapped_coefficients ) const
{
   int ncomp = a_unmapped_coefficients.nComp();
   IntVect ghosts = a_unmapped_coefficients.ghostVect();

   const DisjointBoxLayout& grids = m_geometry.grids();
   LevelData<FluxBox> N(grids, SpaceDim*SpaceDim, ghosts);
   getPointwiseN(N);

   LevelData<FluxBox> Jinverse(grids, 1, ghosts);
   getPointwiseJinverse(Jinverse);

   LevelData<FluxBox> unmapped_coefficients(grids, ncomp, ghosts);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      unmapped_coefficients[dit].copy(a_unmapped_coefficients[dit]);
   }   

   if (m_discretization_order == 4) {
      fillTransverseGhosts(unmapped_coefficients);
   }

   a_mapped_coefficients.define(grids, ncomp, ghosts);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const FluxBox& this_unmapped_coef = unmapped_coefficients[dit];
      FluxBox& this_mapped_coef = a_mapped_coefficients[dit];
      FluxBox& this_N = N[dit];
      FluxBox& this_Jinverse = Jinverse[dit];
         
      for (int dir=0; dir<SpaceDim; ++dir) {
         FArrayBox& this_mapped_coef_dir = this_mapped_coef[dir];
         const Box& box_dir = this_mapped_coef_dir.box();
         
         FORT_COMPUTE_MAPPED_COEFFICIENTS(CHF_BOX(box_dir),
                                          CHF_CONST_FRA(this_unmapped_coef[dir]),
                                          CHF_CONST_FRA(this_N[dir]),
                                          CHF_CONST_FRA1(this_Jinverse[dir],0),
                                          CHF_FRA(this_mapped_coef_dir));
      }
   }

   if (m_discretization_order == 4) {
      // Convert the mapped coefficients from face-centered to face-averaged
      fourthOrderAverage(a_mapped_coefficients);
   }
}



void
MBSolverUtil::fillTransverseGhosts( LevelData<FluxBox>& a_data ) const
{
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const ProblemDomain& block_domain = m_coord_sys_ptr->blockDomainOfBox(grids[dit]);

      for (int dir=0; dir<SpaceDim; ++dir) {
         Box interior = surroundingNodes(grids[dit],dir);

         secondOrderTransExtrapAtDomainBdry(a_data[dit][dir],
                                            dir,
                                            interior,
                                            block_domain);
      }
   }

   a_data.exchange();
}


void
MBSolverUtil::getPointwiseN( LevelData<FluxBox>& a_N ) const
{
   const DisjointBoxLayout& grids = a_N.disjointBoxLayout();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const NewCoordSys* coord_sys = m_coord_sys_ptr->getCoordSys(grids[dit]);
      FluxBox& this_N = a_N[dit];
      for (int dir=0; dir<SpaceDim; ++dir) {
         FArrayBox& this_N_dir = this_N[dir];

         const Box& box = this_N_dir.box();
         FArrayBox xi(box, SpaceDim);
         coord_sys->getCenterMappedCoordinates(xi, box);

#if CH_SPACEDIM == 2
         coord_sys->dXdXi(this_N_dir, xi, 0, 1, 1, box);
         coord_sys->dXdXi(this_N_dir, xi, 1, 1, 0, box);
         coord_sys->dXdXi(this_N_dir, xi, 2, 0, 1, box);
         coord_sys->dXdXi(this_N_dir, xi, 3, 0, 0, box);
         this_N_dir.negate(1,2);
#endif
#if CH_SPACEDIM == 3
         int m = 0;
         for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
               coord_sys->dXdXi(this_N_dir, xi, m++, i, j, box);
            }
         }

         double a[3][3], b[3][3];
         for (BoxIterator bit(box); bit.ok(); ++bit) {
            IntVect iv = bit();
            m = 0;
            for (int i=0; i<3; ++i) {
               for (int j=0; j<3; ++j) {
                  a[i][j] = this_N_dir(iv,m++);
               }
            }

            for(int i=0; i<3; ++i) {
               for(int j=0; j<3; j++) {
                  b[i][j] = (a[(i+1)%3][(j+1)%3] * a[(i+2)%3][(j+2)%3])
                     - (a[(i+1)%3][(j+2)%3] * a[(i+2)%3][(j+1)%3]);
               }
            }

            m = 0;
            for (int i=0; i<3; ++i) {
               for (int j=0; j<3; ++j) {
                  this_N_dir(iv,m++) = b[i][j];
               }
            }
         }
#endif
      }
   }
}




void
MBSolverUtil::getPointwiseJinverse( LevelData<FluxBox>& a_Jinverse ) const
{
   const DisjointBoxLayout& grids = a_Jinverse.disjointBoxLayout();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const NewCoordSys* coord_sys = m_coord_sys_ptr->getCoordSys(grids[dit]);
      FluxBox& this_Jinverse = a_Jinverse[dit];
      for (int dir=0; dir<SpaceDim; ++dir) {
         FArrayBox& this_Jinverse_dir = this_Jinverse[dir];

         const Box& box = this_Jinverse_dir.box();
         FArrayBox xi(box, SpaceDim);
         coord_sys->getCenterMappedCoordinates(xi, box);

         coord_sys->pointwiseJ(this_Jinverse_dir, xi, box);
         this_Jinverse_dir.invert(1.);
      }
   }
}


void
MBSolverUtil::getCellVolumes( LevelData<FArrayBox>& a_volume ) const
{
   const DisjointBoxLayout& grids = a_volume.disjointBoxLayout();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const NewCoordSys* coord_sys = m_coord_sys_ptr->getCoordSys(grids[dit]);
      FArrayBox& this_volume = a_volume[dit];

      const Box& box = this_volume.box();

      // The flux needs two extra ghost cells, one to accomodate the
      // loss of a ghost cell in the fourth-order multiplication by N^T
      // that happens in computeDivergence and one to accomodate the
      // loss of a ghost cell in fourthOrderAverage.
      Box fbox(box);
      fbox.grow(2);
      FluxBox F(fbox, SpaceDim);

      FluxBox N(fbox, coord_sys->getNumN());
      coord_sys->getN(N, fbox);

      F.setVal(0.0);
      for (int dir=0; dir<SpaceDim; dir++) {
         FArrayBox& thisFdir = F[dir];
         const Box& this_box_dir = thisFdir.box();

         FArrayBox xi(this_box_dir, SpaceDim);
         coord_sys->getCenterMappedCoordinates(xi, this_box_dir);

         BoxIterator bit(this_box_dir);
         // this is going to be slow, but we can
         // eventually move this into fortran
         for (bit.begin(); bit.ok(); ++bit) {
            IntVect iv = bit();
            RealVect mappedLoc;
            for (int n=0; n<SpaceDim; ++n) {
               mappedLoc[n] = xi(iv,n);
            }
            RealVect realLoc = coord_sys->realCoord(mappedLoc);

            for (int n=0; n<SpaceDim; ++n) {
               thisFdir(iv,n) = realLoc[n];
            }
         }
      } // end loop over directions

      // convert point values to 4th-order face averages
      fourthOrderAverageFace(F);

      Interval divInterval = this_volume.interval();
      coord_sys->computeDivergence(this_volume, F, N, box, divInterval);

      this_volume /= SpaceDim;
   }
}



#include "NamespaceFooter.H"
