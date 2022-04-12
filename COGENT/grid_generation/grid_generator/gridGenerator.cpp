#include "BoxIterator.H"

#include "BlockMapping.H"
#include "GridOutput.H"

#include "CONSTANTS.H"

#include "CH_HDF5.H"
#include "parstream.H"
#ifdef CH_MPI
#include "CH_Attach.H"
#endif


#include "UsingNamespace.H"
enum SingleNullBlockType {LCORE,RCORE,LCSOL,RCSOL,LSOL,RSOL,LPF,RPF,MCORE,MCSOL,NUM_SINGLE_NULL_BLOCKS};
enum directions {RADIAL_DIR, POLOIDAL_DIR};

void
getUnitSquareGrid(const int  a_n_radial,
                  const int  a_n_poloidal,
                  FArrayBox& a_xi)
{
   Box box(IntVect::Zero, IntVect(a_n_radial,a_n_poloidal));
   RealVect dx(1./a_n_radial,1./a_n_poloidal);
   
   a_xi.define(box,2);
   
   for (BoxIterator bit(box);bit.ok();++bit) {
      IntVect iv = bit();
      for (int dir=0; dir<SpaceDim; ++dir) {
         a_xi(iv,dir) = iv[dir]*dx[dir];
      }
   }
}


int main( int a_argc, char* a_argv[] )
{
   
#ifdef CH_MPI
   // Start MPI
   MPI_Init( &a_argc, &a_argv );
   setChomboMPIErrorHandler();
#endif

   if (a_argc!=24) {
      cout << "Twenty-four parameter input is required" << endl;
      cout << "Refer to README file for usage" << endl;
      return -1;
   }
   
   string DCT_file_name    = a_argv[1];          //name of the DCT file
   string output_file_name = a_argv[2];          //name of the output file
   double psiNormInner     = atof(a_argv[3]);    //normalized flux at the inner core boundary
   double psiNormOuter     = atof(a_argv[4]);    //normalized flux at the outer sol boundary

   int n_psiCORE           = atoi(a_argv[5]);    //number of radial cells in the core region
   int n_thetaCORE         = atoi(a_argv[6]);    //number of poloidal cells in the core region
   int n_thetaInnerLeg     = atoi(a_argv[7]);    //number of poloidal cells in the inner leg
   int n_thetaOuterLeg     = atoi(a_argv[8]);    //number of poloidal cells in the outer leg
   int n_psiPF             = atoi(a_argv[9]);    //number of radial cells in the PF region
   int n_thetaMCORE        = atoi(a_argv[10]);   //number of poloidal cells in the mcore region

   int n_extrpInnerLeg     = atoi(a_argv[11]);   //number of poloidal extrapolation cells in the inner leg
   int n_extrpOuterLeg     = atoi(a_argv[12]);   //number of poloidal extrapolation cells in the outer leg

   int n_psiGhost          = atoi(a_argv[13]);   //number of radial ghost cells
   int n_thetaGhost        = atoi(a_argv[14]);   //number of poloidal ghost cells

   double trans_rad        = atof(a_argv[15]);   //transition radius (set to a negative value for flux-aligned grids)
   double trans_length     = atof(a_argv[16]);   //transition length = distance from X point beyond which grid is locally orthogonal
   double blending_factor  = atof(a_argv[17]);   //blending factor weighting the real and block-aligned flux
   
   RealVect init_guessX;
   RealVect init_guessO;
   
   init_guessX[0]    = atof(a_argv[18]);         //initial (r) guess for Xpt
   init_guessX[1]    = atof(a_argv[19]);         //initial (z) guess for Xpt

   init_guessO[0]    = atof(a_argv[20]);         //initial (r) guess for Opt
   init_guessO[1]    = atof(a_argv[21]);         //initial (z) guess for Opt
   
   bool compute_crit_pts = atoi(a_argv[22]);  //should we compute critical points (if not, rely on the data in DCT file).
   bool output_gridUe = atoi(a_argv[23]);        //should we output data for gridUE generator?
   
   if (n_thetaMCORE > 0) {
      CH_assert(n_thetaMCORE < n_thetaCORE);
   }
   
   /*
    Step 1: Compute critical points and poloidal arc-length mapping in all blocks
    */
   RealVect Xpoint;
   RealVect Opoint;

   BlockMapping block_mapping_tmp(DCT_file_name, LCORE, 0, 0, 0., 0.);
   if (compute_crit_pts) {
      Xpoint = block_mapping_tmp.findCriticalPoint(init_guessX, 0.01);
      Opoint = block_mapping_tmp.findCriticalPoint(init_guessO, 0.01);
      block_mapping_tmp.m_Xpoint  = Xpoint;
      block_mapping_tmp.m_magAxis = Opoint;
   }
   block_mapping_tmp.init(trans_rad);
   
   //Get core (poloidal) arc-length mapping
   Vector<RealVect*> mappingCORE;
   double lengthCORE(0.);
   block_mapping_tmp.getArcLengthMapping(mappingCORE, lengthCORE, LCORE);
   
   //Get outer leg mapping (default length is lengthCORE/4.0)
   block_mapping_tmp.m_block_number = RPF;
   double outer_div_leg_length = lengthCORE/2.0;
   Vector<RealVect*> mappingOuterLeg;
   block_mapping_tmp.getArcLengthMapping(mappingOuterLeg, outer_div_leg_length, RPF);

   //Get inner leg mapping (default length is lengthCORE/4.0)
   block_mapping_tmp.m_block_number = LPF;
   double inner_div_leg_length = lengthCORE/2.0;
   Vector<RealVect*> mappingInnerLeg;
   block_mapping_tmp.getArcLengthMapping(mappingInnerLeg, inner_div_leg_length, LPF);
   
   if (procID()==0) {
      cout<<"Done constructing arc-length mapping"<<endl;
   }

   
   /*
    Step 2: Loop over blocks and assemble the grid
    */
   
   GridOutput output(output_file_name,n_psiGhost,n_thetaGhost,block_mapping_tmp);
   
   double dpsi = (1.0 - psiNormInner)/n_psiCORE; //set the radial cell size
   double psiNorm_smallEnd(1.0);  // value of normalized psi at the separatrix
   double psiNorm_bigEnd;
   int n_psi, n_theta, n_extrp(0);
   
   Vector<FArrayBox*> grid;
   
   for (int block_id(0); block_id<NUM_SINGLE_NULL_BLOCKS; ++block_id) {
      
      // Skip MCORE and MCSOL, which, if being constructed, will be
      // constructed from the LCORE and RCORE blocks determined here
      if (block_id == MCORE || block_id == MCSOL) continue;

      if (block_id == RCORE || block_id == LCORE)  {
         n_psi = n_psiCORE;
         n_theta = n_thetaCORE/2;
         psiNorm_bigEnd = psiNormInner;
      }

      if (block_id == RCSOL || block_id == LCSOL || block_id == RSOL || block_id == LSOL  ) {
         n_psi = int(round((psiNormOuter-1.0)/dpsi));
         n_theta = n_thetaCORE/2;
         psiNorm_bigEnd = psiNormOuter;
      }

      if (block_id == RSOL) {
         psiNorm_bigEnd = psiNormOuter;
         n_theta = n_thetaOuterLeg;
         n_extrp = n_extrpOuterLeg;
      }

      if (block_id == LSOL) {
         psiNorm_bigEnd = psiNormOuter;
         n_theta = n_thetaInnerLeg;
         n_extrp = n_extrpInnerLeg;
      }

      if (block_id == RPF) {
         n_psi = n_psiPF;
         psiNorm_bigEnd = psiNorm_smallEnd - dpsi * n_psiPF;
         n_theta = n_thetaOuterLeg;
         n_extrp = n_extrpOuterLeg;
      }

      if (block_id == LPF) {
         n_psi = n_psiPF;
         psiNorm_bigEnd = psiNorm_smallEnd - dpsi * n_psiPF;
         n_theta = n_thetaInnerLeg;
         n_extrp = n_extrpInnerLeg;
      }

      //Get a grid on the unit square. The argument sets the size.
      FArrayBox xi;
      getUnitSquareGrid(n_psi, n_theta, xi);
      const Box& box(xi.box());

      //Construct the block mapping object
      BlockMapping block_mapping(DCT_file_name, block_id, psiNorm_bigEnd, n_extrp, trans_length, blending_factor);
      if (compute_crit_pts) {
         block_mapping.m_Xpoint  = Xpoint;
         block_mapping.m_magAxis = Opoint;
      }
      block_mapping.init(trans_rad);
      block_mapping.m_lengthCORE = lengthCORE;
      block_mapping.m_nthetaCORE = n_thetaCORE;
      
      //Get the physical coordinates
      FArrayBox physical_coordinates(box,2);
      if (block_id < LSOL) block_mapping.getPhysicalCoordinates(physical_coordinates, mappingCORE);
      else if (block_id == LPF || block_id == LSOL) block_mapping.getPhysicalCoordinates(physical_coordinates, mappingInnerLeg);
      else block_mapping.getPhysicalCoordinates(physical_coordinates, mappingOuterLeg);

      //Assemble the grid
      FArrayBox ordered_coords(box,2);
      output.orderBlock(block_id, physical_coordinates, ordered_coords);
      
      FArrayBox *ordered_coords_ptr = new FArrayBox;
      (*ordered_coords_ptr).define(box,2);
      (*ordered_coords_ptr).copy(ordered_coords);

      grid.push_back(ordered_coords_ptr);

      if (procID()==0) {
         cout<< "Done computing block " << block_id << " coordinates" << endl;
      }

      // Perform search of the minimal cell size
      bool grid_diagnostics = true;
      if (grid_diagnostics) {
         for (int dir=0; dir<2; ++dir) {
            Real size;
            RealVect loc;
            block_mapping.getMinimalCellSize(size, loc, physical_coordinates, dir);
            if (procID()==0) {
               cout<< "Minimal cell size in dir " << dir << " is " << size << " at location " << loc << endl;
            }
         }
      }

   }
   
   /*
    Step 3: Output the grid
   */
   
   if (procID()==0) {
      output.writeGrid(grid, n_thetaMCORE, output_gridUe);
   }
   
  
#ifdef CH_MPI
   CH_TIMER_REPORT();
   MPI_Finalize();
#endif

   return 0;
}


