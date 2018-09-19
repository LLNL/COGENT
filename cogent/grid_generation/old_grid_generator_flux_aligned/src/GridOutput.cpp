#include <fstream>
#include <float.h>
#include "GridOutput.H"

#include "NamespaceHeader.H"

enum SingleNullBlockType {LCORE,RCORE,LCSOL,RCSOL,LSOL,RSOL,LPF,RPF,MCORE,MCSOL,NUM_SINGLE_NULL_BLOCKS};
enum directions {RADIAL_DIR, POLOIDAL_DIR};

GridOutput::GridOutput(const string&         a_mapping_file_name,
                       const int             a_radial_extention,
                       const int             a_poloidal_extention,
                       const BlockMapping&   a_block_mapping)

  :
   m_mapping_file_name(a_mapping_file_name),
   m_ext_rad(a_radial_extention),
   m_ext_pol(a_poloidal_extention),
   m_block_mapping(a_block_mapping)

{
}


GridOutput::~GridOutput()
{
}

void
GridOutput::writeGrid(const Vector<FArrayBox*> & grid,
                      const int                  n_thetaMCORE)
{
   
   /*
    Write extended grid into a file.
    */
   
   //const Box& box1(grid[1]->box());
   const Box& box1((*grid[1]).box());
   
   char file_name[80];
   sprintf(file_name, m_mapping_file_name.c_str());
   FILE* fd;
   fd = fopen(file_name, "w");
   
   char block_file_name[80];
   sprintf(block_file_name, "coords");
   
   char extended_block_file_name[80];
   sprintf(extended_block_file_name, "extended_coords");
   
   int num_blocks = 8;
   int half_mcore = 0;
   if (n_thetaMCORE > 0) {
      CH_assert(n_thetaMCORE%2 == 0);
      half_mcore = n_thetaMCORE / 2;
      num_blocks = 10;
   }

   for (int block_id(0); block_id<num_blocks; ++block_id) {

      if (procID()==0) {
         cout<< "Writing block " << block_id << " coordinates" << endl;
      }

      FArrayBox extended_data;
      Box grown_box;
      
      if (block_id == LCORE) {
         
         const Box& box = grid[block_id]->box();
         grown_box = box;
      
         fprintf(fd, "lcore %d %d %d %d \n", box.size(0)-m_ext_rad, m_ext_rad, box.size(1)-half_mcore, m_ext_pol);
         
         grown_box.growLo(POLOIDAL_DIR,-half_mcore);
         grown_box.grow(POLOIDAL_DIR, m_ext_pol);
         grown_box.growHi(RADIAL_DIR, m_ext_rad);
         extended_data.define(grown_box,2);
         
         for (BoxIterator bit(grown_box);bit.ok();++bit) {
            IntVect iv = bit();
            
            if (iv[POLOIDAL_DIR] < 0 && iv[RADIAL_DIR] < box.size(0)) {
               IntVect ip(iv);
               ip[1] = ((*grid[RCORE]).box()).size(1) - 1 + iv[1];
               extended_data(iv,0) = (*grid[RCORE])(ip,0);
               extended_data(iv,1) = (*grid[RCORE])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] < 0 && iv[RADIAL_DIR] >= box.size(0)) {
               IntVect ip(iv);
               ip[0] = iv[0] - box.size(0) + 1;
               ip[1] = ((*grid[RCSOL]).box()).size(1) - 1 + iv[1];
               extended_data(iv,0) = (*grid[RCSOL])(ip,0);
               extended_data(iv,1) = (*grid[RCSOL])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= 0 && iv[POLOIDAL_DIR] < box.size(1) && iv[RADIAL_DIR] < box.size(0))  {
               IntVect ip(iv);
               extended_data(iv,0) = (*grid[LCORE])(ip,0);
               extended_data(iv,1) = (*grid[LCORE])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= 0 && iv[POLOIDAL_DIR] < box.size(1) && iv[RADIAL_DIR] >= box.size(0))  {
               IntVect ip(iv);
               ip[0] = iv[0] - box.size(0) + 1;
               extended_data(iv,0) = (*grid[LCSOL])(ip,0);
               extended_data(iv,1) = (*grid[LCSOL])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= box.size(1) && iv[RADIAL_DIR] < box.size(0))  {
               IntVect ip(iv);
               ip[1] = iv[1] - box.size(1) + 1;
               extended_data(iv,0) = (*grid[RCORE])(ip,0);
               extended_data(iv,1) = (*grid[RCORE])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= box.size(1) && iv[RADIAL_DIR] >= box.size(0))  {
               IntVect ip(iv);
               ip[0] = iv[0] - box.size(0) + 1;
               ip[1] = iv[1] - box.size(1) + 1;
               extended_data(iv,0) = (*grid[LSOL])(ip,0);
               extended_data(iv,1) = (*grid[LSOL])(ip,1);
            }

            RealVect phys_coord(extended_data(iv,0), extended_data(iv,1));
            RealVect RB = getRB(phys_coord);
            
            fprintf(fd, "%20.12e %20.12e %20.12e %20.12e \n", extended_data(iv,0), extended_data(iv,1), RB[0], RB[1]);
            
         }
      }
      
      else if (block_id == RCORE) {
         
         const Box& box = grid[block_id]->box();
         grown_box = box;

         fprintf(fd, "rcore %d %d %d %d \n", box.size(0)-m_ext_rad, m_ext_rad, box.size(1)-half_mcore, m_ext_pol);
         
         grown_box.growHi(POLOIDAL_DIR,-half_mcore);
         grown_box.grow(POLOIDAL_DIR, m_ext_pol);
         grown_box.growHi(RADIAL_DIR, m_ext_rad);
         extended_data.define(grown_box,2);
         
         for (BoxIterator bit(grown_box);bit.ok();++bit) {
            IntVect iv = bit();
            
            if (iv[POLOIDAL_DIR] < 0 && iv[RADIAL_DIR] < box.size(0)) {
               IntVect ip(iv);
               ip[1] = ((*grid[LCORE]).box()).size(1) - 1 + iv[1];
               extended_data(iv,0) = (*grid[LCORE])(ip,0);
               extended_data(iv,1) = (*grid[LCORE])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] < 0 && iv[RADIAL_DIR] >= box.size(0)) {
               IntVect ip(iv);
               ip[0] = iv[0] - box.size(0) + 1;
               ip[1] = ((*grid[RSOL]).box()).size(1) - 1 + iv[1];
               extended_data(iv,0) = (*grid[RSOL])(ip,0);
               extended_data(iv,1) = (*grid[RSOL])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= 0 && iv[POLOIDAL_DIR] < box.size(1) && iv[RADIAL_DIR] < box.size(0))  {
               IntVect ip(iv);
               extended_data(iv,0) = (*grid[RCORE])(ip,0);
               extended_data(iv,1) = (*grid[RCORE])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= 0 && iv[POLOIDAL_DIR] < box.size(1) && iv[RADIAL_DIR] >= box.size(0))  {
               IntVect ip(iv);
               ip[0] = iv[0] - box.size(0) + 1;
               extended_data(iv,0) = (*grid[RCSOL])(ip,0);
               extended_data(iv,1) = (*grid[RCSOL])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= box.size(1) && iv[RADIAL_DIR] < box.size(0))  {
               IntVect ip(iv);
               ip[1] = iv[1] - box.size(1) + 1;
               extended_data(iv,0) = (*grid[LCORE])(ip,0);
               extended_data(iv,1) = (*grid[LCORE])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= box.size(1) && iv[RADIAL_DIR] >= box.size(0))  {
               IntVect ip(iv);
               ip[0] = iv[0] - box.size(0) + 1;
               ip[1] = iv[1] - box.size(1) + 1;
               extended_data(iv,0) = (*grid[LCSOL])(ip,0);
               extended_data(iv,1) = (*grid[LCSOL])(ip,1);
            }
            
            RealVect phys_coord(extended_data(iv,0), extended_data(iv,1));
            RealVect RB = getRB(phys_coord);
            
            fprintf(fd, "%20.12e %20.12e %20.12e %20.12e \n", extended_data(iv,0), extended_data(iv,1), RB[0], RB[1]);
            
         }
      }
      
      else if (block_id == LCSOL) {
         
         const Box& box = grid[block_id]->box();
         grown_box = box;

         fprintf(fd, "lcsol %d %d %d %d \n", box.size(0)-m_ext_rad, m_ext_rad, box.size(1)-half_mcore, m_ext_pol);
         
         grown_box.growLo(POLOIDAL_DIR,-half_mcore);
         grown_box.grow(POLOIDAL_DIR, m_ext_pol);
         grown_box.growLo(RADIAL_DIR, m_ext_rad);
         extended_data.define(grown_box,2);
         
         for (BoxIterator bit(grown_box);bit.ok();++bit) {
            IntVect iv = bit();
            
            if (iv[POLOIDAL_DIR] < 0 && iv[RADIAL_DIR] >= 0) {
               IntVect ip(iv);
               ip[1] = ((*grid[RCSOL]).box()).size(1) - 1 + iv[1];
               extended_data(iv,0) = (*grid[RCSOL])(ip,0);
               extended_data(iv,1) = (*grid[RCSOL])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] < 0 && iv[RADIAL_DIR] < 0) {
               IntVect ip(iv);
               ip[0] = ((*grid[RCORE]).box()).size(0) - 1 + iv[0];
               ip[1] = ((*grid[RCORE]).box()).size(1) - 1 + iv[1];
               extended_data(iv,0) = (*grid[RCORE])(ip,0);
               extended_data(iv,1) = (*grid[RCORE])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= 0 && iv[POLOIDAL_DIR] < box.size(1) && iv[RADIAL_DIR] >= 0)  {
               IntVect ip(iv);
               extended_data(iv,0) = (*grid[LCSOL])(ip,0);
               extended_data(iv,1) = (*grid[LCSOL])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= 0 && iv[POLOIDAL_DIR] < box.size(1) && iv[RADIAL_DIR] < 0)  {
               IntVect ip(iv);
               ip[0] = ((*grid[LCORE]).box()).size(0) - 1 + iv[0];
               extended_data(iv,0) = (*grid[LCORE])(ip,0);
               extended_data(iv,1) = (*grid[LCORE])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= box.size(1) && iv[RADIAL_DIR] >= 0)  {
               IntVect ip(iv);
               ip[1] = iv[1] - box.size(1) + 1;
               extended_data(iv,0) = (*grid[LSOL])(ip,0);
               extended_data(iv,1) = (*grid[LSOL])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= box.size(1) && iv[RADIAL_DIR] < 0)  {
               IntVect ip(iv);
               ip[0] = ((*grid[LPF]).box()).size(0) - 1 + iv[0];
               ip[1] = iv[1] - box.size(1) + 1;
               extended_data(iv,0) = (*grid[LPF])(ip,0);
               extended_data(iv,1) = (*grid[LPF])(ip,1);
            }

            RealVect phys_coord(extended_data(iv,0), extended_data(iv,1));
            RealVect RB = getRB(phys_coord);

            fprintf(fd, "%20.12e %20.12e %20.12e %20.12e \n", extended_data(iv,0), extended_data(iv,1), RB[0], RB[1]);
            
         }
      }
      
      else if (block_id == RCSOL) {
         
         const Box& box = grid[block_id]->box();
         grown_box = box;

         fprintf(fd, "rcsol %d %d %d %d \n", box.size(0)-m_ext_rad, m_ext_rad, box.size(1)-half_mcore, m_ext_pol);
         
         grown_box.growHi(POLOIDAL_DIR,-half_mcore);
         grown_box.grow(POLOIDAL_DIR, m_ext_pol);
         grown_box.growLo(RADIAL_DIR, m_ext_rad);
         extended_data.define(grown_box,2);
         
         for (BoxIterator bit(grown_box);bit.ok();++bit) {
            IntVect iv = bit();
            
            if (iv[POLOIDAL_DIR] < 0 && iv[RADIAL_DIR] >= 0) {
               IntVect ip(iv);
               ip[1] = ((*grid[RSOL]).box()).size(1) - 1 + iv[1];
               extended_data(iv,0) = (*grid[RSOL])(ip,0);
               extended_data(iv,1) = (*grid[RSOL])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] < 0 && iv[RADIAL_DIR] < 0) {
               IntVect ip(iv);
               ip[0] = ((*grid[RPF]).box()).size(0) - 1 + iv[0];
               ip[1] = ((*grid[RPF]).box()).size(1) - 1 + iv[1];
               extended_data(iv,0) = (*grid[RPF])(ip,0);
               extended_data(iv,1) = (*grid[RPF])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= 0 && iv[POLOIDAL_DIR] < box.size(1) && iv[RADIAL_DIR] >= 0)  {
               IntVect ip(iv);
               extended_data(iv,0) = (*grid[RCSOL])(ip,0);
               extended_data(iv,1) = (*grid[RCSOL])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= 0 && iv[POLOIDAL_DIR] < box.size(1) && iv[RADIAL_DIR] < 0)  {
               IntVect ip(iv);
               ip[0] = ((*grid[RCORE]).box()).size(0) - 1 + iv[0];
               extended_data(iv,0) = (*grid[RCORE])(ip,0);
               extended_data(iv,1) = (*grid[RCORE])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= box.size(1) && iv[RADIAL_DIR] >= 0)  {
               IntVect ip(iv);
               ip[1] = iv[1] - box.size(1) + 1;
               extended_data(iv,0) = (*grid[LCSOL])(ip,0);
               extended_data(iv,1) = (*grid[LCSOL])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= box.size(1) && iv[RADIAL_DIR] < 0)  {
               IntVect ip(iv);
               ip[0] = ((*grid[LCORE]).box()).size(0) - 1 + iv[0];
               ip[1] = iv[1] - box.size(1) + 1;
               extended_data(iv,0) = (*grid[LCORE])(ip,0);
               extended_data(iv,1) = (*grid[LCORE])(ip,1);
            }

            RealVect phys_coord(extended_data(iv,0), extended_data(iv,1));
            RealVect RB = getRB(phys_coord);

            fprintf(fd, "%20.12e %20.12e %20.12e %20.12e \n", extended_data(iv,0), extended_data(iv,1), RB[0], RB[1]);
            
         }
      }
      
      else if (block_id == LSOL) {
         
         const Box& box = grid[block_id]->box();
         grown_box = box;

         fprintf(fd, "lsol %d %d %d %d \n", box.size(0)-m_ext_rad, m_ext_rad, box.size(1)-m_ext_pol, m_ext_pol);
         
         grown_box.growLo(POLOIDAL_DIR, m_ext_pol);
         grown_box.growLo(RADIAL_DIR, m_ext_rad);
         extended_data.define(grown_box,2);
         
         for (BoxIterator bit(grown_box);bit.ok();++bit) {
            IntVect iv = bit();
            
            if (iv[POLOIDAL_DIR] < 0 && iv[RADIAL_DIR] >= 0) {
               IntVect ip(iv);
               ip[1] = ((*grid[LCSOL]).box()).size(1) - 1 + iv[1];
               extended_data(iv,0) = (*grid[LCSOL])(ip,0);
               extended_data(iv,1) = (*grid[LCSOL])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] < 0 && iv[RADIAL_DIR] < 0) {
               IntVect ip(iv);
               ip[0] = ((*grid[LCORE]).box()).size(0) - 1 + iv[0];
               ip[1] = ((*grid[LCORE]).box()).size(1) - 1 + iv[1];
               extended_data(iv,0) = (*grid[LCORE])(ip,0);
               extended_data(iv,1) = (*grid[LCORE])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= 0 && iv[RADIAL_DIR] >= 0)  {
               IntVect ip(iv);
               extended_data(iv,0) = (*grid[LSOL])(ip,0);
               extended_data(iv,1) = (*grid[LSOL])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= 0 && iv[RADIAL_DIR] < 0)  {
               IntVect ip(iv);
               ip[0] = ((*grid[LPF]).box()).size(0) - 1 + iv[0];
               extended_data(iv,0) = (*grid[LPF])(ip,0);
               extended_data(iv,1) = (*grid[LPF])(ip,1);
            }
            
            RealVect phys_coord(extended_data(iv,0), extended_data(iv,1));
            RealVect RB = getRB(phys_coord);

            fprintf(fd, "%20.12e %20.12e %20.12e %20.12e \n", extended_data(iv,0), extended_data(iv,1), RB[0], RB[1]);
            
         }
      }
      
      else if (block_id == RSOL) {
         
         const Box& box = grid[block_id]->box();
         grown_box = box;

         fprintf(fd, "rsol %d %d %d %d \n", box.size(0)-m_ext_rad, m_ext_rad, box.size(1)-m_ext_pol, m_ext_pol);
         
         grown_box.growHi(POLOIDAL_DIR, m_ext_pol);
         grown_box.growLo(RADIAL_DIR, m_ext_rad);
         extended_data.define(grown_box,2);
         
         for (BoxIterator bit(grown_box);bit.ok();++bit) {
            IntVect iv = bit();
            
            if (iv[POLOIDAL_DIR] < box.size(1) && iv[RADIAL_DIR] >= 0) {
               IntVect ip(iv);
               extended_data(iv,0) = (*grid[RSOL])(ip,0);
               extended_data(iv,1) = (*grid[RSOL])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] < box.size(1) && iv[RADIAL_DIR] < 0) {
               IntVect ip(iv);
               ip[0] = ((*grid[RPF]).box()).size(0) - 1 + iv[0];
               extended_data(iv,0) = (*grid[RPF])(ip,0);
               extended_data(iv,1) = (*grid[RPF])(ip,1);            }
            
            if (iv[POLOIDAL_DIR] >= box.size(1) && iv[RADIAL_DIR] >= 0)  {
               IntVect ip(iv);
               ip[1] = iv[1] - box.size(1) + 1;
               extended_data(iv,0) = (*grid[RCSOL])(ip,0);
               extended_data(iv,1) = (*grid[RCSOL])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= box.size(1) && iv[RADIAL_DIR] < 0)  {
               IntVect ip(iv);
               ip[0] = ((*grid[RCORE]).box()).size(0) - 1 + iv[0];
               ip[1] = iv[1] - box.size(1) + 1;
               extended_data(iv,0) = (*grid[RCORE])(ip,0);
               extended_data(iv,1) = (*grid[RCORE])(ip,1);
            }
            
            RealVect phys_coord(extended_data(iv,0), extended_data(iv,1));
            RealVect RB = getRB(phys_coord);

            fprintf(fd, "%20.12e %20.12e %20.12e %20.12e \n", extended_data(iv,0), extended_data(iv,1), RB[0], RB[1]);
            
         }
      }
      
      else if (block_id == LPF) {
         
         const Box& box = grid[block_id]->box();
         grown_box = box;

         fprintf(fd, "lpf %d %d %d %d \n", box.size(0)-m_ext_rad, m_ext_rad, box.size(1)-m_ext_pol, m_ext_pol);
         
         grown_box.growLo(POLOIDAL_DIR, m_ext_pol);
         grown_box.growHi(RADIAL_DIR, m_ext_rad);
         extended_data.define(grown_box,2);
         
         for (BoxIterator bit(grown_box);bit.ok();++bit) {
            IntVect iv = bit();
            
            if (iv[POLOIDAL_DIR] < 0 && iv[RADIAL_DIR] < box.size(0)) {
               IntVect ip(iv);
               ip[1] = ((*grid[RPF]).box()).size(1) - 1 + iv[1];
               extended_data(iv,0) = (*grid[RPF])(ip,0);
               extended_data(iv,1) = (*grid[RPF])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] < 0 && iv[RADIAL_DIR] >= box.size(0) ) {
               IntVect ip(iv);
               ip[0] = iv[0] - box.size(0) + 1;
               ip[1] = ((*grid[LCSOL]).box()).size(1) - 1 + iv[1];
               extended_data(iv,0) = (*grid[LCSOL])(ip,0);
               extended_data(iv,1) = (*grid[LCSOL])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= 0 && iv[RADIAL_DIR] < box.size(0))  {
               IntVect ip(iv);
               extended_data(iv,0) = (*grid[LPF])(ip,0);
               extended_data(iv,1) = (*grid[LPF])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= 0 && iv[RADIAL_DIR] >= box.size(0))  {
               IntVect ip(iv);
               ip[0] = iv[0] - box.size(0) + 1;
               extended_data(iv,0) = (*grid[LSOL])(ip,0);
               extended_data(iv,1) = (*grid[LSOL])(ip,1);
            }
            
            RealVect phys_coord(extended_data(iv,0), extended_data(iv,1));
            RealVect RB = getRB(phys_coord);

            fprintf(fd, "%20.12e %20.12e %20.12e %20.12e \n", extended_data(iv,0), extended_data(iv,1), RB[0], RB[1]);
            
         }
      }
      
      else if (block_id == RPF) {
         
         const Box& box = grid[block_id]->box();
         grown_box = box;

         fprintf(fd, "rpf %d %d %d %d \n", box.size(0)-m_ext_rad, m_ext_rad, box.size(1)-m_ext_pol, m_ext_pol);
         
         grown_box.growHi(POLOIDAL_DIR, m_ext_pol);
         grown_box.growHi(RADIAL_DIR, m_ext_rad);
         extended_data.define(grown_box,2);
         
         for (BoxIterator bit(grown_box);bit.ok();++bit) {
            IntVect iv = bit();
            
            if (iv[POLOIDAL_DIR] < box.size(1) && iv[RADIAL_DIR] < box.size(0)) {
               IntVect ip(iv);
               extended_data(iv,0) = (*grid[RPF])(ip,0);
               extended_data(iv,1) = (*grid[RPF])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] < box.size(1) && iv[RADIAL_DIR] >= box.size(0) ) {
               IntVect ip(iv);
               ip[0] = iv[0] - box.size(0) + 1;
               extended_data(iv,0) = (*grid[RSOL])(ip,0);
               extended_data(iv,1) = (*grid[RSOL])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= box.size(1) && iv[RADIAL_DIR] < box.size(0))  {
               IntVect ip(iv);
               ip[1] = iv[1] - box.size(1) + 1;
               extended_data(iv,0) = (*grid[LPF])(ip,0);
               extended_data(iv,1) = (*grid[LPF])(ip,1);
            }
            
            if (iv[POLOIDAL_DIR] >= box.size(1) && iv[RADIAL_DIR] >= box.size(0))  {
               IntVect ip(iv);
               ip[0] = iv[0] - box.size(0) + 1;
               ip[1] = iv[1] - box.size(1) + 1;
               extended_data(iv,0) = (*grid[RCSOL])(ip,0);
               extended_data(iv,1) = (*grid[RCSOL])(ip,1);
            }
            
            RealVect phys_coord(extended_data(iv,0), extended_data(iv,1));
            RealVect RB = getRB(phys_coord);

            fprintf(fd, "%20.12e %20.12e %20.12e %20.12e \n", extended_data(iv,0), extended_data(iv,1), RB[0], RB[1]);
            
         }
      }
      
      else if (block_id == MCORE) {
         
         const Box& lcore_box = grid[LCORE]->box();

         fprintf(fd, "mcore %d %d %d %d \n", lcore_box.size(0)-m_ext_rad, m_ext_rad, 2*half_mcore+1, m_ext_pol);
         
         const Box box(IntVect(lcore_box.smallEnd(RADIAL_DIR),0),IntVect(lcore_box.bigEnd(RADIAL_DIR),2*half_mcore));
         grown_box = box;

         grown_box.grow(POLOIDAL_DIR, m_ext_pol);
         grown_box.growHi(RADIAL_DIR, m_ext_rad);
         extended_data.define(grown_box,2);

         for (BoxIterator bit(grown_box);bit.ok();++bit) {
            IntVect iv = bit();
            
            if (iv[POLOIDAL_DIR] < half_mcore + 1) {
               if (iv[RADIAL_DIR] < box.size(0)) {
                  IntVect ip(iv);
                  ip[1] = ((*grid[RCORE]).box()).size(1) - 1 + iv[1] - half_mcore;
                  extended_data(iv,0) = (*grid[RCORE])(ip,0);
                  extended_data(iv,1) = (*grid[RCORE])(ip,1);
               }
               else {
                  IntVect ip;
                  ip[0] = iv[0] - box.size(0) + 1;
                  ip[1] = ((*grid[RCSOL]).box()).size(1) - 1 + iv[1] - half_mcore;
                  extended_data(iv,0) = (*grid[RCSOL])(ip,0);
                  extended_data(iv,1) = (*grid[RCSOL])(ip,1);
               }
            }
            else {
               if (iv[RADIAL_DIR] < box.size(0))  {
                  IntVect ip(iv);
                  ip[1] = iv[1] - half_mcore;
                  extended_data(iv,0) = (*grid[LCORE])(ip,0);
                  extended_data(iv,1) = (*grid[LCORE])(ip,1);
               }
               else {
                  IntVect ip;
                  ip[0] = iv[0] - box.size(0) + 1;
                  ip[1] = iv[1] - half_mcore;
                  extended_data(iv,0) = (*grid[LCSOL])(ip,0);
                  extended_data(iv,1) = (*grid[LCSOL])(ip,1);
               }
            }
            
            RealVect phys_coord(extended_data(iv,0), extended_data(iv,1));
            RealVect RB = getRB(phys_coord);

            fprintf(fd, "%20.12e %20.12e %20.12e %20.12e \n", extended_data(iv,0), extended_data(iv,1), RB[0], RB[1]);
            
         }
      }
         
      else if (block_id == MCSOL) {
         
         const Box& lcsol_box = grid[LCSOL]->box();

         fprintf(fd, "mcsol %d %d %d %d \n", lcsol_box.size(0)-m_ext_rad, m_ext_rad, 2*half_mcore+1, m_ext_pol);
         
         const Box box(IntVect(lcsol_box.smallEnd(RADIAL_DIR),0),IntVect(lcsol_box.bigEnd(RADIAL_DIR),2*half_mcore));
         grown_box = box;

         grown_box.grow(POLOIDAL_DIR, m_ext_pol);
         grown_box.growLo(RADIAL_DIR, m_ext_rad);
         extended_data.define(grown_box,2);
         
         for (BoxIterator bit(grown_box);bit.ok();++bit) {
            IntVect iv = bit();
            
            if (iv[POLOIDAL_DIR] < half_mcore + 1) {
               if (iv[RADIAL_DIR] < 0) {
                  IntVect ip;
                  ip[0] = ((*grid[RCORE]).box()).size(0) - 1 + iv[0];
                  ip[1] = ((*grid[RCSOL]).box()).size(1) - 1 + iv[1] - half_mcore;
                  extended_data(iv,0) = (*grid[RCORE])(ip,0);
                  extended_data(iv,1) = (*grid[RCORE])(ip,1);
               }
               else {
                  IntVect ip(iv);
                  ip[1] = ((*grid[RCSOL]).box()).size(1) - 1 + iv[1] - half_mcore;
                  extended_data(iv,0) = (*grid[RCSOL])(ip,0);
                  extended_data(iv,1) = (*grid[RCSOL])(ip,1);
               }
            }
            else {
               if (iv[RADIAL_DIR] < 0)  {
                  IntVect ip(iv);
                  ip[0] = ((*grid[LCORE]).box()).size(0) - 1 + iv[0];
                  ip[1] = iv[1] - half_mcore;
                  extended_data(iv,0) = (*grid[LCORE])(ip,0);
                  extended_data(iv,1) = (*grid[LCORE])(ip,1);
               }
               else {
                  IntVect ip(iv);
                  ip[1] = iv[1] - half_mcore;
                  extended_data(iv,0) = (*grid[LCSOL])(ip,0);
                  extended_data(iv,1) = (*grid[LCSOL])(ip,1);
               }
            }
            
            RealVect phys_coord(extended_data(iv,0), extended_data(iv,1));
            RealVect RB = getRB(phys_coord);

            fprintf(fd, "%20.12e %20.12e %20.12e %20.12e \n", extended_data(iv,0), extended_data(iv,1), RB[0], RB[1]);
         }
      }

      writeBlockCoordinates(extended_block_file_name, extended_data, block_id);

      Box valid_box = grow(grown_box,-IntVect(m_ext_rad,m_ext_pol));
      FArrayBox valid_data(valid_box,2);
      valid_data.copy(extended_data);

      writeBlockCoordinates(block_file_name, valid_data, block_id);
   }
}

void
GridOutput::writeBlockCoordinates(const string& a_file_name,
                                  const FArrayBox& a_physical_coordinates,
                                  const int a_block_number)
{
   char file_name[80];
   sprintf(file_name, "%s%s%d", "output/",a_file_name.c_str(), a_block_number);
   FILE* fd = fopen(file_name, "w");
   
   const Box& box = a_physical_coordinates.box();
   
   fprintf(fd, "%d %d\n", box.size(0), box.size(1));
   
   for (BoxIterator bit(box);bit.ok();++bit) {
      IntVect iv = bit();
      fprintf(fd, "%20.12e %20.12e \n", a_physical_coordinates(iv,0), a_physical_coordinates(iv,1));
   }
   
   fclose(fd);
}

void
GridOutput::orderBlock(const int a_block_id,
                       const FArrayBox & a_coords,
                       FArrayBox& a_ordered_coords)
{
   
   /*
    Order grid consistent with cogent mapping.
    */
   
   const Box& box(a_coords.box());
   for (BoxIterator bit(box);bit.ok();++bit) {
      IntVect iv = bit();
      
      IntVect ip(iv);
      if (a_block_id == RCSOL || a_block_id == LSOL) {
         //do nothing
      }
      
      if (a_block_id == RCORE || a_block_id == LPF) {
         ip[0] = box.size(0) - 1 - iv[0];
      }
      
      if (a_block_id == LCSOL || a_block_id == RSOL) {
         ip[1] = box.size(1) - 1 - iv[1];
      }
      
      if (a_block_id == LCORE || a_block_id == RPF) {
         ip[0] = box.size(0) - 1 - iv[0];
         ip[1] = box.size(1) - 1 - iv[1];
      }
      
      a_ordered_coords(iv,0) = a_coords(ip,0);
      a_ordered_coords(iv,1) = a_coords(ip,1);
   }
}

RealVect
GridOutput::getRB(const RealVect& a_phys_coord)
{

   Vector<Real> mag_data = m_block_mapping.getMagFieldData( a_phys_coord );
   double R = a_phys_coord[0];

   RealVect RB;
   RB[0] = R * mag_data[0]; //R*Br
   RB[1] = R * mag_data[1]; //R*Bz
  
   return RB;
}

#include "NamespaceFooter.H"


