#ifndef _GRIDOUTPUT_H_
#define _GRIDOUTPUT_H_


#include "ParmParse.H"
#include "RealVect.H"
#include "IntVect.H"
#include "FArrayBox.H"
#include "BoxIterator.H"
#include "BlockMapping.H"

#include "NamespaceHeader.H"

class GridOutput
{
  public:

   GridOutput(const string&         mapping_file_name,
              const int             radial_extention,
              const int             poloidal_extention,
              const BlockMapping&   block_mapping);

   ~GridOutput();

   void writeGrid(const Vector<FArrayBox*> & grid,
                  const int                  num_mcore_poloidal_cells);

   void orderBlock(const int block_id,
                   const FArrayBox & coords,
                   FArrayBox& ordered_coords);
   
   void writeBlockCoordinates(const string& a_file_name,
                              const FArrayBox& physical_coordinates,
                              const int block_number);
  private:

   
   RealVect getRB(const RealVect& physical_coordinate);
   
   string m_mapping_file_name;
   int m_ext_rad;
   int m_ext_pol;
   const BlockMapping& m_block_mapping;


};

#include "NamespaceFooter.H"

#endif
