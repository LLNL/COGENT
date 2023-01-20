#include "VelCoordSys.H"
#include "newMappedGridIO.H"

#include "NamespaceHeader.H"

const std::string VelCoordSys::pp_name = "velocity_coord_sys";


VelCoordSys::VelCoordSys( ParmParse&               a_parm_parse,
                          const DisjointBoxLayout& a_grids,
                          const ProblemDomain&     a_domain,
                          const RealVect&          a_dX )
   : CartesianCS(a_dX),
     m_grids(a_grids),
     m_domain(a_domain)
{
   // Get the verbosity flag
   m_verbose = false;
   if (a_parm_parse.contains("verbose")) {
      a_parm_parse.get("verbose", m_verbose);
   }

   m_isDefined = true;
}



VelCoordSys::~VelCoordSys()
{
}



void
VelCoordSys::getMetricTerms( LevelData<FluxBox>& a_N,
                             LevelData<FluxBox>& a_tanGradN ) const
{
   CH_assert(a_N.getBoxes()==a_tanGradN.getBoxes());

   DataIterator dit1 = a_N.dataIterator();
   for (dit1.begin(); dit1.ok(); ++dit1) {
      getN(a_N[dit1], a_N[dit1].box());
   }
   a_N.exchange();

   DataIterator dit2 = a_tanGradN.dataIterator();
   for (dit2.begin(); dit2.ok(); ++dit2) {
      NewFourthOrderCoordSys::computeTangentialGrad(a_tanGradN[dit2], a_N[dit2], a_tanGradN[dit2].box());
   }
   a_tanGradN.exchange();
}



void
VelCoordSys::getCellVolumes( FArrayBox& a_volume ) const
{
  Box box(a_volume.box());
  Box grown_box(grow(box,1));
  FluxBox N(grown_box, getNumN());
  getN(N, grown_box);
  cellVol(a_volume, N, box);
}



void
VelCoordSys::getCellVolumes( LevelData<FArrayBox>& a_volume ) const
{
   DataIterator dit = a_volume.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
     getCellVolumes(a_volume[dit]);
   }
   a_volume.exchange();
}



void
VelCoordSys::getJ( LevelData<FArrayBox>& a_J ) const
{
  DataIterator dit = a_J.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
     getAvgJ(a_J[dit], a_J[dit].box());
  }
  a_J.exchange();
}



void
VelCoordSys::getFaceAreas( FluxBox& a_areas ) const
{
  const Box& box = a_areas.box();

  FluxBox N(box, SpaceDim*SpaceDim);

  getN(N, box);

  magnitudeN(a_areas, N, box);
}



void
VelCoordSys::getFaceAreas( LevelData<FluxBox>& a_areas ) const
{
  DataIterator dit = a_areas.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    getFaceAreas(a_areas[dit]);
  }
}



void VelCoordSys::plotCellData( const string&               a_file_name,  
                                const LevelData<FArrayBox>& a_data,
                                const double&               a_time ) const
{
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   Box domain_box = grids.physDomain().domainBox();
   domain_box.grow(a_data.ghostVect());
   WriteMappedUGHDF5(a_file_name.c_str(), grids, a_data, *this, domain_box, a_time);
}



void VelCoordSys::plotFaceData( const string&             a_file_name,
                                const LevelData<FluxBox>& a_data,
                                const double&             a_time ) const
{
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   LevelData<FArrayBox> data_cell(grids, a_data.nComp(), a_data.ghostVect());

   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      const FluxBox& this_data = a_data[dit];
      FArrayBox& this_data_cell = data_cell[dit];
      this_data_cell.setVal(0.);
      for (int dir=0; dir<SpaceDim; dir++) {
         const FArrayBox& this_data_dir = this_data[dir];
         FArrayBox shift_fab(this_data_dir.box(), this_data_dir.nComp());
         shift_fab.copy(this_data_dir);

         shift_fab.shiftHalf(dir,-1);
         this_data_cell.plus(shift_fab);
         shift_fab.shiftHalf(dir,2);
         this_data_cell.plus(shift_fab);
      }
      this_data_cell.mult(0.5/SpaceDim);
   }

   Box domain_box = grids.physDomain().domainBox();
   domain_box.grow(a_data.ghostVect());
   WriteMappedUGHDF5(a_file_name.c_str(), grids, data_cell, *this, domain_box, a_time);
}



#include "NamespaceFooter.H"
