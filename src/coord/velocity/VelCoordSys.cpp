#include "VelCoordSys.H"

#include "NamespaceHeader.H"

char* VelCoordSys::pp_name = {"velocity_coord_sys"};



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



#include "NamespaceFooter.H"
