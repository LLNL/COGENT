#include "OneDimData.H"

#include <iostream>
#include <typeinfo>
#include <string>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <istream>

using namespace std;

#include "rbf_interp_2d.hpp"
#include "r8lib.hpp"
#include "test_interp_2d.hpp"

#include "DataIterator.H"
#include "Directions.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "FourthOrderUtil.H"
#include "LevelData.H"
#include "MayDay.H"
#include "SNCoreBlockCoordSys.H"
#include "SingleNullBlockCoordSys.H"
#include "SNCoreBlockCoordSysModel.H"
#include "SingleNullBlockCoordSysModel.H"
#include "MagGeom.H"

#include "MultiBlockCoordSys.H"
#include "Vector.H"

#include "NamespaceHeader.H"

inline
const MagBlockCoordSys& getCoordSys( const MultiBlockLevelGeom& a_geometry,
                                     const Box& a_box )
{
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   const int block_number( coord_sys.whichBlock( a_box ) );
   const NewCoordSys* block_coord_sys( coord_sys.getCoordSys( block_number ) );
   return static_cast<const MagBlockCoordSys&>( *block_coord_sys );
}


OneDimData::OneDimData( ParmParse& a_pp,
                const int& a_verbosity )
   : m_verbosity(a_verbosity),
     m_subtype("MBA"),
     m_coord_type("flux"),
     m_coords(NULL),
     m_weights(NULL),
     m_MBA(NULL),
     m_init_lattice_MBA(10)

{
   parseParameters( a_pp );
   if (m_subtype == "RBF") createInterpolationRBF();
   if (m_subtype == "MBA") createInterpolationMBA();
}

OneDimData::~OneDimData()
{
    if (m_coords) delete [] m_coords;
    if (m_weights) delete [] m_weights;
    if (m_MBA) delete m_MBA;
}


void OneDimData::assign( LevelData<FArrayBox>& a_data,
                     const MultiBlockLevelGeom& a_geometry,
                     const Real& a_time,
                     const bool& a_cell_averages ) const
{

   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );

   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
 
      int block_number( coord_sys.whichBlock( grids[dit] ) );

      if (a_cell_averages) {
         setCellAverages( a_data[dit], coord_sys, getCoordSys( a_geometry, grids[dit] ), block_number );
      }
      else {
         setPointwise( a_data[dit], getCoordSys( a_geometry, grids[dit] ) );
      }
   }
   a_data.exchange();
}


void OneDimData::assign( LevelData<FluxBox>& a_data,
                     const MultiBlockLevelGeom& a_geometry,
                     const Real& a_time,
                     const bool& a_cell_averages ) const
{

   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );

   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
 
      int block_number( coord_sys.whichBlock( grids[dit] ) );

      if (a_cell_averages) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            setCellAverages( a_data[dit][dir], coord_sys, getCoordSys( a_geometry, grids[dit] ), block_number );
         }
      }
      else {
         for (int dir=0; dir<SpaceDim; ++dir) {
            setPointwise( a_data[dit][dir], getCoordSys( a_geometry, grids[dit] ) );
         }
      }
   }
   a_data.exchange();
}


void OneDimData::assign( FArrayBox& a_data,
                     const MultiBlockLevelGeom& a_geometry,
                     const Box& a_box, // interior box
                     const Real& a_time,
                     const bool& a_cell_averages ) const
{
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   const int block_number( coord_sys.whichBlock( a_box ) );

   if (a_cell_averages) {
      setCellAverages( a_data, coord_sys, getCoordSys( a_geometry, a_box ), block_number );
   }
   else {
      setPointwise( a_data, getCoordSys( a_geometry, a_box ) );
   }
}


void OneDimData::assign( FluxBox& a_data,
                     const MultiBlockLevelGeom& a_geometry,
                     const Box& a_box, // interior box
                     const Real& a_time,
                     const bool& a_cell_averages ) const
{
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   const int block_number( coord_sys.whichBlock( a_box ) );

   if (a_cell_averages) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         setCellAverages( a_data[dir], coord_sys, getCoordSys( a_geometry, a_box ), block_number );
      }
   }
   else {
      for (int dir=0; dir<SpaceDim; ++dir) {
         setPointwise( a_data[dir], getCoordSys( a_geometry, a_box ) );
      }
   }
}


void OneDimData::assign( LevelData<FArrayBox>& a_data,
                     const MultiBlockLevelGeom& a_geometry,
                     const BoundaryBoxLayout& a_bdry_layout,
                     const Real& a_time ) const
{
   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
   // NB: This is a cheat - there's one too many cells at the (dir,side) face
   // of the boundary box, but it doesn't matter because one-sided difference
   // will be used at that face to construct the cell average.  We do want the
   // extra cell in all other directions.
   LevelData<FArrayBox> data_tmp( grids, a_data.nComp(), IntVect::Unit );
   for (DataIterator dit( grids ); dit.ok(); ++dit) {
      const Box box( a_bdry_layout.interiorBox( dit ) );
      const MagBlockCoordSys& coord_sys( ((MagGeom&)a_geometry).getBlockCoordSys( box ) );

      setPointwise( data_tmp[dit], coord_sys );
   }
   for (DataIterator dit( grids ); dit.ok(); ++dit) {
      Box domain_box( data_tmp[dit].box() );
      domain_box.growDir( a_bdry_layout.dir(), a_bdry_layout.side(), -1 );
      ProblemDomain domain( domain_box );
      fourthOrderAverageCell( data_tmp[dit], domain, grids[dit] );
   }
   data_tmp.copyTo( a_data );
   a_data.exchange();
}


inline
void OneDimData::parseParameters( ParmParse& a_pp )
{
   
   if ( a_pp.contains("data_file")) {
        a_pp.get("data_file", m_data_file);
   }
   else {
       MayDay::Error("OneDimData:: No data file specified");
   }

   a_pp.query( "subtype", m_subtype );
   a_pp.query( "coordinate_type", m_coord_type );
    
   if (m_verbosity) {
      printParameters();
   }
}


inline
void OneDimData::setCellAverages( FArrayBox&        a_data,
                              const MultiBlockCoordSys& a_coords,
                              const MagBlockCoordSys&   a_block_coord,
                              const int               block_number ) const

{
   Box box( a_data.box() );
   Box tmp_box( box );
   tmp_box.grow( IntVect::Unit );
   FArrayBox tmp( tmp_box, a_data.nComp() );

   setPointwise( tmp, a_block_coord );

   fourthOrderAverageCell( tmp, a_block_coord.domain(), box );

   a_data.copy( tmp, box );
}


inline
void OneDimData::setPointwise(FArrayBox&                a_data,
                              const MagBlockCoordSys&   a_coord_sys )  const

{
   
   //Check geometry validity for "flux" or "outer_midplane" types
   if (m_coord_type == "flux" || m_coord_type == "outer_midplane") {
      bool not_sn( typeid(a_coord_sys) != typeid(SNCoreBlockCoordSys) );
      not_sn &= (typeid(a_coord_sys) != typeid(SNCoreBlockCoordSysModel));
      not_sn &= (typeid(a_coord_sys) != typeid(SingleNullBlockCoordSys));
      not_sn &= (typeid(a_coord_sys) != typeid(SingleNullBlockCoordSysModel));
      if ( not_sn ) {
         const std::string msg( "OneDimData: Attempt to use not a single-null geometry with the flux or outer_midplane options. ");
         MayDay::Error( msg.c_str() );
      }
   }

   Box box( a_data.box() );
   FArrayBox cc_phys_coords( box, SpaceDim );
   a_coord_sys.getCellCenteredRealCoords( cc_phys_coords );

   FArrayBox cc_mapped_coords( box, SpaceDim );
   a_coord_sys.getCellCenteredMappedCoords( cc_mapped_coords );
   

   //We need this for the RBF method
   int count = 0;
   int npoints = box.numPts();
   double *box_coords = new double[npoints];
   
   
   BoxIterator bit(a_data.box());
   for (bit.begin(); bit.ok(); ++bit)
   {
      
      double coord;
      
      IntVect iv = bit();
      RealVect loc(iv);
      RealVect phys_coordinate(iv);
      
      //Here assume that dir == RADIAL_DIR, however that can be changed to any direction
      int dir = RADIAL_DIR;

      if (m_coord_type == "physical")  coord = cc_phys_coords(iv, dir);
      if (m_coord_type == "mapped" || m_coord_type == "flux" || m_coord_type == "outer_midplane")  coord = cc_mapped_coords(iv, dir);
   
      for (int dir=0; dir<SpaceDim; dir++) {
         phys_coordinate[dir] = cc_phys_coords(iv, dir);
      }
      
      if (m_coord_type == "flux") {
         if (typeid(a_coord_sys) == typeid(SNCoreBlockCoordSys))  coord = ((const SNCoreBlockCoordSys&)a_coord_sys).getNormMagneticFlux(phys_coordinate);
         if (typeid(a_coord_sys) == typeid(SingleNullBlockCoordSys))  coord = ((const SingleNullBlockCoordSys&)a_coord_sys).getNormMagneticFlux(phys_coordinate);
         
         if (typeid(a_coord_sys) == typeid(SNCoreBlockCoordSysModel))  coord = ((const SNCoreBlockCoordSysModel&)a_coord_sys).getNormMagneticFlux(phys_coordinate);
         if (typeid(a_coord_sys) == typeid(SingleNullBlockCoordSysModel))  coord = ((const SingleNullBlockCoordSysModel&)a_coord_sys).getNormMagneticFlux(phys_coordinate);
      }
      
      if (m_coord_type == "outer_midplane") {
         
         if (typeid(a_coord_sys) == typeid(SNCoreBlockCoordSys)) {
            double fluxNorm  = ((const SNCoreBlockCoordSys&)a_coord_sys).getNormMagneticFlux(phys_coordinate);
            double Rsep  = ((const SNCoreBlockCoordSys&)a_coord_sys).getOuterRsep();
            coord = ((const SNCoreBlockCoordSys&)a_coord_sys).getOuterMidplaneCoord(fluxNorm) - Rsep;
            
         }
         
         if (typeid(a_coord_sys) == typeid(SingleNullBlockCoordSys)) {
            double fluxNorm  = ((const SingleNullBlockCoordSys&)a_coord_sys).getNormMagneticFlux(phys_coordinate);
            double Rsep  = ((const SingleNullBlockCoordSys&)a_coord_sys).getOuterRsep();
            coord = ((const SingleNullBlockCoordSys&)a_coord_sys).getOuterMidplaneCoord(fluxNorm) - Rsep;
         }
         
         if (typeid(a_coord_sys) == typeid(SNCoreBlockCoordSysModel)) {
            double fluxNorm  = ((const SNCoreBlockCoordSysModel&)a_coord_sys).getNormMagneticFlux(phys_coordinate);
            double Rsep  = ((const SNCoreBlockCoordSysModel&)a_coord_sys).getOuterRsep();
            coord = ((const SNCoreBlockCoordSysModel&)a_coord_sys).getOuterMidplaneCoord(fluxNorm) - Rsep;
            
         }
         
         if (typeid(a_coord_sys) == typeid(SingleNullBlockCoordSysModel)) {
            double fluxNorm  = ((const SingleNullBlockCoordSysModel&)a_coord_sys).getNormMagneticFlux(phys_coordinate);
            double Rsep  = ((const SingleNullBlockCoordSysModel&)a_coord_sys).getOuterRsep();
            coord = ((const SingleNullBlockCoordSysModel&)a_coord_sys).getOuterMidplaneCoord(fluxNorm) - Rsep;
         }
         
      }

      if (m_subtype == "MBA") {
         std::array<double,1> this_point = mba::make_array<double>(coord);
         a_data(iv,0) = m_MBA->getInterpolatedValue(this_point);
      }
      
      if (m_subtype == "RBF") {
         box_coords[count] = coord;
         count++;
      }
      
   }
   
   if (m_subtype == "RBF") {
    
   
      //Choise of the radial basis functions are given by phi1, phi2, phi3, and phi4 (see rbf_interp_2d.cpp)
      //and should be consistent with phi in rbf_weights
      double *interpolated_data;
      interpolated_data = rbf_interp ( 1, m_data_npoints, m_coords,
                                   m_r0, phi1, m_weights, npoints, box_coords );

    
      //Copy interpolated data to a_data
      int i=0;
      for (bit.begin(); bit.ok(); ++bit) {
         IntVect iv = bit();
         a_data(iv,0)=interpolated_data[i];
         i++;
      }

      delete [] interpolated_data;
   }

   delete [] box_coords;
   
}


void OneDimData::createInterpolationRBF()
{
    
#ifdef CH_MPI
    if (procID() == 0) {
#endif
     
     //Get number of data points lines (i.e., number of line in a data file)
     int npoints = 0;
     string line;
     ifstream datafile(m_data_file.c_str());
        
     while (std::getline(datafile, line))  ++npoints;

        m_data_npoints = npoints;
#ifdef CH_MPI
    }
    
    MPI_Bcast(&m_data_npoints, 1, MPI_INT, 0, MPI_COMM_WORLD);

#endif
    
    m_coords = new double[m_data_npoints];
    m_weights = new double[m_data_npoints];

#ifdef CH_MPI
    if (procID() == 0) {
#endif
    
     //Reading the table data
     double *data = new double[m_data_npoints];
        
     ifstream inFile;
     inFile.open( m_data_file.c_str() );
   
     for (int i=0; i<m_data_npoints; ++i) {
        inFile >> m_coords[i];
        inFile >> data[i];
     }
     inFile.close();
        
     double max = r8vec_max ( m_data_npoints, m_coords );
     double min = r8vec_min ( m_data_npoints, m_coords );
     double volume = ( max - min ) ;
        
     double e = 1.0 / ( double ) ( 1 );
     m_r0 = pow ( volume / m_data_npoints, e );
        
     //Choise of the radial basis functions are given by phi1, phi2, phi3, and phi4 (see rbf_interp_2d.cpp)
     //and should be consistent with phi in rbf_interp
     m_weights = rbf_weight ( 1, m_data_npoints, m_coords, m_r0, phi1, data );
        
#ifdef CH_MPI
    }

    
    MPI_Bcast(&m_r0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(m_coords, m_data_npoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(m_weights, m_data_npoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#endif
    
}

void OneDimData::createInterpolationMBA()
{
    
#ifdef CH_MPI
    if (procID() == 0) {
#endif
        
        //Get number of data points lines (i.e., number of line in a data file)
        int npoints = 0;
        string line;
        ifstream datafile(m_data_file.c_str());
        
        while (std::getline(datafile, line)) ++npoints;
        
        m_data_npoints = npoints;
#ifdef CH_MPI
    }
    
    MPI_Bcast(&m_data_npoints, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
#endif
    
    double *coords = new double[m_data_npoints];
    double *data = new double[m_data_npoints];
    
#ifdef CH_MPI
    if (procID() == 0) {
#endif
        
        //Reading the table data
        ifstream inFile;
        inFile.open( m_data_file.c_str() );
        
        for (int i=0; i<m_data_npoints; ++i) {
            inFile >> coords[i];
            inFile >> data[i];
        }
        inFile.close();
        
#ifdef CH_MPI
    }
    
    
    MPI_Bcast(coords, m_data_npoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(data, m_data_npoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
#endif
    
    double min = r8vec_min ( m_data_npoints, coords );
    double max = r8vec_max ( m_data_npoints, coords );
   
    //Creating the MBA interpolant
    std::vector< std::array<double,1> > data_points ;
    std::vector<double> data_values;
    
    for(int i = 0; i < m_data_npoints; ++i) {
            data_points.push_back( mba::make_array<double>(coords[i]) );
            data_values.push_back( data[i] );
    }
    
    m_MBA = new mba::cloud<1> (
                    mba::make_array<double>(min),
                    mba::make_array<double>(max),
                    data_points, data_values,
                    mba::make_array<size_t>(m_init_lattice_MBA));
    
    delete [] coords;
    delete [] data;

}

void OneDimData::printParameters() const
{
   if (procID()==0) {
      std::cout << "OneDimData grid function parameters:" << std::endl;
      std::cout << "  data_file: "   << m_data_file   << std::endl;
      std::cout << std::endl;
   }
}

#include "NamespaceFooter.H"
