#include "OneDimData.H"

// Needed by external utils
using namespace std;

#include "rbf_interp_2d.hpp"
#include "r8lib.hpp"
#include "test_interp_2d.hpp"

#include "Directions.H"
#include "SNCoreCoordSys.H"
#include "SingleNullCoordSys.H"
#include "SNCoreBlockCoordSys.H"
#include "SingleNullBlockCoordSys.H"
#include "SNCoreBlockCoordSysModel.H"
#include "SingleNullBlockCoordSysModel.H"

#include "NamespaceHeader.H"



OneDimData::OneDimData( ParmParse& a_pp,
                const int& a_verbosity )
   : GridFunction(a_verbosity),
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

void OneDimData::checkGeometryValidity( const MultiBlockLevelGeom& a_geometry ) const
{
  //const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
  //Presently the capbaility works for any geometry, but we might need this later
}

void OneDimData::setPointwise(FArrayBox&                 a_data,
                              const MultiBlockLevelGeom& a_geometry,
                              const FArrayBox&           a_real_coords,
                              const FArrayBox&           a_normalized_flux,
                              const int                  a_block_number) const
{
   const MagBlockCoordSys& coord_sys = getCoordSys(a_geometry, a_block_number);

   Box box( a_data.box() );
   FArrayBox cc_phys_coords( box, SpaceDim );
   coord_sys.getCellCenteredRealCoords( cc_phys_coords );

   FArrayBox cc_mapped_coords( box, SpaceDim );
   coord_sys.getCellCenteredMappedCoords( cc_mapped_coords );
   

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
         if (typeid(coord_sys) == typeid(SNCoreBlockCoordSys))  coord = ((const SNCoreBlockCoordSys&)coord_sys).getNormMagneticFlux(phys_coordinate);
         if (typeid(coord_sys) == typeid(SingleNullBlockCoordSys))  coord = ((const SingleNullBlockCoordSys&)coord_sys).getNormMagneticFlux(phys_coordinate);
         
         if (typeid(coord_sys) == typeid(SNCoreBlockCoordSysModel))  coord = ((const SNCoreBlockCoordSysModel&)coord_sys).getNormMagneticFlux(phys_coordinate);
         if (typeid(coord_sys) == typeid(SingleNullBlockCoordSysModel))  coord = ((const SingleNullBlockCoordSysModel&)coord_sys).getNormMagneticFlux(phys_coordinate);
      }
      
      if (m_coord_type == "outer_midplane") {
         
         if (typeid(coord_sys) == typeid(SNCoreBlockCoordSys)) {
            double fluxNorm  = ((const SNCoreBlockCoordSys&)coord_sys).getNormMagneticFlux(phys_coordinate);
            double Rsep  = ((const SNCoreBlockCoordSys&)coord_sys).getOuterRsep();
            coord = ((const SNCoreBlockCoordSys&)coord_sys).getOuterMidplaneCoord(fluxNorm) - Rsep;
            
         }
         
         if (typeid(coord_sys) == typeid(SingleNullBlockCoordSys)) {
            double fluxNorm  = ((const SingleNullBlockCoordSys&)coord_sys).getNormMagneticFlux(phys_coordinate);
            double Rsep  = ((const SingleNullBlockCoordSys&)coord_sys).getOuterRsep();
            coord = ((const SingleNullBlockCoordSys&)coord_sys).getOuterMidplaneCoord(fluxNorm) - Rsep;
         }
         
         if (typeid(coord_sys) == typeid(SNCoreBlockCoordSysModel)) {
            double fluxNorm  = ((const SNCoreBlockCoordSysModel&)coord_sys).getNormMagneticFlux(phys_coordinate);
            double Rsep  = ((const SNCoreBlockCoordSysModel&)coord_sys).getOuterRsep();
            coord = ((const SNCoreBlockCoordSysModel&)coord_sys).getOuterMidplaneCoord(fluxNorm) - Rsep;
            
         }
         
         if (typeid(coord_sys) == typeid(SingleNullBlockCoordSysModel)) {
            double fluxNorm  = ((const SingleNullBlockCoordSysModel&)coord_sys).getNormMagneticFlux(phys_coordinate);
            double Rsep  = ((const SingleNullBlockCoordSysModel&)coord_sys).getOuterRsep();
            coord = ((const SingleNullBlockCoordSysModel&)coord_sys).getOuterMidplaneCoord(fluxNorm) - Rsep;
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
        
     if( !datafile.is_open() ){
        stringstream error_msg;
        error_msg << "OneDimData::createInterpolationRBF(): File " << m_data_file.c_str() << " not found";
        MayDay::Error(error_msg.str().c_str());
     }
        
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

        if( !datafile.is_open() ){
           stringstream error_msg;
           error_msg << "OneDimData::createInterpolationMBA(): File " << m_data_file.c_str() << " not found";
           MayDay::Error(error_msg.str().c_str());
        }
        
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
