#include "RZdata.H"

// Needed by external utils
using namespace std;

#include "rbf_interp_2d.hpp"
#include "r8lib.hpp"
#include "test_interp_2d.hpp"
#include "Directions.H"

#include "NamespaceHeader.H"

RZdata::RZdata( ParmParse& a_pp,
                const int& a_verbosity )
   : GridFunction(a_verbosity),
     m_subtype("MBA"),
     m_data_grid(NULL),
     m_weights(NULL),
     m_MBA(NULL),
     m_init_lattice_MBA(10)

{
   parseParameters( a_pp );
   if (m_subtype == "RBF") createInterpolationRBF();
   if (m_subtype == "MBA") createInterpolationMBA();
}

RZdata::~RZdata()
{
    if (m_data_grid) delete [] m_data_grid;
    if (m_weights) delete [] m_weights;
    if (m_MBA) delete m_MBA;
}

inline
void RZdata::parseParameters( ParmParse& a_pp )
{
   
   if ( a_pp.contains("data_file")) {
        a_pp.get("data_file", m_data_file);
   }
   else {
       MayDay::Error("RZdata:: No data file specified");
   }

   a_pp.query( "subtype", m_subtype );
    
   if (m_verbosity) {
      printParameters();
   }
}


void RZdata::setPointwise(FArrayBox&                 a_data,
                          const MultiBlockLevelGeom& a_geometry,
                          const FArrayBox&           a_real_coords,
                          const FArrayBox&           a_normalized_flux,
                          const int                  a_block_number) const
{
   const MagBlockCoordSys& coord_sys = getCoordSys(a_geometry, a_block_number);

   Box box( a_data.box() );
   FArrayBox cell_center_coords( box, SpaceDim );
   coord_sys.getCellCenteredRealCoords( cell_center_coords );

   if (m_subtype == "RBF") {
    
     //Get interpolated data
     int npoints = box.numPts();
     double *box_grid = new double[SpaceDim*npoints];
   
     int i=0;
     BoxIterator bit(box);
     for (bit.begin(); bit.ok(); ++bit) {
        IntVect iv = bit();
        box_grid[0+i*2] = cell_center_coords (iv,0);
        box_grid[1+i*2] = cell_center_coords (iv,1);
        i++;
     }
    
     //Choise of the radial basis functions are given by phi1, phi2, phi3, and phi4 (see rbf_interp_2d.cpp)
     //and should be consistent with phi in rbf_weights
     double *interpolated_data;
     interpolated_data = rbf_interp ( SpaceDim, m_data_npoints, m_data_grid,
                                   m_r0, phi1, m_weights, npoints, box_grid );

    
     //Copy interpolated data to a_data
     i=0;
     for (bit.begin(); bit.ok(); ++bit) {
        IntVect iv = bit();
        a_data(iv,0)=interpolated_data[i];
        i++;
     }

     delete [] box_grid;
     delete [] interpolated_data;
  }

  if (m_subtype == "MBA") {

      BoxIterator bit(box);
      for (bit.begin(); bit.ok(); ++bit) {
          IntVect iv = bit();
          double R = cell_center_coords (iv,0);
          double Z = cell_center_coords (iv,1);
          std::array<double,2> this_point = mba::make_array<double>(R, Z);
          a_data(iv,0) = m_MBA->getInterpolatedValue(this_point);
          
      }
       
  }
    
}


void RZdata::createInterpolationRBF()
{
    
#ifdef CH_MPI
    if (procID() == 0) {
#endif
     
     //Get number of data points lines (i.e., number of line in a data file)
     int npoints = 0;
     string line;
     ifstream datafile(m_data_file.c_str());
        
     while (std::getline(datafile, line))
        ++npoints;

        m_data_npoints = npoints;
#ifdef CH_MPI
    }
    
    MPI_Bcast(&m_data_npoints, 1, MPI_INT, 0, MPI_COMM_WORLD);

#endif
    
    m_data_grid = new double[SpaceDim*m_data_npoints];
    m_weights = new double[m_data_npoints];

#ifdef CH_MPI
    if (procID() == 0) {
#endif
    
     //Reading the table data
     double *R = new double[m_data_npoints];
     double *Z = new double[m_data_npoints];
     double *Data = new double[m_data_npoints];
        
     ifstream inFile;
     inFile.open( m_data_file.c_str() );
   
     for (int i=0; i<m_data_npoints; ++i) {
        inFile >> R[i];
        inFile >> Z[i];
        inFile >> Data[i];
     }
     inFile.close();
        
     for ( int i = 0; i < m_data_npoints; i++ )
     {
        m_data_grid[0+i*2] = R[i];
        m_data_grid[1+i*2] = Z[i];
     }
        
     double rmax = r8vec_max ( m_data_npoints, R );
     double rmin = r8vec_min ( m_data_npoints, R );
     double zmax = r8vec_max ( m_data_npoints, Z );
     double zmin = r8vec_min ( m_data_npoints, Z );
     double volume = ( rmax - rmin ) * ( zmax - zmin );
        
     double e = 1.0 / ( double ) ( SpaceDim );
     m_r0 = pow ( volume / m_data_npoints, e );
        
     //Choise of the radial basis functions are given by phi1, phi2, phi3, and phi4 (see rbf_interp_2d.cpp)
     //and should be consistent with phi in rbf_interp
     m_weights = rbf_weight ( SpaceDim, m_data_npoints, m_data_grid, m_r0, phi1, Data );
        
     delete [] R;
     delete [] Z;
     delete [] Data;

#ifdef CH_MPI
    }

    
    MPI_Bcast(&m_r0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(m_data_grid, SpaceDim*m_data_npoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(m_weights, m_data_npoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#endif
    
}

void RZdata::createInterpolationMBA()
{
    
#ifdef CH_MPI
    if (procID() == 0) {
#endif
        
        //Get number of data points lines (i.e., number of line in a data file)
        int npoints = 0;
        string line;
        ifstream datafile(m_data_file.c_str());
        
        while (std::getline(datafile, line))
            ++npoints;
        
        m_data_npoints = npoints;
#ifdef CH_MPI
    }
    
    MPI_Bcast(&m_data_npoints, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
#endif
    
    double *R = new double[m_data_npoints];
    double *Z = new double[m_data_npoints];
    double *Data = new double[m_data_npoints];
    
#ifdef CH_MPI
    if (procID() == 0) {
#endif
        
        //Reading the table data
        ifstream inFile;
        inFile.open( m_data_file.c_str() );
        
        for (int i=0; i<m_data_npoints; ++i) {
            inFile >> R[i];
            inFile >> Z[i];
            inFile >> Data[i];
        }
        inFile.close();
        
#ifdef CH_MPI
    }
    
    
    MPI_Bcast(R, m_data_npoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(Z, m_data_npoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(Data, m_data_npoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
#endif
    
    double rmax = r8vec_max ( m_data_npoints, R );
    double rmin = r8vec_min ( m_data_npoints, R );
    double zmax = r8vec_max ( m_data_npoints, Z );
    double zmin = r8vec_min ( m_data_npoints, Z );
    
    //Creating the MBA interpolant
    std::vector< std::array<double,2> > data_points ;
    std::vector<double> data_values;
    
    for(int i = 0; i < m_data_npoints; ++i) {
            data_points.push_back( mba::make_array<double>(R[i], Z[i]) );
            data_values.push_back( Data[i] );
    }
    
    m_MBA = new mba::cloud<2> (
                    mba::make_array<double>(rmin, zmin),
                    mba::make_array<double>(rmax, zmax),
                    data_points, data_values,
                    mba::make_array<size_t>(m_init_lattice_MBA, m_init_lattice_MBA));
    
    delete [] R;
    delete [] Z;
    delete [] Data;

}

void RZdata::printParameters() const
{
   if (procID()==0) {
      std::cout << "RZdata grid function parameters:" << std::endl;
      std::cout << "  data_file: "   << m_data_file   << std::endl;
      std::cout << std::endl;
   }
}

#include "NamespaceFooter.H"
