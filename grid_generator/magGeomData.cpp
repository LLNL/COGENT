#define CH_SPACEDIM 2

#include "BoxIterator.H"

#include "FieldData.H"

#include "CH_HDF5.H"
#include "parstream.H"
#ifdef CH_MPI
#include "CH_Attach.H"
#endif


#include "UsingNamespace.H"


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

   string RZ_mapping_file_name;
   string field_file_name;
   string block_name;

   if (a_argc==4) {
      RZ_mapping_file_name = a_argv[1];
      field_file_name = a_argv[2];
      block_name = a_argv[3];
   }
   else {
      cout << "Usage:  magGeomData...ex <RZ_mapping_file> [<field_file>] <blockname>," << endl;
      cout << "        where <RZ_mapping_file> defines the coordinate mapping, and" << endl;
      cout << "        where <field_file> defines the magnetic field, and" << endl;
      cout << "        where <blockname> is one of lcore, rcore, lcsol, rcsol, lsol, rsol, lpf, rpf, mcore, mcsol" << endl;
      return -1;
   }

   // Construct the field data object
   FieldData field_data(RZ_mapping_file_name, field_file_name, block_name);

   /*
     Fill an FArrayBox with physical coordinates.  In the
     new mesh generator code using the Brackbill-Saltzmann
     approach, the Euler equation will solved for these
     coordinates.  Here, however, we simply create a grid
     in the unit square, which is the assumed mapped
     domain used internally in the FieldData object and
     get the corresponding physical coordinates.
   */

   // Get a grid on the unit square.  The argument sets the size.
   FArrayBox xi;
   getUnitSquareGrid(10, 50, xi);

   const Box& box(xi.box());

   // Get the corresponding physical coordinates.
   FArrayBox physical_coordinates(box,2);
   field_data.getPhysicalCoordinates(xi, physical_coordinates);

   // Write a file containing the physical coordinates that
   // can be plotted with the plot_coords.m script using
   // either MatLab or Octave.  The file name is "coordsX",
   // where X is the block number.
   field_data.writePhysicalCoordinates(physical_coordinates);

   FArrayBox field_unit_vector(box,6);

#ifdef USE_DCT_FIELD  

   // Using the discrete cosine transform (DCT) field representation,
   // whose coefficients are contained in "field_file_name", get the field
   // unit vector (first two components) and the derivatives of its components
   // with respect to the physical coordinates (next 4 components) at
   // this set of physical coordinates.  The numbering of the 6 components
   // is explained at the top of FieldData::convertToMappedDerivatives(),
   // which converts the derivatives in physical space to derivatives
   // in mapped space.

   field_data.getFieldUnitVectorFromDCT(physical_coordinates, field_unit_vector);
   field_data.convertToMappedDerivatives(xi, field_unit_vector);

   FArrayBox magnetic_flux(physical_coordinates.box(),1);
   field_data.getMagneticFluxFromDCT(physical_coordinates, magnetic_flux);

   field_data.writeMagneticFlux(physical_coordinates, magnetic_flux);
#else

   // Using the field data contained in the coordinate mapping
   // file "geometry_file_name" get the field unit vector (first two
   // components) and the derivatives of its components with
   // respect to the mapped coordinates (next 4 components) at
   // this set of physical coordinates.  The numbering of the 6 components
   // is explained at the top of FieldData::getFieldUnitVectorFromMapping().

   field_data.getFieldUnitVectorFromMappingFile(physical_coordinates, field_unit_vector);

#endif

   // Write a file containing the field unit vector that
   // can be plotted with the plot_vectors.m script using
   // either MatLab or Octave.  The file name is "vectorsX",
   // where X is the block number.
   field_data.writeVectors(physical_coordinates, field_unit_vector);

#ifdef CH_MPI
   MPI_Finalize();
#endif

   return 0;
}
