#define CH_SPACEDIM 2

#include "BoxIterator.H"

#include "SFLocalField.H"

#include "CH_HDF5.H"
#include "parstream.H"
#ifdef CH_MPI
#include "CH_Attach.H"
#endif


#include "UsingNamespace.H"


inline int checkCommandLineArgs( int a_argc, char* a_argv[] )
{
   // Check for command line files
   if (a_argc<=4) {
      pout() << "Usage:  locSFGeom...ex <Rmaj> <zmax> <zm0> <blockname>," << endl;
      pout() << "        where Rmaj is major radius, zmax is max z extent used of local snowflake, " << endl;
      pout() << "        zm0 is z intercept of bounding flux surface at x=0, and" << endl;
      pout() << "        <blockname> is one of lcore, rcore, lcsol, rcsol, lsol, rsol, lpf, rpf" << endl;
      pout() << "An input line with nice inputs would be locSFGeom...ex 10.0 2.0 0.5 rcore " << endl;
      return -1;
   }
   return 0;
}


void
getUnitSquareGrid(const int  a_n,
                  FArrayBox& a_xi)
{
   Box box(IntVect::Zero, a_n*IntVect::Unit);
   double dx = 1./a_n;

   a_xi.define(box,2);

   for (BoxIterator bit(box);bit.ok();++bit) {
      IntVect iv = bit();
      for (int dir=0; dir<SpaceDim; ++dir) {
         a_xi(iv,dir) = iv[dir]*dx;
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

   int status = checkCommandLineArgs( a_argc, a_argv);

   if (status==0) {

      const double Rmaj(atof(a_argv[1]));
      const double zmax(atof(a_argv[2]));
      const double zm0(atof(a_argv[3]));
      string block_name(a_argv[4]);

      // Construct the field data object
      SFLocalField sf_localfield(Rmaj, zmax, zm0, block_name);

      /*
        Fill an FArrayBox with physical coordinates.  In the
        new mesh generator code using the Brackbill-Saltzmann
        approach, the Euler equation will solved for these
        coordinates.  Here, however, we simply create a grid
        in the unit square, which is the assumed mapped
        domain used internally in the SFLocalField object and
        get the corresponding physical coordinates.
      */

      // Get a grid on the unit square.  The argument sets the size.
      FArrayBox xi;
      getUnitSquareGrid(20, xi);  // doubled from original, but will only plot every other fieldline

      const Box& box(xi.box());

      // Get the corresponding physical coordinates.
      FArrayBox physical_coordinates(box,2);
      sf_localfield.getPhysicalCoordinates(xi, physical_coordinates);

      // Write a file containing the physical coordinates that
      // can be plotted with the plot_coords.m script using
      // either MatLab or Octave.  The file name is "coordsX",
      // where X is the block number.
      sf_localfield.writePhysicalCoordinates(physical_coordinates);

      // At this set of physical coordinates, get the field
      // unit vector (first two components) and the derivatives
      // of its components with respect to the physical coordinates
      // (next 4 components).  The numbering of the 6 components
      // is explained at the top of SFLocalField::getFieldUnitVector().
      FArrayBox field_unit_vector(box,6);
      sf_localfield.getFieldUnitVector(physical_coordinates, field_unit_vector);

      // Write a file containing the field unit vector that
      // can be plotted with the plot_vectors.m script using
      // either MatLab or Octave.  The file name is "vectorsX",
      // where X is the block number.
      sf_localfield.writeVectors(physical_coordinates, field_unit_vector);

   }

#ifdef CH_MPI
   MPI_Finalize();
#endif

   return status;
}
