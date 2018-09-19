#include <fstream>

#include "FArrayBox.H"
#include "RealVect.H"

#include "UsingNamespace.H"

using namespace std;

string block_names[10] = {"lcore", "rcore", "lcsol", "rcsol", "lsol", "rsol", "lpf", "rpf", "mcore", "mcsol"};

extern "C" {
  void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
}


double polynomial_interp( const int      npts,
                          const double*  nodes,
                          const double*  node_values,
                          const double   evaluation_point )
{
   double value = 0.;

   for (int n=0; n<npts; ++n) {
      double numer = node_values[n];
      double denom = 1.;
      for (int m=0; m<npts; ++m) {
         if ( m != n ) {
            numer *= (evaluation_point - nodes[m]);
            denom *= (nodes[n] - nodes[m]);
         }
      }
      value += numer / denom;
   }

   return value;
}


void extrapolate_radial( const int   side,
                         const int   poloidal_index,
                         const int   n_extend,
                         const int   degree,
                         FArrayBox&  phys_coords )
{
   CH_assert(side == -1 || side == 1);

   int npts = degree + 1;
   double* arc_length = new double[npts];
   double* nodes = new double[npts];
   double* R = new double[npts];
   double* Z = new double[npts];

   const Box& box = phys_coords.box();

   int last_valid;
   Box fill_box;

   if ( side == 1 ) {
      last_valid = box.bigEnd(0)-n_extend;
      fill_box.define(IntVect(last_valid+1,poloidal_index),IntVect(box.bigEnd(0),poloidal_index));
   }
   else {
      last_valid = box.smallEnd(0)+n_extend;
      fill_box.define(IntVect(box.smallEnd(0),poloidal_index),IntVect(last_valid-1,poloidal_index));
   }

   arc_length[0] = 0.;
   int m = last_valid;
   nodes[0] = double(m);
   R[0] = phys_coords(IntVect(m,poloidal_index),0);
   Z[0] = phys_coords(IntVect(m,poloidal_index),1);
   for (int n=1; n<npts; ++n) {
      m -= side;
      IntVect iv = IntVect(m,poloidal_index);

      R[n] = phys_coords(iv,0);
      Z[n] = phys_coords(iv,1);

      arc_length[n] = arc_length[n-1]
         + sqrt(pow(R[n] - R[n-1],2) + pow(Z[n] - Z[n-1],2));
      nodes[n] = double(m);
   }

   for (BoxIterator bit(fill_box); bit.ok(); ++bit) {
      IntVect iv = bit();

      // Extrapolate the arc length
      double this_arc_length_coord = polynomial_interp(npts, nodes, arc_length, double(iv[0]));

      // Interpolate the coordinates
      phys_coords(iv,0) = polynomial_interp(npts, arc_length, R, this_arc_length_coord );
      phys_coords(iv,1) = polynomial_interp(npts, arc_length, Z, this_arc_length_coord );
   }

   delete [] arc_length;
   delete [] nodes;
   delete [] R;
   delete [] Z;
}


void extrapolate_poloidal( const int   side,
                           const int   radial_index,
                           const int   n_extend,
                           const int   degree,
                           FArrayBox&  phys_coords )
{
   CH_assert(side == -1 || side == 1);

   int npts = degree + 1;
   double* arc_length = new double[npts];
   double* nodes = new double[npts];
   double* R = new double[npts];
   double* Z = new double[npts];

   const Box& box = phys_coords.box();

   int last_valid;
   Box fill_box;

   if ( side == 1 ) {
      last_valid = box.bigEnd(1)-n_extend;
      fill_box.define(IntVect(radial_index,last_valid+1),IntVect(radial_index,box.bigEnd(1)));
   }
   else {
      last_valid = box.smallEnd(1)+n_extend;
      fill_box.define(IntVect(radial_index,box.smallEnd(1)),IntVect(radial_index,last_valid-1));
   }

   arc_length[0] = 0.;
   int m = last_valid;
   nodes[0] = double(m);
   R[0] = phys_coords(IntVect(radial_index,m),0);
   Z[0] = phys_coords(IntVect(radial_index,m),1);
   for (int n=1; n<npts; ++n) {
      m -= side;
      IntVect iv = IntVect(radial_index,m);

      R[n] = phys_coords(iv,0);
      Z[n] = phys_coords(iv,1);

      arc_length[n] = arc_length[n-1]
         + sqrt(pow(R[n] - R[n-1],2) + pow(Z[n] - Z[n-1],2));
      nodes[n] = double(m);
   }

   for (BoxIterator bit(fill_box); bit.ok(); ++bit) {
      IntVect iv = bit();

      // Extrapolate the arc length
      double this_arc_length_coord = polynomial_interp(npts, nodes, arc_length, double(iv[1]));

      // Interpolate the coordinates
      phys_coords(iv,0) = polynomial_interp(npts, arc_length, R, this_arc_length_coord );
      phys_coords(iv,1) = polynomial_interp(npts, arc_length, Z, this_arc_length_coord );
   }

   delete [] arc_length;
   delete [] nodes;
   delete [] R;
   delete [] Z;
}


void appendPoints( const IntVect&     start,
                   const IntVect&     stop,
                   const IntVect&     skip,
                   vector<IntVect>&   node,
                   vector<RealVect>&  node_data,
                   FArrayBox&         phys_coords)
{
   for (int j=start[1]; j<=stop[1]; j+=skip[1]) {
      for (int i=start[0]; i<=stop[0]; i+=skip[0]) {
         IntVect iv(i,j);
         node.push_back(iv);
         node_data.push_back(RealVect(phys_coords(iv,0),phys_coords(iv,1)));
      }
   }
}


void interpolate_RBF( const vector<IntVect>&   node,
                      const vector<RealVect>&  node_data,
                      FArrayBox&               data )
{
   int npts = node.size();
   CH_assert(node_data.size() == npts);

   double* mat = new double[npts*npts];

   int n = 0;
   for (int col=0; col<npts; ++col) {
      for (int row=0; row<npts; ++row) {
         IntVect row_node = node[row];
         IntVect col_node = node[col];
         mat[n++] = pow( pow(double(row_node[0] - col_node[0]),2)
                       + pow(double(row_node[1] - col_node[1]),2), 1.5);
      }
   }

   int nrhs = 2;
   double* rhs = new double[nrhs*npts];
   int* ipvt = new int[npts];
   int info;

   int m = 0;
   for (vector<RealVect>::const_iterator it=node_data.begin(); it != node_data.end(); ++it) {
      RealVect this_node_data = *it;

      rhs[m]      = this_node_data[0];
      rhs[m+npts] = this_node_data[1];
      m++;
   }

   dgesv_(&npts, &nrhs, mat, &npts, ipvt, rhs, &npts, &info);

   if ( info != 0 ) {
      cout << "dgesv returned = " << info << endl;
   }

   for (BoxIterator bit(data.box()); bit.ok(); ++bit) {
      IntVect iv = bit();
      
      data(iv,0) = 0.;
      data(iv,1) = 0.;

      int n = 0;
      for (vector<IntVect>::const_iterator it=node.begin(); it != node.end(); ++it) {
         IntVect this_node = *it;

         double r = pow( pow(double(iv[0] - this_node[0]),2)
                       + pow(double(iv[1] - this_node[1]),2), 1.5);
         data(iv,0) += rhs[n] * r;
         data(iv,1) += rhs[n+npts] * r;

         n++;
      }
   }

   delete [] ipvt;
   delete [] rhs;
   delete [] mat;
}      


void checkForDuplicates( const vector<IntVect>& node )
{
   int npts = node.size();

   for (int n=0; n<npts-1; ++n) {
      IntVect this_node = node[n];
      for (int m=n+1; m<npts; ++m) {
         IntVect other_node = node[m];
         if ( other_node == this_node ) {
            cout << "Duplicate index: (" << this_node << endl;
         }
      }
   }
}


void interpConformingExtendedGrid( const string&     block_name,
                                   FArrayBox&        phys_coords,
                                   const IntVect&    n,
                                   const IntVect&    n_extend,
                                   const IntVect&    extrap_degree,
                                   const string&     option,
                                   const bool&       plot_blocks )
{
   int nr = n[0];
   int nr_extend = n_extend[0];
   int np = n[1];
   int np_extend = n_extend[1];
   int nr_total = nr + 2*nr_extend;
   int np_total = np + 2*np_extend;
   int radial_index, poloidal_index;
   
   // These parameters enable skipping of points added to the RBF interpolation point set to try
   // to limit the size of the resulting RBF dense matrices.
   IntVect skip = IntVect::Unit;

   vector<IntVect> node;
   vector<RealVect> node_data;

   if ( block_name == "mcore" || block_name == "mcsol" ) {

      // Just copy the valid and ghost points.  Even though we're not doing RBF interpolation, we build
      // the point list anyway in case we want to plot the would be interpolation points.

      appendPoints( IntVect::Zero, IntVect(nr_total-1,np_total-1), skip, node, node_data, phys_coords);
   }
   else {

      // STEP 1: Initialize the RBF interpolation point list with the valid points

      appendPoints( IntVect(nr_extend,np_extend), IntVect(nr+nr_extend-1,np+np_extend-1), skip, node, node_data, phys_coords);

      // STEP 2: Add to the RBF interpolation point list the original data at radial block boundaries not containing the X point.

      if ( option == "exclude_r_upper_p_lower" ) {
         appendPoints( IntVect(nr+nr_extend,np_extend), IntVect(nr_total-1,np_total-1), skip, node, node_data, phys_coords);
      }
      else {
         appendPoints( IntVect(nr+nr_extend,0), IntVect(nr_total-1,np_total-1), skip, node, node_data, phys_coords);
      }

      // STEP 3: Add to the RBF interpolation point list the original ghost data at poloidal block boundaries
      // not containing the X point

      appendPoints( IntVect(0, np+np_extend), IntVect(nr+nr_extend-1,np_total-1), skip, node, node_data, phys_coords);

      // STEP 4: Extrapolate or copy the antennae along physical boundaries at the poloidal
      // block boundary containing the X point and add the points to the PBF point list

      radial_index = nr+nr_extend-1;
      if ( option == "exclude_r_upper_p_lower" ) {
         extrapolate_poloidal(-1, radial_index, np_extend, extrap_degree[1], phys_coords);
      }
      appendPoints( IntVect(radial_index,0), IntVect(radial_index,np_extend-1), skip, node, node_data, phys_coords);

      // STEP 5: Extrapolate the two antennae at the separatrix block boundary

      poloidal_index = np_extend;
      extrapolate_radial(-1, poloidal_index, nr_extend, extrap_degree[0], phys_coords);
      appendPoints( IntVect(0,poloidal_index), IntVect(nr_extend-1,poloidal_index), IntVect::Unit, node, node_data, phys_coords);

      poloidal_index = np+np_extend-1;
      extrapolate_radial(-1, poloidal_index, nr_extend, extrap_degree[0], phys_coords);
      appendPoints( IntVect(0, poloidal_index), IntVect(nr_extend-1,poloidal_index), IntVect::Unit, node, node_data, phys_coords);

      // STEP 6: Extrapolate the remaining antennae at the X point in the poloidal direction
      // Extrapolate an antenna at the low poloidal boundary

      radial_index = nr_extend;
      extrapolate_poloidal(-1, radial_index, np_extend, extrap_degree[1], phys_coords);
      appendPoints( IntVect(radial_index,0), IntVect(radial_index,np_extend-1), IntVect::Unit, node, node_data, phys_coords);

      // STEP 7: All of the RBF interpolation points have now been set. Check the list for duplicates
      // (since the RBF matrices will be singular otherwise) and construct the RBF interpolators.

      checkForDuplicates(node);

      // STEP 8: Interpolate to the full index space

      interpolate_RBF(node, node_data, phys_coords );
   }

   if ( plot_blocks ) {
      FILE* f_block_id = fopen(block_name.c_str(),"w");

      fprintf(f_block_id, "%d %d\n", n[0] + 2*n_extend[0], n[1] + 2*n_extend[1]);
      for (BoxIterator bit(phys_coords.box()); bit.ok(); ++bit) {
         IntVect iv = bit();
         fprintf(f_block_id, "%20.13e %20.13e\n", phys_coords(iv,0), phys_coords(iv,1));
      }

      fclose(f_block_id);

      string file_name = block_name + "_scatter";
      FILE* f_scatter_block_id = fopen(file_name.c_str(),"w");
      
      for (vector<RealVect>::const_iterator it=node_data.begin(); it != node_data.end(); ++it) {
         RealVect this_node_data = *it;
         
         fprintf(f_scatter_block_id, "%20.13e %20.13e\n", this_node_data[0], this_node_data[1]);
      }

      fclose(f_scatter_block_id);
   }
}


void reorderBlock( const string&  block,
                   FArrayBox&     data )
{
   if ( block == "lcore" || block == "lcsol" || block == "rsol" || block == "rpf" ) {

      // Flip the columns

      const Box& box = data.box();
      FArrayBox tmp(box,2);

      for (BoxIterator bit(box); bit.ok(); ++bit) {
         IntVect iv(bit());
         IntVect iv_flip(iv[0],box.bigEnd(1)-iv[1]);

         for (int n=0; n<2; ++n) {
            tmp(iv,n) = data(iv_flip,n);
         }
      }
      data.copy(tmp);
   }

   if ( block == "lcore" || block == "rcore" || block == "lpf" || block == "rpf" )  {

      // Flip the rows

      const Box& box = data.box();
      FArrayBox tmp(box,2);

      for (BoxIterator bit(box); bit.ok(); ++bit) {
         IntVect iv(bit());
         IntVect iv_flip(box.bigEnd(0)-iv[0],iv[1]);

         for (int n=0; n<2; ++n) {
            tmp(iv,n) = data(iv_flip,n);
         }
      }
      data.copy(tmp);
   }
}


void add_block_ghosts( const string&  file_name,
                       const string&  block_name,
                       IntVect&       n,
                       IntVect&       n_extend,
                       string&        option,
                       bool&          plot_blocks,
                       FArrayBox&     phys_coords )
{
   ifstream inFile;

   inFile.open( file_name.c_str() );
   if (!inFile) {
      cout << "add_block_ghosts(): Unable to open input file " << file_name << endl;
      exit(1);
   }

   // Find the block
   bool found_block = false;

   for (int block_num=0; block_num<10; ++block_num) {
      string this_block_name;

      inFile >> this_block_name;
      inFile >> n[0];
      inFile >> n_extend[0];
      inFile >> n[1];
      inFile >> n_extend[1];

      phys_coords.define(Box(IntVect::Zero, IntVect(n[0] + 2*n_extend[0] - 1,n[1] + 2*n_extend[1] - 1)),2);

      for (BoxIterator bit(phys_coords.box()); bit.ok(); ++bit) {
         IntVect iv = bit();
         inFile >> phys_coords(iv,0);
         inFile >> phys_coords(iv,1);
         inFile.ignore(numeric_limits<streamsize>::max(),'\n');
      }

      found_block = (this_block_name == block_name);

      if (found_block) {
         break;
      }
   }

   inFile.close();

   if ( !found_block ) {
      cout << "Block " << block_name << " not found" << endl;
      exit(1);
   }

   cout << "Smoothing ghost cells for " << block_name << endl;

   // Reorder the indices in the blocks containing the X point so that increasing indices
   // correspond to an increasing distance from the X point
   reorderBlock(block_name, phys_coords);        

   IntVect extrap_degree(1,1);

   interpConformingExtendedGrid(block_name, phys_coords, n, n_extend, extrap_degree, option, plot_blocks);

   // Reorder the indices back to original
   reorderBlock(block_name, phys_coords);        
}


int main( int a_argc, char* a_argv[] )
{
   bool plot_blocks = false;
   string option;

   if (a_argc!=3) {
      cout << "Usage: smoothGhosts grid_file mapping_file" << endl;
      return -1;
   }
   
   string input_file_name = a_argv[1];
   string output_file_name = a_argv[2];
   CH_assert(input_file_name != output_file_name);
   
   FILE* f_out_id = fopen(output_file_name.c_str(),"w");

   for (int block_number = 0; block_number<10; ++block_number) {
      string block_name = block_names[block_number];
      IntVect n, n_extend;
      FArrayBox phys_coords;

      add_block_ghosts(input_file_name, block_name, n, n_extend, option, plot_blocks, phys_coords);

      fprintf(f_out_id, "%s %d %d %d %d\n", block_name.c_str(), n[0], n_extend[0], n[1], n_extend[1]);
      for (BoxIterator bit(phys_coords.box()); bit.ok(); ++bit) {
         IntVect iv = bit();
         fprintf(f_out_id, "%20.13e %20.13e %20.13e %20.13e\n", phys_coords(iv,0), phys_coords(iv,1), 0., 0.);
      }
   }

   fclose(f_out_id);

   return 0;
}


