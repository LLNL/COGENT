#include <fstream>
#include <float.h>
#include "BlockMapping.H"
#include "BlockMappingF_F.H"
#include "BoxIterator.H"
#include "IntVect.H"
#include "CONSTANTS.H"

#include "NamespaceHeader.H"

enum SingleNullBlockType {LCORE,RCORE,LCSOL,RCSOL,LSOL,RSOL,LPF,RPF,MCORE,MCSOL,NUM_SINGLE_NULL_BLOCKS};
enum directions {RADIAL_DIR, POLOIDAL_DIR};

inline double signum(const double a_x)
{
   
   double result;
   
   if (a_x>0) {
      result = 1.0;
   }
   else if (a_x<0) {
      result = -1.0;
   }
   else {
      result = 0;
   }

   return result;
}


BlockMapping::BlockMapping(const string& a_field_file_name,
                           const int     a_block_number,
                           const double  a_psi_bigEnd,
                           const int     a_nExtrp,
                           const double  a_trans_length,
                           const double  a_blending_factor)
  : m_block_number(a_block_number),
    m_psi_smallEnd(1.0),
    m_psi_bigEnd(a_psi_bigEnd),
    m_nExtrp(a_nExtrp),
    m_tol(1.0e-6),
    m_dl(0.01),
    m_lengthCORE(1.0),
    m_dealigned_grid(false),
    m_straight_block_bndry(true),
    m_max_iterations(100),
    m_psi_expansion_order(-1),
    m_blending_factor(a_blending_factor),
    m_trans_length(a_trans_length)
{
   readFieldFromDCTFile(a_field_file_name);

}


BlockMapping::~BlockMapping()
{
}

void
BlockMapping::init(const double a_trans_rad)
{
   /* 
    Class initialization
    */
   
   m_physFluxOnAxis = getMagneticFluxFromDCT(m_magAxis);
   m_physFluxOnSep = getMagneticFluxFromDCT(m_Xpoint);
   
   m_sep_dir = getSeparatrixDirAnalytic();
   getCutDir(m_sep_dir, m_core_cut_dir, m_sol_cut_dir);
   
   getPsiLinearCoeff(m_lin_1, m_lin_2);
   
   if (a_trans_rad > 0) {
      m_dealigned_grid = true;
      m_trans_rad = a_trans_rad;
   }
    
}

void
BlockMapping::getPhysicalCoordinates(FArrayBox& a_physical_coordinates,
                                     const Vector<RealVect*> & a_mapping) const
{
   /*
    Compute the block mapping on valid cells 
    For a given nodal mapped coordinate compute its physical coordinate.
    NB:: the box type is cell-centered though
   */
 
   CH_TIME("BlockMapping::getPhysicalCoordinates");
   
   const Box& box(a_physical_coordinates.box());
   
   // Compute block cut coordinates and local directions
   Box cut_box(IntVect(0,0), IntVect(box.size(0)-1,0));
   m_block_cut_coords.define(cut_box,2);
   m_block_cut_dir.define(cut_box,2);
   initializeBlockCut(m_block_cut_coords, m_block_cut_dir);
   
   // To speed-up ray tracing we decompose the poloidal direction
   DisjointBoxLayout grids;
   Vector<Box> grid_boxes;
   int decomp_num;

   Box modified_box(box);
   modified_box.growHi(1, -m_nExtrp);
   
   getPoloidalDisjointBoxLayout(grids, grid_boxes, decomp_num, modified_box);
   
   LevelData<FArrayBox> phys_coords_decomposed(grids, 2, IntVect::Zero);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      // Create a box to store a single radial grid line
      Box poloidal_box(grids[dit]);
      poloidal_box.growHi(0, 1-poloidal_box.size()[RADIAL_DIR]);
      
      for (BoxIterator bit(poloidal_box); bit.ok(); ++bit) {
         
         IntVect iv = bit();
         
         int pol_index = iv[POLOIDAL_DIR];
         
         //Get poloidal coordinate (arc length)
         double length = getPoloidalArcLength(pol_index, "cosine");
      
         //Get a radial grid line for the given arc length coordinate
         Box box_radial_line(IntVect(0,pol_index), IntVect(box.size(RADIAL_DIR)-1,pol_index));
         FArrayBox radialGridLine(box_radial_line,2);
         getRadialGridLine(radialGridLine, a_mapping, length);
         phys_coords_decomposed[dit].copy(radialGridLine);
      }
   }
   
   // Assemble data back to a single FArrayBox for output purposes
   assembleDecomposedData(a_physical_coordinates, phys_coords_decomposed, grid_boxes, decomp_num);
   
   
   // Grid extrapolation (sometimes used to extend the grid in dthe ivertor leg region)
   // Because the number of extrapolated cells is typically small, do not
   // decompose/parallelize this calculation
   if (m_nExtrp>0) {
      Box box_extrp(IntVect(0,box.size(1) - m_nExtrp), IntVect(box.size(0)-1,box.size(1)-1));
      extrapolateGrid(a_physical_coordinates,box_extrp);
   }
}

void
BlockMapping::getPoloidalDisjointBoxLayout(DisjointBoxLayout&   a_grids,
                                           Vector<Box>&         a_boxes,
                                           int&                 a_decomp_num,
                                           const Box&           a_box) const
{
   int n_proc = numProc();
   int n_pol = a_box.size(POLOIDAL_DIR);
   
   a_decomp_num = (n_pol>n_proc) ? n_proc : n_pol;
   int n_loc_pol = (n_pol>n_proc) ? int(n_pol/n_proc) : 1;
   
   IntVect domain_hi_end = a_box.bigEnd();
   IntVect domain_lo_end = a_box.smallEnd();
      
   
   Vector<int> proc_ids;
   
   for(int i = 0; i < a_decomp_num; i++) {
      
      int pol_ind_lo = (domain_lo_end[1] + i * n_loc_pol);
      int pol_ind_hi = (domain_lo_end[1] + (i+1) * n_loc_pol - 1);
          
      if (i == a_decomp_num - 1) {
         pol_ind_hi = domain_hi_end[POLOIDAL_DIR];
      }

      IntVect lo_end(domain_lo_end[RADIAL_DIR], pol_ind_lo);
      IntVect hi_end(domain_hi_end[RADIAL_DIR], pol_ind_hi);
      
      Box box(lo_end, hi_end);
      int procID = i;

      a_boxes.push_back( box );
      proc_ids.push_back( procID );
      
   }
   a_grids.define( a_boxes, proc_ids );
}


void
BlockMapping::assembleDecomposedData(FArrayBox&                    a_data,
                                     const LevelData<FArrayBox>&   a_data_decomp,
                                     const Vector<Box>&            a_grid_boxes,
                                     const int&                    a_decomp_num ) const
{
   /*
    Puts decomposed LevelData<FArrayBox> onto a single FArrayBox visible by all procs
    We assume that each proc contains only one box (we created grids that way)
    Also, the processor IDs assigned to grids correspond to 0 .. a_decomp_num-1
   */
   
   CH_TIME("BlockMapping::assembleDecomposedData");
   
   FArrayBox this_data;
   
   Vector<FArrayBox> all_data(numProc());
   
   // Fill local data
   for (int ivec = 0; ivec < numProc(); ivec++) {
      if (procID() == ivec && ivec < a_decomp_num) {
         this_data.define(a_grid_boxes[ivec], 2);
         
         // Iterate over boxes of decomposed object to find the box needed
         for (DataIterator dit( a_data_decomp.dataIterator() ); dit.ok(); ++dit) {
            if (a_data_decomp[dit].box() == a_grid_boxes[ivec]) {
               this_data.copy(a_data_decomp[dit]);
            }
         }
      }
      else {
         if (procID() == ivec) {
            // Or perhas should use the empty constructor to create invalid FArrayBox
            Box empty_box;
            this_data.define(empty_box,1);
         }
      }
   }

#ifdef CH_MPI
   // Gather data across all processors and put it onto all processors
   for (int ivec = 0; ivec < numProc(); ivec++) {
      gather(all_data, this_data, ivec);
   }
   
   // Create global interpolation FArrayBox objects on all procs
   for (int ivec = 0; ivec < a_decomp_num; ivec++) {
      a_data.copy(all_data[ivec], a_grid_boxes[ivec]);
   }
#else
   a_data.copy(this_data, a_grid_boxes[0]);
#endif
}

void BlockMapping::readFieldFromDCTFile(const string& a_file_name)
{
   ifstream fieldCoefficientsFile;

   int NR, NZ;

#ifdef CH_MPI
   if (procID() == 0) {
#endif

      fieldCoefficientsFile.open( a_file_name.c_str() );

      if (!fieldCoefficientsFile) {
         cout << "Unable to open field DCT coefficients file" << endl;
         exit(1);
      }

     //        cout << "Reading field coefficients from " << a_field_file_name.c_str() << endl;

      fieldCoefficientsFile >> m_Rmin;
      fieldCoefficientsFile >> m_Rmax;
      fieldCoefficientsFile >> m_Zmin;
      fieldCoefficientsFile >> m_Zmax;
      fieldCoefficientsFile >> NR;
      fieldCoefficientsFile >> NZ;
      fieldCoefficientsFile >> m_magAxis[0];
      fieldCoefficientsFile >> m_magAxis[1];
      fieldCoefficientsFile >> m_Xpoint[0];
      fieldCoefficientsFile >> m_Xpoint[1];

#ifdef CH_MPI
   }
#endif

#ifdef CH_MPI
   MPI_Bcast(&m_Rmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast(&m_Rmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast(&m_Zmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast(&m_Zmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast(&NR, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&NZ, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
   Box kbox = Box(IntVect::Zero,IntVect(NR-1,NZ-1));

   m_psi_coefs.define(kbox,1);

#ifdef CH_MPI
   if (procID() == 0) {
#endif
      for (BoxIterator bit(kbox); bit.ok(); ++bit) {
            fieldCoefficientsFile >> m_psi_coefs(bit());
      }
        
#ifdef CH_MPI
   }
#endif // end if MPI

#ifdef CH_MPI
   int num_values = kbox.numPts();
   MPI_Bcast(m_psi_coefs.dataPtr(), num_values, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

#ifdef CH_MPI
   if (procID() == 0) {
#endif
   fieldCoefficientsFile.close();
#ifdef CH_MPI
   }
#endif
   
   
   if ((m_psi_expansion_order != -1) && ((m_psi_expansion_order > NR || m_psi_expansion_order > NZ))) {
      MayDay::Error("readFieldFromDCTFile: Insufficient number of DCT coefficients for requested psi_expansion_order");
   }
   
   if (m_psi_expansion_order == -1) {
      m_psi_expansion_order = NR;
      if (NZ < NR) m_psi_expansion_order = NZ;
   }
 
}


double BlockMapping::getNormMagneticFlux( const RealVect& a_physical_coordinate ) const
{
   /*
    Return the normalized magnetic flux.
    psiNorm = 0 at the magnetic axis
    psiNorm = 1 at the separatrix
    */
   

   double physFlux = getMagneticFluxFromDCT(a_physical_coordinate);

   return (physFlux - m_physFluxOnAxis)/(m_physFluxOnSep - m_physFluxOnAxis);

}

void
BlockMapping::getNormFluxCoordinates(FArrayBox&  a_magFlux_coords) const
{
   /*
    Return the normalized magnetic flux values
    on the grid nodes.
    */
   
   CH_assert(a_magFlux_coords.nComp() == 1);
   const Box& box(a_magFlux_coords.box());
   int n_psi = box.size(0) - 1;
   double dpsi = (m_psi_bigEnd - m_psi_smallEnd)/n_psi;
   
   for (BoxIterator bit(box);bit.ok();++bit) {
      IntVect iv = bit();
      a_magFlux_coords(iv,0) = m_psi_smallEnd + float(iv[0]) * dpsi;
   }
   
}

double
BlockMapping::getMagneticFluxFromDCT( const RealVect& a_physical_coordinate ) const
{
   /*
     Given the input physical coordinate a_physical_coordinate, return the magnetic flux
     at that point by evaluating its DCT expansion.
   */

   int NR = m_psi_coefs.box().size(0);
   int NZ = m_psi_coefs.box().size(1);

   int llen = NR>NZ? NR: NZ;

   double * temp = new double[7*llen];
   double * lambda = temp;
   double * facR = lambda + llen;
   double * facZ = facR + llen;
   double * sinfacR = facZ + llen;
   double * cosfacR = sinfacR + llen;
   double * sinfacZ = cosfacR + llen;
   double * cosfacZ = sinfacZ + llen;

   lambda[0] = 1. / sqrt(2.);
   for (int l=1; l<llen; ++l) {
      lambda[l] = 1.;
   }

   for (int i=0; i<NR; ++i) {
      facR[i] = i * Pi / NR;
   }

   for (int j=0; j<NZ; ++j) {
      facZ[j] = j * Pi / NZ;
   }

   double Rscale = (NR-1)/(m_Rmax - m_Rmin);
   double Zscale = (NZ-1)/(m_Zmax - m_Zmin);

   double psi;

   double sR = (a_physical_coordinate[0] - m_Rmin) * Rscale;
   double sZ = (a_physical_coordinate[1] - m_Zmin) * Zscale;

   for (int i=0; i<NR; ++i) {
      double t = facR[i] * (sR + 0.5);
      sinfacR[i] = sin(t);
      cosfacR[i] = cos(t);
   }

   for (int j=0; j<NZ; ++j) {
      double t = facZ[j] * (sZ + 0.5);
      sinfacZ[j] = sin(t);
      cosfacZ[j] = cos(t);
   }

   // Compute psi

   int derivR = 0;
   int derivZ = 0;

   FORT_DCT_INTERP( CHF_CONST_FRA1(m_psi_coefs,0),
                    CHF_CONST_INT(m_psi_expansion_order),
                    CHF_CONST_INT(derivR),
                    CHF_CONST_INT(derivZ),
                    CHF_R1D(facR,NR),
                    CHF_R1D(facZ,NZ),
                    CHF_R1D(sinfacR,NR),
                    CHF_R1D(cosfacR,NR),
                    CHF_R1D(sinfacZ,NZ),
                    CHF_R1D(cosfacZ,NZ),
                    CHF_R1D(lambda,llen),
                    CHF_REAL(psi) );

   delete[] temp;

   if (m_dealigned_grid) {

      double origFlux = psi - m_physFluxOnSep;
   
      double x_new = m_lin_1[0] * (a_physical_coordinate[0] - m_Xpoint[0])
                   + m_lin_1[1] * (a_physical_coordinate[1] - m_Xpoint[1]);
   
      double y_new = m_lin_2[0] * (a_physical_coordinate[0] - m_Xpoint[0])
                   + m_lin_2[1] * (a_physical_coordinate[1] - m_Xpoint[1]);
   
      double linFlux = m_blending_factor * m_trans_rad * (fabs(x_new) - fabs(y_new));
   
      RealVect dist_to_X(a_physical_coordinate);
      dist_to_X -= m_Xpoint;
      dist_to_X /= m_trans_rad;
      double XtoXpt = sqrt(pow(dist_to_X[0],2) + pow(dist_to_X[1],2));
   
      double blend = tanh(XtoXpt);
      psi = blend * origFlux + (1.0 - blend) * linFlux + m_physFluxOnSep;
   }
   
   return psi;
}

Vector<Real>
BlockMapping::getMagFieldData( const RealVect& a_physical_coordinate ) const
{
   /*
    Given the RealVect a_physical_coordinate of physical coordinates, fill the
    return Vector<Real> with the following components:
    
    Component 0: First component of B (i.e., B_R)
    Component 1: Second component of B (i.e., B_Z)
    Component 2: First component of gradPsi (dPsi/dR)
    Component 3: Second component of gradPsi (dPsi/dZ)
    Component 4: d^2Psi/dR^2
    Component 5: d^2Psi/dRdZ
    Component 6: d^2Psi/dZ^2
    
    where B is the field vector, and Psi is the unnormalized magnetic flux.
    */
   
   Vector<Real> result(7);
   
   double RB[2];
   double gradPsi[2];
   double secondDer[3];
   
   int NR = m_psi_coefs.box().size(0);
   int NZ = m_psi_coefs.box().size(1);
   
   int llen = NR>NZ? NR: NZ;
   
   double * temp = new double[7*llen];
   double * lambda = temp;
   double * facR = lambda + llen;
   double * facZ = facR + llen;
   double * sinfacR = facZ + llen;
   double * cosfacR = sinfacR + llen;
   double * sinfacZ = cosfacR + llen;
   double * cosfacZ = sinfacZ + llen;
   
   lambda[0] = 1. / sqrt(2.);
   for (int l=1; l<llen; ++l) {
      lambda[l] = 1.;
   }
   
   for (int i=0; i<NR; ++i) {
      facR[i] = i * Pi / NR;
   }
   
   for (int j=0; j<NZ; ++j) {
      facZ[j] = j * Pi / NZ;
   }
   
   double Rscale = (NR-1)/(m_Rmax - m_Rmin);
   double Zscale = (NZ-1)/(m_Zmax - m_Zmin);
   
   
   double sR = (a_physical_coordinate[0] - m_Rmin) * Rscale;
   double sZ = (a_physical_coordinate[1] - m_Zmin) * Zscale;
   
   for (int i=0; i<NR; ++i) {
      double t = facR[i] * (sR + 0.5);
      sinfacR[i] = sin(t);
      cosfacR[i] = cos(t);
   }
   
   for (int j=0; j<NZ; ++j) {
      double t = facZ[j] * (sZ + 0.5);
      sinfacZ[j] = sin(t);
      cosfacZ[j] = cos(t);
   }
   
   double value;
   
   // Compute dpsi_dR
   
   int derivR = 1;
   int derivZ = 0;
   
   FORT_DCT_INTERP(  CHF_CONST_FRA1(m_psi_coefs,0),
                   CHF_CONST_INT(m_psi_expansion_order),
                   CHF_CONST_INT(derivR),
                   CHF_CONST_INT(derivZ),
                   CHF_R1D(facR,NR),
                   CHF_R1D(facZ,NZ),
                   CHF_R1D(sinfacR,NR),
                   CHF_R1D(cosfacR,NR),
                   CHF_R1D(sinfacZ,NZ),
                   CHF_R1D(cosfacZ,NZ),
                   CHF_R1D(lambda,llen),
                   CHF_REAL(value) );
   
   value *= Rscale;
   
   RB[1] = value;  // R*B_Z
   gradPsi[0] = value;
   
   
   
   // Compute dpsi_dZ
   
   derivR = 0;
   derivZ = 1;
   
   FORT_DCT_INTERP( CHF_CONST_FRA1(m_psi_coefs,0),
                   CHF_CONST_INT(m_psi_expansion_order),
                   CHF_CONST_INT(derivR),
                   CHF_CONST_INT(derivZ),
                   CHF_R1D(facR,NR),
                   CHF_R1D(facZ,NZ),
                   CHF_R1D(sinfacR,NR),
                   CHF_R1D(cosfacR,NR),
                   CHF_R1D(sinfacZ,NZ),
                   CHF_R1D(cosfacZ,NZ),
                   CHF_R1D(lambda,llen),
                   CHF_REAL(value) );
   
   value *= Zscale;
   
   RB[0] = -value;  // R*B_R
   gradPsi[1] = value;
   
   
   
   // Compute d2psi_dR2
   
   derivR = 2;
   derivZ = 0;
   
   FORT_DCT_INTERP(CHF_CONST_FRA1(m_psi_coefs,0),
                   CHF_CONST_INT(m_psi_expansion_order),
                   CHF_CONST_INT(derivR),
                   CHF_CONST_INT(derivZ),
                   CHF_R1D(facR,NR),
                   CHF_R1D(facZ,NZ),
                   CHF_R1D(sinfacR,NR),
                   CHF_R1D(cosfacR,NR),
                   CHF_R1D(sinfacZ,NZ),
                   CHF_R1D(cosfacZ,NZ),
                   CHF_R1D(lambda,llen),
                   CHF_REAL(value) );
   
   value *= (Rscale * Rscale); // d(R*BZ)/dR
   secondDer[0] = value;
   
   // Compute d2psi_dZ2
   
   derivR = 0;
   derivZ = 2;
   
   FORT_DCT_INTERP(CHF_CONST_FRA1(m_psi_coefs,0),
                   CHF_CONST_INT(m_psi_expansion_order),
                   CHF_CONST_INT(derivR),
                   CHF_CONST_INT(derivZ),
                   CHF_R1D(facR,NR),
                   CHF_R1D(facZ,NZ),
                   CHF_R1D(sinfacR,NR),
                   CHF_R1D(cosfacR,NR),
                   CHF_R1D(sinfacZ,NZ),
                   CHF_R1D(cosfacZ,NZ),
                   CHF_R1D(lambda,llen),
                   CHF_REAL(value) );
   
   value *= (Zscale * Zscale); // -d(R*B_R)/dZ
   secondDer[2] = value;
   
   // Compute d2psi_dRdZ
   
   derivR = 1;
   derivZ = 1;
   
   FORT_DCT_INTERP(CHF_CONST_FRA1(m_psi_coefs,0),
                   CHF_CONST_INT(m_psi_expansion_order),
                   CHF_CONST_INT(derivR),
                   CHF_CONST_INT(derivZ),
                   CHF_R1D(facR,NR),
                   CHF_R1D(facZ,NZ),
                   CHF_R1D(sinfacR,NR),
                   CHF_R1D(cosfacR,NR),
                   CHF_R1D(sinfacZ,NZ),
                   CHF_R1D(cosfacZ,NZ),
                   CHF_R1D(lambda,llen),
                   CHF_REAL(value) );
   
   value *= (Rscale * Zscale);
   secondDer[1] = value;
   


   // Compute psi                                                                                                                                                                                 
   double psi;
   derivR = 0;
   derivZ = 0;

   FORT_DCT_INTERP( CHF_CONST_FRA1(m_psi_coefs,0),
		    CHF_CONST_INT(m_psi_expansion_order),
		    CHF_CONST_INT(derivR),
		    CHF_CONST_INT(derivZ),
		    CHF_R1D(facR,NR),
		    CHF_R1D(facZ,NZ),
		    CHF_R1D(sinfacR,NR),
		    CHF_R1D(cosfacR,NR),
		    CHF_R1D(sinfacZ,NZ),
		    CHF_R1D(cosfacZ,NZ),
		    CHF_R1D(lambda,llen),
		    CHF_REAL(psi) );


   delete[] temp;
   
   double R = a_physical_coordinate[0];
   
   result[0] = RB[0] / R;
   result[1] = RB[1] / R;
   
   result[2] = gradPsi[0];
   result[3] = gradPsi[1];
   
   result[4] = secondDer[0];
   result[5] = secondDer[1];
   result[6] = secondDer[2];

   
   if (m_dealigned_grid) {
      
      //Blend original (DCT) flux and the rectalinear flux function near the X-point

      double origFlux = psi - m_physFluxOnSep;

      double x_new = m_lin_1[0] * (a_physical_coordinate[0] - m_Xpoint[0])
                   + m_lin_1[1] * (a_physical_coordinate[1] - m_Xpoint[1]);
      
      double y_new = m_lin_2[0] * (a_physical_coordinate[0] - m_Xpoint[0])
                   + m_lin_2[1] * (a_physical_coordinate[1] - m_Xpoint[1]);

      double linFlux = m_blending_factor * m_trans_rad * (fabs(x_new) - fabs(y_new));

      double linFlux_x = m_blending_factor * m_trans_rad * (m_lin_1[0] * signum(x_new) - m_lin_2[0] * signum(y_new));
      double linFlux_y = m_blending_factor * m_trans_rad * (m_lin_1[1] * signum(x_new) - m_lin_2[1] * signum(y_new));
      
      RealVect dist_to_X(a_physical_coordinate);
      dist_to_X -= m_Xpoint;
      dist_to_X /= m_trans_rad;
      double XtoXpt = sqrt(pow(dist_to_X[0],2.0) + pow(dist_to_X[1],2.0));
      
      double blend = tanh(XtoXpt);
      
      //Now compute spatial derivatives of the physical (unnormalzied) flux, given by
      //physFlux = blend * origFlux + (1.0 - blend) * linFlux + m_physFluxOnSep;

      double blend_r = 1.0/pow(cosh(XtoXpt),2.0);
      double blend_rr = -2.0 * blend * blend_r;

      double XtoXpt_x = (1.0/m_trans_rad) * dist_to_X[0]/XtoXpt;
      double XtoXpt_y = (1.0/m_trans_rad) * dist_to_X[1]/XtoXpt;

      double XtoXpt_xx = 1.0/pow(m_trans_rad,2.0) * (1.0/XtoXpt - dist_to_X[0]*dist_to_X[0]/pow(XtoXpt,3.0));
      double XtoXpt_yy = 1.0/pow(m_trans_rad,2.0) * (1.0/XtoXpt - dist_to_X[1]*dist_to_X[1]/pow(XtoXpt,3.0));
      double XtoXpt_xy = 1.0/pow(m_trans_rad,2.0) * (- dist_to_X[0]*dist_to_X[1]/pow(XtoXpt,3.0));

    
      result[2] = blend_r * XtoXpt_x * (origFlux - linFlux) + blend * gradPsi[0] + (1.0 - blend) * linFlux_x;
      result[3] = blend_r * XtoXpt_y * (origFlux - linFlux) + blend * gradPsi[1] + (1.0 - blend) * linFlux_y;
                                            
      result[4] = (blend_rr * pow(XtoXpt_x,2.0) + blend_r * XtoXpt_xx) * (origFlux - linFlux)
                  + 2.0 * blend_r * XtoXpt_x * (gradPsi[0] - linFlux_x) + blend * secondDer[0];

      result[6] = (blend_rr * pow(XtoXpt_y,2.0) + blend_r * XtoXpt_yy) * (origFlux - linFlux)
                  + 2.0 * blend_r * XtoXpt_y * (gradPsi[1] - linFlux_y) + blend * secondDer[2];
                                            
      result[5] = (blend_rr * XtoXpt_y * XtoXpt_x + blend_r * XtoXpt_xy) * (origFlux - linFlux)
                  + blend_r * XtoXpt_x * (gradPsi[1] - linFlux_y) + blend_r * XtoXpt_y * ( gradPsi[0] - linFlux_x)
                  + blend * secondDer[1];

      result[0] = -result[3] / R;
      result[1] =  result[2] / R;

   }
   
   return result;
}


double BlockMapping::getLinearMagneticFlux( const RealVect& a_physical_coordinate ) const
{
   /*
    Return the linearized flux near the X-point.
    psi = abs(x_new) - abs(y_new)
    */
   
   double x_new = m_lin_1[0] * (a_physical_coordinate[0] - m_Xpoint[0])
                + m_lin_1[1] * (a_physical_coordinate[1] - m_Xpoint[1]);
   
   double y_new = m_lin_2[0] * (a_physical_coordinate[0] - m_Xpoint[0])
                + m_lin_2[1] * (a_physical_coordinate[1] - m_Xpoint[1]);

   return fabs(x_new) - fabs(y_new);

}

double
BlockMapping::getPoloidalArcLength(const int       a_pol_index,
                                   const string&   a_option) const
{
   /*
    Poloidal mesh packing is done here (e.g., uniform: length = dl * a_pol_index)
    */
    
   double length;
   double dl = m_lengthCORE/double(m_nthetaCORE);

   //Define xi (mapped poloidal coordinate)
   //lenght(xi=0)= 0 and lenght(xi=1)= m_lengthCORE/2
   double xi = a_pol_index /double(m_nthetaCORE/2);
   
   /*
    Uniform: uniform distribution of poloidal points
   */
   if (a_option == "uniform") {
      length = xi * m_lengthCORE/2.0;
   }

   /*
    Cosine: packs points closer (far apart) for positive (negative)
    values of amplitude near X-point
   */
   else if (a_option == "cosine") {
      double amplitude = 0.3;
      
      length = xi * m_lengthCORE/2.0 * (1.0 - amplitude * cos(0.5*PI*xi));
      
      //Use constant-size cells (corresponding to the mimimal cell-size at Xpt) in the divertor legs
      //because often-times the legs are short and we want to pack as many points there (without extrapolation) as we can
      if (m_block_number > RCSOL) {
         length = xi * m_lengthCORE/2.0 * (1.0-amplitude);
      }
   }
   
   /*
    Tanh: packs points closer (or far apart) near the X-point
    relative to the rest of the core.
    d(arcLenght) = d(lenght)/d(pol_index) = C1 * tanh(xi/C2) + C0
    lengh = C0*x + C1*C2*Ln(Cosh(x/C2))
    dxi_ratio = 1 + C1/C0 = dxi(xi=inf)/dxi(xi=0)
   */
   
   else if (a_option == "tanh") {
      double dxi_ratio = 0.25;
      double C2 = 0.05;
      double fac = 1.0 +(dxi_ratio-1.0)*C2*log(cosh(1.0/C2));
      double C0 = m_lengthCORE/fac/2.0;
      
      length = C0*(xi+(dxi_ratio-1.0)*C2*log(cosh(xi/C2)));
   }
   
   else {
      MayDay::Error("getPoloidalArcLength: unknown option");
   }
   
   //Correct for the derection of the poloidal-length increase (always increases away from Xptn)
   if (m_block_number == LCORE || m_block_number == LCSOL) length = m_lengthCORE - length;
   
   return length;
}

void
BlockMapping::getRadialGridLine(FArrayBox& a_radialGridLine,
                                const Vector<RealVect*> & a_mapping,
                                const double a_length) const
{
   //Is the grid line a block cut?
   bool is_cut(false);
   if (fabs(a_length) < m_tol || fabs(a_length - m_lengthCORE) < m_tol ) is_cut=true;
   
   //Get a starting point on a separatrix
   RealVect ref_point = is_cut ? m_Xpoint : interpolateFromArcLengthMapping(a_mapping, a_length);

   FArrayBox psi_coords(a_radialGridLine.box(),1);
   getNormFluxCoordinates(psi_coords);
   
   for (BoxIterator bit(a_radialGridLine.box());bit.ok();++bit) {
      IntVect iv = bit();

      if (iv[RADIAL_DIR]==0) {
         //get a point on the separatrix
         a_radialGridLine(iv,0) = ref_point[0];
         a_radialGridLine(iv,1) = ref_point[1];
      }
      
      else {
         
         IntVect old_index(iv);
         old_index[RADIAL_DIR] = iv[RADIAL_DIR] - 1;
         RealVect point;
         point[0] = a_radialGridLine(old_index,0);
         point[1] = a_radialGridLine(old_index,1);

         if (is_cut) {
            //get points on a block cut
            point[0] = m_block_cut_coords(iv,0);
            point[1] = m_block_cut_coords(iv,1);
         }

         else {
            //get interior points
            RealVect cut_dir;
            IntVect iv_cut(iv[0],0);
            cut_dir[0] = m_block_cut_dir(iv_cut,0);
            cut_dir[1] = m_block_cut_dir(iv_cut,1);
            
            double internal_arc_length(a_length);
            if (m_block_number == LCORE || m_block_number == LCSOL) internal_arc_length = m_lengthCORE - a_length;
            
            //use ray tracing to bring the point closer to the next flux surface
            //if immidiate linear pushToFluxSurface will not be good enough
            radialRayTracing(point, psi_coords(iv,0), internal_arc_length, cut_dir);
            
            RealVect increment_dir = getRadialDir(point, internal_arc_length, cut_dir);
            pushToFluxSurface(point,psi_coords(iv,0),increment_dir);
         }

         a_radialGridLine(iv,0) = point[0];
         a_radialGridLine(iv,1) = point[1];
      }
      
   }
   
}

void
BlockMapping::initializeBlockCut(FArrayBox&  a_block_cut_coords,
                                 FArrayBox&  a_block_cut_dir) const
{
   /*
    Computes coordinates and directional vector
    along a block cut
    */
   
   FArrayBox psi_coords(a_block_cut_coords.box(),1);
   getNormFluxCoordinates(psi_coords);
   
   for (BoxIterator bit(a_block_cut_coords.box());bit.ok();++bit) {
      IntVect iv = bit();

      if (iv[RADIAL_DIR]==0) {
         a_block_cut_coords(iv,0) = m_Xpoint[0];
         a_block_cut_coords(iv,1) = m_Xpoint[1];
      }
      
      else {
         //compute block cut coordinates
         IntVect old_index(iv);
         old_index[RADIAL_DIR] = iv[RADIAL_DIR] - 1;
         RealVect point;
         point[0] = a_block_cut_coords(old_index,0);
         point[1] = a_block_cut_coords(old_index,1);

         if (m_straight_block_bndry || iv[RADIAL_DIR]==1) {
            RealVect cut_coord = getStraightCutCoords(psi_coords(iv,0));
            point[0] = cut_coord[0];
            point[1] = cut_coord[1];
         }
         
         else {
            pushToFluxSurface(point,psi_coords(iv,0));
         }

         a_block_cut_coords(iv,0) = point[0];
         a_block_cut_coords(iv,1) = point[1];
         
         //compute block cut directions
         RealVect this_cut_dir;
         this_cut_dir[0] = point[0] - a_block_cut_coords(old_index,0);
         this_cut_dir[1] = point[1] - a_block_cut_coords(old_index,1);
         this_cut_dir /= vectMag(this_cut_dir);
         this_cut_dir *=signum(psi_coords(iv,0) - psi_coords(old_index,0));
      
         a_block_cut_dir(iv,0) = this_cut_dir[0];
         a_block_cut_dir(iv,1) = this_cut_dir[1];
         
      }
   }
}

RealVect
BlockMapping::getStraightCutCoords( const double a_psi ) const
{
   /*
    Return a physical coordinate corresponding to the intersection
    of the curve psi(R,Z)=a_psi and the straight block cut
   */
   

   bool residual_tolerance_satisfied = false;
   double epsilon = m_tol;
   double f_lo, f_hi, f_half;
   RealVect half_x, tmp_lo_x, tmp_hi_x;
   
   
   tmp_lo_x = m_Xpoint;

   if (m_block_number == LCORE || m_block_number == RCORE ) {
      tmp_hi_x[1] = m_magAxis[1];
      tmp_hi_x[0] = m_Xpoint[0] + (m_magAxis[1] - m_Xpoint[1]) * m_core_cut_dir[0]/m_core_cut_dir[1];
   }
   
   else if (m_block_number == LPF || m_block_number == RPF ) {
      tmp_hi_x[1] = m_Zmin;
      tmp_hi_x[0] = m_Xpoint[0] + (m_Zmin - m_Xpoint[1]) * m_core_cut_dir[0]/m_core_cut_dir[1];
   }

   else if (m_block_number == LCSOL || m_block_number == LSOL ) {
      tmp_hi_x[0] = m_Rmin;
      tmp_hi_x[1] = m_Xpoint[1] + (m_Rmin - m_Xpoint[0]) * m_sol_cut_dir[1]/m_sol_cut_dir[0];
   }

   else if (m_block_number == RCSOL || m_block_number == RSOL ) {
      tmp_hi_x[0] = m_Rmax;
      tmp_hi_x[1] = m_Xpoint[1] + (m_Rmax - m_Xpoint[0]) * m_sol_cut_dir[1]/m_sol_cut_dir[0];
   }
   
   int iter = 0;
   while ( iter++ < m_max_iterations && !residual_tolerance_satisfied ) {
      
      f_lo = getNormMagneticFlux(tmp_lo_x) - a_psi;
      f_hi = getNormMagneticFlux(tmp_hi_x) - a_psi;
      
      half_x = tmp_hi_x + tmp_lo_x;
      half_x *= 0.5;

      f_half = getNormMagneticFlux(half_x) - a_psi;
      
      if (f_lo * f_half < 0) {
         tmp_hi_x = half_x;
      }
      
      else {
         tmp_lo_x = half_x;
      }
      
      RealVect dx(tmp_hi_x);
      dx -= tmp_lo_x;
      double dx_mag = sqrt(pow(dx[0],2) + pow(dx[1],2));
      
      if ( dx_mag < epsilon) {
         residual_tolerance_satisfied = true;
      }
      
      if (iter>m_max_iterations-2) {
        cout<<"getCutCoord::max_number_of_iteration_reached"<<endl;
      }
      
   }
   
   return half_x;
   
}


Vector<Real>
BlockMapping::getSeparatrixDirAnalytic() const
{
   
   /*
    Return the separatrix unit vector directions in the CCW direction
    in RCORE(k_out) -- result[0,1]
    in LCORE(k_in)  -- result[2,3]
   */
   
   
   Vector<Real> mag_data = getMagFieldData(m_Xpoint);

   //Near the Xptn the separatrix is given by z=kx; (there are 2 solutions for k)
   //Assume that psi-psi_at_sep = c*x^2 + b*x*y + a*y^2
   double a = mag_data[6]/2.0;
   double b = mag_data[5];
   double c = mag_data[4]/2.0;
   
   double k1 = (-b + sqrt(b*b - 4.0 * a *c))/(2.0 * a);
   double k2 = (-b - sqrt(b*b - 4.0 * a *c))/(2.0 * a);
   
   double k_out, k_in;
   if (k1>0) {
     k_out = k1;
     k_in = k2;
   }

   else {
     k_out =k2;
     k_in = k1;
   }
   
   Vector<Real> result(4);
   
   result[0] = 1.0/sqrt(k_out*k_out+1.0);
   result[1] = k_out/sqrt(k_out*k_out+1.0);
   
   result[2] = 1.0/sqrt(k_in*k_in+1.0);
   result[3] = k_in/sqrt(k_in*k_in+1.0);

   return result;
}

RealVect
BlockMapping::getRadialDir(const RealVect& a_X, const double a_length, const RealVect& a_cut_dir) const
{
   
   /*
    
    Computes the unit vector in the grid "radial" direction at a_X point.
    This vector is parallel to blok cut at the block cut, is orthogonal
    to magnetic flux function at the separatrix, and is a belnd of these
    two directions anywhere in between.
    
    We define the "transition" direction which runs half-between the cut
    and the separatrix. The transition_angle parameter controls the 
    transition width of the radial direction from one assymptotic to the other.
    
    Trans_length parameters define the characteristic value of arc_length, such that 
    the radial direction is orthogonal to the magentic flux everywhere beyond this arc_length.
   */
   
   double transition_angle = PI/30.0;

   
   //Get the agnle transition direction (half of the angle between the block cut and separatrix)
   //loc_tran_dir is directed outward of Xpnt
   RealVect loc_sep_dir;
   RealVect loc_trans_dir;

   RealVect sep_in;
   sep_in[0] = m_sep_dir[2];
   sep_in[1] = m_sep_dir[3];
   
   RealVect sep_out;
   sep_out[0] = m_sep_dir[0];
   sep_out[1] = m_sep_dir[1];
   
   if (m_block_number == LCORE) {
      loc_sep_dir = sep_in * (-1.0);
      loc_trans_dir = loc_sep_dir;
      loc_trans_dir -= a_cut_dir;
   }

   else if (m_block_number == RCORE) {
      loc_sep_dir = sep_out;
      loc_trans_dir = loc_sep_dir;
      loc_trans_dir -= a_cut_dir;
   }

   else if (m_block_number == LCSOL) {
     loc_sep_dir = sep_in * (-1.0);
     loc_trans_dir = loc_sep_dir;
     loc_trans_dir += a_cut_dir;
   }

   else if (m_block_number == RCSOL) {
      loc_sep_dir = sep_out;
      loc_trans_dir = loc_sep_dir;
      loc_trans_dir += a_cut_dir;
   }
   
   else if (m_block_number == LSOL) {
      loc_sep_dir = sep_out * (-1.0);
      loc_trans_dir = loc_sep_dir;
      loc_trans_dir += a_cut_dir;
   }
   
   else if (m_block_number == RSOL) {
      loc_sep_dir = sep_in;
      loc_trans_dir = loc_sep_dir;
      loc_trans_dir += a_cut_dir;
   }
   
   else if (m_block_number == LPF) {
      loc_sep_dir = sep_out * (-1.0);
      loc_trans_dir = loc_sep_dir;
      loc_trans_dir -= a_cut_dir;
   }

   else { //RPF
      loc_sep_dir = sep_in;
      loc_trans_dir = loc_sep_dir;
      loc_trans_dir -= a_cut_dir;
   }
   
   loc_trans_dir/=vectMag(loc_trans_dir);
   
   //Get the angle distance between the local coord and the transition vector
   RealVect XtoXpt(a_X[0]-m_Xpoint[0],a_X[1]-m_Xpoint[1]);
   RealVect XtoXpt_dir = XtoXpt/vectMag(XtoXpt);
   double angle_dist = acos(vectProd(XtoXpt_dir,loc_trans_dir));

   double proj_to_cut = projection(XtoXpt, a_cut_dir);
   double dist_to_cut = sqrt(pow(vectMag(XtoXpt),2) - pow(proj_to_cut,2));

   //Get the "sign" which defines whether the a_X coord lies closer to the separatrix (positive) or to the cut (negative)
   double sign = vectProd(XtoXpt_dir-loc_trans_dir,loc_sep_dir-loc_trans_dir);
   sign /= fabs(sign);
   

   //Compute the blending coefficient
   double blend_coeff_angle = (1.0 - tanh(sign * angle_dist / transition_angle))/2.0;
   double blend_coeff_dist = exp(-pow(a_length/m_trans_length,4));
   double blend_coeff = blend_coeff_angle * blend_coeff_dist;
  
#if 0
   if (m_dealigned_grid) {
    trans_length = 0.6; 
    blend_coeff = exp(-pow(vectMag(XtoXpt)/m_trans_length,4));
   }
#else
   if (m_dealigned_grid) {
      double x = vectMag(XtoXpt)/m_trans_length;
      double exponent = 2;
      if ( x < 1. ) {
         blend_coeff = exp( - pow(x,exponent) / (1. - pow(x,exponent)) );
      }
      else {
         blend_coeff = 0.;
      }
   }
#endif

   RealVect gradPsiNorm_dir = getPsiNorm(a_X, true);
   RealVect result;
   for (int i=0; i<2; ++i) {
     result[i] =  blend_coeff * a_cut_dir[i] + (1.0-blend_coeff) * gradPsiNorm_dir[i];
   }

   result/=vectMag(result);
   return result;
}

RealVect
BlockMapping::getPsiNorm(const RealVect& a_X, const bool a_unit_vect) const
{
   /*
    Return a unit vector in the gradPsiNorm direction
   */
   
   Vector<Real> mag_field_data = getMagFieldData(a_X);
   RealVect gradPsiNorm;
   gradPsiNorm[0] = mag_field_data[2];
   gradPsiNorm[1] = mag_field_data[3];   
   gradPsiNorm /=(m_physFluxOnSep - m_physFluxOnAxis);

   if (a_unit_vect == true) {
     double gradPsiNormMag = sqrt(pow(gradPsiNorm[0],2) + pow(gradPsiNorm[1],2));

     if (gradPsiNormMag > 1.0e-6) {
       gradPsiNorm /= gradPsiNormMag;
     }
     else {
       gradPsiNorm[0] = 0.0;
       gradPsiNorm[1] = 0.0;
     }
   }

   return gradPsiNorm;
}

RealVect
BlockMapping::getPsiTang(const RealVect& a_X, const int a_block_number) const
{
   /*
    Return a unit vector in the tangential direction (i.e., aligned with psi)
    Top blocks: CCW direction, bottom blocks: CW direction
   */
   
   Vector<Real> mag_field_data = getMagFieldData(a_X);
   RealVect Bfield;
   Bfield[0] = mag_field_data[0];
   Bfield[1] = mag_field_data[1];
   double Bmag = sqrt(pow(Bfield[0],2) + pow(Bfield[1],2));
   
   RealVect dist_to_X(a_X);
   dist_to_X -= m_Xpoint;
   double dist_to_X_mag = sqrt(pow(dist_to_X[0],2) + pow(dist_to_X[1],2));
   
   RealVect tang_vect(Bfield);
   
   //Top blocks
   if (a_block_number < LSOL) {
      
      // Within m_dl distance (which is the finest scale in the arc lenght mapping)
      // we approximate tang_dir with m_sep_dir computed at Xpt
      if (dist_to_X_mag < m_dl) {
         tang_vect[0] = (a_X[0] - m_Xpoint[0]) > 0 ? m_sep_dir[0] : m_sep_dir[2];
         tang_vect[1] = (a_X[0] - m_Xpoint[0]) > 0 ? m_sep_dir[1] : m_sep_dir[3];
      }
      
      else {
         RealVect r_loc(a_X);
         r_loc -= m_magAxis;
         double vect_prod =r_loc[0]*Bfield[1] - r_loc[1]*Bfield[0];
         double sign = (vect_prod > 0) ? 1.0 : -1.0;
         tang_vect /= Bmag;
         tang_vect *= sign;
      }
   }
   
   //Bottom blocks
   if (a_block_number >= LSOL) {
      
      // Within a couple of m_dl we approximate tang_dir with m_sep_dir computed at Xpt
      if (dist_to_X_mag < m_dl) {
         tang_vect[0] = (a_X[0] - m_Xpoint[0]) > 0 ? m_sep_dir[2] : m_sep_dir[0];
         tang_vect[1] = (a_X[0] - m_Xpoint[0]) > 0 ? m_sep_dir[3] : m_sep_dir[1];
      }
      
      else {
         RealVect axis(m_core_cut_dir);
         axis *= 0.1;
         axis += m_Xpoint;
         RealVect r_loc(a_X);
         r_loc -= axis;
         double vect_prod =r_loc[0]*Bfield[1] - r_loc[1]*Bfield[0];
         double sign = (vect_prod > 0) ? -1.0 : 1.0;
         tang_vect /= Bmag;
         tang_vect *= sign;
      }
   }
   
   return tang_vect;
}

void
BlockMapping::radialRayTracing( RealVect&       a_X,
                                const double    a_psi_norm,
                                const double    a_arc_length,
                                const RealVect& a_cut_dir) const
{
   
   /*
    If needed, push a_X in radial direction using ray-trace,
    until it is close enough to a_psi_norm flux surface, such that
    linear interpolation in pushToFluxSurface can be used.
    This may be required in the X-point region.
    */

   double tol = 1.0e-5;
   RealVect start_pt(a_X);
   RealVect end_pt;
   bool point_is_linearly_close = false;
   int n_bisection = 0;
   
   while (!point_is_linearly_close) {
      
      double dpsi = a_psi_norm - getNormMagneticFlux(start_pt);

      if (n_bisection > 0) {
         dpsi = dpsi/pow(2,n_bisection);
      }
      
      RealVect gradPsiNorm = getPsiNorm(start_pt,false);
      double gradPsiNorm_mag = sqrt(pow(gradPsiNorm[0],2)+pow(gradPsiNorm[1],2));
      RealVect gradPsiNorm_dir(gradPsiNorm);
      gradPsiNorm_dir /= gradPsiNorm_mag;

      RealVect start_dir = getRadialDir(start_pt, a_arc_length, a_cut_dir);
      double cosine = vectProd(start_dir, gradPsiNorm_dir);
      
      end_pt[0] = start_pt[0] + (dpsi/cosine/gradPsiNorm_mag) * start_dir[0];
      end_pt[1] = start_pt[1] + (dpsi/cosine/gradPsiNorm_mag) * start_dir[1];
      
      RealVect end_dir = getRadialDir(end_pt, a_arc_length, a_cut_dir);
      
      // close enought, so the linear interpolation in the following up
      // call to pushToFluxSurface should have sufficient accuracy
      if ( vectProd(start_dir,end_dir) > (1.0-tol) && n_bisection == 0 ) {
         point_is_linearly_close = true;
      }

      // keep reducing the step in order to do linear ray-tracing
      else if (vectProd(start_dir,end_dir) < (1.0-tol) ) {
         n_bisection = n_bisection + 1;
      }
      
      // push the point to the next position on the ray
      else {
         start_pt = end_pt;
         n_bisection = 0;
      }
   }
   
   a_X = start_pt;
}

void
BlockMapping::pushToFluxSurface(RealVect& a_X, const double a_psi_norm) const
{
   /*
    Push the vector a_X to the flux surface a_psi_norm along the normal direction
   */
   
   double dpsi = a_psi_norm - getNormMagneticFlux(a_X);

   RealVect gradPsiNorm = getPsiNorm(a_X,false);
   double gradPsiNorm_mag = sqrt(pow(gradPsiNorm[0],2)+pow(gradPsiNorm[1],2));
   RealVect gradPsiNorm_dir(gradPsiNorm);
   gradPsiNorm_dir /= gradPsiNorm_mag;
   
   bool residual_tolerance_satisfied(false);
   int iter=0;

   RealVect X_tmp(a_X);
   
   // start with using fast Newton-like (but with fixed gradient) method
   while ( iter++ < m_max_iterations && !residual_tolerance_satisfied ) {
      
      X_tmp[0] = X_tmp[0] + (dpsi/gradPsiNorm_mag) * gradPsiNorm_dir[0];
      X_tmp[1] = X_tmp[1] + (dpsi/gradPsiNorm_mag) * gradPsiNorm_dir[1];
      
      dpsi = a_psi_norm - getNormMagneticFlux(X_tmp);
      
      if ( fabs(dpsi)< m_tol ) {
         residual_tolerance_satisfied = true;
         a_X = X_tmp;
      }

      // if fixed gradient method failed, switch to
      // more robust, but slower bisection method
      if (iter > 20) {
         
         //reset parameters
         iter = 0;
         RealVect tmp_lo_x = a_X;
         RealVect tmp_hi_x = a_X;
         RealVect half_x = a_X;

         dpsi = a_psi_norm - getNormMagneticFlux(a_X);

         double f_lo = getNormMagneticFlux(tmp_lo_x) - a_psi_norm;
         double f_hi = f_lo;
         double f_half = f_lo;
         
         // find f_hi that bounds the solution
         while (f_hi * f_lo > 0) {
            tmp_hi_x[0] = tmp_hi_x[0] + (dpsi/gradPsiNorm_mag) * gradPsiNorm_dir[0];
            tmp_hi_x[1] = tmp_hi_x[1] + (dpsi/gradPsiNorm_mag) * gradPsiNorm_dir[1];
            f_hi = getNormMagneticFlux(tmp_hi_x) - a_psi_norm;
         }
         
         // find solution
         while ( iter++ < m_max_iterations && !residual_tolerance_satisfied ) {
            
            f_lo = getNormMagneticFlux(tmp_lo_x) - a_psi_norm;
            f_hi = getNormMagneticFlux(tmp_hi_x) - a_psi_norm;
            
            half_x = tmp_hi_x + tmp_lo_x;
            half_x *= 0.5;

            f_half = getNormMagneticFlux(half_x) - a_psi_norm;
            
            if (f_lo * f_half < 0) {
               tmp_hi_x = half_x;
            }
            
            else {
               tmp_lo_x = half_x;
            }
            
            RealVect dx(tmp_hi_x);
            dx -= tmp_lo_x;
            double dx_mag = sqrt(pow(dx[0],2) + pow(dx[1],2));
            
            if ( dx_mag < m_tol) {
               residual_tolerance_satisfied = true;
               a_X = half_x;
            }
            
            if (iter>m_max_iterations) {
              cout<<"pushToFluxSurface_norm::max number of iter reached"<<endl;
            }
         }
      }
   }
}

void
BlockMapping::pushToFluxSurface(RealVect& a_X, const double a_psi_norm, const RealVect& a_dir) const
{
   /*
    Push the vector a_X to the flux surface a_psi_norm along the specified direction
   */
                  
   double dpsi = a_psi_norm - getNormMagneticFlux(a_X);

   RealVect gradPsiNorm = getPsiNorm(a_X,false);
   double gradPsiNorm_mag = sqrt(pow(gradPsiNorm[0],2)+pow(gradPsiNorm[1],2));
   RealVect gradPsiNorm_dir(gradPsiNorm);
   gradPsiNorm_dir /= gradPsiNorm_mag;

   double cosine = vectProd(a_dir,gradPsiNorm_dir);
   
   bool residual_tolerance_satisfied(false);
   int iter=0;
   RealVect X_tmp(a_X);
   
   while ( iter++ < m_max_iterations && !residual_tolerance_satisfied ) {
                     
         X_tmp[0] = X_tmp[0] + (dpsi/cosine/gradPsiNorm_mag) * a_dir[0];
         X_tmp[1] = X_tmp[1] + (dpsi/cosine/gradPsiNorm_mag) * a_dir[1];
         
         dpsi = a_psi_norm - getNormMagneticFlux(X_tmp);
         
         if ( fabs(dpsi)< m_tol ) {
            residual_tolerance_satisfied = true;
            a_X = X_tmp;
         }

         // if fixed gradient method failed, switch to
         // more robust, but slower bisection method
         if (iter > 20) {
            
            //reset parameters
            iter = 0;
            RealVect tmp_lo_x = a_X;
            RealVect tmp_hi_x = a_X;
            RealVect half_x = a_X;

            dpsi = a_psi_norm - getNormMagneticFlux(a_X);

            double f_lo = getNormMagneticFlux(tmp_lo_x) - a_psi_norm;
            double f_hi = f_lo;
            double f_half = f_lo;
            
            // find f_hi that bounds the solution
            while (f_hi * f_lo > 0) {
               tmp_hi_x[0] = tmp_hi_x[0] + (dpsi/cosine/gradPsiNorm_mag) * a_dir[0];
               tmp_hi_x[1] = tmp_hi_x[1] + (dpsi/cosine/gradPsiNorm_mag) * a_dir[1];
               f_hi = getNormMagneticFlux(tmp_hi_x) - a_psi_norm;
            }
            
            // find solution
            while ( iter++ < m_max_iterations && !residual_tolerance_satisfied ) {
               
               f_lo = getNormMagneticFlux(tmp_lo_x) - a_psi_norm;
               f_hi = getNormMagneticFlux(tmp_hi_x) - a_psi_norm;
               
               half_x = tmp_hi_x + tmp_lo_x;
               half_x *= 0.5;

               f_half = getNormMagneticFlux(half_x) - a_psi_norm;
               
               if (f_lo * f_half < 0) {
                  tmp_hi_x = half_x;
               }
               
               else {
                  tmp_lo_x = half_x;
               }
               
               RealVect dx(tmp_hi_x);
               dx -= tmp_lo_x;
               double dx_mag = sqrt(pow(dx[0],2) + pow(dx[1],2));
               
               if ( dx_mag < m_tol) {
                  residual_tolerance_satisfied = true;
                  a_X = half_x;
               }
               
               if (iter>m_max_iterations) {
                 cout<<"pushToFluxSurface_norm::max number of iter reached"<<endl;
               }
            }
         }
    }
}

                  
void
BlockMapping::getCutDir(const Vector<Real>& a_sep_dir,
		     RealVect& a_core_cut_dir,
		     RealVect& a_sol_cut_dir ) const
{

   /*
    Compute unit vectors along the core and sol block cuts.
    Direction is along the increase of psi_norm as seen by top and outer blocks
    NB: a_sep_dir is in the CCW direction along the separatrix
   */
   
   double fac;

   a_core_cut_dir[0] = (-a_sep_dir[0] + a_sep_dir[2]);
   a_core_cut_dir[1] = (-a_sep_dir[1] + a_sep_dir[3]);
   fac = sqrt(pow(a_core_cut_dir[0],2) + pow(a_core_cut_dir[1],2));
   a_core_cut_dir[0] /= fac;
   a_core_cut_dir[1] /= fac;


   a_sol_cut_dir[0] = (a_sep_dir[0] + a_sep_dir[2])/sqrt(2.0);
   a_sol_cut_dir[1] = (a_sep_dir[1] + a_sep_dir[3])/sqrt(2.0);
   fac = sqrt(pow(a_sol_cut_dir[0],2) + pow(a_sol_cut_dir[1],2));
   a_sol_cut_dir[0] /= fac;
   a_sol_cut_dir[1] /= fac;
}


void
BlockMapping::getPsiLinearCoeff(RealVect& a_vec_1,
                                RealVect& a_vec_2 ) const
{
   
   /*
    
    Assuming psi = a*x^2+b*x*y+c*y^2, computes a_vec_1(a1,b1) and a_vec_2(a2,b2)
    such that x_new = a1 * x + b1 *y, y_new = a2 * x + b2 * y, and psi = x_new^2 - y_new^2
    
    */

   Vector<Real> mag_data = getMagFieldData(m_Xpoint);
   
   double c = mag_data[6]/2.0;
   double b = mag_data[5];
   double a = mag_data[4]/2.0;

   a_vec_1[0] = -(a + sqrt(b*b + pow((a - c),2)) - c) * sqrt(b*b + 2.0*c*(-a + sqrt(b*b + pow((a - c),2)) + c))/(2.0 *b*pow((b*b + pow((a - c),2)),(1.0/4.0)));
   a_vec_1[1] = -sqrt(b*b + 2.0*c*(-a + sqrt(b*b + pow((a - c),2)) + c))/(2.0 *pow((b*b + pow((a - c),2)),(1.0/4.0)));
   a_vec_2[0] = -sqrt(2*a*a + b*b - 2.0*a*(sqrt(b*b + pow((a - c),2)) + c))/(2.0 *pow((b*b + pow((a - c),2)),(1.0/4.0)));
   a_vec_2[1] = (a + sqrt(b*b + pow((a - c),2)) - c) * sqrt(2*a*a + b*b - 2.0*a*(sqrt(b*b + pow((a - c),2) ) + c))/(2.0 *b*pow((b*b + pow((a - c),2)),(1.0/4.0)));
}

void
BlockMapping::getArcLengthMapping(Vector<RealVect*> & a_mapping, 
                                  double&             a_block_length,
                                  const int           a_block_id) const
{

   CH_TIME("BlockMapping::getArcLengthMapping");
   
   /*
    Compute arc-length mapping along separatrix
    a_mapping is the array of points (RealVect's), which start at Xptn,
    and then follow the separatrix with an equdistant step of m_dl

    a_block_length is the (poloidal) arc length corresponding to the 
    length of the diverter legs if a_block_id =  LPF or RPF, 
    and the lenght of the entire tokamak for a_block_id = LCORE
   */
   
   bool end_is_reached = false;
   bool first_step = true;
   
   RealVect sep_in;
   sep_in[0] = m_sep_dir[2];
   sep_in[1] = m_sep_dir[3];

   RealVect sep_out;
   sep_out[0] = m_sep_dir[0];
   sep_out[1] = m_sep_dir[1];

   RealVect coord(m_Xpoint);
   RealVect tang_vect;
   double length_curr(0.);
   
   while ( !end_is_reached ) {
      
      RealVect *coord_ptr = new RealVect;
      *coord_ptr = coord;
      a_mapping.push_back(coord_ptr);
      
      if (first_step) {
         if (a_block_id == LCORE)    tang_vect = sep_out;
         else if (a_block_id == LPF) tang_vect = sep_out;
         else if (a_block_id == RPF) tang_vect = sep_in;
         else MayDay::Error("getArcLengthMapping: Block_id should be only LCORE, LPF, or RPF");
      }
      
      else tang_vect = getPsiTang(coord,a_block_id);
      
      //Compensate for the nomencalture direction in LPF
      if (a_block_id == LPF) tang_vect*= -1.0;
      
      tang_vect *= m_dl;
      coord += tang_vect;
      pushToFluxSurface(coord, 1.0);
      
      length_curr += m_dl;
      
      RealVect dist_to_X(coord);
      dist_to_X -= m_Xpoint;
      double dist_to_X_mag = sqrt(pow(dist_to_X[0],2) + pow(dist_to_X[1],2));
      
      //Check to see if we have returned to the Xpt for the core region, or
      //if we have reached the specified block_lenght for the divertor leg regions
      if ((a_block_id == LCORE && !first_step && dist_to_X_mag < m_dl && length_curr > 4.0*m_dl) ||
          (a_block_id != LCORE && (a_block_length < length_curr )) ) {
         
         end_is_reached = true;

         RealVect *coord_ptr = new RealVect;
         *coord_ptr = coord;
         a_mapping.push_back(coord_ptr);
         
         length_curr += dist_to_X_mag;
         
      }
      
      first_step = false;
      
      //Debugging
      //if (a_block_id == LCORE) cout<<"coord "<<coord[0]<<" "<<coord[1]<<" "<<dist_to_X_mag<<endl;
   }
   
   if (a_block_id == LCORE) a_block_length = length_curr;

}
RealVect
BlockMapping::interpolateFromArcLengthMapping(const Vector<RealVect*> & a_mapping,
                                              const double a_length) const
{
   /*
    Linearly interpolate the poloidal arc lenght mapping
    to return the point (RealVect) on the separatrix
    corresponding to the poloidal length a_length
    */
   
   int lo_ind = int(floor(a_length/m_dl));
   int hi_ind = lo_ind+1;

   double weight = double(a_length/m_dl - lo_ind);

   RealVect lo_vect = *a_mapping[lo_ind];
   RealVect hi_vect = *a_mapping[hi_ind];
  
   RealVect result(hi_vect);
   result -= lo_vect;
   result *= weight;
   result += lo_vect;

   pushToFluxSurface(result, 1.0);
   
   return result;

}


RealVect
BlockMapping::findCriticalPoint(const RealVect& a_init_guess,
                                const double    a_init_step ) const
{

  /*
   Finds a local minimum of (dpsi/dR)^2+(dpsi/dZ)^2
   Requires an initial guess (a_init_guess)
  */
   
   
  RealVect X(a_init_guess);
  double step_size(a_init_step);

  int iter = 0;
  bool residual_tolerance_satisfied = false;

  while ( iter++ < m_max_iterations && !residual_tolerance_satisfied ) {

    Vector<Real> mag_data = getMagFieldData(X);
    double residual = pow(mag_data[2],2) + pow(mag_data[3],2);

    RealVect grad;
    grad[0] = mag_data[2]*mag_data[4] + mag_data[3]*mag_data[5];
    grad[1] = mag_data[3]*mag_data[6] + mag_data[2]*mag_data[5];
    grad *= -2.0;

    RealVect incr(grad);
    incr /= vectMag(grad);
    incr *= step_size;
    RealVect tmp_X(X);
    tmp_X += incr;
    
    Vector<Real> tmp_mag_data = getMagFieldData(tmp_X);
    double tmp_residual = pow(tmp_mag_data[2],2) + pow(tmp_mag_data[3],2); 

    if (tmp_residual < residual) X = tmp_X;
    else step_size *= 0.5;

    if (step_size < m_tol) residual_tolerance_satisfied = true;
    if (iter == m_max_iterations - 2) cout<<"findCriticalPoint::max number of iter reached"<<endl;

  }

  return X;
}


void
BlockMapping::extrapolateGrid(FArrayBox& a_physical_coordinates, const Box& a_box) const
{
   /*
    Linearly extrapolate the grid into a_box index space
   */
   
   
   //Box iteration goes from left to right (min psi -> max psi), bottom to top (min theta -> max theta);
   for (BoxIterator bit(a_box);bit.ok();++bit) {

      IntVect iv = bit();
      IntVect index1(iv);
      index1[POLOIDAL_DIR] = iv[POLOIDAL_DIR] - 1;

      RealVect coord1;
      coord1[0] = a_physical_coordinates(index1,0);
      coord1[1] = a_physical_coordinates(index1,1);
      
      IntVect index2(iv);
      index2[POLOIDAL_DIR] = iv[POLOIDAL_DIR] - 2;
      RealVect coord2;
      coord2[0] = a_physical_coordinates(index2,0);
      coord2[1] = a_physical_coordinates(index2,1);
      
      RealVect increment(coord1);
      increment -= coord2;
      
      a_physical_coordinates(iv,0) = a_physical_coordinates(index1,0) + increment[0];
      a_physical_coordinates(iv,1) = a_physical_coordinates(index1,1) + increment[1];
   }
}

double
BlockMapping::vectProd(const RealVect& a_1, const RealVect& a_2) const
{
   
   /* 
    Return a scalar product 
   */
   
   return a_1[0]*a_2[0] + a_1[1]*a_2[1];
   
}

double
BlockMapping::vectMag(const RealVect& a_vect) const
{

   /*
    Return a vector length
    */

   return sqrt(a_vect[0]*a_vect[0] + a_vect[1]*a_vect[1]);
   
}

double
BlockMapping::projection(const RealVect& a_vec, const RealVect& a_targ) const
{
   /*
    Compute the absolute value of the projection of a_vec onto a_targ
   */
   
   double result;
   double a_targ_mag = sqrt(pow(a_targ[0],2) + pow(a_targ[1],2));
   result = fabs(a_vec[0]*a_targ[0] + a_vec[1]*a_targ[1]) / a_targ_mag;
   
   return result;
}

RealVect
BlockMapping::rotation(const RealVect& a_vec, const int a_dir) const
{
   /*
    Rotate vector on PI/2 in either left, CCW, (dir = 1) or right, CW, (dir = -1)
   */
   
   RealVect result;
   
   if (a_dir == 1) {
      result[0] = -a_vec[1];
      result[1] = a_vec[0];
   }
   
   else if (a_dir == -1) {
      result[0] = a_vec[1];
      result[1] = -a_vec[0];
   }
   
   else {
      MayDay::Error("Rotation: Unrecognized rotation dir");
   }
   
   return result;
}


void BlockMapping::getMinimalCellSize(Real&            a_size,
                                      RealVect&        a_location,
                                      const FArrayBox& a_coords,
                                      const int        a_dir) const
{
   const Box& box = a_coords.box();
   Box box_reduced(box);
   
   box_reduced.growHi(a_dir, -1);

   a_size = DBL_MAX;
   
   for (BoxIterator bit(box_reduced);bit.ok();++bit) {
      IntVect iv = bit();
      IntVect iv_hi(iv);
      iv_hi[a_dir] = iv[a_dir] + 1;
      
      RealVect x;
      RealVect x_hi;
     
      for (int dir=0; dir<2; ++dir) {
         x[dir] = a_coords(iv, dir);
         x_hi[dir] = a_coords(iv_hi, dir);
      }
         
      RealVect dx(x_hi);
      dx -= x;
      if (vectMag(dx) < a_size) {
         a_size = vectMag(dx);
         a_location = x;
      }
   }
}

/*
 ******************************************************************
 Some obsolete functions, which might be still useful in the future
 ******************************************************************
 */

#if 0

Vector<Real>
BlockMapping::getSeparatrixDirNumeric() const
{
   //Returns separatrix directions in RCORE(k_out) [0,1] and LCORE(k_in) [2,3] that follow separatrix CCW
   
   Vector<Real> result(4);
   
   double theta, fac;
   RealVect phys_coord;
   
   theta = PI/2.0/80.0;
   phys_coord = poloidaCoordRefPoint(theta, 1.0, RCORE);
   result[0] = phys_coord[0]-m_Xpoint[0];
   result[1] = phys_coord[1]-m_Xpoint[1];
   
   fac = sqrt(pow(result[0],2)+pow(result[1],2));
   
   result[0] /= fac;
   result[1] /= fac;
   
   theta = -PI/2.0/80.0;
   phys_coord = poloidaCoordRefPoint(theta, 1.0, LCORE);
   result[2] = -phys_coord[0]+m_Xpoint[0];
   result[3] = -phys_coord[1]+m_Xpoint[1];
   
   fac = sqrt(pow(result[2],2)+pow(result[3],2));
   
   result[2] /= fac;
   result[3] /= fac;
   
   
   return result;
   
}

RealVect
BlockMapping::poloidaCoordRefPoint(const double a_theta,
                                const double a_psi,
                                const int    a_block_number ) const
{
   //input a_theta: poloidal angle that goes counter clockwise with zero being at the Xpoint (reflected for lower blocks)
   //input a_psi:  normalized magneitc flux
   //output: physical coordinate corresponding to the intersection of the curve a_psi=const and the line a_theta=const
   
   
   //Get polar axis shift w/r/t Z-axis
   RealVect vert_dir(0,-1);
   double cosine = vectProd(m_core_cut_dir,vert_dir);
   double theta_shift = acos(cosine);
   
   //Get polar axis coordinates and shifted_theta
   //NB::Opportunity to introduce poloidal points packing
   double r = 0.8;
   RealVect axis(m_core_cut_dir);
   double theta_sc;
   if (a_block_number < 4 ){
      axis *= -r;
      theta_sc = a_theta - theta_shift;
      //theta_sc = theta_sc - PI/48.8 * cos(theta_sc);
   }
   else {
      axis *= r;
      theta_sc = a_theta + theta_shift;
      //theta_sc =theta_sc - PI/100.0 * cos(theta_sc);
   }
   axis += m_Xpoint;
   
   //Get maximum value for the solution raduis (important such that we don't enter into top/bottom blocks)
   double rmax = r / cos(a_theta);
   //Think more carefully -- TEMPORARY FIX
   //rmax = (rmax > 2*(m_Rmax-m_magAxis[0])) ? 1.5 * (m_Rmax-m_magAxis[0]) : rmax;
   
   if (rmax > 2 * (m_Rmax - m_Xpoint[0]) || cos(a_theta)<0.0 ) {
      rmax = 1.5 * (m_Rmax - m_Xpoint[0]);
   }
   
   //Do the Newton solve
   bool residual_tolerance_satisfied = false;
   double epsilon = m_tol;
   double half_r, tmp_lo_r, tmp_hi_r, f_lo_r, f_hi_r, f_half_r, x_lo, x_hi, half_x, z_lo, z_hi, half_z;
   
   tmp_lo_r = 0.0;
   tmp_hi_r = rmax;
   
   int side = (a_block_number < 4) ?  1 : -1;
   
   int iter = 0;
   while ( iter++ < m_max_iterations && !residual_tolerance_satisfied ) {
      
      x_lo = tmp_lo_r * sin(theta_sc) + axis[0];
      z_lo = axis[1] - side * tmp_lo_r * cos(theta_sc);
      f_lo_r = getNormMagneticFlux(RealVect(x_lo,z_lo)) - a_psi;
      
      x_hi = tmp_hi_r * sin(theta_sc) + axis[0];
      z_hi = axis[1] - side * tmp_hi_r * cos(theta_sc);
      f_hi_r = getNormMagneticFlux(RealVect(x_hi,z_hi)) - a_psi;
      
      half_r = 0.5 * (tmp_lo_r + tmp_hi_r);
      half_x = half_r * sin(theta_sc) + axis[0];
      half_z = axis[1] - side * half_r * cos(theta_sc);
      f_half_r = getNormMagneticFlux(RealVect(half_x,half_z)) - a_psi;
      
      if (f_lo_r * f_half_r < 0) {
         tmp_hi_r = half_r;
      }
      
      else {
         tmp_lo_r = half_r;
      }
      
      double x_test_hi = (half_r+m_tol) * sin(theta_sc) + axis[0];
      double x_test_lo = (half_r-m_tol) * sin(theta_sc) + axis[0];
      double z_test_hi = axis[1] - side * (half_r + m_tol) * cos(theta_sc);
      double z_test_lo = axis[1] - side * (half_r - m_tol) * cos(theta_sc);
      
      if ( (getNormMagneticFlux(RealVect(x_test_hi,z_test_hi)) - a_psi) * (getNormMagneticFlux(RealVect(x_test_lo,z_test_lo)) - a_psi) < 0 ) {
         residual_tolerance_satisfied = true;
      }
      
      if (iter>m_max_iterations-2) {
         cout<<"poloidaCoordRefPoint::max_number_of_iteration_reached"<<endl;
      }
      
   }
   
   RealVect result;
   result[0] = half_r * sin(theta_sc) + axis[0];
   result[1] = axis[1] - side * half_r * cos(theta_sc);
   
   return result;
   
}
#endif

#include "NamespaceFooter.H"


