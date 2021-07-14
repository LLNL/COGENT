#include <fstream>
#include <array>
#include "SingleNullPoloidalUtil.H"
#include "SingleNullPoloidalUtilF_F.H"
#include "HermiteInterp.H"
#include "BSplineInterp.H"
#include "SplineInterp.H"
#include "CONSTANTS.H"

enum interpolation_mode {FUNCTION, RADIAL_DERIVATIVE, POLOIDAL_DERIVATIVE, SECOND_RADIAL_DERIVATIVE, RADIAL_POLOIDAL_DERIVATIVE, SECOND_POLOIDAL_DERIVATIVE};
enum field_vars {RBR_VAR, RBZ_VAR, RB_VAR, UNIT_BR_VAR, UNIT_BPHI_VAR, UNIT_BZ_VAR, NUM_FIELD_VARS};
enum rz_vars {R_VAR, Z_VAR};

#include "NamespaceHeader.H"


SingleNullPoloidalUtil::SingleNullPoloidalUtil( const ParmParse&  a_pp,
                                                bool              a_spectral_field,
                                                double            a_RB_toroidal )
   : m_RB_toroidal(a_RB_toroidal),
     m_RZ_interp(NULL),
     m_field_interp(NULL)
{
   if ( a_spectral_field ) {
      
      ifstream fieldCoefficientsFile;
      
      int NR, NZ;
      
#ifdef CH_MPI
      if (procID() == 0) {
#endif
         string field_coefficients_file;
         if ( a_pp.contains("field_coefficients_file")) {
            a_pp.get("field_coefficients_file", field_coefficients_file);
         }
         else {
            MayDay::Error("SNCoreBlockCoordSys::init(): No field coefficients file specified");
         }
         
         fieldCoefficientsFile.open( field_coefficients_file.c_str() );
         
         if (!fieldCoefficientsFile) {
            cout << "Unable to open field coefficients file" << endl;
            exit(1);
         }
         
         cout << "Reading field coefficients from " << field_coefficients_file << endl;
         
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
      
      MPI_Bcast(&m_Rmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&m_Rmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&m_Zmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&m_Zmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&NR, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&NZ, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&m_magAxis, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&m_Xpoint, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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
      
      if ( a_pp.contains("psi_expansion_order")) {
         a_pp.get("psi_expansion_order", m_psi_expansion_order);
         if (m_psi_expansion_order > NR || m_psi_expansion_order > NZ) {
            MayDay::Error("SingleNullPoloidalUtil::SingleNullPoloidalUtil(): Insufficient number of DCT coefficients for requested psi_expansion_order");
         }
      }
      else {
         m_psi_expansion_order = NR;
         if (NZ < NR) m_psi_expansion_order = NZ;
      }
   }
}


SingleNullPoloidalUtil::~SingleNullPoloidalUtil()
{
   if ( m_RZ_interp ) delete m_RZ_interp;
   if ( m_field_interp ) delete m_field_interp;
}


void
SingleNullPoloidalUtil::setRZInterp( const ParmParse&  a_pp,
                                     const string&     a_interp_method,
                                     const FArrayBox&  a_interp_node_coords,
                                     const FArrayBox&  a_RZ_data )
{
   if ( a_interp_method == "bspline" ) {
      m_RZ_interp = new BSplineInterp(a_pp, a_interp_node_coords, a_RZ_data);
   }
   else if ( a_interp_method == "hermite" ) {
      m_RZ_interp = new HermiteInterp(a_pp, a_interp_node_coords, a_RZ_data);
   }
   else if ( a_interp_method == "spline") {
      m_RZ_interp = new SplineInterp(a_pp, a_interp_node_coords, a_RZ_data);
   }
}
                               


void
SingleNullPoloidalUtil::setFieldInterp( const ParmParse&  a_pp,
                                        const string&     a_interp_method,
                                        const FArrayBox&  a_interp_node_coords,
                                        const FArrayBox&  a_field_data )
{
   if ( a_interp_method == "bspline" ) {
      m_field_interp = new BSplineInterp(a_pp, a_interp_node_coords, a_field_data);
   }
   else if ( a_interp_method == "hermite" ) {
      m_field_interp = new HermiteInterp(a_pp, a_interp_node_coords, a_field_data);
   }
   else if ( a_interp_method == "spline") {
      m_field_interp = new SplineInterp(a_pp, a_interp_node_coords, a_field_data);
   }
}
                               


void
SingleNullPoloidalUtil::computeFieldDataSpectral( const FArrayBox&  a_RZ,
                                                  const int         a_dir,
                                                  FArrayBox&        a_BField,
                                                  FArrayBox&        a_BFieldMag,
                                                  FArrayBox&        a_BFieldDir,
                                                  FArrayBox&        a_gradBFieldMag,
                                                  FArrayBox&        a_curlBFieldDir,
                                                  FArrayBox&        a_BFieldDirdotcurlBFieldDir,
                                                  const bool        a_derived_data_only ) const
{
   // Get box intersection
   Box box = a_BField.box();
   box    &= a_BFieldMag.box();
   box    &= a_BFieldDir.box();
   box    &= a_gradBFieldMag.box();
   box    &= a_curlBFieldDir.box();
   box    &= a_BFieldDirdotcurlBFieldDir.box();
   
   FArrayBox RB(box,2);
   FArrayBox dRBdR(box,2);
   FArrayBox dRBdZ(box,2);
   
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
   
   for (BoxIterator bit(box); bit.ok(); ++bit) {
      IntVect iv = bit();
      
      double sR = (a_RZ(iv,0) - m_Rmin) * Rscale;
      double sZ = (a_RZ(iv,1) - m_Zmin) * Zscale;
      
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
      
      value *= Rscale;
      
      RB(iv,1) = value;  // R*B_Z
      
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
      
      RB(iv,0) = -value;  // R*B_R
      
      // Compute d2psi_dR2
      
      derivR = 2;
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
                      CHF_REAL(value) );
      
      value *= (Rscale * Rscale);
      
      dRBdR(iv,1) = value;  // d(R*BZ)/dR
      
      // Compute d2psi_dZ2
      
      derivR = 0;
      derivZ = 2;
      
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
      
      value *= (Zscale * Zscale);
      
      dRBdZ(iv,0) = -value;  // d(R*B_R)/dZ
      
      // Compute d2psi_dRdZ
      
      derivR = 1;
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
      
      value *= (Rscale * Zscale);
      
      dRBdR(iv,0) = -value;  // d(R*B_R)/dR
      dRBdZ(iv,1) =  value;  // d(R*B_Z)/dZ
   }
   
   delete[] temp;
   
   FORT_GET_FIELD_DATA_POLOIDAL( CHF_BOX(box),
                                 CHF_CONST_FRA(a_RZ),
                                 CHF_CONST_FRA(RB),
                                 CHF_CONST_FRA(dRBdR),
                                 CHF_CONST_FRA(dRBdZ),
                                 CHF_CONST_REAL(m_RB_toroidal),
                                 CHF_FRA(a_BField),
                                 CHF_FRA1(a_BFieldMag,0),
                                 CHF_FRA(a_BFieldDir),
                                 CHF_FRA(a_gradBFieldMag),
                                 CHF_FRA(a_curlBFieldDir),
                                 CHF_FRA1(a_BFieldDirdotcurlBFieldDir,0) );
}


array<double,3>
SingleNullPoloidalUtil::computeBField( const RealVect& a_X ) const
{
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
   
   
   double sR = (a_X[0] - m_Rmin) * Rscale;
   double sZ = (a_X[1] - m_Zmin) * Zscale;
   
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
   
   RealVect RB;
   
   // Compute dpsi_dR
   
   int derivR = 1;
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
                   CHF_REAL(value) );
   
   value *= Rscale;
   
   RB[1] = value;  // R*B_Z
   
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
   
   array<double,3> result;
   
   result[0] = RB[0]/a_X[0];
   result[1] = m_RB_toroidal/ a_X[0];
   result[2] = RB[1]/a_X[0];
   
   delete[] temp;
   
   return result;
}


array<double,3>
SingleNullPoloidalUtil::computeBUnit( const RealVect& a_X ) const
{
   array<double,3> BField = computeBField(a_X);

   double Bmag = 0.;
   for (int n=0; n<3; ++n) {
      Bmag += BField[n]*BField[n];
   }
   Bmag = sqrt(Bmag);

   array<double,3> b;
   for (int n=0; n<3; ++n) {
      b[n] = BField[n] / Bmag;
   }

   return b;
}


void
SingleNullPoloidalUtil::getNodalFieldData( const FArrayBox&  a_RZ,
                                           FArrayBox&        a_A,
                                           FArrayBox&        a_b,
                                           FArrayBox&        a_Bmag) const
{
   const Box& box = a_RZ.box();
   FArrayBox psi(box,1);
   FArrayBox RB(box,2);
   FArrayBox dRBdR(box,2);
   FArrayBox dRBdZ(box,2);
      
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
      
   for (BoxIterator bit(box); bit.ok(); ++bit) {
      IntVect iv = bit();
         
      double sR = (a_RZ(iv,0) - m_Rmin) * Rscale;
      double sZ = (a_RZ(iv,1) - m_Zmin) * Zscale;
         
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
                       CHF_REAL(value) );
         
      psi(iv,0) = value;
         
      // Compute dpsi_dR
         
      derivR = 1;
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
                       CHF_REAL(value) );
         
      value *= Rscale;
         
      RB(iv,1) = value;  // R*B_Z
         
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
         
      RB(iv,0) = -value;  // R*B_R
   }
      
   delete[] temp;
      
   FORT_GET_NODAL_FIELD_DATA_POLOIDAL( CHF_BOX(box),
                                       CHF_CONST_FRA(a_RZ),
                                       CHF_CONST_FRA1(psi,0),
                                       CHF_CONST_FRA(RB),
                                       CHF_CONST_REAL(m_RB_toroidal),
                                       CHF_FRA(a_A),
                                       CHF_FRA(a_b),
                                       CHF_FRA(a_Bmag));
}


void
SingleNullPoloidalUtil::getMagneticFluxFromDCT( const FArrayBox&  a_physical_coordinates,
                                                FArrayBox&        a_magnetic_flux ) const
{
   /*
    Given the FArrayBox a_physical_coordinates of physical coordinates, fill the
    components of a_magnetic_flux with the magnetic flux by evaluating its
    DCT expansion.
    */
   
   const Box& box(a_magnetic_flux.box());
   CH_assert(a_physical_coordinates.box().contains(box));
   
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
   
   for (BoxIterator bit(box); bit.ok(); ++bit) {
      IntVect iv = bit();
      
      double sR = (a_physical_coordinates(iv,0) - m_Rmin) * Rscale;
      double sZ = (a_physical_coordinates(iv,1) - m_Zmin) * Zscale;
      
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
                      CHF_REAL(a_magnetic_flux(iv,0)) );
   }
   
   delete[] temp;
}


double
SingleNullPoloidalUtil::getMagneticFluxFromDCT( const RealVect& a_physical_coordinate ) const
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
   
   return psi;
}


RealVect SingleNullPoloidalUtil::mappedCoordNewton( const RealVect&  a_X,
                                                    const RealVect&  a_xi_initial,
                                                    const IntVect&   a_iv_initial ) const
{
   /*
     Use Newton iteration to evaluate the mapping of physical to computational
     coordinates by inverting the mapping of computational to physical coordinates.
   */
   double tol = 1.e-10;         // Read from input if we find some reason to later
   int max_iter = 40;

   RealVect xi = a_xi_initial;

   RealVect residual;

   double R = m_RZ_interp->interpolate(R_VAR, FUNCTION, xi);
   double Z = m_RZ_interp->interpolate(Z_VAR, FUNCTION, xi);
   double dRdr = m_RZ_interp->interpolate(R_VAR, RADIAL_DERIVATIVE, xi);
   double dRdtheta = m_RZ_interp->interpolate(R_VAR, POLOIDAL_DERIVATIVE, xi);
   double dZdr = m_RZ_interp->interpolate(Z_VAR, RADIAL_DERIVATIVE, xi);
   double dZdtheta = m_RZ_interp->interpolate(Z_VAR, POLOIDAL_DERIVATIVE, xi);

   residual[0] = R - a_X[0];
   residual[1] = Z - a_X[1];

   double Fnorm = residual.vectorLength();

   bool converged = Fnorm <= tol;
   int num_iters = 0;
   int bt_steps;

   RealVect xi_saved[100];
   double Fnorm_saved[100];

   xi_saved[0] = xi;
   Fnorm_saved[0] = Fnorm;

   while ( !converged && num_iters < max_iter) {

      // Invert the Jacobian to get the update
      double a = dRdr;
      double b = dRdtheta;
      double c = dZdr;
      double d = dZdtheta;
      double J = a*d - b*c;

      RealVect delta;
      delta[0] = - ( d * residual[0] - b * residual[1] ) / J;
      delta[1] = - (-c * residual[0] + a * residual[1] ) / J;

      double s = 1.;

      RealVect xi_proposed = xi + delta;

      R = m_RZ_interp->interpolate(R_VAR, FUNCTION, xi_proposed);
      Z = m_RZ_interp->interpolate(Z_VAR, FUNCTION, xi_proposed);
      dRdr = m_RZ_interp->interpolate(R_VAR, RADIAL_DERIVATIVE, xi_proposed);
      dRdtheta = m_RZ_interp->interpolate(R_VAR, POLOIDAL_DERIVATIVE, xi_proposed);
      dZdr = m_RZ_interp->interpolate(Z_VAR, RADIAL_DERIVATIVE, xi_proposed);
      dZdtheta = m_RZ_interp->interpolate(Z_VAR, POLOIDAL_DERIVATIVE, xi_proposed);

      residual[0] = R - a_X[0];
      residual[1] = Z - a_X[1];

      double Fnorm_proposed = residual.vectorLength();

#if 1
      /*
        Backtracking algorithm from R. P. Pawlowski, J. N. Shadid, J. P. Simonis and H. F. Walker,
        "Globalization Techniques for Newton-Krylov Methods and Applications to the Fully Coupled
        Solution of the Navier-Stokes Equations", SIAM Review, 48 (4), 2006, pp. 700-721.
      */

      double t = 1.e-4;
      double eta = 0.;
      bt_steps = 0;
      int bt_steps_max = 15;

      bool backtracking = Fnorm_proposed > (1. - t*(1. - eta)) * Fnorm && s > 0.;

      while ( backtracking ) {

         double theta = 0.1;
         delta *= theta;
         eta = 1. - theta*(1. - eta);

         xi_proposed = xi + delta;

         R = m_RZ_interp->interpolate(R_VAR, FUNCTION, xi_proposed);
         Z = m_RZ_interp->interpolate(Z_VAR, FUNCTION, xi_proposed);
         dRdr = m_RZ_interp->interpolate(R_VAR, RADIAL_DERIVATIVE, xi_proposed);
         dRdtheta = m_RZ_interp->interpolate(R_VAR, POLOIDAL_DERIVATIVE, xi_proposed);
         dZdr = m_RZ_interp->interpolate(Z_VAR, RADIAL_DERIVATIVE, xi_proposed);
         dZdtheta = m_RZ_interp->interpolate(Z_VAR, POLOIDAL_DERIVATIVE, xi_proposed);

         residual[0] = R - a_X[0];
         residual[1] = Z - a_X[1];

         Fnorm_proposed = residual.vectorLength();

         backtracking = (Fnorm_proposed > (1. - t*(1. - eta)) * Fnorm) && (bt_steps < bt_steps_max);

         bt_steps++;
      }
#endif

      // Update the current solution and residual norm
      xi = xi_proposed;
      Fnorm = Fnorm_proposed;

      // Test convergence
      converged = Fnorm <= tol;

      num_iters++;

      xi_saved[num_iters] = xi;
      Fnorm_saved[num_iters] = Fnorm;
   }

#ifdef REPORT_NEWTON_FAILURE
   if ( !converged && num_iters >= max_iter ) {
      cout << "On block (" << m_poloidal_index << "," << m_toroidal_index << "): Newton solve did not converge at " << a_X << ", a_xi_initial = " << a_xi_initial << ", iv_init = " << a_iv_initial << ", Fnorm = " << Fnorm << endl;
#if 0
      for (int i=0; i<num_iters; ++i) {
         cout << i << " xi = " << xi_saved[i] << ", Fnorm = " << Fnorm_saved[i] << endl;
      }
#endif
   }
#endif

   return xi;
}


void
SingleNullPoloidalUtil::computeFieldDataFromMappingFile( const FArrayBox&  a_Xi,
                                                         const int         a_dir,
                                                         FArrayBox&        a_BField,
                                                         FArrayBox&        a_BFieldMag,
                                                         FArrayBox&        a_BFieldDir,
                                                         FArrayBox&        a_gradBFieldMag,
                                                         FArrayBox&        a_curlBFieldDir,
                                                         FArrayBox&        a_BFieldDirdotcurlBFieldDir,
                                                         const bool        a_derived_data_only ) const
{
   CH_assert(m_field_interp != NULL);

   Box box = a_Xi.box();

   FArrayBox RZ_data(box,5);
   m_RZ_interp->interpolate(R_VAR, 0, FUNCTION,            a_Xi, box, RZ_data);
   m_RZ_interp->interpolate(R_VAR, 1, RADIAL_DERIVATIVE,   a_Xi, box, RZ_data);
   m_RZ_interp->interpolate(R_VAR, 2, POLOIDAL_DERIVATIVE, a_Xi, box, RZ_data);
   m_RZ_interp->interpolate(Z_VAR, 3, RADIAL_DERIVATIVE,   a_Xi, box, RZ_data);
   m_RZ_interp->interpolate(Z_VAR, 4, POLOIDAL_DERIVATIVE, a_Xi, box, RZ_data);

   if( !a_derived_data_only ) {

      FArrayBox RBpoloidal_pol(box, SpaceDim);
      m_field_interp->interpolate(RBR_VAR, 0, FUNCTION, a_Xi, box, RBpoloidal_pol);
      m_field_interp->interpolate(RBZ_VAR, 1, FUNCTION, a_Xi, box, RBpoloidal_pol);

      // Construct the field from the mapping data, RBpol and RBtor
      FORT_GET_FIELD_FROM_RBPOL_RBTOR(CHF_BOX(box),
                                      CHF_CONST_REAL(m_RB_toroidal),
                                      CHF_CONST_FRA(RBpoloidal_pol),
                                      CHF_CONST_FRA1(RZ_data,0),
                                      CHF_FRA(a_BField));
   }

   // Compute the field magnitude and unit vector
   FORT_GET_FIELD_MAGNITUDE_AND_UNITVECTOR(CHF_BOX(box),
                                           CHF_CONST_FRA(a_BField),
                                           CHF_FRA1(a_BFieldMag,0),
                                           CHF_FRA(a_BFieldDir));

   // Compute RB = the major radius times the field magnitude
   FArrayBox RB(box, 1);
   RB.copy(RZ_data,0,0,1);
   RB *= a_BFieldMag;
   
   // Differentiate RB and the field unit vector components
   FArrayBox derivative_data(box, 8);
   if( a_derived_data_only ) {

      FArrayBox diffed_data(box, 4);
      BoxIterator bit(box);
      for (bit.begin();bit.ok();++bit) {
         IntVect iv = bit();
         diffed_data(iv,0) = RB(iv);
         for (int m=1; m<4; ++m) {
            diffed_data(iv,m) = a_BFieldDir(iv,m-1);
         }
      }

      ParmParse pp;
      Interp* interp = new HermiteInterp(pp, a_Xi, diffed_data);

      for( int i=0; i<4; ++i ) {
         interp->interpolate(i, 2*i  , RADIAL_DERIVATIVE,   a_Xi, box, derivative_data);
         interp->interpolate(i, 2*i+1, POLOIDAL_DERIVATIVE, a_Xi, box, derivative_data);
      }

      delete interp;
   }
   else {
      m_field_interp->interpolate(RB_VAR,        0, RADIAL_DERIVATIVE,   a_Xi, box, derivative_data);
      m_field_interp->interpolate(RB_VAR,        1, POLOIDAL_DERIVATIVE, a_Xi, box, derivative_data);
      m_field_interp->interpolate(UNIT_BR_VAR,   2, RADIAL_DERIVATIVE,   a_Xi, box, derivative_data);
      m_field_interp->interpolate(UNIT_BR_VAR,   3, POLOIDAL_DERIVATIVE, a_Xi, box, derivative_data);
      m_field_interp->interpolate(UNIT_BPHI_VAR, 4, RADIAL_DERIVATIVE,   a_Xi, box, derivative_data);
      m_field_interp->interpolate(UNIT_BPHI_VAR, 5, POLOIDAL_DERIVATIVE, a_Xi, box, derivative_data);
      m_field_interp->interpolate(UNIT_BZ_VAR,   6, RADIAL_DERIVATIVE,   a_Xi, box, derivative_data);
      m_field_interp->interpolate(UNIT_BZ_VAR,   7, POLOIDAL_DERIVATIVE, a_Xi, box, derivative_data);
   }

   //   int axisymmetric = m_axisymmetric? 1: 0;
   int axisymmetric = 1;

   // Compute the differentiated field quantities: grad(Bmag), curl(b) and b dot curl(b)
   FORT_GET_FIELD_DERIVATIVE_DATA(CHF_BOX(box),
                                  CHF_CONST_INT(axisymmetric),
                                  CHF_CONST_FRA1(RB,0),                      // RB
                                  CHF_CONST_FRA1(derivative_data,0),         // dRBdr
                                  CHF_CONST_FRA1(derivative_data,1),         // dRBdt
                                  CHF_CONST_FRA(a_BFieldDir),                // bunit
                                  CHF_CONST_FRA1(derivative_data,2),         // dbunitRdr
                                  CHF_CONST_FRA1(derivative_data,3),         // dbunitRdt
                                  CHF_CONST_FRA1(derivative_data,4),         // dbunitphidr
                                  CHF_CONST_FRA1(derivative_data,5),         // dbunitphdt
                                  CHF_CONST_FRA1(derivative_data,6),         // dbunitZdr
                                  CHF_CONST_FRA1(derivative_data,7),         // dbunitZdt
                                  CHF_CONST_FRA1(RZ_data,0),                 // R
                                  CHF_CONST_FRA1(RZ_data,1),                 // Rr
                                  CHF_CONST_FRA1(RZ_data,2),                 // Rt
                                  CHF_CONST_FRA1(RZ_data,3),                 // Zr
                                  CHF_CONST_FRA1(RZ_data,4),                 // Zt
                                  CHF_FRA(a_gradBFieldMag),                  // gradB
                                  CHF_FRA(a_curlBFieldDir),                  // curlbunit
                                  CHF_FRA1(a_BFieldDirdotcurlBFieldDir,0));  // bdotcurlbunit
}

#include "NamespaceFooter.H"


