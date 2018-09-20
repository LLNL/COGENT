#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <fstream>
#include <iostream>
#include <sstream>
#include "mappedGridIO.H"
using std::cout;
using std::endl;
using std::ofstream;


#include "parstream.H"
#include "ParmParse.H"
#include "MillerBlockCoordSys.H"
#include "FourthOrderUtil.H"
#include "BoxIterator.H"
#include "computeNorm.H"
#include "FABView.H"

/// Global variables for handling output:
static const char* pgmname = "testMiller" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
static bool verbose = true ;
static bool writePlotFiles = true;

#ifdef CH_USE_DOUBLE
static Real precision = 1.0e-15;
#else
static Real precision = 1.0e-7;
#endif


enum probtypeEnum {constantF = 0,
                   linearF,
                   quadraticF,
                   cubicF,
                   quarticF,
                   quinticF,
                   miller,
                   num_probtype};

//int probtype = constantF;
// int probtype = linearF;
//int probtype = quadraticF;
//int probtype = cubicF;
//int probtype = quarticF;
//int probtype = quinticF;
int probtype = miller;

double inner_radial_bdry;
double outer_radial_bdry;
double R0;
int const_minor_radius;

///
// Parse the standard test options (-v -q -h) and
// app-specific options (-S <domain_size>) out of the command line
///
void
parseTestOptions( int argc ,char* argv[] )
{
  for( int i = 1 ; i < argc ; ++i )
    {
      if( argv[i][0] == '-' ) //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if( strncmp( argv[i] ,"-v" ,3 ) == 0 )
            {
              verbose = true ;
              //argv[i] = "" ;
            }
          else if( strncmp( argv[i] ,"-q" ,3 ) == 0 )
            {
              verbose = false ;
              //argv[i] = "" ;
            }
          else if( strncmp( argv[i] ,"-h" ,3 ) == 0 )
            {
              pout() << "usage: " << pgmname << " [-hqv]" << std::endl ;
              exit( 99 ) ;
            }
        }
    }
  return ;
}

void initF(LevelData<FluxBox>& a_F,
           const ProblemDomain& a_levelDomain,
           const RealVect& a_dxLevel,
           const MillerBlockCoordSys* coordsys)
{
  DataIterator dit = a_F.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisF = a_F[dit];
      for (int dir=0; dir<SpaceDim; dir++)
        {
          // set offest (location of face center relative
          // to dx*intVect
          RealVect offset = a_dxLevel;
          offset *= 0.5;
          offset[dir] = 0.0;

          FArrayBox& thisFdir = thisF[dir];
          // this is the really slow way to do this
          BoxIterator bit(thisFdir.box());
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect loc = iv*a_dxLevel + offset;

              // first shot at things --  set up a constant field
              if (probtype == constantF)
                {
                  D_TERM(
                         thisFdir(iv, 0) = 1;,
                         thisFdir(iv,1) = 1;,
                         thisFdir(iv,2) = 1;)
                }
              else if (probtype == linearF)
                {
                  D_TERM(
                         thisFdir(iv,0) = loc[0];,
                         thisFdir(iv,1) = loc[1];,
                         thisFdir(iv,2) = loc[2];)
                }
              else if (probtype == quadraticF)
                {
                  D_TERM(
                         thisFdir(iv,0) = loc[0]*loc[0];,
                         thisFdir(iv,1) = loc[1]*loc[1];,
                         thisFdir(iv,2) = loc[2]*loc[2];)
                }
              else if (probtype == cubicF)
                {
                  RealVect hiLoc = loc+offset;
                  RealVect loLoc = loc-offset;

                  RealVect hiLoc4 = hiLoc*hiLoc*hiLoc*hiLoc;
                  RealVect loLoc4 = loLoc*loLoc*loLoc*loLoc;

                  RealVect avgVal = 0.25*(hiLoc4 - loLoc4);
                  avgVal /= a_dxLevel;
                  avgVal[dir] = loc[dir]*loc[dir]*loc[dir];

                  D_TERM(
                         thisFdir(iv,0) = avgVal[0];,
                         thisFdir(iv,1) = avgVal[1];,
                         thisFdir(iv,2) = avgVal[2];)
                }
              else if (probtype == quarticF)
                {
                  RealVect hiLoc = loc+offset;
                  RealVect loLoc = loc-offset;

                  RealVect hiLoc5 = hiLoc*hiLoc*hiLoc*hiLoc*hiLoc;
                  RealVect loLoc5 = loLoc*loLoc*loLoc*loLoc*loLoc;

                  RealVect avgVal = 0.2*(hiLoc5 - loLoc5);
                  avgVal /= a_dxLevel;
                  // average on normal face is just the point value
                  avgVal[dir] = loc[dir]*loc[dir]*loc[dir]*loc[dir];

                  D_TERM(
                         thisFdir(iv,0) = avgVal[0];,
                         thisFdir(iv,1) = avgVal[1];,
                         thisFdir(iv,2) = avgVal[2];)

                }
              else if (probtype == miller)
                 {
                    RealVect real_loc= coordsys->realCoord(loc);

                    double r = loc[0] + inner_radial_bdry;
                    double theta = loc[1];
                    double R = real_loc[0];

#if 0
                    D_TERM(
                           thisFdir(iv,0) = 1./R;,
                           thisFdir(iv,1) = 0.;,
                           thisFdir(iv,2) = 0.;)
#else

                    if (dir==0) {
                       double pi = 4.*atan(1.);
                       double half_theta_lo = 0.5*(theta - 0.5*a_dxLevel[1]);
                       double half_theta_hi = 0.5*(theta + 0.5*a_dxLevel[1]);

                       for (int i=0; i<4; ++i) {
                          if (half_theta_hi> 0.5*pi) {
                             half_theta_lo -= pi;
                             half_theta_hi -= pi;
                          }
                       }

                       double sq = sqrt(R0*R0 - r*r);
                       double lo = atan(sq*tan(half_theta_lo)/(R0+r));
                       double hi = atan(sq*tan(half_theta_hi)/(R0+r));
                       double av = 2.*(hi-lo)/(a_dxLevel[1]*sq);

                       D_TERM(
                              thisFdir(iv,0) = av;,
                              thisFdir(iv,1) = 0.;,
                              thisFdir(iv,2) = 0.;)

                    }
                    else if (dir==1) {

                       double ct = cos(theta);
                       double r_lo = r - 0.5*a_dxLevel[0];
                       double r_hi = r + 0.5*a_dxLevel[0];
                          double av;

                          if (fabs(ct) > 1.e-13) {
                          av = log((R0 + r_hi*ct)/(R0 + r_lo*ct)) / (a_dxLevel[0]*ct);
                       }
                       else {
                          av = 1./R0;
                       }

                       D_TERM(
                              thisFdir(iv,0) = av;,
                              thisFdir(iv,1) = 0.;,
                              thisFdir(iv,2) = 0.;)
                    }
#endif
                       }
              else
                {
                  // bad probtype
                  pout() << "Bad probtype == " << probtype << endl;
                  MayDay::Error("Don't know what to do");
                }
            }
        } //end loop over face directions

    } // end loop over boxes

#if 0
  if (probtype == miller) {
     // Convert face-centered values to face averages
     fourthOrderAverage(a_F);
  }
#endif

}



void exactDivF(LevelData<FArrayBox>& a_divF,
               const ProblemDomain& a_levelDomain,
               const RealVect& a_dxLevel)
{
  // set offest (location of cell center relative
  // to dx*intVect
  RealVect offset = a_dxLevel;
  offset *= 0.5;

  //Real cellVol = D_TERM(a_dxLevel[0],*a_dxLevel[1],*a_dxLevel[2]);

  DataIterator dit = a_divF.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisDiv = a_divF[dit];

      // this is the really slow way to do this
      BoxIterator bit(thisDiv.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect loc = iv*a_dxLevel + offset;

          if (probtype == constantF)
            {
              thisDiv(iv, 0) = 0.0;
            }
          else if (probtype == linearF)
            {
              thisDiv(iv, 0) = SpaceDim;
            }
          else if (probtype == quadraticF)
            {
              thisDiv(iv, 0) = 2*(D_TERM(loc[0],+loc[1],+loc[2]));
            }
          else if (probtype == cubicF)
            {
              RealVect hiLoc = loc + 0.5*a_dxLevel;
              RealVect loLoc = loc - 0.5*a_dxLevel;
              //              thisDiv(iv, 0) = 3*(D_TERM(loc[0]*loc[0],+loc[1]*loc[1],+loc[2]*loc[2]));
              thisDiv(iv,0) = D_TERM((pow(hiLoc[0],3)-pow(loLoc[0],3))/a_dxLevel[0],
                                     +(pow(hiLoc[1],3)-pow(loLoc[1],3))/a_dxLevel[1],
                                     +(pow(hiLoc[2],3)-pow(loLoc[2],3))/a_dxLevel[2]);
              //thisDiv(iv,0) /= cellVol;
            }


          else if (probtype == quarticF)
            {
              //thisDiv(iv, 0) = 4*(D_TERM(loc[0]*loc[0]*loc[0],+loc[1]*loc[1]*loc[1],+loc[2]*loc[2]*loc[2]));
              RealVect hiLoc = loc + 0.5*a_dxLevel;
              RealVect loLoc = loc - 0.5*a_dxLevel;
              //              thisDiv(iv, 0) = 3*(D_TERM(loc[0]*loc[0],+loc[1]*loc[1],+loc[2]*loc[2]));
              thisDiv(iv,0) = D_TERM((pow(hiLoc[0],4)-pow(loLoc[0],4))/a_dxLevel[0],
                                     +(pow(hiLoc[1],4)-pow(loLoc[1],4))/a_dxLevel[1],
                                     +(pow(hiLoc[2],4)-pow(loLoc[2],4))/a_dxLevel[2]);
              //thisDiv(iv,0) /= cellVol;

            }
          else if (probtype == miller)
             {
                thisDiv(iv,0) = 0.;
             }
          else
            {
              // bad probtype
              pout() << "Bad probtype == " << probtype << endl;
              MayDay::Error("Don't know what to do");
            }
        }
    } // end loop over boxes
}


int
testMiller(ParmParse& a_pp);

int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;

  ParmParse pp( argc-2, argv+2, NULL, argv[1] );

  if( verbose )
    pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int status = testMiller(pp);

  if( status == 0 )
    cout << indent << "MillerBlockCoordSys test" << " passed." << endl ;
  else
    cout << indent << "MillerBlockCoordSys test" << " failed with return code " << status << endl ;



  if( status == 0 )
    cout << indent << pgmname << " passed." << endl ;
  else
    cout << indent << pgmname << " failed with return code " << status << endl ;


#ifdef CH_MPI
  MPI_Finalize ();
#endif
  return status ;
}

int
testMiller(ParmParse& a_pp)
{
  int returnCode = 0;

  ParmParse gksystemPP("gksystem");

  Vector<int> numCells;
  gksystemPP.getarr("num_cells", numCells, 0, SpaceDim);

  Vector<int> isPeriodic;
  gksystemPP.getarr("is_periodic", isPeriodic, 0, SpaceDim);

  ParmParse millerPP("gksystem.magnetic_geometry_mapping.miller");

  millerPP.get("inner_radial_bdry", inner_radial_bdry);
  millerPP.get("outer_radial_bdry", outer_radial_bdry);

  Vector<Real> origin;
  millerPP.getarr("origin", origin, 0, SpaceDim);
  R0 = origin[0];

  if (millerPP.contains("l_const_minorrad")) {
     millerPP.get("l_const_minorrad", const_minor_radius);
  }
  else {
     const_minor_radius = 0;   // default, the real Miller geometrys
   }

  if (probtype == constantF)
    {
      pout() << "constant F test problem" << endl;
    }
  else if (probtype == linearF)
    {
      pout() << "linear test problem" << endl;
    }
  else if (probtype == quadraticF)
    {
      pout() << "quadratic test problem" << endl;
    }
  else if (probtype == cubicF)
    {
      pout() << "cubic test problem" << endl;
    }
  else if (probtype == quarticF)
    {
      pout() << "quartic test problem" << endl;
    }
  else if (probtype == miller)
    {
      pout() << "miller test problem" << endl;
    }
  else
    {
      pout() << "unknown test problem type" << endl;
    }

  Real TwoPi = 8.*atan(1.);
  Real radial_length = outer_radial_bdry - inner_radial_bdry;
  Real poloidal_length = TwoPi;
  Real toroidal_length = TwoPi;
  RealVect dxVect(D_DECL(radial_length/numCells[0], poloidal_length/numCells[1], toroidal_length/numCells[2]));
  IntVect hiVect(numCells);
  hiVect -= IntVect::Unit;
  Box domainBox(IntVect::Zero, hiVect);
  ProblemDomain baseDomain(domainBox);
  for (int dir=0; dir<SpaceDim; dir++) {
     baseDomain.setPeriodic(dir,isPeriodic[dir]!=0);
  }

  int numLevels = 4;

  Vector<int> nRefVect(numLevels-1, 2);
  RealVect stretch = RealVect::Unit;

  // simplest possible test
  ProblemDomain levelDomain = baseDomain;
  Vector<Box> boxes(1, domainBox);
  Vector<int> procAssign(1,0);

  IntVect ghostVect = 2*IntVect::Unit;

  Real L1Crse, L2Crse, maxCrse;
  Real L1volCrse, L2volCrse,maxVolCrse;

  RealVect dxLevel = dxVect;

  for (int level = 0; level<numLevels; level++)
    {
      DisjointBoxLayout levelGrids(boxes, procAssign, levelDomain);

      MillerBlockCoordSys* levelCoordSys = new MillerBlockCoordSys(millerPP,
                                                                   levelDomain);
      Box domain_box = levelDomain.domainBox();
      cout << "num_cells = " << domain_box.size(0) << " " << domain_box.size(1) << endl;

      levelCoordSys->define(levelDomain, dxLevel);

      if (levelCoordSys == NULL)
        {
          returnCode += 1;
          pout() << "MillerBlockCoordSysFactory::getCoordSys FAILED" << endl;
        }
      else
        {
          FourthOrderCoordSys* FOCS = dynamic_cast<FourthOrderCoordSys*>(levelCoordSys);
          const LevelData<FArrayBox>& cellVol = FOCS->getCellVolumes();
          LevelData<FArrayBox> exactVolume(cellVol.getBoxes(), cellVol.nComp());

#if 0
          // check face metric terms
          {
            const LevelData<FluxBox>& metricTerms = FOCS->getFaceMetricTerms();

            LevelData<FluxBox> metricErrors(levelGrids, metricTerms.nComp(),
                                            metricTerms.ghostVect());
            Real maxErr = 0.0;
            DataIterator dit = metricTerms.dataIterator();
            for (dit.begin(); dit.ok(); ++dit)
              {
                const FluxBox& thisMetricTerms = metricTerms[dit];
                FluxBox& thisError = metricErrors[dit];
                thisError.setVal(0.0);

                for (int faceDir=0; faceDir<SpaceDim; faceDir++)
                  {
                    const FArrayBox& thisMetricDir = thisMetricTerms[faceDir];
                    FArrayBox& thisErrorDir = thisError[faceDir];
                    thisErrorDir.setVal(0.0);
                    for (int dir=0; dir<SpaceDim; dir++)
                      {
                        // exactVal is the area of the face in the dir
                        // component, zero otherwise
                        Real exactVal = D_TERM(dxLevel[0],
                                               *dxLevel[1],
                                               *dxLevel[2]);
                        exactVal /= dxLevel[dir];

                        int normcomp = FOCS->getMetricTermComponent(dir,dir);
                        thisErrorDir.setVal(exactVal, normcomp);
                      }
                    thisErrorDir -= thisMetricDir;

                    if (verbose)
                      pout() << "FaceDir = " << faceDir << endl;

                    for (int comp=0; comp<thisErrorDir.nComp(); comp++)
                      {
                        Real thisMin = thisErrorDir.min(comp);
                        Real thisMax = thisErrorDir.max(comp);
                        if (verbose)
                          {
                            pout () << "comp " << comp
                                    << " thisMax = " << thisMax
                                    << endl;

                            pout () << "comp " << comp
                                    << " thisMin = " << thisMin
                                    << endl;
                          }

                        if (abs(thisMin) > maxErr) maxErr = abs(thisMin);
                        if (abs(thisMax) > maxErr) maxErr = abs(thisMax);
                      }
                  } // end loop over face directions
              } // end loop over grids

            Real eps = precision;
            if (maxErr > eps)
              {
                // fail
                returnCode += 100;
                pout() << "Face metric terms test FAILS"
                       << " -- Max(error) for face metric terms = "
                       << maxErr << endl;
              }

          } // end test of face metric terms
#endif


          // test cell volumes
          // move cellVolErr outside of context so we can put it into
          // the plotfile
          LevelData<FArrayBox> cellVolErr(cellVol.getBoxes(), cellVol.nComp());
          {

            cellVol.copyTo(cellVolErr);

            DataIterator dit=cellVol.dataIterator();

            double rbar = 0.5 * (inner_radial_bdry + outer_radial_bdry);

            if (probtype == miller) {
               for (dit.begin(); dit.ok(); ++dit)
                  {
                     Real TwoPi = 8.*atan(1.);
                     FArrayBox& this_exact_volume = exactVolume[dit];
                     const FArrayBox& this_cell_vol = cellVol[dit];
                     BoxIterator bit(this_exact_volume.box());
                     for (bit.begin(); bit.ok(); ++bit)
                        {
                           IntVect iv = bit();
                           RealVect loc = iv*dxLevel + 0.5*dxLevel;
                           Real rmid = inner_radial_bdry + loc[0];
                           Real theta_mid = loc[1];

                           double dr = dxLevel[0];
                           double dtheta = dxLevel[1];
                           double r_lo = rmid - 0.5*dr;
                           double r_hi = rmid + 0.5*dr;

#if 1
                           if (const_minor_radius == 1) {
                              this_exact_volume(iv,0) = TwoPi*rbar*dr*(R0*dtheta + rbar*2.*cos(theta_mid)*sin(0.5*dtheta));
                           }
                           else {
                              this_exact_volume(iv,0) = TwoPi*(R0*rmid*dr*dtheta
                                                               + (2./3.)*dr*(r_lo*r_lo + r_lo*r_hi + r_hi*r_hi)*cos(theta_mid)*sin(0.5*dtheta));
                           }

                           if(iv[1]==0) cout << this_exact_volume(iv,0) << " " << this_cell_vol(iv,0) << "  diff = " << this_exact_volume(iv) - this_cell_vol(iv) << endl;
#else
                           this_exact_volume(iv,0) = rmid*dr*dtheta;
#endif
                           //                      if(iv[1]==0) cout << this_exact_volume(iv,0) << " " << this_cell_vol(iv,0) << "  diff = " << this_exact_volume(iv) - this_cell_vol(iv) << endl;
                        }
                     cellVolErr[dit]-=this_exact_volume;
                  }
            }
            else {
               Real exactVol = D_TERM(dxLevel[0],*dxLevel[1],*dxLevel[2]);
               for (dit.begin(); dit.ok(); ++dit)
                  {
                     cellVolErr[dit]-=exactVol;
                  }
            }

            DisjointBoxLayout* fineGridPtr = NULL;
            int nRefFine = -1;
            Interval divInterval(0,0);

            Real L1Err = computeNorm(cellVolErr, fineGridPtr,
                                     nRefFine, dxLevel[0],
                                     divInterval, 1);

            Real L2Err = computeNorm(cellVolErr, fineGridPtr,
                                     nRefFine, dxLevel[0],
                                     divInterval, 2);

            Real maxErr = computeNorm(cellVolErr, fineGridPtr,
                                      nRefFine, dxLevel[0],
                                      divInterval, 0);

            pout () << "level " << level << ":" << endl;
            pout()  << "   err(cellVol): L1 = " << L1Err
                    << "   L2 = " << L2Err
                    << "   Max = " << maxErr << endl;

            if (level > 0)
              {
                Real L1Volrate, L2Volrate, maxVolRate;
                L1Volrate = log(Abs(L1volCrse/L1Err))/log(2.0);
                L2Volrate = log(Abs(L2volCrse/L2Err))/log(2.0);
                maxVolRate = log(Abs(maxVolCrse/maxErr))/log(2.0);
                if (L1Err == 0) L1Volrate = 0.0;
                if (L2Err == 0) L2Volrate = 0.0;
                if (maxErr == 0) maxVolRate = 0.0;

                pout() << "         rate:  L1 = " << L1Volrate
                       << "   L2 = " << L2Volrate
                       << "   Max = " << maxVolRate << endl;

                // test convergence rate, and fail if it isn't good enough
                bool fail = false;
                Real baselineRate = 3.75 + SpaceDim;

                if ((L1Volrate < baselineRate) && (L1Err > precision) )
                  {
                    fail = true;
                  }
                if ((L2Volrate < baselineRate) && (L2Err > precision) )
                  {
                    fail = true;
                  }
                if ((maxVolRate < baselineRate) && (maxErr > precision) )
                  {
                    fail = true;
                  }

                if (fail)
                  {
                    pout() << "Cell volume convergence test FAILS!" << endl;
                    returnCode += 1000;
                  }
              }
            L1volCrse = L1Err;
            L2volCrse = L2Err;
            maxVolCrse = maxErr;

          }
#if 0
          {
            // check face-centered 1/J and N/J terms
            const LevelData<FluxBox>& Jinverse = FOCS->getJInverse();
            const LevelData<FluxBox>& NJinverse = FOCS->getNJinverse();
            const LevelData<FluxBox>& metricTerms = FOCS->getFaceMetricTerms();

            LevelData<FluxBox> JInverseErrors(levelGrids, Jinverse.nComp(),
                                              Jinverse.ghostVect());

            LevelData<FluxBox> NJInverseErrors(levelGrids, NJinverse.nComp(),
                                               NJinverse.ghostVect());

            Real maxErr = 0.0;
            Real maxErrNJinverse[SpaceDim];
            D_TERM(
                   maxErrNJinverse[0] = 0;,
                   maxErrNJinverse[1] = 0;,
                   maxErrNJinverse[2] = 0;)

            DataIterator dit = Jinverse.dataIterator();
            for (dit.begin(); dit.ok(); ++dit)
              {
                const FluxBox& thisJinverse = Jinverse[dit];
                FluxBox& thisError = JInverseErrors[dit];
                thisError.setVal(0.0);

                const FluxBox& thisNJinverse = NJinverse[dit];
                FluxBox& thisNJerr = NJInverseErrors[dit];

                const FluxBox& thisMetricTerms = metricTerms[dit];

                // exact value is inverse of product of the stretches
                Real exactJinverseVal = 1.0/(D_TERM(stretch[0],
                                                    *stretch[1],
                                                    *stretch[2]));

                for (int dir=0; dir<SpaceDim; dir++)
                  {
                    const FArrayBox& thisJinverseDir = thisJinverse[dir];

                    FArrayBox& thisErrorDir = thisError[dir];
                    thisErrorDir.setVal(exactJinverseVal);
                    thisErrorDir -= thisJinverseDir;

                    if (verbose)
                      pout() << "Testing 1/J, N/J: FaceDir = " << dir << endl;

                    // don't expect correct values at highest and lowest
                    // faces in FaceDir, so don't count those in the norms
                    Box normBox(thisErrorDir.box());
                    normBox.grow(dir, -1);

                    for (int comp=0; comp<thisErrorDir.nComp(); comp++)
                      {
                        Real thisMin = thisErrorDir.min(normBox, comp);
                        Real thisMax = thisErrorDir.max(normBox, comp);
                        if (verbose)
                          {
                            pout () << "comp " << comp
                                    << " thisMax(1/J) = " << thisMax
                                    << endl;

                            pout () << "comp " << comp
                                    << " thisMin(1/J) = " << thisMin
                                    << endl;
                          }

                        if (abs(thisMin) > maxErr) maxErr = abs(thisMin);
                        if (abs(thisMax) > maxErr) maxErr = abs(thisMax);
                      }

                    const FArrayBox& thisNJinverseDir = thisNJinverse[dir];
                    const FArrayBox& thisNdir = thisMetricTerms[dir];
                    FArrayBox& thisNJerrDir = thisNJerr[dir];

                    Box NJintersectBox(thisNdir.box());
                    NJintersectBox &= thisNJinverseDir.box();

                    // now compute errors for N/J -- since N is a constant
                    // over faces, and we've already checked N, easiest thing
                    // to do here is to simply multiply the exact value of 1/J
                    // with N
                    thisNJerrDir.copy(thisNdir, NJintersectBox);
                    thisNJerrDir *= exactJinverseVal;
                    thisNJerrDir -= thisNJinverseDir;

                    for (int comp=0; comp<thisNJinverseDir.nComp(); comp++)
                      {
                        Real thisMin = thisNJerrDir.min(comp);
                        Real thisMax = thisNJerrDir.max(comp);
                        if (verbose)
                          {
                            pout () << "comp " << comp
                                    << " thisMax(N/J) = " << thisMax
                                    << endl;

                            pout () << "comp " << comp
                                    << " thisMin(N/J) = " << thisMin
                                    << endl;
                          }

                        if (abs(thisMin) > maxErrNJinverse[comp])
                          maxErrNJinverse[comp] = abs(thisMin);
                        if (abs(thisMax) > maxErrNJinverse[comp])
                          maxErrNJinverse[comp] = abs(thisMax);
                      }



                  } // end loop over face directions
              } // end loop over grids

            Real eps = 1.0e-7;
            if (maxErr > eps)
              {
                // fail
                returnCode += 10;
                pout() << "Jinverse test FAILS"
                       << " -- Max(error) for Jinverse terms = "
                       << maxErr << endl;
              }

            for (int comp=0; comp<SpaceDim; comp++)
              {
                if (maxErrNJinverse[comp] > eps)
                  {
                    // fail
                    returnCode += 20;
                    pout() << "NJinverse test FAILS"
                           << " -- Max(error) for NJinverse terms = "
                           << maxErrNJinverse << endl;
                  }
              } // end loop over components

          }   // end test of Jinverse
#endif

          // now test divergence operator
          LevelData<FluxBox> F(levelGrids, SpaceDim, ghostVect);

          initF(F, levelDomain, dxLevel, levelCoordSys);

          LevelData<FArrayBox> divF(levelGrids, 1, IntVect::Zero);

          levelCoordSys->mappedGridDivergence(divF, F);

          // now evaluate error
          {
            LevelData<FArrayBox> divError(levelGrids, 1);
            exactDivF(divError, levelDomain, dxLevel);

            //            Real cellVolume = D_TERM(dxLevel[0],*dxLevel[1],*dxLevel[2]);

            DataIterator dit = divError.dataIterator();
            for (dit.begin(); dit.ok(); ++dit)
              {
                 divF[dit]/=cellVol[dit];
                 divError[dit].minus(divF[dit]);
              }

            DisjointBoxLayout* fineGridPtr = NULL;
            int nRefFine = -1;
            Interval divInterval(0,0);

            Real L1Err = computeNorm(divError, fineGridPtr,
                                     nRefFine, dxLevel[0],
                                     divInterval, 1);

            Real L2Err = computeNorm(divError, fineGridPtr,
                                     nRefFine, dxLevel[0],
                                     divInterval, 2);

            Real maxErr = computeNorm(divError, fineGridPtr,
                                      nRefFine, dxLevel[0],
                                      divInterval, 0);

#if 1
            cout << "   L1(err(div) = " << L1Err
                   << "   L2(err(div) = " << L2Err
                   << "   Max(err(div) = " << maxErr << endl;
#else
            pout() << "   L1(err(div) = " << L1Err
                   << "   L2(err(div) = " << L2Err
                   << "   Max(err(div) = " << maxErr << endl;
#endif

            if (level > 0)
              {
                Real L1rate, L2rate, maxRate;
                L1rate = log(Abs(L1Crse/L1Err))/log(2.0);
                L2rate = log(Abs(L2Crse/L2Err))/log(2.0);
                maxRate = log(Abs(maxCrse/maxErr))/log(2.0);
                if (L1Err == 0) L1rate = 0.0;
                if (L2Err == 0) L2rate = 0.0;
                if (maxErr == 0) maxRate = 0.;
                pout() << "         rate:  L1 = " << L1rate
                       << "   L2 = " << L2rate
                       << "   Max = " << maxRate << endl;

                // check convergence rates for pass/fail
                Real baselineRate = 3.8;
                bool fail = false;
                if ((L1rate < baselineRate) && (L1Err > precision))
                  {
                    fail = true;
                  }
                if ((L2rate < baselineRate) && (L2Err > precision))
                  {
                    fail = true;
                  }
                if ((maxRate < baselineRate) && (maxErr > precision))
                  {
                    fail = true;
                  }

                if (fail)
                  {
                    returnCode += 100000;
                    pout() << "divergence convergence rate test FAILS" << endl;
                  }

              }
            L1Crse = L1Err;
            L2Crse = L2Err;
            maxCrse = maxErr;

            // this is a good place to write a test plot file
            // with the divergence and the error in the divergence
            if (writePlotFiles)
              {
                 //                LevelData<FArrayBox> plotData(levelGrids, 2);
                LevelData<FArrayBox> plotData(levelGrids, 1);
                const LevelData<FArrayBox>& cellVol = FOCS->getCellVolumes();
#if 0
                divF.copyTo(divF.interval(),
                            plotData, divF.interval());
#endif
                cellVol.copyTo(cellVol.interval(),
                               plotData, cellVol.interval());

#if 0
                Interval errorInterval(1,1);
                exactVolume.copyTo(cellVolErr.interval(),
                                  plotData, errorInterval);

#if 0
                cellVolErr.copyTo(cellVolErr.interval(),
                                  plotData, errorInterval);
#endif

#if 0
                divError.copyTo(divError.interval(),
                                plotData, errorInterval);
#endif
#endif
                stringstream fileRootString;
                fileRootString << "MillerData" << level << ".";

                string fileRoot(fileRootString.str().c_str());
                // use unigrid version, since we've only got one level here
                WriteMappedUGHDF5(fileRoot, levelGrids,
                                  plotData, *FOCS);

              } // end if we're writing plot files

          }

          // finally, if desired, write mesh to a file
          bool writeMesh = true;
          if (writeMesh && level == 0)
            {
              char filestr[80];
              sprintf(filestr, "grid%d.out", level);
              ofstream os(filestr, ios::out);

              if (os.fail())
                {
                  pout() << "cannot open grid output file " << os << endl;
                  MayDay::Error();
                }


              const Box& domBox = levelDomain.domainBox();

              // x-direction
              int dir=0;
              {
                IntVect loVect = domBox.smallEnd();
                IntVect hiVect = domBox.bigEnd();
                // account for cell centering by adding the Unit vector
                hiVect += IntVect::Unit;

                int ilo = loVect[dir];
                int ihi = hiVect[dir];
                for (int i=ilo; i< ihi; i++)
                  {
                    loVect[dir] = i;
                    hiVect[dir] = i;

                    RealVect loLoc(loVect);
                    RealVect hiLoc(hiVect);
                    loLoc *= dxLevel;
                    hiLoc *= dxLevel;
                    loLoc /= stretch;
                    hiLoc /= stretch;

                    RealVect xlo = levelCoordSys->realCoord(loLoc);
                    RealVect xhi = levelCoordSys->realCoord(hiLoc);

                    os << D_TERM(xlo[0] << "  ",
                                 << xlo[1] << "  ",
                                 << "   ")  << endl;

                    os << D_TERM(xhi[0] << "  ",
                                 << xhi[1] << "  ",
                                 << "   ")  << endl;

                    os << endl;
                  } // end loop over cells in this direction
              } // end x-direction

              // y-direction
              dir=1;
              {
                IntVect loVect = domBox.smallEnd();
                IntVect hiVect = domBox.bigEnd();
                // account for cell centering by adding the Unit vector
                hiVect += IntVect::Unit;

                int ilo = loVect[dir];
                int ihi = hiVect[dir];
                for (int i=ilo; i< ihi; i++)
                  {
                    loVect[dir] = i;
                    hiVect[dir] = i;

                    RealVect loLoc(loVect);
                    RealVect hiLoc(hiVect);
                    loLoc *= dxLevel;
                    hiLoc *= dxLevel;
                    loLoc /= stretch;
                    hiLoc /= stretch;

                    RealVect xlo = levelCoordSys->realCoord(loLoc);
                    RealVect xhi = levelCoordSys->realCoord(hiLoc);

                    os << D_TERM(xlo[0] << "  ",
                                 << xlo[1] << "  ",
                                 << "   ")  << endl;

                    os << D_TERM(xhi[0] << "  ",
                                 << xhi[1] << "  ",
                                 << "   ")  << endl;

                    os << endl;
                  } // end loop over cells in this direction
              } // end y-direction



              os.close();
            } // end if we're writing out the mesh

          delete levelCoordSys;
          if (level < nRefVect.size())
            {
              levelDomain.refine(nRefVect[level]);
              dxLevel /= nRefVect[level];
              for (int i=0; i<boxes.size(); i++)
                {
                  boxes[i].refine(nRefVect[level]);
                }
            }
        } // end if we returnd a valid CoordSys
    } // end loop over levels

  return returnCode;

}

