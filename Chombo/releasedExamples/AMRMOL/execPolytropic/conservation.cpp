#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <fstream>
using std::ifstream;
using std::ios;

#include "FABView.H"
// this lets us use dumpIVS, other dump functions
#include "DebugDump.H"

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "CH_Timer.H"
#include "memusage.H"
#include "OldTimer.H"

#include "AMR.H"
#include "AMRLevel.H"
#include "AMRLevelConsFactory.H"
#include "AMRLevelCons.H"

#include "MOLPolytropicPhysics.H"

#include "RampIBC.H"
#include "ExplosionIBC.H"
#include "Explosion1dIBC.H"
#include "GaussianIBC.H"
#include "Gaussian1dIBC.H"
#include "Gaussian1dvIBC.H"
#include "GaussianPBC.H"
#include "GaussianSmoothBC.H"
#include "CosinePowerIBC.H"
#include "CosineEachPowerIBC.H"
#include "VortexIBC.H"
#include "ZZSVortexIBC.H"
#include "ShearIBC.H"
#include "TrigIBC.H"
#include "SineIBC.H"
#include "WaveIBC.H"
#include "SourceIBC.H"
#include "ShockTubeIBC.H"
#include "ChannelShockIBC.H"
#include "ChannelModianoIBC.H"
#include "generalFuncs.H"

// #define HALEM_PROC_SPEED
#ifdef HALEM_PROC_SPEED
#include <cstdio>
#include <sys/sysinfo.h>
#include <machine/hal_sysinfo.h>
#endif

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#ifdef CH_Linux
// Should be undefined by default
#define TRAP_FPE
#undef  TRAP_FPE
#endif

#ifdef TRAP_FPE
static void enableFpExceptions();
#endif

OldTimer Everything    ("gov Everything", 0);
OldTimer TimeReadInput ("gov Read Input",   Everything);
OldTimer TimeSetupAMR  ("gov Setup AMR",    Everything);
OldTimer TimeRun       ("gov Run",          Everything);
OldTimer TimeConclude  ("gov Conclude",     Everything);

// Possible pressure relationships for the initial condition
#define PRESSURE_ISENTROPIC 0
#define PRESSURE_CONSTANT   1

// amrRun is a function (as opposed to inline in main()) to get
// around MPI scoping problems
void amrRun();

// setupFixedGrids allows fixed grids to be read in and used in this AMR
// computation example
void setupFixedGrids(Vector<Vector<Box> >& a_amrGrids,
                     const ProblemDomain&  a_domain,
                     int                   a_maxLevel,
                     int                   a_maxGridSize,
                     int                   a_blockFactor,
                     int                   a_verbosity,
                     std::string           a_gridFile);

// setupVortices allows vortex parameters to be read in and used in this AMR
// computation example
void setupVortices(Vector<RealVect>& a_center,
                   Vector<Real>&     a_radius,
                   Vector<Real>&     a_strength,
                   int               a_verbosity,
                   std::string       a_vortexFile);

// One more function for MPI
void dumpmemoryatexit();

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  // Start MPI
  MPI_Init(&a_argc,&a_argv);
#ifdef CH_AIX
  H5dont_atexit();
#endif
  setChomboMPIErrorHandler();
#endif

  int rank, number_procs;

#ifdef CH_MPI
  MPI_Comm_rank(Chombo_MPI::comm, &rank);
  MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
  rank = 0;
  number_procs = 1;
#endif

  if (rank == 0)
  {
    pout() << " number_procs = " << number_procs << endl;
  }

  OldTimer::TimerInit(rank);

  Everything.start();

  // Check for an input file
  char* inFile = NULL;

  if (a_argc > 1)
  {
    inFile = a_argv[1];
  }
  else
  {
    pout() << "Usage:  amrRun...ex <inputfile>" << endl;
    pout() << "No input file specified" << endl;
    return -1;
  }

  // Parse the command line and the input file (if any)
  ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

#ifdef TRAP_FPE
  enableFpExceptions ();
#endif

  // Run amrRun, i.e., do the computation
  amrRun();

  Everything.stop();

#ifndef CH_NTIMER
  Real end_memory = get_memory_usage_from_OS();

  pout() << endl
         << "Everything completed --- "
         << "mem: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << end_memory
         << " MB, time: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << Everything.wc_time()
         << " sec (wall-clock)" << endl << endl;
#endif

#if !defined(CH_NTIMER) && defined(CH_MPI)
  Real avg_memory, min_memory, max_memory;
  gather_memory_from_procs(end_memory, avg_memory, min_memory, max_memory);
#endif

  // OldTimer::TimerSummary();

#ifdef CH_MPI
  // Exit MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}

void amrRun()
{
  // Start timing the reading of the input file
  TimeReadInput.start();

  /*
    AMR inputs
   */

  // Read inputs that are prefixed with "amr."
  ParmParse ppAmr("amr");

  // This determines the amount of diagnostic output generated
  int verbosity = 0;
  ppAmr.query("verbosity",verbosity);
  CH_assert(verbosity >= 0);

  // Maximum AMR level limit
  int maxLevel = 0;
  ppAmr.get("max_level",maxLevel);
  int numReadLevels = Max(maxLevel,1);

  // Refinement ratios between levels
  std::vector<int> refRatios;
  // Note: this requires a refRatio to be defined for the finest level
  // (even though it will never be used)
  ppAmr.getarr("ref_ratio",refRatios,0,numReadLevels+1);

  // Number of coarse time steps from one regridding to the next
  std::vector<int> regridIntervals;
  ppAmr.getarr("regrid_interval",regridIntervals,0,numReadLevels);

  // Minimum dimension of a grid
  int blockFactor = 1;
  ppAmr.get("block_factor",blockFactor);

  // Maximum dimension of a grid
  int maxGridSize = 32;
  ppAmr.get("max_grid_size",maxGridSize);

  Real fillRatio = 0.75;
  ppAmr.get("fill_ratio",fillRatio);

  // Set up checkpointing
  int checkpointInterval = 0;
  ppAmr.query("checkpoint_interval",checkpointInterval);

  // Set up plot file writing
  int plotInterval = 0;
  ppAmr.query("plot_interval",plotInterval);

  // Set up whether to use subcycling in time.
  int useSubcyclingInt = 1;
  ppAmr.query("use_subcycling", useSubcyclingInt);
  bool useSubcycling = (useSubcyclingInt == 1);

  // For MeshRefine. The hyperbolic codes use a grid buffer of 1
  int gridBufferSize = 1;
  ppAmr.query("grid_buffer_size", gridBufferSize);

  // Determine if a fixed or variable time step will be used
  Real fixedDt = -1;
  ppAmr.query("fixed_dt",fixedDt);

  // Limit the time step growth
  Real maxDtGrowth = 1.1;
  ppAmr.get("max_dt_growth",maxDtGrowth);

  // Let the time step grow by this factor above the "maximum" before
  // reducing it
  Real dtToleranceFactor = 1.1;
  ppAmr.get("dt_tolerance_factor",dtToleranceFactor);

  /*
    AMRLevelFactory inputs
   */

  // Read inputs that are prefixed with "amrlevelfactory."
  ParmParse ppAmrlevelfactory("amrlevelfactory");

  // How far to extend refinement from cells newly tagged for refinement
  int tagBufferSize = 3;
  ppAmrlevelfactory.get("tag_buffer_size",tagBufferSize);

  // Threshold that triggers refinement
  Real refineThresh = 0.3;
  ppAmrlevelfactory.get ("refine_thresh",refineThresh);

  // Whether refinement threshold is scaled with dx
  int refinementIsScaledInt = 0;
  ppAmrlevelfactory.query("refinement_is_scaled", refinementIsScaledInt);
  bool refinementIsScaled = (refinementIsScaledInt == 1);

  // Whether to tag on pressure instead of on density
  int tagPressureInt = 0;
  ppAmrlevelfactory.query("tag_pressure", tagPressureInt);
  bool tagPressure = (tagPressureInt == 1);

  // Whether to tag on vorticity instead of on density
  int tagVorticityInt = 0;
  ppAmrlevelfactory.query("tag_vorticity", tagVorticityInt);
  bool tagVorticity = (tagVorticityInt == 1);

  // Use fourth order slopes:  default true
  int inFourthOrderSlopes = 1;
  bool useFourthOrderSlopes;
  ppAmrlevelfactory.query("use_fourth_order_slopes",inFourthOrderSlopes);
  useFourthOrderSlopes = (inFourthOrderSlopes == 1);

  // Do slope limiting:  default true
  int inPrimLimiting = 1;
  bool usePrimLimiting;
  ppAmrlevelfactory.query("use_prim_limiting",inPrimLimiting);
  usePrimLimiting = (inPrimLimiting == 1);

  // This should actually be 1 even if slope limiting is off
  int highOrderLimiterInt = 1;
  ppAmrlevelfactory.query("high_order_limiter", highOrderLimiterInt);
  bool highOrderLimiter = (highOrderLimiterInt == 1);

  // NEW Kreiss-Oliger artificial viscosity
  int inArtVisc = 0;
  ppAmrlevelfactory.query("use_art_visc", inArtVisc);
  bool useArtVisc = (inArtVisc == 1);

  Real ratioArtVisc = 0.;
  if (useArtVisc)
    {
      ppAmrlevelfactory.query("ratio_art_visc", ratioArtVisc);
    }

  int inForwardEuler = 0;
  ppAmrlevelfactory.query("forward_euler", inForwardEuler);
  bool forwardEuler = (inForwardEuler == 1);

  // Initial values are average:  default false
  int inInitialAverage = 0;
  bool initialAverage;
  ppAmrlevelfactory.query("initial_average", inInitialAverage);
  initialAverage = (inInitialAverage == 1);

  // Do slope flattening:  default false
  int inFlattening = 0;
  bool useFlattening;
  ppAmrlevelfactory.query("use_flattening", inFlattening);
  useFlattening = (inFlattening == 1);

  // Avoid PPM:  default false
  int inNoPPM = 0;
  bool noPPM;
  ppAmrlevelfactory.query("no_ppm", inNoPPM);
  noPPM = (inNoPPM == 1);

  // Do deconvolution:  default true
  int inDoDeconvolution = 1;
  bool doDeconvolution;
  ppAmrlevelfactory.query("do_deconvolution", inDoDeconvolution);
  doDeconvolution = (inDoDeconvolution == 1);

  // Do deconvolution:  default true
  int inDoFaceDeconvolution = 1;
  bool doFaceDeconvolution;
  ppAmrlevelfactory.query("do_face_deconvolution", inDoFaceDeconvolution);
  doFaceDeconvolution = (inDoFaceDeconvolution == 1);

  // Apply artificial viscosity based on divergence.
  int inArtificialViscosity = 1;
  bool useArtificialViscosity;
  ppAmrlevelfactory.get("use_artificial_viscosity",inArtificialViscosity);
  useArtificialViscosity = (inArtificialViscosity == 1);

  // Artificial viscosity coefficient/multiplier
  Real artificialViscosity = 0.1;
  ppAmrlevelfactory.get("artificial_viscosity",artificialViscosity);

  // CFL multiplier
  Real cfl = 0.8;
  ppAmrlevelfactory.get("cfl",cfl);

  // Initial CFL multiplier
  Real initialCFL = 0.1;
  ppAmrlevelfactory.get("initial_cfl",initialCFL);

  /*
    Test problem inputs
  */

  // Read inputs that are prefixed with "godunov."
  ParmParse ppRun("run");

  // Stop after this number of steps
  int nstop = 0;
  ppRun.get("max_step",nstop);

  // Stop when the simulation time get here
  Real stopTime = 0.0;
  ppRun.get("max_time",stopTime);

  // Set the physical size of the longest dimension of the domain
  Real domainLength = 1.0;
  ppRun.get("domain_length",domainLength);

  // Set the resolution of the coarsest level
  vector<int> numCells(SpaceDim);
  for (int i = 0; i < SpaceDim; ++i)
  {
    numCells[i] = 0;
  }
  ppRun.getarr("num_cells",numCells,0,SpaceDim);
  // Check that every component of numCells is positive and even
  for (int i = 0; i < SpaceDim; ++i)
    {
      CH_assert(numCells[i] > 0);
      CH_assert(numCells[i] % 2 == 0);
    }

  // Determine which spatial directions are periodic
  vector<int> isPeriodica(SpaceDim,0);
  bool isPeriodic[SpaceDim];
  ppRun.getarr("is_periodic",isPeriodica,0,SpaceDim);
  // convert periodic from int to bool
  for (int dim = 0; dim < SpaceDim; dim++)
  {
    isPeriodic[dim] = (isPeriodica[dim] == 1);
    if (isPeriodic[dim] && verbosity >= 2 && procID() == 0)
    {
      pout() << "Using Periodic BCs in direction: " << dim << endl;
    }
  }

  /*
    Polytropic inputs
   */

  // Read inputs that are prefixed with "polytropic."
  ParmParse ppPolytropic("polytropic");

  // For all gas dynamics
  Real gamma = 1.4;
  ppPolytropic.query("gamma", gamma);

  Real M0sq=-1.;
  ppPolytropic.query("fourth_order_artificial_viscosity_parameter",M0sq);

  // End timing the reading of the input file
  TimeReadInput.stop();

#ifndef CH_NTIMER
  pout() << "Input Read completed --- "
         << "mem: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << get_memory_usage_from_OS()
         << " MB, time: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << TimeReadInput.wc_time()
         << " sec (wall-clock)" << endl;
#endif

  // Start timing AMR solver setup
  TimeSetupAMR.start();

  // Create and define IBC (initial and boundary condition) object
  PhysIBC* ibc;

  // A minimum pressure needed to construct MOLPolytropicPhysics - used in slope
  // flattening
  Real smallPressure;

  // Don't use source term by default
  bool useSourceTerm = false;

  // Source term multiplier
  Real sourceTermScaling = 0.0;

  ProblemDomain probDomain (IntVect::Zero,
                            IntVect(D_DECL6(numCells[0]-1,
                                            numCells[1]-1,
                                            numCells[2]-1,
                                            numCells[3]-1,
                                            numCells[4]-1,
                                            numCells[5]-1)),
                            isPeriodic);

  // Read inputs that are prefixed with "run.ibc."
  ParmParse ppRunIBC("run.ibc");

  // Determine the sample problem specified
  std::string problemString;
  if (ppRunIBC.contains("problem"))
  {
    ppRunIBC.query("problem",problemString);

    // Print some parameters
    if (verbosity >= 2)
    {
      pout() << "problem = " << problemString << endl;
      pout() << "gamma = " << gamma << endl;
    }

    if (problemString == "ramp")
    {
      if (isPeriodic[0] || isPeriodic[1])
      {
        MayDay::Error("Neither x or y boundaries can be periodic");
      }

      // Ramp problem
      Real alpha = 30.0;
      ppRunIBC.get("angle_deg",alpha);

      Real ms = 10.0;
      ppRunIBC.get("shock_mach",ms);

      Real xcorner = 0.1;
      ppRunIBC.get("xcorner",xcorner);

      if (verbosity >= 2)
      {
        pout() << "alpha = " << alpha << endl;
        pout() << "shock_mach = " << ms << endl;
        pout() << "xcorner = " << xcorner << endl;
      }

      // Define IBC for ramp problem
      RampIBC* rampibc = new RampIBC;
      rampibc->setFortranCommon(smallPressure,
                                gamma,
                                alpha,
                                ms,
                                xcorner,
                                artificialViscosity);
      ibc = rampibc;
    }
    else if (problemString == "channelShock")
    {
      if (isPeriodic[0])
      {
        MayDay::Error("The x boundary can't be periodic");
      }

      // Shock problem in a channel
      Real ms = 10.0;
      ppRunIBC.get("shock_mach",ms);

      Real position = 0.1;
      ppRunIBC.get("position",position);

      if (verbosity >= 2)
      {
        pout() << "shock_mach = " << ms << endl;
        pout() << "position = " << position << endl;
      }

      // Define IBC for channel shock problem
      ChannelShockIBC* channelshockibc = new ChannelShockIBC;
      channelshockibc->setFortranCommon(smallPressure,
                                        gamma,
                                        ms,
                                        position,
                                        artificialViscosity);
      ibc = channelshockibc;
    }
    else if (problemString == "channelModiano")
    {
      if (isPeriodic[0])
      {
        MayDay::Error("The x boundary can't be periodic");
      }

      // Modiano problem in a channel
      Real ambientDensity = 1.4;
      ppRunIBC.get("ambient_density",ambientDensity);

      Real deltaDensity = 0.014;
      ppRunIBC.get("delta_density",deltaDensity);

      Real center = 0.5;
      ppRunIBC.get("center",center);

      Real width = 0.5;
      ppRunIBC.get("width",width);

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "delta_density = " << deltaDensity << endl;
        pout() << "center = " << center << endl;
        pout() << "width = " << width << endl;
      }

      // Define IBC for channel Modiano problem
      ChannelModianoIBC* channelmodianoibc = new ChannelModianoIBC;
      channelmodianoibc->setFortranCommon(smallPressure,
                                          gamma,
                                          ambientDensity,
                                          deltaDensity,
                                          center,
                                          width,
                                          artificialViscosity);
      ibc = channelmodianoibc;
    }
    else if (problemString == "explosion")
    {
      // Explosion problem
      Real ms = 10.0;
      ppRunIBC.get("shock_mach",ms);

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppRunIBC.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      Real size = 0.25;
      ppRunIBC.get("initial_size",size);

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppRunIBC.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "shock_mach = " << ms << endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                center[1] << "  " <<,
                                                center[2] << "  " <<,
                                                center[3] << "  " <<,
                                                center[4] << "  " <<,
                                                center[5] << ) endl;
        pout() << "initial_size = " << size << endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                  velocity[1] << "  " <<,
                                                  velocity[2] << "  " <<,
                                                  velocity[3] << "  " <<,
                                                  velocity[4] << "  " <<,
                                                  velocity[5] << ) endl;
      }

      // Define IBC for explosion problem
      ExplosionIBC* explosionibc = new ExplosionIBC;
      explosionibc->setFortranCommon(smallPressure,
                                     gamma,
                                     ms,
                                     center,
                                     size,
                                     velocity,
                                     artificialViscosity);
      ibc = explosionibc;
    }
    else if (problemString == "explosion1d")
    {
      // Explosion1d problem
      Real ms = 10.0;
      ppRunIBC.get("shock_mach",ms);

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppRunIBC.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      Real size = 0.25;
      ppRunIBC.get("initial_size",size);

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppRunIBC.getarr("initial_velocity",velocitypp,0,SpaceDim);
      // use only component 0; set the other components to 0.
      velocity[0] = velocitypp[0];
      for (int i = 1; i < SpaceDim; i++)
      {
        velocity[i] = 0.;
      }

      if (verbosity >= 2)
      {
        pout() << "shock_mach = " << ms << endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                center[1] << "  " <<,
                                                center[2] << "  " <<,
                                                center[3] << "  " <<,
                                                center[4] << "  " <<,
                                                center[5] << ) endl;
        pout() << "initial_size = " << size << endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                  velocity[1] << "  " <<,
                                                  velocity[2] << "  " <<,
                                                  velocity[3] << "  " <<,
                                                  velocity[4] << "  " <<,
                                                  velocity[5] << ) endl;
      }

      // Define IBC for explosion problem
      Explosion1dIBC* explosion1dibc = new Explosion1dIBC;
      explosion1dibc->setFortranCommon(smallPressure,
                                       gamma,
                                       ms,
                                       center,
                                       size,
                                       velocity,
                                       artificialViscosity);
      ibc = explosion1dibc;
    }
    else if (problemString == "shockTube")
    {
      ParmParse ppshocktube("shockTube");
      Real center;
      ppshocktube.get("initial_center",center);

      int direction = 0;
      ppshocktube.query("direction", direction);

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppshocktube.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
        {
          velocity[i] = velocitypp[i];
        }

      // define left and right states
      Real rhoLeft, rhoRight, eLeft, eRight;
      ppshocktube.get("rhoLeft", rhoLeft);
      ppshocktube.get("rhoRight", rhoRight);
      ppshocktube.get("eLeft", eLeft);
      ppshocktube.get("eRight", eRight);

      if (verbosity >= 2)
        {
          pout() << "rhoLeft = " << rhoLeft << endl;
          pout() << "rhoRight = " << rhoRight << endl;
          pout() << "eLeft = " << eLeft << endl;
          pout() << "eRight = " << eRight << endl;

          pout() << "initial_center = " << center << endl;
          pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                    velocity[1] << "  " <<,
                                                    velocity[2] << "  " <<,
                                                    velocity[3] << "  " <<,
                                                    velocity[4] << "  " <<,
                                                    velocity[5] << ) endl;
        }

      ShockTubeIBC* shocktubeibc = new ShockTubeIBC;
      shocktubeibc->setFortranCommon(smallPressure,
                                     gamma,
                                     rhoLeft, rhoRight,
                                     eLeft, eRight,
                                     center,
                                     direction,
                                     velocity,
                                     artificialViscosity);
      ibc = shocktubeibc;
    }
    else if (problemString == "gaussian")
    {
      // Gaussian problem
      Real ambientDensity = 1.4;
      ppRunIBC.get("ambient_density",ambientDensity);

      Real deltaDensity = 0.014;
      ppRunIBC.get("delta_density",deltaDensity);

      int pressure = -1;
      std::string pressureString;
      ppRunIBC.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
      {
        pressure = PRESSURE_ISENTROPIC;
      }
      else if (pressureString == "constant")
      {
        pressure = PRESSURE_CONSTANT;
      }

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString << "\", specified in input file" << endl << endl;
        return;
      }

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppRunIBC.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      Real size = 0.25;
      ppRunIBC.get("initial_size",size);

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppRunIBC.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "delta_density = " << deltaDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                center[1] << "  " <<,
                                                center[2] << "  " <<,
                                                center[3] << "  " <<,
                                                center[4] << "  " <<,
                                                center[5] << ) endl;
        pout() << "initial_size = " << size << endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                  velocity[1] << "  " <<,
                                                  velocity[2] << "  " <<,
                                                  velocity[3] << "  " <<,
                                                  velocity[4] << "  " <<,
                                                  velocity[5] << ) endl;
      }

      // Define IBC for gaussian problem
      GaussianIBC* gaussianibc = new GaussianIBC;
      gaussianibc->setFortranCommon(smallPressure,
                                    gamma,
                                    ambientDensity,
                                    deltaDensity,
                                    pressure,
                                    center,
                                    size,
                                    velocity,
                                    artificialViscosity);
      ibc = gaussianibc;
    }
    else if (problemString == "gaussiansmooth")
    {
      // smoothed Gaussian problem
      Real ambientDensity = 1.4;
      ppRunIBC.get("ambient_density",ambientDensity);

      Real deltaDensity = 0.014;
      ppRunIBC.get("delta_density",deltaDensity);

      int pressure = -1;
      std::string pressureString;
      ppRunIBC.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
      {
        pressure = PRESSURE_ISENTROPIC;
      }
      else if (pressureString == "constant")
      {
        pressure = PRESSURE_CONSTANT;
      }

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString << "\", specified in input file" << endl << endl;
        return;
      }

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppRunIBC.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      // domainWidth = physical width of domain
      RealVect domainWidth = RealVect::Unit;
      int minWidthDir = domainWidth.minDir(true);
      Real widthMin = domainWidth[minWidthDir];

      Real size = 0.25;
      ppRunIBC.get("initial_size",size);

      Real outerFlatBuffer = 0.;
      ppRunIBC.query("initial_buffer", outerFlatBuffer);
      Real radmax = (1.0 - outerFlatBuffer) * 0.5 * widthMin;

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppRunIBC.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "delta_density = " << deltaDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                center[1] << "  " <<,
                                                center[2] << "  " <<,
                                                center[3] << "  " <<,
                                                center[4] << "  " <<,
                                                center[5] << ) endl;
        pout() << "initial_buffer = " << outerFlatBuffer << endl;
        pout() << "initial_size = " << size << endl;
        pout() << "initial_max_radius = " << radmax << endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                  velocity[1] << "  " <<,
                                                  velocity[2] << "  " <<,
                                                  velocity[3] << "  " <<,
                                                  velocity[4] << "  " <<,
                                                  velocity[5] << ) endl;
      }

      // Define IBC for gaussian problem
      GaussianSmoothBC* gaussianibc = new GaussianSmoothBC;
      gaussianibc->setFortranCommon(smallPressure,
                                    gamma,
                                    ambientDensity,
                                    deltaDensity,
                                    pressure,
                                    center,
                                    size,
                                    radmax,
                                    velocity,
                                    artificialViscosity);
      ibc = gaussianibc;
    }
    else if (problemString == "cosinepower")
    {
      // cosine power problem
      Real ambientDensity = 1.4;
      ppRunIBC.get("ambient_density",ambientDensity);

      Real deltaDensity = 0.014;
      ppRunIBC.get("delta_density",deltaDensity);

      int pressure = -1;
      std::string pressureString;
      ppRunIBC.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
      {
        pressure = PRESSURE_ISENTROPIC;
      }
      else if (pressureString == "constant")
      {
        pressure = PRESSURE_CONSTANT;
      }

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString << "\", specified in input file" << endl << endl;
        return;
      }

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppRunIBC.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      // domainWidth = physical width of domain
      RealVect domainWidth = RealVect::Unit;
      int minWidthDir = domainWidth.minDir(true);
      Real widthMin = domainWidth[minWidthDir];

      int cosinePower = 6;
      ppRunIBC.query("cosine_power", cosinePower);

      Real outerFlatBuffer = 0.;
      ppRunIBC.query("initial_buffer", outerFlatBuffer);
      Real radmax = (1.0 - outerFlatBuffer) * 0.5 * widthMin;

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppRunIBC.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "delta_density = " << deltaDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                center[1] << "  " <<,
                                                center[2] << "  " <<,
                                                center[3] << "  " <<,
                                                center[4] << "  " <<,
                                                center[5] << ) endl;
        pout() << "initial_buffer = " << outerFlatBuffer << endl;
        pout() << "initial_max_radius = " << radmax << endl;
        pout() << "initial_cosine_power = " << cosinePower << endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                  velocity[1] << "  " <<,
                                                  velocity[2] << "  " <<,
                                                  velocity[3] << "  " <<,
                                                  velocity[4] << "  " <<,
                                                  velocity[5] << ) endl;
      }

      // Define IBC for cosine-power problem
      CosinePowerIBC* cospoweribc = new CosinePowerIBC;
      cospoweribc->setFortranCommon(smallPressure,
                                    gamma,
                                    ambientDensity,
                                    deltaDensity,
                                    pressure,
                                    center,
                                    radmax,
                                    cosinePower,
                                    velocity,
                                    artificialViscosity);
      ibc = cospoweribc;
    }
    else if (problemString == "cosineeachpower")
    {
      // cosine power problem
      Real ambientDensity = 1.4;
      ppRunIBC.get("ambient_density",ambientDensity);

      Real deltaDensity = 0.014;
      ppRunIBC.get("delta_density",deltaDensity);

      int pressure = -1;
      std::string pressureString;
      ppRunIBC.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
      {
        pressure = PRESSURE_ISENTROPIC;
      }
      else if (pressureString == "constant")
      {
        pressure = PRESSURE_CONSTANT;
      }

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString << "\", specified in input file" << endl << endl;
        return;
      }

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppRunIBC.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      // domainWidth = physical width of domain
      RealVect domainWidth = RealVect::Unit;
      int minWidthDir = domainWidth.minDir(true);
      Real widthMin = domainWidth[minWidthDir];

      int cosinePower = 6;
      ppRunIBC.query("cosine_power", cosinePower);

      Real outerFlatBuffer = 0.;
      ppRunIBC.query("initial_buffer", outerFlatBuffer);
      Real radmax = (1.0 - outerFlatBuffer) * 0.5 * widthMin;

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppRunIBC.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "delta_density = " << deltaDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                center[1] << "  " <<,
                                                center[2] << "  " <<,
                                                center[3] << "  " <<,
                                                center[4] << "  " <<,
                                                center[5] << ) endl;
        pout() << "initial_buffer = " << outerFlatBuffer << endl;
        pout() << "initial_max_radius = " << radmax << endl;
        pout() << "initial_cosine_power = " << cosinePower << endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                  velocity[1] << "  " <<,
                                                  velocity[2] << "  " <<,
                                                  velocity[3] << "  " <<,
                                                  velocity[4] << "  " <<,
                                                  velocity[5] << ) endl;
      }

      // Define IBC for cosine-power problem
      CosineEachPowerIBC* cospoweribc = new CosineEachPowerIBC;
      cospoweribc->setFortranCommon(smallPressure,
                                    gamma,
                                    ambientDensity,
                                    deltaDensity,
                                    pressure,
                                    center,
                                    radmax,
                                    cosinePower,
                                    velocity,
                                    artificialViscosity);
      ibc = cospoweribc;
    }
    else if (problemString == "vortex")
    {
      Real ambientDensity = 1.4;
      ppRunIBC.get("ambient_density",ambientDensity);

      int pressure = -1;
      std::string pressureString;
      ppRunIBC.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
      {
        pressure = PRESSURE_ISENTROPIC;
      }
      else if (pressureString == "constant")
      {
        pressure = PRESSURE_CONSTANT;
      }

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString
               << "\", specified in input file" << endl << endl;
        return;
      }

      std::string vortexFile;
      ppRunIBC.query("vortex_file", vortexFile);

      Vector<RealVect> center;
      Vector<Real> radius;
      Vector<Real> strength;
      setupVortices(center, radius, strength, verbosity, vortexFile);

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        for (int iVort = 0; iVort < center.size(); iVort++)
          {
            RealVect& thisCenter = center[iVort];
            pout() << "vortex_center = " << D_TERM6(thisCenter[0] << "  " <<,
                                                   thisCenter[1] << "  " <<,
                                                   thisCenter[2] << "  " <<,
                                                   thisCenter[3] << "  " <<,
                                                   thisCenter[4] << "  " <<,
                                                   thisCenter[5] << ) endl;
            pout() << "vortex_radius = " << radius[iVort] << endl;
            pout() << "vortex_strength = " << strength[iVort] << endl;
          }
      }

      // Define IBC for vortex problem
      VortexIBC* vortexibc = new VortexIBC;
      vortexibc->setFortranCommon(smallPressure,
                                  gamma,
                                  ambientDensity,
                                  pressure,
                                  center,
                                  radius,
                                  strength,
                                  artificialViscosity);
      ibc = vortexibc;
    }
    else if (problemString == "ZZSvortex")
      { // Vortex of Zhang, Zhang, and Shu
      Real ambientDensity = 1.;
      ppRunIBC.get("ambient_density",ambientDensity);

      int pressure = PRESSURE_ISENTROPIC;

      Vector<RealVect> center(1);
      Vector<Real> radius(1);
      Vector<Real> strength(1);

      center[0] = RealVect(D_DECL6(5., 5., 0., 0., 0., 0.));
      radius[0] = 1.;
      strength[0] = 5.;

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        for (int iVort = 0; iVort < center.size(); iVort++)
          {
            RealVect& thisCenter = center[iVort];
            pout() << "vortex_center = " << D_TERM6(thisCenter[0] << "  " <<,
                                                   thisCenter[1] << "  " <<,
                                                   thisCenter[2] << "  " <<,
                                                   thisCenter[3] << "  " <<,
                                                   thisCenter[4] << "  " <<,
                                                   thisCenter[5] << ) endl;
            pout() << "vortex_radius = " << radius[iVort] << endl;
            pout() << "vortex_strength = " << strength[iVort] << endl;
          }
      }

      // Define IBC for vortex problem
      ZZSVortexIBC* vortexibc = new ZZSVortexIBC;
      vortexibc->setFortranCommon(smallPressure,
                                  gamma,
                                  ambientDensity,
                                  pressure,
                                  center,
                                  radius,
                                  strength,
                                  artificialViscosity);
      ibc = vortexibc;
    }
    else if (problemString == "shear")
    {
      // shear-layer problem
      Real ambientDensity = 1.4;
      ppRunIBC.get("ambient_density",ambientDensity);

      int pressure = -1;
      std::string pressureString;
      ppRunIBC.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
      {
        pressure = PRESSURE_ISENTROPIC;
      }
      else if (pressureString == "constant")
      {
        pressure = PRESSURE_CONSTANT;
      }

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString
               << "\", specified in input file" << endl << endl;
        return;
      }

      Real widthinv = 42.;
      ppRunIBC.get("shear_widthinv", widthinv);

      Real delt = 0.05;
      ppRunIBC.get("shear_delt", delt);

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        pout() << "shear_widthinv = " << widthinv << endl;
        pout() << "shear_delt = " << delt << endl;
      }

      // Define IBC for vortex problem
      ShearIBC* shearibc = new ShearIBC;
      shearibc->setFortranCommon(smallPressure,
                                 gamma,
                                 ambientDensity,
                                 pressure,
                                 widthinv,
                                 delt,
                                 artificialViscosity);
      ibc = shearibc;
    }
    else if (problemString == "trig")
    {
      // trigonometric problem
      Real ambientDensity = 1.4;
      ppRunIBC.get("ambient_density",ambientDensity);

      int pressure = -1;
      std::string pressureString;
      ppRunIBC.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
      {
        pressure = PRESSURE_ISENTROPIC;
      }
      else if (pressureString == "constant")
      {
        pressure = PRESSURE_CONSTANT;
      }

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString
               << "\", specified in input file" << endl << endl;
        return;
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
      }

      // Define IBC for vortex problem
      TrigIBC* trigibc = new TrigIBC;
      trigibc->setFortranCommon(smallPressure,
                                gamma,
                                ambientDensity,
                                pressure,
                                artificialViscosity);
      ibc = trigibc;
    }
    else if (problemString == "gaussian1d")
    {
      // 1D Gaussian problem
      Real ambientDensity = 1.4;
      ppRunIBC.get("ambient_density",ambientDensity);

      Real deltaDensity = 0.014;
      ppRunIBC.get("delta_density",deltaDensity);

      int pressure = -1;
      std::string pressureString;
      ppRunIBC.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
      {
        pressure = PRESSURE_ISENTROPIC;
      }
      else if (pressureString == "constant")
      {
        pressure = PRESSURE_CONSTANT;
      }

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString << "\", specified in input file" << endl << endl;
        return;
      }

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppRunIBC.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      Real size = 0.25;
      ppRunIBC.get("initial_size",size);

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppRunIBC.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "delta_density = " << deltaDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                center[1] << "  " <<,
                                                center[2] << "  " <<,
                                                center[3] << "  " <<,
                                                center[4] << "  " <<,
                                                center[5] << ) endl;
        pout() << "initial_size = " << size << endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                  velocity[1] << "  " <<,
                                                  velocity[2] << "  " <<,
                                                  velocity[3] << "  " <<,
                                                  velocity[4] << "  " <<,
                                                  velocity[5] << ) endl;
      }

      // Define IBC for gaussian problem
      Gaussian1dIBC* gaussianibc = new Gaussian1dIBC;
      gaussianibc->setFortranCommon(smallPressure,
                                    gamma,
                                    ambientDensity,
                                    deltaDensity,
                                    pressure,
                                    center,
                                    size,
                                    velocity,
                                    artificialViscosity);
      ibc = gaussianibc;
    }
    else if (problemString == "gaussian1dv")
    {
      // 1D Gaussian problem
      Real ambientDensity = 1.4;
      ppRunIBC.get("ambient_density",ambientDensity);

      Real deltaDensity = 0.014;
      ppRunIBC.get("delta_density",deltaDensity);

      int pressure = -1;
      std::string pressureString;
      ppRunIBC.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
      {
        pressure = PRESSURE_ISENTROPIC;
      }
      else if (pressureString == "constant")
      {
        pressure = PRESSURE_CONSTANT;
      }

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString << "\", specified in input file" << endl << endl;
        return;
      }

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppRunIBC.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      Real size = 0.25;
      ppRunIBC.get("initial_size",size);

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppRunIBC.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "delta_density = " << deltaDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                center[1] << "  " <<,
                                                center[2] << "  " <<,
                                                center[3] << "  " <<,
                                                center[4] << "  " <<,
                                                center[5] << ) endl;
        pout() << "initial_size = " << size << endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                  velocity[1] << "  " <<,
                                                  velocity[2] << "  " <<,
                                                  velocity[3] << "  " <<,
                                                  velocity[4] << "  " <<,
                                                  velocity[5] << ) endl;
      }

      // Define IBC for gaussian problem
      Gaussian1dvIBC* gaussianibc = new Gaussian1dvIBC;
      gaussianibc->setFortranCommon(smallPressure,
                                    gamma,
                                    ambientDensity,
                                    deltaDensity,
                                    pressure,
                                    center,
                                    size,
                                    velocity,
                                    artificialViscosity);
      ibc = gaussianibc;
    }
    else if (problemString == "gaussianp")
    {
      // Gaussian pressure problem
      Real ambientDensity = 1.4;
      ppRunIBC.get("ambient_density",ambientDensity);

      Real ambientPressure = 1.4;
      ppRunIBC.get("ambient_pressure",ambientPressure);

      Real deltaPressure = 0.014;
      ppRunIBC.get("delta_pressure",deltaPressure);

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppRunIBC.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      Real size = 0.25;
      ppRunIBC.get("initial_size",size);

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppRunIBC.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "ambient_pressure = " << ambientPressure << endl;
        pout() << "delta_pressure = " << deltaPressure << endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                center[1] << "  " <<,
                                                center[2] << "  " <<,
                                                center[3] << "  " <<,
                                                center[4] << "  " <<,
                                                center[5] << ) endl;
        pout() << "initial_size = " << size << endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                  velocity[1] << "  " <<,
                                                  velocity[2] << "  " <<,
                                                  velocity[3] << "  " <<,
                                                  velocity[4] << "  " <<,
                                                  velocity[5] << ) endl;
      }

      // Define IBC for gaussian problem
      GaussianPBC* gaussianpbc = new GaussianPBC;
      gaussianpbc->setFortranCommon(smallPressure,
                                    gamma,
                                    ambientDensity,
                                    ambientPressure,
                                    deltaPressure,
                                    center,
                                    size,
                                    velocity,
                                    artificialViscosity);
      ibc = gaussianpbc;
    }
    else if (problemString == "sine")
    {
      // Sine problem
      Real ambientDensity = 1.4;
      ppRunIBC.get("ambient_density",ambientDensity);

      Real deltaDensity = 0.014;
      ppRunIBC.get("delta_density",deltaDensity);

      int pressure = -1;
      std::string pressureString;
      ppRunIBC.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
      {
        pressure = PRESSURE_ISENTROPIC;
      }
      else if (pressureString == "constant")
      {
        pressure = PRESSURE_CONSTANT;
      }

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString << "\", specified in input file" << endl << endl;
        return;
      }

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppRunIBC.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppRunIBC.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "delta_density = " << deltaDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                center[1] << "  " <<,
                                                center[2] << "  " <<,
                                                center[3] << "  " <<,
                                                center[4] << "  " <<,
                                                center[5] << ) endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                  velocity[1] << "  " <<,
                                                  velocity[2] << "  " <<,
                                                  velocity[3] << "  " <<,
                                                  velocity[4] << "  " <<,
                                                  velocity[5] << ) endl;
      }

      // Define IBC for sine problem
      SineIBC* sineibc = new SineIBC;
      sineibc->setFortranCommon(smallPressure,
                                gamma,
                                ambientDensity,
                                deltaDensity,
                                pressure,
                                center,
                                velocity,
                                artificialViscosity);
      ibc = sineibc;
    }
    else if (problemString == "wave")
    {
      // Plane wave problem
      Real ambientDensity = 1.4;
      ppRunIBC.get("ambient_density",ambientDensity);

      Real deltaDensity = 0.014;
      ppRunIBC.get("delta_density",deltaDensity);

      int pressure = -1;
      std::string pressureString;
      ppRunIBC.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
        {
          pressure = PRESSURE_ISENTROPIC;
        }
      else if (pressureString == "constant")
        {
          pressure = PRESSURE_CONSTANT;
        }

      if (pressure == -1)
        {
          pout() << "Invalid pressure, \"" << pressureString << "\", specified in input file" << endl << endl;
          return;
        }
      // petermc, 7 July 2008:  allow PRESSURE_CONSTANT.
      //      else
      //        {
      //          pressure = PRESSURE_ISENTROPIC;
      //        }

      vector<int> waveNumberpp(SpaceDim,0);
      waveNumberpp[0] = 1;
      IntVect waveNumber;
      ppRunIBC.getarr("wave_number",waveNumberpp,0,SpaceDim);
      int norm2 = 0;
      for (int i = 0; i < SpaceDim; i++)
      {
        waveNumber[i] = waveNumberpp[i];
        norm2 += waveNumber[i]*waveNumber[i];
      }

      if (norm2 == 0)
      {
        pout() << "One component of the wave number must be non-zero" << endl << endl;
        return;
      }

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppRunIBC.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppRunIBC.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "delta_density = " << deltaDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        pout() << "wave_number = " << D_TERM6(waveNumber[0] << "  " <<,
                                             waveNumber[1] << "  " <<,
                                             waveNumber[2] << "  " <<,
                                             waveNumber[3] << "  " <<,
                                             waveNumber[4] << "  " <<,
                                             waveNumber[5] << ) endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                center[1] << "  " <<,
                                                center[2] << "  " <<,
                                                center[3] << "  " <<,
                                                center[4] << "  " <<,
                                                center[5] << ) endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                  velocity[1] << "  " <<,
                                                  velocity[2] << "  " <<,
                                                  velocity[3] << "  " <<,
                                                  velocity[4] << "  " <<,
                                                  velocity[5] << ) endl;
      }

      // Define IBC for plane wave problem
      WaveIBC* waveibc = new WaveIBC;
      waveibc->setFortranCommon(smallPressure,
                                gamma,
                                ambientDensity,
                                deltaDensity,
                                pressure,
                                waveNumber,
                                center,
                                velocity,
                                artificialViscosity);
      ibc = waveibc;
    }
    else if (problemString == "source")
    {
      // Plane wave problem with a source
      Real ambientDensity = 1.4;
      ppRunIBC.get("ambient_density",ambientDensity);

      Real deltaDensity = 0.014;
      ppRunIBC.get("delta_density",deltaDensity);

      int pressure = -1;
      std::string pressureString;
      ppRunIBC.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
        {
          pressure = PRESSURE_ISENTROPIC;
        }
      else if (pressureString == "constant")
        {
          pressure = PRESSURE_CONSTANT;
        }

      useSourceTerm = true;
      ppRunIBC.get("source_term_scaling",sourceTermScaling);

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString << "\", specified in input file" << endl << endl;
        return;
      }
      else
      {
        pressure = PRESSURE_ISENTROPIC;
      }

      vector<int> waveNumberpp(SpaceDim,0);
      waveNumberpp[0] = 1;
      IntVect waveNumber;
      ppRunIBC.getarr("wave_number",waveNumberpp,0,SpaceDim);
      int norm2 = 0;
      for (int i = 0; i < SpaceDim; i++)
      {
        waveNumber[i] = waveNumberpp[i];
        norm2 += waveNumber[i]*waveNumber[i];
      }

      if (norm2 == 0)
      {
        pout() << "One component of the wave number must be non-zero" << endl << endl;
        return;
      }

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppRunIBC.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppRunIBC.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "delta_density = " << deltaDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        pout() << "wave_number = " << D_TERM6(waveNumber[0] << "  " <<,
                                             waveNumber[1] << "  " <<,
                                             waveNumber[2] << "  " <<,
                                             waveNumber[3] << "  " <<,
                                             waveNumber[4] << "  " <<,
                                             waveNumber[5] << ) endl;
        pout() << "initial_center = " << D_TERM6(center[0] << "  " <<,
                                                center[1] << "  " <<,
                                                center[2] << "  " <<,
                                                center[3] << "  " <<,
                                                center[4] << "  " <<,
                                                center[5] << ) endl;
        pout() << "initial_velocity = " << D_TERM6(velocity[0] << "  " <<,
                                                  velocity[1] << "  " <<,
                                                  velocity[2] << "  " <<,
                                                  velocity[3] << "  " <<,
                                                  velocity[4] << "  " <<,
                                                  velocity[5] << ) endl;
      }

      // Define IBC for plane wave with source problem
      SourceIBC* sourceibc = new SourceIBC;
      sourceibc->setFortranCommon(smallPressure,
                                  gamma,
                                  ambientDensity,
                                  deltaDensity,
                                  pressure,
                                  waveNumber,
                                  center,
                                  velocity,
                                  artificialViscosity);
      ibc = sourceibc;
    }
    else
    {
      // The sample problem name given isn't valid
      pout() << "Invalid problem, \"" << problemString << "\", specified in input file" << endl << endl;
      return;
    }
  }
  else
  {
    // A sample problem must be specified
    pout() << "\"run.ibc.problem\" not specified in input file" << endl << endl;
    return;
  }

  if (verbosity >= 2)
  {
    pout() << "verbosity = " << verbosity << endl;

    pout() << "maximum_step = " << nstop << endl;
    pout() << "maximum_time = " << stopTime << endl;
    if (fixedDt > 0)
    {
      pout() << "fixed_dt = " << fixedDt << endl;
    }

    pout() << "number_of_cells = " << D_TERM6(numCells[0] << "  " <<,
                                             numCells[1] << "  " <<,
                                             numCells[2] << "  " <<,
                                             numCells[3] << "  " <<,
                                             numCells[4] << "  " <<,
                                             numCells[5] << ) endl;
    pout() << "is_period = " << D_TERM6(isPeriodic[0] << "  " <<,
                                       isPeriodic[1] << "  " <<,
                                       isPeriodic[2] << "  " <<,
                                       isPeriodic[3] << "  " <<,
                                       isPeriodic[4] << "  " <<,
                                       isPeriodic[5] << ) endl;

    pout() << "maximum_level = " << maxLevel << endl;
    pout() << "refinement_ratio = ";
    for (int i = 0; i < refRatios.size(); ++i)
    {
      pout() << refRatios[i] << " ";
    }
    pout() << endl;

    pout() << "regrid_interval = ";
    for (int i = 0; i < regridIntervals.size(); ++i)
    {
      pout() << regridIntervals[i] << " ";
    }
    pout() << endl;
    pout() << "tag_buffer_size = " << tagBufferSize << endl;

    pout() << "refinement_threshold = " << refineThresh << endl;

    pout() << "blocking_factor = " << blockFactor << endl;
    pout() << "max_grid_size = " << maxGridSize << endl;
    pout() << "fill_ratio = " << fillRatio << endl;

    pout() << "slope_order = "
           << (useFourthOrderSlopes ? "4th" : "2nd") << endl;
    pout() << "use_primitive_slope_limiting = "
           << (usePrimLimiting ? "yes" : "no") << endl;
    pout() << "initial_average = "
           << (initialAverage ? "yes" : "no") << endl;
    pout() << "use_slope_flattening = "
           << (useFlattening ? "yes" : "no") << endl;

    pout() << "use_artificial_viscosity = "
           << (useArtificialViscosity ? "yes" : "no") << endl;
    if (useArtificialViscosity)
    {
      pout() << "artificial_viscosity = " << artificialViscosity << endl;
    }

    pout() << "use_source_term = "
           << (useSourceTerm ? "yes" : "no") << endl;
    if (useSourceTerm)
    {
      pout() << "source_term_scaling = " << sourceTermScaling << endl;
    }

    pout() << "checkpoint_interval = " << checkpointInterval << endl;
    pout() << "plot_interval = " << plotInterval << endl;

    pout() << "CFL = " << cfl << endl;
    pout() << "initial_CFL = " << initialCFL << endl;

    pout() << "maximum_dt_growth = " << maxDtGrowth << endl;
    pout() << "dt_tolerance_factor = " << dtToleranceFactor << endl;
  }

  // Set up the physics for polytropic gas dynamics
  MOLPolytropicPhysics polytropicPhysics(smallPressure);
  polytropicPhysics.setPhysIBC(ibc);
  if (M0sq > 0.)
    polytropicPhysics.setFourthOrderArtificialViscosityParameter(M0sq);
  pout() << "fourth-order artificial viscosity parameters = " << M0sq << endl;
  // Cast to physics base class pointer for technical reasons
  MOLPhysics* molPhysics = static_cast<MOLPhysics*> (&polytropicPhysics);

  int spaceOrder = (useFourthOrderSlopes) ? 4 : 2;

  // Set up the AMRLevel... factory
  AMRLevelConsFactory amrLevelFact;
  amrLevelFact.spaceOrder(spaceOrder);
  amrLevelFact.limitFaceValues(usePrimLimiting); // WAS limitFaceValues
  amrLevelFact.highOrderLimiter(highOrderLimiter); // WAS limitFaceValues
  amrLevelFact.initialAverage(initialAverage);
  amrLevelFact.useFlattening(useFlattening);
  amrLevelFact.noPPM(noPPM);
  amrLevelFact.doDeconvolution(doDeconvolution);
  amrLevelFact.doFaceDeconvolution(doFaceDeconvolution);
  amrLevelFact.useArtificialViscosity(useArtificialViscosity);
  amrLevelFact.artificialViscosity(artificialViscosity);
  amrLevelFact.useArtVisc(useArtVisc);
  amrLevelFact.ratioArtVisc(ratioArtVisc);
  amrLevelFact.forwardEuler(forwardEuler);
  // amrLevelFact.enforceMinVal(enforceMinVal, minVal);
  amrLevelFact.CFL(cfl);
  amrLevelFact.domainLength(domainLength);
  amrLevelFact.refinementThreshold(refineThresh);
  amrLevelFact.refinementIsScaled(refinementIsScaled);
  amrLevelFact.tagPressure(tagPressure);
  amrLevelFact.tagVorticity(tagVorticity);
  amrLevelFact.tagBufferSize(tagBufferSize);
  amrLevelFact.verbosity(verbosity);
  amrLevelFact.initialDtMultiplier(initialCFL);
  // amrLevelFact.IBC(ibc);
  amrLevelFact.molPhysics(molPhysics);

  { // scope of AMR amr;
  AMR amr;

  // Set up the AMR object
  amr.define(maxLevel, refRatios, probDomain, &amrLevelFact);

  if (fixedDt > 0)
  {
    amr.fixedDt(fixedDt);
  }

  // Set grid generation parameters
  amr.maxGridSize(maxGridSize);
  amr.blockFactor(blockFactor);
  amr.fillRatio(fillRatio);

  amr.gridBufferSize(gridBufferSize);

  // Set output parameters
  amr.checkpointInterval(checkpointInterval);
  amr.plotInterval(plotInterval);
  amr.regridIntervals(regridIntervals);
  amr.maxDtGrow(maxDtGrowth);
  amr.dtToleranceFactor(dtToleranceFactor);
  amr.useSubcyclingInTime(useSubcycling);

  // Set up output files
  if (ppAmr.contains("plot_prefix"))
  {
    std::string prefix;
    ppAmr.query("plot_prefix",prefix);
    amr.plotPrefix(prefix);
  }

  if (ppAmr.contains("chk_prefix"))
  {
    std::string prefix;
    ppAmr.query("chk_prefix",prefix);
    amr.checkpointPrefix(prefix);
  }

  amr.verbosity(verbosity);

  // Set up input files
  if (!ppAmr.contains("restart_file"))
  {
    if (!ppAmr.contains("fixed_hierarchy"))
    {
      // initialize from scratch for AMR run
      // initialize hierarchy of levels
      amr.setupForNewAMRRun();
    }
    else
    {
      //      std::string gridFile;
      //      ppAmr.query("fixed_hierarchy",gridFile);

      // initialize from a list of grids in "gridFile"
      int numLevels = maxLevel+1;
      Vector<Vector<Box> > amrGrids(numLevels);
//       setupFixedGrids(amrGrids,
//                       probDomain,
//                       maxLevel,
//                       maxGridSize,
//                       blockFactor,
//                       verbosity,
//                       gridFile);
      GenFuncs::readBoxes(amrGrids, ppAmr, probDomain,
                maxGridSize, blockFactor,
                numLevels, refRatios, verbosity);
      amr.setupForFixedHierarchyRun(amrGrids,1);
    }
  }
  else // read from restart file
  {
    std::string restartFile;
    ppAmr.query("restart_file",restartFile);

#ifdef CH_USE_HDF5
    HDF5Handle handle(restartFile,HDF5Handle::OPEN_RDONLY);
    // read from checkpoint file
    amr.setupForRestart(handle);
    handle.close();
#else
    MayDay::Error("amrRun restart only defined with hdf5");
#endif
    // This section added by petermc, 20 Jan 2010
    if (ppAmr.contains("fixed_hierarchy"))
      {
        int numLevels = maxLevel+1;
        Vector<int> regridInterval(numLevels);
        for (int ilev = 0; ilev < numLevels; ilev++) regridInterval[ilev] = 0;
        amr.regridIntervals(regridInterval);
      }
  }

  // End timing AMR solver setup
  TimeSetupAMR.stop();

#ifndef CH_NTIMER
  pout() << "AMR Setup completed ---- "
         << "mem: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << get_memory_usage_from_OS()
         << " MB, time: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << TimeSetupAMR.wc_time()
         << " sec (wall-clock)" << endl;

  if (verbosity >= 1)
  {
    pout() << endl;
  }
#endif

  // Run and time the computation
  TimeRun.start();
  amr.run(stopTime,nstop);
  TimeRun.stop();

#ifndef CN_NTIMER
  if (verbosity >= 1)
  {
    pout() << endl;
  }

  pout() << "AMR Run completed ------ "
         << "mem: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << get_memory_usage_from_OS()
         << " MB, time: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << TimeRun.wc_time()
         << " sec (wall-clock)" << endl;
#endif

  // Output the last plot file and statistics - time the process
  TimeConclude.start();
  amr.conclude();
  TimeConclude.stop();

#ifndef CH_NTIMER
  pout() << "AMR Conclude completed - "
         << "mem: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << get_memory_usage_from_OS()
         << " MB, time: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << TimeConclude.wc_time()
         << " sec (wall-clock)" << endl;
#endif
  } // scope of amr
}

// setupFixedGrids allows fixed grids to be read in and used in this AMR
// computation example
void setupFixedGrids(Vector<Vector<Box> >& a_amrGrids,
                     const ProblemDomain&  a_domain,
                     int                   a_maxLevel,
                     int                   a_maxGridSize,
                     int                   a_blockFactor,
                     int                   a_verbosity,
                     std::string           a_gridFile)
{
  // Run this task on one processor
  if (procID() == uniqueProc(SerialTask::compute))
  {
    a_amrGrids.push_back(Vector<Box>(1,a_domain.domainBox()));

    // Read in predefined grids
    ifstream is(a_gridFile.c_str(), ios::in);

    if (is.fail())
    {
      MayDay::Error("Cannot open grids file");
    }

    // Format of file:
    //   number of levels, then for each level (starting with level 1):
    //   number of grids on level, list of boxes

    int inNumLevels;
    is >> inNumLevels;

    CH_assert (inNumLevels <= a_maxLevel+1);

    if (a_verbosity >= 3)
    {
      pout() << "numLevels = " << inNumLevels << endl;
    }

    while (is.get() != '\n');

    a_amrGrids.resize(inNumLevels);

    // Check to see if coarsest level needs to be broken up
    domainSplit(a_domain,a_amrGrids[0],a_maxGridSize,a_blockFactor);

    if (a_verbosity >= 3)
    {
      pout() << "level 0: ";
      for (int n = 0; n < a_amrGrids[0].size(); n++)
      {
        pout() << a_amrGrids[0][0] << endl;
      }
    }

    // Now loop over levels, starting with level 1
    int ngrid;
    for (int lev = 1; lev < inNumLevels; lev++)
    {
      is >> ngrid;

      if (a_verbosity >= 3)
      {
        pout() << "level " << lev << " numGrids = " << ngrid << endl;
        pout() << "Grids: ";
      }

      while (is.get() != '\n');

      a_amrGrids[lev].resize(ngrid);

      for (int i = 0; i < ngrid; i++)
      {
        Box bx;
        is >> bx;

        while (is.get() != '\n');

        // Quick check on box size
        Box bxRef(bx);

        if (bxRef.longside() > a_maxGridSize)
        {
          pout() << "Grid " << bx << " too large" << endl;
          MayDay::Error();
        }

        if (a_verbosity >= 3)
        {
          pout() << bx << endl;
        }

        a_amrGrids[lev][i] = bx;
      } // End loop over boxes on this level
    } // End loop over levels
  }

  // Broadcast results to all the processors
  broadcast(a_amrGrids,uniqueProc(SerialTask::compute));
}

// setupVortices allows vortex parameters to be read in and used in this AMR
// computation example
void setupVortices(Vector<RealVect>& a_center,
                   Vector<Real>&     a_radius,
                   Vector<Real>&     a_strength,
                   int               a_verbosity,
                   std::string       a_vortexFile)
{
  // Run this task on one processor
  if (procID() == uniqueProc(SerialTask::compute))
  {
    // Read in predefined grids
    ifstream is(a_vortexFile.c_str(), ios::in);

    if (is.fail())
      {
        MayDay::Error("Cannot open vortex file");
      }

    // Format of file:
    // number of vortices
    // (then for each vortex)
    // x y z of center
    // radius of vortex from x,y,z=0
    // strength

    int inNumVortices;
    // while (is.get() != '\n');
    is >> inNumVortices;

    a_center.resize(inNumVortices);
    a_radius.resize(inNumVortices);
    a_strength.resize(inNumVortices);

    // Now loop over vortices
    for (int iVort = 0; iVort < inNumVortices; iVort++)
      {
        Real x, y, z;
        //        while (is.get() != '\n');
        is >> x;
        is >> y;
        is >> z;
        a_center[iVort] = RealVect(D_DECL6(x, y, z, 0., 0., 0.));

        Real radius;
        //        while (is.get() != '\n');
        is >> radius;
        a_radius[iVort] = radius;

        Real strength;
        //        while (is.get() != '\n');
        is >> strength;
        a_strength[iVort] = strength;
      } // End loop over vortices
  }

  // Broadcast results to all the processors
  D_TERM6(
         broadcast(a_center[0], uniqueProc(SerialTask::compute)); ,
         broadcast(a_center[1], uniqueProc(SerialTask::compute)); ,
         broadcast(a_center[2], uniqueProc(SerialTask::compute)); ,
         broadcast(a_center[3], uniqueProc(SerialTask::compute)); ,
         broadcast(a_center[4], uniqueProc(SerialTask::compute)); ,
         broadcast(a_center[5], uniqueProc(SerialTask::compute)); );
  broadcast(a_radius, uniqueProc(SerialTask::compute));
  broadcast(a_strength, uniqueProc(SerialTask::compute));
}

#ifdef TRAP_FPE
#include <fenv.h>

// FE_INEXACT    inexact result
// FE_DIVBYZERO  division by zero
// FE_UNDERFLOW  result not representable due to underflow
// FE_OVERFLOW   result not representable due to overflow
// FE_INVALID    invalid operation

static void enableFpExceptions ()
{
  if (feclearexcept(FE_ALL_EXCEPT) != 0)
  {
    MayDay::Abort("feclearexcept failed");
  }

  int flags = FE_DIVBYZERO |
              FE_INVALID   |
//              FE_UNDERFLOW |
              FE_OVERFLOW  ;

  if (feenableexcept(flags) == -1)
  {
    MayDay::Abort("feenableexcept failed");
  }
}
#endif
