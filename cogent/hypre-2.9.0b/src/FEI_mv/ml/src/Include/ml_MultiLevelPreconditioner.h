/*BHEADER**********************************************************************
 * Copyright (c) 2008,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of HYPRE.  See file COPYRIGHT for details.
 *
 * HYPRE is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * $Revision: 1.5 $
 ***********************************************************************EHEADER*/




/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

/*!
 * \file ml_MultiLevelPreconditioner.h
 *
 * \class MultiLevelPreconditioner
 *
 * \brief ML black-box preconditioner for Epetra_RowMatrix derived classes.
 *
 * \date Last update do Doxygen: 22-Jul-04
 *
 */

#ifndef ML_MULTILEVELPRECONDITIONER_H
#define ML_MULTILEVELPRECONDITIONER_H

#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) 
// define the following to allow compilation without AztecOO
#ifndef HAVE_ML_AZTECOO
#ifndef AZ_PROC_SIZE
#define AZ_PROC_SIZE 1
#endif
#ifndef AZ_OPTIONS_SIZE
#define AZ_OPTIONS_SIZE 1
#endif
#ifndef AZ_PARAMS_SIZE
#define AZ_PARAMS_SIZE 1
#endif
#ifndef AZ_STATUS_SIZE
#define AZ_STATUS_SIZE 1
#endif
#endif

class Epetra_Map;
class Epetra_BlockMap;
class Epetra_MultiVector;
class Epetra_Comm;
class Epetra_CrsMatrix;
class Epetra_FECrsMatrix;
class Epetra_VbrMatrix;

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"

#define ML_MEM_SIZE      20
#define ML_MEM_INITIAL    0
#define ML_MEM_FINAL      1
#define ML_MEM_SMOOTHER   2
#define ML_MEM_COARSE     3
#define ML_MEM_HIERARCHY  4
#define ML_MEM_PREC_FIRST 5
#define ML_MEM_PREC_OTHER 6
#define ML_MEM_TOT1       7
#define ML_MEM_TOT2       8
#define ML_MEM_INITIAL_MALLOC    10
#define ML_MEM_FINAL_MALLOC      11
#define ML_MEM_SMOOTHER_MALLOC   12
#define ML_MEM_COARSE_MALLOC     13
#define ML_MEM_HIERARCHY_MALLOC  14
#define ML_MEM_PREC_FIRST_MALLOC 15
#define ML_MEM_PREC_OTHER_MALLOC 16
#define ML_MEM_TOT1_MALLOC       17
#define ML_MEM_TOT2_MALLOC       18

#ifdef HAVE_ML_TRIUTILS
#include "Trilinos_Util_CommandLineParser.h"
#endif

#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#ifdef HAVE_ML_AZTECOO
#include "Epetra_MsrMatrix.h"
#endif
#include "Teuchos_ParameterList.hpp"

namespace ML_Epetra
{

  //! Sets default parameters for aggregation-based 2-level domain decomposition preconditioners.
  /*! This function, defined in the namespace ML_Epetra, can be used to set
    default values in a user's defined Teuchos::ParameterList.
    \param ProblemType (In) : a string, whose possible values are:
       - "DD" : defaults for 2-level domain decomposition preconditioners based
       on aggregation;
       - "DD-ML" : 3-level domain decomposition preconditioners, with coarser
       spaces defined by aggregation;
       - "SA" : classical smoothed aggregation preconditioners;
       - "maxwell" : default values for Maxwell.
    \param List (Out) : list which will populated by the default parameters
    \param options (In) : integer array, of size \c AZ_OPTIONS_SIZE, that will be
    populated with suitable values. A pointer to \c options will be stick into
    the parameters list. Note that this array is still required to apply the
    preconditioner! Do not delete options, nor let it go out of scope. The default value is 
    0, meaning that \c SetDefaults() will allocate the array. It is
    responsibility of the user to free this memory.
    \param params (Out) : double array, of size \c AZ_PARAMS_SIZE. See comments
    for \c options.    
   */
  int SetDefaults(string ProblemType, Teuchos::ParameterList & List,
		  int * options = 0, double * params = 0);
  
  //! Sets default parameters for aggregation-based 2-level domain decomposition preconditioners.
  int SetDefaultsDD(Teuchos::ParameterList & List, 
		    int * options = 0, double * params = 0);
  
  //! Sets default parameters for aggregation-based 2-level domain decomposition preconditioners, using LU on each subdomain
  int SetDefaultsDD_LU(Teuchos::ParameterList & List, 
		       int * options = 0, double * params = 0);
  
  //! Sets default parameters for aggregation-based 3-level domain decomposition preconditioners.  
  int SetDefaultsDD_3Levels(Teuchos::ParameterList & List, 
			    int * options = 0, double * params = 0);
  
  //! Sets default parameters for aggregation-based 3-level domain decomposition preconditioners with LU.
  int SetDefaultsDD_3Levels_LU(Teuchos::ParameterList & List, 
			       int * options = 0, double * params = 0);

  //! Sets default parameters for Maxwell's equations.
  int SetDefaultsMaxwell(Teuchos::ParameterList & List, 
			 int * options = 0, double * params = 0);
  
  //! Sets classical smoothed aggregation.
  int SetDefaultsSA(Teuchos::ParameterList & List, 
		    int * options = 0, double * params = 0);

/*!
 
   \brief MultiLevelPreconditioner: a class to define black-box multilevel preconditioners using aggregation methods.

   Class ML_Epetra::MultiLevelPreconditioner defined black-box algebraic
   multilevel preconditioners of matrices defined as Epetra_RowMatrix derived
   objects. The resulting preconditioner can be used in AztecOO, and in any
   other solver that accepts Epetra_Operator derived objects, and apply the
   action of the given Epetra_Operator using ApplyInverse(). 
  
   Please refer to the user's guide for a detailed introduction to
   this class, examples, and description of input parameters.
  
    This file requires ML to be configured with the following options:
    - \c --enable-epetra
    - \c --enable-teuchos
    
    The following option is suggested:
    - \c --enable-amesos
    - \c --enable-ifpack

    Some part of this class needs the following options:
    - \c --enable-aztecoo
    - \c --enable-anasazi
    
    It is important to note that ML is more restrictive than Epetra for
    the definition of maps. It is required that RowMatrixRowMap() is equal 
    to OperatorRangeMap(). This is because ML needs to perform matrix-vector
    product, as well as getrow() functions, on the same data distribution.
    
    Also, for square matrices, OperatorDomainMap() must be as 
    OperatorRangeMap(). 

    Several examples are provided in the \c examples subdirectories:
    - \ref ml_preconditioner_cpp is an introductory 
      example;
    - \ref ml_2level_DD_cpp shows how to
      define a 2-level domain decomposition preconditioner using 
      this class;
    - \ref ml_viz_cpp details how to visualize the aggregates;
    - \ref ml_maxwell_cpp reports how to
      use this class for Maxwell problems.
      
   \note
   Namespace ML_Epetra contains another Epetra_Operator derived class, 
   ML_Epetra::MultiLevelOperator. 
   - you should use MultiLevelOperator
     when your code already defines the required ML objects, with the optimal
     choice of parameters, and you just want to wrap the already defined ML 
     preconditioners for AztecOO problems;
   - you should use MultiLevelPreconditioner
     when you have an Epetra_RowMatrix, and you don't want to code the
     conversion to ML_Operator, the creation of the hierarchy and the
     parameters, simply changing some parameters in a Teuchos::ParameterList.
  
   Defaults parameters can be specified using function SetDefaults().

    \warning The Maxwell interface is still under development. 

    \author Marzio Sala, SNL 9214
*/  
class MultiLevelPreconditioner : public virtual Epetra_Operator {
      
public:  

  //@{ \name Constructors.

  //! Constructs an MultiLevelPreconditioner with default values.

  MultiLevelPreconditioner(const Epetra_RowMatrix & RowMatrix,
                           const bool ComputePrec);

  //! Constructs an MultiLevelPreconditioner. Retrives parameters from \c List.
  
  MultiLevelPreconditioner(const Epetra_RowMatrix & RowMatrix,
			   const Teuchos::ParameterList & List,
			   const bool ComputePrec = true);

  //! Constructs an MultiLevelPreconditioner from an ML_Operator. Retrives parameters from \c List.
  
  MultiLevelPreconditioner(ML_Operator* Operator,
			   const Teuchos::ParameterList& List,
			   const bool ComputePrec = true);
  
  //! Constructs an MultiLevelPreconditioner for Maxwell equations. Retrives parameters from \c List.
  /*! Constructs an MultiLevelPreconditioner for Maxwell equations. The constructor
    requires the edge matrix, the connectivity matrix T, the nodal matrix.
  */
  MultiLevelPreconditioner(const Epetra_RowMatrix& EdgeMatrix,
			   const Epetra_RowMatrix& TMatrix,
			   const Epetra_RowMatrix& NodeMatrix,
			   const Teuchos::ParameterList& List,
			   const bool ComputePrec = true);
  //! Constructs an MultiLevelPreconditioner for Maxwell equations. Retrives parameters from \c List.
  /*! Constructs an MultiLevelPreconditioner for Maxwell equations. The constructor
    requires the edge matrix, the connectivity matrix T, the nodal matrix.
  */
#ifdef HAVE_ML_AZTECOO
MultiLevelPreconditioner(const Epetra_MsrMatrix & EdgeMatrix,
                         ML_Operator * ML_TMatrix,
                         AZ_MATRIX * AZ_NodeMatrix,
                         int       * proc_config,
                         const Teuchos::ParameterList & List,
                         const bool ComputePrec);
#endif

  //@}
  
  //@{ \name Destructor.

  //! Destroys the preconditioner.
  ~MultiLevelPreconditioner() {
    if (IsComputePreconditionerOK_) 
      DestroyPreconditioner(); 
  }

  //@}
  
  //@{ \name Query functions

  //! Prints label associated to this object.
  const char* Label() const{return(Label_);};  
  
  //! Prints unused parameters in the input ParameterList on standard output. */
  void PrintUnused() const
  {
    List_.unused(std::cout);
  }

  //! Prints unused parameters in the input ParameterList on the specified stream.
  void PrintUnused(ostream & os) const
  {
    List_.unused(os);
  }

  //! Prints unused parameters in the input ParameterList to cout on proc \c MyPID. 
  /*! Mispelled parameters are simply ignored. Therefore, it is often the best
   * choice to print out the parameters that have not been used in the
   * construction phase. 
   * - \param MyPID (In) : ID of process that should print the unused parameters.
   */
  void PrintUnused(const int MyPID) const;

  //! Gets a reference to the internally stored parameters' list.
  Teuchos::ParameterList& GetList() 
  {
    return List_;
  }

  // Get a copy of the internally stored output list.
  Teuchos::ParameterList GetOutputList() 
  {
    return OutputList_;
  }

  //! Prints on \c cout the values of the internally stored parameter list for processor \c MyPID
  void PrintList(int MyPID);

  //! Copies \c List into the internally stored parameter list object.
  int SetParameterList(const Teuchos::ParameterList& List);

  //@}
  
  //@{ \name Mathematical functions.

  //! Apply the inverse of the preconditioner to an Epetra_MultiVector (NOT AVAILABLE)
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(-1);}

  //! Apply the preconditioner to an Epetra_MultiVector X, puts the result in Y
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //@}
  
  //@{ \name Atribute access functions


  //! Computes the multilevel hierarchy.
  /*! Computes the multilevel hierarchy. This function retrives the user's defines parameters (as
    specified in the input ParameterList), or takes default values otherwise, and creates the ML
    objects for aggregation and hierarchy. Allocated data can be freed used DestroyPreconditioner(),
    or by the destructor,

    In a Newton-type procedure, several linear systems have to be solved, Often, these systems
    are not too different. In this case, it might be convenient to keep the already 
    computed preconditioner (with hierarchy, coarse solver, smoothers), and use it to
    precondition the next linear system. ML offers a way to determine whether the 
    already available preconditioner is "good enough" for the next linear system. 
    The user should proceed as follows:
    - define \c "adaptive: enable" == \c true
    - solve the first linear system. ML tries to estimate the rate of convergence, and record it;
    - change the values of the linear system matrix (but NOT its structure)
    - compute the new preconditioner as \c ComputePreconditioner(true)
    It is supposed that the pointer to the Epetra_RowMatrix remains constant. Currently,
    it is not possible to modify this pointer (other than creating a new preconditioner)
  */
  
  int ComputePreconditioner(const bool CheckFiltering = false);

  //! Recomputed the preconditioner (not implemented for Maxwell).
  int ReComputePreconditioner();

  //! Print the individual operators in the multigrid hierarchy.
  void Print(const char *whichHierarchy = "main");

  int ComputeAdaptivePreconditioner(int TentativeNullSpaceSize,
                                    double* TentativeNullSpace);

  //! Queries whether multilevel hierarchy has been computed or not.
  int IsPreconditionerComputed()  const
  {
    return(IsComputePreconditionerOK_);
  }

  // following functions are required to derive Epetra_RowMatrix objects.

  //! Sets ownership.
  int SetOwnership(bool ownership){ ownership_ = ownership; return(-1);};

  //! Sets use transpose (not implemented).
  int SetUseTranspose(bool UseTranspose){return(-1);}

  //! Returns the infinity norm (not implemented).
  double NormInf() const {return(0.0);};

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(false);};
  
  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  bool HasNormInf() const{return(false);};
  
  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  const Epetra_Comm& Comm() const{return(*Comm_);};
  
  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map& OperatorDomainMap() const {return(*DomainMap_);};
  
  //! Returns the Epetra_Map object associated with the range of this operator.
  const Epetra_Map& OperatorRangeMap() const {return(*RangeMap_);};
  //@}

  //! Destroys all structures allocated in \c ComputePreconditioner() if the preconditioner has been computed.
  int DestroyPreconditioner();

  //! Returns a reference to the internally stored RowMatrix.
  const Epetra_RowMatrix& RowMatrix() const
  {
    return(*RowMatrix_);
  }

  //! Returns a reference to RowMatrix->Map().
  const Epetra_BlockMap& Map() const
  {
    return(RowMatrix_->Map());
  }

  //! Returns the global number of rows in the matrix.
  int NumGlobalRows() const
  {
    return(RowMatrix_->NumGlobalRows());
  }
  
  //! Returns the global number of columns in the matrix.
  int NumGlobalCols() const
  {
    return(RowMatrix_->NumGlobalCols());
  }
  
  //! Returns the local number of rows in the matrix.
  int NumMyRows() const
  {
    return(RowMatrix_->NumMyRows());
  }
  
  //! Returns the local number of columns in the matrix.
  int NumMyCols() const
  {
    return(RowMatrix_->NumMyCols());
  }
  
  //@}
  //@{ \name debugging and other utilities

  //! Stops the code, waiting for a debugger to attach
  /*! BreakForDebugger() is useful when the user wants to attach to the running
   * process(es). This is a very easy task for serial runs -- just run gdb.
   * Parallel runs may result more problematic. In this case, one can proceed as
   * follows:
   * - define the enviromental variable ML_BREAK_FOR_DEBUGGER (example, in BASH,
   *   \c export \c ML_BREAK_FOR_DEBUGGER=1 )
   * - run the parallel code on a terminal (example, \c mpirun \c -np \c 4 \c
   *   ml_example.exe )
   * - the code will stop in the first call to ComputePreconditioner(). This may
   *   occur in the construction phase. Some information about the ID number of
   *   each process will be shown.
   * - in another terminal, attach to the desired process.
   * - insert one character to let the code continue, and debug as required.
   */
  int BreakForDebugger();

  //! Prints the computational stencil for the specified row and equation (for 2D Cartesian grids only)
  /*! For problems defined on 2D Cartesian grids (with node numbering increasing
   * along the x-axis), this function prints out the stencil in an intelligible
   * form.
   * \param nx (In) : number of nodes along the X-axis
   * \param ny (In) : number of nodes along the Y-axis
   * \param NodeID (In) : (local) ID of node that will be used to print the
   *   stencil. If set to -1, the code will automatically chose an internal node.
   *   Default: -1.
   * \param EquationID (In) : ID of the equation that will be used to print the
   *   stencil (default = 0)  
   */
  int PrintStencil2D(const int nx, const int ny, 
		     int NodeID = -1,
		     const int EquationID = 0);

  //! Cheap analysis of each level matrix.
  int AnalyzeHierarchy(const bool AnalyzeMatrices, 
                       const int PreCycles, const int PostCycles,
                       const int MLCycles);

  //! Analyze the effect of each level's smoother on a random vector.
  int AnalyzeSmoothers(const int NumPreCycles = 1,
                       const int NumPostCycles = 1);

  //! Analyze the effect of the coarse solver on a random vector.
  int AnalyzeCoarse();

  //! Analyze the effect of the ML cycle on a random vector.
  int AnalyzeCycle(const int NumCycles = 1);

  //! Test several smoothers on fine-level matrix.
  int TestSmoothers(Teuchos::ParameterList& InputList,
		    const bool IsSymmetric = false);

  //! Test several smoothers on fine-level matrix using the current parameters.
  int TestSmoothers(const bool IsSymmetric = false) {
    return(TestSmoothers(List_,IsSymmetric));
  }

  //! Returns a pointer to the internally stored ml pointer
  const ML* GetML(const int WhichML= -1) const
  {
    if (WhichML < 0)
      return ml_;
    else if (WhichML == 0)
      return ml_nodes_;
    else
      return(0);
  }

  bool SolvingMaxwell() const
  {
    return SolvingMaxwell_;
  }

  //! Returns a pointer to the internally stored agg pointer
  const ML_Aggregate* GetML_Aggregate() const 
  {
    return agg_;
  }

  //! Generic interface to visualization methods.
  int Visualize(bool VizAggre, bool VizPreSmoother,
		bool VizPostSmoother, bool VizCycle,
		int NumApplPreSmoother, int NumApplPostSmoother,
		int NumCycleSmoother);

  //! Visualizes the shape of the aggregates.
  int VisualizeAggregates();

  //! Visualizes the effect of smoothers on a random vector.
  int VisualizeSmoothers(int NumPrecCycles = 1,
			 int NumPostCycles = 1);

  //! Visualizes the effect of the ML cycle on a random vector.
  int VisualizeCycle(int NumCycles = 1);

//@}

private:

  //! Copy constructor (NOT DEFINED)
  MultiLevelPreconditioner(const MultiLevelPreconditioner & rhs) 
  {};

  //! operator = (NOT DEFINED)
  MultiLevelPreconditioner & operator = (const MultiLevelPreconditioner & rhs)
  {
    return *this;
  };

  //@{ \name Internal setting functions
  //! Initializes object with defauls values.
  int Initialize();

  //! Sets smoothers for non-Maxwell equations.
  int SetSmoothers();

  //! Sets smoothers for Maxwell equations.
  int SetSmoothersMaxwell();

  //! Sets coarse level solvers.
  int SetCoarse();

  //! Sets aggregation schemes.
  int SetAggregation();

  //! Sets preconditioner type (usually, V-cycle).
  int SetPreconditioner();

  //! Sets the null space for non-Maxwell problems.
  int SetNullSpace();

  //! Sets the null space for Maxwell equations.
  int SetNullSpaceMaxwell();

  //! Sets prolongator smoother parameters.
  int SetSmoothingDamping();

  //! Sets damping parameter for classical smoothed aggregation.
  int SetSmoothingDampingClassic();
  
  //! Creates label for this object (printed out by AztecOO)
  int CreateLabel();

#define OLD_AUX
#ifdef OLD_AUX
  int CreateAuxiliaryMatrixCrs(Epetra_FECrsMatrix * & FakeMatrix);

  int CreateAuxiliaryMatrixVbr(Epetra_VbrMatrix * & FakeMatrix);
#endif

  int SetupCoordinates();

  void PrintMem(char *fmt, int size, int, int);

  void PrintMemoryUsage();

  int SetFiltering();

  void RandomAndZero(double *, double *, int);
  
  //! Checks whether the previously computed preconditioner is still valuable for the newly available linear system.
  /*! Used only when \c "adaptive: enable" is \c true, and
   * ComputePreconditioner(true) is called. */
  bool CheckPreconditionerKrylov();

  void VectorNorms(double*, int, double*,double*);

  //@}

  //@{ \name Internal data
  
  //! Pointer to ML_Struct
  ML* ml_;
  //! ML communicator, convenient to have separately from ml_,
  //  ml_nodes_, one or all of which may be null.
  ML_Comm* ml_comm_;
  //! ML_Aggregate, contains aggregate information
  ML_Aggregate* agg_;
  //! Label for this object
  char* Label_;

  //! pointer to linear system matrix
  const Epetra_RowMatrix* RowMatrix_;
  //! specifies whether a hierarchy already exists or not.
  bool IsComputePreconditionerOK_;
  
  //! Number of levels
  int NumLevels_;
  //! Domain Map
  const Epetra_Map* DomainMap_;
  //! Range Map
  const Epetra_Map* RangeMap_;
  //! Epetra communicator object
  const Epetra_Comm* Comm_;
  bool  ownership_;
  //! proc_config for Aztec smoothers
  int   ProcConfig_[AZ_PROC_SIZE];
  //! options for Aztec smoothers
  int   SmootherOptions_[AZ_OPTIONS_SIZE];
  //! params for Aztec smoothers
  double SmootherParams_[AZ_PARAMS_SIZE];
  //! status for Aztec smoothers
  double SmootherStatus_[AZ_STATUS_SIZE];

  //! List containing all input parameters.
  Teuchos::ParameterList List_;
  //! List containing all output parameters
  Teuchos::ParameterList OutputList_;      

  //! Maximum number of levels
  int MaxLevels_;

  //! Number of applications of the ML cycle
  int CycleApplications_;

  //! If \c true, zero starting solution is used in the application of the cycle.
  bool ZeroStartingSolution_;

  //! Integer array used to easily handle ML_INCREASING and ML_DECREASING
  /*! Integer array, of size MaxLevels_, that contain the ML level ID
    for the first logical level, and so on for all levels. The ML level ID
    of logical level L is LevelID_[L].
    In this interface, all levels move from 0 to MaxLevels-1.
    ML's level for interface's level i is LevelID_[i]
  */
  vector<int> LevelID_;

  //! If not NULL, contains the allocated null space vector 
  double* NullSpaceToFree_;              

  //! all cout's have this prefix (default'd in Initialize() )
  string PrintMsg_;
  //! all cerr's have this prefix (default'd in Initialize() )
  char ErrorMsg_[80];
  //! true of information has to be printed on this process
  bool verbose_;
  //! Number of PDE equations.
  int NumPDEEqns_;

  //@}

  //@{ \name Maxwell variables

  //! true if Maxwell equations are used
  bool SolvingMaxwell_;
  //! Main matrix for Maxwell
  const Epetra_RowMatrix* EdgeMatrix_;
  //! aux matrix for Maxwell
  const Epetra_RowMatrix* NodeMatrix_;
  bool CreatedNodeMatrix_;
  //! Auxiliary matrix used in intermediate step
  ML_Operator* ML_Kn_;
  bool CreatedML_Kn_;
  //! T matrix for Maxwell
  const Epetra_RowMatrix* TMatrix_;
  bool CreatedTMatrix_;
  ML_Operator* TMatrixML_;
  ML_Operator* TMatrixTransposeML_;
  ML_Operator** Tmat_array, ** Tmat_trans_array;
  //! ML structures for Maxwell
  ML* ml_nodes_;

  void** nodal_args_,** edge_args_;

  //@}

  //@{ \name Variables for Timing
  //! Number of applications
  int NumApplications_;
  //! CPU time for all applications of the preconditioner
  double ApplicationTime_;
  bool FirstApplication_;
  //@ CPU time for first application
  double FirstApplicationTime_;
  //! Number of construction phases
  int NumConstructions_;
  //! CPU time for construction of the preconditioner.
  double ConstructionTime_;

  //@}
  
  // other stuff for old ML's compatibility
  Epetra_CrsMatrix* RowMatrixAllocated_;

  bool AnalyzeMemory_;
  
  int memory_[ML_MEM_SIZE];

  // filtering stuff
  vector<double> flt_NullSpace_;
  ML* flt_ml_;
  ML_Aggregate* flt_agg_;
  
  // for reuse of preconditioning
  double RateOfConvergence_;

}; // class MultiLevelPreconditioner
 
} // namespace ML_Epetra

#endif /* defined HAVE_ML_EPETRA and HAVE_ML_TEUCHOS */

#endif /* define ML_MULTILEVELPRECONDITIONER_H */
