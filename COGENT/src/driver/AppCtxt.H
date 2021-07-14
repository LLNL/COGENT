#ifndef  _APPCTXT_H_
#define  _APPCTXT_H_

#include "ParmParse.H"
#include "parstream.H"
#include "TiDefinitions.H"
#include "ODEVector.H"
#include "Preconditioner.H"

#include "NamespaceHeader.H"

class AppCtxt 
{
  public:

    AppCtxt() {}
    virtual ~AppCtxt() {}
    
    virtual void initialize( const int, const double ) = 0;
    virtual Real stableDt( const int, const TimeIntegratorType ) = 0;
    
    virtual void writePlotFile(const char*, const int, const double& ) = 0;
    virtual void writeFieldHistory(int, double, bool) = 0;
    virtual void writeCheckpointFile( HDF5Handle&, 
                                      const int, 
                                      const double, 
                                      const double ) = 0;
    virtual void readCheckpointFile(  HDF5Handle&,
                                      int&,
                                      double&, 
                                      double& ) = 0;
    
    virtual void printDiagnostics() = 0;
    
    virtual void preTimeStep  (const ODEVector&, int, Real) = 0;
    virtual void postTimeStep (ODEVector&, int, Real, Real) = 0;
    virtual void postTimeStage( const int,
                                const Real,
                                const Real,
                                ODEVector&,
                                const int ) = 0;
    
    virtual int getVectorSize() = 0;
    
    virtual void copyStateToArray (ODEVector&) const = 0;
    virtual void copyStateFromArray (const ODEVector&) = 0;
    virtual void copyRHSToArray (ODEVector&) const = 0;
    virtual void copyRHSFromArray (const ODEVector&)  = 0;
    
    virtual void preOpEval( const ODEVector&,
                            const Real,
                            const Checkpoint& ) = 0;
    
    virtual void preSolutionOpEval( const ODEVector&,
                                    const Real,
                                    const Checkpoint& ) = 0;
    
    virtual void explicitOp(ODEVector&, 
                            const Real, 
                            const ODEVector& ) = 0;
    
    virtual void explicitOpImEx(ODEVector&,
                                const Real,
                                const ODEVector& ) = 0;
    
    virtual void implicitOpImEx(ODEVector&,
                                const Real,
                                const ODEVector& ) = 0;
    
    virtual void solutionOp(ODEVector&,
                            const Real,
                            const ODEVector& ) = 0;
    
    virtual void explicitPC(ODEVector&,
                            const Real,
                            const ODEVector& ) = 0;
    
    virtual void implicitPCImEx(ODEVector&,
                                const Real,
                                const ODEVector& ) = 0;
    
    virtual void solutionPC(ODEVector&,
                            const Real,
                            const ODEVector& ) = 0;
    
    virtual void 
    defineMultiPhysicsPC( std::vector<Preconditioner<ODEVector,AppCtxt>*>&,
                          std::vector<DOFList>&,
                          const ODEVector&,
                          const std::string&,
                          const std::string&,
                          bool ) = 0;
    
    virtual void 
    updateMultiPhysicsPC( std::vector<Preconditioner<ODEVector,AppCtxt>*>&,
                              const ODEVector&,
                              const Real,
                              const int, 
                              const int,
                              const Real,
                              const bool ) = 0;
    
    virtual void solveSolutionPC( ODEVector&, 
                                  const ODEVector&, 
                                  int a_i = -1 ) = 0;
    
    virtual void solveFluidOpPCImEx(ODEVector&, 
                                    const ODEVector&, 
                                    int a_i = -1 ) = 0;
    
    virtual void 
    solveVlasovOpPCImEx(  Preconditioner<ODEVector,AppCtxt>* const,
                          ODEVector&, 
                          const ODEVector&, 
                          int a_i = -1 ) = 0;
    
    virtual bool isLinear() = 0;
    virtual void printFunctionCounts() = 0;
    
    virtual void defineODEVector( ODEVector& ) const = 0;
    
    virtual void useScales( const bool ) = 0;
    virtual void setScaleTolerance(const Real a_tol) = 0;
    virtual void computeUScale( const ODEVector& ) = 0;
    virtual void computeRHSOpScale( const ODEVector& ) = 0;
    virtual void computeLHSOpScale( const ODEVector& ) = 0;
    virtual void scaleU(ODEVector& )  = 0;
    virtual void scaleRHSOp(ODEVector&) = 0;
    virtual void scaleLHSOp(ODEVector&) = 0;
    virtual void deScaleU(ODEVector&) = 0;
    virtual void deScaleRHSOp(ODEVector&) = 0;
    virtual void deScaleLHSOp(ODEVector&) = 0;
    virtual const ODEVector& uScale() = 0;
    virtual const ODEVector& rhsOpScale() = 0;
    virtual const ODEVector& lhsOpScale() = 0;
    virtual Real scaledVecNormU( const ODEVector& ) = 0;
    virtual Real scaledVecNormRHSOp( const ODEVector& ) = 0;
    virtual Real scaledVecNormLHSOp( const ODEVector& ) = 0;
    
    virtual  bool trivialSolutionOp() const = 0;

    virtual const std::string& tiClass() const = 0;
    virtual const std::string& tiMethod() const = 0;
};

#include "NamespaceFooter.H"

#endif
