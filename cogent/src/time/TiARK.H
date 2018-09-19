#ifndef _TiARK_H_
#define _TiARK_H_

#include <iostream>
#include <string>

#include "ParmParse.H"
#include "parstream.H"
#include "TimeIntegrator.H"
#include "ImplicitStageFunction.H"
#include "ImplicitStageJacobian.H"
#include "ExplicitStageFunction.H"
#include "ExplicitStageJacobian.H"
#include "NewtonSolver.H"

#define _TI_ARK_ "ark"

#include "NamespaceHeader.H"

template <class Vector, class Ops>
class TiARK : public TimeIntegrator<Vector,Ops>
{

  public:

    /// Constructor
    /**
     * Constructor: set m_is_Defined to false.
     */ 
   TiARK<Vector,Ops>() : m_is_Defined(false) {}

    /// Destructor
    /*
     * Clean up allocations
     */
    virtual ~TiARK();

    /// Define the specific ARK method
    /**
     * define the specific ARK method (eg, "1bee", "2a", "3", "4", etc)
     *
     * @param[in] a_name string containing the method name
     */
     virtual void define( ParmParse&, std::string, Vector&, Real );

    /// Advance one time step
    /**
     * Advance one time step.
     *
     * @param[in] a_dt time step
     * @param[in] a_time current simulation time
     * @param[out] a_Y solution
     */
     virtual void advance( const Real& a_time, Vector& a_Y );

    /// Check if method is defined
    /**
     * Returns the value of m_is_Defined
     */
    bool isDefined() const { return m_is_Defined; }
    
    /// Get the operators for the time integrator
    /**
     * get the operators for the time integrator
     */
     virtual Ops& getOperators() { return m_Operators; }

    /// Set the time step size
    /**
     * set the time step size for the time integrator
     *
     * @param[in] a_dt the specified time step size
     */
     virtual void setTimeStepSize( const Real& a_dt ) { m_dt = a_dt; }

    /// Get the time step size
    /**
     * get the time step size of the time integrator
     *
     * @param[out] a_dt the time step size
     */ 
     virtual void getTimeStepSize( Real& a_dt ) const { a_dt = m_dt; }

    /// Set the time step 
    /**
     * set the time step for the time integrator
     *
     * @param[in] a_n the specified time step
     */
     virtual void setTimeStep( const int& a_n ) { m_cur_step = a_n; }

    /// Get the time step
    /**
     * get the time step of the time integrator
     *
     * @param[out] a_n the time step
     */ 
     virtual void getTimeStep( int& a_n ) const { a_n = m_cur_step; }

    /// Set the current simulation time 
    /**
     * set the current simulation time
     *
     * @param[in] a_time the specified simulation time
     */ 
     virtual void setCurrentTime( const Real& a_time ) { m_time = a_time; }

    /// Get the current simulation time
    /*
     * get the current simulation time
     *
     * @param[out] a_time the current simulation time
     */ 
    virtual void getCurrentTime( Real& a_time ) const { a_time = m_time; }

    virtual bool isExplicit() const { return false; }
   
    virtual bool isImEx() const { return true; }

    virtual void printCounts() const
      {
        if (!procID()) {
          cout << "  Time integrator counts:-\n";
          cout << "    Time steps          : " << m_count << "\n";
          if (m_Operators.trivialSolutionOp()) {
            cout << "    Nonlinear iterations: " 
                 << m_ImplicitStageSolver.getCount()
                 << "\n";
            cout << "    Linear iterations   : " 
                 << m_ImplicitStageSolver.getLinearCount()
                 << "\n";
          } else {
            cout << "    Nonlinear iterations: " 
                 << m_ImplicitStageSolver.getCount() + m_ExplicitStageSolver.getCount()
                 << "\n";
            cout << "    Linear iterations   : " 
                 << m_ImplicitStageSolver.getLinearCount() + m_ExplicitStageSolver.getLinearCount()
                 << "\n";
          }
        }
      }

  protected:

  private:
    bool          m_is_Defined;
    bool          m_isLinear;
    bool          m_stagePredictor;
    std::string   m_name;
    int           m_nstages;
    int           m_pinterp;
    Real          *m_Ae, *m_be, *m_ce;
    Real          *m_Ai, *m_bi, *m_ci;
    Real          *m_binterpe, *m_binterpi;
    Vector        m_YStage, m_YPrev,
                  *m_rhsStage_exp, *m_rhsStage_imp, 
                  *m_rhsStage_exp_prev, *m_rhsStage_imp_prev,
                  m_R, m_Z;
    Ops           m_Operators;
    Real          m_time;
    Real          m_dt;
    Real          m_epsJFNK;
    int           m_cur_step,
                  m_count;

    ExplicitStageFunction<Vector,Ops>       *m_Efunction;
    ExplicitStageJacobian<Vector,Ops>       *m_EJacobian;
    ImplicitStageFunction<Vector,Ops>   *m_Ifunction;
    ImplicitStageJacobian<Vector,Ops>   *m_IJacobian;
    NewtonSolver<Vector,Ops>                m_ImplicitStageSolver,
                                            m_ExplicitStageSolver;

    void setCoefficients( int, const Real*, const Real*, const Real*, const Real*,
                          int, const Real*, const Real*);

    void extrapolate( Real, Vector& );
    void parseParameters( ParmParse& );
};


template <class Vector, class Ops>
void TiARK<Vector,Ops>::define(ParmParse& a_pp, std::string a_name, Vector& a_Y, Real a_dt)
{
  m_dt = a_dt;
  m_stagePredictor = false;
  m_epsJFNK = 1e-6;

  /* set Butcher coefficients */
  if (a_name == "1bee") {

    /* 1st order backward Euler represented as an ARK method.
     */

    m_name = a_name;
    m_nstages = 3;

    const Real
      Ae[3][3] = {{0.0,0.0,0.0},
                  {0.0,0.0,0.0},
                  {0.0,0.5,0.0}},
      Ai[3][3] = {{1.0,0.0,0.0},
                  {0.0,0.5,0.0},
                  {0.0,0.5,0.5}};
    setCoefficients(m_nstages,&Ae[0][0],&Ai[0][0],NULL,NULL,1,NULL,NULL);

  } else if (a_name == "2a") {

    /* 2nd order method with A-stable implicit part; 1 explicit, 1 implicit stage. 
     * Ref: Emil Constantinescu 
     */

    m_name = a_name;
    m_nstages = 2;

    const Real
      Ae[2][2] = {{0.0,0.0},
                  {1.0,0.0}},
      Ai[2][2] = {{0.0,0.0},
                  {0.5,0.5}};
    setCoefficients(m_nstages,&Ae[0][0],&Ai[0][0],NULL,NULL,1,NULL,NULL);

  } else if (a_name == "2c") {

    /* 2nd order method with L-stable implicit part; 1 explicit and 2 implicit stages.
     * Ref: Emil Constantinescu
     */

    m_name = a_name;
    m_nstages = 3;

    const Real
      Ae[3][3] = {{0,0,0},
                  {2-1.414213562373095048802,0,0},
                  {0.5,0.5,0}},
      Ai[3][3] = {{0,0,0},
                  {1-1/1.414213562373095048802,1-1/1.414213562373095048802,0},
                  {1/(2*1.414213562373095048802),1/(2*1.414213562373095048802),1-1/1.414213562373095048802}},
      binterpi[3][2] = {{1.0/1.414213562373095048802,-1.0/(2.0*1.414213562373095048802)},
                        {1.0/1.414213562373095048802,-1.0/(2.0*1.414213562373095048802)},
                        {1.0-1.414213562373095048802,1.0/1.414213562373095048802}};
    setCoefficients(m_nstages,&Ae[0][0],&Ai[0][0],NULL,NULL,2,NULL,&binterpi[0][0]);

  } else if (a_name == "2e") {

    /* 2nd order method with L-stable implicit part; 1 explicitand 2 implicit stages.
     * Ref: Emil Constantinescu
     */

    m_name = a_name;
    m_nstages = 3;

    const Real
      Ae[3][3] = {{0,0,0},
                  {2-1.414213562373095048802,0,0},
                  {(3-2*1.414213562373095048802)/6,(3+2*1.414213562373095048802)/6,0}},
      Ai[3][3] = {{0,0,0},
                  {1-1/1.414213562373095048802,1-1/1.414213562373095048802,0},
                  {1/(2*1.414213562373095048802),1/(2*1.414213562373095048802),1-1/1.414213562373095048802}},
      binterpi[3][2] =  {{1.0/1.414213562373095048802,-1.0/(2.0*1.414213562373095048802)},
                         {1.0/1.414213562373095048802,-1.0/(2.0*1.414213562373095048802)},
                         {1.0-1.414213562373095048802,1.0/1.414213562373095048802}};
    setCoefficients(m_nstages,&Ae[0][0],&Ai[0][0],NULL,NULL,2,NULL,&binterpi[0][0]);

  } else if (a_name == "prssp2") {

    /* 2nd order SSP method; 3 implicit stages.
     * Ref: Pareschi & Russo (SSP2-(3,3,2))
     */

    m_name = a_name;
    m_nstages = 3;

    const Real
      Ae[3][3] = {{0,0,0},
                  {0.5,0,0},
                  {0.5,0.5,0}},
      Ai[3][3] = {{0.25,0,0},
                  {0,0.25,0},
                  {1./3,1./3,1./3}};
    setCoefficients(m_nstages,&Ae[0][0],&Ai[0][0],NULL,NULL,0,NULL,NULL);

  } else if (a_name == "3") {

    /* 3rd order method with L-stable implicit part; 1 explicit, 3 implicit stages.
     * Ref: Kennedy & Carpenter
     */

    m_name = a_name;
    m_nstages = 4;

    const Real
      Ae[4][4] = {{0,0,0,0},
                  {1767732205903./2027836641118.,0,0,0},
                  {5535828885825./10492691773637.,788022342437./10882634858940.,0,0},
                  {6485989280629./16251701735622.,-4246266847089./9704473918619.,10755448449292./10357097424841.,0}},
      Ai[4][4] = {{0,0,0,0},
                  {1767732205903./4055673282236.,1767732205903./4055673282236.,0,0},
                  {2746238789719./10658868560708.,-640167445237./6845629431997.,1767732205903./4055673282236.,0},
                  {1471266399579./7840856788654.,-4482444167858./7529755066697.,11266239266428./11593286722821.,1767732205903./4055673282236.}},
      binterpi[4][2] = {{4655552711362./22874653954995., -215264564351./13552729205753.},
                        {-18682724506714./9892148508045.,17870216137069./13817060693119.},
                        {34259539580243./13192909600954.,-28141676662227./17317692491321.},
                        {584795268549./6622622206610.,   2508943948391./7218656332882.}};
    setCoefficients(m_nstages,&Ae[0][0],&Ai[0][0],NULL,NULL,2,NULL,&binterpi[0][0]);

  } else if (a_name == "4") {

    /* Fourth order ARK scheme with L-stable implicit part.
       This method has one explicit stage and four implicit stages.
       References: Kennedy and Carpenter 2003.
    */
    m_name = a_name;
    m_nstages = 6;

    const Real 
      Ae[6][6] = {{0,0,0,0,0,0},
                 {1./2,0,0,0,0,0},
                 {13861./62500.,6889./62500.,0,0,0,0},
                 {-116923316275./2393684061468.,-2731218467317./15368042101831.,9408046702089./11113171139209.,0,0,0},
                 {-451086348788./2902428689909.,-2682348792572./7519795681897.,12662868775082./11960479115383.,3355817975965./11060851509271.,0,0},
                 {647845179188./3216320057751.,73281519250./8382639484533.,552539513391./3454668386233.,3354512671639./8306763924573.,4040./17871.,0}},
      Ai[6][6] = {{0,0,0,0,0,0},
                  {1./4,1./4,0,0,0,0},
                  {8611./62500.,-1743./31250.,1./4,0,0,0},
                  {5012029./34652500.,-654441./2922500.,174375./388108.,1./4,0,0},
                  {15267082809./155376265600.,-71443401./120774400.,730878875./902184768.,2285395./8070912.,1./4,0},
                  {82889./524892.,0,15625./83664.,69875./102672.,-2260./8211,1./4}},
      binterpi[6][3] = {{6943876665148./7220017795957.,-54480133./30881146.,6818779379841./7100303317025.},
                        {0,0,0},
                        {7640104374378./9702883013639.,-11436875./14766696.,2173542590792./12501825683035.},
                        {-20649996744609./7521556579894.,174696575./18121608.,-31592104683404./5083833661969.},
                        {8854892464581./2390941311638.,-12120380./966161.,61146701046299./7138195549469.},
                        {-11397109935349./6675773540249.,3843./706.,-17219254887155./4939391667607.}};
    setCoefficients(m_nstages,&Ae[0][0],&Ai[0][0],NULL,NULL,3,NULL,&binterpi[0][0]);

  } else if (a_name == "5") {

    /* 5th order method with L-stable implicit part. 1 explicit, 5 implicit stages.
     * Ref: Kennedy & Carpenter
     */

    m_name = a_name;
    m_nstages = 8;

    const Real
      Ae[8][8] = {{0,0,0,0,0,0,0,0},
                  {41./100,0,0,0,0,0,0,0},
                  {367902744464./2072280473677.,677623207551./8224143866563.,0,0,0,0,0,0},
                  {1268023523408./10340822734521.,0,1029933939417./13636558850479.,0,0,0,0,0},
                  {14463281900351./6315353703477.,0,66114435211212./5879490589093.,-54053170152839./4284798021562.,0,0,0,0},
                  {14090043504691./34967701212078.,0,15191511035443./11219624916014.,-18461159152457./12425892160975.,-281667163811./9011619295870.,0,0,0},
                  {19230459214898./13134317526959.,0,21275331358303./2942455364971.,-38145345988419./4862620318723.,-1./8,-1./8,0,0},
                  {-19977161125411./11928030595625.,0,-40795976796054./6384907823539.,177454434618887./12078138498510.,782672205425./8267701900261.,-69563011059811./9646580694205.,7356628210526./4942186776405.,0}},
      Ai[8][8] = {{0,0,0,0,0,0,0,0},
                  {41./200.,41./200.,0,0,0,0,0,0},
                  {41./400.,-567603406766./11931857230679.,41./200.,0,0,0,0,0},
                  {683785636431./9252920307686.,0,-110385047103./1367015193373.,41./200.,0,0,0,0},
                  {3016520224154./10081342136671.,0,30586259806659./12414158314087.,-22760509404356./11113319521817.,41./200.,0,0,0},
                  {218866479029./1489978393911.,0,638256894668./5436446318841.,-1179710474555./5321154724896.,-60928119172./8023461067671.,41./200.,0,0},
                  {1020004230633./5715676835656.,0,25762820946817./25263940353407.,-2161375909145./9755907335909.,-211217309593./5846859502534.,-4269925059573./7827059040749.,41./200,0},
                  {-872700587467./9133579230613.,0,0,22348218063261./9555858737531.,-1143369518992./8141816002931.,-39379526789629./19018526304540.,32727382324388./42900044865799.,41./200.}},
      binterpi[8][3] = {{-17674230611817./10670229744614.,  43486358583215./12773830924787., -9257016797708./5021505065439.},
                        {0,  0, 0                            },
                        {0,  0, 0                            },
                        {65168852399939./7868540260826.,  -91478233927265./11067650958493., 26096422576131./11239449250142.},
                        {15494834004392./5936557850923.,  -79368583304911./10890268929626., 92396832856987./20362823103730.},
                        {-99329723586156./26959484932159.,  -12239297817655./9152339842473., 30029262896817./10175596800299.},
                        {-19024464361622./5461577185407.,  115839755401235./10719374521269., -26136350496073./3983972220547.},
                        {-6511271360970./6095937251113.,  5843115559534./2180450260947., -5289405421727./3760307252460. }};
    setCoefficients(m_nstages,&Ae[0][0],&Ai[0][0],NULL,NULL,3,NULL,&binterpi[0][0]);

  } else if (a_name == "rk4") {

    /* explicit RK4, disguised as an ARK method -- for debugging */
    m_name = a_name;
    m_nstages = 4;

    const Real
      A[4][4] = {{0,0,0,0},
                 {0.5,0,0,0},
                 {0,0.5,0,0},
                 {0,0,1.0,0}},
      b[4]    = {1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0};
    setCoefficients(m_nstages,&A[0][0],&A[0][0],&b[0],&b[0],0,NULL,NULL);

  } else {

    /* default: ARK4 */
    if (!procID()) cout << "Error: unknown ARK method specified " << a_name << ". Using default.\n";
    m_name = "4";
    m_nstages = 6;

    const Real 
      Ae[6][6] = {{0,0,0,0,0,0},
                 {1./2,0,0,0,0,0},
                 {13861./62500.,6889./62500.,0,0,0,0},
                 {-116923316275./2393684061468.,-2731218467317./15368042101831.,9408046702089./11113171139209.,0,0,0},
                 {-451086348788./2902428689909.,-2682348792572./7519795681897.,12662868775082./11960479115383.,3355817975965./11060851509271.,0,0},
                 {647845179188./3216320057751.,73281519250./8382639484533.,552539513391./3454668386233.,3354512671639./8306763924573.,4040./17871.,0}},
      Ai[6][6] = {{0,0,0,0,0,0},
                  {1./4,1./4,0,0,0,0},
                  {8611./62500.,-1743./31250.,1./4,0,0,0},
                  {5012029./34652500.,-654441./2922500.,174375./388108.,1./4,0,0},
                  {15267082809./155376265600.,-71443401./120774400.,730878875./902184768.,2285395./8070912.,1./4,0},
                  {82889./524892.,0,15625./83664.,69875./102672.,-2260./8211,1./4}},
      binterpi[6][3] = {{6943876665148./7220017795957.,-54480133./30881146.,6818779379841./7100303317025.},
                        {0,0,0},
                        {7640104374378./9702883013639.,-11436875./14766696.,2173542590792./12501825683035.},
                        {-20649996744609./7521556579894.,174696575./18121608.,-31592104683404./5083833661969.},
                        {8854892464581./2390941311638.,-12120380./966161.,61146701046299./7138195549469.},
                        {-11397109935349./6675773540249.,3843./706.,-17219254887155./4939391667607.}};
    setCoefficients(m_nstages,&Ae[0][0],&Ai[0][0],NULL,NULL,3,NULL,&binterpi[0][0]);

  }

  ParmParse ppARK("ark");
  parseParameters( ppARK );
  if ((m_pinterp == 0) && (m_stagePredictor)) {
    if (!procID()) {
      cout << "ARK method \"" << m_name << "\" does not have dense output. ";
      cout << "Cannot use stage predictors. Set ark.stage_predictor to false ";
      cout << "or remove it from input file (default: false).\n";
    }
    return;
  }

  /* allocate RHS */
  m_rhsStage_exp  = new Vector[m_nstages];
  m_rhsStage_imp  = new Vector[m_nstages];

  /* define the work vectors and operators */
  m_YStage.define(a_Y);
  m_R.define(a_Y);
  m_Z.define(a_Y);
  for (int i=0; i<m_nstages; i++) {
    m_rhsStage_exp[i].define(a_Y);
    m_rhsStage_imp[i].define(a_Y);
  }
  m_Operators.define((*a_Y.getState()), m_dt);
  m_isLinear = m_Operators.isLinear();

  /* define vectors for stage predictors */
  if (m_stagePredictor) {
    m_YPrev.define(a_Y);
    m_rhsStage_exp_prev  = new Vector[m_nstages];
    m_rhsStage_imp_prev  = new Vector[m_nstages];
    for (int i=0; i<m_nstages; i++) {
      m_rhsStage_exp_prev[i].define(a_Y);
      m_rhsStage_imp_prev[i].define(a_Y);
    }
  } else m_rhsStage_exp_prev = m_rhsStage_imp_prev = NULL;

  /* define string prefixes for options parsing */
  std::string isolverOutPrefix = "TiARK.i";
  std::string esolverOutPrefix = "TiARK.e";
  std::string solverOptPrefix = "ark";

  /* define the Newton solvers */
  Function<Vector,Ops>       
    **ImplicitStageSolverFunction = m_ImplicitStageSolver.getFunction();
  LinearFunction<Vector,Ops> 
    **ImplicitStageSolverJacobian = m_ImplicitStageSolver.getJacobian();

  *ImplicitStageSolverFunction = new ImplicitStageFunction<Vector,Ops>;
  *ImplicitStageSolverJacobian = new ImplicitStageJacobian<Vector,Ops>;
  
  m_Ifunction = dynamic_cast<ImplicitStageFunction<Vector,Ops>*>
                  (*ImplicitStageSolverFunction);
  m_IJacobian = dynamic_cast<ImplicitStageJacobian<Vector,Ops>*>
                  (*ImplicitStageSolverJacobian);

  m_Ifunction->define(m_R,m_Operators);
  m_IJacobian->define(m_R,m_Operators,isolverOutPrefix,solverOptPrefix);
  m_Ifunction->setJac(m_IJacobian);
  m_IJacobian->setIsLinear(m_isLinear);
  m_IJacobian->setJFNKEps(m_epsJFNK);

  m_ImplicitStageSolver.define(a_pp,m_R,isolverOutPrefix,solverOptPrefix);
  m_ImplicitStageSolver.setIsLinear(m_isLinear);
  
  if (!m_Operators.trivialSolutionOp()) {
    Function<Vector,Ops>       
      **ExplicitStageSolverFunction = m_ExplicitStageSolver.getFunction();
    LinearFunction<Vector,Ops> 
      **ExplicitStageSolverJacobian = m_ExplicitStageSolver.getJacobian();
  
    *ExplicitStageSolverFunction = new ExplicitStageFunction<Vector,Ops>;
    *ExplicitStageSolverJacobian = new ExplicitStageJacobian<Vector,Ops>;
    
    m_Efunction = dynamic_cast<ExplicitStageFunction<Vector,Ops>*>
                    (*ExplicitStageSolverFunction);
    m_EJacobian = dynamic_cast<ExplicitStageJacobian<Vector,Ops>*>
                    (*ExplicitStageSolverJacobian);
  
    m_Efunction->define(m_R,m_Operators);
    m_EJacobian->define(m_R,m_Operators,esolverOutPrefix,solverOptPrefix);
    m_Efunction->setJac(m_EJacobian);
    m_EJacobian->setIsLinear(m_isLinear);
    m_EJacobian->setJFNKEps(m_epsJFNK);
  
    m_ExplicitStageSolver.define(a_pp,m_R,esolverOutPrefix,solverOptPrefix);
    m_ExplicitStageSolver.setIsLinear(m_isLinear);
  }

  m_count = 0;

  /* done */
  m_is_Defined = true;
  if (!procID()) cout << "Time integration method: ark (" << m_name << ")\n" ;
}

template <class Vector, class Ops>
void TiARK<Vector, Ops>::setCoefficients( int a_nstages,
                                          const Real* a_Ae, const Real* a_Ai, 
                                          const Real* a_be, const Real* a_bi,
                                          int a_pinterp,
                                          const Real* a_binterpe, const Real* a_binterpi
                                        )
{
  CH_assert(!isDefined());
  CH_assert(a_nstages == m_nstages);

  /* allocate Butcher tableaux coefficients 
   * deallocated in destructor */
  m_Ae = new Real[m_nstages*m_nstages];
  m_be = new Real[m_nstages];
  m_ce = new Real[m_nstages];
  m_Ai = new Real[m_nstages*m_nstages];
  m_bi = new Real[m_nstages];
  m_ci = new Real[m_nstages];

  int i, j;
  for (i=0; i<m_nstages*m_nstages; i++) {
    m_Ae[i] = a_Ae[i];
    m_Ai[i] = a_Ai[i];
  }
  if (a_bi) for (i=0; i<m_nstages; i++) m_bi[i] = a_bi[i];
  else      for (i=0; i<m_nstages; i++) m_bi[i] = a_Ai[(m_nstages-1)*m_nstages+i];
  if (a_be) for (i=0; i<m_nstages; i++) m_be[i] = a_be[i];
  else      for (i=0; i<m_nstages; i++) m_be[i] = m_bi[i];

  for (i=0; i<m_nstages; i++) {
    m_ce[i] = m_ci[i] = 0.0;
    for(j=0; j<m_nstages; j++) { 
      m_ce[i] += a_Ae[i*m_nstages+j]; 
      m_ci[i] += a_Ai[i*m_nstages+j]; 
    }
  }

  if (a_pinterp != 0) {
    m_pinterp = a_pinterp;
    m_binterpe = new Real[m_nstages*m_pinterp];
    m_binterpi = new Real[m_nstages*m_pinterp];
    for (i=0; i<m_nstages; i++) {
      for (j=0; j<m_pinterp; j++) {
        m_binterpi[i*m_pinterp+j] = (a_binterpi?a_binterpi[i*m_pinterp+j]:m_be[i]);
        m_binterpe[i*m_pinterp+j] = (a_binterpe?a_binterpe[i*m_pinterp+j]:m_binterpi[i*m_pinterp+j]);
      }
    }
  } else m_binterpe = m_binterpi = NULL;
}

template <class Vector, class Ops>
TiARK<Vector, Ops>::~TiARK()
{
  delete[] m_Ae;
  delete[] m_be;
  delete[] m_ce;
  delete[] m_Ai;
  delete[] m_bi;
  delete[] m_ci;
  if (m_pinterp != 0) {
    delete[] m_binterpe;
    delete[] m_binterpi;
  }
  delete[] m_rhsStage_exp;
  delete[] m_rhsStage_imp;
  if (m_stagePredictor) {
    delete[] m_rhsStage_exp_prev;
    delete[] m_rhsStage_imp_prev;
  }
  delete m_Ifunction;
  delete m_IJacobian;
  if (!m_Operators.trivialSolutionOp()) {
    delete m_Efunction;
    delete m_EJacobian;
  }
}

template <class Vector, class Ops>
void TiARK<Vector, Ops>::advance( const Real& a_time, Vector& a_Y )
{
  //CH_TIMERS("TiARK::advance");
  CH_assert(isDefined());
  CH_assert(m_time == a_time);
  m_Z.copy(a_Y);

  /* Assert explicit-first-stage - post-time-stage function implementation 
   * in GKOps.cpp assumes that YStage(1) = y_n!
  */
  CH_assert(m_Ai[0] == 0.0);

  /* Stage calculations */
  int i, j;
  for (i = 0; i < m_nstages; i++) {

    Real stage_time = m_time+m_ce[i]*m_dt;

    if (m_Operators.trivialSolutionOp()) {
       m_YStage.copy(a_Y);
    }
    else {
       m_Operators.preOpEval(a_Y,m_time);
       m_Operators.solutionOp(m_YStage,m_time,a_Y);
    }

    for (j=0; j<i; j++) m_YStage.increment(m_rhsStage_exp[j],(m_dt*m_Ae[i*m_nstages+j]));
    for (j=0; j<i; j++) m_YStage.increment(m_rhsStage_imp[j],(m_dt*m_Ai[i*m_nstages+j]));

    if (m_Ai[i*m_nstages+i] != 0.0) {
      
      /* implicit stage */
      Real shift = 1.0 / (m_Ai[i*m_nstages+i] * m_dt);

      m_Ifunction->setShift(shift);
      m_Ifunction->setStageTime(stage_time);

      m_IJacobian->setShift(shift);
      m_IJacobian->setStageTime(stage_time);

      /* Build precond matrix from m_Z, which is
       * first stage - a_Y
       * subsequent stages - m_YStage from previous stage */
      m_IJacobian->updatePreCondMat(m_Z);

      /* right-hand side */
      m_R.copy(m_YStage);  
      m_R.scale(shift);

      /* Before the Newton solve begins, m_Z contains the
       * initial guess. If using stage predictor, initial 
       * guess is computed by extrapolating from previous step */
      if (m_stagePredictor && (m_cur_step > 0)) extrapolate(m_ce[i],m_Z);

      m_Z = m_YStage;

      m_ImplicitStageSolver.solve(m_Z,m_R);
      m_YStage.copy(m_Z);

    } else if (!m_Operators.trivialSolutionOp()) {

      m_Efunction->setStageTime(stage_time);

      m_EJacobian->setStageTime(stage_time);

      /* Build precond matrix from m_Z, which is
       * first stage - a_Y
       * subsequent stages - m_YStage from previous stage */
      m_EJacobian->updatePreCondMat(m_Z);

      /* right-hand side */
      m_R.copy(m_YStage);  

      m_ExplicitStageSolver.solve(m_Z, m_R);
      m_YStage.copy(m_Z);

    }

    m_Operators.postTimeStage(m_cur_step,stage_time,m_YStage,i);
    m_Operators.preOpEval(m_YStage,stage_time);
    m_Operators.explicitOpImEx(m_rhsStage_exp[i],stage_time,m_YStage);
    m_Operators.implicitOpImEx(m_rhsStage_imp[i],stage_time,m_YStage);
  }

  /* save stuff for stage predictor */
  if (m_stagePredictor) {
    m_YPrev.copy(a_Y);
    for (i = 0; i< m_nstages; i++) {
      m_rhsStage_exp_prev[i].copy(m_rhsStage_exp[i]);
      m_rhsStage_imp_prev[i].copy(m_rhsStage_imp[i]);
    }
  }

  /* Step completion */
  if (m_Operators.trivialSolutionOp()) {

    for (i = 0; i < m_nstages; i++) {
      a_Y.increment(m_rhsStage_exp[i],(m_dt*m_be[i]));
      a_Y.increment(m_rhsStage_imp[i],(m_dt*m_bi[i]));
    }

  } else {

    m_Efunction->setStageTime(m_time+m_dt);

    m_EJacobian->setStageTime(m_time+m_dt);

    m_EJacobian->updatePreCondMat(a_Y);

    m_Operators.preOpEval(a_Y,(m_time+m_dt));
    m_Operators.solutionOp(m_R,(m_time+m_dt),a_Y);
    for (i = 0; i < m_nstages; i++) {
      m_R.increment(m_rhsStage_exp[i],(m_dt*m_be[i]));
      m_R.increment(m_rhsStage_imp[i],(m_dt*m_bi[i]));
    }

    m_ExplicitStageSolver.solve(a_Y,m_R);
  }

  /* update current time and step number */
  m_cur_step++;
  m_time += m_dt; 
  m_count++;
}

/* Use dense output to extrapolate from previous stage */
template <class Vector, class Ops>
void TiARK<Vector, Ops>::extrapolate(Real c, Vector& a_Y)
{
  Real h, tt, t, *bt,*b;
  int i,j;
  if (!m_binterpi || !m_binterpe) {
    if (!procID()) cout << "Error: ARK method \"" << m_name << "\" does not have an interpolation formula.\n";
    return;
  }

  t = 1.0 + c;
  h = m_dt;
  bt = (Real*) calloc (m_nstages,sizeof(Real));
  b  = (Real*) calloc (m_nstages,sizeof(Real));

  for (i=0; i<m_nstages; i++) bt[i] = b[i] = 0;
  for (j=0,tt=t; j<m_pinterp; j++,tt*=t) {
    for (i=0; i<m_nstages; i++) {
      bt[i] += h * m_binterpi[i*m_pinterp+j] * tt;
      b[i]  += h * m_binterpe[i*m_pinterp+j] * tt;
    }
  }

  a_Y.copy(m_YPrev);
  for (i=0; i<m_nstages; i++) {
    a_Y.increment(m_rhsStage_imp_prev[i],bt[i]);
    a_Y.increment(m_rhsStage_exp_prev[i],b[i] );
  }

  free(bt);
  free(b);
}

template <class Vector, class Ops>
void TiARK<Vector,Ops>::parseParameters( ParmParse& a_pp)
{
  /* use stage predictor based on dense output */
  a_pp.query("stage_predictor", m_stagePredictor); 
  a_pp.query("jfnk_epsilon", m_epsJFNK); 
}

#include "NamespaceFooter.H"

#endif
