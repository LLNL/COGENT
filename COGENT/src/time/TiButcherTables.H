#ifndef _TIBUTCHER_H_
#define _TIBUTCHER_H_

#include <string.h>

#include "NamespaceHeader.H"

namespace TiButcherTables{

  void setCoefficients( Real**            a_A,
                        Real**            a_b,
                        Real**            a_c,
                        Real**            a_bemb,
                        const int         a_nstages,
                        const Real* const A,
                        const Real* const b,
                        const Real* const bemb = NULL )
  {
    (*a_A) = new Real[a_nstages*a_nstages];
    (*a_b) = new Real[a_nstages];
    (*a_c) = new Real[a_nstages];
    if (bemb != NULL) {
      (*a_bemb) = new Real[a_nstages];
    }
  
    int i, j;
    for (i=0; i<a_nstages*a_nstages; i++) (*a_A)[i] = A[i];
    for (i=0; i<a_nstages;           i++) (*a_b)[i] = b[i];
    for (i=0; i<a_nstages; i++) {
      (*a_c)[i] = 0.0; 
      for(j=0; j<a_nstages; j++) (*a_c)[i] += (*a_A)[i*a_nstages+j]; 
    }
    if (bemb != NULL) {
      (*a_bemb) = new Real[a_nstages];
      for (i=0; i<a_nstages; i++) (*a_bemb)[i] = bemb[i];
    }
  }

  void getButcherTableRK( const std::string&  a_name,
                          int&                a_nstages,
                          int&                a_order,
                          int&                a_order_emb,
                          Real**              a_A,
                          Real**              a_b,
                          Real**              a_c,
                          Real**              a_bemb )
  {
    if (a_name == "1fe") {
  
      /* forward Euler */
      a_nstages = 1;
      a_order = 1;
      a_order_emb = 0;
  
      const Real
        A[1][1] = {{0.0}},
        b[1]    = {1.0};
      setCoefficients(a_A, a_b, a_c, a_bemb, a_nstages, &A[0][0], &b[0]);
  
    } else if (a_name == "2") {
  
      /* 2nd order, 2-stage Runge-Kutta */
      a_nstages = 2;
      a_order = 2;
      a_order_emb = 0;
      
      const Real
        A[2][2]     = {{0.0,0.0},
                      {0.5,0.0}},
        b[2]        = {0.0,1.0};
      setCoefficients(a_A, a_b, a_c, a_bemb, a_nstages, &A[0][0], &b[0]); 
  
    } else if (a_name == "2a") {
  
      /* 2nd order, 2-stage Runge-Kutta */
      a_nstages = 2;
      a_order = 2;
      a_order_emb = 0;
      
      const Real
        A[2][2]     = {{0.0,0.0},
                      {1.0,0.0}},
        b[2]        = {0.5,0.5};
      setCoefficients(a_A, a_b, a_c, a_bemb, a_nstages, &A[0][0], &b[0]); 
  
    } else if (a_name == "3") {
  
      /* 3rd order, 3-stage Runge-Kutta */
      a_nstages = 3;
      a_order = 3;
      a_order_emb = 0;
      
      const Real
        A[3][3] = {{0,0,0},
                   {2.0/3.0,0,0},
                   {-1.0/3.0,1.0,0}},
        b[3]    = {0.25,0.5,0.25};
      setCoefficients(a_A, a_b, a_c, a_bemb, a_nstages, &A[0][0], &b[0]); 
  
    } else if (a_name == "4") {
      
      /* 4th order, 4-stage Runge-Kutta */
      a_nstages = 4;
      a_order = 4;
      a_order_emb = 0;
      
      const Real
        A[4][4] = {{0,0,0,0},
                   {0.5,0,0,0},
                   {0,0.5,0,0},
                   {0,0,1.0,0}},
        b[4]    = {1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0};
      setCoefficients(a_A, a_b, a_c, a_bemb, a_nstages, &A[0][0], &b[0]); 
  
    } else if (a_name == "5dp") {
  
      /* 5th order, 7-stage Dormand-Prince */
      a_nstages = 7;
      a_order = 5;
      a_order_emb = 0;
  
      const Real
        A[7][7]   = {{0,0,0,0,0,0,0},
                     {1.0/5.0,0,0,0,0,0,0},
                     {3.0/40.0,9.0/40.0,0,0,0,0,0},
                     {44.0/45.0,-56.0/15.0,32.0/9.0,0,0,0,0},
                     {19372.0/6561.0,-25360.0/2187.0,64448.0/6561.0,-212.0/729.0,0,0,0},
                     {9017.0/3168.0,-355.0/33.0,46732.0/5247.0,49.0/176.0,-5103.0/18656.0,0,0},
                     {35.0/384.0,0,500.0/1113.0,125.0/192.0,-2187.0/6784.0,11.0/84.0,0}},
        b[7]      = {35.0/384.0,0,500.0/1113.0,125.0/192.0,-2187.0/6784.0,11.0/84.0,0};
      setCoefficients(a_A, a_b, a_c, a_bemb, a_nstages, &A[0][0], &b[0]); 
  
    } else {
      
      /* default: RK4 */
      if (!procID()) {
        cout  << "Warning: unknown RK method specified " 
              << a_name << ". Using default.\n";
      }

      a_nstages = 4;
      a_order = 4;
      a_order_emb = 0;
      
      const Real
        A[4][4] = {{0,0,0,0},
                   {0.5,0,0,0},
                   {0,0.5,0,0},
                   {0,0,1.0,0}},
        b[4]    = {1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0};
      setCoefficients(a_A, a_b, a_c, a_bemb, a_nstages, &A[0][0], &b[0]); 
  
    }

    return;
  }

  void setCoefficients( Real**            a_Ai,
                        Real**            a_Ae,
                        Real**            a_bi,
                        Real**            a_be,
                        Real**            a_c,
                        Real**            a_bembi,
                        Real**            a_bembe,
                        Real**            a_binterpi,
                        Real**            a_binterpe,
                        const int         a_nstages,
                        const Real* const Ae,
                        const Real* const Ai,
                        const Real* const be,
                        const Real* const bi,
                        const int         pinterp,
                        const Real* const binterpe,
                        const Real* const binterpi,
                        const Real* const bembe = NULL,
                        const Real* const bembi = NULL)
  {
    (*a_Ae) = new Real[a_nstages*a_nstages];
    (*a_be) = new Real[a_nstages];
    (*a_c) = new Real[a_nstages];
    (*a_Ai) = new Real[a_nstages*a_nstages];
    (*a_bi) = new Real[a_nstages];
  
    int i, j;
    for (i=0; i<a_nstages*a_nstages; i++) {
      (*a_Ae)[i] = Ae[i];
      (*a_Ai)[i] = Ai[i];
    }
    if (bi) for (i=0; i<a_nstages; i++) (*a_bi)[i] = bi[i];
    else    for (i=0; i<a_nstages; i++) (*a_bi)[i] = Ai[(a_nstages-1)*a_nstages+i];
    if (be) for (i=0; i<a_nstages; i++) (*a_be)[i] = be[i];
    else    for (i=0; i<a_nstages; i++) (*a_be)[i] = (*a_bi)[i];
  
    for (i=0; i<a_nstages; i++) {
      (*a_c)[i] = 0.0;
      for(j=0; j<a_nstages; j++) { 
        (*a_c)[i] += (*a_Ai)[i*a_nstages+j]; 
      }
    }
  
    if (pinterp != 0) {
      (*a_binterpe) = new Real[a_nstages*pinterp];
      (*a_binterpi) = new Real[a_nstages*pinterp];
      for (i=0; i<a_nstages; i++) {
        for (j=0; j<pinterp; j++) {
          (*a_binterpi)[i*pinterp+j] = (binterpi ? binterpi[i*pinterp+j] : (*a_be)[i]);
          (*a_binterpe)[i*pinterp+j] = (binterpe ? binterpe[i*pinterp+j] : (*a_binterpi)[i*pinterp+j]);
        }
      }
    } 

    if (bembe != NULL) {
      (*a_bembe) = new Real[a_nstages];
      for (i=0; i<a_nstages; i++) (*a_bembe)[i] = bembe[i];
    }
    if (bembi != NULL) {
      (*a_bembi) = new Real[a_nstages];
      for (i=0; i<a_nstages; i++) (*a_bembi)[i] = bembi[i];
    }

    return;
  }

  void getButcherTableARK(  const std::string&  a_name,
                            int&                a_nstages,
                            int&                a_order,
                            int&                a_order_emb,
                            Real**              a_Ai,
                            Real**              a_Ae,
                            Real**              a_bi,
                            Real**              a_be,
                            Real**              a_c,
                            Real**              a_bi_emb,
                            Real**              a_be_emb,
                            Real**              a_binterpi,
                            Real**              a_binterpe )
  {
    if (a_name == "1bee") {
  
      /* 1st order backward Euler represented as an ARK method. */
  
      a_nstages = 3;
      a_order = 1;
      a_order_emb = 0;

      const Real
        Ae[3][3] = {{0.0,0.0,0.0},
                    {0.0,0.0,0.0},
                    {0.0,0.5,0.0}},
        Ai[3][3] = {{1.0,0.0,0.0},
                    {0.0,0.5,0.0},
                    {0.0,0.5,0.5}};

      setCoefficients(  a_Ai, a_Ae,
                        a_bi, a_be,
                        a_c,
                        a_bi_emb, a_be_emb,
                        a_binterpi, a_binterpe,
                        a_nstages,
                        &Ae[0][0], &Ai[0][0],
                        NULL, NULL,
                        1,
                        NULL, NULL );
  
    } else if (a_name == "2a") {
  
      /* 2nd order method with A-stable implicit part; 
       * 1 explicit, 1 implicit stage. 
       * Ref: Emil Constantinescu 
       */
  
      a_nstages = 2;
      a_order = 2;
      a_order_emb = 0;
  
      const Real
        Ae[2][2] = {{0.0,0.0},
                    {1.0,0.0}},
        Ai[2][2] = {{0.0,0.0},
                    {0.5,0.5}};

      setCoefficients(  a_Ai, a_Ae,
                        a_bi, a_be,
                        a_c,
                        a_bi_emb, a_be_emb,
                        a_binterpi, a_binterpe,
                        a_nstages,
                        &Ae[0][0], &Ai[0][0],
                        NULL, NULL,
                        1,
                        NULL, NULL );
  
    } else if (a_name == "2c") {
  
      /* 2nd order method with L-stable implicit part; 
       * 1 explicit and 2 implicit stages.
       * Ref: Emil Constantinescu
       */
  
      a_nstages = 3;
      a_order = 2;
      a_order_emb = 0;
  
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
  
      setCoefficients(  a_Ai, a_Ae,
                        a_bi, a_be,
                        a_c,
                        a_bi_emb, a_be_emb,
                        a_binterpi, a_binterpe,
                        a_nstages,
                        &Ae[0][0], &Ai[0][0],
                        NULL, NULL,
                        1,
                        NULL, &binterpi[0][0] );
  
    } else if (a_name == "2e") {
  
      /* 2nd order method with L-stable implicit part; 
       * 1 explicitand 2 implicit stages.
       * Ref: Emil Constantinescu
       */
  
      a_nstages = 3;
      a_order = 2;
      a_order_emb = 0;
  
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
  
      setCoefficients(  a_Ai, a_Ae,
                        a_bi, a_be,
                        a_c,
                        a_bi_emb, a_be_emb,
                        a_binterpi, a_binterpe,
                        a_nstages,
                        &Ae[0][0], &Ai[0][0],
                        NULL, NULL,
                        2,
                        NULL, &binterpi[0][0] );
  
    } else if (a_name == "3") {
  
      /* 3rd order method with L-stable implicit part; 
       * 1 explicit, 3 implicit stages.
       * Ref: Kennedy & Carpenter
       */
  
      a_nstages = 4;
      a_order = 3;
      a_order_emb = 0;
  
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
  
      setCoefficients(  a_Ai, a_Ae,
                        a_bi, a_be,
                        a_c,
                        a_bi_emb, a_be_emb,
                        a_binterpi, a_binterpe,
                        a_nstages,
                        &Ae[0][0], &Ai[0][0],
                        NULL, NULL,
                        2,
                        NULL, &binterpi[0][0] );
  
    } else if (a_name == "4") {
  
      /* Fourth order ARK scheme with L-stable implicit part.
         This method has one explicit stage and four implicit stages.
         References: Kennedy and Carpenter 2003.
      */
      a_nstages = 6;
      a_order = 4;
      a_order_emb = 0;
  
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
  
      setCoefficients(  a_Ai, a_Ae,
                        a_bi, a_be,
                        a_c,
                        a_bi_emb, a_be_emb,
                        a_binterpi, a_binterpe,
                        a_nstages,
                        &Ae[0][0], &Ai[0][0],
                        NULL, NULL,
                        1,
                        NULL, &binterpi[0][0] );
  
    } else if (a_name == "5") {
  
      /* 5th order method with L-stable implicit part. 
       * 1 explicit, 5 implicit stages.
       * Ref: Kennedy & Carpenter
       */
  
      a_nstages = 8;
      a_order = 5;
      a_order_emb = 0;
  
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
  
      setCoefficients(  a_Ai, a_Ae,
                        a_bi, a_be,
                        a_c,
                        a_bi_emb, a_be_emb,
                        a_binterpi, a_binterpe,
                        a_nstages,
                        &Ae[0][0], &Ai[0][0],
                        NULL, NULL,
                        1,
                        NULL, &binterpi[0][0] );
  
    } else if (a_name == "rk4") {
  
      /* explicit RK4, disguised as an ARK method -- for debugging */
      a_nstages = 4;
      a_order = 4;
      a_order_emb = 0;
  
      const Real
        A[4][4] = {{0,0,0,0},
                   {0.5,0,0,0},
                   {0,0.5,0,0},
                   {0,0,1.0,0}},
        b[4]    = {1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0};
  
      setCoefficients(  a_Ai, a_Ae,
                        a_bi, a_be,
                        a_c,
                        a_bi_emb, a_be_emb,
                        a_binterpi, a_binterpe,
                        a_nstages,
                        &A[0][0], &A[0][0],
                        &b[0], &b[0],
                        0,
                        NULL, NULL );
  
    } else {
  
      /* default: ARK4 */
      if (!procID()) {
        cout  << "Error: unknown ARK method specified " 
              << a_name << ". Using default.\n";
      }

      a_nstages = 6;
      a_order = 4;
      a_order_emb = 0;
  
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
  
      setCoefficients(  a_Ai, a_Ae,
                        a_bi, a_be,
                        a_c,
                        a_bi_emb, a_be_emb,
                        a_binterpi, a_binterpe,
                        a_nstages,
                        &Ae[0][0], &Ai[0][0],
                        NULL, NULL,
                        1,
                        NULL, &binterpi[0][0] );
  
    }
    return;
  }

}

#include "NamespaceFooter.H"
#endif
