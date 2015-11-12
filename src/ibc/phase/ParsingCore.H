//Parsing class for arbitrary equation
//Written by Wonjae Lee
//
#ifndef _PARSEEQUATION_H_
#define _PARSEEQUATION_H_
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <iostream>
#include <sstream>


namespace ParsingSpace{
    enum DataKind{
        OPERATOR_   = 0, 
        VALUE_      = 1 
    };
    
    enum ErrorCode{
        NOERR_            = 0, 
        ERR_PREPROC_EMAN_ = 1,
        ERR_ITP_          = 2,	
        ERR_PFEV_	      = 3
    };

    enum PriorityType{
        PR1_   = 0, //OPENPAREN_
        PR2_   = 10, //ADD_, SUBTRACT_   : INFIX_
        PR3_   = 20, //MULTI_, DIVIDE_   : INFIX_
        PR4_   = 25, //PPHM_         : INFIX_
        PR5_   = 30, //SIN_, COS_    : PREFIX_
        PR6_   = 31, //PLUS_, MINUS_ : PREFIX_
        PR7_   = 35, //PMHM_         : INFIX_
        PR8_   = 40, //POWER_, NSQRT_    : INFIX_
        PR9_   = 50, //INVERSE_, SQUARE_ : POSTFIX_
        PR10_  = 100 //EP_, EI_          : INFIX_
    };

   enum OperType{
        PREFIX_  = 3,
        MACONST_ = 5,
        INFIX_   = 6,
        POSTFIX_ = 9,
        VAR_     = 7,
        DEFAULT_ = 29 //default type
    };

    enum OperCode{
        NOTDEF_ = 299 ,
        OPENPAREN_    ,
        CLOSEPAREN_   ,
        PLUS_         ,
        MINUS_        , 
        MULTI_        ,
        DIVIDE_       ,
        SIN_          , 
        COS_          ,
        TAN_          ,
        CSC_          ,
        SEC_          ,
        COT_          ,
        ASIN_         ,
        ACOS_         ,
        ATAN_         ,
        ABS_          ,
        LOG_          ,
        LN_           ,
        TENPOWER_     ,
        EXPF_         ,
        SQRT_         ,
        SINH_         ,
        COSH_         ,
        TANH_         ,
        POWER_        ,
        NSQRT_        ,
        INVERSE_      ,
        SQUARE_       ,
        ADD_          ,
        SUBTRACT_     ,
        PMHM_         , //past macro hidden multiplication
        PPHM_         , //past prefix hidden multiplication
        PI__          ,
        X__           ,
        Y__           ,
        Z__           ,
        MU__          ,
        VPAR__        ,
        ANS__         ,
        EP_           ,
        EI_            
    };
}
class OperatorData{
    int	m_preInPost;//operator characteristic {PREFIX_,INFIX_,POSTFIX_} 
    int	m_priority;//operation priority[1,2,3,4,5...]
    int	m_idenCode;//operator code SIN_ COS_..
    public:
    OperatorData(){m_preInPost=0;m_priority=0;m_idenCode=0;};
    int setOperatorData(std::string ops);// cos sin 

    int getPreInPost(){return m_preInPost;}//
    int getPriority(){return m_priority;}//
    int getIdenCode(){return m_idenCode;}

    int operToCode(std::string ops);//
    std::string codeToStr();//
    std::string codeToStr(int opercode);//
};

struct stackData{
    OperatorData	obOperatorData;	
    double		val; 

    struct stackData *down; 
};
typedef struct stackData stackData;

class UserStack{
    stackData	*top;
    int	         cnt;
    public:
    UserStack(){cnt=0;top=NULL;}

    void		dPush(double argdou);
    void		operStrPush(std::string paramStr);

    OperatorData	oPop();
    double		dPop();

    OperatorData	oTop();
    double		dTop();

    int		empty();
    ~UserStack();
};

struct queueData{
    int	   kind;
    struct sQ{
        OperatorData	 obOperatorData;
        std::string	 obCstr;
    }sQv;
    struct queueData *next;
};
typedef struct queueData queueData;

class UserQueue{
    queueData  *head;
    queueData  *rear;
    int         cnt;
    public:
    UserQueue(){head=rear=NULL;cnt=0;}

    void enqueue(queueData	&buffQData);
    void operStrEnqueue(std::string paramStr);
    void valStrEnqueue(std::string paramStr);

    queueData	dequeue();
    int showCnt(){return cnt;}
    ~UserQueue();
};

class ParsingCore{
    bool	isRad; 
    double	result; 
    double  exAnswer;
    UserQueue	postFix; 
    std::string	formula; 
    std::string postStr;
    std::string manipStr;

    double	m_x; 
    double	m_y; 
    double	m_z; 
    double	m_mu; 
    double	m_vpar; 

    public:	

    ParsingCore(){isRad=true;exAnswer=result=0.0;m_x=0.0;m_y=0.0;m_z=0.0;m_vpar=0.0;m_mu=0.0;}
    ParsingCore(const char *pure ){isRad=true;exAnswer=result=0.0;m_x=m_y=m_z=m_vpar=m_mu=0.0;setFormula(pure);preProcess();inToPostFix();postFixEvaluate();}
    ParsingCore(double ax,const char *pure ){isRad=true;exAnswer=result=0.0;m_x=ax;m_y=m_z=m_vpar=m_mu=0.0;setFormula(pure);preProcess();inToPostFix();postFixEvaluate();}
    ParsingCore(double ax,double ay, const char *pure ){isRad=true;exAnswer=result=0.0;m_x=ax;m_y=ay;m_z=m_vpar=m_mu=0.0;setFormula(pure);preProcess();inToPostFix();postFixEvaluate();}
    ParsingCore(double ax,double ay,double az, const char *pure ){isRad=true;exAnswer=result=0.0;m_x=ax;m_y=ay;m_z=az;m_vpar=m_mu=0.0;setFormula(pure);preProcess();inToPostFix();postFixEvaluate();}
    ParsingCore(double a_x,double a_y,double a_z, const char *pure, double a_vpar ){isRad=true;exAnswer=result=0.0;m_x=a_x;m_y=a_y;m_z=a_z;m_vpar=a_vpar;m_mu=0.0;setFormula(pure);preProcess();inToPostFix();postFixEvaluate();}
    ParsingCore(double a_x,double a_y,double a_z, const char *pure, double a_vpar, double a_mu ){isRad=true;exAnswer=result=0.0;m_x=a_x;m_y=a_y;m_z=a_z;m_vpar=a_vpar;m_mu=a_mu;setFormula(pure);preProcess();inToPostFix();postFixEvaluate();}
    ParsingCore(double paramAns,bool radOrDeg){result=0.0;exAnswer=paramAns;isRad=radOrDeg;}

    double calc2d(double a_x, double a_y){m_x=a_x;m_y=a_y;inToPostFix();postFixEvaluate();return result;}
    double calc3d(double a_x, double a_y, double a_z){m_x=a_x;m_y=a_y;m_z=a_z;inToPostFix();postFixEvaluate();return result;}
    double calc4d(double a_x, double a_y, double a_vpar, double a_mu){m_x=a_x;m_y=a_y;m_vpar=a_vpar;m_mu=a_mu;inToPostFix();postFixEvaluate();return result;}
    double calc5d(double a_x, double a_y, double a_z, double a_vpar, double a_mu){m_x=a_x;m_y=a_y;m_z=a_z;m_vpar=a_vpar;m_mu=a_mu;inToPostFix();postFixEvaluate();return result;}

    std::string getFormula(){return formula;}
    std::string getManipStr(){return manipStr;}
    std::string getPostStr(){return postStr;}
    void setFormula(std::string pure){formula=pure;}
    void setFormula(char* pure){formula.assign(pure);}
    int preProcess();
    int inToPostFix();
    int postFixEvaluate();
    double calculate(OperatorData op,double operand1);
    double calculate(double operand1, OperatorData op, double operand2);
    double calculate(double operand2, OperatorData op);
    int noEnterable(std::string currentStr, std::string stackStr);

    double showResult(){return result;}
};
#endif
