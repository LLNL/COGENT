//Parse Equation
#include "ParsingCore.H"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <iostream>
#include <sstream>

//
//extern "C" double wrappedParsingCore(char *inputstring) {
//    std::cout<<" From the wrappedParsingCore "<<std::endl;
//    ParsingCore pcore(inputstring);
//    std::cout<<"ParsingCore.formula = "<<pcore.getFormula()<<std::endl;
//    std::cout<<"ParsingCore.manipStr = "<<pcore.getManipStr()<<std::endl;
//    std::cout<<"ParsingCore.postStr = "<<pcore.getPostStr()<<std::endl;
//    std::cout<<"result = "<<pcore.showResult()<<std::endl;
//    return pcore.showResult();
//}
//extern "C" double wrappedParsingCoreX(char *inputstring, double inputx ) {
//    std::cout<<" From the wrappedParsingCore "<<std::endl;
//    ParsingCore pcore(inputx, inputstring);
//    std::cout<<"ParsingCore.formula = "<<pcore.getFormula()<<std::endl;
//    std::cout<<"ParsingCore.manipStr = "<<pcore.getManipStr()<<std::endl;
//    std::cout<<"ParsingCore.postStr = "<<pcore.getPostStr()<<std::endl;
//    std::cout<<"result = "<<pcore.showResult()<<std::endl;
//    return pcore.showResult();
//}
//extern "C" double wrappedParsingCoreXY(char *inputstring, double inputx, double inputy) {
//    std::cout<<" From the wrappedParsingCore "<<std::endl;
//    ParsingCore pcore(inputx, inputy, inputstring);
//    std::cout<<"ParsingCore.formula = "<<pcore.getFormula()<<std::endl;
//    std::cout<<"ParsingCore.manipStr = "<<pcore.getManipStr()<<std::endl;
//    std::cout<<"ParsingCore.postStr = "<<pcore.getPostStr()<<std::endl;
//    std::cout<<"result = "<<pcore.showResult()<<std::endl;
//    return pcore.showResult();
//}
//extern "C" double wrappedParsingCoreXYZ(char *inputstring, double inputx, double inputy, double inputz) {
//    std::cout<<" From the wrappedParsingCore "<<std::endl;
//    ParsingCore pcore(inputx, inputy, inputz, inputstring);
//    std::cout<<"ParsingCore.formula = "<<pcore.getFormula()<<std::endl;
//    std::cout<<"ParsingCore.manipStr = "<<pcore.getManipStr()<<std::endl;
//    std::cout<<"ParsingCore.postStr = "<<pcore.getPostStr()<<std::endl;
//    std::cout<<"result = "<<pcore.showResult()<<std::endl;
//    return pcore.showResult();
//}
//



int OperatorData::operToCode(std::string ops)
{
    if(ops=="sin"||ops=="SIN"||ops=="Sin"){
        return ParsingSpace::SIN_;
    }
    else if(ops=="cos"||ops=="COS"||ops=="Cos"){
        return ParsingSpace::COS_;
    }
    else if(ops=="tan"||ops=="TAN"||ops=="Tan"){
        return ParsingSpace::TAN_;
    }
    else if(ops=="csc"||ops=="CSC"||ops=="Csc"||ops=="cosec"||ops=="Cosec"||ops=="COSEC"){
        return ParsingSpace::CSC_;
    }
    else if(ops=="sec"||ops=="SEC"||ops=="Sec"){
        return ParsingSpace::SEC_;
    }
    else if(ops=="cot"||ops=="COT"||ops=="Cot"){
        return ParsingSpace::COT_;
    }
    else if(ops=="arcsin"||ops=="ARCSIN"||ops=="Arcsin"||ops=="asin"||ops=="ASIN"||ops=="Asin"){
        return ParsingSpace::ASIN_;
    }
    else if(ops=="arccos"||ops=="ARCCOS"||ops=="Arccos"||ops=="acos"||ops=="ACOS"||ops=="Acos"){
        return ParsingSpace::ACOS_;
    }
    else if(ops=="arctan"||ops=="ARCTAN"||ops=="Arctan"||ops=="atan"||ops=="ATAN"||ops=="Atan"){
        return ParsingSpace::ATAN_;
    }
    else if(ops=="abs"||ops=="ABS"||ops=="Abs"){
        return ParsingSpace::ABS_;
    }
    else if(ops=="log"||ops=="LOG"||ops=="Log"){
        return ParsingSpace::LOG_;
    }
    else if(ops=="ln"||ops=="LN"||ops=="Ln"){
        return ParsingSpace::LN_;
    }
    else if(ops=="tenpower"||ops=="TENPOWER"||ops=="Tenpower"){
        return ParsingSpace::TENPOWER_;
    }
    else if(ops=="exp"||ops=="EXP"||ops=="Exp"){
        return ParsingSpace::EXPF_;
    }
    else if(ops=="sqrt"||ops=="SQRT"||ops=="Sqrt"){
        return ParsingSpace::SQRT_;
    }
    else if(ops=="sinh"||ops=="SINH"||ops=="Sinh"){
        return ParsingSpace::SINH_;
    }
    else if(ops=="cosh"||ops=="COSH"||ops=="Cosh"){
        return ParsingSpace::COSH_;
    }
    else if(ops=="tanh"||ops=="TANH"||ops=="Tanh"){
        return ParsingSpace::TANH_;
    }
    else if(ops=="power"||ops=="POWER"||ops=="Power"){
        return ParsingSpace::POWER_;
    }
    else if(ops=="nsqrt"||ops=="NSQRT"||ops=="Nsqrt"){
        return ParsingSpace::NSQRT_;
    }
    else if(ops=="inverse\""||ops=="INVERSE\""||ops=="Inverse\""){
        return ParsingSpace::INVERSE_;
    }
    else if(ops=="square\""||ops=="SQUARE\""||ops=="Square\""){
        return ParsingSpace::SQUARE_;
    }
    else if(ops=="*"){
        return ParsingSpace::MULTI_;
    }
    else if(ops=="/"){
        return ParsingSpace::DIVIDE_;
    }
    else if(ops=="^"){
        return ParsingSpace::POWER_;
    }
    else if(ops=="+:"){
        return ParsingSpace::PLUS_;
    }
    else if(ops=="-:"){
        return ParsingSpace::MINUS_;
    }
    else if(ops=="+"){
        return ParsingSpace::ADD_;
    }
    else if(ops=="-"){
        return ParsingSpace::SUBTRACT_;
    }
    else if(ops==":"){
        return ParsingSpace::PMHM_;
    }
    else if(ops==";"){
        return ParsingSpace::PPHM_;
    }
    else if(ops=="Ep"){
        return ParsingSpace::EP_;
    }
    else if(ops=="Ei"){
        return ParsingSpace::EI_;
    }
    else if(ops=="("){
        return ParsingSpace::OPENPAREN_;
    }
    else if(ops==")"){
        return ParsingSpace::CLOSEPAREN_;
    }
    else if(ops=="pi\'"||ops=="PI\'"||ops=="Pi\'"){
        return ParsingSpace::PI__;
    }
    else if(ops=="x\'"||ops=="X\'"){
        return ParsingSpace::X__;
    }
    else if(ops=="y\'"||ops=="Y\'"){
        return ParsingSpace::Y__;
    }
    else if(ops=="z\'"||ops=="Z\'"){
        return ParsingSpace::Z__;
    }
    else if(ops=="mu\'"||ops=="MU\'"||ops=="Mu\'"){
        return ParsingSpace::MU__;
    }
    else if(ops=="vpar\'"||ops=="VPAR\'"||ops=="Vpar\'"){
        return ParsingSpace::VPAR__;
    }
    else if(ops=="ans\'"||ops=="ANS\'"||ops=="Ans\'"){
        return ParsingSpace::ANS__;
    }

    else{ 
        return ParsingSpace::NOTDEF_;
    }
}

std::string OperatorData::codeToStr()
{
    std::string  tempStr="notDef";
    switch(m_idenCode){
        case ParsingSpace::OPENPAREN_  :tempStr="(";break;
        case ParsingSpace::CLOSEPAREN_ :tempStr=")";break;
        case ParsingSpace::SIN_        :tempStr="sin"; break;
        case ParsingSpace::COS_        :tempStr="cos"; break;
        case ParsingSpace::TAN_        :tempStr="tan"; break;
        case ParsingSpace::CSC_        :tempStr="csc"; break;
        case ParsingSpace::SEC_        :tempStr="sec"; break;
        case ParsingSpace::COT_        :tempStr="cot"; break;
        case ParsingSpace::ASIN_       :tempStr="arcsin"; break;
        case ParsingSpace::ACOS_       :tempStr="arccos"; break;
        case ParsingSpace::ATAN_       :tempStr="arctan"; break;
        case ParsingSpace::ABS_        :tempStr="abs"; break;
        case ParsingSpace::LOG_        :tempStr="log"; break;
        case ParsingSpace::LN_         :tempStr="ln"; break;
        case ParsingSpace::TENPOWER_   :tempStr="tenpower"; break;
        case ParsingSpace::EXPF_       :tempStr="exp"; break;
        case ParsingSpace::SQRT_       :tempStr="sqrt"; break;
        case ParsingSpace::SINH_       :tempStr="sinh"; break;
        case ParsingSpace::COSH_       :tempStr="cosh"; break;
        case ParsingSpace::TANH_       :tempStr="tanh"; break;
        case ParsingSpace::POWER_      :tempStr="^"; break;
        case ParsingSpace::NSQRT_      :tempStr="nsqrt"; break;
        case ParsingSpace::INVERSE_    :tempStr="inverse\""; break;
        case ParsingSpace::SQUARE_     :tempStr="square\""; break;
        case ParsingSpace::MULTI_      :tempStr="*"; break;
        case ParsingSpace::DIVIDE_     :tempStr="/"; break;
        case ParsingSpace::PLUS_       :tempStr="+:"; break;
        case ParsingSpace::MINUS_      :tempStr="-:"; break;
        case ParsingSpace::ADD_        :tempStr="+"; break;
        case ParsingSpace::SUBTRACT_   :tempStr="-"; break;
        case ParsingSpace::PMHM_       :tempStr=":"; break;
        case ParsingSpace::PPHM_       :tempStr=";"; break;
        case ParsingSpace::EP_         :tempStr="Ep"; break;
        case ParsingSpace::EI_         :tempStr="Ei"; break;
        case ParsingSpace::X__         :tempStr="x\'"; break;
        case ParsingSpace::Y__         :tempStr="y\'"; break;
        case ParsingSpace::Z__         :tempStr="z\'"; break;
        case ParsingSpace::MU__        :tempStr="mu\'"; break;
        case ParsingSpace::VPAR__      :tempStr="vpar\'"; break;
        case ParsingSpace::NOTDEF_     :tempStr="NOTDEF_"; break;
        default                        :tempStr="default"; break;
    }
    return tempStr;
}

std::string OperatorData::codeToStr(int opercode)
{
    std::string  tempStr="notDef";
    switch(opercode){
        case ParsingSpace::OPENPAREN_  :tempStr="(";break;
        case ParsingSpace::CLOSEPAREN_ :tempStr=")";break;
        case ParsingSpace::SIN_        :tempStr="sin"; break;
        case ParsingSpace::COS_        :tempStr="cos"; break;
        case ParsingSpace::TAN_        :tempStr="tan"; break;
        case ParsingSpace::CSC_        :tempStr="csc"; break;
        case ParsingSpace::SEC_        :tempStr="sec"; break;
        case ParsingSpace::COT_        :tempStr="cot"; break;
        case ParsingSpace::ASIN_       :tempStr="arcsin"; break;
        case ParsingSpace::ACOS_       :tempStr="arccos"; break;
        case ParsingSpace::ATAN_       :tempStr="arctan"; break;
        case ParsingSpace::ABS_        :tempStr="abs"; break;
        case ParsingSpace::LOG_        :tempStr="log"; break;
        case ParsingSpace::LN_         :tempStr="ln"; break;
        case ParsingSpace::TENPOWER_   :tempStr="tenpower"; break;
        case ParsingSpace::EXPF_       :tempStr="exp"; break;
        case ParsingSpace::SQRT_       :tempStr="sqrt"; break;
        case ParsingSpace::SINH_       :tempStr="sinh"; break;
        case ParsingSpace::COSH_       :tempStr="cosh"; break;
        case ParsingSpace::TANH_       :tempStr="tanh"; break;
        case ParsingSpace::POWER_      :tempStr="^"; break;
        case ParsingSpace::NSQRT_      :tempStr="nsqrt"; break;
        case ParsingSpace::INVERSE_    :tempStr="inverse\""; break;
        case ParsingSpace::SQUARE_     :tempStr="square\""; break;
        case ParsingSpace::MULTI_      :tempStr="*"; break;
        case ParsingSpace::DIVIDE_     :tempStr="/"; break;
        case ParsingSpace::PLUS_       :tempStr="+:"; break;
        case ParsingSpace::MINUS_      :tempStr="-:"; break;
        case ParsingSpace::ADD_        :tempStr="+"; break;
        case ParsingSpace::SUBTRACT_   :tempStr="-"; break;
        case ParsingSpace::PMHM_       :tempStr=":"; break;
        case ParsingSpace::PPHM_       :tempStr=";"; break;
        case ParsingSpace::EP_         :tempStr="Ep"; break;
        case ParsingSpace::EI_         :tempStr="Ei"; break;
        case ParsingSpace::X__         :tempStr="x\'"; break;
        case ParsingSpace::Y__         :tempStr="y\'"; break;
        case ParsingSpace::Z__         :tempStr="z\'"; break;
        case ParsingSpace::MU__        :tempStr="mu\'"; break;
        case ParsingSpace::VPAR__      :tempStr="vpar\'"; break;
        case ParsingSpace::NOTDEF_     :tempStr="NOTDEF_"; break;
        default                        :tempStr="default"; break;
    }
    return tempStr;
}

int OperatorData::setOperatorData(std::string ops)
{
    //initializing
    m_preInPost=ParsingSpace::DEFAULT_;
    m_priority=ParsingSpace::DEFAULT_;
    m_idenCode=ParsingSpace::DEFAULT_;

    m_idenCode=operToCode(ops);
    switch(m_idenCode){
        case ParsingSpace::SIN_:
            m_preInPost=ParsingSpace::PREFIX_;
            m_priority=ParsingSpace::PR5_;//30
            break;
        case ParsingSpace::COS_:
            m_preInPost=ParsingSpace::PREFIX_;
            m_priority=ParsingSpace::PR5_;//30
            break;
        case ParsingSpace::TAN_:
            m_preInPost=ParsingSpace::PREFIX_;
            m_priority=ParsingSpace::PR5_;//30;
            break;
        case ParsingSpace::CSC_:
            m_preInPost=ParsingSpace::PREFIX_;
            m_priority=ParsingSpace::PR5_;//30;
            break;
        case ParsingSpace::SEC_:
            m_preInPost=ParsingSpace::PREFIX_;
            m_priority=ParsingSpace::PR5_;//30;
            break;
        case ParsingSpace::COT_:
            m_preInPost=ParsingSpace::PREFIX_;
            m_priority=ParsingSpace::PR5_;//30;
            break;
        case ParsingSpace::ASIN_:
            m_preInPost=ParsingSpace::PREFIX_;
            m_priority=ParsingSpace::PR5_;//30;
            break;
        case ParsingSpace::ACOS_:
            m_preInPost=ParsingSpace::PREFIX_;
            m_priority=ParsingSpace::PR5_;//30;
            break;
        case ParsingSpace::ATAN_:
            m_preInPost=ParsingSpace::PREFIX_;
            m_priority=ParsingSpace::PR5_;//30;
            break;
        case ParsingSpace::ABS_:
            m_preInPost=ParsingSpace::PREFIX_;
            m_priority=ParsingSpace::PR5_;//30;
            break;
        case ParsingSpace::LOG_:
            m_preInPost=ParsingSpace::PREFIX_;
            m_priority=ParsingSpace::PR5_;//30;
            break;
        case ParsingSpace::LN_:
            m_preInPost=ParsingSpace::PREFIX_;
            m_priority=ParsingSpace::PR5_;//30;
            break;
        case ParsingSpace::TENPOWER_:
            m_preInPost=ParsingSpace::PREFIX_;
            m_priority=ParsingSpace::PR5_;//30;
            break;
        case ParsingSpace::EXPF_:
            m_preInPost=ParsingSpace::PREFIX_;
            m_priority=ParsingSpace::PR5_;//30;
            break;
        case ParsingSpace::SQRT_:
            m_preInPost=ParsingSpace::PREFIX_;
            m_priority=ParsingSpace::PR5_;//30;
            break;
        case ParsingSpace::SINH_:
            m_preInPost=ParsingSpace::PREFIX_;
            m_priority=ParsingSpace::PR5_;//30;
            break;
        case ParsingSpace::COSH_:
            m_preInPost=ParsingSpace::PREFIX_;
            m_priority=ParsingSpace::PR5_;//30;
            break;
        case ParsingSpace::TANH_:
            m_preInPost=ParsingSpace::PREFIX_;
            m_priority=ParsingSpace::PR5_;//30;
            break;
        case ParsingSpace::POWER_:
            m_preInPost=ParsingSpace::INFIX_;
            m_priority=ParsingSpace::PR8_;//40;
            break;
        case ParsingSpace::NSQRT_:
            m_preInPost=ParsingSpace::INFIX_;
            m_priority=ParsingSpace::PR8_;//40;
            break;
        case ParsingSpace::PLUS_:
            m_preInPost=ParsingSpace::PREFIX_;
            m_priority=ParsingSpace::PR6_;//31;     
            break;
        case ParsingSpace::MINUS_:
            m_preInPost=ParsingSpace::PREFIX_;
            m_priority=ParsingSpace::PR6_;//31;	
            break;
        case ParsingSpace::ADD_:
            m_preInPost=ParsingSpace::INFIX_;
            m_priority=ParsingSpace::PR2_;//10;
            break;
        case ParsingSpace::SUBTRACT_:
            m_preInPost=ParsingSpace::INFIX_;
            m_priority=ParsingSpace::PR2_;//10;
            break;
        case ParsingSpace::OPENPAREN_:
            m_priority=ParsingSpace::PR1_;//0;
            break;
        case ParsingSpace::INVERSE_:
            m_preInPost=ParsingSpace::POSTFIX_;
            m_priority=ParsingSpace::PR9_;//50;
            break;
        case ParsingSpace::SQUARE_:
            m_preInPost=ParsingSpace::POSTFIX_;
            m_priority=ParsingSpace::PR9_;//50;
            break;
        case ParsingSpace::MULTI_:
            m_preInPost=ParsingSpace::INFIX_;
            m_priority=ParsingSpace::PR3_;//20;
            break;
        case ParsingSpace::DIVIDE_:
            m_preInPost=ParsingSpace::INFIX_;
            m_priority=ParsingSpace::PR3_;//20;
            break;
        case ParsingSpace::PI__:
            m_preInPost=ParsingSpace::MACONST_;
            break;
        case ParsingSpace::X__:
            m_preInPost=ParsingSpace::VAR_;
            break;
        case ParsingSpace::Y__:
            m_preInPost=ParsingSpace::VAR_;
            break;
        case ParsingSpace::Z__:
            m_preInPost=ParsingSpace::VAR_;
            break;
        case ParsingSpace::MU__:
            m_preInPost=ParsingSpace::VAR_;
            break;
        case ParsingSpace::VPAR__:
            m_preInPost=ParsingSpace::VAR_;
            break;
        case ParsingSpace::ANS__:
            m_preInPost=ParsingSpace::MACONST_;
            break;
        case ParsingSpace::PMHM_:
            m_preInPost=ParsingSpace::INFIX_;
            m_priority=ParsingSpace::PR7_;//35;
            break;
        case ParsingSpace::PPHM_:
            m_preInPost=ParsingSpace::INFIX_;
            m_priority=ParsingSpace::PR4_;//25;
            break;
        case ParsingSpace::EP_:
            m_preInPost=ParsingSpace::PREFIX_;
            m_priority=ParsingSpace::PR10_;//100;
            break;
        case ParsingSpace::EI_:
            m_preInPost=ParsingSpace::INFIX_;
            m_priority=ParsingSpace::PR10_;//100;
            break;
        default:
            m_idenCode=ParsingSpace::NOTDEF_;
            m_preInPost=ParsingSpace::NOTDEF_;
            return 1;
    }
    return 0;
}

void UserStack::operStrPush(std::string paramStr)
{
    stackData	*p;
    p=new stackData;
    if(p==NULL){
        exit(1);
    }
    p->down=top;
    top=p;
    p->obOperatorData.setOperatorData(paramStr);
    cnt++;
}
UserStack::~UserStack()
{	
    while(cnt){
        stackData *p;
        p=top;
        top=p->down;
        delete p;
        cnt--;
    }
}		

void UserStack::dPush(double argDou)
{
    stackData  *p;
    p=new stackData;
    if(p==NULL){
        exit(1);
    }
    p->down=top;
    top=p;
    p->val=argDou;
    cnt++;
}
OperatorData UserStack::oPop()
{
    stackData  *p;
    OperatorData temp;
    if(!cnt){
        temp.setOperatorData("!");
        return temp;
    }
    temp=top->obOperatorData;
    p=top;
    top=p->down;
    cnt--;
    delete p;
    return temp;
}
double UserStack::dPop()
{
    stackData  *p;
    double	temp=0.0;
    if(!cnt||top==NULL){
        return temp;
    }
    temp=top->val;
    p=top;
    top=p->down;
    cnt--;
    delete p;
    return temp;
}
OperatorData UserStack::oTop()
{
    OperatorData temp;
    if(!cnt||top==NULL){
        temp.setOperatorData("!");
        return temp;
    }
    temp=top->obOperatorData;
    return temp;
}
double UserStack::dTop()
{
    double temp=0.0;
    if(!cnt||top==NULL){
        return temp;
    }
    temp = top->val;
    return temp;
}
int UserStack::empty()
{
    if(cnt==0)
        return 1;
    else
        return 0;
}

void UserQueue::operStrEnqueue(std::string paramStr)
{
    queueData  *p;
    p=new queueData;
    if(p==NULL){
        exit(1);
    }
    if(head==NULL||rear==NULL){
        head=rear=p;
    }
    else{
        rear->next=p;
        rear=p;
    }
    p->kind=ParsingSpace::OPERATOR_;
    p->sQv.obOperatorData.setOperatorData(paramStr);
    p->next=NULL;
    cnt++;
}
void UserQueue::valStrEnqueue(std::string paramStr)
{
    queueData *p;
    p=new queueData;
    if(p==NULL){
        exit(1);
    }
    if(head==NULL||rear==NULL){
        head=rear=p;
    }
    else{
        rear->next=p;
        rear=p;
    }
    p->kind=ParsingSpace::VALUE_;
    p->sQv.obCstr=paramStr;
    p->next=NULL;
    cnt++;
}
void UserQueue::enqueue(queueData &buffQData)
{
    queueData *p;
    p = new queueData;
    if(p==NULL){
        exit(1);
    }
    if(head==NULL||rear==NULL){
        head=rear=p;
    }
    else{
        rear->next=p;
        rear=p;
    }

    *p=buffQData;

    p->next=NULL;
    cnt++;
}
queueData UserQueue::dequeue()
{
    queueData temp;
    queueData *p;

    temp=*head;
    temp.next=NULL;
    p=head;
    head=head->next;
    delete p;
    cnt--;
    return temp;
}
UserQueue::~UserQueue()
{
    while(cnt){
        queueData *p;
        p=head;
        head=head->next;
        delete p;
        cnt--;
    }
}


int ParsingCore::preProcess()
{
    manipStr=formula;
    std::string postFixList[6]={"inverse","INVERSE","Inverse","square","SQUARE","Square"};
    std::string macroList1[6]={"pi","PI","Pi","ans","ANS","Ans"};
    std::string varList1[12]={"x","X","y","Y","z","Z","mu","MU","Mu","vpar","VPAR","Vpar"};
    std::string preFixList[]={"arcsin","ARCSIN","Arcsin","asin","ASIN","Asin",
        "arccos","ARCCOS","Arccos","acos","ACOS","Acos",
        "arctan","ARCTAN","Arctan","atan","ATAN","Atan", 
        "abs","ABS","Abs",
        "sinh","SINH","Sinh","cosh","COSH","Cosh","tanh","TANH","Tanh",
        "sin","SIN","Sin","cos","COS","Cos","tan","TAN","Tan",
        "csc","CSC","Csc","cosec","Cosec","COSEC",
        "sec","SEC","Sec","cot","COT","Cot",
        "log","LOG","Log","ln","LN","Ln","tenpower","TENPOWER","Tenpower",
        "exp","EXP","Exp","sqrt","SQRT","Sqrt"};

    std::string pStr;
    OperatorData serviceOperData;

    int i=-6;
    int j,tempCode=0;
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 0"<<std::endl;
#endif
//debug

    while( (i= manipStr.find(" "))!=-1){
        manipStr.erase(i,1);		
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 1"<<std::endl;
#endif
//debug
    }

    i=-6;
    for(j=0;j<6; j++,i=-6){
//debug
#ifdef DEBUG_PARSER
        std::cout<<"postFixList[j]="<< postFixList[j]<<"  : preprocess 2"<<std::endl;
        std::cout<<"i=manipStr.find(postFixList[j],i+6))="<<  manipStr.find(postFixList[j],i+6)<<"  : preprocess 2"<<std::endl;
#endif
//debug
        while( (i=manipStr.find(postFixList[j],i+6))!=-1){
            if(j<=2){
                manipStr.insert(i+7,"\"");
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 2"<<std::endl;
#endif
//debug
            }
            else{
                manipStr.insert(i+6,"\"");			
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 3"<<std::endl;
#endif
//debug
            }
        }
    }

    i=-2;
    for(j=0;j<6;j++,i=-2){
        while( (i=manipStr.find(macroList1[j],i+2))!=-1){
            if(j<=2){
                manipStr.insert(i+2,"\'");
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 4"<<std::endl;
#endif
//debug
            }
            else{
                manipStr.insert(i+3,"\'");
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 5"<<std::endl;
#endif
//debug
            }
        }
    }

    i=-2;
    for(j=0;j<12;j++,i=-2){
        while( (i=manipStr.find(varList1[j],i+2))!=-1){
            if (j<=1){
                if (i==0){
                    manipStr.insert(i+1,"\'");
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 5.1"<<std::endl;
#endif
//debug
                }
                else if(i==manipStr.length()){
                    manipStr.insert(i+1,"\'");
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 5.2"<<std::endl;
#endif
//debug
                }
                else if(manipStr[i+1]=='p'||manipStr[i+1]=='P'){
                   continue;
                }
                else{
                    manipStr.insert(i+1,"\'");
                }

            }
            else if(j<=5){
                manipStr.insert(i+1,"\'");
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 5.3"<<std::endl;
#endif
//debug
            }
            else if (j<=8){
                manipStr.insert(i+2,"\'");
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 5.4"<<std::endl;
#endif
//debug
            }
            else {
                manipStr.insert(i+4,"\'");
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 5.5"<<std::endl;
#endif
            }
        }
    }

    i=-2;
    for(j=0;j<54;j++,i=-2){
        while( (i=manipStr.find(preFixList[j],i+2))!=-1){
            if( i!=0 && (manipStr[i-1]==')'||isdigit(manipStr[i-1])||manipStr[i-1]=='\''||manipStr[i-1]=='\"') ){
                manipStr.insert(i,";");
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 6"<<std::endl;
#endif
//debug
            }
        }
    }
    i=-2;
    for(j=0;j<6;j++,i=-2){
        while( (i=manipStr.find(macroList1[j],i+2))!=-1){
            if( i!=0 && (manipStr[i-1]==')'||isdigit(manipStr[i-1])||manipStr[i-1]=='\''||manipStr[i-1]=='\"') ){
                manipStr.insert(i,":");
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 7"<<std::endl;
#endif
//debug
            }
        }
    }	
    i=-2;
    while( (i=manipStr.find('(',i+2))!=-1){
        if( i!=0 && (manipStr[i-1]==')'||isdigit(manipStr[i-1])||manipStr[i-1]=='\''||manipStr[i-1]=='\"') ){
            manipStr.insert(i,"*");
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 8"<<std::endl;
#endif
//debug
        }
    }

    for(i=0; i<manipStr.length();){
        if(isalpha(manipStr[i]) ){
        pStr.resize(0);
        pStr+=manipStr[i];
            while( i<= manipStr.length()-2){
                if( isalpha(manipStr[i+1])||manipStr[i+1]=='\''||manipStr[i+1]=='\"' ){
                    pStr+=manipStr[++i];
//debug
#ifdef DEBUG_PARSER
        std::cout<<"pStr="<<pStr<<"  : preprocess 9"<<std::endl;
#endif
//debug
                    if( (tempCode=serviceOperData.operToCode(pStr))==ParsingSpace::NOTDEF_ ){
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 9"<<std::endl;
#endif
//debug
                        continue;
                    }
                    else if(tempCode==ParsingSpace::SIN_||tempCode==ParsingSpace::COS_||tempCode==ParsingSpace::TAN_){
                        if( i<=manipStr.length()-2 && (manipStr[i+1]=='h'||manipStr[i+1]=='H') )
                            pStr+=manipStr[++i];
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 10"<<std::endl;
#endif
//debug
                        break;
                    }
                    else{
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 11"<<std::endl;
#endif
//debug
                        break;
                    }
                }
                else if( manipStr[i]=='E'){
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 11.1"<<std::endl;
#endif
//debug
                    if(i==0){
                        manipStr.insert(++i,"p");
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 12"<<std::endl;
#endif
//debug
                    }
                    else if( isdigit(manipStr[i-1])||manipStr[i-1]=='.' ){
                        manipStr.insert(++i,"i");
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 13"<<std::endl;
#endif
//debug
                    }
                    else{
                        manipStr.insert(++i,"p");
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 14"<<std::endl;
#endif
//debug
                    }

                    break;
                }
                else{
//debug
#ifdef DEBUG_PARSER
        std::cout<<"manipStr="<<manipStr<<"  : preprocess 15"<<std::endl;
#endif
//debug
                    return ParsingSpace::ERR_PREPROC_EMAN_;
                }
            }
            i++;
        }
        else 
            i++;
    }


    return ParsingSpace::NOERR_;
}

int ParsingCore::noEnterable(std::string currentStr,std::string stackStr)
{
    OperatorData oper1,oper2;
    oper1.setOperatorData(currentStr);
    oper2.setOperatorData(stackStr);

    switch(oper2.getPriority() ){
        case ParsingSpace::PR1_ ://0:
            return 0;
        case ParsingSpace::PR2_ ://10:
            if(oper1.getPriority()==ParsingSpace::PR2_){//10 
                return 1;
            }
            else{
                return 0;
            }
        case ParsingSpace::PR3_ ://20:
            if(oper1.getPriority()==ParsingSpace::PR2_ ||oper1.getPriority()==ParsingSpace::PR3_){
                //if(oper1.getPriority()==10 ||oper1.getPriority()==20){
                return 1;
            }
            else{
                return 0;
            }
        case ParsingSpace::PR4_ ://25:
            if(oper1.getPriority()==ParsingSpace::PR2_||oper1.getPriority()==ParsingSpace::PR3_||oper1.getPriority()==ParsingSpace::PR4_){
            //if(oper1.getPriority()==10||oper1.getPriority()==20||oper1.getPriority()==25){
                return 1;
            }
            else{ 
                return 0;
            }
        case ParsingSpace::PR5_ : //30
        case ParsingSpace::PR6_ : //31
            if(oper1.getPriority()==ParsingSpace::PR2_||oper1.getPriority()==ParsingSpace::PR3_||oper1.getPriority()==ParsingSpace::PR4_){
            //if(oper1.getPriority()==10||oper1.getPriority()==20||oper1.getPriority()==25){
                return 1;
            }
            else{
                return 0;
            }

        case ParsingSpace::PR7_ : //35
            if(oper1.getPriority()==ParsingSpace::PR2_ ||oper1.getPriority()==ParsingSpace::PR3_ ||oper1.getPriority()==ParsingSpace::PR4_ ||oper1.getPriority()==ParsingSpace::PR7_ ){
            //if(oper1.getPriority()==10||oper1.getPriority()==20||oper1.getPriority()==25||oper1.getPriority()==35){
                return 1;
            }
            else{ 
                return 0;
            }
        case ParsingSpace::PR8_ : //40
            if(oper1.getPriority()==ParsingSpace::PR2_ ||oper1.getPriority()==ParsingSpace::PR3_  || oper1.getPriority()==ParsingSpace::PR4_ ||oper1.getPriority()==ParsingSpace::PR7_ ||oper1.getPriority()==ParsingSpace::PR8_ ){
                return 1;
            }
            else{
                return 0;
            }
        case ParsingSpace::PR9_ : //50:
            return 1;
        case ParsingSpace::PR10_ ://100:
            if(oper1.getPriority()==ParsingSpace::PR6_){
            //if(oper1.getPriority()==31){
                return 0;
            }
            else {
                return 1;
            }
    }
    return 0;
}

int ParsingCore::inToPostFix()
{
    UserStack		operStk;
    OperatorData	serviceOperData;
    std::string    	parserStr, topParserStr, outParserStr;
    int tempCode=ParsingSpace::NOTDEF_;

    //DEBUG
    postStr.resize(0);//Empty();

//debug
#ifdef DEBUG_PARSER
        std::cout<<"postStr="<<postStr<<" |parsterStr="<<parserStr<<" |topParserStr="<<topParserStr<<" |outParserStr="<<outParserStr<< "  || inToPostFix()1"<<std::endl;
#endif
//debug

    for(int i=0;i<manipStr.length();i++){
        parserStr.resize(0);//Empty();

        parserStr+=manipStr[i];
//debug
#ifdef DEBUG_PARSER
        std::cout<<"postStr="<<postStr<<" |parsterStr="<<parserStr<<" |topParserStr="<<topParserStr<<" |outParserStr="<<outParserStr<< "  || inToPostFix()2"<<"|i="<<i<<std::endl;
#endif
//debug

        if(	parserStr=="("	){
            operStk.operStrPush(parserStr);
//debug
#ifdef DEBUG_PARSER
        std::cout<<"postStr="<<postStr<<" |parsterStr="<<parserStr<<" |topParserStr="<<topParserStr<<" |outParserStr="<<outParserStr<< "  || inToPostFix()3"<<"|i="<<i<<std::endl;
#endif
//debug
        }
        else if(	parserStr==")"	){
            if(!operStk.empty() ){
                parserStr=operStk.oPop().codeToStr();
//debug
#ifdef DEBUG_PARSER
        std::cout<<"postStr="<<postStr<<" |parsterStr="<<parserStr<<" |topParserStr="<<topParserStr<<" |outParserStr="<<outParserStr<< "  || inToPostFix()4"<<"|i="<<i<<std::endl;
#endif
//debug
            }
            else{
                return ParsingSpace::ERR_ITP_;
//debug
#ifdef DEBUG_PARSER
        std::cout<<"postStr="<<postStr<<" |parsterStr="<<parserStr<<" |topParserStr="<<topParserStr<<" |outParserStr="<<outParserStr<< "  || inToPostFix()5"<<"|i="<<i<<std::endl;
#endif
//debug
            }
            while(	parserStr!="("	){
                //DEBUG
                postStr+="{";
                postStr+=parserStr;
                postStr+="}";
                //
                postFix.operStrEnqueue(parserStr);
//debug
#ifdef DEBUG_PARSER
        std::cout<<"postStr="<<postStr<<" |parsterStr="<<parserStr<<" |topParserStr="<<topParserStr<<" |outParserStr="<<outParserStr<< "  || inToPostFix()6"<<"|i="<<i<<std::endl;
#endif
//debug
                if(!operStk.empty() ){
                    parserStr=operStk.oPop().codeToStr();
//debug
#ifdef DEBUG_PARSER
        std::cout<<"postStr="<<postStr<<" |parsterStr="<<parserStr<<" |topParserStr="<<topParserStr<<" |outParserStr="<<outParserStr<< "  || inToPostFix()7"<<"|i="<<i<<std::endl;
#endif
//debug
                }
                else{
//debug
#ifdef DEBUG_PARSER
        std::cout<<"postStr="<<postStr<<" |parsterStr="<<parserStr<<" |topParserStr="<<topParserStr<<" |outParserStr="<<outParserStr<< "  || inToPostFix()7"<<"|i="<<i<<std::endl;
#endif
//debug
                    return ParsingSpace::ERR_ITP_;
                }
            }
        }
        else if(	parserStr=="*" || parserStr=="/" || parserStr==":"|| parserStr==";"|| parserStr=="^"	){
//debug
#ifdef DEBUG_PARSER
        std::cout<<"postStr="<<postStr<<" |parsterStr="<<parserStr<<" |topParserStr="<<topParserStr<<" |outParserStr="<<outParserStr<< "  || inToPostFix()8"<<"|i="<<i<<std::endl;
#endif
//debug
            if(!operStk.empty() ){
                topParserStr=operStk.oTop().codeToStr();
                while(	noEnterable(parserStr, topParserStr)	){
                    outParserStr=operStk.oPop().codeToStr();
                    //DEBUG
                    postStr+="{";
                    postStr+=outParserStr;
                    postStr+="}";
                    //
                    postFix.operStrEnqueue(outParserStr);
                    if( !operStk.empty() ){
                        topParserStr=operStk.oTop().codeToStr();
//debug
#ifdef DEBUG_PARSER
        std::cout<<"postStr="<<postStr<<" |parsterStr="<<parserStr<<" |topParserStr="<<topParserStr<<" |outParserStr="<<outParserStr<< "  || inToPostFix()9"<<"|i="<<i<<std::endl;
#endif
//debug
                    }
                    else{
                        break;
                    }
                }
            }
            operStk.operStrPush(parserStr);
//debug
#ifdef DEBUG_PARSER
        std::cout<<"postStr="<<postStr<<" |parsterStr="<<parserStr<<" |topParserStr="<<topParserStr<<" |outParserStr="<<outParserStr<< "  || inToPostFix()10"<<"|i="<<i<<std::endl;
#endif
//debug
        }
        else if( parserStr=="+"||parserStr=="-"){
            if(i==0){
                parserStr+=':';
                operStk.operStrPush(parserStr);
            }
            else if(manipStr[i-1]==')'||isdigit(manipStr[i-1])||manipStr[i-1]=='\''||manipStr[i-1]=='\"'){
                if( !operStk.empty() ){
                    topParserStr=operStk.oTop().codeToStr();
                    while( noEnterable(parserStr, topParserStr) ){
                        outParserStr=operStk.oPop().codeToStr();
                        //DEBUG
                        postStr+="{";
                        postStr+=outParserStr;
                        postStr+="}";
                        //
                        postFix.operStrEnqueue(outParserStr);
//debug
#ifdef DEBUG_PARSER
        std::cout<<"postStr="<<postStr<<" |parsterStr="<<parserStr<<" |topParserStr="<<topParserStr<<" |outParserStr="<<outParserStr<< "  || inToPostFix()11"<<"|i="<<i<<std::endl;
#endif
//debug
                        if( !operStk.empty() ){
                            topParserStr=operStk.oTop().codeToStr();
                        }
                        else{
                            break;
                        }
                    }
                }
                operStk.operStrPush(parserStr);
            }
            else{
                parserStr+=':';
                if( !operStk.empty() ){
                    topParserStr=operStk.oTop().codeToStr();
                    while( noEnterable(parserStr, topParserStr) ){
                        outParserStr=operStk.oPop().codeToStr();
                        //DEBUG
                        postStr+="{";
                        postStr+=outParserStr;
                        postStr+="}";
                        //
                        postFix.operStrEnqueue(outParserStr);
                        if( !operStk.empty() ){
                            topParserStr=operStk.oTop().codeToStr();
                        }
                        else{
                            break;
                        }
                    }
                }
                operStk.operStrPush(parserStr);
            }
        }
        else if( isalpha(parserStr[0]) ){
            while( i<=manipStr.length()-2 ){
                if( isalpha(manipStr[i+1])||manipStr[i+1]=='\''||manipStr[i+1]=='\"' ){
                    parserStr+=manipStr[++i];
                    if( (tempCode=serviceOperData.operToCode(parserStr) )==ParsingSpace::NOTDEF_){//sinh conh tanh 
                        continue;
                    }
                    else if(tempCode==ParsingSpace::SIN_||tempCode==ParsingSpace::COS_||tempCode==ParsingSpace::TAN_){
                        if( i<=manipStr.length()-2 && (manipStr[i+1]=='h'||manipStr[i+1]=='H') )
                            parserStr+=manipStr[++i];
                        break;
                    }
                    else{
                        break;
                    }
                }
                else{
//debug
#ifdef DEBUG_PARSER
        std::cout<<"postStr="<<postStr<<" |parsterStr="<<parserStr<<" |topParserStr="<<topParserStr<<" |outParserStr="<<outParserStr<< "  || inToPostFix()12"<<"|i="<<i<<std::endl;
#endif
//debug
                    break;
                }
            }
            serviceOperData.setOperatorData(parserStr);
            if( serviceOperData.getPreInPost()==ParsingSpace::MACONST_ ){
                switch(serviceOperData.getIdenCode() ){
                    case ParsingSpace::PI__:
                        parserStr.resize(0);//Empty();
                        parserStr="3.14159265358979323846";
                        //DEBUG
                        postStr+="{";
                        postStr+="PI'";
                        postStr+="}";
                        //
                        postFix.valStrEnqueue(parserStr);
                        break;
                    case ParsingSpace::ANS__:
                        parserStr.resize(0);//Empty();
                        std::stringstream ss;
                        ss << exAnswer;
                        parserStr=ss.str();
                        //parserStr.format("%lf",exAnswer);
                        //DEBUG
                        postStr+="{";
                        postStr+="Ans'";
                        postStr+="}";
                        //
                        postFix.valStrEnqueue(parserStr);
                        break;
                }
            }
            else if( serviceOperData.getPreInPost()==ParsingSpace::VAR_ ){
                        std::stringstream ss;
                switch(serviceOperData.getIdenCode() ){
                    case ParsingSpace::X__:
                        parserStr.resize(0);//Empty();
                        ss << m_x;
                        parserStr=ss.str();
                        //DEBUG
                        postStr+="{";
                        postStr+="x\'";
                        postStr+="[";
                        postStr+=parserStr;
                        postStr+="]";
                        postStr+="}";
                        //
                        postFix.valStrEnqueue(parserStr);
                        break;
                    case ParsingSpace::Y__:
                        parserStr.resize(0);//Empty();
                        ss << m_y;
                        parserStr=ss.str();
                        //DEBUG
                        postStr+="{";
                        postStr+="y\'";
                        postStr+="[";
                        postStr+=parserStr;
                        postStr+="]";
                        postStr+="}";
                        //
                        postFix.valStrEnqueue(parserStr);
                        break;
                    case ParsingSpace::Z__:
                        parserStr.resize(0);//Empty();
                        ss << m_z;
                        parserStr=ss.str();
                        //DEBUG
                        postStr+="{";
                        postStr+="z\'";
                        postStr+="[";
                        postStr+=parserStr;
                        postStr+="]";
                        postStr+="}";
                        //
                        postFix.valStrEnqueue(parserStr);
                        break;
                    case ParsingSpace::MU__:
                        parserStr.resize(0);//Empty();
                        ss << m_mu;
                        parserStr=ss.str();
                        //DEBUG
                        postStr+="{";
                        postStr+="mu\'";
                        postStr+="[";
                        postStr+=parserStr;
                        postStr+="]";
                        postStr+="}";
                        //
                        postFix.valStrEnqueue(parserStr);
                        break;
                    case ParsingSpace::VPAR__:
                        parserStr.resize(0);//Empty();
                        ss << m_vpar;
                        parserStr=ss.str();
                        //DEBUG
                        postStr+="{";
                        postStr+="vpar\'";
                        postStr+="[";
                        postStr+=parserStr;
                        postStr+="]";
                        postStr+="}";
                        //
                        postFix.valStrEnqueue(parserStr);
                        break;
                    case ParsingSpace::ANS__:
                        parserStr.resize(0);//Empty();
                        ss << exAnswer;
                        parserStr=ss.str();
                        //parserStr.format("%lf",exAnswer);
                        //DEBUG
                        postStr+="{";
                        postStr+="Ans'";
                        postStr+="[";
                        postStr+=parserStr;
                        postStr+="]";
                        postStr+="}";
                        //
                        postFix.valStrEnqueue(parserStr);
                        break;
                }
            }
            else if(!operStk.empty() ){
                topParserStr=operStk.oTop().codeToStr();
                while( noEnterable(parserStr,topParserStr) ){
                    outParserStr=operStk.oPop().codeToStr();
                    //DEBUG
                    postStr+="{";
                    postStr+=outParserStr;
                    postStr+="}";
                    //
                    postFix.operStrEnqueue(outParserStr);
                    if( !operStk.empty() ){
                        topParserStr=operStk.oTop().codeToStr();
                    }
                    else{
                        break;
                    }
                }
                operStk.operStrPush(parserStr);
            }
            else{
                operStk.operStrPush(parserStr);
            }
        }
        else{//
           while(manipStr.length()>i+1 && ( manipStr[i+1]=='.'|| ('0'<=manipStr[i+1]&&manipStr[i+1]<='9') ) ){
                parserStr+=manipStr[++i];
            }
            //DEBUG
            postStr+="{";
            postStr+=parserStr;
            postStr+="}";
            //
            postFix.valStrEnqueue(parserStr);
        }
    }
    while( !operStk.empty() ){
        outParserStr=operStk.oPop().codeToStr();
        //DEBUG
        postStr+="{";
        postStr+=outParserStr;
        postStr+="}";
        //
        postFix.operStrEnqueue(outParserStr);
    }
    return ParsingSpace::NOERR_;
}

int ParsingCore::postFixEvaluate()
{
    int exprSize;
    UserStack evaStk;
    double		value=0.0;
    double		operand1,operand2;
    OperatorData  op;
    queueData	bufferQData;

    exprSize=postFix.showCnt();
    for(int i=0;i<exprSize;i++){
        bufferQData=postFix.dequeue();
        if(bufferQData.kind){
            evaStk.dPush(atof(bufferQData.sQv.obCstr.c_str()));
        }
        else{
            op=bufferQData.sQv.obOperatorData;
            if(op.getPreInPost()==ParsingSpace::PREFIX_){
                if(!evaStk.empty()){
                    operand1=evaStk.dPop();
                }
                else{
                    return ParsingSpace::ERR_PFEV_;
                }
                value=calculate(op,operand1);
                evaStk.dPush(value);
            }
            else if(op.getPreInPost()==ParsingSpace::INFIX_){
                if(!evaStk.empty()){
                    operand2=evaStk.dPop();
                }
                else{
                    return ParsingSpace::ERR_PFEV_;
                }
                if(!evaStk.empty() ){
                    operand1=evaStk.dPop();
                }
                else{
                    return ParsingSpace::ERR_PFEV_;
                }
                value=calculate(operand1,op,operand2);
                evaStk.dPush(value);
            }
            else if(op.getPreInPost()==ParsingSpace::POSTFIX_){
                if(!evaStk.empty() ){
                    operand2=evaStk.dPop();
                }
                else{
                    return ParsingSpace::ERR_PFEV_;
                }
                value=calculate(operand2,op);
                evaStk.dPush(value);
            }
            else{
                return ParsingSpace::ERR_PFEV_;
            }
        }
    }
    if(!evaStk.empty()){
        value=evaStk.dPop();
    }
    else{
        return ParsingSpace::ERR_PFEV_;
    }
    result=value;
    return ParsingSpace::NOERR_;
}


double ParsingCore::calculate( OperatorData op, double operand1)
{
    double temp;
    if(isRad==true){
        switch(op.getIdenCode() ){
            case ParsingSpace::SIN_:
                temp=sin(operand1);
                break;
            case ParsingSpace::COS_:
                temp=cos(operand1);
                break;
            case ParsingSpace::TAN_:
                temp=tan(operand1);
                break;
            case ParsingSpace::CSC_:
                temp=1.0/sin(operand1);
                break;
            case ParsingSpace::SEC_:
                temp=1.0/cos(operand1);
                break;
            case ParsingSpace::COT_:
                temp=1.0/tan(operand1);
                break;
            case ParsingSpace::ASIN_:
                temp=asin(operand1);
                break;
            case ParsingSpace::ACOS_:
                temp=acos(operand1);
                break;
            case ParsingSpace::ATAN_:
                temp=atan(operand1);
                break;
            case ParsingSpace::ABS_:
                temp=fabs(operand1);
                break;
            case ParsingSpace::LOG_:
                temp=log10(operand1);
                break;
            case ParsingSpace::LN_:
                temp=log(operand1);
                break;
            case ParsingSpace::TENPOWER_:
                temp=pow(10, operand1);
                break;
            case ParsingSpace::EXPF_:
                temp=exp(operand1);
                break;
            case ParsingSpace::SQRT_:
                temp=sqrt(operand1);
                break;
            case ParsingSpace::SINH_:
                temp=sinh(operand1);
                break;
            case ParsingSpace::COSH_:
                temp=cosh(operand1);
                break;
            case ParsingSpace::TANH_:
                temp=tanh(operand1);
                break;
            case ParsingSpace::PLUS_:
                temp=operand1;
                break;
            case ParsingSpace::MINUS_:
                temp= -operand1;
                break;
            case ParsingSpace::EP_:
                temp=pow(10,operand1);
                break;
            default:
                temp=operand1;
                break;
        }
    }
    else{//degree
        switch(op.getIdenCode() ){
            case ParsingSpace::SIN_:
                operand1=operand1/180.0*3.14159265358979323846;
                temp=sin(operand1);
                break;
            case ParsingSpace::COS_:
                operand1=operand1/180.0*3.14159265358979323846;
                temp=cos(operand1);
                break;
            case ParsingSpace::TAN_:
                operand1=operand1/180.0*3.14159265358979323846;
                temp=tan(operand1);
                break;
            case ParsingSpace::CSC_:
                operand1=operand1/180.0*3.14159265358979323846;
                temp=1.0/sin(operand1);
                break;
            case ParsingSpace::SEC_:
                operand1=operand1/180.0*3.14159265358979323846;
                temp=1.0/cos(operand1);
                break;
            case ParsingSpace::COT_:
                operand1=operand1/180.0*3.14159265358979323846;
                temp=1.0/tan(operand1);
                break;
            case ParsingSpace::ASIN_:
                temp=asin(operand1);
                temp=temp*180.0/3.14159265358979323846;
                break;
            case ParsingSpace::ACOS_:
                temp=acos(operand1);
                temp=temp*180.0/3.14159265358979323846;
                break;
            case ParsingSpace::ATAN_:
                temp=atan(operand1);
                temp=temp*180.0/3.14159265358979323846;
                break;
            case ParsingSpace::ABS_:
                temp=fabs(operand1);
                break;
            case ParsingSpace::LOG_:
                temp=log10(operand1);
                break;
            case ParsingSpace::LN_:
                temp=log(operand1);
                break;
            case ParsingSpace::TENPOWER_:
                temp=pow(10, operand1);
                break;
            case ParsingSpace::EXPF_:
                temp=exp(operand1);
                break;
            case ParsingSpace::SQRT_:
                temp=sqrt(operand1);
                break;
            case ParsingSpace::SINH_:
                temp=sinh(operand1);
                break;
            case ParsingSpace::COSH_:
                temp=cosh(operand1);
                break;
            case ParsingSpace::TANH_:
                temp=tanh(operand1);
                break;
            case ParsingSpace::PLUS_:
                temp=operand1;
                break;
            case ParsingSpace::MINUS_:
                temp= -operand1;
                break;
            case ParsingSpace::EP_:
                temp=pow(10,operand1);
                break;
            default:
                temp=operand1;
                break;
        }
    }

    return temp;
}

double ParsingCore::calculate(double operand1, OperatorData op, double operand2)
{
    double temp=0.0;
    switch(op.getIdenCode() ){
        case ParsingSpace::ADD_:
            temp=operand1+operand2;
            break;
        case ParsingSpace::SUBTRACT_:
            temp=operand1-operand2;
            break;
        case ParsingSpace::MULTI_:
        case ParsingSpace::PMHM_:
        case ParsingSpace::PPHM_:
            temp=operand1*operand2;
            break;
        case ParsingSpace::DIVIDE_:
            temp=operand1/operand2;
            break;
        case ParsingSpace::POWER_:
            temp=pow(operand1,operand2);
            break;
        case ParsingSpace::NSQRT_:
            temp=pow(operand2,1.0/operand2);
            break;
        case ParsingSpace::EI_:
            temp=operand1*pow(10,operand2);
            break;
        default:
            //exit(1);
            ;
    }
    return temp;
}
double ParsingCore::calculate(double operand2, OperatorData op)
{
    double temp=0.0;
    switch(op.getIdenCode() ){
        case ParsingSpace::INVERSE_:
            temp=pow(operand2,-1);
            break;
        case ParsingSpace::SQUARE_:
            temp=pow(operand2,2);
            break;
        default:
            //		exit(1)
            ;
    }
    return temp;
}



