#ifndef _H_QBBASETASKS_INCLUDED
#define _H_QBBASETASKS_INCLUDED
#include "QBMatrix.hpp"
#include "QBExpression.hpp"
#include<limits>
#include<cmath>
class QBResolveTask{
    public:
        void * qbmatr;
        QBResolveTask(){
        }
        void putPtr(void *qb){
            this->qbmatr=qb;
        }
};
class QBIntegrateTask{
    public:
        void * qbmatr;
        QBIntegrateTask(){
        }
        void putPtr(void *qb){
            this->qbmatr=qb;
        }
};
QBIntegrateTask& operator<<(QBIntegrateTask& integrator,QBExpression& qbe){
    QBMatrix * qbm=(QBMatrix*)integrator.qbmatr;
    for(int i=0;i<qbe.b1.size();i++) qbm->add(qbe.b1[i],qbe.b2[i],qbe.coef12[i]);
    for(int i=0;i<qbe.bt.size();i++) qbm->qadd(qbe.bt[i],qbe.bt[i],qbe.coeflin[i]);
    return integrator;
}
void runIntegrator(QBMatrix& qbm,QBExpression& qbe){
    for(int i=0;i<qbe.b1.size();i++) qbm.add(qbe.b1[i],qbe.b2[i],qbe.coef12[i]);
    for(int i=0;i<qbe.bt.size();i++) qbm.qadd(qbe.bt[i],qbe.bt[i],qbe.coeflin[i]);
}
QBExpression&  runResolve(QBMatrix& qbm,QBExpression& qbe){
    QBExpression& ret=*(new QBExpression(qbe.coeflin,qbe.bt));
    ret.const_expr=qbe.const_expr;
    for(int i=0;i<qbe.b1.size();i++) {
        /*unsigned introduced_bit=qbm->popUnassigned();
        qbm->carry(qbe.b1[i],qbe.b2[i],introduced_bit);*/
        ret.bt.push_back(qbm.pcarry(qbe.b1[i],qbe.b2[i]));
        ret.coeflin.push_back(qbe.coef12[i]);
    }
    return ret;
}
QBExpression& createPacked(QBMatrix& qbm,double start,unsigned len,bool invertStart=false,double cexpr=0){
    QBExpression& ret=*(new QBExpression());
    if(invertStart) ret.pushLin(qbm.popUnassigned(),start*(invertStart?-1:1));
    for(unsigned i=1;i<len;i++) ret.pushLin(qbm.popUnassigned(),start/=2.);
    ret.const_expr=cexpr;
    return ret;
}
QBExpression& createRevPacked(QBMatrix& qbm,double end,unsigned len,bool invertStart=false,double cexpr=0){
    QBExpression& ret=*(new QBExpression());
    if(len>1) ret.pushLin(qbm.popUnassigned(),end);
    for(unsigned i=1;i<len-1;i++) ret.pushLin(qbm.popUnassigned(),end*=2.);
    if(invertStart) ret.pushLin(qbm.popUnassigned(),end*2.*(invertStart?-1:1));
    ret.const_expr=cexpr;
    return ret;
}
#define FLIP_FIRST 1
#define FLIP_ALL 2
#define EXTRAPOLATE 4
#define INTVAR 4
QBExpression&  createPacked(QBMatrix& qbm,QBExpression& qbe,unsigned char autoDeg=0){
    //ret.const_expr=qbe.const_expr;
    double posSum=0;
    double negSum=0;
    double smallestab= std::numeric_limits<double>::infinity();
    for(size_t i=0;i<qbe.coeflin.size();i++) {
        if(qbe.coeflin[i]>0) posSum+=qbe.coeflin[i];
        else negSum-=qbe.coeflin[i];
        if(smallestab>((qbe.coeflin[i]>0)?qbe.coeflin[i]:-qbe.coeflin[i])) 
        smallestab=((qbe.coeflin[i]>0)?qbe.coeflin[i]:-qbe.coeflin[i]);
    }
    for(size_t i=0;i<qbe.b1.size();i++) {
        if(qbe.coef12[i]>0) posSum+=qbe.coef12[i];
        else negSum-=qbe.coeflin[i];
        if(smallestab>((qbe.coef12[i]>0)?qbe.coef12[i]:-qbe.coef12[i])) 
        smallestab=((qbe.coef12[i]>0)?qbe.coef12[i]:-qbe.coef12[i]);
    }
    double makcoef=posSum+negSum;
    unsigned projectedLength=(unsigned)round(log2(makcoef/smallestab))+1;
    if(autoDeg==0) return createRevPacked(qbm,smallestab,projectedLength,false,qbe.const_expr);
    if(autoDeg==1) return createRevPacked(qbm,smallestab,projectedLength+1,true,qbe.const_expr);
    if(autoDeg==2) return createRevPacked(qbm,-smallestab,projectedLength,false,qbe.const_expr);
    if(autoDeg==3) return createRevPacked(qbm,-smallestab,projectedLength+1,true,qbe.const_expr);
    if(autoDeg==4){
        if(posSum>0&&negSum==0) return createRevPacked(qbm,smallestab,projectedLength,false,qbe.const_expr);
        if(posSum&&negSum<0) return createRevPacked(qbm,-smallestab,projectedLength,false,qbe.const_expr);
        if(autoDeg&2){
            return createRevPacked(qbm,-smallestab,projectedLength+1,true,qbe.const_expr);
        }else{
            return createRevPacked(qbm,-smallestab,projectedLength+1,false,qbe.const_expr);
        }
    }
    //return NULL;
}
QBExpression&  runPacked(QBMatrix& qbm,QBExpression& qbe,unsigned char autoDeg=4){
    QBExpression& ret=createPacked(qbm,qbe,autoDeg);
    QBExpression& diff=qbe-ret;
    runIntegrator(qbm,diff*diff);
    return ret;
}
QBExpression& operator<<(QBResolveTask& integrator,QBExpression& qbe){
    QBMatrix * qbm=(QBMatrix*)integrator.qbmatr;
    QBExpression& ret=*(new QBExpression(qbe.coeflin,qbe.bt));
    ret.const_expr=qbe.const_expr;
    for(int i=0;i<qbe.b1.size();i++) {
        /*unsigned introduced_bit=qbm->popUnassigned();
        qbm->carry(qbe.b1[i],qbe.b2[i],introduced_bit);*/
        ret.bt.push_back(qbm->pcarry(qbe.b1[i],qbe.b2[i]));
        ret.coeflin.push_back(qbe.coef12[i]);
    }
    return ret;
}
QBExpression& makeVar(QBMatrix& qm,std::vector<double> w){
    std::cout<<w.size()<<std::endl;
    QBExpression& rt=*(new QBExpression(w,qm.popUnassigned(w.size())));
    rt.delinfo=1;
    return rt;
}
std::vector<double>& makeWeights(unsigned size,double sz,unsigned char mode=0){
    std::vector<double>& ret=*(new std::vector<double>(size));
    if(mode&INTVAR){
        unsigned val=1;
        //Break condition because of potential overflow
        for(unsigned i=size-1;i<size;i--) {
            ret[i]=val;
            val<<=1;
        }
    }else{
        for(unsigned i=0;i<size;i++) {
            ret[i]=sz;
            sz/=2.;
        }
    }
    if(mode&FLIP_ALL) for(unsigned i=0;i<size;i++) {
            ret[i]*=-1.;
    }
    if(mode&FLIP_FIRST) ret[0]*=-1.;
    return ret;
}
QBExpression& makeVar(QBMatrix& qm,unsigned size,unsigned char mode=0){
    return makeVar(qm,makeWeights(size,1.,mode));
}
QBExpression& makeVar(QBMatrix& qm,unsigned size,double sz,unsigned char mode=0){
    return makeVar(qm,makeWeights(size,sz,mode));
}
#endif