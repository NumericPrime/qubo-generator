#ifndef _H_QBVARS_INCLUDED
#define _H_QBVARS_INCLUDED
#include<string.h>
class QBVar{
    public:
    unsigned len;
    double * weights;
    int * bits;
    unsigned char delarray=0;
    QBVar(double * _weights,int * _bits,unsigned len,unsigned char delar=0){
        //std::cout<<"constr2"<<std::endl;
        this->len=len;
        this->delarray=delar;
        this->weights=_weights;
        this->bits=_bits;
        //std::cout<<"constre"<<std::endl;
    }
    ~QBVar(){
        if(delarray&1) delete weights;
        if(delarray&2) delete bits;
    }
};
QBVar& operator+(QBVar& v1,QBVar& v2){
    double * nweights=new double[v2.len+v1.len];
    int * nbits=new int[v2.len+v1.len];
    memcpy(nweights,v1.weights,sizeof(double)*v1.len);
    memcpy(nweights+v1.len,v2.weights,sizeof(double)*v2.len);
    memcpy(nbits,v1.bits,sizeof(double)*v1.len);
    memcpy(nbits+v1.len,v2.bits,sizeof(double)*v2.len);
    std::cout<<"start"<<nweights<<" "<<nbits<<" "<<v2.len+v1.len<<std::endl;
    //QBVar rt2(nweights,nbits,v2.len+v1.len,3);
    //std::cout<<"end"<<std::endl;
    QBVar * rt=new QBVar(nweights,nbits,v2.len+v1.len,3);
    std::cout<<"end"<<std::endl;
    return *rt;
}
QBVar& operator*(double fct,QBVar& v2){
    double * nweights=new double[v2.len];
    for(unsigned i=0;i<v2.len;i++) nweights[i]=v2.weights[i]*fct;
    QBVar& ret=*(new QBVar(nweights,v2.bits,v2.len,1|v2.delarray));

    return ret;
}
QBVar& operator*(QBVar& v2,double fct){
    double * nweights=new double[v2.len];
    for(unsigned i=0;i<v2.len;i++) nweights[i]=v2.weights[i]*fct;
    QBVar& ret=*(new QBVar(nweights,v2.bits,v2.len,1|v2.delarray));
    return ret;
}
#endif