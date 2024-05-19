#ifndef _H_QBMATRIX_INCLUDED
#define _H_QBMATRIX_INCLUDED
#include <iostream>
#include <vector>
#include<string.h>
#include "QMisc.hpp"
//using namespace std;
class QBMatrix {
    public:
        unsigned sz;
        double ** data;
        unsigned * unassignedBase;
        unsigned * currentUnassigned;
        double scale=1;
        double secscale=1;
        double scaleStep=10;
        std::unordered_map<std::pair<unsigned,unsigned>,unsigned> * um=NULL;
        QBMatrix(unsigned sz){
            this->sz=sz;
            data = new double*[sz];
            for(unsigned i=0;i<sz;i++) {
                data[i]=new double[sz];
                memset(data[i],0,sz*sizeof(double));
            }
            unassignedBase=new unsigned[sz];
            currentUnassigned=unassignedBase;
            for(unsigned i=0;i<sz;i++) unassignedBase[i]=i;
        }
        QBMatrix& track_carries(){
            um=new std::unordered_map<std::pair<unsigned,unsigned>,unsigned> ();
            return *this;
        }
        ~QBMatrix(){
            for(unsigned i=0;i<sz;i++) delete data[i];
            delete data;
            delete unassignedBase;
            if(um) delete um;
        }
        QBMatrix * add(unsigned i,unsigned j,double val);
        QBMatrix * qadd(unsigned i,unsigned j,double val);
        QBMatrix * carry(unsigned i,unsigned j,unsigned t);
        unsigned pcarry(unsigned i,unsigned j);
        QBMatrix * addquad(double constant_,std::vector<double>& coefs,std::vector<int>& vars);
        unsigned popUnassigned();
        QBMatrix * pushUnassigned(unsigned a);
        std::vector<unsigned>& popUnassigned(unsigned n);
        
};
QBMatrix * QBMatrix::add(unsigned i,unsigned j,double val) {
    if(i<=j) return qadd(i,j,val);
    else return qadd(j,i,val);
}
QBMatrix * QBMatrix::qadd(unsigned i,unsigned j,double val) {
    data[i][j]+=val*scale*secscale;
    return this;
}
QBMatrix * QBMatrix::carry(unsigned i,unsigned j,unsigned t) {
    qadd(t,t,3);
    add(i,t,-2);
    add(j,t,-2);
    return add(i,j,1);
}
unsigned QBMatrix::pcarry(unsigned i,unsigned j) {
    unsigned buff;
    if(!um) {
        buff=popUnassigned();
        carry(i,j,buff);
        return buff;
    }
    if(i>j){
        buff=i;
        i=j;
        j=buff;
    }
    std::pair<unsigned,unsigned> pa(i,j);
    if(um->count(pa)){
        return um->at(pa);
    }else{
        buff=popUnassigned();
        carry(i,j,buff);
        um->insert({pa,buff});
        return buff;
    }
}
unsigned QBMatrix::popUnassigned(){
    return *(currentUnassigned++);
}
QBMatrix * QBMatrix::pushUnassigned(unsigned val){
    *(--currentUnassigned)=val;
    return this;
}
QBMatrix * QBMatrix::addquad(double constant_,std::vector<double>& coefs,std::vector<int>& vars){
    unsigned len=coefs.size();
    for(unsigned i=0;i<len;i++) for(unsigned j=i+1;j<len;j++){
        add(vars[i],vars[j],2*coefs[i]*coefs[j]);
    }
    for(unsigned i=0;i<len;i++) {
        qadd(vars[i],vars[i],coefs[i]*coefs[i]+2.*constant_*coefs[i]);
    }
    return this;
}
std::ostream& operator<<(std::ostream& stream,QBMatrix& qbm){
    for(unsigned i=0;i<qbm.sz;i++) {
        for(unsigned j=0;j<qbm.sz;j++) 
        stream<<qbm.data[i][j]<<" ";
        stream<<std::endl;
        }
    return stream;
}
std::vector<unsigned>& QBMatrix::popUnassigned(unsigned n) {
    std::vector<unsigned>& ret=*(new std::vector<unsigned>(n));
    for(unsigned i=0;i<n;i++) ret[i]=popUnassigned();
    return ret;
}
#endif