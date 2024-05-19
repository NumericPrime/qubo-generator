#ifndef _H_QEXPRESSIONS_INCLUDED
#define _H_QEXPRESSIONS_INCLUDED
//#include"QBVars.hpp"
#include<vector>
class QBExpression{
    public:
    double const_expr;
    std::vector<unsigned> b1;
    std::vector<unsigned> b2;
    std::vector<double> coef12;
    std::vector<unsigned> bt;
    std::vector<double> coeflin;
    unsigned char delinfo;
    QBExpression(unsigned char _delinfo=0){
        this->delinfo=_delinfo;
        const_expr=0;
    }/*
    QBExpression(QBVar& qvar){
        for(std::size_t i=0;i<qvar.len;i++){
            bt.push_back(qvar.bits[i]);
            coeflin.push_back(qvar.weights[i]);
        }
        const_expr=0;
    }*/
    QBExpression(std::vector<double> weights,std::vector<unsigned> bits){
        std::cout<<"Create obj";
        bt=bits;
        coeflin=weights;
        const_expr=0;
    }
    QBExpression(double * _weights,unsigned * bits,unsigned len){
        bt.insert(bt.end(),bits,bits+len);
        coeflin.insert(coeflin.end(),_weights,_weights+len);
        const_expr=0;
    }
    ~QBExpression(){
        if(!delinfo){
        bt.~vector<unsigned>();
        b1.~vector<unsigned>();
        b2.~vector<unsigned>();
        coef12.~vector<double>();
        coeflin.~vector<double>();}
    }
    void pushDual(unsigned _b1,unsigned _b2,double coef){
        if(_b1==_b2) return pushLin(_b1,coef);
        for(std::size_t i=0;i<b1.size();i++) {
            if(b1[i]==_b1&&b2[i]==_b2){
                coef12[i]+=coef;
                return;
            }
            if(b1[i]==_b2&&b2[i]==_b1){
                coef12[i]+=coef;
                return;
            }
        }
        b1.push_back(_b1);
        b2.push_back(_b2);
        coef12.push_back(coef);
    }
    void pushLin(int _b1,double coef){
        for(std::size_t i=0;i<bt.size();i++) {
            if(bt[i]==_b1){
                coeflin[i]+=coef;
                return;
            }
        }
        coeflin.push_back(coef);
        bt.push_back(_b1);
    }
    QBExpression * copy(){
        QBExpression * ret=new QBExpression();
        for(std::size_t i=0;i<b1.size();i++){
            ret->b1.push_back(b1[i]);
            ret->b2.push_back(b1[i]);
            ret->coef12.push_back(coef12[i]);
        }
        ret->const_expr=const_expr;
        for(std::size_t i=0;i<bt.size();i++){
            ret->bt.push_back(bt[i]);
            ret->coeflin.push_back(coeflin[i]);
        }
        return ret;
    }
};
QBExpression& operator+(QBExpression& s1,QBExpression& s2){
    QBExpression * rt=s1.copy();
    QBExpression& rf=*rt;
    for(std::size_t i=0;i<s2.b1.size();i++){
        rf.pushDual(s2.b1[i],s2.b2[i],s2.coef12[i]);
    }
    for(std::size_t i=0;i<s2.bt.size();i++){
        rf.pushLin(s2.bt[i],s2.coeflin[i]);
    }
    rf.const_expr+=s2.const_expr;
    return rf;
}
QBExpression& operator*(QBExpression& s1,QBExpression& s2){
    QBExpression& rf=*(new QBExpression());
    for(std::size_t i=0;i<s1.bt.size();i++){
        for(std::size_t j=0;j<s2.bt.size();j++){
            rf.pushDual(s1.bt[i],s2.bt[j],s1.coeflin[i]*s2.coeflin[j]);
        }
    }
    for(std::size_t i=0;i<s2.bt.size();i++){
        rf.pushLin(s2.bt[i],s2.coeflin[i]*s1.const_expr);
    }
    for(std::size_t i=0;i<s1.bt.size();i++){
        rf.pushLin(s1.bt[i],s1.coeflin[i]*s2.const_expr);
    }
    rf.const_expr=s1.const_expr*s2.const_expr;
    return rf;
}
QBExpression& operator*(double mul,QBExpression& exp){
        QBExpression * ret=new QBExpression();
        for(std::size_t i=0;i<exp.b1.size();i++){
            ret->b1.push_back(exp.b1[i]);
            ret->b2.push_back(exp.b1[i]);
            ret->coef12.push_back(exp.coef12[i]*mul);
        }
        ret->const_expr=exp.const_expr*mul;
        for(std::size_t i=0;i<exp.bt.size();i++){
            ret->bt.push_back(exp.bt[i]);
            ret->coeflin.push_back(exp.coeflin[i]*mul);
        }
        return *ret;
    }
QBExpression& operator*(QBExpression& exp,double mul){
    return mul*exp;
}
QBExpression& operator+(double s1,QBExpression& s2){
    QBExpression& rf=*(s2.copy());
    rf.const_expr+=s1;
    return rf;
}
QBExpression& operator+(QBExpression& exp,double mul){
    return mul+exp;
}
QBExpression& operator-(QBExpression& exp,double mul){
    return (-mul)+exp;
}
QBExpression& operator-(double mul,QBExpression& exp){
    return (-mul)+exp;
}
QBExpression& operator-(QBExpression& s1,QBExpression& s2){
    return s1+(-1*s2);
}
std::ostream& operator<<(std::ostream& ostr,QBExpression& exp){
    ostr<<"Expr : "<<exp.const_expr<<std::endl;
    for(unsigned i=0;i<exp.coeflin.size();i++) ostr<<exp.bt[i]<<" : "<<exp.coeflin[i]<<std::endl;
    for(unsigned i=0;i<exp.coef12.size();i++) ostr<<exp.b1[i]<<" "<<exp.b2[i]<<" : "<<exp.coef12[i]<<std::endl;
    ostr<<"VAR DONE "<<std::endl;
    return ostr;
}
/*
QBExpression& operator+(double s1,QBVar& s2){
    QBExpression& rf=*(new QBExpression(s2));
    rf.const_expr+=s1;
    return rf;
}
QBExpression& operator+(QBVar& exp,double mul){
    return mul+exp;
}
QBExpression& operator*(QBVar& f1,QBVar& f2){
    QBExpression& rf1=*(new QBExpression(f1));
    QBExpression& rf2=*(new QBExpression(f2));
    return rf1*rf2;
}*/
#endif