#include "opt.h"

//##### STEEPEST DESCENT #####

//Default constructor
template <typename M>
SteepestDescent<M>::SteepestDescent() {
    itMax=1000;
    ls=1e-3;
    tol=1e-6;
}

//Constructor with algorithm search parameters
template <typename M>
SteepestDescent<M>::SteepestDescent(int iterationLimit, double lineSearchInc, double convergenceTolerance) {
    itMax=iterationLimit;
    ls=lineSearchInc;
    tol=convergenceTolerance;
}

//Optimisation Algorithm
template <typename M>
VecF<int> SteepestDescent<M>::operator()(M& model, VecF<double> &x) {

    //Initialise optimisation variables
    int it=0; //number of iterations
    double f0, f1, f; //previous, current and generic function value
    double df; //difference between current and previous function evaluations
    VecF<double> g(x); //derivative of function
    VecF<double> d(x); //search direction
    VecF<int> status(2); //optimisation status: convergence status, iterations
    f0=numeric_limits<double>::infinity();

    //Evaluate gradient and check non-zero before commencing descent
    g=model.gradient(x);
    if(vAsum(g)<tol){
        status[0]=1;
        status[1]=0;
        return status;
    }

    //Steepest descent algorithm
    for(int i=0; i<itMax; ++i){

        //line search
        d=-g*ls;
        f1=model.function(x);
        for(;;){
            x+=d;
            f=model.function(x);
            if(f>f1){//passed through minimum, backtrack
                x-=d;
                ++it;
                break;
            }
            else f1=f;
        }

        //check convergence
        df=abs(f1-f0);
        if(abs(f1)<tol || df<tol) break;
        else f0=f1;
        //recalculate gradient
        g=model.gradient(x);
    }

    if(it==itMax-1) status[0]=2;
    else status[0]=0;
    status[1]=it;

    return status;
}


//##### STEEPEST DESCENT WITH ARMIJO BACKTRACKING LINE SEARCH #####

//Default constructor
template <typename M>
SteepestDescentArmijo<M>::SteepestDescentArmijo() {
    itMax=1000;
    tau=0.5;
    tol=1e-6;
}

//Constructor with algorithm search parameters
template <typename M>
SteepestDescentArmijo<M>::SteepestDescentArmijo(int iterationLimit, double lineSearchInc, double convergenceTolerance) {
    itMax=iterationLimit;
    tau=lineSearchInc;
    tol=convergenceTolerance;
}

//Optimisation Algorithm
template <typename M>
VecF<int> SteepestDescentArmijo<M>::operator()(M& model, VecF<double> &x) {

    //Initialise optimisation variables
    int it=0; //number of iterations
    double f0, f1, f; //previous, current and generic function value
    double df; //difference between current and previous function evaluations
    double alpha; //overall line search step size
    double gSq; //modulus square of gradient
    VecF<double> g(x); //derivative of function
    VecF<double> d(x); //coordinates after search direction
    VecF<int> status(2); //optimisation status: convergence status, iterations
    f0=numeric_limits<double>::infinity();

    //Evaluate gradient and check non-zero before commencing descent
    g=model.gradient(x);
    if(vAsum(g)<tol){
        status[0]=1;
        status[1]=0;
        return status;
    }

    //Steepest descent algorithm
    for(int i=0; i<itMax; ++i){

        //backtracking line search
        alpha=1.0;
        gSq=vNormSq(g);
        f=model.function(x);
        for(;;){
            d=x-g*alpha;
            f1=model.function(d);
            f1+=0.5*alpha*gSq;
            if(f1<=f){
                x=d;
                ++it;
                break;
            }
            else alpha*=tau;
        }

        //check convergence
        df=abs(f1-f0);
        if(abs(f1)<tol || df<tol) break;
        else f0=f1;
        //recalculate gradient
        g=model.gradient(x);
    }

    if(it==itMax-1) status[0]=2;
    else status[0]=0;
    status[1]=it;

    return status;
}