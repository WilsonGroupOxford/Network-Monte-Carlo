//Optimisation Algorithms
//Takes model to evaluate function and derivative
#ifndef NL_OPT_H
#define NL_OPT_H

#include <iostream>
#include <limits>
#include "vecf.h"
#include "vec_func.h"

using namespace std;


//Steepest descent algorithm
template <typename M>
class SteepestDescentMultiDim{

private:

    //Data members
    int itMax; //maximum iterations
    double ls; //line search increment
    double tol; //convergence tolerance

public:

    //Constructors
    SteepestDescentMultiDim();
    SteepestDescentMultiDim(int iterationLimit, double lineSearchInc, double convergenceTolerance);

    //Optimisation function
    VecF<int> operator() (M& model, VecF<double>& x);
};


//Steepest descent algorithm with Armijo Search
template <typename M>
class SteepestDescentArmijoMultiDim{

private:

    //Data members
    int itMax; //maximum iterations
    double tau; //armijo factor
    double tol; //convergence tolerance

public:

    //Constructors
    SteepestDescentArmijoMultiDim();
    SteepestDescentArmijoMultiDim(int iterationLimit, double lineSearchInc, double convergenceTolerance);

    //Optimisation function
    VecF<int> operator() (M& model, VecF<double>& x);
};

//Newton root finder
template <typename M>
class Newton{

private:

    //Data members
    int itMax; //maximum iterations
    long double tol; //convergence tolerance

public:

    //Constructors
    Newton();
    Newton(int iterationLimit, long double convergenceTolerance);

    //Root finding
    VecF<int> operator() (M& model, long double& x);
};

//Halley root finder
template <typename M>
class Halley{

private:

    //Data members
    int itMax; //maximum iterations
    long double tol; //convergence tolerance

public:

    //Constructors
    Halley();
    Halley(int iterationLimit, long double convergenceTolerance);

    //Root finding
    VecF<int> operator() (M& model, long double& x);
};

//One dimensional model with function and derivative
class FuncGrad{

public:

    virtual long double function(long double& x)=0;
    virtual long double gradient(long double& x)=0;

};

//One dimensional model with function, derivative and second derivative
class FuncGradHess{

public:

    virtual long double function(long double& x)=0;
    virtual long double gradient(long double& x)=0;
    virtual long double hessian(long double& x)=0;

};

//Multidimensional model with function and derivative
class FuncGradMultiDim{

public:

    virtual double function(VecF<double>& x)=0;
    virtual VecF<double> gradient(VecF<double>& x)=0;

};


#include "opt.tpp"

#endif //NL_OPT_H
