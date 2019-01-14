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
class SteepestDescent{

private:

    //Data members
    int itMax; //maximum iterations
    double ls; //line search increment
    double tol; //convergence tolerance

public:

    //Constructors
    SteepestDescent();
    SteepestDescent(int iterationLimit, double lineSearchInc, double convergenceTolerance);

    //Optimisation function
    VecF<int> operator() (M& model, VecF<double>& x);
};


//Steepest descent algorithm with Armijo Search
template <typename M>
class SteepestDescentArmijo{

private:

    //Data members
    int itMax; //maximum iterations
    double tau; //armijo factor
    double tol; //convergence tolerance

public:

    //Constructors
    SteepestDescentArmijo();
    SteepestDescentArmijo(int iterationLimit, double lineSearchInc, double convergenceTolerance);

    //Optimisation function
    VecF<int> operator() (M& model, VecF<double>& x);
};


//Model with function and derivative
class Model{

public:

    virtual double function(VecF<double>& x)=0;
    virtual VecF<double> gradient(VecF<double>& x)=0;

};

#include "opt.tpp"

#endif //NL_OPT_H
