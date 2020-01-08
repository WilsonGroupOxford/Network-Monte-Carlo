//1D Potential Models Compatible with Optimistaion Routines
#ifndef NL_POT1D_H
#define NL_POT1D_H

#include <iostream>
#include "pot.h"

//Harmonic Bonds, Lennard-Jones Repulsions
class HLJ1D: public BasePotentialModel1D{

public:

    //Constructor
    HLJ1D();

    //Virtual to define
    double bndPotential(double& x0, double& x1, int& p) override;
    double repPotential(double& x0, double& x1, int& p) override;
    double gcnPotential(double& x0) override;
    void bndForce(double& x0, double& x1, double& f0, double &f1, int& p) override;
    void repForce(double& x0, double& x1, double& f0, double &f1, int& p) override;
    void gcnForce(double& x0, double& f0) override;

};

//Harmonic Bonds, Lennard-Jones Repulsions, Periodic Boundary Conditions
class HLJ1DP: public BasePotentialModel1D{

public:

    //Constructor
    HLJ1DP(double periodicX);

    //Periodic Boundary
    double pbx,pbrx; //x cell length and reciprocal
    void wrap(VecF<double>& x); //wrap coordinates inside periodic cell

    //Virtual to define
    double bndPotential(double& x0, double& x1, int& p) override;
    double repPotential(double& x0, double& x1, int& p) override;
    double gcnPotential(double& x0) override;
    void bndForce(double& x0, double& x1, double& f0, double &f1, int& p) override;
    void repForce(double& x0, double& x1, double& f0, double &f1, int& p) override;
    void gcnForce(double& x0, double& f0) override;

};

#endif //NL_POT1D_H
