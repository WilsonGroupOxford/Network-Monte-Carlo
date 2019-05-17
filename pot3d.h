//3D Potential Models Compatible with Optimistaion Routines
#ifndef NL_POT3D_H
#define NL_POT3D_H

#include <iostream>
#include "pot.h"

//Harmonic Bonds, Lennard-Jones Repulsions
class HLJ3D: public BasePotentialModel3D{

public:

    //Constructor
    HLJ3D();

    //Virtual to define
    double bndPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, int& p) override;
    double angPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2, int& p) override;
    double repPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, int& p) override;
    double intPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2, double& x3, double& y3, double& z3, int& p) override;
    double gcnPotential(double& x0, double& y0, double& z0) override;
    void bndForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1,
                          double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1,  int& p) override;
    void angForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2,
                          double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1, double& fx2, double& fy2, double& fz2,  int& p) override;
    void repForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1,
                          double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1,  int& p) override;
    void intForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2, double& x3, double& y3, double& z3,
                          double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1, double& fx2, double& fy2, double& fz2, double& fx3, double& fy3, double& fz3, int& p) override;
    void gcnForce(double& x0, double& y0, double& z0, double& fx0, double& fy0, double& fz0) override;

};

//Harmonic Bonds, Lennard-Jones Repulsions, Periodic Boundary
class HLJ3DP: public BasePotentialModel3D{

public:

    //Constructor
    HLJ3DP(double periodicX, double periodicY, double periodicZ);

    //Periodic Boundary
    double pbx,pby,pbz,pbrx,pbry,pbrz; //cell lengths and reciprocals
    void wrap(VecF<double>& x); //wrap coordinates inside periodic cell

    //Virtual to define
    double bndPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, int& p) override;
    double angPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2, int& p) override;
    double repPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, int& p) override;
    double intPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2, double& x3, double& y3, double& z3, int& p) override;
    double gcnPotential(double& x0, double& y0, double& z0) override;
    void bndForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1,
                  double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1,  int& p) override;
    void angForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2,
                  double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1, double& fx2, double& fy2, double& fz2,  int& p) override;
    void repForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1,
                  double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1,  int& p) override;
    void intForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2, double& x3, double& y3, double& z3,
                  double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1, double& fx2, double& fy2, double& fz2, double& fx3, double& fy3, double& fz3, int& p) override;
    void gcnForce(double& x0, double& y0, double& z0, double& fx0, double& fy0, double& fz0) override;

};

//Harmonic Bonds, Lennard-Jones Repulsions, Constrained to Sphere
class HLJ3DS: public BasePotentialModel3D{

public:

    //Constructor
    HLJ3DS();

    //Virtual to define
    double bndPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, int& p) override;
    double angPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2, int& p) override;
    double repPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, int& p) override;
    double intPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2, double& x3, double& y3, double& z3, int& p) override;
    double gcnPotential(double& x0, double& y0, double& z0) override;
    void bndForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1,
                  double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1,  int& p) override;
    void angForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2,
                  double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1, double& fx2, double& fy2, double& fz2,  int& p) override;
    void repForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1,
                  double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1,  int& p) override;
    void intForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2, double& x3, double& y3, double& z3,
                  double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1, double& fx2, double& fy2, double& fz2, double& fx3, double& fy3, double& fz3, int& p) override;
    void gcnForce(double& x0, double& y0, double& z0, double& fx0, double& fy0, double& fz0) override;

};

//Harmonic Bonds and Angles, Infinite Cost Arc Intersections, Constrained to Sphere
class HI3DS: public BasePotentialModel3D{

public:

    //Constructor
    HI3DS();

    //Virtual to define
    double bndPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, int& p) override;
    double angPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2, int& p) override;
    double repPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, int& p) override;
    double intPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2, double& x3, double& y3, double& z3, int& p) override;
    double gcnPotential(double& x0, double& y0, double& z0) override;
    void bndForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1,
                  double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1,  int& p) override;
    void angForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2,
                  double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1, double& fx2, double& fy2, double& fz2,  int& p) override;
    void repForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1,
                  double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1,  int& p) override;
    void intForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2, double& x3, double& y3, double& z3,
                  double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1, double& fx2, double& fy2, double& fz2, double& fx3, double& fy3, double& fz3, int& p) override;
    void gcnForce(double& x0, double& y0, double& z0, double& fx0, double& fy0, double& fz0) override;
};

//Harmonic Bonds and Restricted Angles, Infinite Cost Arc Intersections, Constrained to Sphere
class HRI3DS: public BasePotentialModel3D{

public:

    //Constructor
    HRI3DS();

    //Virtual to define
    double bndPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, int& p) override;
    double angPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2, int& p) override;
    double repPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, int& p) override;
    double intPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2, double& x3, double& y3, double& z3, int& p) override;
    double gcnPotential(double& x0, double& y0, double& z0) override;
    void bndForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1,
                  double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1,  int& p) override;
    void angForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2,
                  double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1, double& fx2, double& fy2, double& fz2,  int& p) override;
    void repForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1,
                  double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1,  int& p) override;
    void intForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2, double& x3, double& y3, double& z3,
                  double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1, double& fx2, double& fy2, double& fz2, double& fx3, double& fy3, double& fz3, int& p) override;
    void gcnForce(double& x0, double& y0, double& z0, double& fx0, double& fy0, double& fz0) override;
};

//Addition 3D vector functions
inline VecF<double> crossProduct(VecF<double>& v0, VecF<double>& v1);
inline double dotProduct(VecF<double>& v0, VecF<double>& v1);
inline double orthodromicDist(VecF<double>& v0, VecF<double>& v1);

#endif //NL_POT3D_H