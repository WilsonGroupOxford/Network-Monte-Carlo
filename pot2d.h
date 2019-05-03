//2D Potential Models Compatible with Optimistaion Routines
#ifndef NL_POT2D_H
#define NL_POT2D_H

#include <iostream>
#include "pot.h"

//Harmonic Bonds and Angles, Lennard-Jones Repulsions
class HLJ2D: public BasePotentialModel2D{

public:

    //Constructor
    HLJ2D();

    //Virtual to define
    double bndPotential(double& x0, double& y0, double& x1, double& y1, int& p) override;
    double angPotential(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, int& p) override;
    double repPotential(double& x0, double& y0, double& x1, double& y1, int& p) override;
    double intPotential(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, double& x3, double& y3, int& p) override;
    double gcnPotential(double& x0, double& y0) override;
    void bndForce(double& x0, double& y0, double& x1, double& y1,
                          double& fx0, double& fy0, double& fx1, double& fy1, int& p) override;
    void angForce(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2,
                          double& fx0, double& fy0, double& fx1, double& fy1, double& fx2, double& fy2, int& p) override;
    void repForce(double& x0, double& y0, double& x1, double& y1,
                          double& fx0, double& fy0, double& fx1, double& fy1,int& p) override;
    void intForce(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, double& x3, double& y3,
                          double& fx0, double& fy0, double& fx1, double& fy1, double& fx2, double& fy2, double& fx3, double& fy3, int& p) override;
    void gcnForce(double& x0, double& y0, double& fx0, double& fy0) override;

};

//Harmonic Bonds, Lennard-Jones Repulsions, Periodic Boundary Conditions
class HLJ2DP: public BasePotentialModel2D{

public:

    //Constructor
    HLJ2DP(double periodicX, double periodicY);

    //Periodic Boundary
    double pbx,pby,pbrx,pbry; //cell lengths and reciprocals
    void wrap(VecF<double>& x); //wrap coordinates inside periodic cell


    //Virtual to define
    double bndPotential(double& x0, double& y0, double& x1, double& y1, int& p) override;
    double angPotential(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, int& p) override;
    double repPotential(double& x0, double& y0, double& x1, double& y1, int& p) override;
    double intPotential(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, double& x3, double& y3, int& p) override;
    double gcnPotential(double& x0, double& y0) override;
    void bndForce(double& x0, double& y0, double& x1, double& y1,
                  double& fx0, double& fy0, double& fx1, double& fy1, int& p) override;
    void angForce(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2,
                  double& fx0, double& fy0, double& fx1, double& fy1, double& fx2, double& fy2, int& p) override;
    void repForce(double& x0, double& y0, double& x1, double& y1,
                  double& fx0, double& fy0, double& fx1, double& fy1,int& p) override;
    void intForce(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, double& x3, double& y3,
                  double& fx0, double& fy0, double& fx1, double& fy1, double& fx2, double& fy2, double& fx3, double& fy3, int& p) override;
    void gcnForce(double& x0, double& y0, double& fx0, double& fy0) override;

};

//Harmonic Bonds, Lennard-Jones Repulsions, Constrained to Circle
class HLJ2DC: public BasePotentialModel2D{

public:

    //Constructor
    HLJ2DC();

    //Virtual to define
    double bndPotential(double& x0, double& y0, double& x1, double& y1, int& p) override;
    double angPotential(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, int& p) override;
    double repPotential(double& x0, double& y0, double& x1, double& y1, int& p) override;
    double intPotential(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, double& x3, double& y3, int& p) override;
    double gcnPotential(double& x0, double& y0) override;
    void bndForce(double& x0, double& y0, double& x1, double& y1,
                  double& fx0, double& fy0, double& fx1, double& fy1, int& p) override;
    void angForce(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2,
                  double& fx0, double& fy0, double& fx1, double& fy1, double& fx2, double& fy2, int& p) override;
    void repForce(double& x0, double& y0, double& x1, double& y1,
                  double& fx0, double& fy0, double& fx1, double& fy1,int& p) override;
    void intForce(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, double& x3, double& y3,
                  double& fx0, double& fy0, double& fx1, double& fy1, double& fx2, double& fy2, double& fx3, double& fy3, int& p) override;
    void gcnForce(double& x0, double& y0, double& fx0, double& fy0) override;

};

//Harmonic Bonds, Infinite Cost Line Intersections
class HI2D: public BasePotentialModel2D{

public:

    //Constructor
    HI2D();

    //Virtual to define
    double bndPotential(double& x0, double& y0, double& x1, double& y1, int& p) override;
    double angPotential(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, int& p) override;
    double repPotential(double& x0, double& y0, double& x1, double& y1, int& p) override;
    double intPotential(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, double& x3, double& y3, int& p) override;
    double gcnPotential(double& x0, double& y0) override;
    void bndForce(double& x0, double& y0, double& x1, double& y1,
                  double& fx0, double& fy0, double& fx1, double& fy1, int& p) override;
    void angForce(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2,
                  double& fx0, double& fy0, double& fx1, double& fy1, double& fx2, double& fy2, int& p) override;
    void repForce(double& x0, double& y0, double& x1, double& y1,
                  double& fx0, double& fy0, double& fx1, double& fy1,int& p) override;
    void intForce(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, double& x3, double& y3,
                  double& fx0, double& fy0, double& fx1, double& fy1, double& fx2, double& fy2, double& fx3, double& fy3, int& p) override;
    void gcnForce(double& x0, double& y0, double& fx0, double& fy0) override;

};

//Harmonic Bonds, Infinite Cost Line Intersections, Periodic Boundary Conditions
class HI2DP: public BasePotentialModel2D{

public:

    //Constructor
    HI2DP(double periodicX, double periodicY);

    //Periodic Boundary
    double pbx,pby,pbrx,pbry; //cell lengths and reciprocals
    void wrap(VecF<double>& x); //wrap coordinates inside periodic cell

    //Virtual to define
    double bndPotential(double& x0, double& y0, double& x1, double& y1, int& p) override;
    double angPotential(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, int& p) override;
    double repPotential(double& x0, double& y0, double& x1, double& y1, int& p) override;
    double intPotential(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, double& x3, double& y3, int& p) override;
    double gcnPotential(double& x0, double& y0) override;
    void bndForce(double& x0, double& y0, double& x1, double& y1,
                  double& fx0, double& fy0, double& fx1, double& fy1, int& p) override;
    void angForce(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2,
                  double& fx0, double& fy0, double& fx1, double& fy1, double& fx2, double& fy2, int& p) override;
    void repForce(double& x0, double& y0, double& x1, double& y1,
                  double& fx0, double& fy0, double& fx1, double& fy1,int& p) override;
    void intForce(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, double& x3, double& y3,
                  double& fx0, double& fy0, double& fx1, double& fy1, double& fx2, double& fy2, double& fx3, double& fy3, int& p) override;
    void gcnForce(double& x0, double& y0, double& fx0, double& fy0) override;

};

//Harmonic Windowed Bonds and Angles, Harmonic Angles, Infinite Cost Line Intersections
class HWI2D: public BasePotentialModel2D{

public:

    //Constructor
    HWI2D();

    //Virtual to define
    double bndPotential(double& x0, double& y0, double& x1, double& y1, int& p) override;
    double angPotential(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, int& p) override;
    double repPotential(double& x0, double& y0, double& x1, double& y1, int& p) override;
    double intPotential(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, double& x3, double& y3, int& p) override;
    double gcnPotential(double& x0, double& y0) override;
    void bndForce(double& x0, double& y0, double& x1, double& y1,
                  double& fx0, double& fy0, double& fx1, double& fy1, int& p) override;
    void angForce(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2,
                  double& fx0, double& fy0, double& fx1, double& fy1, double& fx2, double& fy2, int& p) override;
    void repForce(double& x0, double& y0, double& x1, double& y1,
                  double& fx0, double& fy0, double& fx1, double& fy1,int& p) override;
    void intForce(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, double& x3, double& y3,
                  double& fx0, double& fy0, double& fx1, double& fy1, double& fx2, double& fy2, double& fx3, double& fy3, int& p) override;
    void gcnForce(double& x0, double& y0, double& fx0, double& fy0) override;

};

//Harmonic Windowed Bonds and Angles, Harmonic Angles, Infinite Cost Line Intersections, Periodic Boundary Conditions
class HWI2DP: public BasePotentialModel2D{

public:

    //Constructor
    HWI2DP(double periodicX, double periodicY);

    //Periodic Boundary
    double pbx,pby,pbrx,pbry; //cell lengths and reciprocals
    void wrap(VecF<double>& x); //wrap coordinates inside periodic cell

    //Virtual to define
    double bndPotential(double& x0, double& y0, double& x1, double& y1, int& p) override;
    double angPotential(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, int& p) override;
    double repPotential(double& x0, double& y0, double& x1, double& y1, int& p) override;
    double intPotential(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, double& x3, double& y3, int& p) override;
    double gcnPotential(double& x0, double& y0) override;
    void bndForce(double& x0, double& y0, double& x1, double& y1,
                  double& fx0, double& fy0, double& fx1, double& fy1, int& p) override;
    void angForce(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2,
                  double& fx0, double& fy0, double& fx1, double& fy1, double& fx2, double& fy2, int& p) override;
    void repForce(double& x0, double& y0, double& x1, double& y1,
                  double& fx0, double& fy0, double& fx1, double& fy1,int& p) override;
    void intForce(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, double& x3, double& y3,
                  double& fx0, double& fy0, double& fx1, double& fy1, double& fx2, double& fy2, double& fx3, double& fy3, int& p) override;
    void gcnForce(double& x0, double& y0, double& fx0, double& fy0) override;

};

//Functions for line intersections, implemeted from "Computational Geometry in C"
inline double signedAreaSqTriangle(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2){
    return x0*y1-y0*x1+y0*x2-x0*y2+x1*y2-x2*y1;
}
inline bool leftTriangle(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2){
    return signedAreaSqTriangle(x0,y0,x1,y1,x2,y2)>0.0;
}
inline bool collinearPoints(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, double threshold=1e-6){
    return fabs(signedAreaSqTriangle(x0,y0,x1,y1,x2,y2))<threshold;
}
inline bool properIntersectionLines(double& x1a, double& y1a, double& x1b, double& y1b, double& x2a, double& y2a, double& x2b, double& y2b){
    if(collinearPoints(x1a,y1a,x1b,y1b,x2a,y2a)
       || collinearPoints(x1a,y1a,x1b,y1b,x2b,y2b)
       || collinearPoints(x2a,y2a,x2b,y2b,x1a,y1a)
       || collinearPoints(x2a,y2a,x2b,y2b,x1b,y1b)) return false;
    return (leftTriangle(x1a,y1a,x1b,y1b,x2a,y2a)^leftTriangle(x1a,y1a,x1b,y1b,x2b,y2b))
           && (leftTriangle(x2a,y2a,x2b,y2b,x1a,y1a)^leftTriangle(x2a,y2a,x2b,y2b,x1b,y1b));
}

#endif //NL_POT2D_H