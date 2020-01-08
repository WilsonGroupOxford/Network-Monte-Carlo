//Base Potential Model Compatible with Optimistaion Routines
#ifndef NL_POT_H
#define NL_POT_H

#include <iostream>
#include "opt.h"

//Base potential model
class BasePotentialModel: public FuncGradMultiDim{

public:

    //Model parameters
    bool useBnds,useAngs,useReps,useIntx,useFixd,useGcns; //flags to include interaction types
    VecF<int> bnds; //ids of atoms in bond(pair)
    VecF<int> angs; //ids of atoms in angle(trio)
    VecF<int> reps; //ids of atoms in repulsion(pair)
    VecF<int> intx; //ids of atoms in intersections(quad)
    VecF<int> fixd; //ids of fixed atoms
    VecF<int> gcns; //ids of atoms constrained by geometry
    VecF<double> bndP; //bond parameters
    VecF<double> angP; //angle parameters
    VecF<double> repP; //rep parameters
    VecF<double> intP; //intx parameters
    VecF<double> gcnP; //gcns parameters


    //Constructor, Setters
    BasePotentialModel();
    void setBonds(VecF<int> bonds, VecF<double> params);
    void setAngles(VecF<int> angles, VecF<double> params);
    void setRepulsions(VecF<int> repulsions, VecF<double> params);
    void setIntersections(VecF<int> intersections, VecF<double> params);
    void setFixedAtoms(VecF<int> fixed);
    void setGeomConstraints(VecF<int> constrained, VecF<double> params);
    void reset();

    //Virtual to define
    double function(VecF<double>& x) override;
    VecF<double> gradient(VecF<double>& x) override;

    //Potentials
    virtual double bndsPotential(VecF<double>& x)=0;
    virtual double angsPotential(VecF<double>& x)=0;
    virtual double repsPotential(VecF<double>& x)=0;
    virtual double intxPotential(VecF<double>& x)=0;
    virtual double gcnsPotential(VecF<double>& x)=0;

    //Forces
    virtual void bndsForce(VecF<double>& f,VecF<double>& x)=0;
    virtual void angsForce(VecF<double>& f,VecF<double>& x)=0;
    virtual void repsForce(VecF<double>& f,VecF<double>& x)=0;
    virtual void intxForce(VecF<double>& f,VecF<double>& x)=0;
    virtual void fixdForce(VecF<double>& f,VecF<double>& x)=0;
    virtual void gcnsForce(VecF<double>& f,VecF<double>& x)=0;
};


//Base model for 1 dimensional systems
class BasePotentialModel1D: public BasePotentialModel{

public:

    //Constructor
    BasePotentialModel1D();

    //Virtual to define
    double bndsPotential(VecF<double>& x) override;
    double angsPotential(VecF<double>& x) override;
    double repsPotential(VecF<double>& x) override;
    double intxPotential(VecF<double>& x) override;
    double gcnsPotential(VecF<double>& x) override;

    void bndsForce(VecF<double>& f,VecF<double>& x) override;
    void angsForce(VecF<double>& f,VecF<double>& x) override;
    void repsForce(VecF<double>& f,VecF<double>& x) override;
    void intxForce(VecF<double>& f,VecF<double>& x) override;
    void fixdForce(VecF<double>& f,VecF<double>& x) override;
    void gcnsForce(VecF<double>& f,VecF<double>& x) override;

    //Individual potentials and forces
    virtual double bndPotential(double& x0, double& x1, int& p)=0;
    virtual double repPotential(double& x0, double& x1, int& p)=0;
    virtual double gcnPotential(double& x0)=0;
    virtual void bndForce(double& x0, double& x1, double& f0, double& f1, int& p)=0;
    virtual void repForce(double& x0, double& x1, double& f0, double& f1, int& p)=0;
    virtual void gcnForce(double& x0, double& f0)=0;
};

//Base model for 2 dimensional systems
class BasePotentialModel2D: public BasePotentialModel{

public:

    //Constructor
    BasePotentialModel2D();

    //Virtual to define
    double bndsPotential(VecF<double>& x) override;
    double angsPotential(VecF<double>& x) override;
    double repsPotential(VecF<double>& x) override;
    double intxPotential(VecF<double>& x) override;
    double gcnsPotential(VecF<double>& x) override;

    void bndsForce(VecF<double>& f,VecF<double>& x) override;
    void angsForce(VecF<double>& f,VecF<double>& x) override;
    void repsForce(VecF<double>& f,VecF<double>& x) override;
    void intxForce(VecF<double>& f,VecF<double>& x) override;
    void fixdForce(VecF<double>& f,VecF<double>& x) override;
    void gcnsForce(VecF<double>& f,VecF<double>& x) override;

    //Individual potentials and forces
    virtual double bndPotential(double& x0, double& y0, double& x1, double& y1, int& p)=0;
    virtual double angPotential(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, int& p)=0;
    virtual double repPotential(double& x0, double& y0, double& x1, double& y1, int& p)=0;
    virtual double intPotential(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, double& x3, double& y3, int& p)=0;
    virtual double gcnPotential(double& x0, double& y0)=0;
    virtual void bndForce(double& x0, double& y0, double& x1, double& y1,
                          double& fx0, double& fy0, double& fx1, double& fy1, int& p)=0;
    virtual void angForce(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2,
                          double& fx0, double& fy0, double& fx1, double& fy1, double& fx2, double& fy2, int& p)=0;
    virtual void repForce(double& x0, double& y0, double& x1, double& y1,
                          double& fx0, double& fy0, double& fx1, double& fy1,int& p)=0;
    virtual void intForce(double& x0, double& y0, double& x1, double& y1, double& x2, double& y2, double& x3, double& y3,
                          double& fx0, double& fy0, double& fx1, double& fy1, double& fx2, double& fy2, double& fx3, double& fy3, int& p)=0;
    virtual void gcnForce(double& x0, double& y0, double& fx0, double& fy0)=0;

};

//Base model for 3 dimensional systems
class BasePotentialModel3D: public BasePotentialModel{

public:

    //Constructor
    BasePotentialModel3D();

    //Virtual to define
    double bndsPotential(VecF<double>& x) override;
    double angsPotential(VecF<double>& x) override;
    double repsPotential(VecF<double>& x) override;
    double intxPotential(VecF<double>& x) override;
    double gcnsPotential(VecF<double>& x) override;

    void bndsForce(VecF<double>& f,VecF<double>& x) override;
    void angsForce(VecF<double>& f,VecF<double>& x) override;
    void repsForce(VecF<double>& f,VecF<double>& x) override;
    void intxForce(VecF<double>& f,VecF<double>& x) override;
    void fixdForce(VecF<double>& f,VecF<double>& x) override;
    void gcnsForce(VecF<double>& f,VecF<double>& x) override;

    //Individual potentials and forces
    virtual double bndPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, int& p)=0;
    virtual double angPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2, int& p)=0;
    virtual double repPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, int& p)=0;
    virtual double intPotential(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2, double& x3, double& y3, double& z3, int& p)=0;
    virtual double gcnPotential(double& x0, double& y0, double& z0)=0;
    virtual void bndForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1,
                              double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1,  int& p)=0;
    virtual void angForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2,
                          double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1, double& fx2, double& fy2, double& fz2,  int& p)=0;
    virtual void repForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1,
                          double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1,  int& p)=0;
    virtual void intForce(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2, double& x3, double& y3, double& z3,
                          double& fx0, double& fy0, double& fz0, double& fx1, double& fy1, double& fz1, double& fx2, double& fy2, double& fz2, double& fx3, double& fy3, double& fz3, int& p)=0;
    virtual void gcnForce(double& x0, double& y0, double& z0, double& fx0, double& fy0, double& fz0)=0;

};
#endif //NL_POT_H
