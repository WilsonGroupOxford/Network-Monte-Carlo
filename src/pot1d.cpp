#include "pot1d.h"

//##### HARMONIC BONDS, LJ REPULSIONS #####
HLJ1D::HLJ1D():BasePotentialModel1D(){};

//Potential
double HLJ1D::bndPotential(double &x0, double &x1, int &p) {
    double r=abs(x1-x0);
    return 0.5*bndP[2*p]*pow((r-bndP[2*p+1]),2);
}

double HLJ1D::repPotential(double &x0, double &x1, int &p) {
    double u=0.0;
    double dx=(x1-x0);
    double r02=repP[2*p];
    double r2=dx*dx;
    if(r2<r02){
        double d2=r02/r2;
        double d12=pow(d2,6);
        double d24=pow(d12,2);
        double ep=repP[2*p+1];
        u+=ep*(d24-2.0*d12)+ep;
    }
    return u;
}

double HLJ1D::gcnPotential(double &x0) {
    return 0.0;
}

//Force
void HLJ1D::bndForce(double &x0, double &x1, double &f0, double &f1, int &p) {
    double dx=x1-x0;
    double r=abs(dx);
    double mag=-bndP[2*p]*(r-bndP[2*p+1])/r;
    dx*=mag;
    f0-=dx;
    f1+=dx;
}

void HLJ1D::repForce(double &x0, double &x1, double &f0, double &f1, int &p) {
    double dx=(x1-x0);
    double r02=repP[2*p];
    double r2=dx*dx;
    if(r2<r02){
        double d2=r02/r2;
        double d12=pow(d2,6);
        double d24=pow(d12,2);
        double ep=repP[2*p+1];
        double mag=24.0*ep*(d24-d12)/r2;
        dx*=mag;
        f0-=dx;
        f1+=dx;
    }
}

void HLJ1D::gcnForce(double &x0, double &f0) {
    return;
}

//##### HARMONIC BONDS, LJ REPULSIONS, PERIODIC BOUNDARY #####
HLJ1DP::HLJ1DP(double periodicX):BasePotentialModel1D(){
    pbx=periodicX;
    pbrx=1.0/pbx;
};

//Potential
double HLJ1DP::bndPotential(double &x0, double &x1, int &p) {
    double dx=x1-x0;
    dx-=pbx*nearbyint(dx*pbrx);
    double r=abs(dx);
    return 0.5*bndP[2*p]*pow((r-bndP[2*p+1]),2);
}

double HLJ1DP::repPotential(double &x0, double &x1, int &p) {
    double u=0.0;
    double dx=(x1-x0);
    dx-=pbx*nearbyint(dx*pbrx);
    double r02=repP[2*p];
    double r2=dx*dx;
    if(r2<r02){
        double d2=r02/r2;
        double d12=pow(d2,6);
        double d24=pow(d12,2);
        double ep=repP[2*p+1];
        u+=ep*(d24-2.0*d12)+ep;
    }
    return u;
}

double HLJ1DP::gcnPotential(double &x0) {
    return 0.0;
}

//Force
void HLJ1DP::bndForce(double &x0, double &x1, double &f0, double &f1, int &p) {
    double dx=x1-x0;
    dx-=pbx*nearbyint(dx*pbrx);
    double r=abs(dx);
    double mag=-bndP[2*p]*(r-bndP[2*p+1])/r;
    dx*=mag;
    f0-=dx;
    f1+=dx;
}

void HLJ1DP::repForce(double &x0, double &x1, double &f0, double &f1, int &p) {
    double dx=(x1-x0);
    dx-=pbx*nearbyint(dx*pbrx);
    double r02=repP[2*p];
    double r2=dx*dx;
    if(r2<r02){
        double d2=r02/r2;
        double d12=pow(d2,6);
        double d24=pow(d12,2);
        double ep=repP[2*p+1];
        double mag=24.0*ep*(d24-d12)/r2;
        dx*=mag;
        f0-=dx;
        f1+=dx;
    }
}

void HLJ1DP::gcnForce(double &x0, double &f0) {
    return;
}

void HLJ1DP::wrap(VecF<double>& x) {
    for(int i=0; i<x.n; ++i) x[i]-=pbx*nearbyint(x[i]*pbrx);
    x+=0.5*pbx;
}