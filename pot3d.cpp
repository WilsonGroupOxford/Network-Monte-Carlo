#include "pot3d.h"

//##### HARMONIC BONDS, LJ REPULSIONS #####
HLJ3D::HLJ3D():BasePotentialModel3D(){};

//Potential
double HLJ3D::bndPotential(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, int &p) {
    double dx=x1-x0;
    double dy=y1-y0;
    double dz=z1-z0;
    double r=sqrt((dx*dx+dy*dy+dz*dz));
    return 0.5*bndP[2*p]*pow((r-bndP[2*p+1]),2);
}

double HLJ3D::angPotential(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &x2, double &y2, double &z2, int &p) {
    return 0.0;
}

double HLJ3D::repPotential(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, int &p) {
    double u=0.0;
    double dx=x1-x0;
    double dy=y1-y0;
    double dz=z1-z0;
    double r02=repP[2*p];
    double r2=(dx*dx+dy*dy+dz*dz);
    if(r2<r02){
        double d2=r02/r2;
        double d12=pow(d2,6);
        double d24=pow(d12,2);
        double ep=repP[2*p+1];
        u+=ep*(d24-2.0*d12)+ep;
    }
    return u;
}

double HLJ3D::intPotential(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &x2,
                           double &y2, double &z2, double &x3, double &y3, double &z3, int &p) {
    return 0.0;
}

double HLJ3D::gcnPotential(double &x0, double &y0, double &z0) {
    return 0.0;
}

//Forces
void HLJ3D::bndForce(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1,
                     double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, int &p) {
    double dx=x1-x0;
    double dy=y1-y0;
    double dz=z1-z0;
    double r=sqrt(dx*dx+dy*dy+dz*dz);
    double m=-bndP[2*p]*(r-bndP[2*p+1])/r;
    dx*=m;
    dy*=m;
    dz*=m;
    fx0-=dx;
    fy0-=dy;
    fz0-=dz;
    fx1+=dx;
    fy1+=dy;
    fz1+=dz;
}

void HLJ3D::angForce(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &x2, double &y2, double &z2,
                     double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, double &fx2, double &fy2, double &fz2, int &p) {
    return;
}

void HLJ3D::repForce(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1,
                     double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, int &p) {
    double dx=x1-x0;
    double dy=y1-y0;
    double dz=z1-z0;
    double r02=repP[2*p];
    double r2=(dx*dx+dy*dy+dz*dz);
    if(r2<r02){
        double d2=r02/r2;
        double d12=pow(d2,6);
        double d24=pow(d12,2);
        double ep=repP[2*p+1];
        double m=24.0*ep*(d24-d12)/r2;
        dx*=m;
        dy*=m;
        dz*=m;
        fx0-=dx;
        fy0-=dy;
        fz0-=dz;
        fx1+=dx;
        fy1+=dy;
        fz1+=dz;
    }
}

void HLJ3D::intForce(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &x2, double &y2,
                     double &z2, double &x3, double &y3, double &z3, double &fx0, double &fy0, double &fz0, double &fx1,
                     double &fy1, double &fz1, double &fx2, double &fy2, double &fz2, double &fx3, double &fy3,
                     double &fz3, int &p) {
    return;
}

void HLJ3D::gcnForce(double &x0, double &y0, double &z0, double &fx0, double &fy0, double &fz0) {
    return;
}

//##### HARMONIC BONDS, LJ REPULSIONS, PERIODIC BOUNDARY #####
HLJ3DP::HLJ3DP(double periodicX, double periodicY, double periodicZ):BasePotentialModel3D(){
    pbx=periodicX;
    pby=periodicY;
    pbz=periodicZ;
    pbrx=1.0/pbx;
    pbry=1.0/pby;
    pbrz=1.0/pbz;
};

//Potential
double HLJ3DP::bndPotential(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, int &p) {
    double dx=x1-x0;
    double dy=y1-y0;
    double dz=z1-z0;
    dx-=pbx*nearbyint(dx*pbrx);
    dy-=pby*nearbyint(dy*pbry);
    dz-=pbz*nearbyint(dz*pbrz);
    double r=sqrt((dx*dx+dy*dy+dz*dz));
    return 0.5*bndP[2*p]*pow((r-bndP[2*p+1]),2);
}

double HLJ3DP::angPotential(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &x2, double &y2, double &z2, int &p) {
    return 0.0;
}

double HLJ3DP::repPotential(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, int &p) {
    double u=0.0;
    double dx=x1-x0;
    double dy=y1-y0;
    double dz=z1-z0;
    dx-=pbx*nearbyint(dx*pbrx);
    dy-=pby*nearbyint(dy*pbry);
    dz-=pbz*nearbyint(dz*pbrz);
    double r02=repP[2*p];
    double r2=(dx*dx+dy*dy+dz*dz);
    if(r2<r02){
        double d2=r02/r2;
        double d12=pow(d2,6);
        double d24=pow(d12,2);
        double ep=repP[2*p+1];
        u+=ep*(d24-2.0*d12)+ep;
    }
    return u;
}

double HLJ3DP::intPotential(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &x2,
                           double &y2, double &z2, double &x3, double &y3, double &z3, int &p) {
    return 0.0;
}

double HLJ3DP::gcnPotential(double &x0, double &y0, double &z0) {
    return 0.0;
}

//Forces
void HLJ3DP::bndForce(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1,
                     double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, int &p) {
    double dx=x1-x0;
    double dy=y1-y0;
    double dz=z1-z0;
    dx-=pbx*nearbyint(dx*pbrx);
    dy-=pby*nearbyint(dy*pbry);
    dz-=pbz*nearbyint(dz*pbrz);
    double r=sqrt(dx*dx+dy*dy+dz*dz);
    double m=-bndP[2*p]*(r-bndP[2*p+1])/r;
    dx*=m;
    dy*=m;
    dz*=m;
    fx0-=dx;
    fy0-=dy;
    fz0-=dz;
    fx1+=dx;
    fy1+=dy;
    fz1+=dz;
}

void HLJ3DP::angForce(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &x2, double &y2, double &z2,
                     double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, double &fx2, double &fy2, double &fz2, int &p) {
    return;
}

void HLJ3DP::repForce(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1,
                     double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, int &p) {
    double dx=x1-x0;
    double dy=y1-y0;
    double dz=z1-z0;
    dx-=pbx*nearbyint(dx*pbrx);
    dy-=pby*nearbyint(dy*pbry);
    dz-=pbz*nearbyint(dz*pbrz);
    double r02=repP[2*p];
    double r2=(dx*dx+dy*dy+dz*dz);
    if(r2<r02){
        double d2=r02/r2;
        double d12=pow(d2,6);
        double d24=pow(d12,2);
        double ep=repP[2*p+1];
        double m=24.0*ep*(d24-d12)/r2;
        dx*=m;
        dy*=m;
        dz*=m;
        fx0-=dx;
        fy0-=dy;
        fz0-=dz;
        fx1+=dx;
        fy1+=dy;
        fz1+=dz;
    }
}

void HLJ3DP::intForce(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &x2, double &y2,
                     double &z2, double &x3, double &y3, double &z3, double &fx0, double &fy0, double &fz0, double &fx1,
                     double &fy1, double &fz1, double &fx2, double &fy2, double &fz2, double &fx3, double &fy3,
                     double &fz3, int &p) {
    return;
}

void HLJ3DP::gcnForce(double &x0, double &y0, double &z0, double &fx0, double &fy0, double &fz0) {
    return;
}

void HLJ3DP::wrap(VecF<double>& x) {
    for(int i=0,j=1,k=2; i<x.n; i+=3, j+=3, k+=3){
        x[i]-=pbx*nearbyint(x[i]*pbrx)-pbx*0.5;
        x[j]-=pby*nearbyint(x[j]*pbry)-pby*0.5;
        x[k]-=pbz*nearbyint(x[k]*pbrz)-pbz*0.5;
    }
}


//##### HARMONIC BONDS, LJ REPULSIONS, CONSTRAINED TO SPHERE #####
HLJ3DS::HLJ3DS():BasePotentialModel3D(){};

//Potential
double HLJ3DS::bndPotential(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, int &p) {
    double dx=x1-x0;
    double dy=y1-y0;
    double dz=z1-z0;
    double r=sqrt((dx*dx+dy*dy+dz*dz));
    return 0.5*bndP[2*p]*pow((r-bndP[2*p+1]),2);
}

double HLJ3DS::angPotential(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &x2, double &y2, double &z2, int &p) {
    return 0.0;
}

double HLJ3DS::repPotential(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, int &p) {
    double u=0.0;
    double dx=x1-x0;
    double dy=y1-y0;
    double dz=z1-z0;
    double r02=repP[2*p];
    double r2=(dx*dx+dy*dy+dz*dz);
    if(r2<r02){
        double d2=r02/r2;
        double d12=pow(d2,6);
        double d24=pow(d12,2);
        double ep=repP[2*p+1];
        u+=ep*(d24-2.0*d12)+ep;
    }
    return u;
}

double HLJ3DS::intPotential(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &x2,
                           double &y2, double &z2, double &x3, double &y3, double &z3, int &p) {
    return 0.0;
}

double HLJ3DS::gcnPotential(double &x0, double &y0, double &z0) {

    double r=sqrt((x0*x0+y0*y0+z0*z0));
    return 0.5*gcnP[0]*pow((r-gcnP[1]),2);
}

//Forces
void HLJ3DS::bndForce(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1,
                     double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, int &p) {
    double dx=x1-x0;
    double dy=y1-y0;
    double dz=z1-z0;
    double r=sqrt(dx*dx+dy*dy+dz*dz);
    double m=-bndP[2*p]*(r-bndP[2*p+1])/r;
    dx*=m;
    dy*=m;
    dz*=m;
    fx0-=dx;
    fy0-=dy;
    fz0-=dz;
    fx1+=dx;
    fy1+=dy;
    fz1+=dz;
}

void HLJ3DS::angForce(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &x2, double &y2, double &z2,
                     double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, double &fx2, double &fy2, double &fz2, int &p) {
    return;
}

void HLJ3DS::repForce(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1,
                     double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, int &p) {
    double dx=x1-x0;
    double dy=y1-y0;
    double dz=z1-z0;
    double r02=repP[2*p];
    double r2=(dx*dx+dy*dy+dz*dz);
    if(r2<r02){
        double d2=r02/r2;
        double d12=pow(d2,6);
        double d24=pow(d12,2);
        double ep=repP[2*p+1];
        double m=24.0*ep*(d24-d12)/r2;
        dx*=m;
        dy*=m;
        dz*=m;
        fx0-=dx;
        fy0-=dy;
        fz0-=dz;
        fx1+=dx;
        fy1+=dy;
        fz1+=dz;
    }
}

void HLJ3DS::intForce(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &x2, double &y2,
                     double &z2, double &x3, double &y3, double &z3, double &fx0, double &fy0, double &fz0, double &fx1,
                     double &fy1, double &fz1, double &fx2, double &fy2, double &fz2, double &fx3, double &fy3,
                     double &fz3, int &p) {
    return;
}

void HLJ3DS::gcnForce(double &x0, double &y0, double &z0, double &fx0, double &fy0, double &fz0) {
    double r=sqrt(x0*x0+y0*y0+z0*z0);
    double m=-gcnP[0]*(r-gcnP[1])/r;
    double dx=m*x0;
    double dy=m*y0;
    double dz=m*z0;
    fx0+=dx;
    fy0+=dy;
    fz0+=dz;
    return;
}






