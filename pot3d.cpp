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



//##### HARMONIC BONDS AND ANGLES, INFINITE COST ARC INTERSECTIONS, CONSTRAINED TO SPHERE #####
HI3DS::HI3DS():BasePotentialModel3D(){};

//Potential
double HI3DS::bndPotential(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, int &p) {
    double dx=x1-x0;
    double dy=y1-y0;
    double dz=z1-z0;
    double r=sqrt((dx*dx+dy*dy+dz*dz));
    return 0.5*bndP[2*p]*pow((r-bndP[2*p+1]),2);
}

double HI3DS::angPotential(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &x2, double &y2, double &z2, int &p) {
    double dx1=x1-x0;
    double dy1=y1-y0;
    double dz1=z1-z0;
    double dx2=x2-x0;
    double dy2=y2-y0;
    double dz2=z2-z0;
    double n1=1/sqrt(dx1*dx1+dy1*dy1+dz1*dz1);
    double n2=1/sqrt(dx2*dx2+dy2*dy2+dz2*dz2);
    dx1*=n1;
    dy1*=n1;
    dz1*=n1;
    dx2*=n2;
    dy2*=n2;
    dz2*=n2;
    double cosTheta=dx1*dx2+dy1*dy2+dz1*dz2;
    return 0.5*angP[2*p]*pow((cosTheta-angP[2*p+1]),2);
}

double HI3DS::repPotential(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, int &p) {
    return 0.0;
}

double HI3DS::intPotential(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &x2,
                            double &y2, double &z2, double &x3, double &y3, double &z3, int &p) {
    VecF<double> p0(3),p1(3),p2(3),p3(3);
    p0[0]=x0;
    p0[1]=y0;
    p0[2]=z0;
    p1[0]=x1;
    p1[1]=y1;
    p1[2]=z1;
    p2[0]=x2;
    p2[1]=y2;
    p2[2]=z2;
    p3[0]=x3;
    p3[1]=y3;
    p3[2]=z3;
    p0/=vNorm(p0);
    p1/=vNorm(p1);
    p2/=vNorm(p2);
    p3/=vNorm(p3);

    VecF<double> v0=crossProduct(p0,p1);
    VecF<double> v1=crossProduct(p2,p3);
    v0/=vNorm(v0);
    v1/=vNorm(v1);

    VecF<double> d=crossProduct(v0,v1);
    VecF<double> s0=d/vNorm(d);
    VecF<double> s1=s0*-1.0;

    bool arc0s0=fabs(orthodromicDist(p0,p1)-orthodromicDist(p0,s0)-orthodromicDist(p1,s0))<1e-8;
    bool arc1s0=fabs(orthodromicDist(p2,p3)-orthodromicDist(p2,s0)-orthodromicDist(p3,s0))<1e-8;
    if(arc0s0 && arc1s0) return numeric_limits<double>::infinity();

    bool arc0s1=fabs(orthodromicDist(p0,p1)-orthodromicDist(p0,s1)-orthodromicDist(p1,s1))<1e-8;
    bool arc1s1=fabs(orthodromicDist(p2,p3)-orthodromicDist(p2,s1)-orthodromicDist(p3,s1))<1e-8;
    if(arc0s1 && arc1s1) return numeric_limits<double>::infinity();

    return 0.0;
}

double HI3DS::gcnPotential(double &x0, double &y0, double &z0) {

    double r=sqrt((x0*x0+y0*y0+z0*z0));
    return 0.5*gcnP[0]*pow((r-gcnP[1]),2);
}

//Forces
void HI3DS::bndForce(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1,
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

void HI3DS::angForce(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &x2, double &y2, double &z2,
                      double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, double &fx2, double &fy2, double &fz2, int &p) {
    double dx1=x1-x0;
    double dy1=y1-y0;
    double dz1=z1-z0;
    double dx2=x2-x0;
    double dy2=y2-y0;
    double dz2=z2-z0;
    double n1=1/sqrt(dx1*dx1+dy1*dy1+dz1*dz1);
    double n2=1/sqrt(dx2*dx2+dy2*dy2+dz2*dz2);
    dx1*=n1;
    dy1*=n1;
    dz1*=n1;
    dx2*=n2;
    dy2*=n2;
    dz2*=n2;
    double cosTheta=dx1*dx2+dy1*dy2+dz1*dz2;
    double ff=angP[2*p]*(cosTheta-angP[2*p+1]);
    double ffx1=-ff*n1*(dx2-cosTheta*dx1);
    double ffy1=-ff*n1*(dy2-cosTheta*dy1);
    double ffz1=-ff*n1*(dz2-cosTheta*dz1);
    double ffx2=-ff*n2*(dx1-cosTheta*dx2);
    double ffy2=-ff*n2*(dy1-cosTheta*dy2);
    double ffz2=-ff*n2*(dz1-cosTheta*dz2);
    double ffx0=-ffx1-ffx2;
    double ffy0=-ffy1-ffy2;
    double ffz0=-ffz1-ffz2;
    fx0+=ffx0;
    fy0+=ffy0;
    fz0+=ffz0;
    fx1+=ffx1;
    fy1+=ffy1;
    fz1+=ffz1;
    fx2+=ffx2;
    fy2+=ffy2;
    fz2+=ffz2;
    return;
}

void HI3DS::repForce(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1,
                      double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, int &p) {
    return;
}

void HI3DS::intForce(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &x2, double &y2,
                      double &z2, double &x3, double &y3, double &z3, double &fx0, double &fy0, double &fz0, double &fx1,
                      double &fy1, double &fz1, double &fx2, double &fy2, double &fz2, double &fx3, double &fy3,
                      double &fz3, int &p) {
    return;
}

void HI3DS::gcnForce(double &x0, double &y0, double &z0, double &fx0, double &fy0, double &fz0) {
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

//##### HARMONIC BONDS, RESTRICTED ANGLES, INFINITE COST ARC INTERSECTIONS, CONSTRAINED TO SPHERE #####
HRI3DS::HRI3DS():BasePotentialModel3D(){};

//Potential
double HRI3DS::bndPotential(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, int &p) {
    double dx=x1-x0;
    double dy=y1-y0;
    double dz=z1-z0;
    double r=sqrt((dx*dx+dy*dy+dz*dz));
    return 0.5*bndP[2*p]*pow((r-bndP[2*p+1]),2);
}

double HRI3DS::angPotential(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &x2, double &y2, double &z2, int &p) {
    double dx1=x1-x0;
    double dy1=y1-y0;
    double dz1=z1-z0;
    double dx2=x2-x0;
    double dy2=y2-y0;
    double dz2=z2-z0;
    double n1=1/sqrt(dx1*dx1+dy1*dy1+dz1*dz1);
    double n2=1/sqrt(dx2*dx2+dy2*dy2+dz2*dz2);
    dx1*=n1;
    dy1*=n1;
    dz1*=n1;
    dx2*=n2;
    dy2*=n2;
    dz2*=n2;
    double cosTheta=dx1*dx2+dy1*dy2+dz1*dz2;
    double sinThetaSq=1.0-cosTheta*cosTheta;
    return 0.5*angP[2*p]*pow((cosTheta-angP[2*p+1]),2)/sinThetaSq;
}

double HRI3DS::repPotential(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, int &p) {
    return 0.0;
}

double HRI3DS::intPotential(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &x2,
                           double &y2, double &z2, double &x3, double &y3, double &z3, int &p) {
    VecF<double> p0(3),p1(3),p2(3),p3(3);
    p0[0]=x0;
    p0[1]=y0;
    p0[2]=z0;
    p1[0]=x1;
    p1[1]=y1;
    p1[2]=z1;
    p2[0]=x2;
    p2[1]=y2;
    p2[2]=z2;
    p3[0]=x3;
    p3[1]=y3;
    p3[2]=z3;
    p0/=vNorm(p0);
    p1/=vNorm(p1);
    p2/=vNorm(p2);
    p3/=vNorm(p3);

    VecF<double> v0=crossProduct(p0,p1);
    VecF<double> v1=crossProduct(p2,p3);
    v0/=vNorm(v0);
    v1/=vNorm(v1);

    VecF<double> d=crossProduct(v0,v1);
    VecF<double> s0=d/vNorm(d);
    VecF<double> s1=s0*-1.0;

    bool arc0s0=fabs(orthodromicDist(p0,p1)-orthodromicDist(p0,s0)-orthodromicDist(p1,s0))<1e-8;
    bool arc1s0=fabs(orthodromicDist(p2,p3)-orthodromicDist(p2,s0)-orthodromicDist(p3,s0))<1e-8;
    if(arc0s0 && arc1s0) return numeric_limits<double>::infinity();

    bool arc0s1=fabs(orthodromicDist(p0,p1)-orthodromicDist(p0,s1)-orthodromicDist(p1,s1))<1e-8;
    bool arc1s1=fabs(orthodromicDist(p2,p3)-orthodromicDist(p2,s1)-orthodromicDist(p3,s1))<1e-8;
    if(arc0s1 && arc1s1) return numeric_limits<double>::infinity();

    return 0.0;
}

double HRI3DS::gcnPotential(double &x0, double &y0, double &z0) {

    double r=sqrt((x0*x0+y0*y0+z0*z0));
    return 0.5*gcnP[0]*pow((r-gcnP[1]),2);
}

//Forces
void HRI3DS::bndForce(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1,
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

void HRI3DS::angForce(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &x2, double &y2, double &z2,
                     double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, double &fx2, double &fy2, double &fz2, int &p) {
    double dx1=x1-x0;
    double dy1=y1-y0;
    double dz1=z1-z0;
    double dx2=x2-x0;
    double dy2=y2-y0;
    double dz2=z2-z0;
    double n1=1/sqrt(dx1*dx1+dy1*dy1+dz1*dz1);
    double n2=1/sqrt(dx2*dx2+dy2*dy2+dz2*dz2);
    dx1*=n1;
    dy1*=n1;
    dz1*=n1;
    dx2*=n2;
    dy2*=n2;
    dz2*=n2;
    double cosTheta=dx1*dx2+dy1*dy2+dz1*dz2;
    double sinThetaSq=1.0-cosTheta*cosTheta;
    double ff=-2*angP[2*p]*(cosTheta-angP[2*p+1])*(1-cosTheta*angP[2*p+1])/pow(sinThetaSq,2);
    double ffx1=-ff*n1*(dx2-cosTheta*dx1);
    double ffy1=-ff*n1*(dy2-cosTheta*dy1);
    double ffz1=-ff*n1*(dz2-cosTheta*dz1);
    double ffx2=-ff*n2*(dx1-cosTheta*dx2);
    double ffy2=-ff*n2*(dy1-cosTheta*dy2);
    double ffz2=-ff*n2*(dz1-cosTheta*dz2);
    double ffx0=-ffx1-ffx2;
    double ffy0=-ffy1-ffy2;
    double ffz0=-ffz1-ffz2;
    fx0+=ffx0;
    fy0+=ffy0;
    fz0+=ffz0;
    fx1+=ffx1;
    fy1+=ffy1;
    fz1+=ffz1;
    fx2+=ffx2;
    fy2+=ffy2;
    fz2+=ffz2;
    return;
}

void HRI3DS::repForce(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1,
                     double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, int &p) {
    return;
}

void HRI3DS::intForce(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &x2, double &y2,
                     double &z2, double &x3, double &y3, double &z3, double &fx0, double &fy0, double &fz0, double &fx1,
                     double &fy1, double &fz1, double &fx2, double &fy2, double &fz2, double &fx3, double &fy3,
                     double &fz3, int &p) {
    return;
}

void HRI3DS::gcnForce(double &x0, double &y0, double &z0, double &fx0, double &fy0, double &fz0) {
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


//Additional 3D vector functions
VecF<double> crossProduct(VecF<double>& v0, VecF<double>& v1){
    VecF<double> v2(3);
    v2[0]=v0[1]*v1[2]-v0[2]*v1[1];
    v2[1]=v1[0]*v0[2]-v1[2]*v0[0];
    v2[2]=v0[0]*v1[1]-v0[1]*v1[0];
    return v2;
}

inline double dotProduct(VecF<double>& v0, VecF<double>& v1){
    return v0[0]*v1[0]+v0[1]*v1[1]+v0[2]*v1[2];
}

inline double orthodromicDist(VecF<double>& v0, VecF<double>& v1){
    double num=vNorm(crossProduct(v0,v1));
    double den=dotProduct(v0,v1);
    return atan(num/den);
}
