#include "pot.h"

//##### BASE POTENTIAL MODEL #####

//Default constructor, no interactions
BasePotentialModel::BasePotentialModel() {
    useBnds=false;
    useAngs=false;
    useReps=false;
    useIntx=false;
    useFixd=false;
    useGcns=false;
}

//Setters: interactions types
void BasePotentialModel::setBonds(VecF<int> bonds, VecF<double> params) {
    bnds=bonds;
    bndP=params;
    useBnds=true;
}

void BasePotentialModel::setAngles(VecF<int> angles, VecF<double> params) {
    angs=angles;
    angP=params;
    useAngs=true;
}

void BasePotentialModel::setRepulsions(VecF<int> repulsions, VecF<double> params) {
    reps=repulsions;
    repP=params;
    useReps=true;
}

void BasePotentialModel::setIntersections(VecF<int> intersections, VecF<double> params) {
    intx=intersections;
    intP=params;
    useIntx=true;
}

void BasePotentialModel::setFixedAtoms(VecF<int> fixed) {
    fixd=fixed;
    useFixd=true;
}

void BasePotentialModel::setGeomConstraints(VecF<int> constrained, VecF<double> params) {
    gcns=constrained;
    gcnP=params;
    useGcns=true;
}

void BasePotentialModel::reset() {
    useBnds=false;
    useAngs=false;
    useReps=false;
    useIntx=false;
    useFixd=false;
    useGcns=false;
}

//Potential energy, the objective function
double BasePotentialModel::function(VecF<double>& x) {
    double u=0.0;

    if(useBnds) u+=bndsPotential(x);
    if(useAngs) u+=angsPotential(x);
    if(useReps) u+=repsPotential(x);
    if(useIntx) u+=intxPotential(x);
    if(useGcns) u+=gcnsPotential(x);

    return u;
}

//Negative force, the gradient
VecF<double> BasePotentialModel::gradient(VecF<double> &x) {
    VecF<double> f(x.n);

    if(useBnds) bndsForce(f,x);
    if(useAngs) angsForce(f,x);
    if(useReps) repsForce(f,x);
    if(useIntx) intxForce(f,x);
    if(useFixd) fixdForce(f,x);
    if(useGcns) gcnsForce(f,x);

    return -f;
}

//##### BASE POTENTIAL MODEL 1D #####

//Default constructor
BasePotentialModel1D::BasePotentialModel1D():BasePotentialModel(){};

//Potential energy
double BasePotentialModel1D::bndsPotential(VecF<double> &x) {
    double u=0.0;
    double x0, x1;
    for(int i=0,j=1,k=0;i<bnds.n;i+=2,j+=2,++k){
        x0=x[bnds[i]];
        x1=x[bnds[j]];
        u+=bndPotential(x0,x1,k);
    }
    return u;
}

double BasePotentialModel1D::angsPotential(VecF<double> &x) {
    return 0.0; //cannot have 1D angle
}

double BasePotentialModel1D::repsPotential(VecF<double> &x) {
    double u=0.0;
    double x0, x1;
    for(int i=0,j=1,k=0;i<reps.n;i+=2,j+=2,++k){
        x0=x[reps[i]];
        x1=x[reps[j]];
        u+=repPotential(x0,x1,k);
    }
    return u;
}

double BasePotentialModel1D::intxPotential(VecF<double> &x) {
    return 0.0; //cannot have 1d intersection
}

double BasePotentialModel1D::gcnsPotential(VecF<double> &x) {
    double u=0.0;
    double x0;
    for(int i=0; i<gcns.n; ++i){
        x0=x[gcns[i]];
        u+=gcnPotential(x0);
    }
    return u;
}

//Forces
void BasePotentialModel1D::bndsForce(VecF<double> &f, VecF<double> &x) {

    int id0, id1;
    for(int i=0,j=1,k=0;i<bnds.n;i+=2,j+=2,++k){
        id0=bnds[i];
        id1=bnds[j];
        bndForce(x[id0],x[id1],f[id0],f[id1],k);
    }
}

void BasePotentialModel1D::angsForce(VecF<double> &f, VecF<double> &x) {
    return; //cannot have 1D angles
}

void BasePotentialModel1D::intxForce(VecF<double> &f, VecF<double> &x) {
    return; //cannot have 1D intersections
}

void BasePotentialModel1D::repsForce(VecF<double> &f, VecF<double> &x) {

    int id0, id1;
    for(int i=0,j=1,k=0;i<reps.n;i+=2,j+=2,++k){
        id0=reps[i];
        id1=reps[j];
        repForce(x[id0],x[id1],f[id0],f[id1],k);
    }
}

void BasePotentialModel1D::fixdForce(VecF<double> &f, VecF<double> &x) {

    for(int i=0; i<fixd.n; ++i){
        f[fixd[i]]=0.0;
    }
}

void BasePotentialModel1D::gcnsForce(VecF<double> &f, VecF<double> &x) {

    int id;
    for(int i=0; i<gcns.n; ++i){
        id=gcns[i];
        gcnForce(x[id],f[id]);
    }

}

//##### BASE POTENTIAL MODEL 2D #####

//Default constructor
BasePotentialModel2D::BasePotentialModel2D():BasePotentialModel(){};

//Potentials
double BasePotentialModel2D::bndsPotential(VecF<double> &x) {
    double u=0.0;
    int id0, id1;
    for(int i=0,j=1,k=0;i<bnds.n;i+=2,j+=2,++k){
        id0=bnds[i];
        id1=bnds[j];
        u+=bndPotential(x[2*id0],x[2*id0+1],x[2*id1],x[2*id1+1],k);
    }
    return u;
}

double BasePotentialModel2D::angsPotential(VecF<double> &x) {
    double u=0.0;
    int id0, id1, id2;
    for(int i=0,j=1,k=2,l=0;i<angs.n;i+=3,j+=3,k+=3,++l){
        id0=angs[i];
        id1=angs[j];
        id2=angs[k];
        u+=angPotential(x[2*id0],x[2*id0+1],x[2*id1],x[2*id1+1],x[2*id2],x[2*id2+1],l);
    }
    return u;
}

double BasePotentialModel2D::repsPotential(VecF<double> &x) {
    double u=0.0;
    int id0, id1;
    for(int i=0,j=1,k=0;i<reps.n;i+=2,j+=2,++k){
        id0=reps[i];
        id1=reps[j];
        u+=repPotential(x[2*id0],x[2*id0+1],x[2*id1],x[2*id1+1],k);
    }
    return u;
}

double BasePotentialModel2D::intxPotential(VecF<double> &x) {
    double u=0.0;
    int id0, id1, id2, id3;
    for(int i=0,j=1,k=2,l=3,m=0;i<intx.n;i+=4,j+=4,k+=4,l+=4,++m){
        id0=intx[i];
        id1=intx[j];
        id2=intx[k];
        id3=intx[l];
        u+=intPotential(x[2*id0],x[2*id0+1],x[2*id1],x[2*id1+1],x[2*id2],x[2*id2+1],x[2*id3],x[2*id3+1],m);
    }
    return u;
}

double BasePotentialModel2D::gcnsPotential(VecF<double> &x) {
    double u=0.0;
    int id;
    for(int i=0; i<gcns.n; ++i){
        id=gcns[i];
        u+=gcnPotential(x[2*id],x[2*id+1]);
    }
    return u;
}

//Forces
void BasePotentialModel2D::bndsForce(VecF<double> &f, VecF<double> &x) {
    int id0, id1;
    for(int i=0,j=1,k=0;i<bnds.n;i+=2,j+=2,++k){
        id0=bnds[i];
        id1=bnds[j];
        bndForce(x[2*id0],x[2*id0+1],x[2*id1],x[2*id1+1],
                 f[2*id0],f[2*id0+1],f[2*id1],f[2*id1+1],k);
    }
}

void BasePotentialModel2D::angsForce(VecF<double> &f, VecF<double> &x){
    int id0, id1, id2;
    for(int i=0,j=1,k=2,l=0;i<angs.n;i+=3,j+=3,k+=3,++l){
        id0=angs[i];
        id1=angs[j];
        id2=angs[k];
        angForce(x[2*id0],x[2*id0+1],x[2*id1],x[2*id1+1],x[2*id2],x[2*id2+1],
                 f[2*id0],f[2*id0+1],f[2*id1],f[2*id1+1],f[2*id2],f[2*id2+1],l);
    }
}

void BasePotentialModel2D::repsForce(VecF<double> &f, VecF<double> &x) {
    int id0, id1;
    for(int i=0,j=1,k=0;i<reps.n;i+=2,j+=2,++k){
        id0=reps[i];
        id1=reps[j];
        repForce(x[2*id0],x[2*id0+1],x[2*id1],x[2*id1+1],
                 f[2*id0],f[2*id0+1],f[2*id1],f[2*id1+1],k);
    }
}

void BasePotentialModel2D::intxForce(VecF<double> &f, VecF<double> &x) {
    int id0, id1, id2, id3;
    for(int i=0,j=1,k=2,l=3,m=0;i<intx.n;i+=4,j+=4,k+=4,l+=4,++m){
        id0=intx[i];
        id1=intx[j];
        id2=intx[k];
        id3=intx[l];
        intForce(x[2*id0],x[2*id0+1],x[2*id1],x[2*id1+1],x[2*id2],x[2*id2+1],x[2*id3],x[2*id3+1],
                 f[2*id0],f[2*id0+1],f[2*id1],f[2*id1+1],f[2*id2],f[2*id2+1],f[2*id3],f[2*id3+1],m);
    }
}

void BasePotentialModel2D::fixdForce(VecF<double> &f, VecF<double> &x) {

    for(int i=0; i<fixd.n; ++i){
        f[2*fixd[i]]=0.0;
        f[2*fixd[i]+1]=0.0;
    }
}

void BasePotentialModel2D::gcnsForce(VecF<double> &f, VecF<double> &x) {
    int id;
    for(int i=0; i<gcns.n; ++i){
        id=gcns[i];
        gcnForce(x[2*id],x[2*id+1],f[2*id],f[2*id+1]);
    }
}

//##### BASE POTENTIAL MODEL 3D #####

//Default constructor
BasePotentialModel3D::BasePotentialModel3D():BasePotentialModel(){};

//Potentials
double BasePotentialModel3D::bndsPotential(VecF<double> &x) {
    double u=0.0;
    int id0, id1;
    for(int i=0,j=1,k=0;i<bnds.n;i+=2,j+=2,++k){
        id0=bnds[i];
        id1=bnds[j];
        u+=bndPotential(x[3*id0],x[3*id0+1],x[3*id0+2],x[3*id1],x[3*id1+1],x[3*id1+2],k);
    }
    return u;
}

double BasePotentialModel3D::angsPotential(VecF<double> &x) {
    double u=0.0;
    int id0, id1, id2;
    for(int i=0,j=1,k=2,l=0;i<angs.n;i+=3,j+=3,k+=3,++l){
        id0=angs[i];
        id1=angs[j];
        id2=angs[k];
        u+=angPotential(x[3*id0],x[3*id0+1],x[3*id0+2],x[3*id1],x[3*id1+1],x[3*id1+2],x[3*id2],x[3*id2+1],x[3*id2+2],l);
    }
    return u;
}

double BasePotentialModel3D::repsPotential(VecF<double> &x) {
    double u=0.0;
    int id0, id1;
    for(int i=0,j=1,k=0;i<reps.n;i+=2,j+=2,++k){
        id0=reps[i];
        id1=reps[j];
        u+=repPotential(x[3*id0],x[3*id0+1],x[3*id0+2],x[3*id1],x[3*id1+1],x[3*id1+2],k);
    }
    return u;
}

double BasePotentialModel3D::intxPotential(VecF<double> &x) {
    double u=0.0;
    int id0, id1, id2, id3;
    for(int i=0,j=1,k=2,l=3,m=0;i<intx.n;i+=4,j+=4,k+=4,l+=4,++m){
        id0=intx[i];
        id1=intx[j];
        id2=intx[k];
        id3=intx[l];
        u+=intPotential(x[3*id0],x[3*id0+1],x[3*id0+2],x[3*id1],x[3*id1+1],x[3*id1+2],x[3*id2],x[3*id2+1],x[3*id2+2],x[3*id3],x[3*id3+1],x[3*id3+2],m);
    }
    return u;
}

double BasePotentialModel3D::gcnsPotential(VecF<double> &x) {
    double u=0.0;
    int id;
    for(int i=0; i<gcns.n; ++i){
        id=gcns[i];
        u+=gcnPotential(x[3*id],x[3*id+1],x[3*id+2]);
    }
    return u;
}

//Forces
void BasePotentialModel3D::bndsForce(VecF<double> &f, VecF<double> &x) {
    int id0, id1;
    for(int i=0,j=1,k=0;i<bnds.n;i+=2,j+=2,++k){
        id0=bnds[i];
        id1=bnds[j];
        bndForce(x[3*id0],x[3*id0+1],x[3*id0+2],x[3*id1],x[3*id1+1],x[3*id1+2],
                 f[3*id0],f[3*id0+1],f[3*id0+2],f[3*id1],f[3*id1+1],f[3*id1+2],k);
    }
}

void BasePotentialModel3D::angsForce(VecF<double> &f, VecF<double> &x){
    int id0, id1, id2;
    for(int i=0,j=1,k=2,l=0;i<angs.n;i+=3,j+=3,k+=3,++l){
        id0=angs[i];
        id1=angs[j];
        id2=angs[k];
        angForce(x[3*id0],x[3*id0+1],x[3*id0+2],x[3*id1],x[3*id1+1],x[3*id1+2],x[3*id2],x[3*id2+1],x[3*id2+2],
                 f[3*id0],f[3*id0+1],f[3*id0+2],f[3*id1],f[3*id1+1],f[3*id1+2],f[3*id2],f[3*id2+1],f[3*id2+2],l);
    }
}

void BasePotentialModel3D::repsForce(VecF<double> &f, VecF<double> &x) {
    int id0, id1;
    for(int i=0,j=1,k=0;i<reps.n;i+=2,j+=2,++k){
        id0=reps[i];
        id1=reps[j];
        repForce(x[3*id0],x[3*id0+1],x[3*id0+2],x[3*id1],x[3*id1+1],x[3*id1+2],
                 f[3*id0],f[3*id0+1],f[3*id0+2],f[3*id1],f[3*id1+1],f[3*id1+2],k);
    }
}

void BasePotentialModel3D::intxForce(VecF<double> &f, VecF<double> &x) {
    int id0, id1, id2, id3;
    for(int i=0,j=1,k=2,l=3,m=0;i<intx.n;i+=4,j+=4,k+=4,l+=4,++m){
        id0=intx[i];
        id1=intx[j];
        id2=intx[k];
        id3=intx[l];
        intForce(x[3*id0],x[3*id0+1],x[3*id0+2],x[3*id1],x[3*id1+1],x[3*id1+2],x[3*id2],x[3*id2+1],x[3*id2+2],x[3*id3],x[3*id3+1],x[3*id3+2],
                 f[3*id0],f[3*id0+1],f[3*id0+2],f[3*id1],f[3*id1+1],f[3*id1+2],f[3*id2],f[3*id2+1],f[3*id2+2],f[3*id3],f[3*id3+1],f[3*id3+2],m);
    }
}

void BasePotentialModel3D::fixdForce(VecF<double> &f, VecF<double> &x) {
    for(int i=0; i<fixd.n; ++i){
        f[3*fixd[i]]=0.0;
        f[3*fixd[i]+1]=0.0;
        f[3*fixd[i]+2]=0.0;
    }
}

void BasePotentialModel3D::gcnsForce(VecF<double> &f, VecF<double> &x) {
    int id=0;
    for(int i=0; i<gcns.n; ++i){
        id=gcns[i];
        gcnForce(x[3*id],x[3*id+1],x[3*id+2],f[3*id],f[3*id+1],f[3*id+2]);
    }
}