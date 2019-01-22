#include "linked_network.h"

//Default constructor
LinkedNetwork::LinkedNetwork() {
    //
}

//Construct with defined A network, B generated as the dual of A
LinkedNetwork::LinkedNetwork(int nodesA, string latticeA, int maxACnxs, int maxBCnxs, int minCnxs) {

    //Initialise lattices
    if(latticeA=="square" || latticeA=="triangular"){
        networkA=Network(nodesA,latticeA,maxACnxs);
        networkB=networkA.constructDual(maxBCnxs);
    }
    else if(latticeA=="hexagonal"){
        networkB=Network(nodesA,"triangular",maxBCnxs);
        networkA=networkB.constructDual(maxACnxs);
        rescale(sqrt(3.0));
    }
    else if(latticeA=="cubic"){
        networkA=Network(nodesA,"cubic",maxACnxs);
        networkB=networkA.constructDual(maxBCnxs);
    }
    else if(latticeA=="inv_cubic"){
        networkB=Network(nodesA,"cubic",maxBCnxs);
        networkA=networkB.constructDual(maxACnxs);
    }
    else if(latticeA=="geodesic"){
        networkA=Network(nodesA,"geodesic",maxACnxs);
        networkB=networkA.constructDual(maxBCnxs);
    }
    else if(latticeA=="goldberg"){
        networkB=Network(nodesA,"geodesic",maxBCnxs);
        networkA=networkB.constructDual(maxACnxs);
    }
    minNodeCnxs=minCnxs;
}

//Construct by loading networks from files
LinkedNetwork::LinkedNetwork(string prefix) {
    networkA=Network(prefix+"_A");
    networkB=Network(prefix+"_B");
}

//Set up potential model with single angle and bond parameter set
void LinkedNetwork::initialisePotentialModel(double ak, double bk, double ck, int convexity) {

    //Make copy of lattice A coordinates
    if(networkA.geometryCode=="2DE"){
        crdsA=VecF<double>(2*networkA.nodes.n);
        for(int i=0; i<networkA.nodes.n; ++i){
            crdsA[2*i]=networkA.nodes[i].crd[0];
            crdsA[2*i+1]=networkA.nodes[i].crd[1];
        }
    }
    else{
        crdsA=VecF<double>(3*networkA.nodes.n);
        for(int i=0; i<networkA.nodes.n; ++i){
            crdsA[3*i]=networkA.nodes[i].crd[0];
            crdsA[3*i+1]=networkA.nodes[i].crd[1];
            crdsA[3*i+2]=networkA.nodes[i].crd[2];
        }
    }

    //Initialise potential model parameters
    //Angle parameters
    potParamsA=VecF<double>(4); //for 3 and 4 coordinate
    potParamsA[0]=ak;
    potParamsA[1]=sqrt(3.0);
    potParamsA[2]=ak;
    potParamsA[3]=sqrt(2.0);

    //Bond parameters
    potParamsB=VecF<double>(2);
    potParamsB[0]=bk;
    potParamsB[1]=1.0;

    //Geometry constraint parameters
    potParamsC=VecF<double>(2);
    potParamsC[0]=ck; //k, r0 updated through optimal projection

    //Line intersection parameters
    potParamsD=VecF<int>(2);
    potParamsD[0]=1;
    potParamsD[1]=convexity;

}

//Set up geometry optimisation parameters
void LinkedNetwork::initialiseGeometryOpt(int iterations, double tau, double tolerance, int localExtent) {

    goptParamsA=VecF<int>(2);
    goptParamsA[0]=iterations;
    goptParamsA[1]=localExtent;
    goptParamsB=VecF<double>(2);
    goptParamsB[0]=tau;
    goptParamsB[1]=tolerance;
}

//Set up monte carlo and random number generators
void LinkedNetwork::initialiseMonteCarlo(double temperature, int seed) {

    double energy=globalPotentialEnergy(potParamsD[0],potParamsD[1]);
    mc=Metropolis(seed,temperature,energy);
    mtGen.seed(seed);
}

//Rescale lattice dimensions
void LinkedNetwork::rescale(double scaleFactor) {
    networkA.rescale(scaleFactor);
    networkB.rescale(scaleFactor);
}

//Project lattice onto different geometry
void LinkedNetwork::project(string projType, double param) {
    if(projType=="sphere"){
        networkA.project(projType,param);
        networkB.project(projType,param);
    }
}

//Project lattice onto different geometry with optimal parameters
void LinkedNetwork::optimalProjection(string projType) {

    if(projType=="sphere"){
        //Geometry optimise changing sphere radius until hit minimum in energy

        /* Find initial radius to nearest 1
         * 1) may be multiple minima so try all values in range
         * 2) find lowest value and take limits as radii either side */
        int searchLim=20;
        VecF<double> saveCrdsA=crdsA;
        VecF<double> energies(20);
        energies[0]=numeric_limits<double>::infinity();
        double radius=1.0;
        for(int i=1; i<searchLim; ++i){
            networkA.project(projType,radius);
            potParamsC[1]=radius;
            for(int i=0; i<networkA.nodes.n; ++i){
                crdsA[3*i]=networkA.nodes[i].crd[0];
                crdsA[3*i+1]=networkA.nodes[i].crd[1];
                crdsA[3*i+2]=networkA.nodes[i].crd[2];
            }
            globalGeometryOptimisation(false,false);
            energies[i]=globalPotentialEnergy(false,false);
            cout<<radius<<" "<<energies[i]<<endl;
            radius+=1.0;
            crdsA=saveCrdsA;
        }
        int id0;
        double e0=energies[0],e1;
        double lowerLim, upperLim, minRadius, minEnergy;
        for(int i=1; i<searchLim; ++i){
            if(energies[i]<e0){
                e0=energies[i];
                id0=i;
            }
        }
        if(id0==searchLim-1) throw "Initial spherical minimisation reached search limit";
        else{
            lowerLim=id0-1.0;
            minRadius=id0;
            upperLim=id0+1.0;
            networkA.project(projType,minRadius);
            potParamsC[1]=minRadius;
            for(int i=0; i<networkA.nodes.n; ++i){
                crdsA[3*i]=networkA.nodes[i].crd[0];
                crdsA[3*i+1]=networkA.nodes[i].crd[1];
                crdsA[3*i+2]=networkA.nodes[i].crd[2];
            }
            saveCrdsA=crdsA;
            globalGeometryOptimisation(false,false);
            minEnergy=globalPotentialEnergy(false,false);
        }

        /* Refine minimimum
         * 1) could have multiple minima if small - hence reset coordinates as can cause instability
         * 2) search between lower and upper limits until pass through minimum
         * 3) decrease search increment and search again */
        double radiusInc=0.1;
        for(int i=0; i<3; ++i){
            radius=lowerLim;
            e0=numeric_limits<double>::infinity();
            for(;;){
                networkA.project(projType,radius);
                for(int i=0; i<networkA.nodes.n; ++i){
                    crdsA[3*i]=networkA.nodes[i].crd[0];
                    crdsA[3*i+1]=networkA.nodes[i].crd[1];
                    crdsA[3*i+2]=networkA.nodes[i].crd[2];
                }
                potParamsC[1]=radius;
                globalGeometryOptimisation(false,false);
                e1=globalPotentialEnergy(false,false);
                cout<<radius<<" "<<e1<<" "<<minEnergy<<endl;
                if(e1>e0 && e0<=minEnergy){
                    lowerLim=radius-2*radiusInc;
                    minRadius=radius-radiusInc;
                    upperLim=radius;
                    minEnergy=e0;
                    break;
                }
                else e0=e1;
                radius+=radiusInc;
                if(radius>upperLim) throw "Could not find minimum in initial spherical minimisation";
            }
            radiusInc/=10.0;
        }

        //Optimise with minimum radius
        potParamsC[1]=minRadius;
        globalGeometryOptimisation(false,false);
    }
}

//Select nodes forming a random edge in lattice A, and linked nodes in lattice B, only for 3/4 coordinate nodes
int LinkedNetwork::randomCnx34(int& a, int& b, int& u, int& v, mt19937& gen) {

    //pick random node and one of its random connections
    uniform_int_distribution<int> randomNode(0,networkA.nodes.n-1);
    int n0=randomNode(gen);
    int cnd0=networkA.nodes[n0].netCnxs.n;
    uniform_int_distribution<int> randomCnx(0,cnd0-1);
    int n1=networkA.nodes[n0].netCnxs[randomCnx(gen)];
    int cnd1=networkA.nodes[n1].netCnxs.n;

    //find connection type and assign a,b
    int cnxType;
    if(cnd0==3 && cnd1==3){
        cnxType=33;
        a=n0;
        b=n1;
    }
    else if (cnd0==4 && cnd1==4){
        cnxType=44;
        a=n0;
        b=n1;
    }
    else if (cnd0==3 && cnd1==4){
        cnxType=34;
        a=n0;
        b=n1;
    }
    else if (cnd0==4 && cnd1==3){
        cnxType=34;
        a=n1;
        b=n0;
    }
    else throw "Error in random connection - incorrect coordinations";

    //get nodes in dual in random orientation
    uniform_int_distribution<int> randomDirection(0,1);
    VecR<int> common=vCommonValues(networkA.nodes[a].dualCnxs,networkA.nodes[b].dualCnxs);
    if(common.n!=2) throw "Error in random connection - incorrect dual ids";
    int randIndex=randomDirection(gen);
    u=common[randIndex];
    v=common[1-randIndex];

    return cnxType;
}

//Generate all ids of nodes in lattices A and B needed for switch move, only for 3/4 coordinate nodes
int LinkedNetwork::generateSwitchIds34(int cnxType, VecF<int> &switchIdsA, VecF<int> &switchIdsB, int a, int b, int u, int v) {

    //lots of error checking to remove any potential pathological cases
    if(a==b || u==v){
        cout<<"Note: skip in switch generation"<<endl;
        return 1;
    }

    if(cnxType==33){
        /* Switch connectivities in lattice and dual
         * 3-3 coordination connection
         * a,b,c,d,e,f are nodes in lattice A
         * u,v,w,x are nodes in lattice B
         * a-b, a-c, a-e
         * b-a, b-d, b-f
         * a-b share u-v
         * c-a-b-d share u
         * e-a-b-f share v
         * c-a-e share w
         * d-b-f share x
         * u-v, u-w, u-x
         * v-u, v-w, v-x */
        int errorFlag=0;
        int c,d,e,f;
        int w,x;
        VecR<int> common,common1;
        common=vCommonValues(networkA.nodes[a].netCnxs,networkB.nodes[u].dualCnxs);
        common.delValue(b);
        if(common.n!=1) errorFlag=1;
        c=common[0];
        common=vCommonValues(networkA.nodes[b].netCnxs,networkB.nodes[u].dualCnxs);
        common.delValue(a);
        if(common.n!=1) errorFlag=2;
        d=common[0];
        common=vCommonValues(networkA.nodes[a].netCnxs,networkB.nodes[v].dualCnxs);
        common.delValue(b);
        if(common.n!=1) errorFlag=3;
        e=common[0];
        common=vCommonValues(networkA.nodes[b].netCnxs,networkB.nodes[v].dualCnxs);
        common.delValue(a);
        if(common.n!=1) errorFlag=4;
        f=common[0];

        common=vCommonValues(networkA.nodes[a].dualCnxs,networkA.nodes[c].dualCnxs);
        common1=vCommonValues(networkA.nodes[a].dualCnxs,networkA.nodes[e].dualCnxs);
        common.delValue(u);
        common1.delValue(v);
        if(common.n!=1 || common1.n!=1 || common[0]!=common1[0]) errorFlag=5;
        w=common[0];
        common=vCommonValues(networkA.nodes[b].dualCnxs,networkA.nodes[d].dualCnxs);
        common1=vCommonValues(networkA.nodes[b].dualCnxs,networkA.nodes[f].dualCnxs);
        common.delValue(u);
        common1.delValue(v);
        if(common.n!=1 || common1.n!=1 || common[0]!=common1[0]) errorFlag=5;
        x=common[0];

        if(errorFlag!=0){
            cout<<"Note: skip in switch generation"<<endl;
            return 1;
        }

        //check move will not violate dual connectivity limits
        if(networkB.nodes[u].netCnxs.n==minNodeCnxs || networkB.nodes[v].netCnxs.n==minNodeCnxs
            || networkB.nodes[w].netCnxs.n==networkB.nodes[w].netCnxs.nMax
            || networkB.nodes[x].netCnxs.n==networkB.nodes[x].netCnxs.nMax) return 1;
        else {
            switchIdsA = VecF<int>(6);
            switchIdsB = VecF<int>(4);
            switchIdsA[0] = a;
            switchIdsA[1] = b;
            switchIdsA[2] = c;
            switchIdsA[3] = d;
            switchIdsA[4] = e;
            switchIdsA[5] = f;
            switchIdsB[0] = u;
            switchIdsB[1] = v;
            switchIdsB[2] = w;
            switchIdsB[3] = x;
            return 0;
        }
    }
    else if(cnxType==34 || cnxType==44) throw "Not yet implemented!";
    return 0;


//    /* Switch connectivities in lattice and dual
//     * a,b,c,d,e,f are nodes in lattice A
//     * u,v,w,x,y,z are nodes in lattice B
//     * a-b share u-v
//     * a-c share w-u
//     * b-d share u-y
//     * a-e share x-v
//     * b-f share v,z
//     * if 3-coordinate lattice set x==v, y==u
//     * works for 3 or 4 coordinate lattice */
//
//    //find all a-f and u-z
//    int c,d,e,f;
//    int w,x,y,z;
//    VecR<int> common;
//    if(u==-1 || v==-1) {//optional whether u,v, supplied
//        common = vCommonValues(networkA.nodes[a].dualCnxs, networkA.nodes[b].dualCnxs);
//        u = common[0];
//        v = common[1];
//    }
//    common=vCommonValues(networkB.nodes[u].dualCnxs,networkA.nodes[a].netCnxs);
//    common.delValue(b);
//    c=common[0];
//    common=vCommonValues(networkB.nodes[u].dualCnxs,networkA.nodes[b].netCnxs);
//    common.delValue(a);
//    d=common[0];
//    common=vCommonValues(networkB.nodes[v].dualCnxs,networkA.nodes[a].netCnxs);
//    common.delValue(b);
//    e=common[0];
//    common=vCommonValues(networkB.nodes[v].dualCnxs,networkA.nodes[b].netCnxs);
//    common.delValue(a);
//    f=common[0];
//    common=vCommonValues(networkA.nodes[a].dualCnxs,networkA.nodes[c].dualCnxs);
//    common.delValue(u);
//    w=common[0];
//    common=vCommonValues(networkA.nodes[b].dualCnxs,networkA.nodes[f].dualCnxs);
//    common.delValue(v);
//    z=common[0];
//    common=vCommonValues(networkA.nodes[a].dualCnxs,networkA.nodes[e].dualCnxs);
//    common.delValue(v);
//    x=common[0];
//    if(x==w) x=v;
//    common=vCommonValues(networkA.nodes[b].dualCnxs,networkA.nodes[d].dualCnxs);
//    common.delValue(u);
//    y=common[0];
//    if(y==z) y=u;
//
//    //check move will not violate dual connectivity limits
//    if(networkB.nodes[u].netCnxs.n==minNodeCnxs || networkB.nodes[v].netCnxs.n==minNodeCnxs
//       || networkB.nodes[w].netCnxs.n==networkB.nodes[w].netCnxs.nMax
//       || networkB.nodes[z].netCnxs.n==networkB.nodes[z].netCnxs.nMax) return 1;
//    else{
//        switchIdsA=VecF<int>(6);
//        switchIdsB=VecF<int>(6);
//        switchIdsA[0]=a;
//        switchIdsA[1]=b;
//        switchIdsA[2]=c;
//        switchIdsA[3]=d;
//        switchIdsA[4]=e;
//        switchIdsA[5]=f;
//        switchIdsB[0]=u;
//        switchIdsB[1]=v;
//        switchIdsB[2]=w;
//        switchIdsB[3]=x;
//        switchIdsB[4]=y;
//        switchIdsB[5]=z;
//        return 0;
//    }
}

//Switch connectivities in lattice between 2x3 coordinate nodes
void LinkedNetwork::switchCnx33(VecF<int> switchIdsA, VecF<int> switchIdsB) {

    //unpck parameters
    int a,b,c,d,e,f;
    int u,v,w,x;
    a=switchIdsA[0];
    b=switchIdsA[1];
    c=switchIdsA[2];
    d=switchIdsA[3];
    e=switchIdsA[4];
    f=switchIdsA[5];
    u=switchIdsB[0];
    v=switchIdsB[1];
    w=switchIdsB[2];
    x=switchIdsB[3];

    //Apply changes to descriptors due to breaking connections
    //For network A node distribution and edge distribution will remain unchanged

    //For network B node and edge distribution will change
    int nu, nv, nw, nx;
    nu=networkB.nodes[u].netCnxs.n;
    nv=networkB.nodes[v].netCnxs.n;
    nw=networkB.nodes[w].netCnxs.n;
    nx=networkB.nodes[x].netCnxs.n;
    --networkB.nodeDistribution[nu];
    --networkB.nodeDistribution[nv];
    --networkB.nodeDistribution[nw];
    --networkB.nodeDistribution[nx];
    for(int i=0; i<nu; ++i){
        int id=networkB.nodes[u].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nu][nCnx];
        if(id!=v && id!=w && id!=x) --networkB.edgeDistribution[nCnx][nu]; //prevent double counting
    }
    for(int i=0; i<nv; ++i){
        int id=networkB.nodes[v].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nv][nCnx];
        if(id!=u && id!=w && id!=x) --networkB.edgeDistribution[nCnx][nv]; //prevent double counting
    }
    for(int i=0; i<nw; ++i){
        int id=networkB.nodes[w].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nw][nCnx];
        if(id!=u && id!=v && id!=x) --networkB.edgeDistribution[nCnx][nw]; //prevent double counting
    }
    for(int i=0; i<nx; ++i){
        int id=networkB.nodes[x].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nx][nCnx];
        if(id!=u && id!=v && id!=w) --networkB.edgeDistribution[nCnx][nx]; //prevent double counting
    }

    //A-A connectivities
    //swap a->e/d, b->d/e, d->b/a, e->a/b
    networkA.nodes[a].netCnxs.swapValue(e,d);
    networkA.nodes[b].netCnxs.swapValue(d,e);
    networkA.nodes[d].netCnxs.swapValue(b,a);
    networkA.nodes[e].netCnxs.swapValue(a,b);

    //A-B connectvities
    //swap a->v/x, b->u/w
    networkA.nodes[a].dualCnxs.swapValue(v,x);
    networkA.nodes[b].dualCnxs.swapValue(u,w);

    //B-B connectivities
    //have to account for the fact that two nodes may connect multiple times
    //break insert w:u-(x)-v, x:u-(w)-v
    networkB.nodes[u].netCnxs.swapValue(v,-1,w,x);
    networkB.nodes[v].netCnxs.swapValue(u,-1,w,x);
    networkB.nodes[u].netCnxs.delValue(-1);
    networkB.nodes[v].netCnxs.delValue(-1);
    networkB.nodes[w].netCnxs.insertValue(x,u,v);
    networkB.nodes[x].netCnxs.insertValue(w,u,v);

    //B-A connectivities
    //break u->b, v->a, insert w:a-(b)-e, z:b-(a)-d
    networkB.nodes[u].dualCnxs.delValue(b);
    networkB.nodes[v].dualCnxs.delValue(a);
    networkB.nodes[w].dualCnxs.insertValue(b,a,e);
    networkB.nodes[x].dualCnxs.insertValue(a,b,d);

    //Apply changes to descriptors due to making connections
    //Network B
    nu=networkB.nodes[u].netCnxs.n;
    nv=networkB.nodes[v].netCnxs.n;
    nw=networkB.nodes[w].netCnxs.n;
    nx=networkB.nodes[x].netCnxs.n;
    ++networkB.nodeDistribution[nu];
    ++networkB.nodeDistribution[nv];
    ++networkB.nodeDistribution[nw];
    ++networkB.nodeDistribution[nx];
    for(int i=0; i<nu; ++i){
        int id=networkB.nodes[u].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nu][nCnx];
        if(id!=v && id!=w && id!=x) ++networkB.edgeDistribution[nCnx][nu]; //prevent double counting
    }
    for(int i=0; i<nv; ++i){
        int id=networkB.nodes[v].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nv][nCnx];
        if(id!=u && id!=w && id!=x) ++networkB.edgeDistribution[nCnx][nv]; //prevent double counting
    }
    for(int i=0; i<nw; ++i){
        int id=networkB.nodes[w].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nw][nCnx];
        if(id!=u && id!=v && id!=x) ++networkB.edgeDistribution[nCnx][nw]; //prevent double counting
    }
    for(int i=0; i<nx; ++i){
        int id=networkB.nodes[x].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nx][nCnx];
        if(id!=u && id!=v && id!=w) ++networkB.edgeDistribution[nCnx][nx]; //prevent double counting
    }
}

//Swtich connectivities in lattice
void LinkedNetwork::switchCnx(VecF<int> switchIdsA, VecF<int> switchIdsB) {

    //unpack parameters
    int a,b,c,d,e,f;
    int u,v,w,x,y,z;
    a=switchIdsA[0];
    b=switchIdsA[1];
    c=switchIdsA[2];
    d=switchIdsA[3];
    e=switchIdsA[4];
    f=switchIdsA[5];
    u=switchIdsB[0];
    v=switchIdsB[1];
    w=switchIdsB[2];
    x=switchIdsB[3];
    y=switchIdsB[4];
    z=switchIdsB[5];

    //Apply changes to descriptors due to breaking connections
    //For network A node distribution will remain unchanged on remaking connections
    //Edge distribution will only change if there is a mixed coordination system
    int na, nb, nc, nf;
    na=networkA.nodes[a].netCnxs.n;
    nb=networkA.nodes[b].netCnxs.n;
    nc=networkA.nodes[c].netCnxs.n;
    nf=networkA.nodes[f].netCnxs.n;
    --networkA.edgeDistribution[na][nc];
    --networkA.edgeDistribution[nc][na];
    --networkA.edgeDistribution[nb][nf];
    --networkA.edgeDistribution[nf][nb];
    //For network B node and edge distribution will change
    int nu, nv, nw, nz;
    nu=networkB.nodes[u].netCnxs.n;
    nv=networkB.nodes[v].netCnxs.n;
    nw=networkB.nodes[w].netCnxs.n;
    nz=networkB.nodes[z].netCnxs.n;
    --networkB.nodeDistribution[nu];
    --networkB.nodeDistribution[nv];
    --networkB.nodeDistribution[nw];
    --networkB.nodeDistribution[nz];
    for(int i=0; i<nu; ++i){
        int id=networkB.nodes[u].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nu][nCnx];
        if(id!=v && id!=w && id!=z) --networkB.edgeDistribution[nCnx][nu]; //prevent double counting
    }
    for(int i=0; i<nv; ++i){
        int id=networkB.nodes[v].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nv][nCnx];
        if(id!=u && id!=w && id!=z) --networkB.edgeDistribution[nCnx][nv]; //prevent double counting
    }
    for(int i=0; i<nw; ++i){
        int id=networkB.nodes[w].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nw][nCnx];
        if(id!=u && id!=v && id!=z) --networkB.edgeDistribution[nCnx][nw]; //prevent double counting
    }
    for(int i=0; i<nz; ++i){
        int id=networkB.nodes[z].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nz][nCnx];
        if(id!=u && id!=v && id!=w) --networkB.edgeDistribution[nCnx][nz]; //prevent double counting
    }

    //A-A connectivities
    //break a->c, b->f, insert a: b-(f)-e, b: a-(c)-d
    networkA.nodes[a].netCnxs.delValue(c);
    networkA.nodes[b].netCnxs.delValue(f);
    networkA.nodes[a].netCnxs.insertValue(f,b,e);
    networkA.nodes[b].netCnxs.insertValue(c,a,d);
    //swap c->a/b, f->b/a
    networkA.nodes[c].netCnxs.swapValue(a,b);
    networkA.nodes[f].netCnxs.swapValue(b,a);

    //A-B connectvities
    //swap a: u->z, b: v->w
    networkA.nodes[a].dualCnxs.swapValue(u,z);
    networkA.nodes[b].dualCnxs.swapValue(v,w);

    //B-B connectivities
    //break u<->v, insert w: u-(z)-x, z: y-(w)-v
    networkB.nodes[u].netCnxs.delValue(v);
    networkB.nodes[v].netCnxs.delValue(u);
    networkB.nodes[w].netCnxs.insertValue(z,u,x);
    networkB.nodes[z].netCnxs.insertValue(w,y,v);

    //B-A connectivities
    //break u-a, v-b, insert w: a-(b)-c, z: b-(a)-f
    networkB.nodes[u].dualCnxs.delValue(a);
    networkB.nodes[v].dualCnxs.delValue(b);
    networkB.nodes[w].dualCnxs.insertValue(b,a,c);
    networkB.nodes[z].dualCnxs.insertValue(a,b,f);

    //Apply changes to descriptors due to making connections
    //Network A
    na=networkA.nodes[a].netCnxs.n;
    nb=networkA.nodes[b].netCnxs.n;
    nc=networkA.nodes[c].netCnxs.n;
    nf=networkA.nodes[f].netCnxs.n;
    ++networkA.edgeDistribution[na][nf];
    ++networkA.edgeDistribution[nf][na];
    ++networkA.edgeDistribution[nb][nc];
    ++networkA.edgeDistribution[nc][nb];
    //Network B
    nu=networkB.nodes[u].netCnxs.n;
    nv=networkB.nodes[v].netCnxs.n;
    nw=networkB.nodes[w].netCnxs.n;
    nz=networkB.nodes[z].netCnxs.n;
    ++networkB.nodeDistribution[nu];
    ++networkB.nodeDistribution[nv];
    ++networkB.nodeDistribution[nw];
    ++networkB.nodeDistribution[nz];
    for(int i=0; i<nu; ++i){
        int id=networkB.nodes[u].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nu][nCnx];
        if(id!=v && id!=w && id!=z) ++networkB.edgeDistribution[nCnx][nu]; //prevent double counting
    }
    for(int i=0; i<nv; ++i){
        int id=networkB.nodes[v].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nv][nCnx];
        if(id!=u && id!=w && id!=z) ++networkB.edgeDistribution[nCnx][nv]; //prevent double counting
    }
    for(int i=0; i<nw; ++i){
        int id=networkB.nodes[w].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nw][nCnx];
        if(id!=u && id!=v && id!=z) ++networkB.edgeDistribution[nCnx][nw]; //prevent double counting
    }
    for(int i=0; i<nz; ++i){
        int id=networkB.nodes[z].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nz][nCnx];
        if(id!=u && id!=v && id!=w) ++networkB.edgeDistribution[nCnx][nz]; //prevent double counting
    }
}

//Single monte carlo switching move
VecF<int> LinkedNetwork::monteCarloSwitchMove(double& energy) {

    /* Single MC switch move
     * 1) select random connection
     * 2) switch connection
     * 3) optimise and evaluate energy
     * 4) accept or reject */

    //Select valid random connection - that will not violate connection limits
    int a,b,u,v;
    VecF<int> switchIdsA, switchIdsB;
    int validMove;
    int cnxType;
    for(int i=0; i<networkA.nodes.n*networkA.nodes.n; ++i){//catch in case cannot find any valid moves
        cnxType=randomCnx34(a,b,u,v,mtGen);
        validMove=generateSwitchIds34(cnxType,switchIdsA,switchIdsB,a,b,u,v);
        if(validMove==0) break;
    }
    if(validMove==1) throw "Cannot find any valid switch moves";

    //Save current state
    double saveEnergy=energy;
    VecF<double> saveCrdsA=crdsA;
    VecF<int> saveNodeDistA=networkA.nodeDistribution;
    VecF<int> saveNodeDistB=networkB.nodeDistribution;
    VecF< VecF<int> > saveEdgeDistA=networkA.edgeDistribution;
    VecF< VecF<int> > saveEdgeDistB=networkB.edgeDistribution;
    VecF<Node> saveNodesA(switchIdsA.n), saveNodesB(switchIdsB.n);
    for(int i=0; i<saveNodesA.n; ++i) saveNodesA[i]=networkA.nodes[switchIdsA[i]];
    for(int i=0; i<saveNodesB.n; ++i) saveNodesB[i]=networkB.nodes[switchIdsB[i]];

    //Switch and geometry optimise
    VecF<int> optStatus;
    if(cnxType==33) switchCnx33(switchIdsA,switchIdsB);
    else throw "Not yet implemented!";
    localGeometryOptimisation(a,b,1,false,false); //bond switch atoms only
    optStatus=localGeometryOptimisation(a,b,goptParamsA[1],potParamsD[0],potParamsD[1]); //wider area

    //Accept or reject
    energy=globalPotentialEnergy(potParamsD[0],potParamsD[1]);
    int accept=mc.acceptanceCriterion(energy);
//    accept=1;
    if(accept==0){
        energy=saveEnergy;
        crdsA=saveCrdsA;
        networkA.nodeDistribution=saveNodeDistA;
        networkA.edgeDistribution=saveEdgeDistA;
        networkB.nodeDistribution=saveNodeDistB;
        networkB.edgeDistribution=saveEdgeDistB;
        for(int i=0; i<saveNodesA.n; ++i) networkA.nodes[switchIdsA[i]]=saveNodesA[i];
        for(int i=0; i<saveNodesB.n; ++i) networkB.nodes[switchIdsB[i]]=saveNodesB[i];
    }

    /* Status report
     * [0] accepted/rejected 1/0
     * [1] optimisation code 0=successful 1=successful(zero force) 2=unsuccessful(it limit) 3=unsuccessful(intersection)
     * [2] optimisation iterations */
    VecF<int> status(3);
    status[0]=accept;
    status[1]=optStatus[0];
    status[2]=optStatus[1];

    return status;
}


//Calculate potential energy of entire system
double LinkedNetwork::globalPotentialEnergy(bool useIntx, bool keepConvex) {

    /* Potential model
     * Bonds as harmonic
     * Angles as approximated harmonic
     * Local line intersections */

    //Set up potential model
    //Bond and angle harmonics
    VecR<int> bonds(0,networkA.nodes.n*100);
    VecR<double> bondParams(0,networkA.nodes.n*100);
    for(int i=0; i<networkA.nodes.n; ++i){
        generateHarmonics(i,bonds,bondParams);
    }
    //Intersections
    VecR<int> intersections(0,networkA.nodes.n*100);
    if(useIntx){
        for(int i=0; i<networkA.nodes.n; ++i){
            generateIntersections(i,intersections,keepConvex);
        }
    }

    //Assign to fixed size arrays
    VecF<int> bnds(bonds.n),intx(intersections.n);
    VecF<double> bndP(bondParams.n),intxP; //placeholder for intersections
    for(int i=0; i<bnds.n; ++i) bnds[i]=bonds[i];
    for(int i=0; i<bndP.n; ++i) bndP[i]=bondParams[i];
    for(int i=0; i<intersections.n; ++i) intx[i]=intersections[i];

    //Potential model based on geometry code
    double potEnergy;
    if(networkA.geometryCode=="2DE"){
        HI2DP potModel(networkA.pb[0],networkA.pb[1]);
        potModel.setBonds(bnds,bndP);
        if(useIntx) potModel.setIntersections(intx,intxP);
        potEnergy=potModel.function(crdsA);
    }
    else if(networkA.geometryCode=="2DS"){
        VecF<int> constrained(networkA.nodes.n);
        for(int i=0; i<networkA.nodes.n; ++i) constrained[i]=i;
        HI3DS potModel;
        potModel.setBonds(bnds,bndP);
        potModel.setGeomConstraints(constrained,potParamsC);
        if(useIntx) potModel.setIntersections(intx,intxP);
        potEnergy=potModel.function(crdsA);
    }

    return potEnergy;
}

//Geometry optimise entire system
void LinkedNetwork::globalGeometryOptimisation(bool useIntx, bool keepConvex) {

    /* Potential model
     * Bonds as harmonic
     * Angles as approximated harmonic
     * Local line intersections */

    //Set up potential model
    //Bond and angle harmonics
    VecR<int> bonds(0,networkA.nodes.n*100);
    VecR<double> bondParams(0,networkA.nodes.n*100);
    for(int i=0; i<networkA.nodes.n; ++i){
        generateHarmonics(i,bonds,bondParams);
    }
    //Intersections
    VecR<int> intersections(0,networkA.nodes.n*100);
    if(useIntx){
        for(int i=0; i<networkA.nodes.n; ++i){
            generateIntersections(i,intersections,keepConvex);
        }
    }

    //Assign to fixed size arrays
    VecF<int> bnds(bonds.n),intx(intersections.n);
    VecF<double> bndP(bondParams.n),intxP; //placeholder for intersections
    for(int i=0; i<bnds.n; ++i) bnds[i]=bonds[i];
    for(int i=0; i<bndP.n; ++i) bndP[i]=bondParams[i];
    for(int i=0; i<intersections.n; ++i) intx[i]=intersections[i];

    //Geometry optimise
    if(networkA.geometryCode=="2DE"){
        HI2DP potModel(networkA.pb[0],networkA.pb[1]);
        potModel.setBonds(bnds,bndP);
        if(useIntx) potModel.setIntersections(intx,intxP);
        if(potModel.function(crdsA)<numeric_limits<double>::infinity()){//only optimise if no line intersections
            SteepestDescentArmijo<HI2DP> optimiser(goptParamsA[0],goptParamsB[0],goptParamsB[1]);
            optimiser(potModel,crdsA);
        }
    }
    else if(networkA.geometryCode=="2DS"){
        VecF<int> constrained(networkA.nodes.n);
        for(int i=0; i<networkA.nodes.n; ++i) constrained[i]=i;
        HI3DS potModel;
        potModel.setBonds(bnds,bndP);
        potModel.setGeomConstraints(constrained,potParamsC);
        if(useIntx) potModel.setIntersections(intx,intxP);
        if(potModel.function(crdsA)<numeric_limits<double>::infinity()){//only optimise if no line intersections
            SteepestDescentArmijo<HI3DS> optimiser(goptParamsA[0], goptParamsB[0], goptParamsB[1]);
            optimiser(potModel, crdsA);
        }
    }

}

//Geometry optimise subsection of system by only including interactions in a specified range
VecF<int> LinkedNetwork::localGeometryOptimisation(int centreA, int centreB, int extent, bool useIntx, bool keepConvex) {

    /* Find three regions
     * 1) local (full interactions)
     * 2) fixed inner (interactions with local and fixed but immobile)
     * 3) fixed outer (interactions with fixed inner but immobile) */
    VecR<int> local, fixedInner, fixedOuter;
    networkA.findLocalRegion(centreA,centreB,extent,local,fixedInner,fixedOuter);

    //Harmonics
    VecR<int> bonds(0,local.n*100);
    VecR<double> bondParams(0,local.n*100);
    for(int i=0; i<local.n; ++i){
        generateHarmonics(local[i],bonds,bondParams);
    }
    for(int i=0; i<fixedInner.n; ++i){
        generateHarmonics(fixedInner[i],bonds,bondParams);
    }

    //Intersections
    VecR<int> intersections(0,local.n*1000);
    if(useIntx) {
        for(int i=0; i<local.n; ++i){
            generateIntersections(local[i],intersections,keepConvex);
        }
        for(int i=0; i<fixedInner.n; ++i){
            generateIntersections(fixedInner[i],intersections,keepConvex);
        }
    }

    //Assign to fixed size arrays
    VecF<int> bnds(bonds.n),fixd(fixedInner.n+fixedOuter.n),intx(intersections.n);
    VecF<double> bndP(bondParams.n),intxP; //placeholder for intersections
    for(int i=0; i<bnds.n; ++i) bnds[i]=bonds[i];
    for(int i=0; i<bndP.n; ++i) bndP[i]=bondParams[i];
    for(int i=0; i<fixedInner.n; ++i) fixd[i]=fixedInner[i];
    for(int i=0; i<fixedOuter.n; ++i) fixd[i+fixedInner.n]=fixedOuter[i];
    for(int i=0; i<intersections.n; ++i) intx[i]=intersections[i];

    //Geometry optimise
    VecF<int> optStatus(2);
    if(networkA.geometryCode=="2DE"){
        HI2DP potModel(networkA.pb[0],networkA.pb[1]);
        potModel.setBonds(bnds,bndP);
        potModel.setFixedAtoms(fixd);
        if(useIntx) potModel.setIntersections(intx,intxP);
        if(potModel.function(crdsA)<numeric_limits<double>::infinity()){//optimise if no line intersections
            SteepestDescentArmijo<HI2DP> optimiser(goptParamsA[0],goptParamsB[0],goptParamsB[1]);
            optStatus=optimiser(potModel,crdsA);
        }
        else{
            optStatus[0]=3;
            optStatus[1]=0;
        }
    }
    else if(networkA.geometryCode=="2DS"){
        VecF<int> constrained(networkA.nodes.n);
        for(int i=0; i<networkA.nodes.n; ++i) constrained[i]=i;
        HI3DS potModel;
        potModel.setBonds(bnds,bndP);
        potModel.setFixedAtoms(fixd);
        potModel.setGeomConstraints(constrained,potParamsC);
        if(useIntx) potModel.setIntersections(intx,intxP);
        if(potModel.function(crdsA)<numeric_limits<double>::infinity()) {//optimise if no line intersections
            SteepestDescentArmijo<HI3DS> optimiser(goptParamsA[0], goptParamsB[0], goptParamsB[1]);
            optStatus = optimiser(potModel, crdsA);
        }
        else{
            optStatus[0]=3;
            optStatus[1]=0;
        }
    }
    return optStatus;
}


//Generate harmonic potentials corresponding to bonds and angles associated with a given node in lattice A
void LinkedNetwork::generateHarmonics(int id, VecR<int>& bonds, VecR<double>& bondParams) {

    //Potential parameters
    double bk=potParamsB[0], br=potParamsB[1]; //bond k and r0
    double ak, ar; //angle k and r0

    //Harmonics
    int cnd; //coordination
    int id0 = id, id1, id2;
    cnd = networkA.nodes[id0].netCnxs.n;
    if (cnd == 3) {
        ak = potParamsA[0];
        ar = potParamsA[1];
    } else if (cnd == 4) {
        ak = potParamsA[2];
        ar = potParamsA[3];
    }
    for (int i = 0; i < cnd; ++i) {
        id1 = networkA.nodes[id0].netCnxs[i];
        id2 = networkA.nodes[id0].netCnxs[(i + 1) % cnd];
        //bonds
        if (id0 < id1) {//prevent double counting
            bonds.addValue(id0);
            bonds.addValue(id1);
            bondParams.addValue(bk);
            bondParams.addValue(br);
        }
        //angles
        bonds.addValue(id1);
        bondParams.addValue(ak);
        bonds.addValue(id2);
        bondParams.addValue(ar);
    }
}

//Generate intersections associated with a given node from lattice A
void LinkedNetwork::generateIntersections(int id, VecR<int> &intersections, bool keepConvex) {

    //Maintain convexity
    int cnd=networkA.nodes[id].netCnxs.n;
    int id0,id1,id2;
    if(keepConvex) {
        for (int i = 0; i < cnd; ++i) {
            id0 = networkA.nodes[id].netCnxs[i];
            id1 = networkA.nodes[id].netCnxs[(i + 1) % cnd];
            id2 = networkA.nodes[id].netCnxs[(i + 2) % cnd];
            intersections.addValue(id0);
            intersections.addValue(id1);
            intersections.addValue(id);
            intersections.addValue(id2);
        }
    }

    /* Create unique intersections with local environment
     * treat 3/4 coordinate separately to optimise*/
    int idA0,idA1,idA2,idA3,idA4,idB;
    int l0a,l0b,l1a,l1b,l2a,l2b;
    idA0=id;
    cnd=networkA.nodes[idA0].netCnxs.n;
    if(cnd==3) {
        for (int j = 0; j < cnd; ++j) {
            idA1 = networkA.nodes[idA0].netCnxs[j];
            idA2 = networkA.nodes[idA0].netCnxs[(j + 1) % cnd];
            idB = vCommonValues(networkA.nodes[idA1].dualCnxs, networkA.nodes[idA2].dualCnxs)[0];
            VecR<int> edgeIds = networkB.nodes[idB].dualCnxs;
            VecR<int> edges(0, edgeIds.n * 2);
            for (int k = 0; k<edgeIds.n; ++k) {
                idA3 = edgeIds[k];
                idA4 = edgeIds[(k + 1) % edgeIds.n];
                if (idA3 != idA0 && idA4 != idA0) {
                    edges.addValue(idA3);
                    edges.addValue(idA4);
                }
            }
            idA3 = networkA.nodes[idA0].netCnxs[(j + 2) % cnd];
            if(idA0<idA3){
                l0a=idA0;
                l0b=idA3;
            }
            else{
                l0a=idA3;
                l0b=idA0;
            }
            //prevent overlap of edges
            for(int k=0,l=1; k<edges.n; k+=2,l+=2){
                if(edges[k]<edges[l]){
                    l1a=edges[k];
                    l1b=edges[l];
                }
                else{
                    l1a=edges[l];
                    l1b=edges[k];
                }
                if(true==true){
//                if(l0a<l1a) {//prevent double counting
                    intersections.addValue(l0a);
                    intersections.addValue(l0b);
                    intersections.addValue(l1a);
                    intersections.addValue(l1b);
                }
            }
        }
    }
    else if(cnd==4) {
        for (int j = 0; j < cnd; ++j) {
            idA1 = networkA.nodes[idA0].netCnxs[j];
            idA2 = networkA.nodes[idA0].netCnxs[(j + 1) % cnd];
            idB = vCommonValues(networkA.nodes[idA1].dualCnxs, networkA.nodes[idA2].dualCnxs)[0];
            VecR<int> edgeIds = networkB.nodes[idB].dualCnxs;
            VecR<int> edges(0, edgeIds.n * 2);
            for (int k = 0; k<edgeIds.n; ++k) {
                idA3 = edgeIds[k];
                idA4 = edgeIds[(k + 1) % edgeIds.n];
                if (idA3 != idA0 && idA4 != idA0) {
                    edges.addValue(idA3);
                    edges.addValue(idA4);
                }
            }
            idA3 = networkA.nodes[idA0].netCnxs[(j + 2) % cnd];
            idA4 = networkA.nodes[idA0].netCnxs[(j + 3) % cnd];
            if(idA0<idA3){
                l0a=idA0;
                l0b=idA3;
            }
            else{
                l0a=idA3;
                l0b=idA0;
            }
            if(idA0<idA4){
                l1a=idA0;
                l1b=idA4;
            }
            else{
                l1a=idA4;
                l1b=idA0;
            }
            //prevent overlap of edges
            for(int k=0,l=1; k<edges.n; k+=2,l+=2){
                if(edges[k]<edges[l]){
                    l2a=edges[k];
                    l2b=edges[l];
                }
                else{
                    l2a=edges[l];
                    l2b=edges[k];
                }
                if(l0a<l2a) {
                    intersections.addValue(l0a);
                    intersections.addValue(l0b);
                    intersections.addValue(l2a);
                    intersections.addValue(l2b);
                }
                if(l1a<l2a) {
                    intersections.addValue(l1a);
                    intersections.addValue(l1b);
                    intersections.addValue(l2a);
                    intersections.addValue(l2b);
                }
            }
        }
    }
}

//Update networks with geometry optimised coordinates
void LinkedNetwork::syncCoordinates() {

    //Sync A coordinates
    if(networkA.geometryCode=="2DE"){
        for(int i=0; i<networkA.nodes.n; ++i){
            networkA.nodes[i].crd[0]=crdsA[2*i];
            networkA.nodes[i].crd[1]=crdsA[2*i+1];
        }
    }
    else{
        for(int i=0; i<networkA.nodes.n; ++i){
            networkA.nodes[i].crd[0]=crdsA[3*i];
            networkA.nodes[i].crd[1]=crdsA[3*i+1];
            networkA.nodes[i].crd[2]=crdsA[3*i+2];
        }
    }

    //Sync B coordinates
    if(networkB.geometryCode=="2DE"){
        for(int i=0; i<networkB.nodes.n; ++i){
            VecF<double> x(networkB.nodes[i].dualCnxs.n);
            VecF<double> y(networkB.nodes[i].dualCnxs.n);
            for (int j = 0; j < networkB.nodes[i].dualCnxs.n; ++j) {
                x[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[0];
                y[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[1];
            }
            VecF<double> origin(2);
            origin[0]=x[0];
            origin[1]=y[0];
            x-=origin[0];
            y-=origin[1];
            for(int j=0; j<x.n; ++j) x[j]-=networkB.pb[0]*nearbyint(x[j]*networkB.rpb[0]);
            for(int j=0; j<y.n; ++j) y[j]-=networkB.pb[1]*nearbyint(y[j]*networkB.rpb[1]);
            VecF<double> c(2);
            c[0]=origin[0]+vMean(x);
            c[1]=origin[1]+vMean(y);
            networkB.nodes[i].crd=c;
        }
    }
    else{
        for (int i = 0; i < networkB.nodes.n; ++i) {
            VecF<double> x(networkB.nodes[i].dualCnxs.n);
            VecF<double> y(networkB.nodes[i].dualCnxs.n);
            VecF<double> z(networkB.nodes[i].dualCnxs.n);
            for (int j = 0; j < networkB.nodes[i].dualCnxs.n; ++j) {
                x[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[0];
                y[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[1];
                z[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[2];
            }
            VecF<double> origin(3);
            origin[0] = x[0];
            origin[1] = y[0];
            origin[2] = z[0];
            x -= origin[0];
            y -= origin[1];
            z -= origin[2];
            for (int j = 0; j < x.n; ++j) x[j] -= networkB.pb[0] * nearbyint(x[j] * networkB.rpb[0]);
            for (int j = 0; j < y.n; ++j) y[j] -= networkB.pb[1] * nearbyint(y[j] * networkB.rpb[1]);
            for (int j = 0; j < z.n; ++j) z[j] -= networkB.pb[2] * nearbyint(z[j] * networkB.rpb[2]);
            VecF<double> c(3);
            c[0] = origin[0] + vMean(x);
            c[1] = origin[1] + vMean(y);
            c[2] = origin[2] + vMean(z);
            networkB.nodes[i].crd = c;
        }
    }
}

//Get normalised probability distribution of nodes of each size in given lattice
VecF<double> LinkedNetwork::getNodeDistribution(string lattice) {

    if(lattice=="A") return networkA.getNodeDistribution();
    else return networkB.getNodeDistribution();
}

//Get Aboav-Weaire fitting parameters
VecF<double> LinkedNetwork::getAboavWeaire(string lattice) {

    if(lattice=="A") return networkA.aboavWeaireParams();
    else return networkB.aboavWeaireParams();
}

//Get assortativity
double LinkedNetwork::getAssortativity(string lattice) {

    if(lattice=="A") return networkA.assortativity();
    else return networkB.assortativity();

}

//Get alpha estimate
double LinkedNetwork::getAboavWeaireEstimate(string lattice) {

    if(lattice=="A") return networkA.aboavWeaireEstimate();
    else return networkB.aboavWeaireEstimate();
}

//Get entropy
VecF<double> LinkedNetwork::getEntropy(string lattice) {

    if(lattice=="A") return networkA.entropy();
    else return networkB.entropy();
}

//Check linked networks for consistency
bool LinkedNetwork::checkConsistency() {

    bool checkCnx=checkCnxConsistency();
    bool checkDesc=checkDescriptorConsistency();
    bool check=checkCnx*checkDesc;

    return check;
}

//Check linked networks have mutual network and dual connections
bool LinkedNetwork::checkCnxConsistency() {

    //check number of network connections is equal to number of dual connections
    bool netDualEquality=true;
    for(int i=0; i<networkA.nodes.n; ++i){
        if(networkA.nodes[i].netCnxs.n!=networkA.nodes[i].dualCnxs.n) netDualEquality=false;
    }
    for(int i=0; i<networkB.nodes.n; ++i){
        if(networkB.nodes[i].netCnxs.n!=networkB.nodes[i].dualCnxs.n) netDualEquality=false;
    }

    //check mutual network connections
    bool mutualNetCnx=true;
    int id0, id1;
    for(int i=0; i<networkA.nodes.n; ++i){
        id0=i;
        for(int j=0; j<networkA.nodes[i].netCnxs.n; ++j){
            id1=networkA.nodes[i].netCnxs[j];
            mutualNetCnx=vContains(networkA.nodes[id1].netCnxs,id0);
        }
    }
    for(int i=0; i<networkB.nodes.n; ++i){
        id0=i;
        for(int j=0; j<networkB.nodes[i].netCnxs.n; ++j){
            id1=networkB.nodes[i].netCnxs[j];
            mutualNetCnx=vContains(networkB.nodes[id1].netCnxs,id0);
        }
    }

    //check mutual dual connections
    bool mutualDualCnx=true;
    for(int i=0; i<networkA.nodes.n; ++i){
        id0=i;
        for(int j=0; j<networkA.nodes[i].dualCnxs.n; ++j){
            id1=networkA.nodes[i].dualCnxs[j];
            mutualDualCnx=vContains(networkB.nodes[id1].dualCnxs,id0);
        }
    }
    for(int i=0; i<networkB.nodes.n; ++i){
        id0=i;
        for(int j=0; j<networkB.nodes[i].dualCnxs.n; ++j){
            id1=networkB.nodes[i].dualCnxs[j];
            mutualDualCnx=vContains(networkA.nodes[id1].dualCnxs,id0);
        }
    }

    //check network connections are neighbours by lying on same ring (some highly strained cases could give a false positive)
    bool nbNetCnx=true;
    for(int i=0; i<networkA.nodes.n; ++i){
        int nCnxs=networkA.nodes[i].netCnxs.n;
        for(int j=0; j<nCnxs; ++j){
            id0=networkA.nodes[i].netCnxs[j];
            id1=networkA.nodes[i].netCnxs[(j+1)%nCnxs];
            VecR<int> common=vCommonValues(networkA.nodes[id0].dualCnxs,networkA.nodes[id1].dualCnxs);
            if(common.n==0) nbNetCnx=false;
        }
    }
    for(int i=0; i<networkB.nodes.n; ++i){
        int nCnxs=networkB.nodes[i].netCnxs.n;
        for(int j=0; j<nCnxs; ++j){
            id0=networkB.nodes[i].netCnxs[j];
            id1=networkB.nodes[i].netCnxs[(j+1)%nCnxs];
            VecR<int> common=vCommonValues(networkB.nodes[id0].dualCnxs,networkB.nodes[id1].dualCnxs);
            if(common.n==0){
                int debug=1;
                nbNetCnx=false;
            }
        }
    }

    //check dual connections are neighbours by lying on same ring (some highly strained cases could give a false positive)
    bool nbDualCnx=true;
    for(int i=0; i<networkA.nodes.n; ++i){
        int nCnxs=networkA.nodes[i].dualCnxs.n;
        for(int j=0; j<nCnxs; ++j){
            id0=networkA.nodes[i].dualCnxs[j];
            id1=networkA.nodes[i].dualCnxs[(j+1)%nCnxs];
            VecR<int> common=vCommonValues(networkB.nodes[id0].dualCnxs,networkB.nodes[id1].dualCnxs);
            common.delValue(i);
            if(common.n==0) nbDualCnx=false;
        }
    }
    for(int i=0; i<networkB.nodes.n; ++i){
        int nCnxs=networkB.nodes[i].dualCnxs.n;
        for(int j=0; j<nCnxs; ++j){
            id0=networkB.nodes[i].dualCnxs[j];
            id1=networkB.nodes[i].dualCnxs[(j+1)%nCnxs];
            VecR<int> common=vCommonValues(networkA.nodes[id0].dualCnxs,networkA.nodes[id1].dualCnxs);
            common.delValue(i);
            if(common.n==0) nbDualCnx=false;
        }
    }

    //overall flag
    bool consistent=netDualEquality*mutualNetCnx*mutualDualCnx*nbNetCnx*nbDualCnx;

    return consistent;
}

//Check linked networks have accurate descriptors
bool LinkedNetwork::checkDescriptorConsistency() {

    VecF<int> checkNodeA(networkA.nodeDistribution);
    VecF<int> checkNodeB(networkB.nodeDistribution);
    VecF< VecF<int> > checkEdgeA(networkA.edgeDistribution);
    VecF< VecF<int> > checkEdgeB(networkB.edgeDistribution);

    //Check node distribution
    bool nodeA, nodeB;
    checkNodeA=0;
    checkNodeB=0;
    for(int i=0; i<networkA.nodes.n; ++i){
        int n=networkA.nodes[i].netCnxs.n;
        ++checkNodeA[n];
    }
    for(int i=0; i<networkB.nodes.n; ++i){
        int n=networkB.nodes[i].netCnxs.n;
        ++checkNodeB[n];
    }
    nodeA=checkNodeA==networkA.nodeDistribution;
    nodeB=checkNodeB==networkB.nodeDistribution;


    //Check edge distribution
    bool edgeA=true, edgeB=true;
    for(int i=0; i<checkEdgeA.n; ++i) checkEdgeA[i]=0;
    for(int i=0; i<checkEdgeB.n; ++i) checkEdgeB[i]=0;
    for(int i=0; i<networkA.nodes.n; ++i){
        int m=networkA.nodes[i].netCnxs.n;
        for(int j=0; j<m; ++j){
            int n=networkA.nodes[networkA.nodes[i].netCnxs[j]].netCnxs.n;
            ++checkEdgeA[m][n];
        }
    }
    for(int i=0; i<networkB.nodes.n; ++i){
        int m=networkB.nodes[i].netCnxs.n;
        for(int j=0; j<m; ++j){
            int n=networkB.nodes[networkB.nodes[i].netCnxs[j]].netCnxs.n;
            ++checkEdgeB[m][n];
        }
    }
    for(int i=0; i<checkEdgeA.n; ++i){
        if(!(checkEdgeA[i]==networkA.edgeDistribution[i])){
            edgeA=false;
            break;
        }
    };
    for(int i=0; i<checkEdgeB.n; ++i){
        if(!(checkEdgeB[i]==networkB.edgeDistribution[i])){
            edgeB=false;
            break;
        }
    };

    //Overall flag
    bool consistent=nodeA*nodeB*edgeA*edgeB;

    return consistent;
}

//Wrap coordinates of lattice A if periodic
void LinkedNetwork::wrapCoordinates() {

    if(networkA.geometryCode=="2DE"){
        HI2DP potModel(networkA.pb[0],networkA.pb[1]);
        potModel.wrap(crdsA);
    }
}

//Write xyz file format of networks
void LinkedNetwork::writeXYZ(string prefix) {
    networkA.writeXYZ(prefix+"_A","O");
    networkB.writeXYZ(prefix+"_B","N");
}

//Write networks in format that can be loaded and visualised
void LinkedNetwork::write(string prefix) {
    networkA.write(prefix+"_A");
    networkB.write(prefix+"_B");
}
