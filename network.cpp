#include "network.h"

//Default constructor
Network::Network() {
    nodes=VecR<Node>(0,1);
}

//Construct with number of nodes and connectivity limits
Network::Network(int nNodes, int maxCnxs) {
    nodes=VecR<Node>(0,nNodes);
    for(int i=0; i<nNodes; ++i){
        Node node(i,maxCnxs,maxCnxs);
        nodes.addValue(node);
    }
}

//Construct with choice of default lattice
Network::Network(int nNodes, string lattice, int maxCnxs) {
    if(lattice=="square") initialiseSquareLattice(sqrt(nNodes),maxCnxs);
    else if(lattice=="triangular") initialiseTriangularLattice(sqrt(nNodes),maxCnxs);
    else if(lattice=="cubic") initialiseCubicLattice(nNodes,maxCnxs);
    else if(lattice=="geodesic") initialiseGeodesicLattice(nNodes,maxCnxs);
    initialiseDescriptors(maxCnxs);
}

//Initialise node and edge statistics
void Network::initialiseDescriptors(int maxCnxs) {

    //Set sizes of vectors and matrices
    nodeDistribution=VecF<int>(maxCnxs+1);
    edgeDistribution=VecF< VecF<int> >(maxCnxs+1);
    for(int i=0; i<maxCnxs+1; ++i) edgeDistribution[i]=VecF<int>(maxCnxs+1);

    //Count number of each node type and add to vector
    for(int i=0; i<nodes.n; ++i) ++nodeDistribution[nodes[i].netCnxs.n];

    //Double count number of each edge type and add to vector
    for(int i=0; i<nodes.n; ++i){
        for(int j=0; j<nodes[i].netCnxs.n; ++j){
            ++edgeDistribution[nodes[i].netCnxs.n][nodes[j].netCnxs.n];
        }
    }
}

//Construct by loading from files
Network::Network(string prefix) {

    //Initialise variables with aux file information
    string line;
    istringstream ss("");
    ifstream auxFile(prefix+"_aux.dat", ios::in);
    int nNodes, maxNetCnxs, maxDualCnxs;
    getline(auxFile,line);
    istringstream(line)>>nNodes;
    getline(auxFile,line);
    ss.str(line);
    ss>>maxNetCnxs;
    ss>>maxDualCnxs;
    getline(auxFile,line);
    istringstream(line)>>geometryCode;
    nodes=VecR<Node>(0,nNodes);
    if(geometryCode=="2DE"){
        pb=VecF<double>(2);
        rpb=VecF<double>(2);
    }
    else if(geometryCode=="3DE" || geometryCode=="2DS"){
        pb=VecF<double>(3);
        rpb=VecF<double>(3);
    }
    getline(auxFile,line);
    ss.str(line);
    ss>>pb[0];
    ss>>pb[1];
    if(geometryCode=="3DE" || geometryCode=="2DS") ss>>pb[2];
    getline(auxFile,line);
    ss.str(line);
    ss>>rpb[0];
    ss>>rpb[1];
    if(geometryCode=="3DE" || geometryCode=="2DS") ss>>rpb[2];
    for(int i=0; i<nodes.nMax; ++i){
        Node node(i,maxNetCnxs,maxDualCnxs);
        nodes.addValue(node);
    }
    auxFile.close();

    //Read coordinates
    ifstream crdFile(prefix+"_crds.dat", ios::in);
    if(geometryCode=="2DE"){
        VecF<double> crd(2);
        for(int i=0; i<nodes.n; ++i){
            getline(crdFile,line);
            ss.str(line);
            ss>>crd[0];
            ss>>crd[1];
            nodes[i].crd=crd;
        }
    }
    else if(geometryCode=="3DE" || geometryCode=="2DS"){
        VecF<double> crd(3);
        for(int i=0; i<nodes.n; ++i){
            getline(crdFile,line);
            ss.str(line);
            ss>>crd[0];
            ss>>crd[1];
            ss>>crd[2];
            nodes[i].crd=crd;
        }
    }
    crdFile.close();

    //Read network connections
    int cnx;
    ifstream netFile(prefix+"_net.dat", ios::in);
    for(int i=0; i<nodes.n; ++i){
        getline(netFile,line);
        ss.str(line);
        while(ss>>cnx) nodes[i].netCnxs.addValue(cnx);
    }
    netFile.close();

    //Read dual connections
    ifstream dualFile(prefix+"_dual.dat", ios::in);
    for(int i=0; i<nodes.n; ++i){
        getline(dualFile,line);
        ss.str(line);
        while(ss>>cnx) nodes[i].dualCnxs.addValue(cnx);
    }
    dualFile.close();

    //Set up descriptors
    initialiseDescriptors(maxNetCnxs);
}

//Initialise square lattice of periodic 4-coordinate nodes
void Network::initialiseSquareLattice(int dim, int& maxCnxs) {
    geometryCode="2DE"; //2D euclidean
    int dimSq=dim*dim;
    nodes=VecR<Node>(0,dimSq);

    //make 4 coordinate nodes
    if(maxCnxs<4) maxCnxs=4; //need at least 4 connections
    for(int i=0; i<nodes.nMax; ++i){
        Node node(i,maxCnxs,maxCnxs);
        nodes.addValue(node);
    }

    //assign coordinates in layers, with unit bond lengths
    pb=VecF<double>(2);
    rpb=VecF<double>(2);
    pb=dim;
    rpb=1.0/dim;
    VecF<double> c(2);
    for(int y=0; y<dim; ++y){
        c[1]=0.5+y;
        for(int x=0; x<dim; ++x){
            c[0]=0.5+x;
            nodes[x+y*dim].crd=c;
        }
    }

    //make connections to nodes in clockwise order
    int id=0, cnx;
    for(int y=0; y<dim; ++y){
        for(int x=0; x<dim; ++x){
            cnx=y*dim+(id+dim-1)%dim;
            nodes[id].netCnxs.addValue(cnx);
            cnx=(id+dim)%dimSq;
            nodes[id].netCnxs.addValue(cnx);
            cnx=y*dim+(id+1)%dim;
            nodes[id].netCnxs.addValue(cnx);
            cnx=(id+dimSq-dim)%dimSq;
            nodes[id].netCnxs.addValue(cnx);
            ++id;
        }
    }

    //make connections to dual nodes in clockwise order
    id=0;
    for(int y=0; y<dim; ++y){
        for(int x=0; x<dim; ++x){
            cnx=id;
            nodes[id].dualCnxs.addValue(cnx);
            cnx=y*dim+(id+1)%dim;
            nodes[id].dualCnxs.addValue(cnx);
            cnx=((y*dim+(id+1)%dim)+dimSq-dim)%dimSq;
            nodes[id].dualCnxs.addValue(cnx);
            cnx=(id+dimSq-dim)%dimSq;
            nodes[id].dualCnxs.addValue(cnx);
            ++id;
        }
    }
}

//Initialise triangular lattice of periodic 6-coordinate nodes
void Network::initialiseTriangularLattice(int dim, int& maxCnxs) {
    geometryCode="2DE"; //2D euclidean
    int dimSq=dim*dim;
    nodes=VecR<Node>(0,dimSq);

    //make 6 coordinate nodes
    if(maxCnxs<6) maxCnxs=6; //need at least 6 connections
    for(int i=0; i<nodes.nMax; ++i){
        Node node(i,maxCnxs,maxCnxs);
        nodes.addValue(node);
    }

    //assign coordinates in layers, with unit bond lengths
    pb=VecF<double>(2);
    rpb=VecF<double>(2);
    pb[0]=dim;
    pb[1]=dim*sqrt(3)*0.5;
    rpb=1.0/dim;
    VecF<double> c(2);
    double dy=sqrt(3.0)*0.5;
    for(int y=0; y<dim; ++y){
        c[1]=0.5*dy+y*dy;
        for(int x=0; x<dim; ++x){
            c[0]=0.5*(y%2)+x;
            nodes[x+y*dim].crd=c;
        }
    }

    //make connections to nodes in clockwise order
    int id=0, cnx;
    for(int y=0; y<dim; ++y){
        for(int x=0; x<dim; ++x){
            cnx=y*dim+(id+dim-1)%dim;
            nodes[id].netCnxs.addValue(cnx);
            if(y%2==0){
                cnx=(cnx+dim)%dimSq;
                nodes[id].netCnxs.addValue(cnx);
                cnx=(id+dim)%dimSq;
                nodes[id].netCnxs.addValue(cnx);
            }
            else{
                cnx=(id+dim)%dimSq;
                nodes[id].netCnxs.addValue(cnx);
                cnx=((y+1)*dim+(id+1)%dim)%dimSq;
                nodes[id].netCnxs.addValue(cnx);
            }
            cnx=y*dim+(id+1)%dim;
            nodes[id].netCnxs.addValue(cnx);
            if(y%2==0){
                cnx=(id+dimSq-dim)%dimSq;
                nodes[id].netCnxs.addValue(cnx);
                cnx=(dimSq+(y-1)*dim+(id+dim-1)%dim)%dimSq;
                nodes[id].netCnxs.addValue(cnx);
            }
            else{
                cnx=(dimSq+(y-1)*dim+(id+dim+1)%dim)%dimSq;
                nodes[id].netCnxs.addValue(cnx);
                cnx=(dimSq+(y-1)*dim+(id+dim)%dim)%dimSq;
                nodes[id].netCnxs.addValue(cnx);
            }
            ++id;
        }
    }

    //make connections to dual nodes in clockwise order
    id=0;
    int dimSq2=2*dimSq;
    for(int y=0; y<dim; ++y){
        for(int x=0; x<dim; ++x){
            cnx=(2*y*dim+(id+dim-1)%dim);
            nodes[id].dualCnxs.addValue(cnx);
            cnx=((2*y+1)*dim+(id)%dim);
            nodes[id].dualCnxs.addValue(cnx);
            cnx=(2*y*dim+id%dim);
            nodes[id].dualCnxs.addValue(cnx);
            if(y%2==0){
                cnx=(2*y*dim+id%dim-dim+dimSq2)%dimSq2;
                nodes[id].dualCnxs.addValue(cnx);
                cnx=(2*(y-1)*dim+(id+dim-1)%dim+dimSq2)%dimSq2;
                nodes[id].dualCnxs.addValue(cnx);
                cnx=(2*(y-1)*dim+(id+dim-1)%dim+dim+dimSq2)%dimSq2;
                nodes[id].dualCnxs.addValue(cnx);
            }
            else{
                cnx=(2*y*dim+(id+1)%dim-dim+dimSq2)%dimSq2;
                nodes[id].dualCnxs.addValue(cnx);
                cnx=(2*(y-1)*dim+id%dim+dimSq2)%dimSq2;
                nodes[id].dualCnxs.addValue(cnx);
                cnx=(2*(y-1)*dim+id%dim+dim+dimSq2)%dimSq2;
                nodes[id].dualCnxs.addValue(cnx);
            }
            ++id;
        }
    }
}

//Initialise cubic lattice of 3/4-coordinate nodes i.e. square lattice on a cube
void Network::initialiseCubicLattice(int nNodes, int& maxCnxs) {
    geometryCode = "3DE"; //3D euclidean

    /*The lattice will be constructed and numbered as follows:
     * consider cube as a net
     * long section consisting of 4 squares
     * short section consisting of 2 squares */
    nodes=VecR<Node>(0,nNodes);
    int l=(12+sqrt(144-24*(8-nNodes)))/12; //number of nodes on a given edge including vertices
    int dimL0 = (l-1)*4; //number of nodes in long section x direction
    int dimL1 = l; //number of nodes in long section y direction

    //make 4 coordinate nodes
    if(maxCnxs<4) maxCnxs=4; //need at least 4 connections
    for(int i=0; i<nodes.nMax; ++i){
        Node node(i,maxCnxs,maxCnxs);
        nodes.addValue(node);
    }

    //assign coordinates with unit bond lengths
    pb=VecF<double>(3);
    rpb=VecF<double>(3);
    pb=l*100; //not a periodic system so set box size as arbitrarily large
    rpb=1.0/l;
    VecF<double> c(3);
    VecF<double> dir(3);
    int id=0;
    for(int z=0; z<dimL1; ++z){
        for(int i=0; i<dimL0; ++i){
            if((id/(l-1))%4==0){
                dir[0]=1;
                dir[1]=0;
            }
            else if((id/(l-1))%4==1){
                dir[0]=0;
                dir[1]=1;
            }
            else if((id/(l-1))%4==2){
                dir[0]=-1;
                dir[1]=0;
            }
            else if((id/(l-1))%4==3){
                dir[0]=0;
                dir[1]=-1;
            }
            c[2]=z;
            nodes[id].crd=c;
            c+=dir;
            ++id;
        }
    }
    for(int z=0; z<2; ++z){
        c[2]=z*(l-1);
        for(int y=0; y<(l-2); ++y){
            c[1]=(y+1);
            for(int x=0; x<(l-2); ++x){
                c[0]=(x+1);
                nodes[id].crd=c;
                ++id;
            }
        }
    }

    /*make node connections by brute force
     * search all other nodes to find which are at distance 1
     * arrange in order*/
    for(int i=0; i<nNodes; ++i){
        VecR<int> nb(0,4);
        for(int j=0; j<nNodes; ++j){
            if(fabs(vNormSq(nodes[i].crd-nodes[j].crd)-1.0)<1e-2) nb.addValue(j);
        }
        if(nb.n==3){//if 3 neighbours already ordered
            for(int j=0; j<3; ++j){
                nodes[i].netCnxs.addValue(nb[j]);
            }
        }
        else if(nb.n==4){//determine order
            int pos0=-1, pos2=-1;
            for(int j=0; j<3; ++j){
                for(int k=j+1; k<4; ++k){
                    if(fabs(vNormSq(nodes[nb[j]].crd-nodes[nb[k]].crd)-4.0)<1e-2){
                        pos0=j;
                        pos2=k;
                        break;
                    };
                }
                if(pos0>=0) break;
            }
            int pos1, pos3;
            for(int j=0; j<4; ++j){
                if(pos0!=j && pos2!=j){
                    pos1=j;
                    break;
                }
            }
            for(int j=0; j<4; ++j){
                if(pos0!=j && pos2!=j && pos1!=j) pos3=j;
            }
            nodes[i].netCnxs.addValue(nb[pos0]);
            nodes[i].netCnxs.addValue(nb[pos1]);
            nodes[i].netCnxs.addValue(nb[pos2]);
            nodes[i].netCnxs.addValue(nb[pos3]);
            if(pos0<0 || pos1<0 || pos2<0 || pos3<0) throw("Error in cubic network connection generation");
        }
        else throw("Error in cubic network connection generation");
    }

    /*make dual connections by brute force
     * generate temporary dual network with coordinates
     * search nodes to find those at distance root(2)/2 */
    Network dualTemp(6*(l-1)*(l-1),1);
    //make coordinates
    c=0.5;
    c[2]=0.0;
    id=0;
    for(int i=0; i<l-1; ++i){
        c[0]=0.5;
        for(int j=0; j<l-1; ++j){
            dualTemp.nodes[id].crd=c;
            c[0]+=1.0;
            ++id;
        }
        c[1]+=1.0;
    }
    c=0.5;
    c[2]=(l-1);
    for(int i=0; i<l-1; ++i){
        c[0]=0.5;
        for(int j=0; j<l-1; ++j){
            dualTemp.nodes[id].crd=c;
            c[0]+=1.0;
            ++id;
        }
        c[1]+=1.0;
    }
    c=0.5;
    c[1]=0;
    for(int i=0; i<l-1; ++i){
        c[0]=0.5;
        for(int j=0; j<l-1; ++j){
            dualTemp.nodes[id].crd=c;
            c[0]+=1.0;
            ++id;
        }
        c[2]+=1.0;
    }
    c=0.5;
    c[1]=(l-1);
    for(int i=0; i<l-1; ++i){
        c[0]=0.5;
        for(int j=0; j<l-1; ++j){
            dualTemp.nodes[id].crd=c;
            c[0]+=1.0;
            ++id;
        }
        c[2]+=1.0;
    }
    c=0.5;
    c[0]=0;
    for(int i=0; i<l-1; ++i){
        c[1]=0.5;
        for(int j=0; j<l-1; ++j){
            dualTemp.nodes[id].crd=c;
            c[1]+=1.0;
            ++id;
        }
        c[2]+=1.0;
    }
    c=0.5;
    c[0]=(l-1);
    for(int i=0; i<l-1; ++i){
        c[1]=0.5;
        for(int j=0; j<l-1; ++j){
            dualTemp.nodes[id].crd=c;
            c[1]+=1.0;
            ++id;
        }
        c[2]+=1.0;
    }
    //find connections
    for(int i=0; i<nNodes; ++i){
        VecR<int> nb(0,4);
        for(int j=0; j<dualTemp.nodes.n; ++j){
            if(fabs(vNormSq(nodes[i].crd-dualTemp.nodes[j].crd)-0.5)<1e-2) nb.addValue(j);
        }
        if(nb.n==3){//if 3 neighbours already ordered
            for(int j=0; j<3; ++j){
                nodes[i].dualCnxs.addValue(nb[j]);
            }
        }
        else if(nb.n==4){//determine order
            int pos0=-1, pos2=-1;
            for(int j=0; j<3; ++j){
                for(int k=j+1; k<4; ++k){
                    if(fabs(vNormSq(dualTemp.nodes[nb[j]].crd-dualTemp.nodes[nb[k]].crd)-2.0)<1e-2){
                        pos0=j;
                        pos2=k;
                        break;
                    }
                    else if(fabs(vNormSq(dualTemp.nodes[nb[j]].crd-dualTemp.nodes[nb[k]].crd)-1.5)<1e-2){
                        pos0=j;
                        pos2=k;
                        break;
                    }
                }
                if(pos0>=0) break;
            }
            int pos1, pos3;
            for(int j=0; j<4; ++j){
                if(pos0!=j && pos2!=j){
                    pos1=j;
                    break;
                }
            }
            for(int j=0; j<4; ++j){
                if(pos0!=j && pos2!=j && pos1!=j) pos3=j;
            }
            nodes[i].dualCnxs.addValue(nb[pos0]);
            nodes[i].dualCnxs.addValue(nb[pos1]);
            nodes[i].dualCnxs.addValue(nb[pos2]);
            nodes[i].dualCnxs.addValue(nb[pos3]);
            if(pos0<0 || pos1<0 || pos2<0 || pos3<0) throw("Error in cubic network dual connection generation");
        }
        else throw("Error in cubic network dual connection generation");
    }

    //recentre coordinates on origin
    VecF<double> com(3);
    for(int i=0; i<nNodes; ++i) com+=nodes[i].crd;
    com/=nNodes;
    for(int i=0; i<nNodes; ++i) nodes[i].crd-=com;
}

//Initialise geodesic lattice of 5/6-coordinate nodes i.e. triangular lattice on a icosahedron
void Network::initialiseGeodesicLattice(int nNodes, int& maxCnxs) {
    geometryCode = "3DE"; //3D euclidean

    /*The lattice will be constructed and numbered as follows:
     * consider icosahedron as top (5 faces), middle (10 faces) and bottom (5 faces)
     * vertices lie on 3 orthogonal golden rectangles
     * divide each face into triangular lattices */
    nodes = VecR<Node>(0, nNodes);
    int l = (20+sqrt(400-40*(12-nNodes)))/20; //number of nodes on a given edge including vertices

    //make 4 coordinate nodes
    if (maxCnxs < 6) maxCnxs = 6; //need at least 6 connections
    for (int i = 0; i < nodes.nMax; ++i) {
        Node node(i, maxCnxs, maxCnxs);
        nodes.addValue(node);
    }

    /* assign coordinates with unit bond lengths
     * make vertex nodes
     * make edge nodes
     * make face nodes */
    pb = VecF<double>(3);
    rpb = VecF<double>(3);
    pb = l * 100; //not a periodic system so set box size as arbitrarily large
    rpb = 1.0 / l;
    //make vertex coordinates
    VecF<double> c(3);
    c[1]=0.5*(l-1);
    c[2]=0.25*(1.0+sqrt(5))*(l-1);
    VecF<double> xt(c), xmt(c), xmb(c), xb(c); //vertex coordinates of rectangle in y-z plane (top,middle(top/bottom),bottom)
    xb=-xb;
    xmt[1]=-xmt[1];
    xmb[2]=-xmb[2];
    c=vCyclicPermutation(c);
    VecF<double> ymta(c), ymtb(c), ymba(c), ymbb(c); //vertex coordinates of rectangle in x-z plane
    ymta[0]=-ymta[0];
    ymba=-ymba;
    ymbb[2]=-ymbb[2];
    c=vCyclicPermutation(c);
    VecF<double> zmta(c), zmtb(c), zmba(c), zmbb(c); //vertex coordinates of rectangle in x-y plane
    zmta[0]=-zmta[0];
    zmba=-zmba;
    zmbb[1]=-zmbb[1];
    //assign vertex coordinates to nodes
    nodes[0].crd=xt;
    nodes[1].crd=xmt;
    nodes[2].crd=xmb;
    nodes[3].crd=xb;
    nodes[4].crd=ymta;
    nodes[5].crd=ymtb;
    nodes[6].crd=ymba;
    nodes[7].crd=ymbb;
    nodes[8].crd=zmta;
    nodes[9].crd=zmtb;
    nodes[10].crd=zmba;
    nodes[11].crd=zmbb;
    int id=12;
    //group vertices into layers
    VecF< VecF<double> > layer1(5), layer2(5);
    layer1[0]=xmt;
    layer1[1]=ymta;
    layer1[2]=zmta;
    layer1[3]=zmtb;
    layer1[4]=ymtb;
    layer2[0]=zmba;
    layer2[1]=ymba;
    layer2[2]=xmb;
    layer2[3]=ymbb;
    layer2[4]=zmbb;
    //make edge coordinates
    VecF<double> pA(3), pB(3), pC(3); //three points
    VecF<double> vecAB(3), vecAC(3); //vectors between points
    for(int i=0; i<30; ++i){
        if(i<5){
            pA=xt;
            pB=layer1[i];
        }
        else if(i<10){
            pA=layer1[i%5];
            pB=layer1[(i+1)%5];
        }
        else if(i<15){
            pA=xb;
            pB=layer2[i%5];
        }
        else if(i<20){
            pA=layer2[i%5];
            pB=layer2[(i+1)%5];
        }
        else if(i<25){
            pA=layer1[i%5];
            pB=layer2[i%5];
        }
        else{
            pA=layer1[(i+1)%5];
            pB=layer2[i%5];
        }

        vecAB=(pB-pA)/(l-1);
        for(int j=1; j<l-1; ++j) {
            nodes[id].crd = pA + vecAB * j;
            ++id;
        }
    }
    //make face coordinates
    for(int i=0; i<20; ++i){
        if(i<5){
            pA=xt;
            pB=layer1[i%5];
            pC=layer1[(i+1)%5];
        }
        else if(i<10){
            pA=xb;
            pB=layer2[i%5];
            pC=layer2[(i+1)%5];
        }
        else if(i<15){
            pA=layer1[(i+1)%5];
            pB=layer2[i%5];
            pC=layer2[(i+1)%5];
        }
        else if(i<20){
            pA=layer2[i%5];
            pB=layer1[i%5];
            pC=layer1[(i+1)%5];
        }

        vecAB=(pB-pA)/(l-1);
        vecAC=(pC-pA)/(l-1);
        for(int j=1; j<l-2; ++j){
            for(int k=1; k<l-1-j; ++k){
                nodes[id].crd=pA+vecAB*j+vecAC*k;
                ++id;
            }
        }
    }

    /*make node connections by brute force
     * search all other nodes to find which are at distance 1
     * arrange in order*/
    for(int i=0; i<nNodes; ++i) {
        //find neighhbours
        VecR<int> nb(0, 6);
        for (int j = 0; j < nNodes; ++j) {
            if (fabs(vNormSq(nodes[i].crd - nodes[j].crd) - 1.0) < 1e-2) nb.addValue(j);
        }
        //order neighbours
        VecR<int> orderedNb(0, nb.n);
        orderedNb.addValue(nb[0]);
        for (int j = 1; j < nb.n; ++j) {
            int nb0 = -1;
            for (int k = 0; k < nb.n; ++k) {
                if (fabs(vNormSq(nodes[orderedNb[j - 1]].crd - nodes[nb[k]].crd) - 1.0) < 1e-2) {
                    nb0 = k;
                    break;
                }
            }
            int nb1 = -1;
            for (int k = nb0 + 1; k < nb.n; ++k) {
                if (fabs(vNormSq(nodes[orderedNb[j - 1]].crd - nodes[nb[k]].crd) - 1.0) < 1e-2) {
                    nb1 = k;
                    break;
                }
            }
            if (nb0 == -1 || nb1 == -1) throw ("Error in geodesic network connection generation");
            if (!vContains(orderedNb, nb[nb0])) orderedNb.addValue(nb[nb0]);
            else if (!vContains(orderedNb, nb[nb1])) orderedNb.addValue(nb[nb1]);
            else throw ("Error in geodesic network connection generation");
        }
        for (int j = 0; j < orderedNb.n; ++j) nodes[i].netCnxs.addValue(orderedNb[j]);
    }

    /*make dual connections
     * find unique triangles and assigning id
     * attribute ids to nodes */
    map<string,int> triIds;
    int triId=0;
    //find unique triangles by looping over ordered network connections
    for(int i=0; i<nNodes; ++i){
        VecF<int> tri(3);
        string triCode;
        int n=nodes[i].netCnxs.n;
        for(int j=0; j<n; ++j){
            tri[0]=i;
            tri[1]=nodes[i].netCnxs[j];
            tri[2]=nodes[i].netCnxs[(j+1)%n];
            tri=vSort(tri);
            triCode="#"+to_string(tri[0])+"#"+to_string(tri[1])+"#"+to_string(tri[2]);
            if(triIds.count(triCode)==0){
                triIds[triCode]=triId;
                ++triId;
            }
        }
    }
    //add dual ids depending on triangles present - will be ordered
    for(int i=0; i<nNodes; ++i){
        VecF<int> tri(3);
        string triCode;
        int n=nodes[i].netCnxs.n;
        for(int j=0; j<n; ++j){
            tri[0]=i;
            tri[1]=nodes[i].netCnxs[j];
            tri[2]=nodes[i].netCnxs[(j+1)%n];
            tri=vSort(tri);
            triCode="#"+to_string(tri[0])+"#"+to_string(tri[1])+"#"+to_string(tri[2]);
            nodes[i].dualCnxs.addValue(triIds.at(triCode));
        }
    }

}

//Construct dual network from network
Network Network::constructDual(int maxCnxs) {
    //Doesn't need to be very optimised as only used during initialisation

    //find number of unique dual nodes and initialise
    int nNodes=-1;
    for(int i=0; i<nodes.n; ++i){
        for(int j=0; j<nodes[i].dualCnxs.n; ++j){
            if(nodes[i].dualCnxs[j]>nNodes) nNodes=nodes[i].dualCnxs[j];
        }
    }
    ++nNodes;
    Network dualNetwork(nNodes,maxCnxs);

    //add unordered dual connections
    int id;
    for(int i=0; i<nodes.n; ++i){
        for(int j=0; j<nodes[i].dualCnxs.n; ++j){
            id=nodes[i].dualCnxs[j];
            dualNetwork.nodes[id].dualCnxs.addValue(i);
        }
    }
    //order
    for(int i=0; i<dualNetwork.nodes.n; ++i){
        VecR<int> unordered=dualNetwork.nodes[i].dualCnxs;
        VecR<int> ordered(0,unordered.nMax);
        ordered.addValue(unordered[0]);
        for(int j=1; j<unordered.n; ++j){
            VecR<int> common=vCommonValues(nodes[ordered[j-1]].netCnxs,unordered);
            if(!vContains(ordered,common[0])) ordered.addValue(common[0]);
            else ordered.addValue(common[1]);
        }
        dualNetwork.nodes[i].dualCnxs=ordered;
    }

    //add unordered network connections
    for(int i=0; i<dualNetwork.nodes.n; ++i){
        VecR<int> nonUnique(0,100);
        for(int j=0; j<dualNetwork.nodes[i].dualCnxs.n; ++j){
            for(int k=0; k<nodes[dualNetwork.nodes[i].dualCnxs[j]].dualCnxs.n; ++k){
                nonUnique.addValue(nodes[dualNetwork.nodes[i].dualCnxs[j]].dualCnxs[k]);
            }
        }
        VecR<int> unique=vUnique(nonUnique);
        unique.delValue(i);
        for(int j=0; j<unique.n; ++j){//need to share an edge to be neighbours
            VecR<int> common=vCommonValues(dualNetwork.nodes[i].dualCnxs,dualNetwork.nodes[unique[j]].dualCnxs);
            if(common.n==2) dualNetwork.nodes[i].netCnxs.addValue(unique[j]);
        }
    }
    //order
    for(int i=0; i<dualNetwork.nodes.n; ++i){
        VecR<int> unordered=dualNetwork.nodes[i].netCnxs;
        VecR<int> ordered(0,unordered.nMax);
        VecR<int> dualC=dualNetwork.nodes[i].dualCnxs;
        VecR<int> common0=vCommonValues(unordered,nodes[dualC[0]].dualCnxs);
        VecR<int> common1=vCommonValues(unordered,nodes[dualC[1]].dualCnxs);
        VecR<int> common2=vCommonValues(common0,common1);
        ordered.addValue(common2[0]);
        for(int j=2; j<dualC.n; ++j){
            common0=common1;
            common1=vCommonValues(unordered,nodes[dualC[j]].dualCnxs);
            common2=vCommonValues(common0,common1);
            ordered.addValue(common2[0]);
        }
        common0=common1;
        common1=vCommonValues(unordered,nodes[dualC[0]].dualCnxs);
        common2=vCommonValues(common0,common1);
        ordered.addValue(common2[0]);
        dualNetwork.nodes[i].netCnxs=ordered;
    }


    //make coordinate at centre of dual connnections
    if(geometryCode=="2DE") {
        for (int i = 0; i < dualNetwork.nodes.n; ++i) {
            VecF<double> x(dualNetwork.nodes[i].dualCnxs.n);
            VecF<double> y(dualNetwork.nodes[i].dualCnxs.n);
            for (int j = 0; j < dualNetwork.nodes[i].dualCnxs.n; ++j) {
                x[j] = nodes[dualNetwork.nodes[i].dualCnxs[j]].crd[0];
                y[j] = nodes[dualNetwork.nodes[i].dualCnxs[j]].crd[1];
            }
            VecF<double> origin(2);
            origin[0]=x[0];
            origin[1]=y[0];
            x-=origin[0];
            y-=origin[1];
            for(int j=0; j<x.n; ++j) x[j]-=pb[0]*nearbyint(x[j]*rpb[0]);
            for(int j=0; j<y.n; ++j) y[j]-=pb[1]*nearbyint(y[j]*rpb[1]);
            VecF<double> c(2);
            c[0]=origin[0]+vMean(x);
            c[1]=origin[1]+vMean(y);
            dualNetwork.nodes[i].crd=c;
        }
    }
    else if(geometryCode=="3DE") {
        for (int i = 0; i < dualNetwork.nodes.n; ++i) {
            VecF<double> x(dualNetwork.nodes[i].dualCnxs.n);
            VecF<double> y(dualNetwork.nodes[i].dualCnxs.n);
            VecF<double> z(dualNetwork.nodes[i].dualCnxs.n);
            for (int j = 0; j < dualNetwork.nodes[i].dualCnxs.n; ++j) {
                x[j] = nodes[dualNetwork.nodes[i].dualCnxs[j]].crd[0];
                y[j] = nodes[dualNetwork.nodes[i].dualCnxs[j]].crd[1];
                z[j] = nodes[dualNetwork.nodes[i].dualCnxs[j]].crd[2];
            }
            VecF<double> origin(3);
            origin[0]=x[0];
            origin[1]=y[0];
            origin[2]=z[0];
            x-=origin[0];
            y-=origin[1];
            z-=origin[2];
            for(int j=0; j<x.n; ++j) x[j]-=pb[0]*nearbyint(x[j]*rpb[0]);
            for(int j=0; j<y.n; ++j) y[j]-=pb[1]*nearbyint(y[j]*rpb[1]);
            for(int j=0; j<z.n; ++j) z[j]-=pb[2]*nearbyint(z[j]*rpb[2]);
            VecF<double> c(3);
            c[0]=origin[0]+vMean(x);
            c[1]=origin[1]+vMean(y);
            c[2]=origin[2]+vMean(z);
            dualNetwork.nodes[i].crd=c;
        }
    }

    //set remaining parameters
    dualNetwork.pb=pb;
    dualNetwork.rpb=rpb;
    dualNetwork.geometryCode=geometryCode;
    dualNetwork.initialiseDescriptors(maxCnxs);

    return dualNetwork;
}

//Rescale coordinates and lattice dimensions
void Network::rescale(double scaleFactor) {
    pb*=scaleFactor;
    rpb/=scaleFactor;
    for(int i=0; i<nodes.n; ++i) nodes[i].crd*=scaleFactor;
}

//Project lattice onto different geometry
void Network::project(string projType, double param) {
    if(projType=="sphere"){
        if(geometryCode=="3DE"){
            VecF<double> vec(3);
            for(int i=0; i<nodes.n; ++i){
                vec=nodes[i].crd;
                vec/=vNorm(vec);
                vec*=param;
                nodes[i].crd=vec;
            }
            geometryCode="2DS";
        }
    }
}

//Find local region of lattice, nodes in a given range of central nodes
void Network::findLocalRegion(int a, int b, int extent, VecR<int> &local, VecR<int> &fixedInner, VecR<int> &fixedOuter) {

    VecR<int> localNodes(0,10000), fixedInnerNodes(0,10000), fixedOuterNodes(0,10000);
    VecR<int> layer0(0,2),layer1(0,100),layer2;
    VecR<int> common;
    layer0.addValue(a);
    layer0.addValue(b);
    for(int i=0; i<layer0.n; ++i) localNodes.addValue(layer0[i]);
    for(int i=0; i<layer0.n; ++i){
        for(int j=0; j<nodes[layer0[i]].netCnxs.n; ++j){
            layer1.addValue(nodes[layer0[i]].netCnxs[j]);
        }
    }
    layer1=vUnique(layer1);
    for(int i=0; i<layer0.n; ++i) layer1.delValue(layer0[i]);
    for(int i=0; i<layer1.n; ++i) localNodes.addValue(layer1[i]);

    for(int i=0; i<extent+1; ++i){
        layer2=VecR<int>(0,1000);
        for(int j=0; j<layer1.n; ++j){
            for(int k=0; k<nodes[layer1[j]].netCnxs.n; ++k){
                layer2.addValue(nodes[layer1[j]].netCnxs[k]);
            }
        }
        if(layer2.n==0) break;
        layer2=vUnique(layer2);
        common=vCommonValues(layer2,layer1);
        for(int j=0; j<common.n; ++j) layer2.delValue(common[j]);
        common=vCommonValues(layer2,layer0);
        for(int j=0; j<common.n; ++j) layer2.delValue(common[j]);
        if(i<extent-1){
            for(int j=0; j<layer2.n; ++j) localNodes.addValue(layer2[j]);
        }
        else if(i==extent-1){
            for(int j=0; j<layer2.n; ++j) fixedInnerNodes.addValue(layer2[j]);
        }
        else if(i==extent){
            for(int j=0; j<layer2.n; ++j) fixedOuterNodes.addValue(layer2[j]);
        }
        layer0=layer1;
        layer1=layer2;
    }

    local=VecR<int>(localNodes.n);
    fixedInner=VecR<int>(fixedInnerNodes.n);
    fixedOuter=VecR<int>(fixedOuterNodes.n);
    for(int i=0; i<local.n; ++i) local[i]=localNodes[i];
    for(int i=0; i<fixedInner.n; ++i) fixedInner[i]=fixedInnerNodes[i];
    for(int i=0; i<fixedOuter.n; ++i) fixedOuter[i]=fixedOuterNodes[i];
}

//Get proportion of nodes of each size
VecF<double> Network::getNodeDistribution() {

    VecF<double> normalisedDist(nodeDistribution.n);
    for(int i=0; i<nodeDistribution.n; ++i) normalisedDist[i]=nodeDistribution[i];
    normalisedDist/=vSum(normalisedDist);
    return normalisedDist;
}

//Get proportion of edges with a node of each size at each end
VecF< VecF<double> > Network::getEdgeDistribution() {

    VecF< VecF<double> > normalisedDist(edgeDistribution.n);
    double sum=0.0;
    for(int i=0; i<edgeDistribution.n; ++i){
        normalisedDist[i]=VecF<double>(edgeDistribution[i].n);
        for(int j=0; j<edgeDistribution[i].n; ++j) normalisedDist[i][j]=edgeDistribution[i][j];
        sum+=vSum(edgeDistribution[i]);
    }
    for(int i=0; i<edgeDistribution.n; ++i) normalisedDist[i]/=sum;

    return normalisedDist;
}

//Calculate Aboav-Weaire fitting parameters
VecF<double> Network::aboavWeaireParams() {

    /* Aboav-Weaire's law: nm_n=<n>^2+mu+<n>(1-alpha)(n-<n>)
     * calculate mean node size, <n>^2
     * calculate mean node size about a node of size n
     * perform linear fit */

    //mean from node distribution
    double mean=0.0;
    for(int i=0; i<nodeDistribution.n; ++i) mean+=i*nodeDistribution[i];
    mean/=vSum(nodeDistribution);

    //find x,y only for sizes which are present
    VecR<double> x(0,edgeDistribution.n);
    VecR<double> y(0,edgeDistribution.n);
    for(int i=0; i<edgeDistribution.n; ++i){
        int num=vSum(edgeDistribution[i]);
        if(num>0){
            double mn=0.0;
            for(int j=0; j<edgeDistribution[i].n; ++j) mn+=j*edgeDistribution[i][j];
            mn/=num;
            x.addValue(mean*(i-mean));
            y.addValue(i*mn);
        }
    }

    //linear fit if more than one data point, return: alpha, mu, rsq
    VecF<double> aw(3);
    if(x.n>1){
        VecR<double> fit=vLinearRegression(x,y);
        aw[0]=1.0-fit[0];
        aw[1]=fit[1]-mean*mean;
        aw[2]=fit[2];
    }

    return aw;
}

//Calculate network assortativity through the degree correlation coefficient
double Network::assortativity() {

    /* definitions:
     * 1) e_ij degree correlation matrix, prob of finding nodes with degree i,j at end of random link
     * 2) q_k prob of finding node with degree k at end of random link, q_k=kp_k/<k> */
    VecF< VecF<double> > e=getEdgeDistribution();
    VecF<double> q=getNodeDistribution();
    double mean=0.0;
    for(int i=0; i<q.n; ++i) mean+=i*q[i];
    for(int i=0; i<q.n; ++i) q[i]=i*q[i]/mean;

    /* degree correlation coefficient:
     * r=sum_jk jk(e_jk-q_j*q_k)/sig^2
     * sig^2 = sum_k k^2*q_k - (sum_k k*q_k)^2
     * bounded between -1 (perfectly dissortative) and 1 (perfectly assortative) with 0 as neutral */
    double r=0.0;
    for(int j=0; j<e.n; ++j){
        for(int k=0; k<e.n; ++k) r+=j*k*(e[j][k]-q[j]*q[k]);
    }
    double sigSq,a=0.0,b=0.0; //dummy variables
    for(int k=0; k<q.n; ++k){
        a+=k*k*q[k];
        b+=k*q[k];
    }
    sigSq=a-b*b;
    r/=sigSq;

    return r;
}

//Estimate alpha parameter from degree correlation coefficient
double Network::aboavWeaireEstimate() {

    /* Can derive by substituting aw law into equation for r */
    VecF<double> p=getNodeDistribution();
    VecF<double> k(p.n);
    for(int i=0; i<p.n; ++i) k[i]=i;
    double n,n2,n3,nSq;
    n=vSum(k*p);
    n2=vSum(k*k*p);
    n3=vSum(k*k*k*p);
    nSq=n*n;
    double alpha;
    double r=assortativity();
    double mu=n2-nSq;
    alpha=(-r*(n*n3-n2*n2)-mu*mu)/(nSq*mu);

    return alpha;
}

//Calculate entropy of node and edge distribution
VecF<double> Network::entropy() {

    double s0=0.0, s1=0.0;
    VecF<double> p=getNodeDistribution();
    VecF< VecF<double> > e=getEdgeDistribution();

    for(int i=0; i<p.n; ++i){
        if(p[i]>0.0) s0-=p[i]*log(p[i]);
    }

    for(int i=0; i<e.n; ++i){
        for(int j=0; j<e[i].n; ++j){
            if(e[i][j]>0.0) s1-=e[i][j]*log(e[i][j]);
        }
    }

    VecF<double> entropy(2);
    entropy[0]=s0;
    entropy[1]=s1;

    return entropy;
}

//Write network in format which can be loaded
void Network::write(string prefix) {

    //auxilary information
    ofstream auxFile(prefix + "_aux.dat", ios::in | ios::trunc);
    auxFile << fixed << showpoint << setprecision(1);
    auxFile << setw(10) << left << nodes.n <<endl;
    auxFile << setw(10) << left << nodes[0].netCnxs.nMax << setw(10) << left << nodes[0].dualCnxs.nMax << endl;
    auxFile << setw(10) << left << geometryCode <<endl;
    auxFile << fixed << showpoint << setprecision(6);
    for(int i=0; i<pb.n; ++i) auxFile << setw(20) << left << pb[i];
    auxFile<<endl;
    for(int i=0; i<rpb.n; ++i) auxFile << setw(20) << left << rpb[i];
    auxFile <<endl;
    auxFile.close();

    //coordinates
    ofstream crdFile(prefix + "_crds.dat", ios::in | ios::trunc);
    crdFile << fixed << showpoint << setprecision(6);
    for(int i=0; i<nodes.n; ++i){
        for(int j=0; j<nodes[i].crd.n; ++j){
            crdFile<< setw(20) << left << nodes[i].crd[j];
        }
        crdFile<<endl;
    }
    crdFile.close();

    //network connections
    ofstream netFile(prefix + "_net.dat", ios::in | ios::trunc);
    netFile << fixed << showpoint << setprecision(1);
    for(int i=0; i<nodes.n; ++i){
        for(int j=0; j<nodes[i].netCnxs.n; ++j){
            netFile<< setw(20) << left << nodes[i].netCnxs[j];
        }
        netFile<<endl;
    }
    netFile.close();

    //dual connections
    ofstream dualFile(prefix + "_dual.dat", ios::in | ios::trunc);
    dualFile << fixed << showpoint << setprecision(1);
    for(int i=0; i<nodes.n; ++i){
        for(int j=0; j<nodes[i].dualCnxs.n; ++j){
            dualFile<< setw(20) << left << nodes[i].dualCnxs[j];
        }
        dualFile<<endl;
    }
    dualFile.close();
}

//Write xyz file format of network
void Network::writeXYZ(string prefix, string element) {
    ofstream xyzFile(prefix + ".xyz", ios::in | ios::trunc);
    if (geometryCode == "2DE") {
        xyzFile << nodes.n << endl;
        xyzFile << "# Element X Y Z  " << endl;
        xyzFile << fixed << showpoint << setprecision(6);
        for (int i = 0; i < nodes.n; ++i) {
            xyzFile << setw(10) << left << element + to_string(i);
            xyzFile << setw(20) << left << nodes[i].crd[0];
            xyzFile << setw(20) << left << nodes[i].crd[1];
            xyzFile << setw(20) << left << "0.0" << endl;
        }
    }
    else if (geometryCode == "2DS") {
        xyzFile << nodes.n << endl;
        xyzFile << "# Element X Y Z  " << endl;
        xyzFile << fixed << showpoint << setprecision(6);
        for (int i = 0; i < nodes.n; ++i) {
            xyzFile << setw(10) << left << element + to_string(i);
            xyzFile << setw(20) << left << nodes[i].crd[0];
            xyzFile << setw(20) << left << nodes[i].crd[1];
            xyzFile << setw(20) << left << nodes[i].crd[2]<<endl;
        }
    }
    else if (geometryCode == "3DE") {
        xyzFile << nodes.n << endl;
        xyzFile << "# Element X Y Z  " << endl;
        xyzFile << fixed << showpoint << setprecision(6);
        for (int i = 0; i < nodes.n; ++i) {
            xyzFile << setw(10) << left << element + to_string(i);
            xyzFile << setw(20) << left << nodes[i].crd[0];
            xyzFile << setw(20) << left << nodes[i].crd[1];
            xyzFile << setw(20) << left << nodes[i].crd[2]<<endl;
        }
    }
    xyzFile.close();
}
