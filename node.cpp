#include "node.h"

//Default constructor
Node::Node() {
    id=-1;
    crd=VecF<double>(1);
    netCnxs=VecR<int>(1);
    dualCnxs=VecR<int>(1);
    auxCnxs=VecR<int>(1);
}

//Construct with maximum number of connections, initialise with 0
Node::Node(int nodeId, int maxNetCnxs, int maxDualCnxs, int maxAuxCnxs) {
    id=nodeId;
    netCnxs=VecR<int>(0,maxNetCnxs);
    dualCnxs=VecR<int>(0,maxDualCnxs);
    auxCnxs=VecR<int>(0,maxAuxCnxs);
}

//Copy constructor
Node::Node(const Node& source) {
    id=source.id;
    crd=source.crd;
    netCnxs=source.netCnxs;
    dualCnxs=source.dualCnxs;
    auxCnxs=source.auxCnxs;
}

//Assigment operator
Node& Node::operator=(const Node& source) {
    if(this != &source){
        id=source.id;
        crd=source.crd;
        netCnxs=source.netCnxs;
        dualCnxs=source.dualCnxs;
        auxCnxs=source.auxCnxs;
    }
    return *this;
}

