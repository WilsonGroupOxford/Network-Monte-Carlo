//Node in network containing coordinate, connections to nodes in network and connections to nodes in dual

#ifndef NL_NODE_H
#define NL_NODE_H

#include <iostream>
#include "vecf.h"
#include "vecr.h"
#include "vec_func.h"

using namespace std;

//Node in network, aka vertex in graph
class Node{

public:

    //Data members
    int id; //unique id used by other nodes for connections
    VecF<double> crd; //coordinate
    VecR<int> netCnxs; //connections to nodes in network
    VecR<int> dualCnxs; //connections to nodes in dual

    //Constructors, copy constructor, assignment operator
    Node();
    Node(int nodeId, int maxNetCnxs, int maxDualCnxs);
    Node(const Node& source);
    Node& operator = (const Node& source);
};

#endif //NL_NODE_H
