//Network contains nodes and topological information

#ifndef NL_NETWORK_H
#define NL_NETWORK_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include "node.h"

using namespace std;

class Network {

private:

    //Default lattices
    void initialiseSquareLattice(int dim, int& maxCnxs); //4 coordinate nodes, forming periodic square lattice
    void initialiseTriangularLattice(int dim, int& maxCnxs); //6 coordinate nodes, forming periodic square lattice
    void initialiseSnubSquareLattice(int dim, int& maxCnxs); //5 coordinate nodes, forming periodic snub square lattice
    void initialiseMixedTSLattice(int dim, int& maxCnxs, double mixProportion); //4&6 coordinate nodes, forming mixed triangular and square lattice
    void initialiseCubicLattice(int nNodes, int& maxCnxs); //3&4 coordinate nodes, forming cubic lattice
    void initialiseGeodesicLattice(int nNodes, int& maxCnxs); //5&6 coordinate nodes, forming geodesic lattice
    void initialiseDescriptors(int maxCnxs); //node descriptors

public:

    //Data members
    VecF<double> pb, rpb; //(reciprocal) periodic boundary
    VecR<Node> nodes; //list of nodes
    string geometryCode; //geometry of system
    VecF<int> nodeDistribution; //number of each type of node
    VecF< VecF<int> > edgeDistribution; //number of each type of edge

    //Constructors
    Network();
    Network(int nNodes, int maxCnxs);
    Network(int nNodes, string lattice, int maxCnxs, double mixProportion=0.0); //construct with default lattice
    Network(string prefix); //construct by loading from files

    //Member Functions
    Network constructDual(int maxCnxs); //make dual graph
    void generateAuxConnections(Network dualNetwork, int auxType); //generate auxilary connections
    void rescale(double scaleFactor); //rescale coordinates
    void project(string projType, double param); //project lattice onto different geometry
    void findLocalRegion(int a, int b, int extent, VecR<int>& local, VecR<int>& fixedInner, VecR<int>& fixedOuter);
    VecF<double> getNodeDistribution(); //proportion of each type of node
    VecF< VecF<double> > getEdgeDistribution(); //proportion of each type of node
    VecF<double> aboavWeaireParams(); //calculate Aboav-Weaire parameters
    double assortativity(); //calculate network assortativity
    double aboavWeaireEstimate(); //estimate aw alpha parameter
    VecF<double> entropy(); //calculate entropy of node and edge distribution
    double cluster(int nodeCnd); //get cluster statistics for given node coordination
    void write(string prefix);
    void writeXYZ(string prefix, string element);
};


#endif //NL_NETWORK_H
