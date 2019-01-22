//Linked reciprocal networks - network and dual pair

#ifndef NL_LINKED_NETWORK_H
#define NL_LINKED_NETWORK_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <random>
#include "network.h"
#include "pot2d.h"
#include "pot3d.h"
#include "opt.h"
#include "monte_carlo.h"

using namespace std;

class LinkedNetwork {
private:
public:

    //Data members
    Network networkA, networkB; //two reciprocal networks
    VecF<double> crdsA; //copy of coordinates in network A (for efficient geometry optimisation)
    mt19937 mtGen; //mersenne twister random number generator
    Metropolis mc; //monte carlo metropolis condition
    VecF<double> potParamsA,potParamsB,potParamsC; //potetial model parameters (angles, bonds, constraints)
    VecF<int> potParamsD; //potential model intersection parameters
    VecF<int> goptParamsA; //geometry optimisation parameters
    VecF<double> goptParamsB; //geometry optimisation parameters

    //Additional data members
    int minNodeCnxs;

    //Constructors
    LinkedNetwork();
    LinkedNetwork(int nodesA, string latticeA, int maxACnxs, int maxBCnxs, int minCnxs); //construct with starting A lattice
    LinkedNetwork(string prefix);  //construct by loading from files

    //Member Functions
    void initialisePotentialModel(double ak, double bk, double ck=0.0, int convexity=1); //set up potential model
    void initialiseGeometryOpt(int iterations, double tau, double tolerance, int localExtent); //set up geometry optimsiation parameters
    void initialiseMonteCarlo(double temperature, int seed=0); //set up monte carlo
    void rescale(double scaleFactor); //rescale lattice dimensions
    void project(string projType, double param); //project lattice onto different geometry
    void optimalProjection(string projType); //project lattice onto different geometry with optimal parameters
    void randomCnx(int& a, int& b, int& u, int& v, mt19937& gen); //choose nodes forming random edge in lattice A, and corresponding nodes in lattice B
    int generateSwitchIds(VecF<int>& switchIdsA, VecF<int>& switchIdsB, int a, int b, int u=-1, int v=-1); //get all ids of nodes in lattice A and B needed for switch move
    void switchCnx(VecF<int> switchIdsA, VecF<int> switchIdsB); //switch connectivities in lattice
    VecF<int> monteCarloSwitchMove(double& energy); //monte carlo switching move
    double globalPotentialEnergy(bool useIntx, bool keepConvex); //calculate potential energy of entire system
    void globalGeometryOptimisation(bool useIntx, bool keepConvex); //geometry optimise entire system
    VecF<int> localGeometryOptimisation(int centreA, int centreB, int extent, bool useIntx, bool keepConvex); //geometry optimise subsection of system
    void generateHarmonics(int id, VecR<int>& bonds, VecR<double>& bondParams); //generate harmonic interactions
    void generateIntersections(int id, VecR<int>& intersections, bool keepConvex); //generate intersection interactions
    void wrapCoordinates(); //wrap coordinates if periodic
    void syncCoordinates(); //update geometry optimised coordinates to networks
    VecF<double> getNodeDistribution(string lattice); //get proportion of nodes of each size
    VecF<double> getAboavWeaire(string lattice); //get aboav-weaire parameters
    double getAssortativity(string lattice); //get network assortativity
    double getAboavWeaireEstimate(string lattice); //get estimate of aw alpha parameter from assortativity
    VecF<double> getEntropy(string lattice); //get node and edge distribution entropy
    bool checkConsistency(); //check networks are consistent
    bool checkCnxConsistency(); //check for mutual connections
    bool checkDescriptorConsistency(); //check descriptors are accurate

    //Write Functions
    void writeXYZ(string prefix);
    void write(string prefix);
};


#endif //NL_LINKED_NETWORK_H
