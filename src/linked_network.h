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
    Metropolis mc,mcCost; //monte carlo metropolis condition
    VecF<double> potParamsA,potParamsB,potParamsC; //potential model parameters (angles, bonds, constraints)
    VecF<int> potParamsD; //potential model intersection/convex parameters
    VecF<int> goptParamsA; //geometry optimisation parameters
    VecF<double> goptParamsB; //geometry optimisation parameters
    VecF<double> costParams; //cost function parameters

    //Additional data members
    int minACnxs,maxACnxs,minBCnxs,maxBCnxs;

    //Constructors
    LinkedNetwork();
    LinkedNetwork(int nodesA, string latticeA, int minA, int maxA, int minB, int maxB); //construct with starting A lattice
    LinkedNetwork(string prefix);  //construct by loading from files
    bool checkRingNodesUnique();
    //Member Functions
    void initialisePotentialModel(double ak, double bk, double ck=0.0, int convex=0); //set up potential model
    void initialiseGeometryOpt(int iterations, double tau, double tolerance, int localExtent); //set up geometry optimsiation parameters
    void initialiseMonteCarlo(double temperature, int seed=0, bool globalOpt=false); //set up monte carlo
    void initialiseCostFunction(double temperature, int seed, double pk, double rk); //set up cost function
    void makeCrystal(string crystalCode, string lattice); //perform defined moves to make specific crystal
    void rescale(double scaleFactor); //rescale lattice dimensions
    void project(string projType, double param); //project lattice onto different geometry
    void optimalProjection(string projType); //project lattice onto different geometry with optimal parameters
    int pickRandomCnx34(int &a, int &b, int &u, int &v, mt19937 &gen); //choose nodes forming random edge in lattice A, and corresponding nodes in lattice B
    int pickRandomCnx(int& a, int& b, int& u, int& v, mt19937& gen); //choose nodes forming random edge in lattice A, and corresponding nodes in lattice B
    int generateSwitchIds34(int cnxType, VecF<int>& switchIdsA, VecF<int>& switchIdsB, int a, int b, int u, int v); //get all ids of nodes in lattice A and B needed for switch move
    int generateMixIds34(int cnxType, VecF<int>& mixIdsA, VecF<int>& mixIdsB, int a, int b, int u, int v); //get all ids of nodes in lattice A and B needed for mix move
    int generateMixIds(int cnxType, VecF<int>& mixIdsA, VecF<int>& mixIdsB, int a, int b, int u, int v); //get all ids of nodes in lattice A and B needed for mix move
    int findAssociatedNodeAB(int idA, int idB, int idDel); //
    int findAssociatedNodeAA(int idA, int idB, int idDel); //
    void switchCnx33(VecF<int> switchIdsA, VecF<int> switchIdsB); //switch connectivities in lattice between 2x3 coordinate nodes
    void switchCnx44(VecF<int> switchIdsA, VecF<int> switchIdsB); //switch connectivities in lattice between 2x4 coordinate nodes
    void switchCnx43(VecF<int> switchIdsA, VecF<int> switchIdsB); //switch connectivities in lattice between 4 and 3 coordinate nodes
    void mixCnx34(VecF<int> mixIdsA, VecF<int> mixIdsB); //mix connectivities in lattice between 4 and 3 coordinate nodes
    void mixCnx(VecF<int> mixIdsA, VecF<int> mixIdsB); //mix connectivities in lattice between 4 and 3 coordinate nodes
    bool mixCheckEdges(int id); //prevent edges being part of three rings
    bool convexRearrangement(int cnxType, VecF<int> switchIdsA, VecF<int> switchIdsB); //rearrange nodes after switchi to maintain convexity
    VecF<int> monteCarloSwitchMove(double& energy); //monte carlo switching move
    VecF<int> monteCarloCostSwitchMove(double& cost, double& energy, double pTarget, double rTarget); //monte carlo switching move with cost function
    VecF<int> monteCarloMixMove(double& energy); //monte carlo mixing move
    double costFunction(double& pTarget, double& rTarget); //cost function based on ring statistics and assortative mixing
    double globalPotentialEnergy(bool useIntx, bool restrict); //calculate potential energy of entire system
    void globalGeometryOptimisation(bool useIntx, bool restrict); //geometry optimise entire system
    VecF<int> localGeometryOptimisation(int centreA, int centreB, int extent, bool useIntx, bool restrict); //geometry optimise subsection of system
    void generateHarmonics(int id, VecR<int>& bonds, VecR<double>& bondParams, VecR<int>& angles, VecR<double>& angleParams); //generate harmonic interactions
    void generateRingIntersections(int rId, VecR<int>& intersections); //generate ring intersection interactions
    void generateConvexIntersections(int nId, VecR<int>& intersections); //generate convex ring intersection interactions
    void wrapCoordinates(); //wrap coordinates if periodic
    void syncCoordinates(); //update geometry optimised coordinates to networks
    VecF<double> getNodeDistribution(string lattice); //get proportion of nodes of each size
    VecF< VecF<int> > getEdgeDistribution(string lattice); //get unnormalised proportion of node connections
    VecF<double> getAboavWeaire(string lattice); //get aboav-weaire parameters
    double getAssortativity(string lattice); //get network assortativity
    double getAboavWeaireEstimate(string lattice); //get estimate of aw alpha parameter from assortativity
    VecF<double> getEntropy(string lattice); //get node and edge distribution entropy
    VecF<double> getOptimisationGeometry(VecF<double> &lenHist, VecF<double> &angHist); //get bond/angle mean and standard deviation
    void getRingAreas(VecF<double> &areaSum, VecF<double> &areaSqSum); //get sum of areas and squared areas of each ring size
    double getMaxCluster(string lattice, int nodeCnd); //get cluster statistics for given node coordination
    VecF<int> getMaxClusters(string lattice, int minCnd, int maxCnd); //get cluster statistics for node coordinations
    bool checkConsistency(); //check networks are consistent
    bool checkCnxConsistency(); //check for mutual connections
    bool checkDescriptorConsistency(); //check descriptors are accurate
    bool checkConvexity(); //check all angles are convex
    bool checkConvexity(int id); //check angles are convex around given node


    //Write Functions
    void writeXYZ(string prefix);
    void write(string prefix);
};


#endif //NL_LINKED_NETWORK_H
