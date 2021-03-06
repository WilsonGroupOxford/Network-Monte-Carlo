#include <iostream>
#include <sstream>
#include "outputfile.h"
#include "linked_network.h"

using namespace std;

double rescale_box(LinkedNetwork& link_net, const double scale_x, const double scale_y) {
    //! Rescale the box with in two directions.
    //!
    /*!
     * \param link_net a linked network to rescale
     * \param scale_x the scale factor in the x direction
     * \param scale_y the scale factor in the y direction
     *
     */
    // Rescale the periodic boundaries for networks A and B
    link_net.networkA.pb[0] *= scale_x;
    link_net.networkA.pb[1] *= scale_y;
    link_net.networkA.rpb[0] = 1.0 / link_net.networkA.pb[0];
    link_net.networkA.rpb[1] = 1.0 / link_net.networkA.pb[1];
    
    link_net.networkB.pb[0] *= scale_x;
    link_net.networkB.pb[1] *= scale_y;
    link_net.networkB.rpb[0] = 1.0 / link_net.networkB.pb[0];
    link_net.networkB.rpb[1] = 1.0 / link_net.networkB.pb[1];
    
    // Change crdsA and then force them into the connected networks.
    for(int i=0; i<link_net.networkA.nodes.n; ++i){
        link_net.crdsA[2*i] *= scale_x;
        link_net.crdsA[2*i+1] *= scale_y;       
    }
    link_net.syncCoordinates();
       
    // Force the network to re-calculate its energy so we don't drop out of sync
    link_net.globalGeometryOptimisation(false, false);
    auto new_energy = link_net.globalPotentialEnergy(false, false);
    link_net.mc.setEnergy(new_energy);
    return new_energy;   
}

int main(){

    //Set up logfile
    Logfile logfile("./netmc.log");
    logfile.datetime("Simulation begun at: ");
    logfile.write("NETwork Monte Carlo");
    logfile.write("Written By: David OM, Wilson Group, 2019");
    logfile.separator();

    //Read input file
    logfile.write("Reading input file");
    ++logfile.currIndent;
    ifstream inputFile("./netmc.inpt", ios::in);
    if(!inputFile.good()) logfile.criticalError("Cannot find input file 'netmc.inpt' in current directory");
    //I/O
    string skip,line,prefixOut; //prefix for output files
    getline(inputFile,skip);
    getline(inputFile,line);
    istringstream(line)>>prefixOut;
    getline(inputFile,skip);
    getline(inputFile,skip);
    logfile.write("I/O read");
    //Network properties
    int nRings,minRingSize,maxRingSize;
    string lattice;
    string crystal;
    string moveType;
    int minCnd,maxCnd;
    getline(inputFile,line);
    istringstream(line)>>nRings;
    getline(inputFile,line);
    istringstream(line)>>minRingSize;
    getline(inputFile,line);
    istringstream(line)>>maxRingSize;
    getline(inputFile,line);
    istringstream(line)>>minCnd;
    getline(inputFile,line);
    istringstream(line)>>maxCnd;
    getline(inputFile,line);
    istringstream(line)>>lattice;
    getline(inputFile,line);
    istringstream(line)>>crystal;
    getline(inputFile,skip);
    getline(inputFile,skip);
    logfile.write("Network properties read");
    //Monte carlo
    string runType;
    int randomSeed;
    getline(inputFile,line);
    istringstream(line)>>runType;
    getline(inputFile,line);
    istringstream(line)>>moveType;
    getline(inputFile,line);
    istringstream(line)>>randomSeed;
    getline(inputFile,skip);
    getline(inputFile,skip);
    logfile.write("Monte Carlo parameters read");
    //Energy search
    int mcSteps,equilSteps;
    double mcStartT,mcEndT,mcIncT,mcThermT;
    getline(inputFile,line);
    istringstream(line)>>mcStartT;
    getline(inputFile,line);
    istringstream(line)>>mcEndT;
    getline(inputFile,line);
    istringstream(line)>>mcIncT;
    getline(inputFile,line);
    istringstream(line)>>mcThermT;
    getline(inputFile,line);
    istringstream(line)>>mcSteps;
    getline(inputFile,line);
    istringstream(line)>>equilSteps;
    getline(inputFile,skip);
    getline(inputFile,skip);
    logfile.write("Energy search parameters read");
    //Cost search
    int costEquilSteps,costInitSteps,costSteps;
    double costLowLimP,costUpLimP,costIncP,costLowLimR,costUpLimR,costIncR,costT;
    getline(inputFile,line);
    istringstream(line)>>costEquilSteps;
    getline(inputFile,line);
    istringstream(line)>>costLowLimP;
    getline(inputFile,line);
    istringstream(line)>>costUpLimP;
    getline(inputFile,line);
    istringstream(line)>>costIncP;
    getline(inputFile,line);
    istringstream(line)>>costLowLimR;
    getline(inputFile,line);
    istringstream(line)>>costUpLimR;
    getline(inputFile,line);
    istringstream(line)>>costIncR;
    getline(inputFile,line);
    istringstream(line)>>costInitSteps;
    getline(inputFile,line);
    istringstream(line)>>costSteps;
    getline(inputFile,line);
    istringstream(line)>>costT;
    getline(inputFile,skip);
    getline(inputFile,skip);
    logfile.write("Cost search parameters read");
    //Cost function
    double costPK,costRK;
    getline(inputFile,line);
    istringstream(line)>>costPK;
    getline(inputFile,line);
    istringstream(line)>>costRK;
    getline(inputFile,skip);
    getline(inputFile,skip);
    logfile.write("Cost function parameters read");
    //Potential model
    int potConvex;
    double potAK,potBK,potCK;
    getline(inputFile,line);
    istringstream(line)>>potBK;
    getline(inputFile,line);
    istringstream(line)>>potAK;
    getline(inputFile,line);
    istringstream(line)>>potCK;
    getline(inputFile,line);
    istringstream(line)>>potConvex;
    getline(inputFile,skip);
    getline(inputFile,skip);
    logfile.write("Potential parameters read");
    //Geometry optimisation
    int goptLocalIt, goptGlobalIt, goptLocalSize;
    double goptTau, goptTol;
    getline(inputFile,line);
    istringstream(line)>>goptLocalIt;
    getline(inputFile,line);
    istringstream(line)>>goptGlobalIt;
    getline(inputFile,line);
    istringstream(line)>>goptTau;
    getline(inputFile,line);
    istringstream(line)>>goptTol;
    getline(inputFile,line);
    istringstream(line)>>goptLocalSize;
    getline(inputFile,skip);
    getline(inputFile,skip);
    logfile.write("Geometry optimisation parameters read");
    //Analysis
    int analysisFreq,writeStructures,structureFreq;
    getline(inputFile,line);
    istringstream(line)>>analysisFreq;
    getline(inputFile,line);
    istringstream(line)>>writeStructures;
    getline(inputFile,line);
    istringstream(line)>>structureFreq;
    getline(inputFile,line);
    logfile.write("Write parameters read");
    getline(inputFile, line);
    double scale_factor = 1.0;
    std::istringstream(line) >> scale_factor;
    double stretchFacX = 1.0;
    double stretchFacY = 1.0;    
    getline(inputFile, line);   
    std::istringstream(line) >> stretchFacX >> stretchFacY; 
    

    int stretchFreq = 1;
    getline(inputFile, line);
    std::istringstream(line) >> stretchFreq;
    
    inputFile.close();
    --logfile.currIndent;
    logfile.write("Input file closed");
    logfile.separator();

    //Initialise network
    logfile.write("Initialising network");
    ++logfile.currIndent;
    logfile.write("Lattice type:", lattice);
    logfile.write("Number of rings:",nRings);
    logfile.write("Min ring size:",minRingSize);
    logfile.write("Max ring size:",maxRingSize);
    logfile.write("Min node coordination:",minCnd);
    logfile.write("Max node coordination:",maxCnd);
    logfile.write("Monte Carlo run type:",runType);
    logfile.write("Monte Carlo move type:",moveType);
    if(runType=="energy") {
        logfile.write("Monte carlo initial log temperature:", mcStartT);
        logfile.write("Monte carlo final log temperature:", mcEndT);
        logfile.write("Monte carlo log temperature increment:", mcIncT);
        logfile.write("Monte carlo steps per increment:", mcSteps);
        logfile.write("Equilibrium steps:", equilSteps);
        logfile.write("Stretching Frequency:", stretchFreq);
        logfile.write("Stretching Factor X:", stretchFacX);
        logfile.write("Stretching Factor Y:", stretchFacY);
        logfile.write("Total Stretch X:", std::pow(stretchFacX, mcSteps/stretchFreq));
        logfile.write("Total Stretch Y:", std::pow(stretchFacY, mcSteps/stretchFreq));
    }
    else if(runType=="cost"){
        logfile.write("Monte carlo p lower limit:", costLowLimP);
        logfile.write("Monte carlo p upper limit:", costUpLimP);
        logfile.write("Monte carlo r lower limit:", costLowLimR);
        logfile.write("Monte carlo r upper limit:", costLowLimR);
        logfile.write("Monte carlo effective temperature:", costT);
        logfile.write("Monte carlo steps per increment:", costSteps);
    }
    bool mixedLattice=false; //whether mixed coordination lattice
    if(lattice.substr(0,3)=="mix" || lattice.substr(0,5)=="cairo" || lattice=="alt_square" || moveType=="mix") mixedLattice=true;
    LinkedNetwork network(nRings,lattice,minCnd,maxCnd,minRingSize,maxRingSize);
    network.write("./output_files/pre");
    network.initialisePotentialModel(potAK,potBK,potCK,potConvex);
    network.initialiseGeometryOpt(goptLocalIt,goptTau,goptTol,goptLocalSize);
    network.initialiseMonteCarlo(pow(10,mcStartT),randomSeed,mixedLattice);
    network.initialiseCostFunction(costT,randomSeed,costPK,costRK);
    if(crystal!="default") network.makeCrystal(crystal,lattice);
    if(lattice=="goldberg" || lattice=="inv_cubic") network.optimalProjection("sphere");
    if(mixedLattice){//get statistics on number of 3/4 coordinate nodes for mixed lattice
        double mixA=network.networkA.nodes.n,mixA3=0,mixA4=0;
        for(int i=0; i<mixA; ++i){
            if(network.networkA.nodes[i].netCnxs.n==3) ++mixA3;
            else if(network.networkA.nodes[i].netCnxs.n==4) ++mixA4;
        }
        VecF<double> ringStats = network.getNodeDistribution("B");
        double mixAv = 0.0;
        for(int i=0; i<=maxRingSize; ++i) mixAv += i*ringStats[i];
        logfile.write("Mixed lattice total rings:",network.networkB.nodes.n);
        logfile.write("Mixed average ring size:",mixAv);
        logfile.write("Mixed lattice total nodes:",mixA);
        logfile.write("Mixed lattice 3 coordinate nodes:",mixA3/mixA);
        logfile.write("Mixed lattice 4 coordinate nodes:",mixA4/mixA);
    }
    --logfile.currIndent;
    logfile.write("Network initialised");
    logfile.separator();

    //Initialise output files
    logfile.write("Initialising analysis output files");
    ++logfile.currIndent;
    OutputFile outRingStats(prefixOut+"_ringstats.out");
    OutputFile outCorr(prefixOut+"_correlations.out");
    OutputFile outEnergy(prefixOut+"_energy.out");
    OutputFile outEntropy(prefixOut+"_entropy.out");
    OutputFile outTemperature(prefixOut+"_temperature.out");
    OutputFile outGeometry(prefixOut+"_geometry.out");
    OutputFile outEmatrix(prefixOut+"_ematrix.out");
    OutputFile outGeomHist(prefixOut+"_geomhist.out");
    OutputFile outAreas(prefixOut+"_areas.out");
    OutputFile outClusterA(prefixOut+"_cluster_a.out");
    OutputFile outClusterB(prefixOut+"_cluster_b.out");
    OutputFile outCndStats(prefixOut+"_cndstats.out");
    OutputFile outPb(prefixOut+"_pb.out");
    outGeometry.initVariables(6,4,60,20);
    outAreas.initVariables(6,4,60,30);
    outEmatrix.initVariables(1,4,60,int(log10(nRings*12))+2);
    outGeomHist.initVariables(6,4,60,20);
    outClusterB.initVariables(1,4,60,10);
    logfile.write("Ring statistics file created");
    logfile.write("Correlations file created");
    logfile.write("Energy file created");
    logfile.write("Entropy file created");
    logfile.write("Temperature file created");
    logfile.write("Geometry file created");
    logfile.write("Geometry histogram file created");
    logfile.write("Edge distribution file created");
    logfile.write("Coordination statistics file created");
    --logfile.currIndent;
    logfile.write("Files initialised");
    logfile.separator();

    //Initialise total analysis variables - only update in main simulation not equilibration
    VecF<double> lenHist(10000),angHist(10000);
    lenHist=0.0;
    angHist=0.0;

    //Run monte carlo
    logfile.write("Box size before rescaling: " + std::to_string(network.networkA.pb[0]) + " x  "+ std::to_string(network.networkA.pb[1]));
    network.rescale(scale_factor);
    logfile.write("Box size after rescaling: " + std::to_string(network.networkA.pb[0]) + " x  "+ std::to_string(network.networkA.pb[1]));
    std::cout<< "Initial energy is " << network.mc.getEnergy()<< "\n";
    int accepted=0,optIterations=0;
    VecF<int> optCodes(5);
    optCodes=0;
    int trackFreq=100;
    VecF<int> moveStatus;
    if(runType=="energy") {//energy run
        //Run monte carlo thermalisation
        logfile.write("Running Monte Carlo thermalisation");
        ++logfile.currIndent;
        double energy=network.mc.getEnergy();
        double mcT = pow(10, mcThermT);
        network.mc.setTemperature(mcT);
        for (int i = 1; i <= equilSteps; ++i) {
            if(!mixedLattice) moveStatus = network.monteCarloSwitchMove(energy);
            else moveStatus = network.monteCarloMixMove(energy);
            accepted += moveStatus[0];
            optCodes[moveStatus[1]] += 1;
            optIterations += moveStatus[2];
//        cout << i << endl;
            if (i % trackFreq == 0) {
                double dt = logfile.timeElapsed();
                std::string track =
                        std::to_string(accepted) + "/" + std::to_string(i) + " moves accepted/completed in " + to_string(dt) +
                        " seconds";
                logfile.write(track);
                std::cout << track << std::endl;
            }
            if (i % analysisFreq == 0) {
                VecF<double> ringStats = network.getNodeDistribution("B");
                double r = network.getAssortativity("B");
                double aEst = network.getAboavWeaireEstimate("B");
                VecF<double> aw = network.getAboavWeaire("B");
                VecF<double> s = network.getEntropy("B");
                VecF<double> corr(6);
                VecF<double> a(maxRingSize+1),aSq(maxRingSize+1);
                double rr = network.getAssortativity("A");
                corr[0] = r;
                corr[1] = aEst;
                corr[2] = aw[0];
                corr[3] = aw[1];
                corr[4] = aw[2];
                corr[5] = rr;
                VecF<double> emptyL,emptyA; //dummy histograms
                VecF<double> geomStats = network.getOptimisationGeometry(emptyL,emptyA);
                VecF< VecF<int> > edgeDist = network.getEdgeDistribution("B");
//                VecF<int> clusters = network.getMaxClusters("B",minRingSize,maxRingSize);
                VecF<double> cndStats = network.getNodeDistribution("A");
                network.getRingAreas(a,aSq);
                outRingStats.writeRowVector(ringStats);
                outCorr.writeRowVector(corr);
                outEnergy.write(energy);
                outEntropy.writeRowVector(s);
                outTemperature.write(mcT);
                outGeometry.writeRowVector(geomStats);
                outAreas.writeRowVector(a);
                outAreas.writeRowVector(aSq);
                for(int j=0; j<edgeDist.n; ++j) outEmatrix.writeRowVector(edgeDist[j]);
//                outClusterB.writeRowVector(clusters);
                outCndStats.writeRowVector(cndStats);
                if(mixedLattice) {
                    VecF<double> cluster(3);
                    cluster[0] = network.getMaxCluster("A", 3);
                    cluster[1] = network.getMaxCluster("A", 4);
                    cluster[2] = network.getAssortativity("A");
                    outClusterA.writeRowVector(cluster);
                }
            }
        }
        --logfile.currIndent;
        logfile.write("Monte Carlo equilibration complete");
        logfile.separator();

        //Perform monte carlo simulation
        logfile.write("Running Monte Carlo simulation");
        ++logfile.currIndent;
        
        int nT = std::max(static_cast<int>((mcEndT - mcStartT) / mcIncT), 0);
        int totalsteps = 0;
        for (int t = 0; t <= nT; ++t) {
            mcT = pow(10, mcStartT + t * mcIncT);
            network.mc.setTemperature(mcT);
            logfile.write("Temperature:", mcT);
            ++logfile.currIndent;
            for (int i = 1; i <= mcSteps; ++i) {
                totalsteps++;
                if (totalsteps % stretchFreq == 0) {
                    energy = rescale_box(network, stretchFacX, stretchFacY);
                }
                
                if(!mixedLattice) {
                    moveStatus = network.monteCarloSwitchMove(energy);
                } else {
                    moveStatus = network.monteCarloMixMove(energy);
                }
                accepted += moveStatus[0];
                optCodes[moveStatus[1]] += 1;
                optIterations += moveStatus[2];
//            cout << i << endl;
                if (i % trackFreq == 0) {
                    double dt = logfile.timeElapsed();
                    string track =
                            to_string(accepted) + "/" + to_string(totalsteps) + " moves accepted/completed in " + to_string(dt) +
                            " seconds. Energy is " + std::to_string(energy);
                    logfile.write(track);
                    std::cout << track << std::endl;
                }
                if (i % analysisFreq == 0) {
                    VecF<double> ringStats = network.getNodeDistribution("B");
                    double r = network.getAssortativity("B");
                    double aEst = network.getAboavWeaireEstimate("B");
                    VecF<double> aw = network.getAboavWeaire("B");
                    VecF<double> s = network.getEntropy("B");
                    VecF<double> corr(6);
                    VecF<double> a(maxRingSize+1),aSq(maxRingSize+1);
                    double rr = network.getAssortativity("A");
                    corr[0] = r;
                    corr[1] = aEst;
                    corr[2] = aw[0];
                    corr[3] = aw[1];
                    corr[4] = aw[2];
                    corr[5] = rr;
                    VecF<double> geomStats = network.getOptimisationGeometry(lenHist,angHist);
                    VecF< VecF<int> > edgeDist = network.getEdgeDistribution("B");
//                    VecF<int> clusters = network.getMaxClusters("B",minRingSize,maxRingSize);
                    VecF<double> cndStats = network.getNodeDistribution("A");
                    network.getRingAreas(a,aSq);
                    outRingStats.writeRowVector(ringStats);
                    outCorr.writeRowVector(corr);
                    outEnergy.write(energy);
                    outEntropy.writeRowVector(s);
                    outTemperature.write(mcT);
                    outPb.write(network.networkA.pb[0], network.networkA.pb[1]);
                    outGeometry.writeRowVector(geomStats);
                    outAreas.writeRowVector(a);
                    outAreas.writeRowVector(aSq);
//                    outClusterB.writeRowVector(clusters);
                    for(int j=0; j<edgeDist.n; ++j) outEmatrix.writeRowVector(edgeDist[j]);
                    outCndStats.writeRowVector(cndStats);
                    if(mixedLattice) {
                        VecF<double> cluster(3);
                        cluster[0] = network.getMaxCluster("A", 3);
                        cluster[1] = network.getMaxCluster("A", 4);
                        cluster[2] = network.getAssortativity("A");
                        outClusterA.writeRowVector(cluster);
                    }
                }
//                cout<<i<<" "<<network.checkConsistency()<<endl;
                if (i%structureFreq==0 && writeStructures==1) {
                    network.syncCoordinates();
                    network.write(prefixOut+"_t"+to_string(t)+"_"+to_string(i));
                }
            }
            --logfile.currIndent;
        }
        --logfile.currIndent;
        logfile.write("Monte Carlo simulation complete");
        logfile.separator();
    }
    else if(runType=="cost"){

        //randomise lattice
        double energy;
        if(costEquilSteps>0){
            network.mc.setTemperature(1e10);
            logfile.write("Running Monte Carlo randomisation");
            ++logfile.currIndent;
            for (int i = 1; i <=costEquilSteps; ++i) {
                moveStatus = network.monteCarloSwitchMove(energy);
                if (i % trackFreq == 0) {
                    double dt = logfile.timeElapsed();
                    string track =
                            to_string(i) + " randomisation moves completed in " + to_string(dt) +" seconds";
                    logfile.write(track);
                    cout << "r" << " " << i << endl;
                }
            }
            --logfile.currIndent;
            logfile.write("Monte Carlo randomisation complete");
            logfile.separator();
        }

        //equilibrate
        logfile.write("Running Monte Carlo equilibration");
        ++logfile.currIndent;
        double cost;
        double costP=costLowLimP;
        double costR=costLowLimR;
        for(int i=1; i<=costInitSteps; ++i){
            moveStatus = network.monteCarloCostSwitchMove(cost,energy,costP,costR);
            accepted += moveStatus[0];
            optCodes[moveStatus[1]] += 1;
            optIterations += moveStatus[2];
//            cout << i << endl;
            if (i % trackFreq == 0) {
                double dt = logfile.timeElapsed();
                string track =
                        to_string(accepted) + "/" + to_string(i) + " moves accepted/completed in " + to_string(dt) +
                        " seconds";
                logfile.write(track);
                cout << "e" << " " << i << " " << accepted << endl;
            }
            if (i % analysisFreq == 0) {
                VecF<double> ringStats = network.getNodeDistribution("B");
                double r = network.getAssortativity("B");
                double aEst = network.getAboavWeaireEstimate("B");
                VecF<double> aw = network.getAboavWeaire("B");
                VecF<double> s = network.getEntropy("B");
                VecF<double> corr(5);
                VecF<double> a(maxRingSize+1),aSq(maxRingSize+1);
                corr[0] = r;
                corr[1] = aEst;
                corr[2] = aw[0];
                corr[3] = aw[1];
                corr[4] = aw[2];
                VecF<double> emptyL,emptyA; //dummy histograms
                VecF<double> geomStats = network.getOptimisationGeometry(emptyL,emptyA);
                VecF< VecF<int> > edgeDist = network.getEdgeDistribution("B");
                VecF<int> clusters = network.getMaxClusters("B",minRingSize,maxRingSize);
                network.getRingAreas(a,aSq);
                outRingStats.writeRowVector(ringStats);
                outCorr.writeRowVector(corr);
                outEnergy.write(energy);
                outEntropy.writeRowVector(s);
                outTemperature.write(costT);
                outGeometry.writeRowVector(geomStats);
                outAreas.writeRowVector(a);
                outAreas.writeRowVector(aSq);
                for(int j=0; j<edgeDist.n; ++j) outEmatrix.writeRowVector(edgeDist[j]);
                outClusterB.writeRowVector(clusters);
                if(mixedLattice) {
                    VecF<double> cluster(3);
                    cluster[0] = network.getMaxCluster("A", 3);
                    cluster[1] = network.getMaxCluster("A", 4);
                    cluster[2] = network.getAssortativity("A");
                    outClusterA.writeRowVector(cluster);
                }
            }
        }
        --logfile.currIndent;
        logfile.write("Monte Carlo equilibration complete");
        logfile.separator();

        //search
        logfile.write("Running Monte Carlo simulation");
        ++logfile.currIndent;
        int nR = ceil((costUpLimR - costLowLimR) / costIncR);
        int nP = ceil((costUpLimP - costLowLimP) / costIncP);
        for(int r=0; r<=nR; ++r){
            costR=costLowLimR+costIncR*r;
            logfile.write("Searching r =",costR);
            ++logfile.currIndent;
            for (int p = 0; p <= nP; ++p) {
                if(r%2==0) costP=costLowLimP+costIncP*p;
                else costP=costUpLimP-costIncP*p;
                network.mcCost.setEnergy(numeric_limits<double>::infinity());
                logfile.write("Searching p =",costP);
                ++logfile.currIndent;
                for(int i=1; i<=costSteps; ++i){
                    moveStatus = network.monteCarloCostSwitchMove(cost,energy,costP,costR);
                    accepted += moveStatus[0];
                    optCodes[moveStatus[1]] += 1;
                    optIterations += moveStatus[2];
//                    cout << i << endl;
                    if (i % trackFreq == 0) {
                        double dt = logfile.timeElapsed();
                        string track =
                                to_string(accepted) + "/" + to_string(i) + " moves accepted/completed in " + to_string(dt) +
                                " seconds";
                        logfile.write(track);
                        cout << costR << " " << costP << " " << i << " " << accepted << endl;
                    }
                    if (i % analysisFreq == 0) {
                        VecF<double> ringStats = network.getNodeDistribution("B");
                        double r = network.getAssortativity("B");
                        double aEst = network.getAboavWeaireEstimate("B");
                        VecF<double> aw = network.getAboavWeaire("B");
                        VecF<double> s = network.getEntropy("B");
                        VecF<double> corr(5);
                        VecF<double> a(maxRingSize+1),aSq(maxRingSize+1);
                        corr[0] = r;
                        corr[1] = aEst;
                        corr[2] = aw[0];
                        corr[3] = aw[1];
                        corr[4] = aw[2];
                        VecF<double> geomStats = network.getOptimisationGeometry(lenHist,angHist);
                        VecF< VecF<int> > edgeDist = network.getEdgeDistribution("B");
                        VecF<int> clusters = network.getMaxClusters("B",minRingSize,maxRingSize);
                        network.getRingAreas(a,aSq);
                        outRingStats.writeRowVector(ringStats);
                        outCorr.writeRowVector(corr);
                        outEnergy.write(energy);
                        outPb.write(network.networkA.pb[0], network.networkA.pb[1]);
                        outEntropy.writeRowVector(s);
                        outTemperature.write(costT);
                        outGeometry.writeRowVector(geomStats);
                        outAreas.writeRowVector(a);
                        outAreas.writeRowVector(aSq);
                        for(int j=0; j<edgeDist.n; ++j) outEmatrix.writeRowVector(edgeDist[j]);
                        outClusterB.writeRowVector(clusters);
                        if(mixedLattice) {
                            VecF<double> cluster(3);
                            cluster[0] = network.getMaxCluster("A", 3);
                            cluster[1] = network.getMaxCluster("A", 4);
                            cluster[2] = network.getAssortativity("A");
                            outClusterA.writeRowVector(cluster);
                        }
                    }
                    if (i%structureFreq==0 && writeStructures==0) {
                        network.syncCoordinates();
                        network.write(prefixOut+"_t"+to_string(r)+"_"+to_string(i));
                    }
                }
                --logfile.currIndent;
            }
            --logfile.currIndent;
        }

        --logfile.currIndent;
        logfile.write("Monte Carlo simulation complete");
        logfile.separator();
    }

    //Write total analysis
    for(int i=0; i<10000; ++i){
        VecF<double> hist(4);
        hist[0]=i*4.0/10000.0;
        hist[1]=lenHist[i];
        hist[2]=i*2*M_PI/10000.0;
        hist[3]=angHist[i];
        outGeomHist.writeRowVector(hist);
    }

    //Check network
    logfile.write("Simulation diagnostics");
    ++logfile.currIndent;
    bool consistent=network.checkConsistency();
    bool convex=network.checkConvexity();
    logfile.write("Network consistent:",consistent);
    logfile.write("Rings convex:",convex);
    logfile.write("Monte Carlo acceptance:",(double)accepted/mcSteps);
    logfile.write("Geometry optimisation codes:");
    ++logfile.currIndent;
    logfile.write("Converged: ",optCodes[0]);
    logfile.write("Converged (zero force): ",optCodes[1]);
    logfile.write("Unconverged: ",optCodes[2]);
    logfile.write("Failed (overlapping): ",optCodes[3]);
    logfile.write("Failed (initially non-convex): ",optCodes[4]);
    --logfile.currIndent;
    logfile.write("Geometry optimisation average iterations:",optIterations/vSum(optCodes));
    --logfile.currIndent;
    logfile.write("Diagnostics complete");
    logfile.separator();

    //Write files
    logfile.write("Writing files");
    ++logfile.currIndent;
    network.wrapCoordinates();
    network.syncCoordinates();
    network.write(prefixOut);
    logfile.write("Writing network files");
    --logfile.currIndent;
    logfile.write("Files written");
    logfile.separator();

    //Close files
    logfile.datetime("Simulation complete at: ");
    return 0;
}
