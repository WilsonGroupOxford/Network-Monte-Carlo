#include <iostream>
#include <sstream>
#include "outputfile.h"
#include "linked_network.h"

using namespace std;

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
    getline(inputFile,line);
    istringstream(line)>>nRings;
    getline(inputFile,line);
    istringstream(line)>>minRingSize;
    getline(inputFile,line);
    istringstream(line)>>maxRingSize;
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
    istringstream(line)>>randomSeed;
    getline(inputFile,skip);
    getline(inputFile,skip);
    logfile.write("Monte Carlo parameters read");
    //Energy search
    int mcSteps,equilSteps;
    double mcStartT,mcEndT,mcIncT;
    getline(inputFile,line);
    istringstream(line)>>mcStartT;
    getline(inputFile,line);
    istringstream(line)>>mcEndT;
    getline(inputFile,line);
    istringstream(line)>>mcIncT;
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
    double potAK,potBK,potCK;
    getline(inputFile,line);
    istringstream(line)>>potBK;
    getline(inputFile,line);
    istringstream(line)>>potAK;
    getline(inputFile,line);
    istringstream(line)>>potCK;
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
    int analysisFreq;
    getline(inputFile,line);
    istringstream(line)>>analysisFreq;
    logfile.write("Analysis parameters read");
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
    logfile.write("Monte Carlo run type:",runType);
    if(runType=="energy") {
        logfile.write("Monte carlo initial log temperature:", mcStartT);
        logfile.write("Monte carlo final log temperature:", mcEndT);
        logfile.write("Monte carlo log temperature increment:", mcIncT);
        logfile.write("Monte carlo steps per increment:", mcSteps);
        logfile.write("Equilibrium steps:", equilSteps);
    }
    else if(runType=="cost"){
        logfile.write("Monte carlo p lower limit:", costLowLimP);
        logfile.write("Monte carlo p upper limit:", costUpLimP);
        logfile.write("Monte carlo r lower limit:", costLowLimR);
        logfile.write("Monte carlo r upper limit:", costLowLimR);
        logfile.write("Monte carlo effective temperature:", costT);
        logfile.write("Monte carlo steps per increment:", costSteps);
    }
    LinkedNetwork network(nRings,lattice,4,maxRingSize,minRingSize);
    network.initialisePotentialModel(potAK,potBK,potCK);
    network.initialiseGeometryOpt(goptLocalIt,goptTau,goptTol,goptLocalSize);
    network.initialiseMonteCarlo(pow(10,mcStartT),randomSeed);
    network.initialiseCostFunction(costT,randomSeed,costPK,costRK);
    if(crystal!="default") network.makeCrystal(crystal,lattice);
    if(lattice=="goldberg" || lattice=="inv_cubic") network.optimalProjection("sphere");
    bool mixedLattice; //whether mixed coordination lattice
    if(lattice.substr(0,3)=="mix"){//get statistics on number of 3/4 coordinate nodes for mixed lattice
        mixedLattice=true;
        double mixA=network.networkA.nodes.n,mixA3=0,mixA4=0;
        for(int i=0; i<mixA; ++i){
            if(network.networkA.nodes[i].netCnxs.n==3) ++mixA3;
            else if(network.networkA.nodes[i].netCnxs.n==4) ++mixA4;
        }
        logfile.write("Mixed lattice total rings:",network.networkB.nodes.n);
        logfile.write("Mixed lattice total nodes:",mixA);
        logfile.write("Mixed lattice 3 coordinate nodes:",mixA3/mixA);
        logfile.write("Mixed lattice 4 coordinate nodes:",mixA4/mixA);
    }
    else mixedLattice=false;
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
    outGeometry.initVariables(6,4,60,20);
    outEmatrix.initVariables(1,4,60,int(log10(nRings*12))+2);
    outGeomHist.initVariables(6,4,60,20);
    logfile.write("Ring statistics file created");
    logfile.write("Correlations file created");
    logfile.write("Energy file created");
    logfile.write("Entropy file created");
    logfile.write("Temperature file created");
    logfile.write("Geometry file created");
    logfile.write("Geometry histogram file created");
    logfile.write("Edge distribution file created");
    --logfile.currIndent;
    logfile.write("Files initialised");
    logfile.separator();

    //Initialise total analysis variables - only update in main simulation not equilibration
    VecF<double> lenHist(10000),angHist(10000);
    lenHist=0.0;
    angHist=0.0;

    //Run monte carlo
    int accepted=0,optIterations=0;
    VecF<int> optCodes(4);
    optCodes=0;
    int trackFreq=100;
    VecF<int> moveStatus;
    if(runType=="energy") {//energy run
        //Run monte carlo equilibrium
        logfile.write("Running Monte Carlo equilibration");
        ++logfile.currIndent;
        double energy=network.mc.getEnergy();
        double mcT = pow(10, mcStartT);
        for (int i = 1; i <= equilSteps; ++i) {
            if(!mixedLattice) moveStatus = network.monteCarloSwitchMove(energy);
            else moveStatus = network.monteCarloMixMove(energy);
            accepted += moveStatus[0];
            optCodes[moveStatus[1]] += 1;
            optIterations += moveStatus[2];
//        cout << i << endl;
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
                corr[0] = r;
                corr[1] = aEst;
                corr[2] = aw[0];
                corr[3] = aw[1];
                corr[4] = aw[2];
                VecF<double> emptyL,emptyA; //dummy histograms
                VecF<double> geomStats = network.getOptimisationGeometry(emptyL,emptyA);
                VecF< VecF<int> > edgeDist = network.getEdgeDistribution("B");
                outRingStats.writeRowVector(ringStats);
                outCorr.writeRowVector(corr);
                outEnergy.write(energy);
                outEntropy.writeRowVector(s);
                outTemperature.write(mcT);
                outGeometry.writeRowVector(geomStats);
                for(int j=0; j<edgeDist.n; ++j) outEmatrix.writeRowVector(edgeDist[j]);
            }
        }
        --logfile.currIndent;
        logfile.write("Monte Carlo equilibration complete");
        logfile.separator();

        //Perform monte carlo simulation
        logfile.write("Running Monte Carlo simulation");
        ++logfile.currIndent;
        int nT = (mcEndT - mcStartT) / mcIncT;
        for (int t = 0; t <= nT; ++t) {
            mcT = pow(10, mcStartT + t * mcIncT);
            network.mc.setTemperature(mcT);
            logfile.write("Temperature:", mcT);
            ++logfile.currIndent;
            for (int i = 1; i <= mcSteps; ++i) {
                if(!mixedLattice) moveStatus = network.monteCarloSwitchMove(energy);
                else moveStatus = network.monteCarloMixMove(energy);
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
                    cout << t << " " << i << " " << accepted << endl;
                }
                if (i % analysisFreq == 0) {
                    VecF<double> ringStats = network.getNodeDistribution("B");
                    double r = network.getAssortativity("B");
                    double aEst = network.getAboavWeaireEstimate("B");
                    VecF<double> aw = network.getAboavWeaire("B");
                    VecF<double> s = network.getEntropy("B");
                    VecF<double> corr(5);
                    corr[0] = r;
                    corr[1] = aEst;
                    corr[2] = aw[0];
                    corr[3] = aw[1];
                    corr[4] = aw[2];
                    VecF<double> geomStats = network.getOptimisationGeometry(lenHist,angHist);
                    VecF< VecF<int> > edgeDist = network.getEdgeDistribution("B");
                    outRingStats.writeRowVector(ringStats);
                    outCorr.writeRowVector(corr);
                    outEnergy.write(energy);
                    outEntropy.writeRowVector(s);
                    outTemperature.write(mcT);
                    outGeometry.writeRowVector(geomStats);
                    for(int j=0; j<edgeDist.n; ++j) outEmatrix.writeRowVector(edgeDist[j]);
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
                corr[0] = r;
                corr[1] = aEst;
                corr[2] = aw[0];
                corr[3] = aw[1];
                corr[4] = aw[2];
                VecF<double> emptyL,emptyA; //dummy histograms
                VecF<double> geomStats = network.getOptimisationGeometry(emptyL,emptyA);
                VecF< VecF<int> > edgeDist = network.getEdgeDistribution("B");
                outRingStats.writeRowVector(ringStats);
                outCorr.writeRowVector(corr);
                outEnergy.write(energy);
                outEntropy.writeRowVector(s);
                outTemperature.write(costT);
                outGeometry.writeRowVector(geomStats);
                for(int j=0; j<edgeDist.n; ++j) outEmatrix.writeRowVector(edgeDist[j]);
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
                        corr[0] = r;
                        corr[1] = aEst;
                        corr[2] = aw[0];
                        corr[3] = aw[1];
                        corr[4] = aw[2];
                        VecF<double> geomStats = network.getOptimisationGeometry(lenHist,angHist);
                        VecF< VecF<int> > edgeDist = network.getEdgeDistribution("B");
                        outRingStats.writeRowVector(ringStats);
                        outCorr.writeRowVector(corr);
                        outEnergy.write(energy);
                        outEntropy.writeRowVector(s);
                        outTemperature.write(costT);
                        outGeometry.writeRowVector(geomStats);
                        for(int j=0; j<edgeDist.n; ++j) outEmatrix.writeRowVector(edgeDist[j]);
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