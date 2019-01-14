#include <iostream>
#include <sstream>
#include "outputfile.h"

using namespace std;

int main(){

    //Set up logfile
    Logfile logfile("./netmc.log");
    logfile.datetime("Simulation begun at: ");
    logfile.write("NETwork Monte Carlo");
    logfile.write("Written By: David OM, Wilson Group, 2019");
    logfile.separator();

    //Read input file
    ifstream inputFile("./netmc.inpt", ios::in);
    if(!inputFile.good()) logfile.criticalError("Cannot find input file 'netmc.inpt' in current directory");
    //I/O
    string skip,line,prefixOut; //prefix for output files
    getline(inputFile,skip);
    getline(inputFile,line);
    istringstream(line)>>prefixOut;
    getline(inputFile,skip);
    getline(inputFile,skip);
    //Network properties
    int nRings,minRingSize,maxRingSize;
    string lattice;
    getline(inputFile,line);
    istringstream(line)>>nRings;
    getline(inputFile,line);
    istringstream(line)>>minRingSize;
    getline(inputFile,line);
    istringstream(line)>>maxRingSize;
    getline(inputFile,line);
    istringstream(line)>>lattice;
    getline(inputFile,skip);
    getline(inputFile,skip);
    //Monte carlo
    int randomSeed;
    double mcTemperature;
    getline(inputFile,line);
    istringstream(line)>>randomSeed;
    getline(inputFile,line);
    istringstream(line)>>mcTemperature;
    getline(inputFile,skip);
    getline(inputFile,skip);
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
    inputFile.close();





    //Close files
    logfile.datetime("Simulation complete at: ");
    return 0;
}