#include "monte_carlo.h"

//##### Metropolis-Hastings #####

//Default constructor
Metropolis::Metropolis() {
    mtGen.seed(0);
    rand01=uniform_real_distribution<double>(0.0,1.0);
    rTemperature=1.0;
    energyPrev=0.0;
}

//Construct with random seed, temperature and initial energy
Metropolis::Metropolis(int seed, double temperature, double energy) {
    mtGen.seed(seed);
    rand01=uniform_real_distribution<double>(0.0,1.0);
    if(temperature<=0.0) throw "Cannot initialise Metropolis algorithm with zero temperature";
    rTemperature=1.0/temperature;
    energyPrev=energy;
}

//Evaluate Metropolis condition, whether to accept or reject move
int Metropolis::acceptanceCriterion(double energy) {

    /* Metropolis algorithm efficiently samples Boltzmann distribution
     * 1) if move downhill in energy accept
     * 2) if move uphill accept with probabilitiy min[1,e^-de/t] */

    double deltaE=energy-energyPrev;
    if(deltaE<0.0){
        energyPrev=energy;
        return 1;
    }
    else{
        double probability=exp(-deltaE*rTemperature);
        double randomNum=rand01(mtGen);
        if(randomNum<probability){
            energyPrev=energy;
            return 1;
        }
        else return 0;
    }
}

