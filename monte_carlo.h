//Monte Carlo Methods
#ifndef NL_MONTECARLO_H
#define NL_MONTECARLO_H

#include <iostream>
#include <random>

using namespace std;

////Metropolis-Hastings algorithm
class Metropolis{

private:

    //Data members
    mt19937 mtGen; //mersenne twister generator
    uniform_real_distribution<double> rand01; //uniform distribution
    double rTemperature; //reciprocal temperature
    double energyPrev; //previous energy

public:

    //Constructors
    Metropolis();
    Metropolis(int seed, double temperature, double energy=0.0);

    //Member functions
    int acceptanceCriterion(double energy);
};

#endif //NL_MONTECARLO_H
