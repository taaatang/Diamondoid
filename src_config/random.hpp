#ifndef __RANDOM_H__
#define __RANDOM_H__

#include <random>
#include <functional>
std::mt19937 generatorI;
template<size_t N>
std::uniform_int_distribution<int> distributionI(0,N);

int diceI(int N){
    std::uniform_int_distribution<int> distributionI(0,N);
    return distributionI(generatorI);
}

std::mt19937 generatorD;
std::uniform_real_distribution<double> distributionD(0.0,1.0);
auto diceD = std::bind(distributionD,generatorD);


#endif // __RANDOM_H__