#ifndef __RANDOM_H__
#define __RANDOM_H__

#include <random>
#include <functional>
#include <chrono>

std::mt19937_64 generatorI;
template<size_t N>
std::uniform_int_distribution<int> distributionI(0,N);

int diceI(int N){
    std::uniform_int_distribution<int> distributionI(0,N);
    return distributionI(generatorI);
}

std::mt19937_64 generatorD;
std::uniform_real_distribution<double> distributionD(0.0,1.0);
double diceD( ) {
    return distributionD(generatorD);
}
// auto diceD = std::bind(distributionD,generatorD);

void randomSeed(bool opt) {
    if (opt) {
        unsigned seedI = std::chrono::system_clock::now().time_since_epoch().count();
        generatorI.seed(seedI);
        unsigned seedD = std::chrono::system_clock::now().time_since_epoch().count();
        generatorD.seed(seedD);
    }
}

std::vector<int> FisherYatesShuffle(int n){
    std::vector<int> v(n,0);
    for(int i = 0; i < n; i++) v[i] = i;
    for(int i = 0; i < n - 1; i++){
        int j = diceI(n-i-1) + i;
        int tmp = v[i]; v[i] = v[j]; v[j] = tmp;
    }
    return v;
}

template<typename T>
void FisherYatesShuffle(std::vector<T>& v){
    int n = v.size();
    for(int i = 0; i < n - 1; i++){
        int j = diceI(n-i-1)+i;
        T tmp = v[i]; v[i] = v[j]; v[j] = tmp;
    }
}

#endif // __RANDOM_H__