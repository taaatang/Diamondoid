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