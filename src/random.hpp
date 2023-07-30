#pragma once

#include <random>
#include <functional>
#include <chrono>

int diceI(int N);

double diceD();

void randomSeed(bool opt);

std::vector<int> FisherYatesShuffle(int n);

template<typename T>
void FisherYatesShuffle(std::vector<T>& v){
    int n = v.size();
    for(int i = 0; i < n - 1; i++){
        int j = diceI(n-i-1)+i;
        T tmp = v[i]; v[i] = v[j]; v[j] = tmp;
    }
}