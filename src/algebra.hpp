#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__

#include <iostream>
#include <array>
#include <cmath>

template <size_t N>
using arr = std::array<int,N>;

template <size_t N>
std::ostream& operator<<(std::ostream& os, const arr<N>& arr_out){
    os<<"[ ";
    for(int i = 0; i < N; i++) os << arr_out[i] << " ";
    os<<"]\n";
    return os;
}

template <size_t N>
arr<N> operator+(const arr<N>& lhs, const arr<N>& rhs){
    arr<N> res;
    for(size_t i = 0; i < N; i++) res[i] = lhs[i] + rhs[i];
    return res;
}

template <size_t N>
arr<N> operator-(const arr<N>& lhs, const arr<N>& rhs){
    arr<N> res;
    for(int i = 0; i < N; i++) res[i] = lhs[i] - rhs[i];
    return res;
}

template <size_t N>
double operator*(const arr<N>& lhs, const arr<N>& rhs){
    double res = 0.0;
    for(int i = 0; i < N; i++) res += lhs[i] * rhs[i];
    return res;
}

template <size_t N>
arr<N> ElProd(const arr<N>& lhs, const arr<N>& rhs){
    arr<N> res;
    for(int i=0; i<N; i++){
        res[i] = lhs[i] * rhs[i];
    }
    return res;
}

template <size_t N>
double norm(const arr<N>& arr_in){
    return std::sqrt(arr_in * arr_in);
}

template <size_t N>
arr<N> operator*(const arr<N>& lhs, int rhs){
    arr<N> res;
    for(int i = 0; i < N; i++) res[i] = lhs[i] * rhs;
    return res;
}

template <size_t N>
arr<N> operator*(int lhs, const arr<N>& rhs){
    arr<N> res;
    for(size_t i = 0; i < N; i++) res[i] = rhs[i] * lhs;
    return res;
}


template <size_t N>
arr<N> operator/(const arr<N>& lhs, int rhs){
    arr<N> res;
    for(int i = 0; i < N; i++) res[i] = lhs[i] / rhs;
    return res;
}
#endif // __ALGEBRA_H__