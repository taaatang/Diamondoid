#pragma once

#include <iostream>
#include <array>
#include <cmath>

template <size_t N>
using arr = std::array<int,N>;

template <size_t N>
std::ostream& operator<<(std::ostream& os, const arr<N>& arr_out){
    os<<"[ ";
    for(size_t i = 0; i < N; i++) os << arr_out[i] << " ";
    os<<"]";
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
    for(size_t i = 0; i < N; i++) res[i] = lhs[i] - rhs[i];
    return res;
}

template <size_t N>
double operator*(const arr<N>& lhs, const arr<N>& rhs){
    double res = 0.0;
    for(size_t i = 0; i < N; i++) res += lhs[i] * rhs[i];
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

template <size_t N>
void min(arr<N> &minimum, const arr<N> &arr_in) {
    for (size_t i = 0; i < N; ++i) {
        if (arr_in[i] < minimum[i]) {
            minimum[i] = arr_in[i];
        }
    }
}

template <size_t N>
void max(arr<N> &maximum, const arr<N> &arr_in) {
    for (size_t i = 0; i < N; ++i) {
        if (arr_in[i] > maximum[i]) {
            maximum[i] = arr_in[i];
        }
    }
}

template <size_t N>
double getVolume(const arr<N> &box) {
    if (N == 0) {
        return 0.0;
    }
    double v = 1.0;
    for (size_t i = 0; i < N; ++i) {
        v *= double(box[i]);
    }
    return v;
}