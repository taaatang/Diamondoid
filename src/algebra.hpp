#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__

#include <array>
#include <iostream>

template <size_t N>
using arr = std::array<int,N>;

template <size_t N>
std::ostream& operator<<(std::ostream& os, const arr<N>& a){
    os<<"[ ";
    for(int i = 0; i < N; i++) os<<a[i]<<" ";
    os<<"]\n";
    return os;
}

template <size_t N>
arr<N> operator+(const arr<N>& a, const arr<N>& b){
    arr<N> c;
    for(int i = 0; i < N; i++) c[i] = a[i] + b[i];
    return c;
}

template <size_t N>
arr<N> operator-(const arr<N>& a, const arr<N>& b){
    arr<N> c;
    for(int i = 0; i < N; i++) c[i] = a[i] - b[i];
    return c;
}

template <size_t N>
double operator*(const arr<N>& a, const arr<N>& b){
    double c = 0.0;
    for(int i = 0; i < N; i++) c += a[i] * b[i];
    return c;
}

template <size_t N>
arr<N> ElProd(const arr<N>& a, const arr<N>&b){
    arr<N> c;
    for(int i=0; i<N; i++){
        c[i] = a[i] * b[i];
    }
    return c;
}

template <size_t N>
double norm(const arr<N>& a){
    return std::sqrt(a*a);
}

template <size_t N>
arr<N> operator*(const arr<N>& a, int b){
    arr<N> c;
    for(int i = 0; i < N; i++) c[i] = a[i] * b;
    return c;
}

template <size_t N>
arr<N> operator*(int b, const arr<N>& a){
    arr<N> c;
    for(int i = 0; i < N; i++) c[i] = a[i] * b;
    return c;
}


template <size_t N>
arr<N> operator/(const arr<N>& a, int b){
    arr<N> c;
    for(int i = 0; i < N; i++) c[i] = a[i] / b;
    return c;
}
#endif // __ALGEBRA_H__