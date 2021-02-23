#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__

#include "global.hpp"

template <size_t N>
using arr = std::array<double,N>;

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
arr<N> operator*(const arr<N>& a, double b){
    arr<N> c;
    for(int i = 0; i < N; i++) c[i] = a[i] * b;
    return c;
}

template <size_t N>
arr<N> operator/(const arr<N>& a, double b){
    arr<N> c;
    for(int i = 0; i < N; i++) c[i] = a[i] / b;
    return c;
}

inline Vec4d QuatMul(Vec4d q2, Vec4d q3){
    Vec4d q1;
    q1[0] = q2[3]*q3[0] - q2[2]*q3[1] + q2[1]*q3[2] + q2[0]*q3[3];
    q1[1] = q2[2]*q3[0] + q2[3]*q3[1] - q2[0]*q3[2] + q2[1]*q3[3];
    q1[2] = -q2[1]*q3[0] + q2[0]*q3[1] + q2[3]*q3[2] + q2[2]*q3[3];
    q1[3] = -q2[0]*q3[0] - q2[1]*q3[1] - q2[2]*q3[2] + q2[3]*q3[3];
    return q1;
}

inline Vec3d Cross(Vec3d v1, Vec3d v2){
    Vec3d v;
    v[0] = v1[1]*v2[2] - v1[2]*v2[1];
    v[1] = v1[2]*v2[0] - v1[0]*v2[2];
    v[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return v;
}

inline Vec3d MV(const RMat& R, const Vec3d& v){
    Vec3d r = {0.0, 0.0, 0.0};
    for(int j=0; j < 3; j++){
        for(int i=0; i<3; i++)
        r[i] += R[i+3*j]*v[j];
    }
    return r;
}

inline bool isDiag(const RMat& mat){
    double tol = 1e-10;
    for(int j=0; j<3; j++){
        for(int i=0; i<3; i++){
            if(i!=j){
                if(std::abs(mat[j*3+i])>tol)return false;
            }
        }
    }
    return true;
}

inline Vec3d getDiag(const RMat& mat){
    Vec3d d = {mat[0], mat[4], mat[8]};
    return d;
}

#endif // __ALGEBRA_H__