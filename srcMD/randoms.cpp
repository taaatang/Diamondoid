#include "randoms.hpp"

int randSeed = 123;
double RandR(){
    randSeed = (randSeed * IMUL + IADD) & MASK;
    return (randSeed * SCALE);
}

Vec3d VRand(){
    double s{2.};
    Vec3d v;
    while(s>1.){
        v[0] = 2. * RandR() - 1.;
        v[1] = 2. * RandR() - 1.;
        s = v[0]*v[0] + v[1]*v[1];
    }
    v[2] = 1. - 2.*s;
    s = 2.*std::sqrt(1.-s);
    v[0] *= s;
    v[1] *= s;
    return v;
}