#ifndef __RANDOMS_H__
#define __RANDOMS_H__

#include "global.hpp"
#define IADD 453806245
#define IMUL 314159269
#define MASK 2147483647
#define SCALE 0.4656612873E-9
extern int randSeed;
double RandR();
Vec3d VRand();
#endif // __RANDOMS_H__