#ifndef __MOLECULE_H__
#define __MOLECULE_H__

#include "global.hpp"
#include "config.hpp"
/* Units
    signma = 1
    epsilon = 1
    m = 1
    t0 = sqrt(m*sigma^2/epsilon)
    kB = 1
*/
const std::array<Vec3d,10> Adamantane = {Vec3d{0., 0.5, 0.5},\
                                       Vec3d{1., 0.5, 0.5},\
                                       Vec3d{0.5, 0.5, 0.},\
                                       Vec3d{0.5, 0.5, 1.},\
                                       Vec3d{0.5, 0., 0.5},\
                                       Vec3d{0.5, 1., 0.5},\
                                       Vec3d{0.75, 0.25, 0.75},\
                                       Vec3d{0.25, 0.25, 0.25},\
                                       Vec3d{0.25, 0.75, 0.75},\
                                       Vec3d{0.75, 0.75, 0.25}};
const Vec3d ucell = {2.0,2.0,2.0};
// Ix,Iy,Iz
Vec3d mInert;
// LJ interaction sites of a molecule
std::array<Vec3d,siteNum> LJSites = Adamantane;
std::array<Vec3d,siteNum> MassSites = Adamantane;

struct FSite{
    Vec3d f,r;
};

struct Molecule{
    Vec3d r, rv, ra, ra1, ra2, r0, rv0;
    Vec4d q, qv, qa, qa1, qa2, q0, qv0;
    Vec3d torq;
};

std::array<FSite,MolNum*siteNum> fsites;
std::array<Molecule,MolNum> mols;
double ETransK, ERotK, EPot;
// predictor-corrector coefficients with k = 4
constexpr Vec3d pr = {19./24., -10./24., 3./24.};
constexpr Vec3d pv = {27./24., -22./24., 7./24.};
constexpr Vec3d cr = {3./24., 10./24., -1./24.};
constexpr Vec3d cv = {7./24., 6./24., -1./24.};

#endif // __MOLECULE_H__