#ifndef __MD_H__
#define __MD_H__
#include "math.h"
#include "../utils/utils.hpp"
#include "global.hpp"
#include "algebra.hpp"
#include "randoms.hpp"


RMat ComputeInert();
void ConstructMol();
void InitCoords();
void InitVels();
void InitAccels();
void InitAngCoords();
Vec4d EulerToQuat(const Vec3d& ang);
void InitAngVels();
void InitAngAccels();
void Initialization();
void InitMeasure();
// LJ Potential
double PairPotential(Vec3d r1, Vec3d r2);
// LJ Force
Vec3d PairForce(Vec3d r1, Vec3d r2);
RMat BuildRotMat(const Vec4d& q, bool transpose);
void GenSiteCoords();
void ComputeForces();
// build rotation matrix
void ComputeTorqs();
// Compute angular veclocity of a molecule
Vec3d ComputeW(int n);
// compute accelaration of quat: q''.
void ComputeQa();
void PR(int n);
void PRV(int n);
void CR(int n);
void CRV(int n);
void PQ(int n);
void PQV(int n);
void CQ(int n);
void CQV(int n);
void PredictorStep();
void CorrectorStep();
void PredictorStepQ();
void CorrectorStepQ();
// PBC(assumes molecules do not move too far in a single step)
void ApplyBoundaryCond();
void AdjustQuat();
// compute potential energy
double ComputeV();
// compute translational kinetic energy
double ComputeTK();
// compute rotational kinetic energy
double ComputeRK();
// Temperature Control
void ApplyThermostat();
void AdjustTemp();
void ChangeTemp(double T);
void EvalProps();
void save(bool isapp);
void saveRdf();
void SingleStep();
void EvalRdf();
bool NotFinished();
bool NotFinishedMeasure();
bool isMeasureRdf();
void printInfo();

#endif // __MD_H__