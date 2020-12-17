#ifndef __CONFIG_H__
#define __CONFIG_H__
#include <cmath>
#include <string>
#include "global.hpp"
const Vec3d region = {10., 10., 10.}; // box for MD simulation with PBC
const int MolNum = 50;
const int siteNum = 10;
const double dt = 0.01;
const double dt2 = dt * dt;
const int tstepNum = 10000;
int stepCount = 0;
const double temp = 1.0;
const int stepAdjustTemp = 20;
const int stepPrintInfo = 100;
const double velMag = std::sqrt(3.*(1.-1./MolNum)*temp);
const double rc = 6.0;
const double rc2 = rc * rc;

std::string DataPath = "Data";
std::string ETKFile = DataPath + "/ETK";
std::string ERKFile = DataPath + "/ERK";
std::string EVFile = DataPath + "/EV";
std::string PosFile = DataPath + "/Pos";
#endif // __CONFIG_H__