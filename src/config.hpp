#ifndef __CONFIG_H__
#define __CONFIG_H__
#include <cmath>
#include <string>
#include "global.hpp"
const Vec3d region = {12., 12., 12.}; // box for MD simulation with PBC
const int MolNum = 27;
const int siteNum = 18;
const double dt = 0.001;
const double dt2 = dt * dt;
const int tstepNum = 15000;
const int tmeasureNum = 6000;
int stepCount = 0;
double temp = 0.45;
double tempI=0.5, tempF=0.1;
double dtemp;
const int stepChangeTemp = 100;
const int stepAdjustTemp = 10;
const int stepPrintInfo = 100;
double velMag = std::sqrt(3.*(1.-1./MolNum)*temp);
const double rc = 12.0;
const double rc2 = rc * rc;

bool isMeasure;
double rangeRdf = 10.0;
int sizeHistRdf = 100;
double RdfDeltaR;
const int stepRdf = 50;
const int limitRdf = 100;
int countRdf  = 0;

std::string DataPath = "Data/Triamantane_8";
std::string ETKFile = DataPath + "/ETK";
std::string ERKFile = DataPath + "/ERK";
std::string EVFile = DataPath + "/EV";
std::string PosFile = DataPath + "/Pos";
std::string HistRdfFile = DataPath + "/HistRdf";
std::string HistRdfDeltaRFile = DataPath + "/RdfDeltaR";
#endif // __CONFIG_H__