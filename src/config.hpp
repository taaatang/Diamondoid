#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "paras.hpp"

Parameters para(".", {"input.txt"});

std::string dataDir = para.maps.at("data dir");
std::string molName = para.maps.at("molecule name");
std::string runid = para.maps.at("runid");
//int molNum = para.mapi.at("molecule num");
int atomNum = para.mapi.at("atom num");

bool isRandseed = para.mapi.at("random seed");

int HEAT_STEP = para.mapi.at("heat step");
int ANNEAL_STEP = para.mapi.at("anneal step");
int COOL_STEP = para.mapi.at("cool step");
int MEASURE_STEP = para.mapi.at("measure step");
int TAdjust_STEP = para.mapi.at("TAdjust step");
double INIT_TEMP = para.mapd.at("init T");
double FINAL_TEMP = para.mapd.at("final T");
int TOT_STEP = HEAT_STEP + ANNEAL_STEP + COOL_STEP;
int PRINT_STEP = (HEAT_STEP + ANNEAL_STEP + COOL_STEP)/10;

double PRESSURE = para.mapd.at("pressure");

double JUMP_LIMIT = para.mapd.at("jump limit");

int roundCast(double positiveNum) {
    return int(positiveNum + 0.5);
}

int getMolNum(const std::string& mol, int atomN) {
    if (mol == "Ada") {
        return roundCast(atomN/10.0);
    } else if (mol == "Dia") {
        return roundCast(atomN/14.0);
    } else if (mol == "Tria") {
        return roundCast(atomN/18.0);
    } else if (mol == "Tetra") {
        return roundCast(atomN/22.0);
    } else if (mol == "Penta") {
        return roundCast(atomN/26.0);
    } else if (mol == "Penta1212") {
        return roundCast(atomN/26.0);
    } else {
        std::cout << "molecule " << mol << "mantane not defined!\n";
        exit(1);
    }
}

#endif // __CONFIG_H__
