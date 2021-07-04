#include "paras.hpp"

Parameters para("/Users/tatang/Documents/work/projects/PPP", {"input.txt"});

std::string dataDir = para.maps.at("data dir");
std::string molName = para.maps.at("molecule name");
std::string runid = para.maps.at("runid");
int molNum = para.mapi.at("molecule num");

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

double JUMP_LIMIT = para.mapd.at("jump limit");
