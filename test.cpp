#include <iostream>
#include <vector>
#include <string>
#include "src/utils.hpp"
#include "src_config/random.hpp"
#include "src_config/molecule.hpp"
#include "src_config/cluster.hpp"

const int MAX_MC_STEP = 10000;
constexpr int MEASURE_STEP = MAX_MC_STEP/100;
const int N_mol = 512;
int main(){
    Timer timer;
    std::string dir = "data/Adamantane_"+std::to_string(N_mol);
    mkdir_fs(dir);

    timer.tik();
    std::vector<Adamantane> Adas(N_mol,Adamantane(0));
    Cluster clus(40,40,40);
    clus.init(Adas);
    timer.tok();
    std::cout<<"Initialization time:"<<timer.elapse()/1000.0<<"s.\n\n";

    timer.tik();
    std::vector<int> bondNum;
    int stepCount = 0;
    while(stepCount < MAX_MC_STEP){
        stepCount++;
        clus.singleStep();
        if(stepCount%MEASURE_STEP==0){
            int mstep = stepCount / MEASURE_STEP;
            std::string outfile = dir + "/step_"+std::to_string(mstep)+".dat";
            auto bnum = clus.countBond();
            bondNum.push_back(bnum);
            clus.saveCoords(outfile);
            std::cout<<"Step: "<<stepCount<<", Total Bonds: "<<bnum<<".\n";
        }
    }
    std::ofstream outf;
    save<int>(bondNum.data(),bondNum.size(),&outf,dir+"/bondNum.dat");
    timer.tok();
    std::cout<<"\nMonte Carlo time:"<<timer.elapse()/1000.0<<"s.\n";

    // clus.saveCoords(outfile);
    return 0;
}