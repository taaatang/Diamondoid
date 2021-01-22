#include <iostream>
#include <vector>
#include <string>
#include "src/utils.hpp"
#include "src_config/random.hpp"
#include "src_config/molecule.hpp"
#include "src_config/cluster.hpp"

const int N_mol = 16;
std::string dir = "data";
int main(){
    mkdir_fs(dir);
    std::string outfile = dir + "/AdamantaneCoords_"+std::to_string(N_mol)+".dat";
    std::vector<Adamantane> Adas(N_mol,Adamantane(0));
    Cluster clus(40,40,40);
    clus.init(Adas);
    int stepCount = 0;
    while(stepCount < MAX_MC_STEP){
        stepCount++;
        clus.singleStep();
        if(stepCount%10==0)std::cout<<"Step: "<<stepCount<<", Total Bonds: "<<clus.countBond()<<".\n\n";
    }
    clus.saveCoords(outfile);
    return 0;
}