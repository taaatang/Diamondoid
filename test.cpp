#include <iostream>
#include <vector>
#include <string>
#include "src/utils.hpp"
#include "src_config/random.hpp"
#include "src_config/molecule.hpp"
#include "src_config/cluster.hpp"

const int MAX_MC_STEP = 10000;
constexpr int MEASURE_STEP = MAX_MC_STEP/100;
const int N_mol = 16;
int main(){
    Timer timer;
    std::string dir = "data/Pentamantane_"+std::to_string(N_mol);
    mkdir_fs(dir);

    std::vector<Atom> Mol;
    timer.tik();
    std::vector<Pentamantane> Mols(N_mol,Pentamantane(0,&Mol));
    std::vector<Molecule*> Molsptr; for(auto& mol:Mols)Molsptr.push_back(&mol);
    Cluster clus(40,40,40);
    clus.init(Molsptr);
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

    timer.tik();
    std::cout<<"\nBegin evaluate surface positions...\n";
    clus.computePos();
    clus.evalPos();
    std::string outfile = dir + "/surface";
    clus.saveSurface(outfile);
    std::cout<<"Surface evaluation time:"<<timer.elapse()/1000.0<<"s.\n";
    timer.tok();
    

    // clus.saveCoords(outfile);
    return 0;
}