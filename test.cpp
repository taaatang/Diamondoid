#include <iostream>
#include <vector>
#include <string>
#include "src/utils.hpp"
#include "src_config/random.hpp"
#include "src_config/molecule.hpp"
#include "src_config/cluster.hpp"


const int MAX_MC_STEP = 100000;
constexpr int MEASURE_STEP = MAX_MC_STEP/1000;

const double INIT_TEMP = 2.0;
const double FINAL_TEMP = 0.1;
const double TEMP = 0.1;
const int TAdjust_STEP = 100;
double dT = (FINAL_TEMP - INIT_TEMP)/MAX_MC_STEP * TAdjust_STEP;

const int N_mol = 32;

int main(){
    Timer timer;
    std::string dir = "data/Adamantane_"+std::to_string(N_mol);
    mkdir_fs(dir);

    std::vector<Atom> Mol;
    timer.tik();
    std::vector<Adamantane> Mols(N_mol,Adamantane(0,&Mol));
    std::vector<Molecule*> Molsptr; for(auto& mol:Mols)Molsptr.push_back(&mol);

    Cluster clus(40,40,40);
    double temp = INIT_TEMP;
    clus.setTemp(temp);
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
        if(stepCount%TAdjust_STEP==0) {
            temp += dT;
            clus.setTemp(temp);
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

    timer.tik();
    clus.Latt.clustering();
    timer.tok();
    std::cout<<"Clustering time:"<<timer.elapse()/1000.0<<"s.\n";
    std::cout<<"Cluster tot size:"<<clus.Latt.latt.size()<<".\n";
    

    // clus.saveCoords(outfile);
    return 0;
}