#include <iostream>
#include <vector>
#include <string>
#include "utils/utils.hpp"
#include "src/random.hpp"
#include "src/molecule.hpp"
#include "src/cluster.hpp"
#include "src/config.hpp"

int main() {
    Timer timer;
    std::string dir = dataDir + "/" + molName + "mantane_" + std::to_string(molNum) + "/run" + runid;
    mkdir_fs(dir);
    para.print(dir + "/para.txt");
    randomSeed(isRandseed);

    // temperature
    auto temp = [=](int step) {
        if (step < HEAT_STEP) return INIT_TEMP;
        else if (step > HEAT_STEP + ANNEAL_STEP) return FINAL_TEMP;
        else return INIT_TEMP + (FINAL_TEMP - INIT_TEMP) / ANNEAL_STEP * (step-HEAT_STEP); 
        // else return INIT_TEMP + (FINAL_TEMP - INIT_TEMP) * (1.0 - std::exp((double)(step-HEAT_STEP)/COOL_STEP))/(1.0 - std::exp(1.0));
    };

    std::vector<Atom> Mol;
    timer.tik();
    std::vector<Molecule*> Molsptr; 
    for (int i = 0; i < molNum; ++i) {
        if (molName == "Ada") {
            Molsptr.push_back(new Adamantane(0, &Mol));
        } else if (molName == "Tria") {
            Molsptr.push_back(new Triamantane(0, &Mol));
        } else if (molName == "Tetra") {
            Molsptr.push_back(new Tetramantane(0, &Mol));
        } else if (molName == "Penta") {
            Molsptr.push_back(new Pentamantane(0, &Mol));
        } else {
            std::cout<<"molecule "<<molName<<"mantane not defined!\n";
            exit(1);
        }
    }
    std::cout<<molNum<<" "<<molName<<"mantanes are created!\n";

    // INIT
    Cluster clus(40,40,40);
    clus.setTemp(temp(0));
    clus.init(Molsptr);
    timer.tok();
    std::cout<<"Initialization time:"<<timer.elapse()/1000.0<<"s.\n\n";

    // MC
    timer.tik();
    std::vector<int> bondNum;
    for (int stepCount = 1; stepCount <= TOT_STEP; ++stepCount){

        clus.singleStep();

        if(stepCount%MEASURE_STEP==0){
            int mstep = stepCount / MEASURE_STEP;
            std::string outfile = dir + "/step_"+std::to_string(mstep)+".dat";
            auto bnum = clus.countBond();
            bondNum.push_back(bnum);
            clus.saveCoords(outfile);

            clus.computePos();
            clus.evalPos();
            clus.saveSurface(dir + "/step_"+std::to_string(mstep) + "_surface");

            clus.clustering();
            clus.saveVacancy(dir + "/step_"+std::to_string(mstep) + "_vacancy.dat");
        }

        if(stepCount%PRINT_STEP==0){
            auto bnum = clus.countBond();
            std::cout<<"Step: "<<stepCount<<" / "<<TOT_STEP<<", T:"<<temp(stepCount)<<", Total Bonds: "<<bnum<<".\n";
        }

        if(stepCount%TAdjust_STEP==0) {
            clus.setTemp(temp(stepCount));
        }

    }

    // save measured bond num
    std::ofstream outf;
    save<int>(bondNum.data(),bondNum.size(),&outf,dir+"/bondNum.dat");
    timer.tok();
    std::cout<<"\nMonte Carlo time:"<<timer.elapse()/1000.0<<"s.\n";

    // evaluate compatible/uncompatible surface positions
    // timer.tik();
    // std::cout<<"\nBegin evaluate surface positions...\n";
    // clus.computePos();
    // clus.evalPos();
    // clus.saveSurface(dir + "/surface");
    // std::cout<<"Surface evaluation time:"<<timer.elapse()/1000.0<<"s.\n";
    // timer.tok();

    // clustering occupied/unoccupied sites. discovering vacancy
    // timer.tik();
    // clus.clustering();
    // clus.saveVacancy(dir + "/vacancy.dat");
    // timer.tok();
    // std::cout<<"Clustering time:"<<timer.elapse()/1000.0<<"s.\n";
    // clus.saveCoords(outfile);

    for (auto& p : Molsptr) {
        delete p;
    }
    return 0;
}