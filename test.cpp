#include <iostream>
#include "src/random.hpp"
#include "src/paras.hpp"
using namespace std;

int main ( ) {
    Parameters para("/Users/tatang/Documents/work/projects/PPP", {"input.txt"});
    randomSeed(para.mapi.at("random seed"));
    cout<<"test diceI:";
    for (int i = 0; i < 10; ++i) {
        cout<<diceI(10)<<" ";
    }
    cout<<"\n";

    cout<<"test diceD:";
    for (int i = 0; i < 10; ++i) {
        cout<<diceD()<<" ";
    }
    cout<<"\n"; 
    return 0;
}