#!/usr/bin/env bash

# env=env.sh
# set environment
cwd=$(pwd)
# cd ~
# source ./${env}
# cd ${cwd}

molname=Penta1212
molnum=12
Ti=2.0
Tf=0.1
JobDir=Job/${molname}/${molnum}
keyword1=RUNID
keyword2=MOLNAME
keyword3=MOLNUM
keyword4=T_INIT
keyword5=T_FINAL
app=main.out
script=run_sherlock.sh
input=input.txt

for i in $(seq 0 1 9)
    do
        appDir=${JobDir}/run_${i}
        mkdir -p ${appDir}
        inputDir=${appDir}
        mkdir -p ${inputDir}

        sed -e "s/${keyword1}/${i}/g;s/${keyword2}/${molname}/g;s/${keyword3}/${molnum}/g;s/${keyword4}/${Ti}/g;s/${keyword5}/${Tf}/g" ${input} > ${inputDir}/${input}
        sed -e "s/${keyword1}/${i}/g;s/${keyword2}/${molname}/g;s/${keyword3}/${molnum}/g" ${script} > ${appDir}/${script}
        cd ${appDir}
        rm -rf job.out
        rm -rf job.err
        sbatch ${script}
        cd ${cwd}
    done
