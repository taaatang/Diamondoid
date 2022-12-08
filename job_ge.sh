#!/usr/bin/env bash

# env=env.sh
# set environment
cwd=$(pwd)
# cd ~
# source ./${env}
# cd ${cwd}

for atomnum in 1000 5000 10000
do
for molname in Ada Dia Tria Tetra Penta Penta1212
   do
       echo "${molname} with ${atomnum} atoms"
       #atomnum=5000
       Ti=0.12
       Tf=0.024
       JobDir=Job/Longtime/Random/D6/${molname}/${atomnum}
       keyword1=RUNID
       keyword2=MOLNAME
       keyword3=ATOMNUM
       keyword4=T_INIT
       keyword5=T_FINAL
       app=diamondMC
       script=run_sherlock.sh
       input=input.txt

       for i in $(seq 0 1 9)
           do
               appDir=${JobDir}/run_${i}
               mkdir -p ${appDir}
               inputDir=${appDir}
               mkdir -p ${inputDir}

               sed -e "s/${keyword1}/${i}/g;s/${keyword2}/${molname}/g;s/${keyword3}/${atomnum}/g;s/${keyword4}/${Ti}/g;s/${keyword5}/${Tf}/g" ${input} > ${inputDir}/${input}
               sed -e "s/${keyword1}/${i}/g;s/${keyword2}/${molname}/g;s/${keyword3}/${atomnum}/g" ${script} > ${appDir}/${script}
               cd ${appDir}
               rm -rf job.out
               rm -rf job.err
               sbatch ${script}
               cd ${cwd}
           done
   done
done
