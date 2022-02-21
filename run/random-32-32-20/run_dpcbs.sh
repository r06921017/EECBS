#!/bin/bash
echo "Running Experiments: DPCBS"

name="random-32-32-20"
map="/home/rdaneel/mapf_benchmark/mapf-map/$name.map"
scen1="even"
scen="/home/rdaneel/mapf_benchmark/scen-$scen1/$name-$scen1"
output="/home/rdaneel/my_exp/$name/NEECBS/$name-$scen1"
sid=0
time=60
w=1.05

for n in $(seq 80 20 100)
do
    for i in $(seq 1 1 25)
    do
        echo "$n agents on instance $name-$scen1-$i  w=$w"
        ../../build/eecbs -m $map -a $scen-$i.scen -k $n -o $output-$n-$w-$sid-DPCBS.csv -t $time --suboptimality $w -s 1 --highLevelSolver "DPS" --inadmissibleH "Global" --heuristics "WDG" --prioritizingConflicts true --bypass true --rectangleReasoning true --corridorReasoning true --targetReasoning true
    done
done