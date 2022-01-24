#!/bin/bash
echo "Running Experiments: NEECBS"

name="random-32-32-20"
map="/home/rdaneel/mapf_benchmark/mapf-map/$name.map"
scen1="even"
scen="/home/rdaneel/mapf_benchmark/scen-$scen1/$name-$scen1"
output="/home/rdaneel/my_exp/$name/NEECBS/$name-$scen1"
sid=0
time=60
w=1.05
b=25

for n in $(seq 20 20 100)
do
    for i in $(seq 1 1 25)
    do
        echo "$n agents on instance $name-$scen1-$i  w=$w  b=$b sid=$sid"
        ../../build/eecbs -m $map -a $scen-$i.scen -k $n -o $output-$n-$w-$sid-$b-NEECBS_WDG_NONE.csv -t $time --suboptimality $w -b $b -s 1 --highLevelSolver "EES" --inadmissibleH "Global" --heuristics "WDG" --prioritizingConflicts false --bypass true --rectangleReasoning false --corridorReasoning false --targetReasoning false
    done
done