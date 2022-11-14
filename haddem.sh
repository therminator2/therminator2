#!/bin/bash

for Hubb in {70..100..5}
do
    for delt in {0..3}
    do  
        dirnameIn=H${Hubb}_D${delt}_file*
        dirnameOut=H${Hubb}D${delt}
        
        hadd ${dirnameOut}/event000.root  ${dirnameIn}/sr/event000.root
    done
done