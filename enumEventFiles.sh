#!/bin/bash

lustrDir="/u/kjedrzej/lustre/hades/user/kjedrzej/HubbDeltLowT"

Harr=(80)
Earr=(2 3) 
Darr=(0)
NOFiles=1000
declare -a outDirName
declare -a inDirName
declare -a ofn
declare -a iter

for H in "${Harr[@]}"
do
    for E in "${Earr[@]}"
    do
        for D in "${Darr[@]}"
        do
            outDirName+=("${lustrDir}/H${H}E${E}D${D}femto")
            if [ ! -d "${outDirName[-1]}" ];
            then
                mkdir ${outDirName[-1]}
            fi
            for (( i=1; i <= ${NOFiles}; i++ ))
            do
                iter+=($((i - 1)))
                inDirName+=("${lustrDir}/H${H}_E${E}_D${D}_file${i}")
                ofn+=("${outDirName[-1]}/event$(printf "%03d" ${iter[-1]}).root")
                if [ -f "${ofn[-1]}" ];
                then
                    echo "${ofn[-1]} exists"
                    continue
                fi
                cp -v "${inDirName[-1]}/sr/event000.root" ${ofn[-1]}
            done
        done
    done
done

unset outDirName
unset inDirName
unset ofn
unset iter
