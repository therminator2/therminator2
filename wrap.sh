#!/bin/bash
jobscript=$1
jobarrayFile=$2
pathoutputlog=$3
arrayoffset=$4

singularity exec \
-B /cvmfs/hadessoft.gsi.de/install/debian8/install:/cvmfs/hades.gsi.de/install \
-B /cvmfs/hadessoft.gsi.de/install/debian8/oracle:/cvmfs/it.gsi.de/oracle \
-B /lustre \
/cvmfs/vae.gsi.de/debian8/containers/user_container-production.sif  ${jobscript} ${jobarrayFile} ${pathoutputlog} ${arrayoffset}

