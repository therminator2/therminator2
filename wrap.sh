#!/bin/bash
jobscript=$1
jobarrayFile=$2
pathoutputlog=$3
arrayoffset=$4

singularity exec \
-B /cvmfs/hadessoft.gsi.de/install/debian10:/cvmfs/hadessoft.gsi.de/install/debian10 \
-B /cvmfs/hadessoft.gsi.de/install/debian10/oracle:/cvmfs/it.gsi.de/oracle \
-B /etc/slurm \
-B /var/run/munge \
-B /var/spool/slurm \
-B /var/lib/sss/pipes/nss \
-B /cvmfs/vae.gsi.de \
-B /lustre \
/cvmfs/vae.gsi.de/debian10/containers/user_container-production.sif  ${jobscript} ${jobarrayFile} ${pathoutputlog} ${arrayoffset}

