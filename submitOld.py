#!python

import os
from shutil import copytree, copy
import numpy as np

submissiondir="/lustre/hades/user/kjedrzej/submit/HubbDeltTest1"
outputdir="/lustre/hades/user/kjedrzej/HubbDeltTest1"
jobarrayfile="jobarray.dat"

try:
    os.remove(jobarrayfile)
    os.mkdir(submissiondir)
    os.mkdir(outputdir)
except Exception, e:
    print e
    

file = open(jobarrayfile,"w")
for i in range(7, 40):
   # tau = i * 0.1
    tau = i
    for j in range(110, 130):
        rho_max = j * 0.1
        dir = "{}/tau_{}_rhomax_{}_para".format(outputdir,tau,rho_max)
        print(dir)
        try:
            os.mkdir(dir)
            os.mkdir(dir + "/fomodel")
            os.mkdir(dir + "/share")
            copy("therm2_events",dir)
            copy("events.ini",dir)
           # replace_command = "perl replace.pl -t fomodel/blastwave_template -d fomodel/blastwave.ini -u {}".format(tau)
            replace_command = "perl replace.pl -t fomodel/SR_template -d fomodel/SR.ini -u {} -m {}".format(tau, rho_max)
            print replace_command
            os.system(replace_command)
            copy("fomodel/blastwave.ini",dir + "/fomodel")
            copy("share/particles.data",dir + "/share")
            copy("share/decays.data",dir + "/share")
        except Exception, e:
            print e
        file.write(dir+"\n")

file.close()
copy(jobarrayfile, submissiondir)
copy("jobScript_SL.sh", submissiondir)

command="sbatch --array=1-441 --mem=1000 --time=0-2:00:00 -D {0}  --output={1}/out/slurm-%A_%a.out -- {0}/jobScript_SL.sh {0}/jobarray.dat {1}".format(submissiondir,outputdir)
print command
os.system(command)

