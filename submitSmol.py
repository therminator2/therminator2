#!python

import os
from shutil import copytree, copy
import numpy as np

scriptdir = "/lustre/hades/user/kjedrzej/defallRoot6Hydra2-5.0.sh"
submissiondir="/lustre/nyx/hades/user/kjedrzej/submit/HubbDeltTest1"
outputdir="/lustre/nyx/hades/user/kjedrzej/HubbDeltTest1"
jobarrayfile="jobarray.dat"

try:
    os.remove(jobarrayfile)
    os.mkdir(submissiondir)
    os.mkdir(outputdir)
except Exception, e:
    print e
    

file = open(jobarrayfile,"w")
for run in range(0,100):
    Hubb = 0.01
    newH = str(Hubb*1000)
    newH = newH.replace('.','')
    delta = 0.0
    newD = int(delta*11)
    dir = "{}/H{}_D{}_file{}".format(outputdir,newH,newD,run)
    print(dir)
    try:
        os.mkdir(dir)
        os.mkdir(dir + "/fomodel")
        os.mkdir(dir + "/share")
        copy("therm2_events",dir)
        copy("events.ini",dir)
        replace_command = "perl replace.pl -t fomodel/SR_template -d fomodel/SR.ini -u {} -m {}".format(Hubb, delta)
        print replace_command
        os.system(replace_command)
        copy("radius_Np_124.xml",dir)
        copy("fomodel/SR.ini",dir + "/fomodel")
        copy("share/particles.data",dir + "/share")
        copy("share/decays.data",dir + "/share")
    except Exception, e:
        print e
    file.write(dir+"\n")

file.close()
copy(jobarrayfile, submissiondir)
copy("jobScript_SL.sh", submissiondir)
chmod_command = "chmod 777 " + submissiondir + "/jobScript_SL.sh"
print chmod_command
os.system(chmod_command)

command="sbatch --array=1-101 --mem=1000 --time=0-3:00:00 -D {0}  --output={1}/out/slurm-%A_%a.out -- {0}/jobScript_SL.sh {0}/jobarray.dat {1} {2}".format(submissiondir,outputdir,scriptdir)
print command
os.system(command)

