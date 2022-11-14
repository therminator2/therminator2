#!python

import os
from shutil import copytree, copy

submissiondir="/lustre/hebe/hades/user/sudol/submit/Therminator3/4/"
outputdir="/lustre/hebe/hades/user/sudol/blabla/therminator_PML_minusError_NEW/"
jobarrayfile="jobarray.dat"
#
try:
    os.remove(jobarrayfile)
    os.mkdir(submissiondir)
    os.mkdir(outputdir)
except Exception, e:
    print e
#
#
file = open(jobarrayfile,"w")
for i in range(1,100):
    dir = "{}/part_{}".format(outputdir,i)
    try:
        os.mkdir(dir)
        os.mkdir(dir + "/fomodel")
        os.mkdir(dir + "/share")
        copy("therm2_events",dir)
    except Exception, e:
        print e
    try:
        replace_command = "perl replace.pl -t events_template -d events.ini -s {}".format(i)
       # replace_command = "perl replace.pl -t fomodel/blastwave_template -d fomodel/blastwave.ini -u {} -m {}".format(tau, rho_max)
        print replace_command
        os.system(replace_command)
        copy("events.ini",dir)
        copy("fomodel/SR.ini",dir + "/fomodel")
        copy("share/particles.data",dir + "/share")
        copy("share/decays.data",dir + "/share")
    except Exception, e:
        print e
    file.write(dir+"\n")

file.close()
copy(jobarrayfile, submissiondir)
copy("jobScript_SL.sh", submissiondir)

command="sbatch --array=1-100 --mem=1000 --time=0-6:00:00 -D {0}  --output={1}/out/slurm-%A_%a.out {0}/jobScript_SL.sh {0}/jobarray.dat {1}".format(submissiondir,outputdir)
print command
os.system(command)

