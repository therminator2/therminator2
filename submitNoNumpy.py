#!python

import os
from shutil import copytree, copy, rmtree

#----------------set input/output dirs here--------------#
unigendir="/lustre/hades/user/harabasz/unigen_new/UniGen-2.3/set.sh"
scriptdir = "/cvmfs/hadessoft.gsi.de/install/debian10/root-6.24.02/bin/thisroot.sh"
submissiondir="/lustre/hades/user/harabasz/submit/HubbDeltHigT_200steps05_noDoubleMomUpdate"
outputdir="/lustre/hades/user/harabasz/HubbDeltHigT_200steps05_noDoubleMomUpdate"
jobarrayfile="jobarray.dat"

#--------------------------------------------------------#

def gen_primes():   #Sieve of Erastothenes
    D = {}
    q = 2
    
    while True:
        if q not in D:
            yield q
            D[q * q] = [q]
        else:
            for p in D[q]:
                D.setdefault(p + q, []).append(p)
            del D[q]
        
        q += 1

def replaceSR(hVal,eVal,dVal):
    fin = open("fomodel/SR_template", "rt")
    fout = open("fomodel/SR.ini", "wt")
    word = []

    for line in fin:
        if "HUBB" in line:
            fout.write(line.replace('HUBB', str(hVal)))
        elif "EPSI" in line:
            fout.write(line.replace('EPSI', str(eVal)))
        elif "DELT" in line:
            fout.write(line.replace('DELT', str(dVal)))
        elif "RadiusFile" in line:
            fout.write(line)
            word = line.split()
        else:
            fout.write(line)
    fin.close()
    fout.close()
    return word[2]  #returns radius file name

def replaceEV(sVal):
    fin = open("events_template", "rt")
    fout = open("events.ini", "wt")

    for line in fin:
        if "SEED" in line:
            fout.write(line.replace('SEED', str(sVal)))
        else:
            fout.write(line)
    fin.close()
    fout.close()

try:
    rmtree(submissiondir)
    os.mkdir(submissiondir)
    os.makedirs(outputdir+"/out", exist_ok=True)
    os.remove(jobarrayfile)
except Exception as e:
    print(e)
    
mygeneraotr = gen_primes()

file = open(jobarrayfile,"w")

#--------------change parameters here-------------------#

HubbList = [0.0225]
epsilonList = [0.0]
deltalist = [0.4]

NoFiles = 200

#-------------------------------------------------------#

for Hubb in HubbList:
    for epsilon in epsilonList:
        for delta in deltalist:
            for i,prime in enumerate(mygeneraotr):
                if i > NoFiles-1:
                    break
                newH = str(Hubb*1000)
                newH = newH.replace('.','')
                newE = abs(int(epsilon*11)) #one would think that using *10 instead would do the trick... but it just won't cut it for some reason
                newD = int(delta*11)
                dir = "{}/H{}_E{}_D{}_file{}".format(outputdir,newH,newE,newD,i+1)
                print(dir)
                
                try:
                    os.mkdir(dir)
                    os.mkdir(dir + "/fomodel")
                    os.mkdir(dir + "/share")
                    copy("therm2_events",dir)
                    radiusFile = replaceSR(Hubb,epsilon,delta)
                    replaceEV(prime)
                    print("Hubble: " + str(Hubb) + "\t Epsilon: " + str(epsilon) + "\t Delta: " + str(delta) + "\t Seed: " + str(prime))
                    copy("events.ini",dir)
                    copy(radiusFile,dir)
                    copy("fomodel/SR.ini",dir + "/fomodel")
                    copy("share/particles.data",dir + "/share")
                    copy("share/decays.data",dir + "/share")
                    file.write(f"{unigendir} {scriptdir} {dir}\n")
                except Exception as e:
                    print(e)
                
file.close()
copy(jobarrayfile, submissiondir)
copy("jobScript_SL.sh", submissiondir)
copy("wrap.sh",submissiondir)
chmod_command = "chmod 777 " + submissiondir + "/jobScript_SL.sh"
print(chmod_command)
os.system(chmod_command)

arrElem = NoFiles * len(HubbList) * len(epsilonList) * len(deltalist)

#command="sbatch --array=1-{3} --mem=2000 --time=0-48:00:00 --partition=long -D {0}  --output={1}/out/slurm-%A_%a.out -- {0}/wrap.sh {0}/jobScript_SL.sh {0}/jobarray.dat {1} {2}".format(submissiondir,outputdir,scriptdir,arrElem+1)
command="sbatch --array=1-{3} --mem=2000 --time=0-48:00:00 --partition=long -D {0}  --output={1}/out/slurm-%A_%a.out -- {0}/jobScript_SL.sh {0}/jobarray.dat {1} {2}".format(submissiondir,outputdir,scriptdir,arrElem+1)
print(command)
os.system(command)

