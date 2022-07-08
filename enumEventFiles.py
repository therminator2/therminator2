#!python

import os
from shutil import copytree, copy

lustrDir = "/u/kjedrzej/lustre/hades/user/kjedrzej/HubbDeltMotornenko"

Harr = [0.0225]#0.015,0.0155,0.016,0.0165,0.017
Earr = [0.0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6]
Darr = [0.4]#0.2,0.3,0.4,0.5,0.6
NoFiles = 1000
append = True

for Hubb in Harr:
    for epsilon in Earr:
        for delta in Darr:
            newH = str(Hubb*1000)
            newH = newH.replace('.','')
            newE = abs(int(epsilon*11)) 
            newD = int(delta*11)

            outDirName = lustrDir + "/H{}E{}D{}femto".format(newH,newE,newD)
            try:
                os.mkdir(outDirName)
            except Exception as e:
                print(e)
            for i,file in enumerate(range(0,NoFiles)):
                inDirName = lustrDir + "/H{}_E{}_D{}_file{}".format(newH,newE,newD,file+1)
                ofn = outDirName + "/event{num:03d}.root".format(num=file)
                if os.path.exists(ofn) and append:
                    print(ofn + "\t exists")
                    continue
                try:
                    if i == 0:
                        copy(inDirName+"/events.ini",outDirName)
                        copy(inDirName+"/sr/SR.ini",outDirName)
                        #copy(inDirName+"/sr/fmult*",outDirName)
                    copy(inDirName + "/sr/event000.root", ofn)
                    print(ofn)
                except Exception as e:
                    print(e)

                
