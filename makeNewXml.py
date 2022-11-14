#!python

fin = open("radius_Np_124.xml","rt")
fout = open("radius_const.xml","wt")

for i,line in enumerate(fin):
    if i > 19 and i < 8483:
        oldVal = float(line)
        fout.write(str(9.86242)+"\n")
    else:
        fout.write(line)

fin.close()
fout.close()