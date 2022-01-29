#!/usr/bin/env python3
import os

rhoTA=[6.66]
rhoBA=[7.00]
rhoSA=[6.50,6.51]
counter=0
for rhoS in rhoSA:
    rhoSS=str(rhoS)
    for rhoB in rhoBA:
        rhoBS=str(rhoB)
        for rhoT in rhoTA:
            rhoTS=str(rhoT)
            counter+=1
            JobN=str(counter)
            #We locate the values and replace them
            #input file
            fin = open("Template/Template.f90", "rt")
            #Make directory for storage
            Dir="Job"+JobN+"_RhoT"+rhoTS+"_RhoB"+rhoBS+"_RhoS"+rhoSS
            os.mkdir(Dir)
            #output file to write the result to
            fout = open(Dir+"/fulltime_2009_05_15_Fortran.f90", "wt")
            #for each line in the input file
            for line in fin:
                #read replace the string and write to output file
                line = line.replace('value1', rhoTS)
                line = line.replace('value2', rhoBS)
                line = line.replace('value3', rhoSS)
                fout.write(line)
            #close input and output files
            fin.close()
            fout.close()
print(counter)
