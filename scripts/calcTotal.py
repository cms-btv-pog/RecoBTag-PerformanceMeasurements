#! /usr/bin/env python

import sys
import os
import math
import optparse

    
if __name__ =='__main__':

    if len(sys.argv) < 2:
        print "\n need to run: calcTotal.py total.txt syst1.txt syst2.txt syst3.txt ...\n"
        sys.exit()
        
    outfilename = sys.argv[1]

    listoffiles = []
    
    for i in range(2,len(sys.argv)):

        listoffiles.append( sys.argv[i] )
        print "adding file "+sys.argv[i]
        
    #listoffiles = ['syst_awayjet.txt','syst_closure.txt','syst_gluon.txt','syst_mupt.txt','syst_ptrel.txt']

    total = []
    bins_min = []
    bins_max = []
    ibin = 0
    countfiles = 0
    for ifile in listoffiles:

        tmpfile = open(ifile)
        ibin = 0
        
        for line in tmpfile:

            tmplist = line.split()

            bin_min = tmplist[0]
            bin_max = tmplist[1]
            error = float(tmplist[2])

            if countfiles==0:
                bins_min.append( bin_min )
                bins_max.append( bin_max )
            if len(total) == ibin:
                total.append( 0 )
            toterror2 = error*error + total[ibin]*total[ibin]

            total[ibin] = float(toterror2)
            ibin += 1
        countfiles += 1

    outfile = open(outfilename,"w")
    
    for ii in range(0, len(bins_min)):

        grandtotal = math.sqrt( total[ii] )  

        #print bins_min[ii] + " " + bins_max[ii] + " " + str( grandtotal )
        outfile.write(bins_min[ii] + " " + bins_max[ii] + " " + str( grandtotal ) +"\n")
    outfile.close()
    
