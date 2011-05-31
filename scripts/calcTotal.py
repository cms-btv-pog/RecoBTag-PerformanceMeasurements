#! /usr/bin/env python

import sys
import os
import math
import optparse

    
if __name__ =='__main__':

    listoffiles = ['syst_awayjet.txt','syst_closure.txt','syst_gluon.txt','syst_mupt.txt','syst_ptrel.txt']

    total = []
    bins = []
    ibin = 0
    countfiles = 0
    for ifile in listoffiles:

        tmpfile = open(ifile)
        ibin = 0
        
        for line in tmpfile:

            tmplist = line.split()

            bin = tmplist[0]
            error = float(tmplist[1])

            if countfiles==0: bins.append( bin )
            if len(total) == ibin:
                total.append( 0 )
            toterror2 = error*error + total[ibin]*total[ibin]

            total[ibin] = float(toterror2)
            ibin += 1
        countfiles += 1

    for ii in range(0, len(bins)):

        grandtotal = math.sqrt( total[ii] )  

        print bins[ii] + " " + str( grandtotal )

