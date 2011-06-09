#! /usr/bin/env python

import sys
import os
import math
import optparse

    
if __name__ =='__main__':

    if len(sys.argv) < 2:
        print "\n need to run: calcTotal.py total.tex away:syst1.txt ptrel:syst2.txt gluon:syst3.txt ...\n"
        sys.exit()
        
    outfilename = sys.argv[1]

    listoffiles = []
    listofsyst = []
    for i in range(2,len(sys.argv)):
        tmplist = sys.argv[i].split(":")
        listofsyst.append( tmplist[0] )
        listoffiles.append( tmplist[1] )
        print "adding file "+tmplist[1] +" for systematic "+tmplist[0]
        
    #listoffiles = ['syst_awayjet.txt','syst_closure.txt','syst_gluon.txt','syst_mupt.txt','syst_ptrel.txt']

    total = []
    bins_min = []
    bins_max = []
    ibin = 0
    countfiles = 0
    rowtable = {}

    wtable = {}
    
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

            if ibin == 0:
                rowtable[listofsyst[countfiles]] = []
                
            rowtable[listofsyst[countfiles]].append( error)
            
            if bin_min == "50":
                wtable[listofsyst[countfiles]] = error
            if bin_min == "60" or bin_min == "70":
                wtable[listofsyst[countfiles]] += error
            
            toterror2 = error*error + total[ibin]*total[ibin]

            total[ibin] = float(toterror2)
            ibin += 1
        countfiles += 1

    outfile = open(outfilename,"w")

    wline = ""
    wtotal = 0.
    for ii in range(0, len(bins_min)):

        #grandtotal = math.sqrt( total[ii] )  
        grandtotal = 0.
        #print bins_min[ii] + " " + bins_max[ii] + " " + str( grandtotal )
        tmpline = ''
        
        for asyst in listofsyst:
            #print asyst
            avect = rowtable[asyst]
            #print avect
            tmpline += str( avect[ii]*100. ) + " & "
            grandtotal += avect[ii]*avect[ii]
                        
        #print tmpline
        grandtotal = math.sqrt( grandtotal )
        aline = bins_min[ii] + " - " + bins_max[ii] + " & " + tmpline + str( round(grandtotal*100.,1) ) + " \\\\ \n"
        print aline
        #outfile.write(bins_min[ii] + " " + bins_max[ii] + " " + str( grandtotal ) +"\n")
        outfile.write( aline )

        if ii>= 3 and ii<=5:
            
            for asyst in listofsyst:
                tmpline += str(wtable[asyst]/3.)+" & "
                wtotal += wtable[asyst]*wtable[asyst]
            wline = tmpline    
        #if bins_min[ii] == "50" or bins_min[ii]== "60" or bins_min[ii]=="70":
    #print "weighted bin:"
    #print "50 - 80 & " + wline + " & " + str(round( math.sqrt(wtotal)/3.,1)) + " \\\\ "
    
    outfile.close()
