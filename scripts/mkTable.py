#! /usr/bin/env python

import sys
import os
import math
import optparse

    
if __name__ =='__main__':

    ## command line options
    parser = optparse.OptionParser ("Usage: %prog --command file-reference.txt [--commandMC file-mc.txt --command1 file-syst1.txt --command2 file-syst2.txt --systname --outlatex syst.tex]")
    # required parameters
    cmdGroup = optparse.OptionGroup (parser, "Command Options ")
    cmdGroup.add_option ('--reference', dest='command', action='store',
                         help = 'reference file')
    cmdGroup.add_option ('--mc', dest='commandMC', action='store',
                         help = 'reference MC file')
    cmdGroup.add_option ('--syst1', dest='command1', action='store',
                         help = 'systematic file')
    cmdGroup.add_option ('--syst2', dest='command2', action='store',
                         help = 'systematic file')
    cmdGroup.add_option ('--systname', dest='systname', action='store',
                         help = 'name of systematic output file')
    cmdGroup.add_option ('--outlatex', dest='outlatex', action='store',
                         help = 'name of the latex output file')
    
    parser.add_option_group (cmdGroup)
    (options, args) = parser.parse_args()
    #if len (args) < 1:
    #    raise RuntimeError, "must provide a path to a directory where production files will be stored."
    #if not options.command:
    #    raise RuntimeError, "Exactly one command option must be specified"

    #print options.command
    #print options.command1

    IsLongTable = False
    
    listoffiles = []
    listofnames = []
    listofbins = []
    if options.command:
        listoffiles.append(options.command)
        listofnames.append('reference')
    if options.commandMC:
        listoffiles.append(options.commandMC)
        listofnames.append('MC')
    if options.command1:
        listoffiles.append(options.command1)
        listofnames.append('syst1')
        
    if options.command2:
        listoffiles.append(options.command2)
        listofnames.append('syst2')
        IsLongTable = True

    bigtable = {} # per bin pt
    rowtable = {}
    iif = 0
    for ifile in listoffiles:
        
        
        samplename = listofnames[iif]
        iif += 1
        print "reading file "+ifile
        
        txtfile = open(ifile)
        gotkeyword = False
        akeyword = "[Binned]"
        bin_pt = '0'
        for line in txtfile:

            if line.find('Generate Group')!=-1:
                print "done reading file"
                break

            gotMCeff = False
            if line.find(akeyword)!=-1:
                gotkeyword = True

            if gotkeyword:
                #print line
                

                if line.find("-- Bin ")!=-1:

                    bin_mc_beff = bin_mc_befferr = -1
                    bin_s8_beff = bin_s8_befferr = -1
                    
                    line = txtfile.next()
                    line = txtfile.next()
                    line = txtfile.next()
                    #print line
                    tmplist = line.split()
                
                    if len(tmplist)>0:
                        bin_pt = tmplist[0]
                        #print bin_pt

                    getMCeff = True
                if line.find("eff_tag_b")!=-1:
                    tmplist = line.split()
                    #print tmplist
                    if getMCeff:
                        bin_mc_beff = tmplist[2]
                        bin_mc_befferr = tmplist[4]
                        getMCeff=False
                        #print bin_mc_beff
                        #print bin_mc_befferr
                        rowtable['mc_beff']=bin_mc_beff
                        rowtable['mc_beff_err']=bin_mc_befferr
                    else:
                        bin_s8_beff = tmplist[2]
                        bin_s8_befferr = tmplist[4]
                        #print bin_s8_beff
                        #print bin_s8_befferr
                        rowtable['s8_beff']=bin_s8_beff
                        rowtable['s8_beff_err']=bin_s8_befferr
                        
                        bigtable[bin_pt+"_"+samplename] = rowtable
                        if iif==1: listofbins.append(bin_pt+"_"+samplename)
                        print "storing row for "+ bin_pt+"_"+samplename
                        print rowtable
                        #print bigtable[bin_pt+"_"+samplename]
                        rowtable = {}
        gotkeywork = False
        
#print bigtable["35.0000_reference"]

print bigtable.keys()

print "Now print table in file: "+options.outlatex +" \n\n"
thetable = '''
\\begin{table}[h]
\\begin{centering}
'''
if IsLongTable:
    thetable +='''
    \\begin{tabular}{|c|c|c|c|c|c|c|} \hline 
    jet $p_T$ & Nominal Data & low & $\Delta$ & high & $\Delta$ & $\Delta_{max}$  \\\ \hline
    '''
else:
    thetable +='''
    \\begin{tabular}{|c|c|c|c|} \hline
    jet $p_T$ & Nominal Data & systematic & $\Delta$ \\\ \hline"
    '''
    
outputfile = None
if options.systname:
    outputfile = open(options.systname,"w")

for irow in listofbins:

    aline = ''
    sp = ' & '
    bin_pt = irow.split("_")[0]
    bin_pt_str = str( float(bin_pt) )
    aline += bin_pt_str + sp
    niisample = 0
    delta = [0,0,0]
    maxdelta = 0
    #keepline = True
    for isample in listofnames:

        #print "bin is "+bin_pt+" sample is "+isample

        if not bigtable.has_key(bin_pt+"_"+isample):
            #keppline = False
            aline += sp + sp
            continue
        tmprow = bigtable[bin_pt+"_"+isample]
        #print tmprow
        if niisample ==0 :
            delta[0] = float(tmprow['s8_beff'])
            #print str(delta)
            aline += tmprow['s8_beff'] + ' $\pm$ ' + tmprow['s8_beff_err'] +sp
        else:    
            if isample!="MC":
                aline += tmprow['s8_beff'] + ' $\pm$ ' + tmprow['s8_beff_err'] +sp
                #print delta
                delta[niisample] = round(math.fabs(delta[0] - float(tmprow['s8_beff']))/delta[0],3)
            else:
                aline += tmprow['mc_beff'] + ' $\pm$ ' + tmprow['mc_beff_err'] +sp
                delta[niisample] = round(math.fabs(delta[0] - float(tmprow['mc_beff']))/delta[0],3)

            if IsLongTable:
                aline += str(delta[niisample]) + sp
            else:
                aline += str(delta[niisample]) 
        niisample += 1

    if niisample>1 and IsLongTable:
        if delta[1]>=delta[2]: maxdelta = delta[1]
        else: maxdelta = delta[2]
        aline += str(maxdelta) + ' \\\\'
    else:
        aline += ' \\\\'

    #print aline
    thetable += aline +'\n'
    
    aline = ''
    if outputfile:
        outputfile.write(bin_pt+" "+str(maxdelta)+"\n")

thetable += '''
\hline
\end{tabular}
%\caption{Systematic uncertainties}\label{tab:tab}
\end{centering}
\end{table}
'''

outputfile.close()

outtexfile = open(options.outlatex,"w")
outtexfile.write(thetable)
outtexfile.close()


