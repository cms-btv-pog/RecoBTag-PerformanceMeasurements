import optparse
import os,sys
import pickle
import math
from rounding import *
import getpass,socket

def ShowSFbSummary(inF,outF):

    #read results from the pickle file
    cachefile = open(inF,'r')
    fitInfo=pickle.load(cachefile)
    fitResults=pickle.load(cachefile)
    mcTruth=pickle.load(cachefile)
    cachefile.close()

    taggerDef=fitInfo['taggerDef']

    #init output file
    out = open(outF,'w')
    out.write('\\documentclass[10pt,a4paper]{report}\n')
    out.write('\\usepackage{graphicx}\n')
    out.write('\\usepackage{color}\n')
    out.write('\\begin{document}\n')
    out.write('\\title{Efficiency measurements for %s using $t\\bar{t}$ events}\n'%taggerDef[0])
    out.write('\\author{%s@%s}\n'%(getpass.getuser(),socket.gethostname()))
    out.write('\\date{\\today}\n')
    out.write('\\maketitle\n')

    #iterate over the results
    for islice in fitResults[0]:

        tableHeader=['Op. point']
        tableTruth=['$\\varepsilon_{\\rm b}^{\\rm MC}$']
        tableObs=['$\\varepsilon_{\\rm b}^{\\rm obs}$']
        tableSF=['${\\rm SF}_{\\rm b}$']
        tableSyst={}

        nPreTagObs    = fitResults[0][islice][''][0]
        nPreTagObsUnc = fitResults[0][islice][''][1]
        nPreTagExp    = mcTruth[0][islice][0]
        nPreTagExpUnc = mcTruth[0][islice][1]
        
        for iop in xrange(1,len(taggerDef)-2):

            tableHeader.append('$>$%3.3f'%taggerDef[iop+1])

            nTagExp    = mcTruth[iop][islice][0]
            nTagExpUnc = mcTruth[iop][islice][1]
            bEffExp=nTagExp/nPreTagExp
            bEffExpUnc=math.sqrt((nTagExp*nPreTagExpUnc)**2+(nTagExpUnc*nPreTagExp)**2)/(nPreTagExp**2)
            tableTruth.append('%s'%toLatexRounded(bEffExp,bEffExpUnc))

            nTagObs    = fitResults[iop][islice][''][0]
            nTagObsUnc = fitResults[iop][islice][''][1]
            bEffObs=nTagObs/nPreTagObs
            bEffObsUnc=math.sqrt((nTagObs*nPreTagObsUnc)**2+(nTagObsUnc*nPreTagObs)**2)/(nPreTagObs**2)
            tableObs.append('%s'%toLatexRounded(bEffObs,bEffObsUnc))

            sfb          = bEffObs/bEffExp
            sfbStatUnc   = (bEffObsUnc*bEffExp)/(bEffExp**2)

            sfbMCStatUnc = (bEffObs*bEffExpUnc)/(bEffExp**2)
            sfbSystUnc   =  sfbMCStatUnc**2
            for var in fitResults[iop][islice]:
                if len(var)==0 : continue
                if var.endswith('dn') : continue
                syst=var[:-2]
                nTagObsUp    = fitResults[iop][islice][syst+'up'][0]
                #sfbObsUp     = ((nTagObsUp-nTagObs)/nPreTagObs)/bEffExp
                sfbObsUp     = ((nTagObsUp-nTagExp)/nPreTagExp)/bEffExp

                nTagObsDn    = fitResults[iop][islice][syst+'dn'][0]
                #sfbObsDn     = ((nTagObsDn-nTagObs)/nPreTagObs)/bEffExp
                sfbObsDn     = ((nTagObsDn-nTagExp)/nPreTagExp)/bEffExp

                if not syst in tableSyst: tableSyst[syst]=['~~~{\small \it %s}'%syst]
                tableSyst[syst].append('${\small %.1g / %.1g }$'%(sfbObsUp,sfbObsDn))
                
                sfbSystUnc += (0.5*(math.fabs(sfbObsUp)+math.fabs(sfbObsDn)))**2

            sfbSystUnc=math.sqrt(sfbSystUnc)
            tableSF.append('%s'%(toLatexRounded(sfb,(sfbStatUnc,sfbSystUnc))))


        table='\\begin{table}[h]\n'
        table+='\\caption{Efficiency measurement,'
        table += ' $%3.1f<%s<%3.1f$.}\n'%(fitInfo['slicebins'][islice][0],fitInfo['slicevar'],fitInfo['slicebins'][islice][1])
        table+='\\begin{center}\n'
        table+='\\begin{tabular}[h]{l%s}\n'%('c'*(len(tableHeader)-1))
        for row in [tableHeader,tableTruth,tableObs,tableSF]:
            table+='\\hline\n'
            for icol in xrange(0,len(row)):
                if icol>0 and icol<len(row): table += '&'
                table += row[icol]
            table+='\\\\\n'
        table+='\\hline\n'
        table+='\\multicolumn{%d}{l}{\\it Systematic uncertainties}\\\\\n'%len(tableHeader)
        for syst in tableSyst:
              for icol in xrange(0,len(tableSyst[syst])):
                if icol>0 and icol<len(tableSyst[syst]): table += '&'
                table += tableSyst[syst][icol]
              table+='\\\\\n'
        table+='\\hline'
        table+='\\end{tabular}\n'
        table+='\\end{center}\n'
        table+='\\end{table}\n'
        out.write(table)

    out.write('\\end{document}\n')
    out.close()
    print 'Output report available @ %s'%outF


"""
steer the script
"""
def main():
    
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',              dest='input',              help='input pickle file',   default=None,  type='string')
    parser.add_option('-o', '--output',             dest='output',             help='output TeX file',     default=None,  type='string')
    (opt, args) = parser.parse_args()
    
    if opt.output is None: opt.output=opt.input+'.tex'

    ShowSFbSummary(inF=opt.input,outF=opt.output)

    #all done here
    exit(0)



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())

