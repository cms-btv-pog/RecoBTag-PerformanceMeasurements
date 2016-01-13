import optparse
import os,sys
import pickle
import math
from rounding import *
import getpass,socket
import ROOT
from Templated_btagEffFitter import SLICEVARTITLES

LUMI=2444.
SUMMARYCTR=0
COLORS=[1,ROOT.kAzure+9,ROOT.kGreen-5,ROOT.kOrange-1] #,ROOT.kMagenta+2]
MARKERS=[20,24,22,25]


def addToSummaryGr(summaryGrs,val,valUnc,xmin,xmax,taggerName,summaryType,uncType,sliceNb,iop):
    sliceType='diff'
    key='%s_%s_%s_%d_%s'%(taggerName,summaryType,uncType,iop,sliceType)
    xcen=0.5*(xmin+xmax)
    dx=(xcen-xmin) if 'total' in uncType else 0.
    np=summaryGrs[key].GetN()
    if 'exp' in key and 'inc' in key:
        summaryGrs[key].SetPoint(np,xmin,val)
        summaryGrs[key].SetPointError(np,0,valUnc)
        summaryGrs[key].SetPoint(np+1,xmax,val)
        summaryGrs[key].SetPointError(np+1,0,valUnc)
        
    else:
        summaryGrs[key].SetPoint(np,xcen,val)
        summaryGrs[key].SetPointError(np,dx,valUnc)
        
"""
Parse pickle file
"""
def buildSFbSummary(inF,title,outDir):

    global SUMMARYCTR
    SUMMARYCTR+=1

    #read results from the pickle file
    cachefile = open(inF,'r')
    fitInfo=pickle.load(cachefile)
    effExpected=pickle.load(cachefile)
    effObserved=pickle.load(cachefile)
    sfbMeasurement=pickle.load(cachefile)
    systUncs=pickle.load(cachefile)    
    cachefile.close()

    taggerDef=fitInfo['taggerDef']
    taggerName=fitInfo['tagger']
    nOPs=len(taggerDef)-2

    # see https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagCalibration#Example_code_in_Python
    sliceVarName=fitInfo['slicevar']
    btvCalib = ROOT.BTagCalibration(title) if sliceVarName=='jetpt' else None

    #prepare graphs
    summaryGrs={}
    for summaryType in ['eff','sfb']:
        for uncType in ['stat','total']:
            for iop in xrange(1,nOPs):
                if not iop in summaryGrs: summaryGrs[iop]={}
                key=(summaryType,uncType)
                summaryGrs[iop][key]=ROOT.TGraphErrors()
                summaryGrs[iop][key].SetName('%s_pass%d_%s_%s_%d'%(taggerName,iop,summaryType,uncType,SUMMARYCTR))
                summaryGrs[iop][key].SetTitle(title)
                summaryGrs[iop][key].SetMarkerColor(COLORS[(SUMMARYCTR-1)%4])
                summaryGrs[iop][key].SetLineColor(COLORS[(SUMMARYCTR-1)%4])
                if uncType=='total':
                    summaryGrs[iop][key].SetFillColor(COLORS[(SUMMARYCTR-1)%4])
                    summaryGrs[iop][key].SetFillStyle(3001+SUMMARYCTR%4)
                    summaryGrs[iop][key].SetMarkerStyle(1)
                else:
                    summaryGrs[iop][key].SetFillColor(0)
                    summaryGrs[iop][key].SetFillStyle(0)
                    summaryGrs[iop][key].SetLineWidth(2)
                    summaryGrs[iop][key].SetMarkerSize(0.6)
                    summaryGrs[iop][key].SetMarkerStyle(MARKERS[(SUMMARYCTR-1)%4])


    #iterate over the results
    table,plotsToInclude=[],[]
    for islice in effExpected[1]:

        tablePerOp={}
        sliceVarMin, sliceVarMax = fitInfo['slicebins'][islice][0], fitInfo['slicebins'][islice][1]
        sliceVarMean = 0.5*(sliceVarMax+sliceVarMin)
        sliceVarDx   = 0.5*(sliceVarMax-sliceVarMin)
        sliceVarMean = sliceVarMean+((SUMMARYCTR-1)%4-2)*sliceVarDx*0.1
        for iop in xrange(1,len(taggerDef)-2):

            systTable=[]
            effExp,effExpUnc = effExpected[iop][islice]
            effObs,effObsUnc = effObserved[iop][islice][0],effObserved[iop][islice][1]
            sfb,sfbStatUnc   = sfbMeasurement[iop][islice]
            sfbSystUnc       = sfbStatUnc**2
            for syst in systUncs[iop][islice]:
                if len(syst)==0 : continue
                if syst.endswith('dn') : continue
                syst=syst[:-2]
                sfbUncUp = systUncs[iop][islice][syst+'up']
                sfbUncDn = systUncs[iop][islice][syst+'dn']
                sfbSystUnc += (0.5*(math.fabs(sfbUncUp)+math.fabs(sfbUncDn)))**2
                systTable.append( ('~~~{\\small \\it %s}'%syst,'${\small %.1g / %.1g }$'%(sfbUncUp,sfbUncDn)) )
            sfbSystUnc=math.sqrt(sfbSystUnc)
            sfbTotalUnc=math.sqrt(sfbStatUnc**2+sfbSystUnc**2)
            
            #report
            if sliceVarName=='jetpt':
                btvCalibParams = ROOT.BTagEntry.Parameters(iop-1, title, 'central', 0, -2.4, 2.4, sliceVarMin,sliceVarMax,0,1)
                entry = ROOT.BTagEntry(str(sfb),btvCalibParams)
                btvCalib.addEntry(entry)
                btvCalibParams = ROOT.BTagEntry.Parameters(iop-1, title, 'up_total', 0, -2.4, 2.4, sliceVarMin,sliceVarMax,0,1)
                entry = ROOT.BTagEntry(str(sfb+sfbTotalUnc),btvCalibParams)
                btvCalib.addEntry(entry)
                btvCalibParams = ROOT.BTagEntry.Parameters(iop-1, title, 'down_total', 0, -2.4, 2.4, sliceVarMin,sliceVarMax,0,1)
                entry = ROOT.BTagEntry(str(sfb-sfbTotalUnc),btvCalibParams)
                btvCalib.addEntry(entry)
                btvCalibParams = ROOT.BTagEntry.Parameters(iop-1, title, 'up_statistics', 0, -2.4, 2.4, sliceVarMin,sliceVarMax,0,1)
                entry = ROOT.BTagEntry(str(sfb+sfbStatUnc),btvCalibParams)
                btvCalib.addEntry(entry)
                btvCalibParams = ROOT.BTagEntry.Parameters(iop-1, title, 'down_statistics', 0, -2.4, 2.4, sliceVarMin,sliceVarMax,0,1)
                entry = ROOT.BTagEntry(str(sfb-sfbStatUnc),btvCalibParams)
                btvCalib.addEntry(entry)


            #fill table rows
            tablePerOp[iop]=[('$\\varepsilon_{\\rm b}^{\\rm MC}$','%s'%toLatexRounded(effExp,effExpUnc)),
                             ('$\\varepsilon_{\\rm b}^{\\rm obs}$','%s'%toLatexRounded(effObs,effObsUnc)),
                             ('${\\rm SF}_{\\rm b}$','%s'%toLatexRounded(sfb,[sfbStatUnc,sfbSystUnc]))] + systTable

            #add points to graphs
            key=('eff','stat')
            np=summaryGrs[iop][key].GetN()
            summaryGrs[iop][key].SetPoint(np,sliceVarMean,effObs)
            summaryGrs[iop][key].SetPointError(np,sliceVarDx*0.1,effObsUnc)

            key=('eff','total')
            summaryGrs[iop][key].SetPoint(np,sliceVarMean,effObs)
            summaryGrs[iop][key].SetPointError(np,sliceVarDx*0.1,effObs*sfbTotalUnc/sfb)

            key=('sfb','stat')
            summaryGrs[iop][key].SetPoint(np,sliceVarMean,sfb)
            summaryGrs[iop][key].SetPointError(np,sliceVarDx*0.1,sfbStatUnc)

            key=('sfb','total')
            summaryGrs[iop][key].SetPoint(np,sliceVarMean,sfb)
            summaryGrs[iop][key].SetPointError(np,sliceVarDx*0.1,sfbTotalUnc)

        sliceVarRangeText='$%3.0f<%s<%3.0f$'%(sliceVarMin,fitInfo['slicevar'],sliceVarMax)
        table.append( (sliceVarRangeText, tablePerOp) )
        plotsToInclude.append( (
                'Result of the fit to the %s disciminator for events with %s. The left (right) panel shows the pass (fail) category defined for events with %s$>$%s.' 
                % (title,sliceVarRangeText,taggerDef[0],taggerDef[iop+1]),
                '%s/%s_%d_slice%d.pdf'
                % (os.path.dirname(inF),taggerName,iop,islice)
                ) )

    #dump to file
    if btvCalib:
        with open('%s/%s_calib.csv'%(outDir,title), 'w') as f : f.write(btvCalib.makeCSV())

    return sliceVarName,summaryGrs,table,plotsToInclude


"""
"""
def produceSummaryFigures(sliceVar,summaryGr,output):

    plotsToInclude=[]

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)

    #show results in canvas
    c=ROOT.TCanvas('c','c',500,500)
    c.SetRightMargin(0)
    c.SetLeftMargin(0)
    c.SetBottomMargin(0)
    c.SetTopMargin(0)

    c.cd()
    p1=ROOT.TPad('cexc','cexc',0.0,0.45,1.0,1.0)
    p1.SetRightMargin(0.02)
    p1.SetTopMargin(0.01)
    p1.SetLeftMargin(0.12)
    p1.SetBottomMargin(0.01)
    p1.Draw()

    c.cd()
    key=('sfb','diff')
    p2=ROOT.TPad('cexcsf','cexcsf',0.0,0,1.0,0.45)
    p2.SetRightMargin(0.02)
    p2.SetTopMargin(0.01)
    p2.SetLeftMargin(0.12)
    p2.SetBottomMargin(0.2)
    p2.Draw()

    firstKey=summaryGr.keys()[0]
    for iop in summaryGr[firstKey]:

        #efficiency pad
        p1.cd()
        baseDrawOpt='a'
        leg=ROOT.TLegend(0.65,0.8,0.95,0.95)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.05)
        leg.SetHeader("#it{Methods}")

        

        for uncType in ['total','stat']:
            drawOpt='2' if uncType=='total' else 'p'
            
            for key in summaryGr:                   
                summaryGr[key][iop][('eff',uncType)].Draw(baseDrawOpt+drawOpt)
                summaryGr[key][iop][('eff',uncType)].GetYaxis().SetTitle('Efficiency')
                summaryGr[key][iop][('eff',uncType)].GetYaxis().SetRangeUser(0.34,1.08)
                summaryGr[key][iop][('eff',uncType)].GetXaxis().SetTitleSize(0.)
                summaryGr[key][iop][('eff',uncType)].GetXaxis().SetLabelSize(0.0)
                summaryGr[key][iop][('eff',uncType)].GetYaxis().SetTitleOffset(0.8)
                summaryGr[key][iop][('eff',uncType)].GetYaxis().SetTitleSize(0.06)
                summaryGr[key][iop][('eff',uncType)].GetYaxis().SetLabelSize(0.06)
                if uncType=='stat': leg.AddEntry(summaryGr[key][iop][('eff',uncType)],key,'p') 
                baseDrawOpt=''

        leg.SetNColumns(2)
        leg.Draw()
        txt=ROOT.TLatex()
        txt.SetNDC(True)
        txt.SetTextFont(43)
        txt.SetTextSize(14)
        txt.SetTextAlign(12)
        if LUMI<100:
            txt.DrawLatex(0.18,0.95,'#bf{CMS} #it{Preliminary} %3.1f pb^{-1} (13 TeV)' % LUMI)
        else:
            txt.DrawLatex(0.18,0.95,'#bf{CMS} #it{Preliminary} %3.1f fb^{-1} (13 TeV)' % (LUMI/1000.))

        p2.cd()
        baseDrawOpt='a'
        for uncType in ['total','stat']:
            drawOpt='2' if uncType=='total' else 'p'            
            for key in summaryGr:
                summaryGr[key][iop][('sfb',uncType)].Draw(baseDrawOpt+drawOpt)
                summaryGr[key][iop][('sfb',uncType)].GetYaxis().SetTitle('Scale factor')
                summaryGr[key][iop][('sfb',uncType)].GetYaxis().SetRangeUser(0.74,1.16)
                summaryGr[key][iop][('sfb',uncType)].GetYaxis().SetTitleOffset(0.7)
                summaryGr[key][iop][('sfb',uncType)].GetYaxis().SetTitleSize(0.08)
                summaryGr[key][iop][('sfb',uncType)].GetYaxis().SetLabelSize(0.08)
                summaryGr[key][iop][('sfb',uncType)].GetXaxis().SetLabelSize(0.08)
                summaryGr[key][iop][('sfb',uncType)].GetXaxis().SetTitleSize(0.08)
                summaryGr[key][iop][('sfb',uncType)].GetXaxis().SetTitle(SLICEVARTITLES[sliceVar])
                baseDrawOpt=''
        c.cd()
        c.Modified()
        c.Update()
        for ext in ['png','pdf']: c.SaveAs('%s/EfficiencySummary_%s_%d.%s'% (output,sliceVar,iop,ext))
        plotsToInclude.append( ('Efficiency measurements as function of %s for the %d-th working point' % (SLICEVARTITLES[sliceVar],iop),
                                '%s/EfficiencySummary_%s_%d.pdf'% (output,sliceVar,iop) ) )

    return plotsToInclude
                    

"""
"""
def createReport(tableCollection,plotsToInclude,output):

    #create output file
    out = open('sfb_report.tex','w')
    out.write('\\documentclass[10pt,a4paper]{article}\n')
    out.write('\\usepackage{graphicx}\n')
    out.write('\\usepackage{color}\n')
    out.write('\\setcounter{tocdepth}{3}\n')
    out.write('\\begin{document}\n')
    out.write('\\title{b-tagging efficiency measurements using $t\\bar{t}$ dilepton events}\n')
    out.write('\\author{%s@%s}\n'%(getpass.getuser(),socket.gethostname()))
    out.write('\\date{\\today}\n')
    out.write('\\maketitle\n')
    out.write('\\tableofcontents\n')
    out.write('\\section{Summary tables}\n')
    out.write('\\label{sec:tables}\n')
    out.write('The following tables report the scale factors and efficiencies measured by different methods for different slice variables.\n')

    for sliceVar in tableCollection:
        out.write('\\subsection{Results as function of %s}\n'%SLICEVARTITLES[sliceVar])
        firstMethod=tableCollection[sliceVar].keys()[0]

        out.write('\\clearpage\n')
        out.write('\\subsubsection{Tables}\n')
        firstMethod=tableCollection[sliceVar].keys()[0]
        for iop in tableCollection[sliceVar][firstMethod][0][1]:

            for i in xrange(0,len(tableCollection[sliceVar][firstMethod])):
                sliceDef=tableCollection[sliceVar][firstMethod][i][0]
                
                out.write('\\begin{table}[h]\n')
                out.write('\\caption{Efficiency measurement for events with %s working point number %d.}\n'%(sliceDef,iop))
                out.write('\\begin{center}\n')
                out.write('\\begin{tabular}[h]{l%s}\n'%('c'*len(tableCollection[sliceVar])))

                table=['Method &']
                for method in tableCollection[sliceVar]:
                    
                    table[0] += ' %20s &'%method
                    rowCtr=1
                    for row,val in tableCollection[sliceVar][method][i][1][iop]:
                        if len(table)<=rowCtr: table.append( ' %20s &' % row )
                        table[rowCtr] += ' %20s &' % val
                        rowCtr+=1

                #dump table to file
                out.write('\\hline\n')
                for row in table: out.write( row[:-1]+'\\\\\n' )
                out.write('\\hline\n')

                out.write('\\end{tabular}\n')
                out.write('\\end{center}\n')
                out.write('\\end{table}\n')
                out.write('\n\n')

        #dump plots as well
        out.write('\\clearpage\n')
        out.write('\\subsubsection{Figures}\n')        
        firstMethod=tableCollection[sliceVar].keys()[0]
        if sliceVar in plotsToInclude:
            iplot=1
            for caption,plot in plotsToInclude[sliceVar]:
                if iplot%4==0 : out.write('\\clearpage\n')
                out.write('\n')
                out.write('\\begin{figure}[!htbp]')
                out.write('\\centering\n')
                print caption,plot
                out.write('\\includegraphics[width=0.99\\textwidth]{%s}\n' % plot)
                out.write('\\caption{%s}\n'%caption)
                out.write('\\end{figure}\n')
                out.write('\n')
                iplot+=1
    out.write('\\end{document}\n')
    
    #close TeX file
    out.close()
    
    #compile and move to output
    os.system('pdflatex sfb_report.tex')
    os.system('mv -v sfb_report.* %s'%output)
    



"""
steer the script
"""
def main():
    
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',              dest='input',              help='input pickle files (csv list)',   default=None,  type='string')
    parser.add_option('-o', '--output',             dest='output',             help='output directory',                default=None,  type='string')
    (opt, args) = parser.parse_args()

    #prepare output
    os.system('mkdir -p %s'%opt.output)

    os.system('ln -s $CMSSW_RELEASE_BASE/src/RecoBTag/PerformanceDB/test/BTagCalibrationStandalone.cc')
    os.system('ln -s $CMSSW_RELEASE_BASE/src/RecoBTag/PerformanceDB/test/BTagCalibrationStandalone.h')
    ROOT.gROOT.ProcessLine('.L BTagCalibrationStandalone.cc+')

    allInputs=opt.input.split(',')
    print allInputs
    summaryTable,summaryGr,plotsToInclude={},{},{}
    for line in allInputs:
        title,inF=line.split(':')        
        sliceVar, effGrs, tables,plots = buildSFbSummary(inF,title,opt.output)
        if not sliceVar in summaryGr: 
            summaryGr[sliceVar]={}
            summaryTable[sliceVar]={}
            plotsToInclude[sliceVar]=[]
             
        summaryGr[sliceVar][title]=effGrs
        summaryTable[sliceVar][title]=tables
        plotsToInclude[sliceVar]+=plots

    for sliceVar in summaryGr:
        plotsToInclude[sliceVar] += produceSummaryFigures(sliceVar,summaryGr[sliceVar],opt.output)

    createReport(summaryTable,plotsToInclude,opt.output)

    #ShowSFbSummary(inF=opt.input,outF=opt.output)



    #all done here
    exit(0)



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())

