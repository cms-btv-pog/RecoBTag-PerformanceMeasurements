import optparse
import os,sys
import pickle
import math
from rounding import *
import getpass,socket
import ROOT


def initSummaryGrs(taggerName,nop):
    summaryGrs={}
    for sliceType in ['inc','diff']:
        for summaryType in ['expEff','obsEff','sfb']:
            for uncType in ['stat','total']:
                for iop in xrange(0,nop):
                    summaryGrs['%s_%s_%s_%d_%s'%(taggerName,summaryType,uncType,iop,sliceType)]=ROOT.TGraphErrors()
    return summaryGrs

def addToSummaryGr(summaryGrs,val,valUnc,xmin,xmax,taggerName,summaryType,uncType,sliceNb,iop):
    sliceType='inc' if sliceNb==0 else 'diff'
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
        
def ShowSFbSummary(inF,outF):

    #read results from the pickle file
    cachefile = open(inF,'r')
    fitInfo=pickle.load(cachefile)
    fitResults=pickle.load(cachefile)
    mcTruth=pickle.load(cachefile)
    cachefile.close()

    taggerDef=fitInfo['taggerDef']
    taggerName=fitInfo['tagger']

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

    #prepare graphs
    summaryGrs=initSummaryGrs(taggerName=taggerName,nop=len(taggerDef)-2)

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
            addToSummaryGr(summaryGrs,
                           bEffExp,bEffExpUnc,
                           fitInfo['slicebins'][islice][0],fitInfo['slicebins'][islice][1],
                           taggerName,'expEff','stat',islice,iop)

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
            sfbTotalUnc=math.sqrt(sfbStatUnc**2+sfbSystUnc**2)
            
            tableSF.append('%s'%(toLatexRounded(sfb,(sfbStatUnc,sfbSystUnc))))
            addToSummaryGr(summaryGrs,
                           sfb,sfbStatUnc,
                           fitInfo['slicebins'][islice][0],fitInfo['slicebins'][islice][1],
                           taggerName,'sfb','stat',islice,iop)
            addToSummaryGr(summaryGrs,
                           sfb,sfbTotalUnc,
                           fitInfo['slicebins'][islice][0],fitInfo['slicebins'][islice][1],
                           taggerName,'sfb','total',islice,iop)
            addToSummaryGr(summaryGrs,
                           bEffObs,bEffObsUnc,
                           fitInfo['slicebins'][islice][0],fitInfo['slicebins'][islice][1],
                           taggerName,'obsEff','stat',islice,iop)
            addToSummaryGr(summaryGrs,
                           bEffObs,bEffObs*sfbTotalUnc/sfb,
                           fitInfo['slicebins'][islice][0],fitInfo['slicebins'][islice][1],
                           taggerName,'obsEff','total',islice,iop)


        #dump table
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


    outDir=os.path.dirname(outF)
    if len(outDir)==0 : outDir='./'
    showSummary(taggerDef,taggerName,fitInfo['slicevar'],summaryGrsMap=[(fitInfo['var'],summaryGrs),],outDir=outDir)

    #include figures in the tex file
    for iop in xrange(1,len(taggerDef)-2):
        out.write('\\begin{figure}[htp]')
        out.write('\\centering\n')
        out.write('\\includegraphics[width=0.8\\textwidth]{%s/%s_%d.pdf}' % (outDir,taggerDef[1],iop))
        out.write('\\caption{\n')
        out.write('Summary of the results for the measurement of the b-tagging efficiency of %s algorithm $>$%3.2f. \n' % (taggerDef[0],taggerDef[1+iop]) )
        out.write('The top panels show the observed b-tagging efficiency compared to the one expected in simulation. \n')
        out.write('While the left panel shows the inclusive results, the right panel shows the results as function of the %s. \n' % fitInfo['slicevar'] )
        out.write('The bottom panels show the SF$_{\\rm b}$ fit to data. \n')
        out.write('The error bars (hatched areas) represent the statistical (total) uncertainty of the measurement.\n')
        out.write('}\n')
        out.write('\\label{fig:%ssummary%d}' % (taggerDef[1],iop))
        out.write('\\end{figure}\n')

    out.write('\\end{document}\n')
    out.close()
    print 'Output report available @ %s'%outF
   

"""
"""
def showSummary(taggerDef,taggerName,sliceVar,summaryGrsMap,outDir):

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)

    #show results in canvas
    c=ROOT.TCanvas('c','c',800,500)
    c.SetRightMargin(0)
    c.SetLeftMargin(0)
    c.SetBottomMargin(0)
    c.SetTopMargin(0)

    colors=[ROOT.kAzure+9,ROOT.kOrange-1,ROOT.kGreen-5]
    markers=[20,21,22]

    c.cd()
    clabel=ROOT.TPad('clabel','clabel',0,0.95,1.0,1.0)
    clabel.SetRightMargin(0.0)
    clabel.SetTopMargin(0.0)
    clabel.SetLeftMargin(0.0)
    clabel.SetBottomMargin(0.0)
    clabel.Draw()

    subpads={}
    c.cd()
    key=('Eff','inc')
    subpads[key]=ROOT.TPad('cinc','cinc',0,0.45,0.2,0.95)
    subpads[key].SetRightMargin(0.02)
    subpads[key].SetTopMargin(0.01)
    subpads[key].SetLeftMargin(0.3)
    subpads[key].SetBottomMargin(0.01)
    subpads[key].Draw()

    c.cd()
    key=('sfb','inc')
    subpads[key]=ROOT.TPad('cincsf','cincsf',0,0,0.2,0.45)
    subpads[key].SetRightMargin(0.02)
    subpads[key].SetTopMargin(0.01)
    subpads[key].SetLeftMargin(0.3)
    subpads[key].SetBottomMargin(0.12)
    subpads[key].Draw()

    c.cd()
    key=('Eff','diff')
    subpads[key]=ROOT.TPad('cexc','cexc',0.2,0.45,1.0,0.95)
    subpads[key].SetRightMargin(0.02)
    subpads[key].SetTopMargin(0.01)
    subpads[key].SetLeftMargin(0.02)
    subpads[key].SetBottomMargin(0.01)
    subpads[key].Draw()

    c.cd()
    key=('sfb','diff')
    subpads[key]=ROOT.TPad('cexcsf','cexcsf',0.2,0,1.0,0.45)
    subpads[key].SetRightMargin(0.02)
    subpads[key].SetTopMargin(0.01)
    subpads[key].SetLeftMargin(0.02)
    subpads[key].SetBottomMargin(0.12)
    subpads[key].Draw()
    
    for iop in xrange(1,len(taggerDef)-2):

        for key in subpads:
            summaryType,sliceType=key[0],key[1]
            subpads[key].cd()
            subpads[key].Clear()
            ivar=0
            for uncType in ['total','stat']:
                ivar=0
                drawOpt='a2' if uncType=='total' else 'p'
                for ivar in xrange(0,len(summaryGrsMap)):
                    var,summaryGrs = summaryGrsMap[ivar][0], summaryGrsMap[ivar][1]

                    grKey='%s_obs%s_%s_%d_%s'%(taggerName,summaryType,uncType,iop,sliceType)
                    if 'sfb' in summaryType:
                        grKey='%s_%s_%s_%d_%s'%(taggerName,summaryType,uncType,iop,sliceType)
                    gr=summaryGrs[grKey]
                
                    gr.SetTitle(var)
                    gr.SetMarkerColor(colors[ivar]+1)
                    gr.SetLineColor(colors[ivar]+1)
                    if uncType=='total': gr.SetFillColor(colors[ivar]+1)
                    else : gr.SetFillColor(0)
                    gr.SetMarkerStyle(markers[ivar])
                    
                    gr.SetFillStyle(3001)
                    gr.Draw(drawOpt)
                    if 'Eff' in summaryType :
                        gr.GetYaxis().SetTitle('Efficiency')
                        gr.GetYaxis().SetRangeUser(0.12,1.5)
                    else :
                        gr.GetYaxis().SetTitle('Scale factor')
                        gr.GetYaxis().SetRangeUser(0.64,1.26)
                    if 'inc' in sliceType:
                        gr.GetYaxis().SetTitleOffset(1.2)
                        gr.GetYaxis().SetTitleSize(0.1)
                        gr.GetYaxis().SetLabelSize(0.1)
                        gr.GetXaxis().SetTitleSize(0.)
                        gr.GetXaxis().SetLabelSize(0.)
                        gr.GetXaxis().SetNdivisions(0)
                    else:                        
                        gr.GetYaxis().SetTitleSize(0.)
                        gr.GetYaxis().SetLabelSize(0.)
                        if 'Eff' in summaryType:
                            gr.GetXaxis().SetTitleSize(0.)
                            gr.GetXaxis().SetLabelSize(0.)
                        else:
                            gr.GetXaxis().SetTitleSize(0.07)
                            gr.GetXaxis().SetLabelSize(0.07)
                            gr.GetXaxis().SetTitle(sliceVar)
                            gr.GetXaxis().SetTitleOffset(0.5)

                    if uncType=='total' : drawOpt='2'
                    ivar+=1
                    
            if 'Eff' in summaryType:
                gr=summaryGrs['%s_expEff_stat_%d_%s'%(taggerName,iop,sliceType)]
                gr.SetMarkerColor(1)
                gr.SetMarkerStyle(1)
                gr.SetLineColor(1)
                gr.SetLineWidth(2)
                gr.Draw('cX')
        c.Modified()
        c.Update()

        clabel.cd()
        clabel.Clear()
        txt=ROOT.TLatex()
        txt.SetNDC(True)
        txt.SetTextFont(43)
        txt.SetTextSize(16)
        txt.SetTextAlign(12)
        lumi=41.6
        if lumi<100:
            txt.DrawLatex(0.05,0.5,'#bf{CMS} #it{Preliminary} %3.1f pb^{-1} (13 TeV)' % lumi)
        else:
            txt.DrawLatex(0.05,0.5,'#bf{CMS} #it{Preliminary} %3.1f fb^{-1} (13 TeV)' % (lumi/1000.))
        txt.DrawLatex(0.7,0.5,'[#scale[0.7]{%s} > %3.3f]' % (taggerDef[0],taggerDef[iop+1]))

        c.cd()
        c.Modified()
        c.Update()
        for ext in ['png','pdf']: c.SaveAs('%s/%s_%d.%s'% (outDir,taggerDef[1],iop,ext))
        


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

