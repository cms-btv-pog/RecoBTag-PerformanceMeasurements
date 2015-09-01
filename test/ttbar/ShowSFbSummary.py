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
    dx=xcen-xmin if 'total' in uncType else 0.
    np=summaryGrs[key].GetN()
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
        
        if islice==0 : print 'pretag: ',nPreTagExp,nPreTagObs

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

            if islice==0 : print iop,nTagExp,nTagObs
           
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
                           bEffObs,bEffObsUnc*sfbTotalUnc/sfb,
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

    showSummary(taggerDef,taggerName,fitInfo['slicevar'],summaryGrsMap=[(fitInfo['var'],summaryGrs),])

    out.write('\\end{document}\n')
    out.close()
    print 'Output report available @ %s'%outF
   

"""
"""
def showSummary(taggerDef,taggerName,sliceVar,summaryGrsMap):

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
    
    c.cd()
    cinc=ROOT.TPad('cinc','cinc',0,0.45,0.2,0.95)
    cinc.SetRightMargin(0.02)
    cinc.SetTopMargin(0.01)
    cinc.SetLeftMargin(0.3)
    cinc.SetBottomMargin(0.01)
    cinc.Draw()

    c.cd()
    cincsf=ROOT.TPad('cincsf','cincsf',0,0,0.2,0.45)
    cincsf.SetRightMargin(0.02)
    cincsf.SetTopMargin(0.01)
    cincsf.SetLeftMargin(0.3)
    cincsf.SetBottomMargin(0.12)
    cincsf.Draw()

    c.cd()
    cexc=ROOT.TPad('cexc','cexc',0.2,0.45,1.0,0.95)
    cexc.SetRightMargin(0.02)
    cexc.SetTopMargin(0.01)
    cexc.SetLeftMargin(0.02)
    cexc.SetBottomMargin(0.01)
    cexc.Draw()

    c.cd()
    cexcsf=ROOT.TPad('cexcsf','cexcsf',0.2,0,1.0,0.45)
    cexcsf.SetRightMargin(0.02)
    cexcsf.SetTopMargin(0.01)
    cexcsf.SetLeftMargin(0.02)
    cexcsf.SetBottomMargin(0.12)
    cexcsf.Draw()
    
    for iop in xrange(0,len(taggerDef)-2):

        #
        #inclusive efficiency measurements
        #
        cinc.cd()
        cinc.Clear()
        ivar=0

        for uncType in ['total','stat']:
            ivar=0
            drawOpt='a2' if uncType=='total' else 'p'
            for ivar in xrange(0,len(summaryGrsMap)):
                var,summaryGrs = summaryGrsMap[ivar][0], summaryGrsMap[ivar][1]
                gr=summaryGrs['%s_obsEff_%s_%d_inc'%(taggerName,uncType,iop)]
                
                gr.SetTitle(var)
                gr.SetMarkerColor(colors[ivar]+1)
                gr.SetLineColor(colors[ivar]+1)
                if uncType=='total': gr.SetFillColor(colors[ivar]+1)
                else : gr.SetFillColor(0)
                gr.SetMarkerStyle(markers[ivar])
                gr.SetMarkerStyle(1)
                gr.SetFillStyle(3001)
                gr.Draw(drawOpt)
                gr.GetYaxis().SetRangeUser(0.12,1.5)
                gr.GetYaxis().SetTitle('Efficiency')
                gr.GetYaxis().SetTitleOffset(1.2)
                gr.GetYaxis().SetTitleSize(0.1)
                gr.GetYaxis().SetLabelSize(0.1)
                gr.GetXaxis().SetTitleSize(0.)
                gr.GetXaxis().SetLabelSize(0.)
                gr.GetXaxis().SetNdivisions(0)
                if uncType=='total' : drawOpt='2'
                ivar+=1
        gr=summaryGrs['%s_expEff_stat_%d_inc'%(taggerName,iop)]
        gr.SetMarkerColor(1)
        gr.SetMarkerStyle(1)
        gr.SetLineColor(1)
        gr.SetLineWidth(2)
        gr.Draw('c')
        c.Modified()
        c.Update()


        #
        #inclusive scale factors
        #
#        grIncSFObsStatColl=convertToScaleFactor(obs=grIncObsStatColl,ref=grIncExpColl)
#        grIncSFObsTotalColl=convertToScaleFactor(obs=grIncObsTotalColl,ref=grIncExpColl)
#        cincsf.cd()
#        cincsf.Clear()
#        ivar=0
#        drawOpt='a2'
#        for var in grIncSFObsTotalColl:
#            grIncSFObsTotalColl[var][iop].Draw(drawOpt)
#            grIncSFObsTotalColl[var][iop].GetYaxis().SetRangeUser(0.74,1.26)
#            grIncSFObsTotalColl[var][iop].GetYaxis().SetTitle('Data-MC scale factor')
#            grIncSFObsTotalColl[var][iop].GetYaxis().SetTitleOffset(1.2)
#            grIncSFObsTotalColl[var][iop].GetYaxis().SetTitleSize(0.1)
#            grIncSFObsTotalColl[var][iop].GetYaxis().SetLabelSize(0.1)
#            grIncSFObsTotalColl[var][iop].GetXaxis().SetTitleSize(0.)
#            grIncSFObsTotalColl[var][iop].GetXaxis().SetLabelSize(0.)
#            grIncSFObsTotalColl[var][iop].GetXaxis().SetNdivisions(0)
#            drawOpt='2'
#            ivar+=1
#
#        ivar=0
#        for var in grIncSFObsStatColl:
#            grIncSFObsStatColl[var][iop].Draw('p')
#            grIncSFObsStatColl[var][iop].SetTitle(var)
#            ivar+=1
#
#        #
#        # differential efficiency measurements
#        #
#        cexc.cd()
#        cexc.Clear()
#        leg = ROOT.TLegend(0.6, 0.95,0.95,0.8)
#        leg.SetBorderSize(0)
#        leg.SetFillStyle(0)
#        leg.SetTextFont(43)
#        leg.SetTextSize(16)
#        if len(fileList)>2 : leg.SetNColumns(2)
#        leg.AddEntry(grExpColl[iop],grExpColl[iop].GetTitle(),"l")
#
#        ivar=0
#        drawOpt='a2'
#        for var in grObsTotalColl:
#            grObsTotalColl[var][iop].SetTitle(var)
#            grObsTotalColl[var][iop].SetMarkerColor(colors[ivar]+1)
#            grObsTotalColl[var][iop].SetLineColor(colors[ivar]+1)
#            grObsTotalColl[var][iop].SetFillColor(colors[ivar]+1)
#            grObsTotalColl[var][iop].SetMarkerStyle(1)
#            grObsTotalColl[var][iop].SetFillStyle(3001)
#            grObsTotalColl[var][iop].Draw(drawOpt)
#            grObsTotalColl[var][iop].GetYaxis().SetTitleSize(0.)
#            grObsTotalColl[var][iop].GetYaxis().SetLabelSize(0.)
#            grObsTotalColl[var][iop].GetYaxis().SetRangeUser(0.0,1.5)
#            drawOpt='2'
#            ivar+=1
#
#        ivar=0
#        for var in grObsStatColl:
#            grObsStatColl[var][iop].SetMarkerColor(colors[ivar])
#            grObsStatColl[var][iop].SetLineColor(colors[ivar])
#            grObsStatColl[var][iop].SetFillColor(0)
#            grObsStatColl[var][iop].SetMarkerStyle(markers[ivar])
#            grObsStatColl[var][iop].SetFillStyle(0)
#            grObsStatColl[var][iop].Draw('p')
#            grObsStatColl[var][iop].SetTitle(var)
#            leg.AddEntry(grObsStatColl[var][iop],grObsStatColl[var][iop].GetTitle(),'p')
#            ivar+=1
#
#        grExpColl[iop].SetMarkerColor(1)
#        grExpColl[iop].SetMarkerStyle(1)
#        grExpColl[iop].SetLineColor(1)
#        grExpColl[iop].SetLineWidth(2)
#        grExpColl[iop].Draw('c')
#
#        leg.Draw()
#
#
#        #
#        # differential scale factor measurements
#        #
#        grSFObsStatColl=convertToScaleFactor(obs=grObsStatColl,ref=grExpColl)
#        grSFObsTotalColl=convertToScaleFactor(obs=grObsTotalColl,ref=grExpColl)
#        cexcsf.cd()
#        cexcsf.Clear()
#        ivar=0
#        drawOpt='a2'
#        for var in grSFObsTotalColl:
#            grSFObsTotalColl[var][iop].Draw(drawOpt)
#            grSFObsTotalColl[var][iop].GetYaxis().SetTitleSize(0.)
#            grSFObsTotalColl[var][iop].GetYaxis().SetLabelSize(0.)
#            grSFObsTotalColl[var][iop].GetYaxis().SetRangeUser(0.74,1.26)
#            grSFObsTotalColl[var][iop].GetXaxis().SetTitle(SLICEVAR)
#            grSFObsTotalColl[var][iop].GetXaxis().SetTitleSize(0.08)
#            grSFObsTotalColl[var][iop].GetXaxis().SetLabelSize(0.08)
#            drawOpt='2'
#            ivar+=1
#
#        ivar=0
#        for var in grSFObsStatColl:
#            grSFObsStatColl[var][iop].Draw('p')
#            grSFObsStatColl[var][iop].SetTitle(var)
#            ivar+=1
#
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
        txt.DrawLatex(0.6,0.5,'[#scale[0.7]{%s} > %3.3f]' % (taggerDef[0],taggerDef[iop+2]))

        c.cd()
        c.Modified()
        c.Update()
        #for ext in ['png','pdf']: c.SaveAs('%s/%s_%d.%s'% (outDir,tagger,iop,ext))
        raw_input()


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

