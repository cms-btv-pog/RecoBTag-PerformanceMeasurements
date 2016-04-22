import optparse
import os,sys
import json
import ROOT
import math

from rounding import *

SYSTCOLORS=[ROOT.kMagenta, ROOT.kRed+1, ROOT.kMagenta+2, ROOT.kAzure+7, ROOT.kMagenta-9, ROOT.kBlue-7]

"""
A wrapper to store data and MC histograms for comparison
"""
class Plot(object):

    def __init__(self,name):
        self.name = name
        self.mc = {}
        self.mcsyst = {}
        self.dataH = None
        self.data = None
        self._garbageList = []
        self.plotformats = ['pdf','png']
        self.savelog = False
        self.ratiorange = (0.46,1.54)

    def add(self, h, title, color, isData,isSyst):
        h.SetTitle(title)
        if isData:
            try:
                self.dataH.Add(h)
            except:
                self.dataH=h
                self.dataH.SetDirectory(0)
                self.dataH.SetMarkerStyle(20)
                self.dataH.SetMarkerSize(1.4)
                self.dataH.SetMarkerColor(color)
                self.dataH.SetLineColor(ROOT.kBlack)
                self.dataH.SetLineWidth(2)
                self.dataH.SetFillColor(0)
                self.dataH.SetFillStyle(0)
                self._garbageList.append(h)
        else:
            if isSyst:
                try:
                    self.mcsyst[title].Add(h)
                except:
                    self.mcsyst[title]=h
                    self.mcsyst[title].SetName('%s_%s' % (self.mcsyst[title].GetName(), title ) )
                    self.mcsyst[title].SetDirectory(0)
                    self.mcsyst[title].SetMarkerStyle(1)
                    self.mcsyst[title].SetMarkerColor(color)
                    self.mcsyst[title].SetLineColor(ROOT.kBlack)
                    self.mcsyst[title].SetLineWidth(1)
                    self.mcsyst[title].SetFillColor(color)
                    self.mcsyst[title].SetFillStyle(1001)
                    self._garbageList.append(h)
            else:
                try:
                    self.mc[title].Add(h)
                except:
                    self.mc[title]=h
                    self.mc[title].SetName('%s_%s' % (self.mc[title].GetName(), title ) )
                    self.mc[title].SetDirectory(0)
                    self.mc[title].SetMarkerStyle(1)
                    self.mc[title].SetMarkerColor(color)
                    self.mc[title].SetLineColor(ROOT.kBlack)
                    self.mc[title].SetLineWidth(1)
                    self.mc[title].SetFillColor(color)
                    self.mc[title].SetFillStyle(1001)
                    self._garbageList.append(h)

    def finalize(self):
        self.data = convertToPoissonErrorGr(self.dataH)

    def appendTo(self,outUrl):
        outF = ROOT.TFile.Open(outUrl,'UPDATE')
        if not outF.cd(self.name):
            outDir = outF.mkdir(self.name)
            outDir.cd()
        for m in self.mc :
            self.mc[m].Write(self.mc[m].GetName(), ROOT.TObject.kOverwrite)
        if self.dataH :
            self.dataH.Write(self.dataH.GetName(), ROOT.TObject.kOverwrite)
        if self.data :
            self.data.Write(self.data.GetName(), ROOT.TObject.kOverwrite)
        outF.Close()

    def reset(self):
        for o in self._garbageList:
            try:
                o.Delete()
            except:
                pass

    def show(self, outDir,lumi,noScale=False,saveTeX=False):

        if len(self.mc)==0:
            print '%s has no MC!' % self.name
            return

        if self.mc.values()[0].InheritsFrom('TH2') :
            print 'Skipping TH2'
            return

        c = ROOT.TCanvas('c','c',500,500)
        c.SetBottomMargin(0.0)
        c.SetLeftMargin(0.0)
        c.SetTopMargin(0)
        c.SetRightMargin(0.00)


        #holds the main plot
        c.cd()
        p1 = ROOT.TPad('p1','p1',0.0,0.2,1.0,1.0)
        p1.Draw()
        p1.SetRightMargin(0.02)
        p1.SetLeftMargin(0.15)
        p1.SetTopMargin(0.07)
        p1.SetBottomMargin(0.01)
        #p1.SetGridx(True)
        self._garbageList.append(p1)
        p1.cd()

        # legend
        leg = ROOT.TLegend(0.5, 0.85-0.03*max(len(self.mc)-2,0), 0.98, 0.9)        
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(43)
        leg.SetTextSize(16)
        nlegCols = 0

        if self.dataH is not None:
            if self.data is None: self.finalize()
            leg.AddEntry( self.data, self.data.GetTitle(),'p')
            nlegCols += 1
        for h in self.mc:
            if not noScale : self.mc[h].Scale(lumi)
            leg.AddEntry(self.mc[h], self.mc[h].GetTitle(), 'f')
            nlegCols += 1
        if nlegCols ==0 :
            print '%s is empty'%self.name
            return
        leg.SetNColumns(ROOT.TMath.Min(nlegCols/2,3))

        # Build the stack to plot from all backgrounds
        totalMC,nominalTTbar = None,None
        stack = ROOT.THStack('mc','mc')

        #sorted(self.mc, key=lambda histo: histo.Integral())
        from collections import OrderedDict
        self.mc = OrderedDict(sorted(self.mc.items(), key=lambda t: t[1].Integral()))

        for h in self.mc:
            stack.Add(self.mc[h],'hist')            
            try:
                totalMC.Add(self.mc[h])
            except:
                totalMC = self.mc[h].Clone('totalmc')
                self._garbageList.append(totalMC)
                totalMC.SetDirectory(0)
            if h=='t#bar{t}': 
                nominalTTbar= self.mc[h].Clone('nomttbar')
                self._garbageList.append(nominalTTbar)
                nominalTTbar.SetDirectory(0)

        systVariations=[]
        for hname in self.mcsyst:
            systvarH=totalMC.Clone('syst%d'%len(systVariations))
            self._garbageList.append(systvarH)
            systvarH.SetDirectory(0)
            systvarH.Add(nominalTTbar,-1)
            self.mcsyst[hname].Scale(lumi)
            systvarH.Add(self.mcsyst[hname])
            htitle=hname.replace('t#bar{t} ','')
            systvarH.SetTitle(htitle)
            systVariations.append( systvarH )
        totalMC.SetTitle('default')
        systVariations.append(totalMC)

        frame = totalMC.Clone('frame') if totalMC is not None else self.dataH.Clone('frame')
        frame.Reset('ICE')
        if totalMC:
            maxY = totalMC.GetMaximum() 
        if self.dataH:
            if maxY<self.dataH.GetMaximum():
                maxY=self.dataH.GetMaximum()
        frame.GetYaxis().SetRangeUser(0.1,maxY*1.3)
        frame.SetDirectory(0)
        frame.Reset('ICE')
        frame.GetYaxis().SetTitle( '%s / (%.g)' % (frame.GetYaxis().GetTitle(),frame.GetXaxis().GetBinWidth(1)) )
        self._garbageList.append(frame)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetLabelSize(0.045)
        frame.GetXaxis().SetTitleSize(0.05)
        frame.GetXaxis().SetLabelSize(0.045)
        frame.GetYaxis().SetNoExponent()
        frame.Draw()
        frame.GetYaxis().SetTitleOffset(1.3)

        if totalMC is not None   : stack.Draw('hist same')
        if self.data is not None : self.data.Draw('p')

        leg.Draw()
        txt=ROOT.TLatex()
        txt.SetNDC(True)
        txt.SetTextFont(43)
        txt.SetTextSize(16)
        txt.SetTextAlign(12)
        if lumi<100:
            txt.DrawLatex(0.18,0.95,'#bf{CMS} #it{preliminary} %3.1f pb^{-1} (13 TeV)' % (lumi) )
        else:
            txt.DrawLatex(0.2,0.88,'#bf{CMS}')
            txt.DrawLatex(0.2,0.83,'#it{Preliminary}')
            txt.DrawLatex(0.67,0.97,'%3.1f fb^{-1} (13 TeV, 25ns)'%(lumi/1000.))
            #txt.DrawLatex(0.18,0.95,'#bf{CMS} #it{preliminary} %3.1f fb^{-1} (13 TeV)' % (lumi/1000.) )

        #holds the ratio
        c.cd()
        p2 = ROOT.TPad('p2','p2',0.0,0.0,1.0,0.2)
        p2.Draw()
        p2.SetRightMargin(0.02)
        p2.SetLeftMargin(0.15)
        p2.SetTopMargin(0.01)
        p2.SetBottomMargin(0.4)

        #p2.SetGridx(True)
        #p2.SetGridy(True)
        self._garbageList.append(p2)
        p2.cd()
        ratioframe=frame.Clone('ratioframe')
        ratioframe.GetYaxis().SetTitle('Data/MC')
        ratioframe.GetYaxis().SetRangeUser(self.ratiorange[0], self.ratiorange[1])
        self._garbageList.append(frame)
        ratioframe.GetYaxis().SetNdivisions(5)
        ratioframe.GetYaxis().SetLabelSize(0.18)        
        ratioframe.GetYaxis().SetTitleSize(0.2)
        ratioframe.GetYaxis().SetTitleOffset(0.3)
        ratioframe.GetXaxis().SetLabelSize(0.18)
        ratioframe.GetXaxis().SetTitleSize(0.2)
        ratioframe.GetXaxis().SetTitleOffset(0.8)
        ratioframe.Draw()

        leg2=None
        if len(systVariations)>1:
            leg2 = ROOT.TLegend(0.15,0.85,0.95,0.75)
            leg2.SetBorderSize(0)
            leg2.SetFillStyle(1001)
            leg2.SetFillColor(0)
            leg2.SetTextFont(43)
            leg2.SetTextSize(10)

        allGrs=[]
        try:            
            for igr in xrange(0,len(systVariations)):
                ratio=self.dataH.Clone('ratio')
                ratio.SetDirectory(0)
                self._garbageList.append(ratio)
                ratio.Divide(systVariations[igr])
                allGrs.append( ROOT.TGraphAsymmErrors(ratio) )
                allGrs[-1].SetTitle(systVariations[igr].GetTitle())
                allGrs[-1].SetName('ratio%d'%igr)
                if igr==len(systVariations)-1:
                    allGrs[-1].SetMarkerStyle(self.data.GetMarkerStyle())
                    allGrs[-1].SetMarkerSize(self.data.GetMarkerSize())
                    allGrs[-1].SetMarkerColor(self.data.GetMarkerColor())
                    allGrs[-1].SetLineColor(self.data.GetLineColor())
                    allGrs[-1].SetLineWidth(self.data.GetLineWidth())
                    allGrs[-1].Draw('p')
                    if leg2 : leg2.AddEntry(allGrs[-1],allGrs[-1].GetTitle(),'p')
                else:
                    allGrs[-1].SetMarkerStyle(1)
                    allGrs[-1].SetMarkerSize(1)
                    allGrs[-1].SetMarkerColor(SYSTCOLORS[igr])
                    allGrs[-1].SetLineColor(SYSTCOLORS[igr])
                    allGrs[-1].SetLineStyle(1+igr%2)
                    allGrs[-1].SetLineWidth(2)
                    allGrs[-1].Draw('lX')
                    if leg2 : leg2.AddEntry(allGrs[-1],allGrs[-1].GetTitle(),'l')
        except:
            pass

        if leg2:
            leg2.SetNColumns(len(allGrs))
            leg2.Draw()

        #all done
        c.cd()
        c.Modified()
        c.Update()

        #save
        for ext in self.plotformats : c.SaveAs(os.path.join(outDir, self.name+'.'+ext))
        if self.savelog:
            p1.cd()
            p1.SetLogy()
            c.cd()
            c.Modified()
            c.Update()
            for ext in self.plotformats : c.SaveAs(os.path.join(outDir, self.name+'_log.'+ext))

        if saveTeX : self.convertToTeX(outDir=outDir)


    def convertToTeX(self, outDir):
        if len(self.mc)==0:
            print '%s is empty' % self.name
            return

        f = open(outDir+'/'+self.name+'.dat','w')
        f.write('------------------------------------------\n')
        f.write("Process".ljust(20),)
        f.write("Events after each cut\n")
        f.write('------------------------------------------\n')

        tot ={}
        err = {}
        f.write(' '.ljust(20),)
        try:
            for xbin in xrange(1,self.mc.values()[0].GetXaxis().GetNbins()+1):
                pcut=self.mc.values()[0].GetXaxis().GetBinLabel(xbin)
                f.write(pcut.ljust(40),)
                tot[xbin]=0
                err[xbin]=0
        except:
            pass
        f.write('\n')
        f.write('------------------------------------------\n')

        for pname in self.mc:
            h = self.mc[pname]
            f.write(pname.ljust(20),)

            for xbin in xrange(1,h.GetXaxis().GetNbins()+1):
                itot=h.GetBinContent(xbin)
                ierr=h.GetBinError(xbin)
                pval=' & %s'%toLatexRounded(itot,ierr)
                f.write(pval.ljust(40),)
                tot[xbin] = tot[xbin]+itot
                err[xbin] = err[xbin]+ierr*ierr
            f.write('\n')

        f.write('------------------------------------------\n')
        f.write('Total'.ljust(20),)
        for xbin in tot:
            pval=' & %s'%toLatexRounded(tot[xbin],math.sqrt(err[xbin]))
            f.write(pval.ljust(40),)
        f.write('\n')

        if self.dataH is None: return
        f.write('------------------------------------------\n')
        f.write('Data'.ljust(20),)
        for xbin in xrange(1,self.dataH.GetXaxis().GetNbins()+1):
            itot=self.dataH.GetBinContent(xbin)
            pval=' & %d'%itot
            f.write(pval.ljust(40))
        f.write('\n')
        f.write('------------------------------------------\n')
        f.close()



"""
converts a histogram to a graph with Poisson error bars
"""
def convertToPoissonErrorGr(h):

    htype=h.ClassName()
    if htype.find('TH1')<0 : return None

    #check https://twiki.cern.ch/twiki/bin/view/CMS/PoissonErrorBars
    alpha = 1 - 0.6827;
    grpois = ROOT.TGraphAsymmErrors(h);
    for i in xrange(0,grpois.GetN()+1) :
        N = grpois.GetY()[i]
        if N<200 :
            L = 0
            if N>0 : L = ROOT.Math.gamma_quantile(alpha/2,N,1.)
            U = ROOT.Math.gamma_quantile_c(alpha/2,N+1,1)
            grpois.SetPointEYlow(i, N-L)
            grpois.SetPointEYhigh(i, U-N)
        else:
            grpois.SetPointEYlow(i, math.sqrt(N))
            grpois.SetPointEYhigh(i,math.sqrt(N))
    return grpois


"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-j', '--json',        dest='json'  ,      help='json with list of files',        default=None,    type='string')
    parser.add_option(      '--systJson',    dest='systJson',    help='json with list of systematics',  default=None,    type='string')
    parser.add_option('-i', '--inDir',       dest='inDir' ,      help='input directory',                default=None,    type='string')
    parser.add_option(      '--saveLog',     dest='saveLog' ,    help='save log versions of the plots', default=False,   action='store_true')
    parser.add_option(      '--silent',      dest='silent' ,     help='only dump to ROOT file',         default=False,   action='store_true')
    parser.add_option(      '--saveTeX',     dest='saveTeX' ,    help='save as tex file as well',       default=False,   action='store_true')
    parser.add_option(      '--rebin',       dest='rebin',       help='rebin factor',                   default=1,       type=int)
    parser.add_option('-l', '--lumi',        dest='lumi' ,       help='lumi to print out',              default=41.6,    type=float)
    parser.add_option(      '--only',        dest='only',        help='plot only these (csv)',          default='',      type='string')
    (opt, args) = parser.parse_args()

    #read list of samples
    jsonFile = open(opt.json,'r')
    samplesList=json.load(jsonFile,encoding='utf-8').items()
    jsonFile.close()
    systSamplesList=None
    if opt.systJson:
        jsonFile = open(opt.systJson,'r')
        systSamplesList=json.load(jsonFile,encoding='utf-8').items()
        jsonFile.close()

    onlyList=opt.only.split(',')

    #read plots 
    plots={}
    for slist,isSyst in [(samplesList,False),(systSamplesList,True)]:
        if slist is None : continue
        for tag,sample in slist: 

            if isSyst and not 't#bar{t}' in sample[3] : continue

            inDir=opt.inDir
            if isSyst : inDir += '/syst'
            fIn=ROOT.TFile.Open('%s/%s.root' % ( inDir, tag) )
            try:
                for tkey in fIn.GetListOfKeys():
                    key=tkey.GetName()
                    keep=False if len(onlyList)>0 else True
                    for tag in onlyList :
                        if tag in key: 
                            keep=True
                    if not keep: continue
                    obj=fIn.Get(key)
                    if not obj.InheritsFrom('TH1') : continue
                    if not key in plots : plots[key]=Plot(key)
                    if opt.rebin>1:  obj.Rebin(opt.rebin)
                    plots[key].add(h=obj,title=sample[3],color=sample[4],isData=sample[1],isSyst=isSyst)
            except:
                print 'Skipping %s'%tag


    #show plots
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)
    outDir=opt.inDir+'/plots'
    os.system('mkdir -p %s' % outDir)
    for p in plots : 
        if opt.saveLog    : plots[p].savelog=True
        if not opt.silent : plots[p].show(outDir=outDir,lumi=opt.lumi,saveTeX=opt.saveTeX)
        plots[p].appendTo(outDir+'/plotter.root')
        plots[p].reset()

    print '-'*50
    print 'Plots and summary ROOT file can be found in %s' % outDir
    print '-'*50

        
"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())

