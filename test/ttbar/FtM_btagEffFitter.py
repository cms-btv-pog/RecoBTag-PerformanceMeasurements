import optparse
import os,sys
import json
import commands
import ROOT
import pickle
from plotter import Plot

CHANNELS      = [-11*11,-13*13,-11*13]
#CHANNELS      = [-11*13]
JETMULTCATEGS = [2,3,4]
SLICEBINS     = [(20,320),(20,60),(60,120),(120,320)]
SLICEVAR      = 'jetpt'
SYSTVARS      = ['','jesup','jesdn','jerup','jerdn','trigdn','trigup','seldn','selup','qcdscaledn','qcdscaleup','hdampdn','hdampup']

"""
Project trees from files to build the templates
"""
def prepareTemplates(tagger,taggerDef,inDir,outDir):
    
    print '...starting %s'%tagger

    histos={}

    nOPs=len(taggerDef)-2
    nSliceCategs=(len(SLICEBINS)-1)**2+1

    #MC efficiency
    for key in ['b','c','l']:        
        for i in xrange(1,nOPs+1):            
            name='%s_%s_pass%d'%(key,tagger,i-1)
            histos[name]=ROOT.TH1F(name,';%s slice bin;Events'%SLICEVAR,len(SLICEBINS),0,len(SLICEBINS))
            for xbin in xrange(0,len(SLICEBINS)):
                label='%d-%d'%(SLICEBINS[xbin][0],SLICEBINS[xbin][1])
                histos[name].GetXaxis().SetBinLabel(xbin+1,label)

    #flavour categories
    flavourCombinationsBinMap=['l_{1}l_{2}','l_{1}c_{2}','l_{1}b_{2}','c_{1}l_{2}','c_{1}c_{2}','c_{1}b_{2}','b_{1}l_{2}','b_{1}c_{2}','b_{1}b_{2}']
    jetCategsBinMap=[]
    j1slice,j2slice=1,1
    for islice in xrange(0,nSliceCategs):
        j1Cuts,j2Cuts=SLICEBINS[0],SLICEBINS[0]
        if islice>0:
            if j2slice==len(SLICEBINS):
                j2slice=1
                j1slice+=1
            j1Cuts,j2Cuts=SLICEBINS[j1slice],SLICEBINS[j2slice]
            j2slice+=1
        if j1Cuts[0]<j2Cuts[0]:continue
        jetCategsBinMap.append( (j1Cuts,j2Cuts) )
    histos['flavcategs']=ROOT.TH2F('flavcategs',';Slice category;Flavour combination',
                                   len(jetCategsBinMap),0,len(jetCategsBinMap),9,0,9)
    for xbin in xrange(0,len(jetCategsBinMap)):
        j1Cuts,j2Cuts=jetCategsBinMap[xbin][0],jetCategsBinMap[xbin][1]
        label='(%d-%d),(%d-%d)'%(j1Cuts[0],j1Cuts[1],j2Cuts[0],j2Cuts[1])
        histos['flavcategs'].GetXaxis().SetBinLabel(xbin+1,label)
    for ybin in xrange(0,len(flavourCombinationsBinMap)):
        histos['flavcategs'].GetYaxis().SetBinLabel(ybin+1,flavourCombinationsBinMap[ybin])


    #tag counting in categories
    tagCountingBinMap=[]
    nJetMultCategs=len(JETMULTCATEGS)
    j1slice,j2slice=1,1
    for islice in xrange(0,nSliceCategs):

        j1Cuts,j2Cuts=SLICEBINS[0],SLICEBINS[0]
        if islice>0:
            if j2slice==len(SLICEBINS):
                j2slice=1
                j1slice+=1
            j1Cuts,j2Cuts=SLICEBINS[j1slice],SLICEBINS[j2slice]
            j2slice+=1

        if j1Cuts[0]<j2Cuts[0]:continue

        for ij in xrange(0,nJetMultCategs):
            jmult=JETMULTCATEGS[ij]
            
            for bmult in xrange(0,3):
                tagCountingBinMap.append( (bmult,jmult,j1Cuts,j2Cuts) )

    print len(tagCountingBinMap),' bins for tag counting...have fun with that'

    for key in ['data','hh','hl','ll']:
        
        for i in xrange(1,nOPs):            
            name='%s_%s_pass%d'%(key,tagger,i)
            histos[name]=ROOT.TH1F(name,';%s b-tag multiplicity;Events'%taggerDef[0],len(tagCountingBinMap),0,len(tagCountingBinMap))

            curJetMult=tagCountingBinMap[0][1]
            curJ1Cut=tagCountingBinMap[0][2]
            curJ2Cut=tagCountingBinMap[0][3]
            for xbin in xrange(1,len(tagCountingBinMap)+1):
                bmult=tagCountingBinMap[xbin-1][0]
                jmult=tagCountingBinMap[xbin-1][1]
                j1cut=tagCountingBinMap[xbin-1][2]
                j2cut=tagCountingBinMap[xbin-1][3]
                printJetMult=False
                if xbin==1 or jmult!=curJetMult:
                    printJetMult=True
                    curJetMult=jmult
                printJetCuts=False
                if xbin==1 or j1cut!=curJ1Cut or j2cut!=curJ2Cut:
                    printJetCuts=True
                    curJ1Cut=j1cut
                    curJ2Cut=j2cut

                label='%dt'%bmult
                if printJetMult : label += ',%dj'%jmult
                if printJetCuts : label='#splitline{%s}{(%d-%d),(%d-%d)}'%(label,j1cut[0],j1cut[1],j2cut[0],j2cut[1])
                histos[name].GetXaxis().SetBinLabel(xbin,label)

    #add files to the corresponding chains
    files = [ f for f in os.listdir(inDir) if '.root' in f ]
    chains={'mc':ROOT.TChain('ftm'),'data':ROOT.TChain('ftm')}
    for f in files: 
        key = 'mc' if 'MC' in f else 'data'
        chains[key].Add(inDir+'/'+f)

    #fill histos
    for key in chains:
        totalEntries=chains[key].GetEntries()
        for i in xrange(0,totalEntries):

            if i%100==0 : sys.stdout.write('\r [ %d/100 ] done for %s' %(int(float(100.*i)/float(totalEntries)),key) )

            chains[key].GetEntry(i)

            #require matching channel
            if not chains[key].ttbar_chan in CHANNELS : continue

            #require at least two jets
            if not chains[key].jetmult in JETMULTCATEGS : continue
            
            #event weight
            weight=chains[key].weight[0]

            ntags=[0]*nOPs
            nheavy=0
            flavourCombLabel=''
            for ij in xrange(0,2):
                    
                #tagger value
                taggerVal = getattr(chains[key],tagger)[ij]

                #count tags
                passTagWgts=[False]*nOPs
                for iop in xrange(1,nOPs):
                    if taggerVal<taggerDef[iop+1]: continue
                    passTagWgts[iop-1]=True
                    ntags[iop-1]+=1

                #MC truth
                flavName='l'
                if abs(chains[key].flavour[ij])==5 :
                    nheavy +=1 
                    flavName='b'
                if abs(chains[key].flavour[ij])==4: 
                    nheavy+=1
                    flavName='c'

                #MC truth for the efficiency as function of the slicing variable
                flavourCombLabel += '%s_{%d}'%(flavName,ij+1)
                jetSliceVarVal=getattr(chains[key],SLICEVAR)[ij]
                for ijcat in xrange(0,len(SLICEBINS)):
                    if jetSliceVarVal<SLICEBINS[ijcat][0] : continue
                    if jetSliceVarVal>SLICEBINS[ijcat][1] : continue 
                    histos['%s_%s_pass0'%(flavName,tagger)].Fill(ijcat,weight)
                    for iop in xrange(1,nOPs):
                        if not passTagWgts[iop-1] : continue
                        name='%s_%s_pass%d'%(flavName,tagger,iop)
                        histos[name].Fill(ijcat,weight)
                        
            #MC truth for the jet flavour combination vs jet slicing category
            if key !='data':
                for iflavComb in xrange(0,len(flavourCombinationsBinMap)):
                    if flavourCombLabel!= flavourCombinationsBinMap[iflavComb] : continue

                    for ijetCateg in xrange(0,len(jetCategsBinMap)):
                        j1Cuts=jetCategsBinMap[ijetCateg][0]
                        j1SliceVarVal=getattr(chains[key],SLICEVAR)[0]
                        if j1SliceVarVal<j1Cuts[0] or j1SliceVarVal>j1Cuts[1] : continue
                    
                        j2Cuts=jetCategsBinMap[ijetCateg][1]
                        j2SliceVarVal=getattr(chains[key],SLICEVAR)[1]
                        if j2SliceVarVal<j2Cuts[0] or j2SliceVarVal>j2Cuts[1] : continue
                
                        histos['flavcategs'].Fill(ijetCateg,iflavComb,weight)

            #tag counting histograms
            flavCat=key
            if key != 'data' :
                flavCat='hh'
                if nheavy==1: flavCat='hl'
                if nheavy==0: flavCat='ll'

            for ibin in xrange(0,len(tagCountingBinMap)):
                
                if chains[key].jetmult != tagCountingBinMap[ibin][1] : continue

                j1Cuts=tagCountingBinMap[ibin][2]
                j1SliceVarVal=getattr(chains[key],SLICEVAR)[0]
                if j1SliceVarVal<j1Cuts[0] or j1SliceVarVal>j1Cuts[1] : continue

                j2Cuts=tagCountingBinMap[ibin][3]
                j2SliceVarVal=getattr(chains[key],SLICEVAR)[1]
                if j2SliceVarVal<j2Cuts[0] or j2SliceVarVal>j2Cuts[1] : continue

                for iop in xrange(1,nOPs):
                    if ntags[iop-1]!=tagCountingBinMap[ibin][0] : continue
                    name='%s_%s_pass%d'%(flavCat,tagger,iop)
                    histos[name].Fill(ibin,weight)
                    
    #save templates to file
    fOut=ROOT.TFile.Open('%s/FtM/%s.root'%(outDir,tagger),'RECREATE')
    for key in histos : histos[key].Write()
    fOut.Close()


"""
Wrapper to be used when run in parallel
"""
def runPrepareTemplatesPacked(args):
    tagger, taggerDef, inDir, outDir = args
    try:
        return prepareTemplates(tagger=tagger,
                                taggerDef=taggerDef,                      
                                inDir=inDir,
                                outDir=outDir)
    except :
        print 50*'<'
        print "  Problem found (%s) baling out of this task" % sys.exc_info()[1]
        print 50*'<'
        return False

"""
Use the templates to prepare the workspace
"""
def prepareWorkspace(tagger,taggerDef,inDir):
    inF=ROOT.TFile.Open('%s/%s.root'%(inDir,tagger))

    colors = [ROOT.kGray,    ROOT.kAzure+7,  ROOT.kGreen, 
              ROOT.kGreen+1, ROOT.kOrange+8, ROOT.kMagenta+2,
              ROOT.kYellow-3, ROOT.kYellow-5, 0]

    #
    #MC EFFICIENCIES
    #
    effGrs={}
    for flav in ['b','c','l']:
        preTag=inF.Get('%s_%s_pass0' % (flav,tagger) )
        if not flav in effGrs: effGrs[flav]=[]
        for iop in xrange(1,len(taggerDef)-2):
            postTag=inF.Get('%s_%s_pass%d' % (flav,tagger,iop) )
            effGrs[flav].append( postTag.Clone() )
            effGrs[flav][-1].Sumw2()
            effGrs[flav][-1].SetTitle('%s>%3.2f' % (tagger,taggerDef[iop+2] ))
            effGrs[flav][-1].SetMarkerStyle(20+iop)
            effGrs[flav][-1].SetMarkerColor(colors[iop])
            effGrs[flav][-1].SetLineColor(colors[iop])
            effGrs[flav][-1].SetFillStyle(0)
            effGrs[flav][-1].Divide(preTag)
        
    #
    #FLAVOUR COMPOSITION
    #
    flavcategs=inF.Get('flavcategs')    
    catFracHistos=[]
    for xbin in xrange(1,flavcategs.GetNbinsX()+1):
        xlabel     = flavcategs.GetXaxis().GetBinLabel(xbin)
        totalInCat = flavcategs.Integral(xbin,xbin,1,flavcategs.GetNbinsY())

        for ybin in xrange(1,flavcategs.GetNbinsY()+1):
            ylabel=flavcategs.GetYaxis().GetBinLabel(ybin)

            #init histo if not yet available
            if len(catFracHistos)<ybin:
                catFracHistos.append( ROOT.TH1F('catfrac%d'%ybin,
                                                '%s;%s;Fraction' %(ylabel,flavcategs.GetTitle()),
                                                flavcategs.GetNbinsX(),flavcategs.GetXaxis().GetXmin(),flavcategs.GetXaxis().GetXmax()))
                catFracHistos[-1].SetTitle(ylabel)
                catFracHistos[-1].SetDirectory(0)
                catFracHistos[-1].Sumw2()
                catFracHistos[-1].SetLineColor(1)
                catFracHistos[-1].SetFillColor(colors[ybin-1])
                catFracHistos[-1].SetMarkerColor(colors[ybin-1])
                catFracHistos[-1].SetFillStyle(1001)
                for binCtr in xrange(1, catFracHistos[-1].GetNbinsX()+1):
                    catFracHistos[-1].GetXaxis().SetBinLabel(binCtr,flavcategs.GetXaxis().GetBinLabel(xbin))
            catFracHistos[ybin-1].SetBinContent(xbin,flavcategs.GetBinContent(xbin,ybin)/totalInCat)
            catFracHistos[ybin-1].SetBinError(xbin,flavcategs.GetBinError(xbin,ybin)/totalInCat)       
    


    #SHOW PLOTS
    ceff=ROOT.TCanvas('ceff','ceff',500,500)
    ceff.SetTopMargin(0)
    ceff.SetBottomMargin(0)
    ceff.SetLeftMargin(0)
    ceff.SetRightMargin(0)

    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(43)
    txt.SetTextSize(16)
    txt.SetTextAlign(12)

    ceff.cd()
    p1=ROOT.TPad('p1','p1',0.,0.,1.0,0.33)
    p1.SetBottomMargin(0.15)
    p1.SetTopMargin(0.01)
    p1.SetLeftMargin(0.12)
    p1.SetRightMargin(0.05)
    p1.Draw()
    p1.cd()
    for i in xrange(0,len(effGrs['b'])):
        drawOpt='E1X0' if i==0 else 'E1X0same'
        effGrs['b'][i].Draw(drawOpt)
        effGrs['b'][i].GetXaxis().SetTitleSize(0.9)
        effGrs['b'][i].GetXaxis().SetLabelSize(0.08)
        effGrs['b'][i].GetYaxis().SetRangeUser(0.12,0.96)
        effGrs['b'][i].GetYaxis().SetTitleSize(0.09)
        effGrs['b'][i].GetYaxis().SetLabelSize(0.08)
        effGrs['b'][i].GetYaxis().SetTitle('Efficiency')
        effGrs['b'][i].GetYaxis().SetTitleOffset(0.6)
    txt.DrawLatex(0.85,0.93,'#bf{[b]}')

    ceff.cd()
    p2=ROOT.TPad('p2','p2',0.,0.33,1.0,0.66)
    p2.SetBottomMargin(0.01)
    p2.SetTopMargin(0.01)
    p2.SetLeftMargin(0.12)
    p2.SetRightMargin(0.05)
    p2.Draw()
    p2.cd()
    for i in xrange(0,len(effGrs['c'])):
        drawOpt='E1X0' if i==0 else 'E1X0same'
        effGrs['c'][i].Draw(drawOpt)
        effGrs['c'][i].GetYaxis().SetRangeUser(0.12,0.96)
        effGrs['c'][i].GetYaxis().SetTitleSize(0.09)
        effGrs['c'][i].GetYaxis().SetLabelSize(0.08)
        effGrs['c'][i].GetYaxis().SetTitle('Efficiency')
        effGrs['c'][i].GetYaxis().SetTitleOffset(0.6)
    txt.DrawLatex(0.85,0.93,'#bf{[c]}')

    ceff.cd()
    p3=ROOT.TPad('p3','p3',0.,0.66,1.0,1.0)
    p3.SetBottomMargin(0.01)
    p3.SetTopMargin(0.02)
    p3.SetLeftMargin(0.12)
    p3.SetRightMargin(0.05)
    p3.Draw()
    p3.cd()
    leg=ROOT.TLegend(0.2,0.75,0.8,0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetNColumns(4)
    leg.SetTextSize(0.06)
    for i in xrange(0,len(effGrs['l'])):
        drawOpt='E1X0' if i==0 else 'E1X0same'
        effGrs['l'][i].Draw(drawOpt)
        effGrs['l'][i].GetYaxis().SetRangeUser(0.12,0.96)
        effGrs['l'][i].GetYaxis().SetTitleSize(0.09)
        effGrs['l'][i].GetYaxis().SetLabelSize(0.08)
        effGrs['l'][i].GetYaxis().SetTitle('Efficiency')
        effGrs['l'][i].GetYaxis().SetTitleOffset(0.6)
        leg.AddEntry(effGrs['l'][i],effGrs['l'][i].GetTitle(),'p')
    leg.Draw()
    txt.DrawLatex(0.9,0.93,'#bf{[l]}')
   
    txt.DrawLatex(0.2,0.93,'#bf{CMS} #it{Simulation}')

    ceff.cd()
    ceff.Modified()
    ceff.Update()
    raw_input()


    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.02)
    c.SetRightMargin(0.1)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.15)
    leg=ROOT.TLegend(0.2,0.75,0.8,0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetNColumns(4)
    leg.SetTextSize(0.04)

    stack=ROOT.THStack()
    for h in catFracHistos:        
        stack.Add(h,'h')
        leg.AddEntry(h,h.GetTitle(),'f')
    stack.Draw('hist')  
    stack.GetYaxis().SetTitle('Fraction')
    stack.GetXaxis().SetTitle('%s category'%SLICEVAR)
    stack.GetXaxis().SetTitleOffset(2)
    stack.GetYaxis().SetRangeUser(0,1)
    leg.Draw()
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(43)
    txt.SetTextSize(16)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.2,0.93,'#bf{CMS} #it{Simulation}')

    raw_input()
        



"""
steer the script
"""
def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    #ROOT.gROOT.SetBatch(True)
    
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-t', '--taggers',            dest='taggers'  ,          help='json with list of taggers',    default=None,        type='string')
    parser.add_option('-i', '--inDir',              dest='inDir',              help='input directory with files',   default=None,        type='string')
    parser.add_option('-l', '--lumi',               dest='lumi' ,              help='lumi to print out',            default=41.6,        type=float)
    parser.add_option('-n', '--njobs',              dest='njobs',              help='# jobs to run in parallel',    default=0,           type='int')
    parser.add_option('-o', '--outDir',             dest='outDir',             help='output directory',             default='analysis',  type='string')
    parser.add_option(      '--recycleTemplates',   dest='recycleTemplates',   help='do not regenerate templates',  default=False,       action='store_true')
    (opt, args) = parser.parse_args()
    
    #read list of samples
    taggersFile = open(opt.taggers,'r')
    taggersList=json.load(taggersFile,encoding='utf-8').items()
    taggersFile.close()
   
    #re-create templates 
    if not opt.recycleTemplates:
   
        task_list=[]
        os.system('mkdir -p %s/FtM'%(opt.outDir))
        for tagger,taggerDef in taggersList:
            task_list.append((tagger,taggerDef,opt.inDir,opt.outDir))
    
        print '%s jobs to run in %d parallel threads' % (len(task_list), opt.njobs)
        if opt.njobs == 0:
            for tagger,taggerDef,inDir,outDir in task_list:
                prepareTemplates(tagger=tagger,
                                 taggerDef=taggerDef,
                                 inDir=inDir,
                                 outDir=outDir)
        else:
            from multiprocessing import Pool
            pool = Pool(opt.njobs)
            pool.map(runPrepareTemplatesPacked, task_list)

    #prepare workspace
    for tagger,taggerDef in taggersList:
        prepareWorkspace(tagger=tagger,taggerDef=taggerDef,inDir=opt.outDir+'/FtM')

    #all done here
    exit(0)


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
