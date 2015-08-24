import optparse
import os,sys
import json
import commands
import ROOT
import pickle
from plotter import Plot

VARSTOFIT  = [('kindisc',-1,1),('jpTagger',0,2),('close_mlj',0,200)]
SLICEBINS  = [(20,320),(20,60),(60,120),(120,320)]
SLICEVAR   = 'jetpt'
SYSTVARS   = ['','jesup','jesdn','jerup','jerdn','puup','pudn','trigdn','trigup','seldn','selup','qcdscaledn','qcdscaleup','hdampdn','hdampup']

"""
Project trees from files to build the templates
"""
def prepareTemplates(tagger,taggerDef,var,varRange,inDir,outDir):
    
    print '...starting %s for %s'%(tagger,var)

    #prepare templates
    baseHisto=ROOT.TH1F(var,';Discriminator;Jets',50,varRange[0],varRange[1])
    histos={}
    for flav in ['b','c','other','data']:
        nOPs=len(taggerDef)-2
        for i in xrange(0,nOPs):
            for islice in xrange(0,len(SLICEBINS)):
                for systVar in SYSTVARS:
                    if flav=='data' and len(systVar)>0 : continue
                    key='%s_pass%d_slice%d%s' % (flav, i, islice,systVar)
                    histos[key]=baseHisto.Clone(key)
                    histos[key].SetDirectory(0)
                    histos[key].Sumw2(0)
    baseHisto.Delete()

    #add files to the corresponding chains
    files = [ f for f in os.listdir(inDir) if '.root' in f ]
    chains={'mc':ROOT.TChain('kin'),'data':ROOT.TChain('kin')}
    for f in files: 
        key = 'mc' if 'MC' in f else 'data'
        chains[key].Add(inDir+'/'+f)

    #fill histos
    for key in chains:
        for i in xrange(0,chains[key].GetEntries()):
            chains[key].GetEntry(i)

            for systVar in SYSTVARS:
                
                if key=='data' and len(systVar)>0 : continue

                wgtIdx, systIdx = 0, 0
                if systVar=='jesup'      : wgtIdx, systIdx = 1,  1
                if systVar=='jesdn'      : wgtIdx, systIdx = 2,  2
                if systVar=='jerup'      : wgtIdx, systIdx = 3,  3
                if systVar=='jerdn'      : wgtIdx, systIdx = 4,  4
                if systVar=='pudn'       : wgtIdx, systIdx = 5,  0
                if systVar=='puup'       : wgtIdx, systIdx = 6,  0
                if systVar=='trigdn'     : wgtIdx, systIdx = 7,  0
                if systVar=='trigup'     : wgtIdx, systIdx = 8,  0
                if systVar=='seldn'      : wgtIdx, systIdx = 9,  0
                if systVar=='selup'      : wgtIdx, systIdx = 10, 0
                if systVar=='qcdscaledn' : wgtIdx, systIdx = 11, 0
                if systVar=='qcdscaleup' : wgtIdx, systIdx = 12, 0
                if systVar=='hdampdn'    : wgtIdx, systIdx = 13, 0
                if systVar=='hdampup'    : wgtIdx, systIdx = 14, 0

                #event weight
                weight      = chains[key].weight[wgtIdx]

                #no need to proceed if event is not selected
                if weight==0: continue

                #variable to slice on, variable to be fit, and tagger to apply
                sliceVarVal = getattr(chains[key],SLICEVAR)[systIdx] if SLICEVAR=='jetpt' else getattr(chains[key],SLICEVAR) 
                varVal      = getattr(chains[key],var)               if var=='jpTagger'   else getattr(chains[key],var)[systIdx]
                taggerVal   = getattr(chains[key],tagger+'Tagger')

                #determine categories
                passSlice=[]
                for islice in xrange(0,len(SLICEBINS)):
                    if sliceVarVal<=SLICEBINS[islice][0] or sliceVarVal>SLICEBINS[islice][1] : continue
                    passSlice.append(islice)

                #assign flavour
                flav='other'
                if chains[key].flavour==3: flav='b'
                if chains[key].flavour==2: flav='c'
                if key=='data' : flav='data'

                #fill the histos
                for islice in passSlice:
                    hkey='%s_pass0_slice%d%s'%(flav,islice,systVar)
                    histos[hkey].Fill(varVal,weight)
                    for i in xrange(2,len(taggerDef)-1):
                        if taggerVal< taggerDef[i] : continue
                        hkey='%s_pass%d_slice%d%s'%(flav,i-1,islice,systVar)
                        histos[hkey].Fill(varVal,weight)

    #save templates to file
    fOut=ROOT.TFile.Open('%s/%s_templates/%s.root'%(outDir,var,tagger),'RECREATE')
    for key in histos : histos[key].Write()
    fOut.Close()


"""
Wrapper to be used when run in parallel
"""
def runPrepareTemplatesPacked(args):
    tagger, taggerDef, var, varRange, inDir, outDir = args
    try:
        return prepareTemplates(tagger=tagger,
                                taggerDef=taggerDef,
                                var=var,
                                varRange=varRange,
                                inDir=inDir,
                                outDir=outDir)
    except :
        print 50*'<'
        print "  Problem found (%s) baling out of this task" % sys.exc_info()[1]
        print 50*'<'
        return False

"""
run the fits
"""
def runSFFits(var,tagger,taggerDef,outDir,lumi):

    flavTemplates=[ ['b'], ['c','other'] ]
    if var=='jpTagger' : flavTemplates=[ ['b'], ['c'], ['other'] ]

    #custom fraction fitter class
    ttFracFitter=ROOT.TTbarFracFitter();

    nOPs=len(taggerDef)-2
    bobs={}
    for iop in xrange(0,nOPs):

        bobs[iop]={}
        for islice in xrange(0,len(SLICEBINS)):

            totalExp=0
            data,pseudoData=None,None
            bobs[iop][islice]={}
            for syst in SYSTVARS:

                #build flavour templates
                baseName='pass%d_slice%d%s' % (iop,islice,syst)
                mc=ROOT.TObjArray( len(flavTemplates) )
                fIn=ROOT.TFile.Open('%s/%s_templates/%s.root' % (outDir, var, tagger) )
                for flavGroup in flavTemplates :

                    title='+'.join(flavGroup)
                    name='%s_%s'%(title,baseName)

                    #build combined template
                    h=None
                    for flav in flavGroup:
                        hname='%s_%s'%(flav,baseName)
                        if h is None:
                            h=fIn.Get(hname).Clone(name)        
                            h.SetDirectory(0)
                            h.SetTitle(title)
                        else:
                            h.Add( fIn.Get(hname) )

                    #check for empty bins...
                    for xbin in xrange(1,h.GetXaxis().GetNbins()+1):
                        cts=h.GetBinContent(xbin)
                        if cts>0: continue
                        h.SetBinContent(xbin,1e-4)

                    #add for template fit
                    mc.Add(h)

                    #add to pseudo-data
                    if syst=='':
                        if data is None:
                            data=fIn.Get('data_%s'%baseName).Clone()
                            data.SetDirectory(0)
                            data.SetTitle('data')
                        if pseudoData is None:
                            pseudoData=h.Clone( name.replace(title,'total') )
                            pseudoData.SetDirectory(0)
                            pseudoData.Reset('ICE')
                        compExp=int(h.Integral())
                        totalExp += h.Integral()                        
                    
                        #full stats
                        pseudoData.Add(h)

                fIn.Close()

                #fit
                saveResultIn = ROOT.TString('%s/%s_templates/%s_%s'%(outDir,var,tagger,baseName) if syst=='' else '')
                res=ttFracFitter.fit(mc,data,0,saveResultIn)
                bobs[iop][islice][syst]=(res.nObs,res.nObsUnc)


    #dump to pickle
    cache = '%s/%s_templates/.%s_fits.pck' % (outDir,var,tagger)
    cachefile = open(cache,'w')
    fitInfo={'var':var,'tagger':tagger,'taggerDef':taggerDef}
    pickle.dump(fitInfo, cachefile,pickle.HIGHEST_PROTOCOL)
    pickle.dump(bobs, cachefile,pickle.HIGHEST_PROTOCOL)
    cachefile.close()
    print 'Fit results have been stored in %s'%cache
    
"""
compares results obtained from different methods
"""
def showEfficiencyResults(fileList,outDir,lumi):

    grIncStatColl={}
    grIncTotalColl={}
    grStatColl={}
    grTotalColl={}

    for ifile in xrange(0,len(fileList)):

        f=fileList[ifile]
        cachefile = open(f,'r')
        fitInfo=pickle.load(cachefile)
        fitResults=pickle.load(cachefile)
        cachefile.close()
    
        var       = fitInfo['var']
        tagger    = fitInfo['tagger']
        taggerDef = fitInfo['taggerDef']
        nOPs=len(taggerDef)-2
        grIncStat, grIncTotal = {}, {}
        grStat,    grTotal = {}, {}
        for iop in xrange(1,nOPs):
            for islice in fitResults[iop]:

                if fitResults[0][islice][''][0]==0: continue
 
                #statistical uncertainty : unc #tag / # pre-tag
                eff={
                    ''     : fitResults[iop][islice][''][0]/fitResults[0][islice][''][0],
                    'stat' : fitResults[iop][islice][''][1]/fitResults[0][islice][''][0]
                    }

                #alternative variations
                for systVar in fitResults[iop][islice]:
                    
                    if systVar=='' : continue

                    #symmetrize differences for up/dn variations
                    newEff = fitResults[iop][islice][systVar][0]/fitResults[0][islice][systVar][0]
                    key    = systVar.replace('up','')
                    key    = key.replace('dn','')
                    newEff = newEff-eff['']
                    if key in eff : eff[key]=0.5*(ROOT.TMath.Abs(eff[key])+ROOT.TMath.Abs(newEff))
                    else          : eff[key]=newEff
                
                #total uncertainty (stat+syst)
                totalUnc=0
                for systVar in eff:
                    if systVar=='' : continue
                    totalUnc += eff[systVar]**2
                eff['total']=ROOT.TMath.Sqrt(totalUnc)

                xmin=SLICEBINS[islice][0]+ifile*4
                xmax=SLICEBINS[islice][1]+ifile*4

                if islice==0 : 
                    grIncStat[iop]=ROOT.TGraphErrors()
                    grIncStat[iop].SetName('%s_%d_inc_stat'%(var,iop))
                    grIncStat[iop].SetPoint(0,0+ifile*4,eff[''])
                    grIncStat[iop].SetPointError(0,5,eff['stat'])

                    grIncTotal[iop] = ROOT.TGraphErrors()
                    grIncTotal[iop].SetName('%s_%d_inc_total'%(var,iop))
                    grIncTotal[iop].SetPoint(0,0+ifile*4,eff[''])
                    grIncTotal[iop].SetPointError(0,5,eff['total'])

                    #init differential plots also
                    grStat[iop]=ROOT.TGraphErrors()
                    grStat[iop].SetName('%s_%d_stat'%(var,iop))
                    grTotal[iop] = ROOT.TGraphErrors()
                    grTotal[iop].SetName('%s_%d_total'%(var,iop))

                else:
                    np=grStat[iop].GetN()
                    grStat[iop].SetPoint(np,0.5*(xmax+xmin),eff[''])
                    grStat[iop].SetPointError(np,0.5*(xmax-xmin),eff['stat'])
                    grTotal[iop].SetPoint(np,0.5*(xmax+xmin),eff[''])
                    grTotal[iop].SetPointError(np,0.5*(xmax-xmin),eff['total'])

        #add collection
        if not (var in grStatColl):
            grIncStatColl[var]={}
            grIncTotalColl[var]={}
            grStatColl[var]={}
            grTotalColl[var]={}
        grIncStatColl[var]  = grIncStat
        grIncTotalColl[var] = grIncTotal
        grStatColl[var]  = grStat
        grTotalColl[var] = grTotal



    #show results in canvas
    c=ROOT.TCanvas('c','c',800,500)
    c.SetRightMargin(0)
    c.SetLeftMargin(0)
    c.SetBottomMargin(0)
    c.SetTopMargin(0)

    colors=[ROOT.kAzure+9,ROOT.kOrange-1,ROOT.kGreen-5]
    markers=[20,21,22]

    clabel=ROOT.TPad('clabel','clabel',0,0.95,1.0,1.0)
    clabel.SetRightMargin(0.0)
    clabel.SetTopMargin(0.0)
    clabel.SetLeftMargin(0.0)
    clabel.SetBottomMargin(0.0)

    cinc=ROOT.TPad('cinc','cinc',0,0,0.2,0.95)
    cinc.SetRightMargin(0.02)
    cinc.SetTopMargin(0.01)
    cinc.SetLeftMargin(0.3)
    cinc.SetBottomMargin(0.12)

    cexc=ROOT.TPad('cexc','cexc',0.2,0,1.0,0.95)
    cexc.SetRightMargin(0.02)
    cexc.SetTopMargin(0.01)
    cexc.SetLeftMargin(0.02)
    cexc.SetBottomMargin(0.12)

    c.cd()
    clabel.Draw()
    c.cd()
    cinc.Draw()
    c.cd()
    cexc.Draw()
    c.cd()
    
    firstvar=grIncStatColl.keys()[0]
    for iop in grIncStatColl[firstvar]:

        cinc.cd()
        cinc.Clear()
        ivar=0
        drawOpt='a2'
        for var in grIncTotalColl:
            grIncTotalColl[var][iop].SetTitle(var)
            grIncTotalColl[var][iop].SetMarkerColor(colors[ivar]+1)
            grIncTotalColl[var][iop].SetLineColor(colors[ivar]+1)
            grIncTotalColl[var][iop].SetFillColor(colors[ivar]+1)
            grIncTotalColl[var][iop].SetMarkerStyle(1)
            grIncTotalColl[var][iop].SetFillStyle(3001)
            grIncTotalColl[var][iop].Draw(drawOpt)
            grIncTotalColl[var][iop].GetYaxis().SetRangeUser(0.0,1.5)
            grIncTotalColl[var][iop].GetYaxis().SetTitle('Efficiency')
            grIncTotalColl[var][iop].GetYaxis().SetTitleSize(0.1)
            grIncTotalColl[var][iop].GetYaxis().SetLabelSize(0.1)
            grIncTotalColl[var][iop].GetXaxis().SetTitleSize(0.)
            grIncTotalColl[var][iop].GetXaxis().SetLabelSize(0.)
            grIncTotalColl[var][iop].GetXaxis().SetNdivisions(0)
            grIncTotalColl[var][iop].GetXaxis().SetTitle(SLICEVAR)
            drawOpt='2'
            ivar+=1

        ivar=0
        for var in grIncStatColl:
            grIncStatColl[var][iop].SetMarkerColor(colors[ivar])
            grIncStatColl[var][iop].SetLineColor(colors[ivar])
            grIncStatColl[var][iop].SetFillColor(0)
            grIncStatColl[var][iop].SetMarkerStyle(markers[ivar])
            grIncStatColl[var][iop].SetFillStyle(0)
            grIncStatColl[var][iop].Draw('p')
            grIncStatColl[var][iop].SetTitle(var)
            ivar+=1

        cexc.cd()
        cexc.Clear()
        leg = ROOT.TLegend(0.6, 0.95,0.95,0.8)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(43)
        leg.SetTextSize(16)
        ivar=0
        drawOpt='a2'
        for var in grTotalColl:
            grTotalColl[var][iop].SetTitle(var)
            grTotalColl[var][iop].SetMarkerColor(colors[ivar]+1)
            grTotalColl[var][iop].SetLineColor(colors[ivar]+1)
            grTotalColl[var][iop].SetFillColor(colors[ivar]+1)
            grTotalColl[var][iop].SetMarkerStyle(1)
            grTotalColl[var][iop].SetFillStyle(3001)
            grTotalColl[var][iop].Draw(drawOpt)
            grTotalColl[var][iop].GetYaxis().SetTitleSize(0.)
            grTotalColl[var][iop].GetYaxis().SetLabelSize(0.)
            grTotalColl[var][iop].GetYaxis().SetRangeUser(0.0,1.5)
            drawOpt='2'
            ivar+=1

        ivar=0
        for var in grStatColl:
            grStatColl[var][iop].SetMarkerColor(colors[ivar])
            grStatColl[var][iop].SetLineColor(colors[ivar])
            grStatColl[var][iop].SetFillColor(0)
            grStatColl[var][iop].SetMarkerStyle(markers[ivar])
            grStatColl[var][iop].SetFillStyle(0)
            grStatColl[var][iop].Draw('p')
            grStatColl[var][iop].SetTitle(var)
            leg.AddEntry(grStatColl[var][iop],grStatColl[var][iop].GetTitle(),'p')
            ivar+=1
        leg.Draw()

        clabel.cd()
        clabel.Clear()
        txt=ROOT.TLatex()
        txt.SetNDC(True)
        txt.SetTextFont(43)
        txt.SetTextSize(16)
        txt.SetTextAlign(12)
        if lumi<100:
            txt.DrawLatex(0.05,0.5,'#bf{CMS} #it{Preliminary} %3.1f pb^{-1} (13 TeV)' % lumi)
        else:
            txt.DrawLatex(0.05,0.5,'#bf{CMS} #it{Preliminary} %3.1f fb^{-1} (13 TeV)' % (lumi/1000.))
        txt.DrawLatex(0.8,0.5,'%s > %3.3f #scale[0.8]{#it{inclusive}}' % (tagger,taggerDef[iop+2]))

        c.cd()
        c.Modified()
        c.Update()
        for ext in ['png','pdf']: c.SaveAs('%s/%s_%d.%s'% (outDir,tagger,iop,ext))



"""
steer the script
"""
def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)
    
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-t', '--taggers',            dest='taggers'  ,          help='json with list of taggers',    default=None,        type='string')
    parser.add_option('-i', '--inDir',              dest='inDir',              help='input directory with files',   default=None,        type='string')
    parser.add_option('-v', '--var',                dest='var',                help='templated variable',           default='kindisc',   type='string')
    parser.add_option(      '--recycleTemplates',   dest='recycleTemplates',   help='recycleTemplates',             default=False,       action='store_true')
    parser.add_option(      '--show',               dest='show',               help='show eff. results',            default=None,        type='string')
    parser.add_option('-l', '--lumi',               dest='lumi' ,              help='lumi to print out',            default=41.6,        type=float)
    parser.add_option('-n', '--njobs',              dest='njobs',              help='# jobs to run in parallel',    default=0,           type='int')
    parser.add_option('-o', '--outDir',             dest='outDir',             help='output directory',             default='analysis',  type='string')
    (opt, args) = parser.parse_args()
    
    #if show result only
    if opt.show: 
        showEfficiencyResults(fileList=opt.show.split(','),outDir=opt.outDir,lumi=opt.lumi)
        exit(0)                 

    #read list of samples
    taggersFile = open(opt.taggers,'r')
    taggersList=json.load(taggersFile,encoding='utf-8').items()
    taggersFile.close()
    

    #re-create templates
    if not opt.recycleTemplates:
        task_list=[]
        for var,varMin,varMax in VARSTOFIT:
            os.system('mkdir -p %s/%s_templates'%(opt.outDir,var))
            for tagger,taggerDef in taggersList:
                if var==tagger : continue
                task_list.append((tagger,taggerDef,var,(varMin,varMax),opt.inDir,opt.outDir))
        #task_list=list(set(task_list))
        print '%s jobs to run in %d parallel threads' % (len(task_list), opt.njobs)
        #run the analysis jobs                                                                                                                                                          
        if opt.njobs == 0:
            for tagger,taggerDef,var,varRange,inDir,outDir in task_list:
                prepareTemplates(tagger=tagger,
                                 taggerDef=taggerDef,
                                 var=var,
                                 varRange=varRange,
                                 inDir=inDir,
                                 outDir=outDir)
        else:
            from multiprocessing import Pool
            pool = Pool(opt.njobs)
            pool.map(runPrepareTemplatesPacked, task_list)

    #run the fits
    ROOT.gSystem.CompileMacro("TTbarSFbFitTools.cc","fk","libTTbarSFbFitTools")
    ROOT.gSystem.Load("libTTbarSFbFitTools.so")
    for var,_,_ in VARSTOFIT:
        for tagger,taggerDef in taggersList:
            runSFFits(var,tagger,taggerDef,opt.outDir,opt.lumi)

    #all done here
    exit(0)



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
