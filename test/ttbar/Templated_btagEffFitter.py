import optparse
import os,sys
import json
import commands
import ROOT
import pickle
from plotter import Plot
from runTTbarAnalysis import LUMI

VARSTOFIT  = [('kindisc',-1,1),('jpTagger',0,2),('close_mlj',0,200)]
SLICEBINS  = [(20,320),(20,30),(30,50),(50,80),(80,120),(120,210),(210,320)]
SLICEVAR   = 'jetpt'
SYSTVARS   = ['','jesup','jesdn','jerup','jerdn','puup','pudn']

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
                    if flav=='data' and var!='' : continue
                    key='%s_pass%d_slice%d%s' % (flav, i, islice,systVar)
                    histos[key]=baseHisto.Clone(key)
                    histos[key].SetDirectory(0)
                    histos[key].Sumw2(0)
    baseHisto.Delete()

    #add files to the corresponding chains
    files = [ f for f in os.listdir(inDir) if '.root' in f ]
    chains={'mc':ROOT.TChain('kin')}
    for f in files: 
        key = 'mc' if 'MC' in f else 'data'
        chains[key].Add(inDir+'/'+f)

    #fill histos
    for key in chains:
        for i in xrange(0,chains[key].GetEntries()):
            chains[key].GetEntry(i)

            for systVar in SYSTVARS:
                
                wgtIdx, systIdx = 0, 0
                if systVar=='jesup' : wgtIdx, systIdx = 1, 1
                if systVar=='jesdn' : wgtIdx, systIdx = 2, 2
                if systVar=='jerup' : wgtIdx, systIdx = 3, 3
                if systVar=='jerdn' : wgtIdx, systIdx = 4, 4
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
def runSFFits(var,tagger,taggerDef,outDir):

    flavTemplates=[ ['b'], ['c','other'] ]
    if var=='jpTagger' : flavTemplates=[ ['b'], ['c'], ['other'] ]

    nOPs=len(taggerDef)-2
    bobs={}
    for iop in xrange(0,nOPs):

        bobs[iop]={}
        for islice in xrange(0,len(SLICEBINS)):

            totalExp=0
            pseudoData=None
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
                        if pseudoData is None:
                            pseudoData=h.Clone( name.replace(title,'total') )
                            pseudoData.SetDirectory(0)
                            pseudoData.Reset('ICE')
                        compExp=int(h.Integral())
                        totalExp += h.Integral()                        
                    
                        #full stats
                        pseudoData.Add(h)
                    
                        #toy
                        #for iev in xrange(0,ROOT.gRandom.Poisson(compExp)) : pseudoData.Fill(h.GetRandom())

                fIn.Close()

                #fix pseudoData uncertainties (not sure if needed)
                if syst=='' :
                    for xbin in xrange(1,pseudoData.GetXaxis().GetNbins()):
                        pseudoData.SetBinError(xbin,ROOT.TMath.Sqrt(pseudoData.GetBinContent(xbin)))

                #fraction fit
                fracFitter=ROOT.TFractionFitter(pseudoData,mc)
                #ROOT.SetOwnership( fracFitter, 0 ) 
                for ipar in xrange(0,mc.GetEntriesFast()) : 
                    expfrac=mc.At(ipar).Integral()/pseudoData.Integral()                    
                    fracFitter.Constrain(ipar,expfrac*0.5,1.0)
                fitRes=fracFitter.Fit()

                totalObs=pseudoData.Integral()
                bfrac,bfracErr  = ROOT.Double(0), ROOT.Double(0)
                fracFitter.GetResult(0,bfrac,bfracErr)
                bobs[iop][islice][syst]=(totalObs*bfrac,totalObs*bfracErr)

                if syst!='' : continue
                showResult(data=pseudoData,
                           mc=mc,
                           fracFitter=fracFitter,
                           name='%s_%s' % (tagger,baseName),
                           outDir='%s/%s_templates' % (outDir,var) )

                #fracFitter.Delete()


    #dump to pickle
    cache = '%s/%s_templates/.%s_fits.pck' % (outDir,var,tagger)
    cachefile = open(cache,'w')
    fitInfo={'var':var,'tagger':tagger,'taggerDef':taggerDef}
    pickle.dump(fitInfo, cachefile,pickle.HIGHEST_PROTOCOL)
    pickle.dump(bobs, cachefile,pickle.HIGHEST_PROTOCOL)
    cachefile.close()
    print 'Fit results have been stored in %s'%cache
    
"""
"""
def showResult(data,mc,fracFitter,name,outDir):
    
    pl=Plot(name)
    pl.add(data,'Data',1,True)
    totalObs=data.Integral()

    colors=[920,865,797]
    for i in xrange(0,mc.GetEntriesFast()):
        frac,fracErr  = ROOT.Double(0), ROOT.Double(0)
        fracFitter.GetResult(i,frac,fracErr)
        h=mc.At(i).Clone('flav_%d'%i)
        sfactor=totalObs*frac/h.Integral()
        h.Scale(sfactor)
        title='%s : %3.0f #pm %3.0f (SF=%3.2f)' % (h.GetTitle(),totalObs*frac,totalObs*fracErr,sfactor)
        h.SetTitle(title)
        pl.add(h,title,colors[i],False)
    pl.finalize()
    pl.show(outDir)
    #pl.reset()


"""
compares results obtained from different methods
"""
def showEfficiencyResults(fileList,outDir):

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
        grStat, grTotal = {}, {}
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
                    xmin, xmax = -40 +ifile*4, -35+ifile*4
                    grStat[iop]=ROOT.TGraphErrors()
                    grStat[iop].SetName('%s_%d_stat'%(var,iop))
                    grTotal[iop] = ROOT.TGraphErrors()
                    grTotal[iop].SetName('%s_%d_total'%(var,iop))
                    print var,tagger,'>',taggerDef[iop+2],eff

                #add point in corresponding slice
                np=grStat[iop].GetN()
                grStat[iop].SetPoint(np,0.5*(xmax+xmin),eff[''])
                grStat[iop].SetPointError(np,0.5*(xmax-xmin),eff['stat'])
                grTotal[iop].SetPoint(np,0.5*(xmax+xmin),eff[''])
                grTotal[iop].SetPointError(np,0.5*(xmax-xmin),eff['total'])

        #add collection
        if not (var in grStatColl):
            grStatColl[var]={}
            grTotalColl[var]={}
        grStatColl[var]  = grStat
        grTotalColl[var] = grTotal

    c=ROOT.TCanvas('c','c',800,500)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.025)
    c.SetLeftMargin(0.1)
    c.SetBottomMargin(0.12)
    colors=[ROOT.kAzure+9,ROOT.kOrange-1,ROOT.kGreen-5]
    markers=[20,21,22]
    firstvar=grStatColl.keys()[0]
    for iop in grStatColl[firstvar]:

        c.Clear()

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
            grTotalColl[var][iop].GetYaxis().SetRangeUser(0.0,1.5)
            grTotalColl[var][iop].GetYaxis().SetTitle('Efficiency')
            grTotalColl[var][iop].GetXaxis().SetTitle(SLICEVAR)
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

        line=ROOT.TLine(0,0,0,1.5)
        line.Draw()

        leg.Draw()
        txt=ROOT.TLatex()
        txt.SetNDC(True)
        txt.SetTextFont(43)
        txt.SetTextSize(16)
        txt.SetTextAlign(12)
        txt.DrawLatex(0.12,0.95,'#bf{CMS} #it{Preliminary} %3.1f fb^{-1} (13 TeV)' % (LUMI/1000.) )
        txt.DrawLatex(0.12,0.9,'%s > %3.3f' % (tagger,taggerDef[iop+2]))
        txt.DrawLatex(0.12,0.8,'#scale[0.8]{#it{inclusive}}')
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
    parser.add_option('-n', '--njobs',              dest='njobs',              help='# jobs to run in parallel',    default=0,           type='int')
    parser.add_option('-o', '--outDir',             dest='outDir',             help='output directory',             default='analysis',  type='string')
    (opt, args) = parser.parse_args()
    
    #if show result only
    if opt.show: 
        showEfficiencyResults(fileList=opt.show.split(','),outDir=opt.outDir)
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
    for var,_,_ in VARSTOFIT:
        for tagger,taggerDef in taggersList:
            runSFFits(var,tagger,taggerDef,opt.outDir)

    #all done here
    exit(0)



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
