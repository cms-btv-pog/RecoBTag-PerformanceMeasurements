import optparse
import os,sys
import json
import commands
import ROOT
import pickle
from plotter import Plot

VARSTOFIT  = [('kindisc',-1,1),('close_mlj',0,250)]  #,('jpTagger',0,2)
SLICEBINS  = [(20,320),(20,60),(60,120),(120,320)]
SLICEVAR   = 'jetpt'
SYSTVARS   = ['','jesup','jesdn','jerup','jerdn','trigdn','trigup','seldn','selup','qcdscaledn','qcdscaleup','hdampdn','hdampup']

"""
Project trees from files to build the templates
"""
def prepareTemplates(tagger,taggerDef,var,varRange,channeList,inDir,outDir):
    
    print '...starting %s for %s'%(tagger,var)

    histos={}

    nOPs=len(taggerDef)-2

    #prepare templates
    baseHisto=ROOT.TH1F(var,';Discriminator;Jets',50,varRange[0],varRange[1])
    for flav in ['b','c','other','data']:
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
                
                if not in channelList : continue

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
                sliceVarVal = getattr(chains[key],SLICEVAR)
                varVal      = getattr(chains[key],var)      if var=='jpTagger'   else getattr(chains[key],var)[systIdx]
                taggerVal   = getattr(chains[key],tagger)

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
                    for iop in xrange(2,len(taggerDef)-1):
                        if taggerVal< taggerDef[iop] : continue
                        hkey='%s_pass%d_slice%d%s'%(flav,iop-1,islice,systVar)
                        histos[hkey].Fill(varVal,weight)

    #save templates to file
    fOut=ROOT.TFile.Open('%s/%s_templates/%s.root'%(outDir,var,tagger),'RECREATE')
    for key in histos : histos[key].Write()
    fOut.Close()


"""
Wrapper to be used when run in parallel
"""
def runPrepareTemplatesPacked(args):
    tagger, taggerDef, var, varRange, channelList, inDir, outDir = args
    try:
        return prepareTemplates(tagger=tagger,
                                taggerDef=taggerDef,
                                var=var,
                                varRange=varRange,
                                channelList,
                                inDir=inDir,
                                outDir=outDir)
    except :
        print 50*'<'
        print "  Problem found (%s) baling out of this task" % sys.exc_info()[1]
        print 50*'<'
        return False

"""
Leave no bins with 0 counts
If negative values (negative weights in MC) set to minimum 
"""
def checkTemplate(h,minVal=1e-5):
    for xbin in xrange(1,h.GetNbinsX()+1):
        y=h.GetBinContent(xbin)
        if y<=0: h.SetBinContent(xbin,minVal)

"""
run the fits
"""
def runSFFits(var,tagger,taggerDef,outDir,lumi):

    flavTemplates=[ ['b'], ['c','other'] ]
    if var=='jpTagger' : flavTemplates=[ ['b'], ['c'], ['other'] ]

    nOPs=len(taggerDef)-2
    bobs={}
    bexp={}
    for iop in xrange(0,nOPs):

        bobs[iop]={}
        bexp[iop]={}
        for islice in xrange(0,len(SLICEBINS)):

            totalExp=0
            data,pseudoData=None,None
            bobs[iop][islice]={}
            baseNameNominal = 'pass%d_slice%d' % (iop,islice)
            for syst in SYSTVARS:

                #build flavour templates
                baseName = baseNameNominal+syst

                mc=ROOT.TObjArray( len(flavTemplates) )
                fIn=ROOT.TFile.Open('%s/%s_templates/%s.root' % (outDir, var, tagger) )
                for flavGroup in flavTemplates :

                    title='+'.join(flavGroup)
                    name='%s_%s'%(title,baseName)

                    #build combined template
                    h=None
                    for flav in flavGroup:
                        hname='%s_%s'%(flav,baseName)
                        ihisto=fIn.Get(hname)
                        checkTemplate(ihisto)

                        hnomname='%s_%s'%(flav,baseNameNominal)
                        inomhisto=fIn.Get(hnomname)
                        checkTemplate(inomhisto)

                        #only interested in shape variations
                        ihisto.Scale(inomhisto.Integral()/ihisto.Integral())

                        if h is None:
                            h=ihisto.Clone(name)        
                            h.SetDirectory(0)
                            h.SetTitle(title)
                        else:
                            h.Add( ihisto )

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
                if len(syst)==0:
                    res=ttFracFitter.fit(mc,data,0,saveResultIn)
                    bobs[iop][islice]['']=(res.nObs,res.nObsUnc)
                    bexp[iop][islice]=(res.nExp,res.nExpUnc)
                    res=ttFracFitter.fit(mc,pseudoData,0,'')
                    bobs[iop][islice]['closureup']=(res.nObs,res.nObsUnc)
                    bobs[iop][islice]['closuredn']=(res.nObs,res.nObsUnc)
                else:
                    res=ttFracFitter.fit(mc,pseudoData,0,saveResultIn)
                    bobs[iop][islice][syst]=(res.nObs,res.nObsUnc)



    #dump to pickle
    cache = '%s/%s_templates/.%s_fits.pck' % (outDir,var,tagger)
    cachefile = open(cache,'w')
    fitInfo={'var':var,'tagger':tagger,'taggerDef':taggerDef,'slicevar':SLICEVAR,'slicebins':SLICEBINS}
    pickle.dump(fitInfo, cachefile,pickle.HIGHEST_PROTOCOL)
    pickle.dump(bobs, cachefile, pickle.HIGHEST_PROTOCOL)
    pickle.dump(bexp, cachefile, pickle.HIGHEST_PROTOCOL)
    cachefile.close()
    print 'Fit results have been stored in %s'%cache
    
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
    parser.add_option('-n', '--njobs',              dest='njobs',              help='# jobs to run in parallel',    default=0,           type='int')
    parser.add_option('-o', '--outDir',             dest='outDir',             help='output directory',             default='analysis',  type='string')
    parser.add_option(      '--channels',           dest='channels',           help='channels to use',              default='-121,-143,-169',  type='string')
    (opt, args) = parser.parse_args()
    
    #read list of samples
    taggersFile = open(opt.taggers,'r')
    taggersList=json.load(taggersFile,encoding='utf-8').items()
    taggersFile.close()

    #channels to filter
    channelList=[ int(k) for k in opt.channels.split(',') ]

    #re-create templates
    if not opt.recycleTemplates:
        task_list=[]
        for var,varMin,varMax in VARSTOFIT:
            os.system('mkdir -p %s/%s_templates'%(opt.outDir,var))
            for tagger,taggerDef in taggersList:
                if var==tagger : continue
                task_list.append((tagger,taggerDef,var,(varMin,varMax),channelList,opt.inDir,opt.outDir))
        #task_list=list(set(task_list))
        print '%s jobs to run in %d parallel threads' % (len(task_list), opt.njobs)
        #run the analysis jobs                                                                                                                                                          
        if opt.njobs == 0:
            for tagger,taggerDef,var,varRange, channelList,inDir,outDir in task_list:
                prepareTemplates(tagger=tagger,
                                 taggerDef=taggerDef,
                                 var=var,
                                 varRange=varRange,
                                 channelList=channelList,
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
