import optparse
import os,sys
import json
import commands
import ROOT
import pickle
import math
from array import array
from storeTools import getEOSlslist
from systTools import getTriggerEfficiency,getLeptonSelectionEfficiencyScaleFactor,getJetEnergyScales, getJetResolutionScales

CHANNELS={-11*11:'ll', -13*13:'ll', -11*13:'emu'}

"""
Perform the analysis on a single file
"""
def runTTbarAnalysis(inFile, outFile, wgt, tmvaWgts=None):

    from ROOT import TTbarEventAnalysis
    evAnalysis=TTbarEventAnalysis()

    #MC specifics
    if 'MC13TeV_DYJetsToLL' in inFile or 'MC13TeV_WJets' in inFile : evAnalysis.setUseOnlySignOfGenWeight(True)
    if 'MC13TeV_TTJets' in inFile: evAnalysis.setReadTTJetsGenWeights(True)

    #Data specifics
    if 'Data' in inFile:
        evAnalysis.setApplyTriggerEff(False)
        evAnalysis.setApplyLepSelEff(False)
        if 'MuonEG'   in inFile : 
            evAnalysis.addTriggerBit(0,-11*13)
            evAnalysis.addTriggerBit(1,-11*13)
        if 'DoubleEG' in inFile : 
            evAnalysis.addTriggerBit(2,-11*11)
        if 'DoubleMu' in inFile :
            evAnalysis.addTriggerBit(3,-13*13)

    for v in ['close_mlj','close_ptrel','close_dphi','close_deta','far_mlj','far_ptrel','far_dphi','far_deta']: evAnalysis.addVarForTMVA(ROOT.TString(v))    
    if not (tmvaWgts is None) : evAnalysis.setTMVAWeightsFile(tmvaWgts)
    evAnalysis.prepareOutput(ROOT.TString(outFile))
    evAnalysis.processFile(ROOT.TString(inFile),wgt)
    evAnalysis.finalizeOutput()

"""
Wrapper to be used when run in parallel
"""
def runTTbarAnalysisPacked(args):
    inFile, outFile, wgt, tmvaWgts = args
    try:
        return runTTbarAnalysis(inFile=inFile, outFile=outFile, wgt=wgt, tmvaWgts=tmvaWgts)
    except :
        print 50*'<'
        print "  Problem  (%s) with %s continuing without"%(sys.exc_info()[1],inFile)
        print 50*'<'
        return False


"""
steer the script
"""
from rounding import *
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-j', '--json',        dest='json'  ,      help='json with list of files',      default=None,        type='string')
    parser.add_option('-i', '--inDir',       dest='inDir',       help='input directory with files',   default=None,        type='string')
    parser.add_option('-o', '--outDir',      dest='outDir',      help='output directory',             default='analysis',  type='string')
    parser.add_option(      '--only',        dest='only',        help='process only matching (csv)',  default='',          type='string')
    parser.add_option(      '--tmvaWgts',    dest='tmvaWgts',    help='tmva weights',             default=None,  type='string')
    parser.add_option('-n', '--njobs',       dest='njobs',       help='# jobs to run in parallel',    default=0,           type='int')
    (opt, args) = parser.parse_args()

    #compile c++ wrapper to run over trees 
    ROOT.gSystem.CompileMacro("TTbarEventAnalysis.cc","fkgd","libTTbarEventAnalysis");
    ROOT.gSystem.Load("libTTbarEventAnalysis.so")

    #read list of samples
    print opt.json
    jsonFile = open(opt.json,'r')
    samplesList=json.load(jsonFile,encoding='utf-8').items()
    jsonFile.close()

    #prepare output
    if len(opt.outDir) : os.system('mkdir -p %s' % opt.outDir)

    #only list
    onlyList=opt.only.split(',')    

    #read normalization
    xsecWgts, integLumi = {}, {}
    cache='%s/.xsecweights.pck'%opt.outDir
    try:
        cachefile = open(cache, 'r')
        xsecWgts  = pickle.load(cachefile)
        integLumi = pickle.load(cachefile)
        cachefile.close()        
        print 'Normalization read from cache (%s)' % cache
    except:
        print 'Computing original number of events and storing in cache, this may take a while if it\'s the first time'
        for tag,sample in samplesList: 

            if sample[1]==1 : 
                xsecWgts[tag]=1.0
                continue

            input_list=getEOSlslist(directory=opt.inDir+'/'+tag)            
            xsec=sample[0]            
            norigEvents=0
            for f in input_list:
                fIn=ROOT.TFile.Open(f)
                norigEvents+=fIn.Get('allEvents/hEventCount').GetBinContent(1)
                fIn.Close()
            xsecWgts[tag]  = xsec/norigEvents  if norigEvents>0 else 0
            integLumi[tag] = norigEvents/xsec  if norigEvents>0 else 0
            print '... %s cross section=%f pb #orig events=%d lumi=%3.2f/fb' % (tag,xsec,norigEvents,integLumi[tag]/1000.)

        #dump to file
        cachefile=open(cache,'w')
        pickle.dump(xsecWgts, cachefile, pickle.HIGHEST_PROTOCOL)
        pickle.dump(integLumi, cachefile, pickle.HIGHEST_PROTOCOL)
        cachefile.close()
        print 'Produced normalization cache (%s)'%cache

    #create the analysis jobs
    runTags = []
    task_list = []
    for tag,_ in samplesList:

        #check if in list
        if len(onlyList)>0:
            veto=True
            for selTag in onlyList:
                if selTag in tag: veto=False
            if veto : continue

        runTags.append(tag)
        input_list=getEOSlslist(directory=opt.inDir+'/'+tag)
        wgt = xsecWgts[tag]
        for nf in xrange(0,len(input_list)) : 
            outF='%s/%s_%d.root'%(opt.outDir,tag,nf)
            task_list.append( (input_list[nf],outF,wgt, opt.tmvaWgts) )

    task_list=list(set(task_list))
    print '%s jobs to run in %d parallel threads' % (len(task_list), opt.njobs)

    #run the analysis jobs
    if opt.njobs == 0:
        for inFile, outFile,wgt, tmvaWgts in task_list: 
            runTTbarAnalysis(inFile=inFile, outFile=outFile, wgt=wgt, tmvaWgts=tmvaWgts)
    else:
        from multiprocessing import Pool
        pool = Pool(opt.njobs)
        pool.map(runTTbarAnalysisPacked, task_list)

    #merge the outputs
    for tag in runTags:
        os.system('hadd -f %s/%s.root %s/%s_*.root' % (opt.outDir,tag,opt.outDir,tag) )
        os.system('rm %s/%s_*.root' % (opt.outDir,tag) )
    print 'Analysis results are available in %s' % opt.outDir

    #all done here
    exit(0)



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
