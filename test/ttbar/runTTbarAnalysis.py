import optparse
import os,sys
import json
import commands
import ROOT
import pickle
import math
from array import array
from storeTools import *
from rounding import *

CHANNELS={-11*11:'ll', -13*13:'ll', -11*13:'emu'}

#configuration
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-j', '--json',        dest='json'  ,      help='json with list of files',      default=None,        type='string')
parser.add_option('-i', '--inDir',       dest='inDir',       help='input directory with files',   default=None,        type='string')
parser.add_option('-o', '--outDir',      dest='outDir',      help='output directory',             default='analysis',  type='string')
parser.add_option(      '--only',        dest='only',        help='process only matching (csv)',  default='',          type='string')
parser.add_option(      '--tmvaWgts',    dest='tmvaWgts',    help='tmva weights',                 default=None,        type='string')
parser.add_option(      '--dyScale',     dest='dyScale',     help='DY scale factor',              default=None,        type='string')
parser.add_option('-n', '--njobs',       dest='njobs',       help='# jobs to run in parallel',    default=0,           type='int')
parser.add_option('--doTwoTag',          dest='doTwoTag',    help='Run 2TagCount (default:  False)', default=0,        type='int')
parser.add_option('--keepGoodFilesBySize',     dest='keepGoodFilesBySize',       help='Review output file sizes before running and skip good ones',    default=1,           type='int')
parser.add_option('--inspectOutput',     dest='inspectOutput',       help='Review output file sizes and re-run baddies',    default=1,           type='int')
(opt, args) = parser.parse_args()


"""
Perform the analysis on a single file
"""
def runTTbarAnalysis(inFile, outFile, wgt, tmvaWgts=None,isData=False):

    from ROOT import TTbarEventAnalysis
    evAnalysis=TTbarEventAnalysis()


    print "opt.doTwoTag",opt.doTwoTag
    evAnalysis.setTwoTagCount(opt.doTwoTag)

    #MC specifics
    if 'TTJets' in inFile: evAnalysis.setReadTTJetsGenWeights(True)

    if isData:
        if 'MuonEG'   in inFile : 
            evAnalysis.addTriggerBit(0,-11*13)
            evAnalysis.addTriggerBit(1,-11*13)
        if 'DoubleElectron' in inFile : 
            evAnalysis.addTriggerBit(2,-11*11)
        if 'DoubleMuon' in inFile :
            evAnalysis.addTriggerBit(3,-13*13)
    else:
            evAnalysis.addTriggerBit(0,-11*13)
            evAnalysis.addTriggerBit(1,-11*13)
            evAnalysis.addTriggerBit(2,-11*11)
            evAnalysis.addTriggerBit(3,-13*13)


    for v in ['close_mlj[0]', 'close_dphi', 'close_deta', 'close_lj2ll_dphi', 'close_lj2ll_deta',
              'far_mlj',      'far_dphi',   'far_deta',   'far_lj2ll_dphi',   'far_lj2ll_deta',
              'j2ll_dphi',    'j2ll_deta']: 
        evAnalysis.addVarForTMVA(ROOT.TString(v))    
    if not (tmvaWgts is None) : evAnalysis.setTMVAWeightsBaseDir(tmvaWgts)
    evAnalysis.prepareOutput(ROOT.TString(outFile))
    evAnalysis.processFile(ROOT.TString(inFile),wgt,isData)
    evAnalysis.finalizeOutput()

"""
Wrapper to be used when run in parallel
"""
def runTTbarAnalysisPacked(args):
    inFile, outFile, wgt, tmvaWgts,isData = args
    try:
        return runTTbarAnalysis(inFile=inFile, outFile=outFile, wgt=wgt, tmvaWgts=tmvaWgts,isData=isData)
    except :
        print 50*'<'
        print "  Problem  (%s) with %s continuing without"%(sys.exc_info()[1],inFile)
        print 50*'<'
        return False


def RunJobs(task_list,nJobs):
    if nJobs == 0:
        for inFile, outFile,wgt, tmvaWgts,isData in task_list: 
            runTTbarAnalysis(inFile=inFile, outFile=outFile, wgt=wgt, tmvaWgts=tmvaWgts, isData=isData)
    else:
        from multiprocessing import Pool
        pool = Pool(nJobs)
        pool.map(runTTbarAnalysisPacked, task_list)

def InspectFiles(outputDir,runTags):
    zero=0
    tiny=0
    norm=0
    
    okFiles=[]
    bustedFiles=[]
    for tag in runTags:
        os.system('ls -l %s/%s_*.root > %s.txt' % (outputDir,tag,tag) )
        sizeFile=open('%s.txt' % (tag))
        linesofls=sizeFile.readlines()
        sizeFile.close()
        os.system('rm %s.txt' % (tag))

        for line in linesofls:
            try:
                items=line.split()
                sizeOfFile=int(items[4])
                fileName=str(items[8])
                if sizeOfFile==0:
                    zero=zero+1
                    bustedFiles.append(fileName)
                elif sizeOfFile<2000:
                    tiny=tiny+1
                    bustedFiles.append(fileName)
                else:
                    norm=norm+1
                    okFiles.append(str(items[8]))
            except:
                print "can't read",line
        
        print tag,"zero, tiny, norm",zero,tiny,norm
        
        print "len bustedFiles",len(bustedFiles)

    return okFiles,bustedFiles


"""
steer the script
"""
def main():

    #compile c++ wrapper to run over trees 
    ROOT.gSystem.Load("libJetMETCorrectionsObjects.so")
    ROOT.gSystem.CompileMacro("TTbarEventAnalysis.cc","fkgd","libTTbarEventAnalysis");
    ROOT.gSystem.Load("libTTbarEventAnalysis.so")
    
    #read list of samples
    jsonFile = open(opt.json,'r')
    samplesList=json.load(jsonFile,encoding='utf-8').items()
    jsonFile.close()

    #prepare output
    if len(opt.outDir) : os.system('mkdir -p %s' % opt.outDir)

    #only list
    onlyList=opt.only.split(',')    

    #read normalization
    xsecWgts, integLumi = {}, {}
    cache='%s/src/RecoBTag/PerformanceMeasurements/test/ttbar/data/.xsecweights.pck'%os.environ['CMSSW_BASE']
    try:
        cachefile = open(cache, 'r')
        xsecWgts  = pickle.load(cachefile)
        integLumi = pickle.load(cachefile)
        cachefile.close()        
        print 'Normalization read from cache (%s)' % cache

        for tag,sample in samplesList:
            if not tag in xsecWgts:
                raise KeyError

    except:
        print '(Re-)Computing original number of events and storing in cache, this may take a while if it\'s the first time'
        print 'Current cache contains already %d processes'%len(xsecWgts)
        xsecWgts, integLumi = produceNormalizationCache(samplesList=samplesList,inDir=opt.inDir,cache=cache, xsecWgts=xsecWgts, integLumi=integLumi)


    for key in xsecWgts:
        print key,xsecWgts[key]

    print
    for key in integLumi:
        print key,integLumi[key]


    #DY scale factor
    if opt.dyScale:
        cachefile=open(opt.dyScale,'r')
        dySF=pickle.load(cachefile)
        cachefile.close()
        for tag in xsecWgts:
            if not 'DY' in tag: continue
            print tag,xsecWgts[tag].GetBinContent(1),' -> ',
            xsecWgts[tag].Scale(dySF[0])
            print xsecWgts[tag].GetBinContent(1)

    
    if opt.keepGoodFilesBySize:
        listOfSampleTags=[]
        for tag,sample in samplesList:
            listOfSampleTags.append(tag)
        goodFiles,badFiles=InspectFiles(opt.outDir,listOfSampleTags)
    else:
        goodFiles=[]

    #create the analysis jobs
    runTags = []
    task_list = []
    for tag,sample in samplesList:
        print tag
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
            if outF in goodFiles:
                continue
            task_list.append( (input_list[nf], outF, wgt, opt.tmvaWgts, sample[1]) )

    task_list=list(set(task_list))
    print '%s jobs to run in %d parallel threads' % (len(task_list), opt.njobs)

    #run the analysis jobs
    RunJobs(task_list,opt.njobs)


    if opt.inspectOutput:
        nLoop=0
        re_task_list=[]
        while(nLoop==0 or len(re_task_list)>0):
            print "nLoop",nLoop
            if nLoop>0:
                print "going at it again",len(re_task_list)
                RunJobs(re_task_list,opt.njobs)
                re_task_list=[]

                goodFiles,badFiles=InspectFiles(opt.outDir,runTags)
                
                for task in task_list:
                    print "task",task
                    outFileNameFromTask=str(task[1])
                    if outFileNameFromTask in bustedFiles:
                        re_task_list.append(task) 
            nLoop=nLoop+1


    #merge the outputs
    for tag in runTags:
        os.system('hadd -f %s/%s.root %s/%s_*.root' % (opt.outDir,tag,opt.outDir,tag) )

    for tag in runTags:
        os.system('rm %s/%s_*.root' % (opt.outDir,tag) )
    print 'Analysis results are available in %s' % opt.outDir

    #all done here
    exit(0)



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
