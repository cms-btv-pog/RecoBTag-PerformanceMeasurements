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

    tagCountingBinMap=[]
    nJetMultCategs=len(JETMULTCATEGS)
    nSliceCategs=(len(SLICEBINS)-1)**2+1
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

    #tag counting
    nOPs=len(taggerDef)-2
    for key in ['data','hh','hl','ll']:
        for i in xrange(1,nOPs+1):            
            name='%s_%s_pass%d'%(key,tagger,i-1)
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

            #count tags
            for iop in xrange(1,nOPs):

                ntags,nheavy=0,0
                for ij in xrange(0,2):
                    taggerVal = getattr(chains[key],tagger)[ij]
                    if taggerVal>taggerDef[iop+1]:ntags+=1
                    if abs(chains[key].flavour[ij])==5 or abs(chains[key].flavour[ij])==4: nheavy+=1

                name='%s_%s_pass%d'%(key,tagger,iop)
                if key != 'data' :
                    flavCat='hh'
                    if nheavy==1: flavCat='hl'
                    if nheavy==0: flavCat='ll'
                    name='%s_%s_pass%d'%(flavCat,tagger,iop)

                for ibin in xrange(0,len(tagCountingBinMap)):

                    if ntags != tagCountingBinMap[ibin][0] : continue
                    if chains[key].jetmult != tagCountingBinMap[ibin][1] : continue
                    j1SliceVarVal=getattr(chains[key],SLICEVAR)[0]
                    j2SliceVarVal=getattr(chains[key],SLICEVAR)[1]

                    leadSliceVarVal,trailSliceVarVal=j1SliceVarVal,j2SliceVarVal
                    if j2SliceVarVal>j1SliceVarVal:
                        leadSliceVarVal,trailSliceVarVal=j2SliceVarVal,j1SliceVarVal

                    if leadSliceVarVal  < tagCountingBinMap[ibin][2][0] or leadSliceVarVal  > tagCountingBinMap[ibin][2][1] : continue
                    if trailSliceVarVal < tagCountingBinMap[ibin][3][0] or trailSliceVarVal > tagCountingBinMap[ibin][3][1] : continue
                    
                    histos[name].Fill(ibin,weight)
                    
                    #flav='other'
                    #if abs(chains[key].flavour[ij])==5: flav='b'
                    #if abs(chains[key].flavour[ij])==4: flav='c'
                    #if key=='data' : flav='data'

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
    (opt, args) = parser.parse_args()
    
    #read list of samples
    taggersFile = open(opt.taggers,'r')
    taggersList=json.load(taggersFile,encoding='utf-8').items()
    taggersFile.close()
    
    #re-create templates
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

    #all done here
    exit(0)



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
