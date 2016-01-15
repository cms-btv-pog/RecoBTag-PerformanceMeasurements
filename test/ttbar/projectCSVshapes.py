import optparse
import os,sys
import json
import commands
import ROOT
import pickle
from plotter import Plot
from array import array

from Templated_btagEffFitter import SLICEBINS,SLICEVAR
    
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
    parser.add_option('-s', '--sliceVar',           dest='sliceVar',           help='slicing variable',             default='jetpt',     type='string')
    parser.add_option('-i', '--inDir',              dest='inDir',              help='input directory with files',   default=None,        type='string')
    parser.add_option('-o', '--outDir',             dest='outDir',             help='output directory',             default='analysis',  type='string')
    parser.add_option(      '--channels',           dest='channels',           help='channels to use',              default='-121,-143,-169',  type='string')
    (opt, args) = parser.parse_args()
    
    global SLICEVAR
    SLICEVAR=opt.sliceVar

    #read list of samples
    taggersFile = open(opt.taggers,'r')
    taggersList=json.load(taggersFile,encoding='utf-8').items()
    taggersFile.close()

    #channels to filter
    channelList=[ int(k) for k in opt.channels.split(',') ]

    files=os.listdir(opt.inDir)
    chain=ROOT.TChain('kin')
    for f in files:
        if 'MC13TeV' in f:
            chain.AddFile(opt.inDir+'/'+f)

    #channels to filter                                                                                                                                                                                 
    channelList=[ int(k) for k in opt.channels.split(',') ]
    channelCond=''
    if len(channelList)>0:
        channelCond='('
        for ch in channelList:
            channelCond += 'ttbar_chan==%d ||'%ch
        channelCond=channelCond[:-2]
        channelCond +=')'
        
    histos=[]
    for islice in xrange(0,len(SLICEBINS[SLICEVAR])):
        cond='%s>%f && %s<%f'%(SLICEVAR,SLICEBINS[SLICEVAR][islice][0],SLICEVAR,SLICEBINS[SLICEVAR][islice][1])
        if channelCond !='' : cond += ' && ' + channelCond
        for tagger,taggerDef in taggersList:
            cond += '&& %s>0. && %s<%f' % (tagger,tagger,taggerDef[-1])
            for flavCond,flavH in [('abs(flavour)==5',ROOT.TH1F('%s_b_slice%d'%(tagger,islice),    ';%s;PDF'%taggerDef[0],50,0,taggerDef[-1])),
                                   ('abs(flavour)!=5',ROOT.TH1F('%s_other_slice%d'%(tagger,islice),';%s;PDF'%taggerDef[0],50,0,taggerDef[-1]))]:
                flavH.Sumw2()
                chain.Draw('%s>>%s'%(tagger,flavH.GetName()),
                           'weight[0]*(%s && %s)'%(cond,flavCond),
                           'goff')
                integ=flavH.Integral()
                if integ>0 : flavH.Scale(1./integ)
                flavH.SetDirectory(0)
                histos.append(flavH)

    fOut=ROOT.TFile('%s/taggersExpected.root'%opt.outDir,'RECREATE')
    for h in histos: h.Write()
    fOut.Close()

    #all done here
    exit(0)


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
