import optparse
import os,sys
import ROOT
import pickle
from rounding import *

def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',       dest='input',       help='input ROOT file with plots',   default=None,        type='string')
    (opt, args) = parser.parse_args()


    #
    # build templates
    #
    inF=ROOT.TFile.Open(opt.input)
    
    #control region
    crPlot='zll_met'
    dyTemplate=inF.Get('%s/%s'%(crPlot,crPlot))
    for key in inF.Get(crPlot).GetListOfKeys():
        keyName=key.GetName()
        if keyName==crPlot or 'DY' in keyName or 'Graph' in keyName: continue
        h=inF.Get('%s/%s'%(crPlot,keyName))        
        dyTemplate.Add( h,-1 )
    
    #signal region
    srPlot='ll_met'
    data=inF.Get('%s/%s'%(srPlot,srPlot))
    otherProc=None
    for key in inF.Get(srPlot).GetListOfKeys():
        keyName=key.GetName()
        if keyName==srPlot or 'Graph' in keyName : continue
        h=inF.Get('%s/%s'%(srPlot,keyName))
        if'DY' in keyName:
            dyTemplate.Scale(h.Integral()/dyTemplate.Integral())
        else:
            if otherProc is None:
                otherProc=h.Clone('%s_others'%crPlot)
            else:
                otherProc.Add(h)
    
    data.SetDirectory(0)
    otherProc.SetDirectory(0)
    otherProc.SetTitle('Non DY')
    dyTemplate.SetDirectory(0)
    dyTemplate.SetTitle('DY')
    inF.Close()


    #compile c++ wrapper to fit
    ROOT.gSystem.CompileMacro("TTbarSFbFitTools.cc","fk","libTTbarSFbFitTools")
    ROOT.gSystem.Load("libTTbarSFbFitTools.so")
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)

    #fit
    ttFracFitter=ROOT.TTbarFracFitter()
    saveResultIn = ROOT.TString('%s/dyFit'%os.path.dirname(opt.input))
    templates=ROOT.TObjArray()
    templates.Add(dyTemplate)
    templates.Add(otherProc)
    res=ttFracFitter.fit(templates,data,0,saveResultIn)

    dySF=(res.sf,res.sfUnc)    
    print '-'*50
    print 'DY fit from data'
    print '-'*50
    print 'Expected:',toLatexRounded(res.nExp,res.nExpUnc)
    print 'Observed:',toLatexRounded(res.nObs,res.nObsUnc)
    print 'Scale factor:',toLatexRounded(res.sf,res.sfUnc)
    print '-'*50

    #dump to file
    cache='%s/.dyScaleFactor.pck' % os.path.dirname(opt.input)
    cachefile=open(cache,'w')
    pickle.dump(dySF, cachefile, pickle.HIGHEST_PROTOCOL)
    cachefile.close()
    print 'Produced normalization cache for DY (%s)'%cache

    #all done here
    exit(0)



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
