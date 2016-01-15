import ROOT
from array import array
import sys
import optparse
import json
import pickle

"""
projects ttbar tree to estimate pt dependent efficiency
"""
def getEfficiencyCurves(tagger,taggerDef,tree):

    print 'Projecting tagging efficiency from ',tree.GetEntries(),' events for ',tagger,' tagger'

    effgrs={}
    ptBins = [0,20,25,30,35,40,45,50,60,70,80,90,100,120,140,160,180,200,250,300,400,500,600,800,1000]
    preTagH=ROOT.TH1F('preTagH',';p_{T} [GeV/c];',len(ptBins)-1,array('d',ptBins))
    preTagH.Sumw2()
    tagH=preTagH.Clone('tagH')

    nOPs=len(taggerDef)-2

    for flav,cond in [('b',   'abs(flavour)==5'),
                      ('c',   'abs(flavour)==4'),
                      ('udsg','abs(flavour)!=5 && abs(flavour)!=4')]:
        print '\t computing for',flav,cond

        if not flav in effgrs: effgrs[flav]={}

        for iop in xrange(1,nOPs):
            
            csvWP='%s>%f'%(tagger,taggerDef[iop+1])

            preTagH.Reset('ICE')
            tagH.Reset('ICE')
            tree.Draw("jetpt >> preTagH",'weight[0]*(%s)'%cond,'goff')
            tree.Draw('jetpt >> tagH',   'weight[0]*(%s && %s)'%(cond,csvWP),'goff')
            
            effgrs[flav][iop]=ROOT.TGraphAsymmErrors()
            effgrs[flav][iop].SetName('%s_%s_pass%d'%(flav,tagger,iop))
            effgrs[flav][iop].Divide(tagH,preTagH,'norm')

    return effgrs
            

"""
steer the script
"""
def main():
    
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-t', '--taggers',            dest='taggers'  ,          help='json with list of taggers',    default=None,               type='string')
    parser.add_option('-i', '--infile',             dest='infile',             help='input file',                   default=None,               type='string')
    parser.add_option('-o', '--output',             dest='output',             help='output file',                  default='data/.expEff.pck', type='string')
    (opt, args) = parser.parse_args()

    #read list of samples
    taggersFile = open(opt.taggers,'r')
    taggersList=json.load(taggersFile,encoding='utf-8').items()
    taggersFile.close()

    #project the efficiencies
    effGrs={}
    fIn=ROOT.TFile.Open(opt.infile)
    ftm=fIn.Get('ftm')
    for tagger,taggerDef in taggersList: 
        effGrs[tagger]=getEfficiencyCurves(tagger=tagger,taggerDef=taggerDef,tree=ftm)
    fIn.Close()

    #save all in a cache
    cachefile = open(opt.output,'w')
    pickle.dump(effGrs, cachefile,pickle.HIGHEST_PROTOCOL)
    cachefile.close()
    print 'Expected tagger efficiencies have been saved to %s'%opt.output


    #all done here
    exit(0)


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())

