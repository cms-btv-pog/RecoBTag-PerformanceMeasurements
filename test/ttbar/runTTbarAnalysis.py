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
def runTTbarAnalysis(inFile, outFile, wgt, taggers, tmvaWgts=None):

    applyMETFilters=False
    applyTriggerEff=True
    applyLepSelEff=True
    applyPUWgt=1    # 0=no, 1=from vertices, 2=from gen PU
    readTTJetsGenWeights=False
    triggerBits=[]
    allowedChannels=[]

    #MC specifics
    useOnlySignOfGenWeight=False
    if 'MC13TeV_DYJetsToLL' in inFile or 'MC13TeV_WJets' in inFile : useOnlySignOfGenWeight=True
    if 'MC13TeV_TTJets' in inFile: 
        readTTJetsGenWeights=True
        print 'Will read generator level weights for ', inFile

    #Data specifics
    if 'Data' in inFile:
        #applyMETFilters=True
        applyPUWgt=0
        applyTriggerEff=False
        applyLepSelEff=False
        if 'DoubleMu' in inFile : 
            allowedChannels.append(-13*13)
            triggerBits.append(3)
        if 'MuonEG'   in inFile : 
            allowedChannels.append(-11*13)
            triggerBits.append(0)
            triggerBits.append(1)        
        if 'DoubleEG' in inFile : 
            allowedChannels.append(-11*11)
            triggerBits.append(2)

    #prepare output and histograms
    outF=ROOT.TFile.Open(outFile,'RECREATE')
    baseHistos={
        'npvinc' : ROOT.TH1F('npvinc', ';N_{PV,good}-N_{HS};Events',              50, 0, 50),
        'npv'    : ROOT.TH1F('npv',    ';N_{PV,good}-N_{HS};Events',              50, 0, 50),
        'rho'    : ROOT.TH1F('rho',    ';#rho [GeV];Events',                      50, 0, 30),
        'mll'    : ROOT.TH1F('mll',    ';Dilepton invariant mass [GeV];Events',   20, 0, 300),
        'mllinc' : ROOT.TH1F('mllinc',    ';Dilepton invariant mass [GeV];Events',   20, 0, 300),
        'met'    : ROOT.TH1F('met',    ';Missing transverse energy [GeV];Events', 15, 0, 300),
        'njets'  : ROOT.TH1F('njets',  ';Jet multiplicity;Events;',               6,  2, 8),
        'leadjpt': ROOT.TH1F('leadjpt',';Leading jet p_{T} [GeV];Events;',        14,30,300),
        'leadlpt': ROOT.TH1F('leadlpt',';Leading lepton p_{T} [GeV];Events;',     9,20,200),
        'trailjpt': ROOT.TH1F('trailjpt',';Trailing jet p_{T} [GeV];Events;',     14,30,300),
        'traillpt': ROOT.TH1F('traillpt',';Trailing lepton p_{T} [GeV];Events;',  9,20,200),
        'evsel': ROOT.TH1F('evsel',';Event selection;Events;', 5,0,5)
        }
    baseHistos['evsel'].GetXaxis().SetBinLabel(1,'Z,#geq 2j')    
    baseHistos['evsel'].GetXaxis().SetBinLabel(2,'#geq 2j')
    baseHistos['evsel'].GetXaxis().SetBinLabel(3,'=2j')
    baseHistos['evsel'].GetXaxis().SetBinLabel(4,'=3j')
    baseHistos['evsel'].GetXaxis().SetBinLabel(5,'#geq4j')

    #init KIN tree
    kinTree=ROOT.TTree('kin','kin training')
    kinTree.SetDirectory(outF)
    kinVars={'flavour'    : (array( 'i', [ 0 ]),   ';Jet flavour;Jets',                   4,  0, 4 ),
             'jetpt'      : (array( 'f', [ 0.]*5), ';Transverse momentum [GeV]; Jets',   50,  0, 200),
             'jeteta'     : (array( 'f', [ 0.]),   ';Pseudo-rapidity; Jets',              25, 0, 2.5),
             'close_mlj'  : (array( 'f', [ 0.]*5), ';M(lepton,jet) [GeV]; Jets',          50, 0, 250),
             'close_deta' : (array( 'f', [ 0.]), ';#Delta#eta(lepton,jet); Jets',       50, 0, 4),
             'close_dphi' : (array( 'f', [ 0.]), ';#Delta#phi(lepton,jet) [rad]; Jets', 50, 0, 3.15),
             'close_ptrel': (array( 'f', [ 0.]), ';p_{T}^{rel}(lepton,jet) [GeV];Jets', 50, 0, 1),
             'far_mlj'    : (array( 'f', [ 0.]), ';M(lepton,jet) [GeV]; Jets',          50, 0, 250),
             'far_deta'   : (array( 'f', [ 0.]), ';#Delta#eta(lepton,jet); Jets',       50, 0, 4),
             'far_dphi'   : (array( 'f', [ 0.]), ';#Delta#phi(lepton,jet) [rad]; Jets', 50, 0, 3.15),
             'far_ptrel'  : (array( 'f', [ 0.]), ';p_{T}^{rel}(lepton,jet) [GeV];Jets', 50, 0, 1),
             'kindisc'    : (array( 'f', [ 0.]*5), ';Kinematics discriminator;Jets',      50, -1, 1),}

    #b-tagging variables (read from json)
    taggersFile = open(taggers,'r')
    taggersList = json.load(taggersFile,encoding='utf-8').items()
    taggersFile.close()
    for tagger, taggerDef in taggersList :
        title, firstOP, lastOP = taggerDef[0], taggerDef[2], taggerDef[-1]
        kinVars[str(tagger+'Tagger')]     = (array('f',[0.]), ';%s;Jets'%title,             50, ROOT.TMath.Min(firstOP,0), lastOP )
    kinVars['weight'] = (array('f',[0.]*15),)

    #add to tree and init histograms
    for v in kinVars.iterkeys():
        vtype = 'I' if v=='flavour' else 'F'
        if v=='kindisc' or v=='close_mlj' or v=='jetpt':
            kinTree.Branch( v, kinVars[v][0], '%s[5]/%s' % (v,vtype) )
        elif v=='weight':
            kinTree.Branch( v, kinVars[v][0], '%s[15]/%s' % (v,vtype) )
        else:
            kinTree.Branch( v, kinVars[v][0], '%s/%s' % (v,vtype) )
        try:
            baseHistos[v] = ROOT.TH1F(v,kinVars[v][1],kinVars[v][2],kinVars[v][3],kinVars[v][4])
        except:
            pass
    baseHistos['flavour'].GetXaxis().SetBinLabel(1,'unmatched')
    baseHistos['flavour'].GetXaxis().SetBinLabel(2,'udsg')
    baseHistos['flavour'].GetXaxis().SetBinLabel(3,'c')
    baseHistos['flavour'].GetXaxis().SetBinLabel(4,'b')

    #init a tmva reader
    tmvaReader=None if tmvaWgts is None else ROOT.TMVA.Reader()
    tmvaVars={}
    ivar=0
    for v in ['close_mlj[0]', 'close_ptrel', 'close_dphi', 'close_deta',
              'far_mlj',      'far_ptrel',   'far_dphi',   'far_deta']:
        tmvaVars[v]=array('f',[0.])
        if not tmvaReader is None : tmvaReader.AddVariable(v,tmvaVars[v])
        ivar+=1
    if not tmvaReader is None : tmvaReader.BookMVA('BDT',tmvaWgts)

    #replicate histos per channel
    histos={}
    for ch in CHANNELS.itervalues():
        for key in baseHistos:
            tag='%s_%s' % (ch,key)
            histos[tag]=baseHistos[key].Clone(tag)
            histos[tag].Sumw2()
            histos[tag].SetDirectory(outF)

    #loop over events
    inF=ROOT.TFile.Open(inFile)
    tree=inF.Get("btagana/ttree")
    nentries=tree.GetEntriesFast()
    print '....opened %s -> analysing %d events -> %s' % (inFile,nentries,outFile)
    for i in xrange(0,nentries): 
        tree.GetEntry(i)

        #progress bar
        if i%100==0 : print '\r[ %3d/100 ] to completion' % int(100.*i/nentries),

        #generator level weights
        genWgt=1.0 if tree.ttbar_nw==0 else tree.ttbar_w[0]
        if useOnlySignOfGenWeight:
            genWgt=1.0 
            if tree.ttbar_w[0]<0: genWgt=-1
        qcdScaleLo,qcdScaleHi,hdampLo,hdampHi = 1.0, 1.0, 1.0, 1.0
        if readTTJetsGenWeights and tree.ttbar_nw>17:
            qcdScaleLo=tree.ttbar_w[9]
            qcdScaleHi=tree.ttbar_w[5]
            hdampLo=tree.ttbar_w[tree.ttbar_nw-17]
            hdampHi=tree.ttbar_w[tree.ttbar_nw-9]

        #pileup weight
        puWgtLo, puWgtNom, puWgtHi = 1.0, 1.0, 1.0

        #
        # MET filters
        #
        if applyMETFilters:
            if tree.ttbar_metfilterWord!=0 : continue
        
        #
        # TRIGGER
        #
        hasTrigger=True if len(triggerBits)==0 else False
        for bit in triggerBits:
            hasTrigger |= ((tree.ttbar_trigWord>>bit) & 1)
        if not hasTrigger : continue

        #
        # CHANNEL
        #
        ch=''
        try:
            ch=CHANNELS[tree.ttbar_chan]
        except:
            continue
        if len(allowedChannels)>0:
            if not tree.ttbar_chan in allowedChannels:
                continue

        #
        # LEPTONS
        #
        if tree.ttbar_nl<2 : continue
        lp4=[]
        for il in xrange(0,tree.ttbar_nl):
            lp4.append( ROOT.TLorentzVector() )
            lp4[-1].SetPtEtaPhiM(tree.ttbar_lpt[il],tree.ttbar_leta[il],tree.ttbar_lphi[il],0.)
        mll=(lp4[0]+lp4[1]).M()

        #trigger efficiency weight
        trigWgtLo, trigWgtNom, trigWgtHi = 1.0, 1.0, 1.0
        if applyTriggerEff :
            eff,effUnc=getTriggerEfficiency(pt1=lp4[0].Pt(),eta1=lp4[0].Eta(),
                                            pt2=lp4[1].Pt(),eta2=lp4[1].Eta(),
                                            ch=tree.ttbar_chan)
            trigWgtLo=eff-effUnc
            trigWgtNom=eff
            trigWgtHi=eff+effUnc

        #lepton selection efficiency
        lepSelEffLo, lepSelEffNom, lepSelEffHi = 1.0, 1.0, 1.0
        if applyLepSelEff:
            for il in xrange(0,2) :
                lepSF, lepSFUnc = getLeptonSelectionEfficiencyScaleFactor(pt=lp4[il].Pt(),
                                                                          eta=lp4[il].Eta(),
                                                                          pdgId=tree.ttbar_lid[il])
                lepSelEffLo  *= (lepSF-lepSFUnc)
                lepSelEffNom *= lepSF
                lepSelEffHi  *= (lepSF+lepSFUnc)
                
        #nominal event weight
        evWgt = wgt*puWgtNom*trigWgtNom*lepSelEffNom*genWgt
        histos[ch+'_npvinc'].Fill(tree.nPV-1,evWgt)
        
        #
        # JET/MET SELECTION
        #
        jetCount=[0]*5
        jetDiff=[ROOT.TLorentzVector(0,0,0,0)]*5
        selJets={}
        leadJPt,trailJPt=0.,0.
        for ij in xrange(0,tree.nJet):

                #update jet energy scale/resolution
                jp4=ROOT.TLorentzVector()
                jp4.SetPtEtaPhiM(tree.Jet_pt[ij],tree.Jet_eta[ij],tree.Jet_phi[ij],0.)
                jrawsf=1./tree.Jet_jes[ij]
                jarea=tree.Jet_area[ij]
                genjpt=tree.Jet_genpt[ij]

                #update JES+JER for this jet
                jesLo, jesNom, jesHi = getJetEnergyScales(jp4.Pt(), eta=jp4.Eta(), rawsf=jrawsf,area=jarea,rho=tree.ttbar_rho)
                jerLo, jerNom, jerHi = getJetResolutionScales(pt=jesNom*jp4.Pt(), eta=jp4.Eta(), genpt=genjpt)
                oldjp4=ROOT.TLorentzVector(jp4)
                jp4 = jp4*jesNom*jerNom

                #cross clean first, as precaution
                minDRlj=9999.
                for il in xrange(0,2) : minDRlj = ROOT.TMath.Min( minDRlj, lp4[il].DeltaR(jp4) )
                if minDRlj<0.4 : continue

                if jp4.Pt()>leadJPt:
                    trailJPt=leadJPt
                    leadJPt=jp4.Pt()
                elif jp4.Pt()>trailJPt:
                    trailJPt=jp4.Pt()

                #apply energy shifts according to systematic variation
                canBeSelected=False
                selJets[ij]=[]
                for iSystVar in xrange(0,5):
                    varjp4=ROOT.TLorentzVector(jp4)
                    if iSystVar==1 : varjp4 *= jesLo/jesNom
                    if iSystVar==2 : varjp4 *= jesHi/jesNom
                    if iSystVar==3 : varjp4 *= jerLo/jerNom
                    if iSystVar==4 : varjp4 *= jerHi/jerNom

                    jetDiff[iSystVar] += (varjp4-oldjp4)
                    selJets[ij].append(varjp4)

                    #select in fiducial region
                    if varjp4.Pt()<30 or math.fabs(varjp4.Eta())>2.5 : continue
                    canBeSelected=True
                    jetCount[iSystVar]+=1                

                #if never in kinematic range forget this jet
                if not canBeSelected : del selJets[ij]

        #base selection
        zCand   = True if 'll' in ch and ROOT.TMath.Abs(mll-91)<15 else False        
        passMet = True if 'emu' in ch or tree.ttbar_metpt>40 else False
        passJets = True if len(selJets)>=2 else False
        
        #n-1 plots
        if passMet and passJets : 
            histos[ch+'_mllinc'].Fill(mll,evWgt)
            if zCand  : histos[ch+'_evsel'].Fill(0,evWgt)
        if not zCand and passJets: 
            histos[ch+'_met'].Fill(tree.ttbar_metpt,evWgt)
        #select event
        if zCand or not passJets or not passMet: continue

        #plots in control region
        histos[ch+'_evsel'].Fill(1,evWgt)
        if len(selJets)-2==0:   histos[ch+'_evsel'].Fill(2,evWgt)
        elif len(selJets)-2==1: histos[ch+'_evsel'].Fill(3,evWgt)
        else:                   histos[ch+'_evsel'].Fill(4,evWgt)

        histos[ch+'_rho'].Fill(tree.ttbar_rho,evWgt) 
        histos[ch+'_npv'].Fill(tree.nPV-1,evWgt)       
        histos[ch+'_njets'].Fill(len(selJets),evWgt)
        histos[ch+'_leadjpt'].Fill(leadJPt,evWgt)
        histos[ch+'_leadlpt'].Fill(lp4[0].Pt(),evWgt)
        histos[ch+'_trailjpt'].Fill(trailJPt,evWgt)
        histos[ch+'_traillpt'].Fill(lp4[1].Pt(),evWgt)
        histos[ch+'_mll'].Fill(mll,evWgt)

        #update weights for systematics
        for iSystVar in xrange(0,len(jetCount)):
            selWeight=1 if jetCount[iSystVar]>=2 else 0
            kinVars['weight'][0][iSystVar]=evWgt*selWeight
        kinVars['weight'][0][5] = evWgt*puWgtLo/puWgtNom
        kinVars['weight'][0][6] = evWgt*puWgtHi/puWgtNom
        kinVars['weight'][0][7] = evWgt*trigWgtLo/trigWgtNom
        kinVars['weight'][0][8] = evWgt*trigWgtHi/trigWgtNom
        kinVars['weight'][0][9] = evWgt*lepSelEffLo/lepSelEffNom
        kinVars['weight'][0][10]= evWgt*lepSelEffHi/lepSelEffNom
        kinVars['weight'][0][11]= evWgt*qcdScaleLo/genWgt
        kinVars['weight'][0][12]= evWgt*qcdScaleHi/genWgt
        kinVars['weight'][0][13]= evWgt*hdampLo/genWgt
        kinVars['weight'][0][14]= evWgt*hdampHi/genWgt

        #save jet tree for later analysis
        for ij in selJets:

            #jet flavour
            flav=0
            if ROOT.TMath.Abs(tree.Jet_flavour[ij]) in [1,2,3,21]: flav=1
            if ROOT.TMath.Abs(tree.Jet_flavour[ij])==4 : flav=2
            if ROOT.TMath.Abs(tree.Jet_flavour[ij])==5 : flav=3
            kinVars['jeteta'][0][0]     = ROOT.TMath.Abs(tree.Jet_eta[ij])
            kinVars['flavour'][0][0]    = flav

            #value of b-taggers for this jet
            for tag,tagDef in taggersList :
                tagVal=getattr(tree,tagDef[1])[ij]
                kinVars[tag+'Tagger'][0][0]=tagVal
                histos['%s_%sTagger' % (ch,tag) ].Fill(tagVal,evWgt)

            for iSystVar in xrange(0,len(selJets[ij])):

                #prepare variables for MVA
                ljkin=[]
                jp4=selJets[ij][iSystVar]
                if jp4.Pt()<30 or math.fabs(jp4.Eta())>2.5 : continue
                for il in xrange(0,2) : 
                    dr=lp4[il].DeltaR(jp4)
                    dphi=ROOT.TMath.Abs(lp4[il].DeltaPhi(jp4))
                    deta=ROOT.TMath.Abs(lp4[il].Eta()-jp4.Eta())
                    ptrel=ROOT.Math.VectorUtil.Perp(lp4[il].Vect(),jp4.Vect().Unit())/lp4[il].P();
                    ljkin.append( ( dr, (lp4[il]+jp4).M(), ptrel, dphi, deta ) )                
                ljkin=sorted(ljkin, key=lambda pair: pair[0])

                #evaluate MVA
                tmvaVars['close_mlj[0]'][0]        = ljkin[0][1]
                tmvaVars['close_ptrel'][0]         = ljkin[0][2]
                tmvaVars['close_dphi'][0]          = ljkin[0][3]
                tmvaVars['close_deta'][0]          = ljkin[0][4]
                tmvaVars['far_mlj'][0]             = ljkin[1][1]
                tmvaVars['far_ptrel'][0]           = ljkin[1][2]
                tmvaVars['far_dphi'][0]            = ljkin[1][3]
                tmvaVars['far_deta'][0]            = ljkin[1][4]                
                kinVars['close_mlj'][0][iSystVar]  = tmvaVars['close_mlj[0]'][0]
                kinVars['kindisc'][0][iSystVar]    = tmvaReader.EvaluateMVA('BDT') if not tmvaReader is None else 0
                kinVars['jetpt'][0][iSystVar]      = jp4.Pt()
                
                if iSystVar>0 : continue
                kinVars['close_ptrel'][0][0]   = tmvaVars['close_ptrel'][0]
                kinVars['close_dphi'][0][0]    = tmvaVars['close_dphi'][0]
                kinVars['close_deta'][0][0]    = tmvaVars['close_deta'][0]
                kinVars['far_mlj'][0][0]       = tmvaVars['far_mlj'][0]
                kinVars['far_ptrel'][0][0]     = tmvaVars['far_ptrel'][0]
                kinVars['far_dphi'][0][0]      = tmvaVars['far_dphi'][0]
                kinVars['far_deta'][0][0]      = tmvaVars['far_deta'][0]
              
            #fill histos, when available
            for v in kinVars :
                try:
                    histos[ch+'_'+v].Fill(kinVars[v][0][0],evWgt)
                except:
                    pass

            #fill tree
            kinTree.Fill()

    #all done with this file
    inF.Close()
    
    #dump results to file
    outF.cd()
    for h in histos: histos[h].Write()
    kinTree.Write()                            
    outF.Close()


"""
Wrapper to be used when run in parallel
"""
def runTTbarAnalysisPacked(args):
    inFile, outFile, wgt, taggers, tmvaWgts = args
    try:
        return runTTbarAnalysis(inFile=inFile, outFile=outFile, wgt=wgt, taggers=taggers, tmvaWgts=tmvaWgts)
    except :
        print 50*'<'
        print "  Problem  (%s) with %s continuing without"%(sys.exc_info()[1],inFile)
        print 50*'<'
        return False


"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-t', '--taggers',     dest='taggers'  ,   help='json with list of taggers',    default=None,        type='string')
    parser.add_option('-j', '--json',        dest='json'  ,      help='json with list of files',      default=None,        type='string')
    parser.add_option('-i', '--inDir',       dest='inDir',       help='input directory with files',   default=None,        type='string')
    parser.add_option('-o', '--outDir',      dest='outDir',      help='output directory',             default='analysis',  type='string')
    parser.add_option(      '--only',        dest='only',        help='process only matching (csv)',  default='',          type='string')
    parser.add_option(      '--tmvaWgts',    dest='tmvaWgts',    help='tmva weights',             default=None,  type='string')
    parser.add_option('-n', '--njobs',       dest='njobs',       help='# jobs to run in parallel',    default=0,           type='int')
    (opt, args) = parser.parse_args()

    
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
            task_list.append( (input_list[nf],outF,wgt, opt.taggers, opt.tmvaWgts) )

    task_list=list(set(task_list))
    print '%s jobs to run in %d parallel threads' % (len(task_list), opt.njobs)

    #run the analysis jobs
    if opt.njobs == 0:
        for inFile, outFile,wgt, taggers, tmvaWgts in task_list: 
            runTTbarAnalysis(inFile=inFile, outFile=outFile, wgt=wgt, taggers=opt.taggers, tmvaWgts=tmvaWgts)
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
