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
SYSTVARS   = ['','jesup','jesdn','jerup','jerdn','trigdn','trigup','seldn','selup','qcdscaledn','qcdscaleup','hdampdn','hdampup']

"""
Project trees from files to build the templates
"""
def prepareTemplates(tagger,taggerDef,var,varRange,inDir,outDir):
    
    print '...starting %s for %s'%(tagger,var)

    histos={}

    nOPs=len(taggerDef)-2

    #tag counting
    tagCounts=[0]*(nOPs-1)
    for key in ['data','mc']:
        for i in xrange(1,nOPs):
            name='%s_%s_pass%d'%(key,tagger,i)
            histos[name]=ROOT.TH1F(name,';b-tag multiplicity;Events',5,0,5)

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
        chains[key].GetEntry(0)
        curEventInfo=[chains[key].EventInfo[0],chains[key].EventInfo[1],chains[key].EventInfo[2]]
        curEventWgt=1.0

        for i in xrange(0,chains[key].GetEntries()):

            chains[key].GetEntry(i)

            #do tag counting
            iEventInfo=[chains[key].EventInfo[0],chains[key].EventInfo[1],chains[key].EventInfo[2]]
            if iEventInfo!=curEventInfo:
                for iop in xrange(0,len(tagCounts)):
                    name='%s_%s_pass%d'%(key,tagger,iop+1)
                    histos[name].Fill(tagCounts[iop],chains[key].weight[0])
                    tagCounts[iop]=0
                curEventInfo=iEventInfo
            
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

                #tag counting (only nominal)
                if len(systVar)==0:
                    for iop in xrange(2,len(taggerDef)-1):
                        if taggerVal> taggerDef[iop] : 
                            tagCounts[iop-2]=tagCounts[iop-2]+1

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

    #custom fraction fitter class
    ttFracFitter=ROOT.TTbarFracFitter();

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
Massive conversion to scale factors given a MC truth reference
"""
def convertToScaleFactor(obs,ref):

    sf={}
    for var in obs:
        sf[var]={}
        for iop in obs[var]:

            newName=obs[var][iop].GetName()+'_sf'
            sf[var][iop]=obs[var][iop].Clone(newName)
            
            xobs, yobs = ROOT.Double(0), ROOT.Double(0)                
            xexp, yexp = ROOT.Double(0), ROOT.Double(0)                
            for ip in xrange(0,obs[var][iop].GetN()):
                obs[var][iop].GetPoint(ip,xobs,yobs)
                xobsUnc=obs[var][iop].GetErrorX(ip)
                yobsUnc=obs[var][iop].GetErrorY(ip)
                
                ref[iop].GetPoint(ip,xexp,yexp)
                xexpUnc=ref[iop].GetErrorX(ip)
                yexpUnc=ref[iop].GetErrorY(ip)
                
                sfVal, sfValUnc = 0., 0.
                if yexp>0:
                    sfVal=yobs/yexp
                    sfValUnc= ROOT.TMath.Sqrt( ROOT.TMath.Power(yobsUnc*yexp,2)+ROOT.TMath.Power(yobs*yexpUnc,2) )/ROOT.TMath.Power(yexp,2)
                sf[var][iop].SetPoint(ip,xobs,sfVal)
                sf[var][iop].SetPointError(ip,xobsUnc,sfValUnc)

    return sf

"""
compares results obtained from different methods
"""
def showEfficiencyResults(fileList,outDir,lumi):

    grIncExpColl,     grExpColl={},{}            # expected inclusive, differential
    grIncObsStatColl, grIncObsTotalColl={},{}    # observed inclusive
    grObsStatColl,    grObsTotalColl={},{}       # observed differential

    for ifile in xrange(0,len(fileList)):

        f=fileList[ifile]
        cachefile = open(f,'r')
        fitInfo=pickle.load(cachefile)
        fitResults=pickle.load(cachefile)
        mcTruth=pickle.load(cachefile)
        cachefile.close()
    
        var       = fitInfo['var']
        tagger    = fitInfo['tagger']
        taggerDef = fitInfo['taggerDef']
        nOPs=len(taggerDef)-2
        grIncObsStat, grIncObsTotal = {}, {}
        grObsStat,    grObsTotal    = {}, {}
        for iop in xrange(1,nOPs):
            for islice in fitResults[iop]:

                if fitResults[0][islice][''][0]==0: continue
 
                #statistical uncertainty : unc #tag / # pre-tag
                eff={
                    ''     : fitResults[iop][islice][''][0]/fitResults[0][islice][''][0],
                    'stat' : fitResults[iop][islice][''][1]/fitResults[0][islice][''][0]
                    }
                mceff={
                    ''     : mcTruth[iop][islice][0]/mcTruth[0][islice][0],
                    'stat' : mcTruth[iop][islice][1]/mcTruth[0][islice][0]
                    }

                #alternative variations
                for systVar in fitResults[iop][islice]:
                    
                    if systVar=='' : continue
                    if systVar=='closure': continue

                    #symmetrize differences for up/dn variations
                    newEff = fitResults[iop][islice][systVar][0]/fitResults[0][islice][systVar][0]
                    key    = systVar.replace('up','')
                    key    = key.replace('dn','')
                    #newEff = newEff-eff['']
                    newEff = newEff-mceff['']
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

                    if ifile==0:
                        grIncExpColl[iop]=ROOT.TGraphErrors()
                        grIncExpColl[iop].SetName('%s_%d_inc_exp'%(var,iop))
                        grIncExpColl[iop].SetTitle('simulation')
                        grIncExpColl[iop].SetPoint(0,-4,mceff[''])
                        grIncExpColl[iop].SetPointError(0,0,mceff['stat'])
                        grIncExpColl[iop].SetPoint(1,4*len(fileList),mceff[''])
                        grIncExpColl[iop].SetPointError(1,0,mceff['stat'])

                        grExpColl[iop]=ROOT.TGraphErrors()
                        grExpColl[iop].SetName('%s_%d_exp'%(var,iop))
                        grExpColl[iop].SetTitle('simulation')

                    grIncObsStat[iop]=ROOT.TGraphErrors()
                    grIncObsStat[iop].SetName('%s_%d_inc_stat'%(var,iop))
                    grIncObsStat[iop].SetPoint(0,0+ifile*4,eff[''])
                    #grIncObsStat[iop].SetPointError(0,5,eff['stat'])
                    grIncObsStat[iop].SetPointError(0,0,eff['stat'])

                    grIncObsTotal[iop] = ROOT.TGraphErrors()
                    grIncObsTotal[iop].SetName('%s_%d_inc_total'%(var,iop))
                    grIncObsTotal[iop].SetPoint(0,0+ifile*4,eff[''])
                    grIncObsTotal[iop].SetPointError(0,5,eff['total'])
                    
                    #init differential plots also
                    grObsStat[iop]=ROOT.TGraphErrors()
                    grObsStat[iop].SetName('%s_%d_stat'%(var,iop))
                    grObsTotal[iop] = ROOT.TGraphErrors()
                    grObsTotal[iop].SetName('%s_%d_total'%(var,iop))

                else:

                    xcen=0.5*(xmax+xmin)
                    dx=(xmax-xmin)*0.5

                    np=grObsStat[iop].GetN()
                    grObsStat[iop].SetPoint(np,xcen,eff[''])
                    #grObsStat[iop].SetPointError(np,dx,eff['stat'])
                    grObsStat[iop].SetPointError(np,0,eff['stat'])
                    grObsTotal[iop].SetPoint(np,xcen,eff[''])
                    grObsTotal[iop].SetPointError(np,dx,eff['total'])

                    if ifile==0:
                        np=grExpColl[iop].GetN()
                        grExpColl[iop].SetPoint(np,xcen,mceff[''])
                        grExpColl[iop].SetPointError(np,0,mceff['stat'])
                        #grExpColl[iop].SetPoint(np+1,xmin+0.6*dx,mceff[''])
                        #grExpColl[iop].SetPointError(np+1,0,mceff['stat'])

        #add collection
        if not (var in grObsStatColl):
            grIncObsStatColl[var], grIncObsTotalColl[var] = {}, {}
            grObsStatColl[var], grObsTotalColl[var] = {}, {}

        grIncObsStatColl[var]  = grIncObsStat
        grIncObsTotalColl[var] = grIncObsTotal
        grObsStatColl[var]     = grObsStat
        grObsTotalColl[var]    = grObsTotal

    #show results in canvas
    c=ROOT.TCanvas('c','c',800,500)
    c.SetRightMargin(0)
    c.SetLeftMargin(0)
    c.SetBottomMargin(0)
    c.SetTopMargin(0)

    colors=[ROOT.kAzure+9,ROOT.kOrange-1,ROOT.kGreen-5]
    markers=[20,21,22]

    c.cd()
    clabel=ROOT.TPad('clabel','clabel',0,0.95,1.0,1.0)
    clabel.SetRightMargin(0.0)
    clabel.SetTopMargin(0.0)
    clabel.SetLeftMargin(0.0)
    clabel.SetBottomMargin(0.0)
    clabel.Draw()
    
    c.cd()
    cinc=ROOT.TPad('cinc','cinc',0,0.45,0.2,0.95)
    cinc.SetRightMargin(0.02)
    cinc.SetTopMargin(0.01)
    cinc.SetLeftMargin(0.3)
    cinc.SetBottomMargin(0.01)
    cinc.Draw()

    c.cd()
    cincsf=ROOT.TPad('cincsf','cincsf',0,0,0.2,0.45)
    cincsf.SetRightMargin(0.02)
    cincsf.SetTopMargin(0.01)
    cincsf.SetLeftMargin(0.3)
    cincsf.SetBottomMargin(0.12)
    cincsf.Draw()

    c.cd()
    cexc=ROOT.TPad('cexc','cexc',0.2,0.45,1.0,0.95)
    cexc.SetRightMargin(0.02)
    cexc.SetTopMargin(0.01)
    cexc.SetLeftMargin(0.02)
    cexc.SetBottomMargin(0.01)
    cexc.Draw()

    c.cd()
    cexcsf=ROOT.TPad('cexcsf','cexcsf',0.2,0,1.0,0.45)
    cexcsf.SetRightMargin(0.02)
    cexcsf.SetTopMargin(0.01)
    cexcsf.SetLeftMargin(0.02)
    cexcsf.SetBottomMargin(0.12)
    cexcsf.Draw()
    
    firstvar=grIncObsStatColl.keys()[0]
    for iop in grIncObsStatColl[firstvar]:

        #
        #inclusive efficiency measurements
        #
        cinc.cd()
        cinc.Clear()
        ivar=0
        drawOpt='a2'
        for var in grIncObsTotalColl:
            grIncObsTotalColl[var][iop].SetTitle(var)
            grIncObsTotalColl[var][iop].SetMarkerColor(colors[ivar]+1)
            grIncObsTotalColl[var][iop].SetLineColor(colors[ivar]+1)
            grIncObsTotalColl[var][iop].SetFillColor(colors[ivar]+1)
            grIncObsTotalColl[var][iop].SetMarkerStyle(1)
            grIncObsTotalColl[var][iop].SetFillStyle(3001)
            grIncObsTotalColl[var][iop].Draw(drawOpt)
            grIncObsTotalColl[var][iop].GetYaxis().SetRangeUser(0.12,1.5)
            grIncObsTotalColl[var][iop].GetYaxis().SetTitle('Efficiency')
            grIncObsTotalColl[var][iop].GetYaxis().SetTitleOffset(1.2)
            grIncObsTotalColl[var][iop].GetYaxis().SetTitleSize(0.1)
            grIncObsTotalColl[var][iop].GetYaxis().SetLabelSize(0.1)
            grIncObsTotalColl[var][iop].GetXaxis().SetTitleSize(0.)
            grIncObsTotalColl[var][iop].GetXaxis().SetLabelSize(0.)
            grIncObsTotalColl[var][iop].GetXaxis().SetNdivisions(0)
            drawOpt='2'
            ivar+=1

        ivar=0
        for var in grIncObsStatColl:
            grIncObsStatColl[var][iop].SetMarkerColor(colors[ivar])
            grIncObsStatColl[var][iop].SetLineColor(colors[ivar])
            grIncObsStatColl[var][iop].SetFillColor(0)
            grIncObsStatColl[var][iop].SetMarkerStyle(markers[ivar])
            grIncObsStatColl[var][iop].SetFillStyle(0)
            grIncObsStatColl[var][iop].Draw('p')
            grIncObsStatColl[var][iop].SetTitle(var)
            ivar+=1

        grIncExpColl[iop].SetMarkerColor(1)
        grIncExpColl[iop].SetMarkerStyle(1)
        grIncExpColl[iop].SetLineColor(1)
        grIncExpColl[iop].SetLineWidth(2)
        grIncExpColl[iop].Draw('c')

        #
        #inclusive scale factors
        #
        grIncSFObsStatColl=convertToScaleFactor(obs=grIncObsStatColl,ref=grIncExpColl)
        grIncSFObsTotalColl=convertToScaleFactor(obs=grIncObsTotalColl,ref=grIncExpColl)
        cincsf.cd()
        cincsf.Clear()
        ivar=0
        drawOpt='a2'
        for var in grIncSFObsTotalColl:
            grIncSFObsTotalColl[var][iop].Draw(drawOpt)
            grIncSFObsTotalColl[var][iop].GetYaxis().SetRangeUser(0.74,1.26)
            grIncSFObsTotalColl[var][iop].GetYaxis().SetTitle('Data-MC scale factor')
            grIncSFObsTotalColl[var][iop].GetYaxis().SetTitleOffset(1.2)
            grIncSFObsTotalColl[var][iop].GetYaxis().SetTitleSize(0.1)
            grIncSFObsTotalColl[var][iop].GetYaxis().SetLabelSize(0.1)
            grIncSFObsTotalColl[var][iop].GetXaxis().SetTitleSize(0.)
            grIncSFObsTotalColl[var][iop].GetXaxis().SetLabelSize(0.)
            grIncSFObsTotalColl[var][iop].GetXaxis().SetNdivisions(0)
            drawOpt='2'
            ivar+=1

        ivar=0
        for var in grIncSFObsStatColl:
            grIncSFObsStatColl[var][iop].Draw('p')
            grIncSFObsStatColl[var][iop].SetTitle(var)
            ivar+=1

        #
        # differential efficiency measurements
        #
        cexc.cd()
        cexc.Clear()
        leg = ROOT.TLegend(0.6, 0.95,0.95,0.8)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(43)
        leg.SetTextSize(16)
        if len(fileList)>2 : leg.SetNColumns(2)
        leg.AddEntry(grExpColl[iop],grExpColl[iop].GetTitle(),"l")

        ivar=0
        drawOpt='a2'
        for var in grObsTotalColl:
            grObsTotalColl[var][iop].SetTitle(var)
            grObsTotalColl[var][iop].SetMarkerColor(colors[ivar]+1)
            grObsTotalColl[var][iop].SetLineColor(colors[ivar]+1)
            grObsTotalColl[var][iop].SetFillColor(colors[ivar]+1)
            grObsTotalColl[var][iop].SetMarkerStyle(1)
            grObsTotalColl[var][iop].SetFillStyle(3001)
            grObsTotalColl[var][iop].Draw(drawOpt)
            grObsTotalColl[var][iop].GetYaxis().SetTitleSize(0.)
            grObsTotalColl[var][iop].GetYaxis().SetLabelSize(0.)
            grObsTotalColl[var][iop].GetYaxis().SetRangeUser(0.0,1.5)
            drawOpt='2'
            ivar+=1

        ivar=0
        for var in grObsStatColl:
            grObsStatColl[var][iop].SetMarkerColor(colors[ivar])
            grObsStatColl[var][iop].SetLineColor(colors[ivar])
            grObsStatColl[var][iop].SetFillColor(0)
            grObsStatColl[var][iop].SetMarkerStyle(markers[ivar])
            grObsStatColl[var][iop].SetFillStyle(0)
            grObsStatColl[var][iop].Draw('p')
            grObsStatColl[var][iop].SetTitle(var)
            leg.AddEntry(grObsStatColl[var][iop],grObsStatColl[var][iop].GetTitle(),'p')
            ivar+=1

        grExpColl[iop].SetMarkerColor(1)
        grExpColl[iop].SetMarkerStyle(1)
        grExpColl[iop].SetLineColor(1)
        grExpColl[iop].SetLineWidth(2)
        grExpColl[iop].Draw('c')

        leg.Draw()


        #
        # differential scale factor measurements
        #
        grSFObsStatColl=convertToScaleFactor(obs=grObsStatColl,ref=grExpColl)
        grSFObsTotalColl=convertToScaleFactor(obs=grObsTotalColl,ref=grExpColl)
        cexcsf.cd()
        cexcsf.Clear()
        ivar=0
        drawOpt='a2'
        for var in grSFObsTotalColl:
            grSFObsTotalColl[var][iop].Draw(drawOpt)
            grSFObsTotalColl[var][iop].GetYaxis().SetTitleSize(0.)
            grSFObsTotalColl[var][iop].GetYaxis().SetLabelSize(0.)
            grSFObsTotalColl[var][iop].GetYaxis().SetRangeUser(0.74,1.26)
            grSFObsTotalColl[var][iop].GetXaxis().SetTitle(SLICEVAR)
            grSFObsTotalColl[var][iop].GetXaxis().SetTitleSize(0.08)
            grSFObsTotalColl[var][iop].GetXaxis().SetLabelSize(0.08)
            drawOpt='2'
            ivar+=1

        ivar=0
        for var in grSFObsStatColl:
            grSFObsStatColl[var][iop].Draw('p')
            grSFObsStatColl[var][iop].SetTitle(var)
            ivar+=1

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
        txt.DrawLatex(0.85,0.5,'[%s > %3.3f]' % (tagger,taggerDef[iop+1]))

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
