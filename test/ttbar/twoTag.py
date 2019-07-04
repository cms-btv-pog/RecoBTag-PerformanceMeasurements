import ROOT
import sys
import math

fileName=sys.argv[1]
tFile=ROOT.TFile.Open(fileName)

lowerPtBinEdges=[20,30,50,70,100,200,300]

def ReturnPtLabel(iPT):
    ptBinLabel=""
    if(iPT==-1):
        ptBinLabel="Inclusive"
    else:
        ptBinLabel+=str(lowerPtBinEdges[iPT])+"to"
        if(iPT==len(lowerPtBinEdges)-1):
            ptBinLabel+="Inf"
        else:
            ptBinLabel+=str(lowerPtBinEdges[iPT+1])
    return ptBinLabel


tags=["twoTags_deepCSVL","only2_twoTags_deepCSVL","twoTags_deepCSVM","only2_twoTags_deepCSVM","twoTags_deepCSVT","only2_twoTags_deepCSVT","twoTags_deepFlavourL","only2_twoTags_deepFlavourL","twoTags_deepFlavourM","only2_twoTags_deepFlavourM","twoTags_deepFlavourT","only2_twoTags_deepFlavourT"]
tags.sort()
samples=["t#bar{t} SL","VV","tW","t#bar{t} DL"]
samples=["t#bar{t} SL","t#bar{t} DL","t#bar{t} AH","ZZ","WZ","WW","tW","t#bar{t} DL","DY"]
systs=["","_fsrConHi","_fsrConLo","_fsrDefHi","_fsrDefLo","_fsrRedHi","_fsrRedLo","_isrConHi","_isrConLo","_isrDefHi","_isrDefLo","_isrRedHi","_isrRedLo","_qcdScaleHi","_qcdScaleLo"]
systs=[""]


statError={}


#can=ROOT.TCanvas("can","",700,700)
results={}
tgraphs={}

can=ROOT.TCanvas("can","",800,800)
outputFile=ROOT.TFile("twoTagSFGraphs.root","RECREATE")


for syst in systs:
    #print syst
    for tag in tags:
        print tag,
        iGraph=0
        if syst=="":
            tgraphs[tag]=ROOT.TGraphErrors()
            tgraphs[tag].SetTitle(tag+";Jet P_{T};Scale Factor")
        if tag.find("twoTags_") == 0:
            print "     ",
        for ptBin in range(-1,len(lowerPtBinEdges)):
            ptLabel=ReturnPtLabel(ptBin)
            basename="emu_"+tag+syst+ptLabel
            print basename
            for sample in samples:
                histName=basename+"/"+basename+"_"+sample
                print histName
                hist=tFile.Get(histName)
                if sample==samples[0]:
                    histBkg=hist.Clone("histBkg")
                else:
                    histBkg.Add(hist)
                #print histName
                #for iBin in range(15):
                #    print iBin,hist.GetBinContent(iBin)
            
            tot={}
            tot["2bmc"]=0
            tot["othermc"]=0
            tot["data"]=0
            tot["2bmc_2btag"]=0
            tot["othermc_2btag"]=0
            tot["data_2btag"]=0
            for iBin in range(15):
                #print iBin,histBkg.GetBinContent(iBin)
                if iBin%4==1:
                    tot["2bmc"]=tot["2bmc"]+histBkg.GetBinContent(iBin)
                    if iBin>9:
                        tot["2bmc_2btag"]=tot["2bmc_2btag"]+histBkg.GetBinContent(iBin)
                else:
                    tot["othermc"]=tot["othermc"]+histBkg.GetBinContent(iBin)
                    if iBin>9:
                        tot["othermc_2btag"]=tot["othermc_2btag"]+histBkg.GetBinContent(iBin)
            twoTagDataHist=basename+"/"+basename
            hist=tFile.Get(twoTagDataHist)
            for iBin in range(15):
                #print iBin,hist.GetBinContent(iBin),histBkg.GetBinContent(iBin)
                tot["data"]=tot["data"]+hist.GetBinContent(iBin)
                if iBin>9:
                    tot["data_2btag"]=tot["data_2btag"]+hist.GetBinContent(iBin)
            
            print tot
            
            print histBkg.Integral(),hist.Integral(),
            if histBkg.Integral()>0:
                hist.Integral()/histBkg.Integral()
            else:
                print "division by zero avoided"
                continue
            
            for name in tot:
                if name.find("mc")!=-1:
                    tot[name]=tot[name]/histBkg.Integral()
                else:
                    tot[name]=tot[name]/hist.Integral()
                    
            
            print tot
            try:
                eff=math.sqrt((tot["data_2btag"]-tot["othermc_2btag"])/tot["2bmc"])
                effmc=math.sqrt(tot["2bmc_2btag"]/tot["2bmc"])
                
                effUp=math.sqrt((tot["data_2btag"]-tot["othermc_2btag"]*1.5)/tot["2bmc"])
                effDown=math.sqrt((tot["data_2btag"]-tot["othermc_2btag"]*0.5)/tot["2bmc"])
                results[tag+syst+ptLabel]=[tot["2bmc_2btag"],tot["2bmc"],tot["data_2btag"],tot["othermc_2btag"],eff,effmc,eff/effmc,effUp/effmc,effDown/effmc]
                
                if syst=="" and ptBin>-1:
                    pt1=lowerPtBinEdges[ptBin]
                    pt2=1000
                    if ptBin!=lowerPtBinEdges[-1]:
                        pt2=lowerPtBinEdges[ptBin+1]
                    tgraphs[tag].SetPoint(iGraph,(pt1+pt2)/2, eff/effmc)
                    tgraphs[tag].SetPointError(iGraph,(pt2-pt1)/2, abs(effUp-effDown)/effmc)
                    iGraph=iGraph+1
                    
            except:
                print "couldn't do this one"
                continue
 

            
            if syst=="":
                print hist.Integral(),histBkg.Integral(),histBkg.GetEffectiveEntries()
                statError[tag+"data_2btag"] = math.sqrt( tot["data_2btag"] * (1 - tot["data_2btag"]) / hist.Integral())
                statError[tag+"othermc_2btag"] = math.sqrt( tot["othermc_2btag"] * (1 - tot["othermc_2btag"]) / histBkg.GetEffectiveEntries())
                statError[tag+"2bmc"] = math.sqrt( tot["2bmc"] * (1 - tot["2bmc"]) /  histBkg.GetEffectiveEntries())
                statError[tag+"eff"] = math.sqrt(math.pow(eff*statError[tag+"2bmc"],2)  + (math.pow(statError[tag+"data_2btag"],2)+math.pow(statError[tag+"othermc_2btag"],2)  ) / math.pow(eff,2)    )/(2*tot["2bmc"])
        tgraphs[tag].Draw("APE")
        can.Update()
        can.SaveAs(tag+".png")
        tgraphs[tag].SetName(tag)
        outputFile.cd()
        tgraphs[tag].Write()

outputFile.Close()
print "tag,",
for syst in systs:
    if syst=="":
        print "data, staterror, nonbbUP, nonbbDOWN,",
    else:
        print syst+",",
print

sortedResults=results.keys()
sortedResults.sort()
for result in sortedResults:
    print result,
    for syst in systs:
        if syst=="":
            print results[result][6],",",results[result][7],",",results[result][8],",",
            #print results[tag+syst+ptLabel][6],",",statError[tag+"eff"],",",results[tag+syst+ptLabel][7],",",results[tag+syst+ptLabel][8],",",
        else:
            print results[result][6],",",
    print
