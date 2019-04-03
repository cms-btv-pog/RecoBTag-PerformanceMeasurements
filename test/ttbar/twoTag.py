import ROOT
import sys
import math

fileName=sys.argv[1]
tFile=ROOT.TFile.Open(fileName)


tags=["twoTags_deepCSVL","only2_twoTags_deepCSVL","twoTags_deepCSVM","only2_twoTags_deepCSVM","twoTags_deepCSVT","only2_twoTags_deepCSVT","twoTags_deepFlavourL","only2_twoTags_deepFlavourL","twoTags_deepFlavourM","only2_twoTags_deepFlavourM","twoTags_deepFlavourT","only2_twoTags_deepFlavourT","twoTags_CSVv2L","only2_twoTags_CSVv2L","twoTags_CSVv2M","only2_twoTags_CSVv2M","twoTags_CSVv2T","only2_twoTags_CSVv2T"]
tags.sort()
samples=["t#bar{t} SL","VV","tW","t#bar{t} DL"]
samples=["t#bar{t} SL","t#bar{t} DL","t#bar{t} AH","ZZ","WZ","WW","tW","t#bar{t} DL","DY"]
systs=[""]
systs=["","_fsrConHi","_fsrConLo","_fsrDefHi","_fsrDefLo","_fsrRedHi","_fsrRedLo","_isrConHi","_isrConLo","_isrDefHi","_isrDefLo","_isrRedHi","_isrRedLo","_qcdScaleHi","_qcdScaleLo"]

# KEY: TH1F  emu_only2_twoTags_CSVv2L_t#bar{t} AH;1  t#bar{t} AH
# KEY: TH1F  emu_only2_twoTags_CSVv2L_other;1  other
# KEY: TH1F  emu_only2_twoTags_CSVv2L_t#bar{t} SL;1  t#bar{t} SL
# KEY: TH1F  emu_only2_twoTags_CSVv2L_tW;1 tW
# KEY: TH1F  emu_only2_twoTags_CSVv2L_t#bar{t} DL;1  t#bar{t} DL


statError={}
#samples=["emu_only2_twoTags/emu_only2_twoTags_t#bar{t} SL","emu_only2_twoTags/emu_only2_twoTags_VV","emu_only2_twoTags/emu_only2_twoTags_tW","emu_only2_twoTags/emu_only2_twoTags_t#bar{t} DL","emu_only2_twoTags/emu_only2_twoTags_DY"]
#twoTagDataHist="emu_only2_twoTags/emu_only2_twoTags"


#can=ROOT.TCanvas("can","",700,700)
results={}

for syst in systs:
    #print syst
    for tag in tags:
        print tag,
        if tag.find("twoTags_") == 0:
            print "     ",
        for sample in samples:
            histName="emu_"+tag+syst+"/emu_"+tag+syst+"_"+sample
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
        twoTagDataHist="emu_"+tag+syst+"/emu_"+tag+syst
        hist=tFile.Get(twoTagDataHist)
        for iBin in range(15):
            #print iBin,hist.GetBinContent(iBin),histBkg.GetBinContent(iBin)
            tot["data"]=tot["data"]+hist.GetBinContent(iBin)
            if iBin>9:
                tot["data_2btag"]=tot["data_2btag"]+hist.GetBinContent(iBin)
        
        #print tot
        
        #print histBkg.Integral(),hist.Integral(),hist.Integral()/histBkg.Integral()
        
        for name in tot:
            if name.find("mc")!=-1:
                tot[name]=tot[name]/histBkg.Integral()
            else:
                tot[name]=tot[name]/hist.Integral()
                
        
        #print tot
        
        eff=math.sqrt((tot["data_2btag"]-tot["othermc_2btag"])/tot["2bmc"])
        effmc=math.sqrt(tot["2bmc_2btag"]/tot["2bmc"])
        
        effUp=math.sqrt((tot["data_2btag"]-tot["othermc_2btag"]*1.5)/tot["2bmc"])
        effDown=math.sqrt((tot["data_2btag"]-tot["othermc_2btag"]*0.5)/tot["2bmc"])
        
        #print tot["2bmc_2btag"],tot["2bmc"],tot["data_2btag"],tot["othermc_2btag"],eff,effmc,eff/effmc,effUp/effmc,effDown/effmc
        results[tag+syst]=[tot["2bmc_2btag"],tot["2bmc"],tot["data_2btag"],tot["othermc_2btag"],eff,effmc,eff/effmc,effUp/effmc,effDown/effmc]

        
        if syst=="":
            print hist.Integral(),histBkg.Integral(),histBkg.GetEffectiveEntries()
            statError[tag+"data_2btag"] = math.sqrt( tot["data_2btag"] * (1 - tot["data_2btag"]) / hist.Integral())
            statError[tag+"othermc_2btag"] = math.sqrt( tot["othermc_2btag"] * (1 - tot["othermc_2btag"]) / histBkg.GetEffectiveEntries())
            statError[tag+"2bmc"] = math.sqrt( tot["2bmc"] * (1 - tot["2bmc"]) /  histBkg.GetEffectiveEntries())
            statError[tag+"eff"] = math.sqrt(math.pow(eff*statError[tag+"2bmc"],2)  + (math.pow(statError[tag+"data_2btag"],2)+math.pow(statError[tag+"othermc_2btag"],2)  ) / math.pow(eff,2)    )/(2*tot["2bmc"])

raw_input()
print "tag,",
for syst in systs:
    if syst=="":
        print "data, staterror, nonbbUP, nonbbDOWN,",
    else:
        print syst+",",
print
for tag in tags:
    print tag+",",
    for syst in systs:
        if syst=="":
            print results[tag+syst][6],",",statError[tag+"eff"],",",results[tag+syst][7],",",results[tag+syst][8],",",
        else:
            print results[tag+syst][6],",",
    print
