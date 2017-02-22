brilcalc trg --prescale --hltpath "HLT_BTagMu_DiJet20_Mu5_v*"  | grep Weights/Prescales/HLT_BTagMu > BTagMu_DiJet20_Mu5.txt
brilcalc trg --prescale --hltpath "HLT_BTagMu_DiJet40_Mu5_v*"  | grep HLT_BTagMu > BTagMu_DiJet40_Mu5.txt
brilcalc trg --prescale --hltpath "HLT_BTagMu_DiJet70_Mu5_v*"  | grep HLT_BTagMu > BTagMu_DiJet70_Mu5.txt
brilcalc trg --prescale --hltpath "HLT_BTagMu_DiJet110_Mu5_v*" | grep HLT_BTagMu > BTagMu_DiJet110_Mu5.txt
brilcalc trg --prescale --hltpath "HLT_BTagMu_Jet300_Mu5_v*"   | grep HLT_BTagMu > BTagMu_Jet300_Mu5.txt

awk '{print $2 "\t"  $4 "\t" $8}' BTagMu_DiJet20_Mu5.txt  > Prescales_DiJet20.txt
awk '{print $2 "\t"  $4 "\t" $8}' BTagMu_DiJet40_Mu5.txt  > Prescales_DiJet40.txt
awk '{print $2 "\t"  $4 "\t" $8}' BTagMu_DiJet70_Mu5.txt  > Prescales_DiJet70.txt
awk '{print $2 "\t"  $4 "\t" $8}' BTagMu_DiJet110_Mu5.txt > Prescales_DiJet110.txt
awk '{print $2 "\t"  $4 "\t" $8}' BTagMu_Jet300_Mu5.txt   > Prescales_Jet300.txt

rm BTagMu*_Mu5.txt 

brilcalc trg --prescale --hltpath "HLT_PFJet40_v*"  | grep HLT_PFJet > JetHT_PFJet40.txt
brilcalc trg --prescale --hltpath "HLT_PFJet60_v*"  | grep HLT_PFJet > JetHT_PFJet60.txt
brilcalc trg --prescale --hltpath "HLT_PFJet80_v*"  | grep HLT_PFJet > JetHT_PFJet80.txt
brilcalc trg --prescale --hltpath "HLT_PFJet260_v*" | grep HLT_PFJet > JetHT_PFJet260.txt

awk '{print $2 "\t"  $4 "\t" $8}' JetHT_PFJet40.txt  > Prescales_PFJet40.txt
awk '{print $2 "\t"  $4 "\t" $8}' JetHT_PFJet60.txt  > Prescales_PFJet60.txt
awk '{print $2 "\t"  $4 "\t" $8}' JetHT_PFJet80.txt  > Prescales_PFJet80.txt
awk '{print $2 "\t"  $4 "\t" $8}' JetHT_PFJet260.txt > Prescales_PFJet260.txt

rm JetHT_PFJet*.txt 
