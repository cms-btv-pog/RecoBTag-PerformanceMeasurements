to calculate the bb scale factors

1) submit jobs -> ./crabConfig.sh
2) make pu weights -> root pu.c+
3) make skims -> root makeskims.c+
4) combine SingleMuon skims -> hadd ./skims/SingleMuon.root ./skims/SingleMuon*.root
5) make data/mc comparisons -> root datamccomparisons.c+
6) make scale factors -> root sf.c+
7) combine the different wp sf to one file -> root sfaddfiles.root+

see results at:
https://twiki.cern.ch/CMS/TTbar_bbsf
