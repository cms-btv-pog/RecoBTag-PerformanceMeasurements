#Loop over all branches in the two output trees and get Mean and RMS. Can be used to compare it to another (previous) output file to see impact of changes
# Use it like: python compare_output.py <file1.root> -b
# Use it like: python compare_output.py <file2.root> -b
# Check diff with: diff <file1.log> <file2.log>

from ROOT import TFile, TTree, TH1

import sys

def compare_output():
  fileName = sys.argv[1]
  outfileName = fileName[:-4] + "log"
  output = open(outfileName, "w")
  file = TFile(fileName, 'READ')
  tree = file.Get('btagana/ttree')
  if(tree):
    for branch in tree.GetListOfBranches():
      branchname = branch.GetName()
      tree.Draw(branchname + '>>hist')
      mean = file.Get('hist').GetMean()
      rms = file.Get('hist').GetRMS()
      output.write('{}: {} {}\n'.format(branchname, mean, rms))

  tree = file.Get('btaganaFatJets/ttree')
  if(tree):
    for branch in tree.GetListOfBranches():
      branchname = branch.GetName()
      tree.Draw(branchname + '>>hist')
      mean = file.Get('hist').GetMean()
      rms = file.Get('hist').GetRMS()
      output.write('{}: {} {}\n'.format(branchname, mean, rms))

  output.close()
  print('Output written to {}'.format(outfileName))


if __name__ == '__main__':
  compare_output()

