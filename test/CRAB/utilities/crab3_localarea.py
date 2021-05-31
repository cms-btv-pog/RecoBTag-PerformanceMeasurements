#!/usr/bin/env python
"""
 rename unnecessary data/ directories to data_OLD/
 (so that they are not included into the crab3 tarball
  and the latter does not exceed its max-allowed size)
"""
from __future__ import print_function

import argparse, os

def EXE(cmd, verbose=False):
    if verbose: print('> '+cmd)

    ret = os.system(cmd)
    if ret: raise SystemExit(ext)

#### main
if __name__ == '__main__':
    ### args
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('--mask-data', dest='mask_data', action='store_true', default=False,
                        help='rename data/ sub-directories to data_OLD/ (such that crab3 ignores them)')

    parser.add_argument('--restore-data', dest='restore_data', action='store_true', default=False,
                        help='rename data_OLD/ sub-directories, if any, back to data/')

    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False,
                        help='enable verbose mode')

    opts, opts_unknown = parser.parse_known_args()
    ###

    INPUT_DIRS = [

      'TopAnalysis/Configuration/analysis/common/data',
      'TopAnalysis/Configuration/analysis/diLeptonic/data',
      'TopAnalysis/Configuration/analysis/ttBSM/data',
      'TopAnalysis/Configuration/analysis/ttH/data',
      'TopAnalysis/Configuration/analysis/ttbb/data',
#      'TopAnalysis/ZTopUtils/data',

      'TTH/CommonClassifier/data',
      'TTH/CommonClassifier/crab',

      'TTH/MEIntegratorStandalone/data',
      'TTH/RecoLikelihoodReconstruction/data',
    ]

    if opts.mask_data or opts.restore_data:

       if 'CMSSW_BASE' not in os.environ:
          raise SystemExit(1)

       INPUT_DIRS = [os.path.abspath(os.environ['CMSSW_BASE']+'/src/'+_tmp) for _tmp in list(set(INPUT_DIRS))]

       for i_dir in INPUT_DIRS:

           inp_dir = i_dir
           out_dir = i_dir+'_OLD'

           if opts.restore_data: inp_dir, out_dir = out_dir, inp_dir

           if not os.path.isdir(inp_dir):
              print(' > target input directory does not exist: '+inp_dir)
              continue

           if os.path.exists(out_dir):
              print(' > target output directory already exists: '+out_dir)
              continue

           EXE('mv '+inp_dir+' '+out_dir, verbose=opts.verbose)
