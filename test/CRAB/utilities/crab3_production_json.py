#!/usr/bin/env python
import os, argparse, json, copy

from common import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='-- create production JSON for crab3 tasks --')
    parser.add_argument('-d', '--datasets-jsons', dest='datasets_jsons',
                        nargs='+', type=str, required=True,
                        help='paths to .json files containing input datasets')
    parser.add_argument('-p', '--production-json', dest='production_json',
                        action='store', default='', required=True,
                        help='path to output production .json file')
    parser.add_argument('-o', '--output-dir', dest='output_dir',
                        action='store', default='', required=True,
                        help='path to target directory for the final output NTuples')
    parser.add_argument('-s', '--storage-dir', dest='storage_dir',
                        action='store', default='', required=True,
                        help='path to storage area for crab3 outputs (Tier-2)')
    parser.add_argument('-c', '--cfg-cmsRun', dest='cfg_cmsRun',
                        action='store', default='ntuple_cfg.py', required=True,
                        help='path to cmsRun configuration file (pset) to be executed via crab3')
    parser.add_argument('-a', '--automatic', dest='automatic',
                        action='store_true', default=False,
                        help='enable crab3\'s automatic job-splitting')
    parser.add_argument('-v', '--verbose', dest='verbose',
                        action='store_true', default=False,
                        help='show verbose printouts during execution')
    opts, opts_unknown = parser.parse_known_args()
    ###

    log_prx = os.path.basename(__file__)+' -- '

    # production_json
    if os.path.exists(opts.production_json):
       KILL(log_prx+'target path to output production .json file already exists [-p]: '+opts.production_json)

    # storage_dir
    if '/store/' not in opts.storage_dir:
       KILL(log_prx+'invalid path to storage area for crab3 outputs (Tier-2), does not contain "/store/" [-s]: '+opts.storage_dir)

    elif not opts.storage_dir.startswith('/store/'):

       opts.storage_dir = opts.storage_dir[opts.storage_dir.find('/store/'):]

       if opts.verbose: WARNING(log_prx+'modified path to storage area (Tier-2) to '+opts.storage_dir)

    # cfg_cmsRun
    if not os.path.isfile(opts.cfg_cmsRun):
       KILL(log_prx+'invalid path to cmsRun configuration file (pset) to be executed via crab3 [-c]: '+opts.cfg_cmsRun)

    # datasets_jsons
    if len(opts.datasets_jsons) == 0:
       KILL(log_prx+'empty list of paths to .json files containing input datasets [-d]')

    ### construction of output json
    datasets_json = {}

    for i_dset_json in opts.datasets_jsons:
        if not os.path.isfile(i_dset_json):
           KILL(log_prx+'invalid path to .json file containing input datasets [-d]: '+i_dset_json)

        i_dset_dict = json.load(open(i_dset_json))

        for i_dset_key in i_dset_dict:

            # copy and update dataset dictionary
            datasets_json[i_dset_key] = copy.deepcopy(i_dset_dict[i_dset_key])

            datasets_json[i_dset_key]['OutputPrePath'] = os.path.abspath(opts.output_dir)+'/'+i_dset_key

            if 'crab3' in datasets_json[i_dset_key]:
               datasets_json[i_dset_key]['crab3']['Data.outLFNDirBase'] = os.path.abspath(opts.storage_dir)+'/'+i_dset_key
               datasets_json[i_dset_key]['crab3']['JobType.psetName']   = os.path.abspath(opts.cfg_cmsRun)

               if 'Data.lumiMask' in datasets_json[i_dset_key]['crab3']:
                  datasets_json[i_dset_key]['crab3']['Data.lumiMask'] = os.path.abspath(os.path.expandvars(datasets_json[i_dset_key]['crab3']['Data.lumiMask']))

                  if not os.path.isfile(datasets_json[i_dset_key]['crab3']['Data.lumiMask']):
                     KILL(log_prx+'argument of "Data.lumiMask" is not a valid file [dataset-key="'+i_dset_key+'"]')

               if opts.automatic:
                  datasets_json[i_dset_key]['crab3']['Data.splitting'] = 'Automatic'

                  if opts.verbose: WARNING(log_prx+'changed crab3 configuration parameter "Data.splitting" to "Automatic" [key='+i_dset_key+']')

            elif opts.verbose:
               WARNING(log_prx+'input dictionary for "'+i_dset_key+'" does not contain a key named "crab3" (where are the crab3 configuration parameters?): file='+i_dset_json)
    ### ---------------------------

    # output production json
    json.dump(datasets_json, open(opts.production_json, 'w'), indent=2, sort_keys=True)

    # disable write permissions
    EXE('chmod 444 '+opts.production_json)

    print colored_text('[output]', ['1', '92']), opts.production_json
