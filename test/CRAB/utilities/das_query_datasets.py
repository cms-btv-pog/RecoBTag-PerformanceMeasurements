#!/usr/bin/env python
"""das_query_datasets.py: helper to list datasets on DAS"""

import argparse, os

from common import *

def das_output_lines(cmd):
    _tmp_out = get_output(cmd)[0]

    _tmp_out_ls = _tmp_out.split('\n')
    _tmp_out_ls = [l for l in _tmp_out_ls if (l != '' and not l.startswith('Showing'))]
#    _tmp_out_ls = sorted(list(set(_tmp_out_ls)))

    return _tmp_out_ls
#---

def das_output(cmd):
    _tmp0_ls = das_output_lines(cmd)

    if len(_tmp0_ls) != 1:
        print cmd, _tmp0_ls
        raise SystemExit

    return _tmp0_ls[0]
#---

#### main
if __name__ == '__main__':
    ### args
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-d', '--dataset', dest='dataset',
                        action='store', default='',
                        help='dataset for DAS query, accepts wildcards (example: /SingleMuon/*/MINIAOD)')
    parser.add_argument('-i', '--input', dest='input',
                        action='store', default='',
                        help='path to input text file with list of datasets for DAS query, accepts wildcards (example: /SingleMuon/*/MINIAOD)')
    parser.add_argument('-s', '--status', dest='status',
                        action='store', default='valid',
                        help='restrict query to datasets in status=X (e.g. \'*\', valid, production)')
    parser.add_argument('--das-client', dest='das_client',
                        action='store', default='/cvmfs/cms.cern.ch/common/dasgoclient',
                        help='path to das-client script (default: dasgoclient)')
    parser.add_argument('--das-options', dest='das_opts',
#                        action='store', default='limit=0,key=~/.globus/userkey.pem,cert=~/.globus/usercert.pem',
                        action='store', default='limit=0',
                        help='comma-sep list of DAS query options')
    opts, opts_unknown = parser.parse_known_args()
    ###

    ### --- configuration

    log_prx = os.path.basename(__file__)+' -- '

    VERBOSE = False

    DSET_INPUT_LS = None

    if opts.dataset and opts.input:
       KILL(log_prx+'logic error: trying to use both dataset [-d '+opts.dataset+'] and list of datasets via input file [-i '+opts.input+']')

    elif opts.dataset and not opts.input:
       DSET_INPUT_LS = [opts.dataset]

    elif not opts.dataset and opts.input:
       if not os.path.isfile(opts.input):
          KILL(log_prx+'target path to input file with list of datasets not found [-i]: '+opts.input)

       with open(opts.input) as ifile:
          ifile_cont = ifile.readlines()
          DSET_INPUT_LS = [i_line.strip().replace(' ','') for i_line in ifile_cont]
          DSET_INPUT_LS = [i_line for i_line in DSET_INPUT_LS if i_line != '']

    else:
       KILL(log_prx+'undefined dataset for DAS query [-d], undefined path to input file with list of datasets [-i]')

    if not os.path.isfile(opts.das_client):
       KILL(log_prx+'invalid path to "das_client" script: '+DAS_CLIENT)

    DAS_CLIENT = opts.das_client

    ### -----------------

    ### implementation --

    VOMS_PROXY_INFO = get_output('which voms-proxy-info')[0]
    if not VOMS_PROXY_INFO: KILL(log_prx+'executable "voms-proxy-info" not available')

    (voms_proxy_info_stdout, voms_proxy_info_stderr) = get_output(VOMS_PROXY_INFO, permissive=True)

    exe_voms_proxy_init = (('timeleft  : 00:00:00' in voms_proxy_info_stdout) or ('Proxy not found' in voms_proxy_info_stderr))

    if exe_voms_proxy_init:

       print '\n >>> WARNING -- '+log_prx+'attempting to run '+DAS_CLIENT+' without user credentials/X509 proxy,',
       print 'creating proxy via "voms-proxy-init -voms cms -rfc"\n'

       EXE('voms-proxy-init -voms cms -rfc')
       EXE('sleep 1')

    DAS_OPTS = ['--'+das_opt for das_opt in opts.das_opts.split(',')]

    DAS_OPTS_STR = ' '.join(DAS_OPTS)

    DSET_LS = []

    for DSET in DSET_INPUT_LS:

        DAS_QUERY = 'dataset dataset='+DSET
        if opts.status != '': DAS_QUERY += ' status='+opts.status

        if VERBOSE: print DAS_CLIENT+' '+DAS_OPTS_STR+' --query=\"'+DAS_QUERY+'\"'

        DAS_OUTPUTS = das_output_lines(DAS_CLIENT+' '+DAS_OPTS_STR+' --query=\"'+DAS_QUERY+'\"')

        # printout
        print ''
        print ' query = "'+DSET+'"'
        print ''

        for i_DAS_OUTPUT in DAS_OUTPUTS: print i_DAS_OUTPUT

        print ''
        print '-'*50

    print ''

    ### -----------------
