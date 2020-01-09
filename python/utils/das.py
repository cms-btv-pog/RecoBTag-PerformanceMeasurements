#!/usr/bin/env python
"""
utilities related to the CMS Data Aggregation System (DAS)
"""
from __future__ import print_function

import json

from RecoBTag.PerformanceMeasurements.utils.common import *

def load_dataset_data(das_name, max_files=-1, max_events=-1, verbose=False):

    if verbose:
       print(das_name)

    dataset_split = das_name.split('/')
    if len(dataset_split) != 4:
       KILL('load_dataset_data -- invalid data-set name (format is incorrect, check slashes): '+das_name)

    dset_data = {'DAS': str(das_name), 'files': []}

    dataset_files = command_output_lines('dasgoclient --query "file dataset='+str(das_name)+'"')
    dataset_files = [_tmp for _tmp in dataset_files if _tmp != '']
    dataset_files = sorted(list(set(dataset_files)))

    if len(dataset_files) == 0:
       KILL('load_dataset_data -- empty list of input files for dataset: '+str(das_name))

    if max_files > 0:
       dataset_files = dataset_files[:max_files]

    totEvents, breakLoop = 0, False

    for i_file_idx, i_file in enumerate(dataset_files):

        if verbose:
           print('  file', i_file_idx, '/', len(dataset_files))

        i_file_nevents = command_output_lines('dasgoclient --query "file='+str(i_file)+' | grep file.nevents"')
        i_file_nevents = [_tmp.replace(' ', '') for _tmp in i_file_nevents]
        i_file_nevents = [_tmp for _tmp in i_file_nevents if _tmp != '']
        i_file_nevents = sorted(list(set(i_file_nevents)))

        if len(i_file_nevents) != 1:
           KILL('AAA')

        i_file_nevents = i_file_nevents[0]

        if not is_int(i_file_nevents):
           KILL('BBB')

        i_file_nevents = int(i_file_nevents)

        totEvents += i_file_nevents
        if (max_events > 0) and (totEvents >= max_events):
           breakLoop = True

        if verbose:
           print(i_file)
           print(i_file_nevents)

        i_file_parents2 = []

        i_file_parents1 = command_output_lines('dasgoclient --query "parent file='+str(i_file)+'"')
        i_file_parents1 = [_tmp for _tmp in i_file_parents1 if _tmp != '']
        i_file_parents1 = sorted(list(set(i_file_parents1)))
        for i_file_aodf in i_file_parents1:
            i_file_parents2_tmp = command_output_lines('dasgoclient --query "parent file='+str(i_file_aodf)+'"')
            i_file_parents2_tmp = [_tmp.replace(' ', '') for _tmp in i_file_parents2_tmp]
            i_file_parents2_tmp = [_tmp for _tmp in i_file_parents2_tmp if _tmp != '']
            i_file_parents2_tmp = sorted(list(set(i_file_parents2_tmp)))
            i_file_parents2 += i_file_parents2_tmp

        i_file_parents2 = sorted(list(set(i_file_parents2)))

        if verbose:
           for _tmp in i_file_parents2:
               print(' '*5, _tmp)

        dset_data['files'] += [{
          'file': i_file,
          'nevents': i_file_nevents,
          'parentFiles_1': i_file_parents1,
          'parentFiles_2': i_file_parents2,
        }]

        if breakLoop:
           break

    del totEvents, breakLoop

    # consistency checks on dataset data
    assert_dataset_data(dset_data=dset_data, verbose=verbose)

    return dset_data

def assert_dataset_data(dset_data, verbose=False):

    if not isinstance(dset_data, dict):
       KILL('assert_dataset_data -- invalid content of dataset .json file [-d]: '+str(file_path))

    if 'DAS' not in dset_data:
       KILL('assert_dataset_data -- 11')
    elif not isinstance(dset_data['DAS'], basestring):
       KILL('assert_dataset_data -- 12')

    if 'files' not in dset_data:
       KILL('assert_dataset_data -- 21')
    elif not isinstance(dset_data['files'], list):
       KILL('assert_dataset_data -- 22')

    for i_ent in dset_data['files']:

        if not isinstance(i_ent, dict):
           KILL('assert_dataset_data -- 31')

        if 'file' not in i_ent:
           KILL('assert_dataset_data -- 41 '+str(i_ent))
        elif not isinstance(i_ent['file'], basestring):
           KILL('assert_dataset_data -- 42 '+str(i_ent))

        if 'nevents' not in i_ent:
           KILL('assert_dataset_data -- 41 '+str(i_ent))
        elif not isinstance(i_ent['nevents'], int):
           KILL('assert_dataset_data -- 42 '+str(i_ent))
        elif i_ent['nevents'] < 0:
           KILL('assert_dataset_data -- 43 '+str(i_ent))

        for _tmp in ['parentFiles_1', 'parentFiles_2']:
            if _tmp not in i_ent:
               KILL('assert_dataset_data -- 51 '+str(i_ent))
            elif not isinstance(i_ent[_tmp], list):
               KILL('assert_dataset_data -- 52 '+str(i_ent))
            else:
               for _tmp2 in i_ent[_tmp]:
                   if not isinstance(_tmp2, basestring):
                      KILL('assert_dataset_data -- 53 '+str(i_ent))

def skim_das_jsondump(file_path, max_files=-1, max_events=-1, verbose=False):

    dset_data = json.load(open(file_path))

    # consistency checks on dataset data
    assert_dataset_data(dset_data, verbose=verbose)

    # skim input list based on max_files and max_events
    if max_files > 0:
       dset_data['files'] = dset_data['files'][:max_files]

    if max_events > 0:
       lastIndex, totEvents = 0, 0
       for _tmp in dset_data['files']:
           lastIndex += 1
           totEvents += _tmp['nevents']
           if totEvents >= max_events:
              break
       if lastIndex != len(dset_data['files']):
          dset_data['files'] = dset_data['files'][:lastIndex]
       del lastIndex, totEvents

    return dset_data
