#!/usr/bin/env python
import json

datasets = json.load(open('datasets.json'))

for i_dset in sorted(datasets.keys()): print i_dset
