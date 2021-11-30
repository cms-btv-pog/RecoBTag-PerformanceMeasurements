import json

from datasets import *

dictt = dataset_dict()

data = {}

for key in dictt:

    data[key] = {
      "Data.inputDataset"  : dictt[key][1],
      "Data.splitting"     : dictt[key][0][0],
      "Data.unitsPerJob"   : dictt[key][0][1],
      "JobType.pyCfgParams": dictt[key][2].split(','),
    }

with open('tmpcopy.json', 'w') as fp: json.dump(data, fp, indent=2, sort_keys=True)
