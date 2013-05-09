----------------------------------------------------------------------------
1) Creating CRAB Jobs:
----------------------------------------------------------------------------

For automatic creation of CRAB jobs for several different datasets, the createCrabJobs.py script can be used.
For a brief description of the input parameters, execute

./createCrabJobs.py --help

which will produce the following output

Usage: createCrabJobs.py [options]
Example: ./createCrabJobs.py -w CRAB_Jobs -d datasetList.txt -c runBTagAnalyzer_cfg.py -t crab_template.cfg

Options:
  -h, --help            show this help message and exit
  -w MAIN_WORKDIR, --main_workdir=MAIN_WORKDIR
                        Main working directory
  -d DATASET_LIST, --dataset_list=DATASET_LIST
                        Text file containing a list of datasets to be
                        processed
  -c CMSSW_CFG, --cmssw_cfg=CMSSW_CFG
                        CMSSW configuration file
  -t CRAB_CFG_TEMPLATE, --crab_cfg_template=CRAB_CFG_TEMPLATE
                        CRAB configuration file template
  -n, --no_creation     Create the necessary configuration files and skip the
                        job creation (This parameter is optional)

The -w (--main_workdir) parameter defines the name of the main working directory, the -d (--dataset_list)
parameter specifies the location of a text file containing a list of datasets to be processed, the -c (--cmssw_cfg)
parameter specifies the location of a common CMSSW configuration file, and the -t (--crab_cfg_template) parameter
specifies the location of a CRAB configuration file template. With the input parameters provided, the createCrabJobs.py
script creates the following file and directory structure

MAIN_WORKDIR/ (Main working directory)
    |
    ------> datasetList.txt (Text file defining a list of datasets to be processed)
    |
    ------> cfg_files/ (Directory containing CMSSW and CRAB configuration files)
    |          |
    |          ------> CMSSW_cfg.py (Common CMSSW configuration file)
    |          |
    |          ------> DATASET1_crab.cfg (CRAB configuration file for DATASET1)
    |          |
    |          ------> DATASET1_lumi_mask.txt (Lumi mask file for DATASET1)
    |          |
    |          ------> DATASET2_crab.cfg (CRAB configuration file for DATASET2)
    |          |
    |          ------> DATASET2_lumi_mask.txt (Lumi mask file for DATASET2)
    |          |
    |          ...
    |
    ------> DATASET1/ (CRAB working directory for DATASET1)
    |
    ------> DATASET2/ (CRAB working directory for DATASET2)
    |
    ...

The dataset list file contains the following entries

# Dataset_name                                                                                                                                                                                  Unit(lumi:1,event:0)   Total_units   Number_of_jobs   Run_selection   Lumi_mask   Scheduler   Use_server             DBS_instance   Publication_name   PyCfg_parameters (optional)
/BprimeBprimeToBHBHinc_M-600_TuneZ2star_8TeV-madgraph/cjenkins-BprimeBprimeToBHBHinc_M600_TuneZ2star_8TeVmadgraphSum12_DR53X_PU_S10_START53_V7A_TLBSMv2-fe5dcf8cf2a24180bf030f68a7d97dda/USER                      0            -1                1               -           -      condor            0   cms_dbs_ph_analysis_02                  -   wantSummary=True
/BprimeBprimeToBHBHinc_M-800_TuneZ2star_8TeV-madgraph/cjenkins-BprimeBprimeToBHBHinc_M800_TuneZ2star_8TeVmadgraphSum12_DR53X_PU_S10_START53_V7A_TLBSMv2-fe5dcf8cf2a24180bf030f68a7d97dda/USER                      0            -1                1               -           -      condor            0   cms_dbs_ph_analysis_02                  -   wantSummary=True

These entries are used to define the working directories and CRAB configuration files for different datasets. Note that
the PyCfg_parameters entry is optional and can be left out. However, all other entries must be defined (if you don't want
to apply them, just set them to -). The lumi mask file can be specified as an absolute path or as a relative path with
respect to the dataset list file.


----------------------------------------------------------------------------
2) Other Operations with CRAB Jobs:
----------------------------------------------------------------------------

Once the CRAB jobs for different datasets have been created, you can manually submit them

crab -submit -c MAIN_WORKDIR/DATASET1/

check their status

crab -status -c MAIN_WORKDIR/DATASET1/

and finally get their output.

crab -getoutput -c MAIN_WORKDIR/DATASET1/

To obtain the final report (including the output lumi summary file), execute one final time

crab -status -c MAIN_WORKDIR/DATASET1/

followed by

crab -report -c MAIN_WORKDIR/DATASET1/

For publication of the output EDM files, execute

crab -publish -c MAIN_WORKDIR/DATASET1/

However, keep in mind that the publication step can be relatively long. It is therefore recommended to start the publication
step inside a 'screen' session (see, for example, http://news.softpedia.com/news/GNU-Screen-Tutorial-44274.shtml).

If you have many datasets in your dataset list file, you will probably find it more convenient to use the checkCrabJobs.py
script provided. This script enables you to perform the following operations on your CRAB jobs:

    status
    submit
    getoutput
    report
    publish
    kill

Note that the checkCrabJobs.py script can only perform the above operations on all jobs of a given dataset. To get output,
kill, or (re)submit a job range, a manual procedure described above has to be followed. For a brief description of
available input parameters, execute

./checkCrabJobs.py --help

which will produce the following output

Usage: checkCrabJobs.py [options]
Example: ./checkCrabJobs.py -w CRAB_Jobs

Options:
  -h, --help            show this help message and exit
  -w MAIN_WORKDIR, --main_workdir=MAIN_WORKDIR
                        Main working directory
  -s, --submit          Submit CRAB jobs (This parameter is optional)
  -g, --getoutput       Get output from CRAB jobs (This parameter is optional)
  -r, --report          Get CRAB report (This parameter is optional)
  -p, --publish         Publish the output of CRAB jobs (This parameter is
                        optional)
  -k, --kill            Kill CRAB jobs (This parameter is optional)

Note that only one optional parameter at a time is allowed. With no optional parameters specified, the CRAB -status command is called.
