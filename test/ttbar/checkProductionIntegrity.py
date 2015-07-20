import os
import sys
import optparse
from storeTools import getEOSlslist

"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inDir',       dest='inDir',       help='input directory with files',      default=None,   type='string')
    parser.add_option('-c', '--cleanup',     dest='cleanup',     help='removes original crab directory', default =False, action='store_true')
    parser.add_option(      '--nocheck',     dest='nocheck',     help='do not prompt user',              default=False,  action='store_true')
    (opt, args) = parser.parse_args()

    dset_list=getEOSlslist(directory=opt.inDir,prepend='')
    for dset in dset_list:
        dsetname=dset.split('/')[-1]

        pub_list=getEOSlslist(directory=dset,prepend='')
        if len(pub_list)!=1 or  not 'crab' in pub_list[0]:
            print 'Ambiguity found @ <publication-name> for <primary-dataset>=%s , bailing out'%dsetname
            continue
        pub=pub_list[0].split('/crab_')[-1]

        time_list=getEOSlslist(directory=pub_list[0],prepend='')
        if len(time_list)!=1:
            print 'Ambiguity found @ <time-stamp> for <primary-dataset>=%s , bailing out'%dsetname
            continue
        time_stamp=time_list[0].split('/')[-1]

        out_list=[]
        count_list=getEOSlslist(directory=time_list[0],prepend='')
        for count in count_list: out_list += getEOSlslist(directory=count,prepend='')
        file_list=[x for x in out_list if '.root' in x]

        newDir='%s/%s' % (opt.inDir,pub)        
        print '<primary-dataset>=%s <publication-name>=crab_%s <time-stamp>=%s has %d files' % (dsetname,pub,time_stamp,len(file_list) )
        if not opt.nocheck:
            choice = raw_input('Will move to %s current output directory. [y/n] ?' % newDir ).lower()
            if not 'y' in choice : continue

        os.system('cmsMkdir %s' % newDir)
        for f in file_list : os.system('cmsStage -f %s %s/' % (f, newDir) )
        if not opt.nocheck and opt.cleanup : 
            choce = raw_input('Will remove output directory. [y/n] ?').lower()
            if 'y' in choice: os.system('cmsRm -r %s' % dset)

        print 'Crab outputs may now be found in %s' % newDir

    print '-'*50
    print 'All done. In case errors were found check that the crab output structure is '
    print '<outLFNDirBase>/<primary-dataset>/<publication-name>/<time-stamp>/<counter>/<file-name>'
    print '-'*50
        


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())

