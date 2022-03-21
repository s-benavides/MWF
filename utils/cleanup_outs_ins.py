import os,glob,pathlib
import numpy as np

dirs='./'

skip = False

tf_paths = list(pathlib.Path().glob(dirs+'/Lx*[!.zstd]/'))
print(tf_paths)

for path in tf_paths:
#	if (('O50B40T90' in str(path)) or ('O10B40T90' in str(path))):
#		print("Skipping  %s" % os.path.split(path)[0][:-4])
    if skip:
        print("Skipping")
        pass
    else:
        print("Working on %s, state" % path)
        # Find last output
        snapshots = list(pathlib.Path().glob(dirs+'/'+str(path)+'/state*.dat'))
#        print(snapshots)
        last_out = np.max(np.array([int((str(snap).split('.')[0]).split('/')[-1][5:]) for snap in snapshots]))
        last_out ="{:0>4s}".format(str(last_out))
        print("Last out: %s" % last_out)

        dels = list(set(pathlib.Path().glob(dirs+'/'+str(path)+'/state*.dat'))-
                    set(pathlib.Path().glob(dirs+'/'+str(path)+'/state'+str(last_out)+'.cdf.dat')))
#        print(dels)
        print(" --- Deleting --- ")
        # Fields:
        count = 0
        for ffile in dels:
            os.remove(ffile)
            count+=1
#            print(ffile)
        if count>0:
            print("Done! Removed %s out files" % count)

        print("Working on %s, turb" % path)
        # Find last output
        snapshots = list(pathlib.Path().glob(dirs+'/'+str(path)+'/turb*.dat'))
#        print(snapshots)
        last_out = np.max(np.array([int((str(snap).split('.')[0]).split('/')[-1][4:]) for snap in snapshots]))
        last_out ="{:0>4s}".format(str(last_out))
        print("Last out: %s" % last_out)

        dels = list(set(pathlib.Path().glob(dirs+'/'+str(path)+'/turb*.dat'))-
                    set(pathlib.Path().glob(dirs+'/'+str(path)+'/turb'+str(last_out)+'.cdf.dat')))
#        print(dels)
        print(" --- Deleting --- ")
        # Fields:
        count = 0
        for ffile in dels:
            os.remove(ffile)
            count+=1
#            print(ffile)
        if count>0:
            print("Done! Removed %s out files" % count)
