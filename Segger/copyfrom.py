from os import listdir
from os.path import isfile, join
import shutil


path = "."
fromPath = "/Users/greg/Dropbox/_mol/Segger/"

for fname in listdir(path) :
    print fname,
    fpathTo = join(path, fname)
    fpathFrom = join(fromPath, fname)
    print " -< ", fromPath,
    if isfile( fpathTo ) and isfile ( fpathFrom ) :
        print " - "
        shutil.copyfile ( fpathFrom, fpathTo  )
    else :
        print "?"
