from os import listdir
from os.path import isfile, join
import shutil


path = "."
fromPath = "/Users/greg/Dropbox/_mol/Segger/"

for fname in listdir(path) :
    print fname,
    if '__init__' in fname :
        print " - not copying"
        continue
    fpathTo = join(path, fname)
    fpathFrom = join(fromPath, fname)
    print " -< ", fromPath,
    if isfile( fpathTo ) and isfile ( fpathFrom ) :
        shutil.copyfile ( fpathFrom, fpathTo  )
        print " -ok- "
    else :
        print "?"
