


from os import listdir
from os.path import isfile, join, splitext
import shutil
import sys

print sys.argv[1]
fromPath = sys.argv[1]

path = "."
print ""

copied, notCopied = [], []
print "Copying..."
print "-----------"
for fname in listdir(fromPath) :
    if '__init__' in fname :
        notCopied.append (fname)
        continue
    if ".pyc" == splitext(fname)[1] :
        continue
    fpathTo = join(path, fname)
    fpathFrom = join(fromPath, fname)
    #print " -< ", fromPath,
    if isfile( fpathTo ) and isfile ( fpathFrom ) :
        shutil.copyfile ( fpathFrom, fpathTo  )
        #print fname
        copied.append ( fname )
    else :
        notCopied.append (fname)

print ""


print "Copied:"
print "-----------"
copied.sort()
for f in copied :
    print f
print ""

print "Not copied:"
print "-----------"
notCopied.sort()
for f in notCopied :
    if splitext(f)[1] != ".pyc" :
        print f
print ""

print "Param:"
print "-----------"
for f in listdir ( fromPath + "/_param" ) :

    fname, fext = splitext (f)
    if fext == ".pdb" :
        print f,
        shutil.copy2 ( fromPath + "/_param/" + f, "./_param/" + f )

print ""
print ""
