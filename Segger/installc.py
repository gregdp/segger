import sys, os, shutil


if len(sys.argv) < 2 :
    exit(0)
    
fromp = sys.argv[1]
top = os.path.split(__file__)[0]

for l in os.listdir ( top ) :

    print l,
    
    fromf = fromp + "/" + l
    tof = top + "/" + l
    try :
        shutil.copyfile(fromf, tof)
        print " - copied"
    except :
        print " - x"
