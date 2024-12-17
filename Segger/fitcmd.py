

import Segger
import Segger.fit

#import fit
print "hi"

#execfile ( "/Users/greg/Dropbox/_mol/Segger/scripts/fit.py" )
Segger.fit.fitAll ( numProc=4, metric="ccm", resolution=2.0, nrot=8, step=10, numOut=30 )

#fit.fitAll ( numProc=6 )
#fit.fitAll ( numProc=1 )
