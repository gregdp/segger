

dmapPath = "/Users/greg/Box Sync/_data/Ribozyme2b/Con2-2.2/Con2-2.2A_sh.mrc"
half1Path = "/Users/greg/Box Sync/_data/Ribozyme2b/Con2-2.2/Con2-2.2A_sh_half_A.ccp4" # can be None to not use half map A
half2Path = "/Users/greg/Box Sync/_data/Ribozyme2b/Con2-2.2/Con2-2.2A_sh_half_B.ccp4" # can be None to not use half map B

if 0 :
    half1Path = None
    half2Path = None

molPath = "/Users/greg/Box Sync/_data/Ribozyme2b/Con2-2.2/Con2hf-coot-2.pdb"
outMolPath = "/Users/greg/Box Sync/_data/Ribozyme2b/Con2-2.2/Con2hf-coot-2_with_ions_and_water.pdb"

thrSigma=3.0 # contour level above which water/ions are placed, so contour level is map_avg + thrSigma * map_stdev
useQ=True # use Q-score in placing water/ions
minQ=0.7 # only place water/ion if Q is above this value
sigQ=0.6 # sigma to use in calculating Q, use 0.6 for 1.5A or lower maps, 0.4 for 1.0 to 1.5A
toChain='' # leave empty to auto-select chain based on nearest atom

minDistI, maxDistI = 1.8, 2.5 # ion distances min - max
minDistW, maxDistW = 2.5, 3.4 # water distances min - max

print ""
print "SWIM"
print ""

from VolumeViewer import open_volume_file

def openMap ( mapPath ) :
    openedMap = None
    try :
        openedMap = open_volume_file ( mapPath, 'ccp4')[0]
        #print " - opened map: %s" % mapPath
    except :
        print " - could not open map: %s" % mapPath
        exit(1)
    return openedMap

dmap = openMap ( dmapPath )

hMapA, hMapB = None, None
if half1Path :
    hMapA = openMap ( half1Path )
if half2Path :
    hMapB = openMap ( half2Path )

from chimera import openModels
mol = openModels.open ( molPath, type='PDB' )[0]
print " - opened model: %s" % molPath
print ""

M = dmap.data.full_matrix()
from numpy import std, average
sdev = std(M)
avg = average(M)
mapThreshold = avg + thrSigma * sdev

print " - for sigma of %.2f, threshold is %.4f" % (thrSigma, mapThreshold)

from Segger import regions
smod = regions.Segmentation(dmap.name, dmap)
smod.calculate_watershed_regions ( dmap, mapThreshold )

print " - got %d regions" % len ( smod.regions )

nearAtoms = [at for at in mol.atoms if not at.element.name == "H"]

from Segger import SWIM
addW, addI = SWIM.goSWIM ( dmap, smod, mol, nearAtoms, toChain, minDistI, maxDistI, minDistW, maxDistW, hMapA, hMapB, useQ, minQ, sigQ )

print "Added %d waters, %d ions" % ( len(addW), len(addI) )

print " - saving to: %s" % outMolPath
import chimera
chimera.PDBio().writePDBfile ( [mol], outMolPath )





# that's it
