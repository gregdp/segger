


# How to run:
# This script should be copied and the parameters changed as needed
#  - or this script can be used as placed in the Chimera/Contents/Resources/share/ folder
#    with default parameters
# - command to use:
#   [path to Chimera exec] --nogui --silent --nostatus [path to map] [path to this script]
#   e.g. ~/_mol/Chimera.app/Contents/MacOS/chimera --nogui --silent --nostatus
#        ~/_data/emd_5001.map seggercmd.py


# The grouping mode paramater:
#  - can be either 'smoothing' or 'connectivity'
#      'smoothing' tends to work better at lower resolutions (4A and lower)
#      'connectivity' tends to work better at higher resolutions (4A and better)
groupingMode = "smoothing"

# Minimum region size in number of voxels (can be turned to A^3 by dividing by
# step size)
#  - as below when set to 1, all regions will be kept; for a value of 5, only
#    regions that have at least 5 voxels would be kept after the first
#    segmentation step
#    (1 means no regions are removed)
minRegionSize = 1

# Minimum contact voxels -
#  (0 means no regions are removed)
minContactVoxels = 0

# when to stop the grouping process (either by)
stopAtNumberOfRegions = 1

# map threshold - only include voxels with map value above this values
# if not sure, leave as None; 3sigma above mean will be used
mapThreshold = None

# parameters for smoothing and grouping:
numSmoothingSteps = 4
smoothingStepSize = 3

# parameters for grouping by connectivity
numConnectivitySteps = 10


# parameters for outputting of regions
outputLargestN = 5
outputRegions = [] # a list of regions to output

# other options
options = {}
options["outputMapSize"] = "box" # can be 'same', 'cube', 'box'
options["borderWidth"] = 4 # border in voxels if mapSize is not "same"
options["mask01"] = False # border in voxels if mapSize is not "same"



import Segger
import chimera
import VolumeViewer
import regions
import numpy

maps = chimera.openModels.list (modelTypes = [VolumeViewer.volume.Volume])

if len(maps) == 0 :
    print " - no maps opened; specify a map on the command line"

dmap = maps[0]

if mapThreshold == None :
    #maxM = numpy.max(M)
    #minM = numpy.min(M)

    M = dmap.data.full_matrix()
    mapThreshold = numpy.average(M)+numpy.std(M)*3.0
    print " - threshold: %f" % mapThreshold


smod = regions.Segmentation(dmap.name, dmap)
smod.calculate_watershed_regions ( dmap, mapThreshold )

if minRegionSize > 1 :
    print "\n - removing regions below %d voxels" % minRegionSize
    smod.remove_small_regions(minRegionSize)

if minContactVoxels > 0 :
    print "\n - removing regions with connection smaller than %d voxels" % minContactVoxels
    smod.remove_contact_regions(minContactVoxels)


if groupingMode == 'smoothing' :
    print "\n - grouping by smoothing"
    smod.smooth_and_group(numSmoothingSteps, smoothingStepSize, stopAtNumberOfRegions)

elif groupingMode == 'connectivity' :
    print "\n - grouping by connectivity"
    self.GroupByCons ( smod, task )


regions = smod.grouped_regions()
print " - %d regions" % len(regions)
outputN = min ( outputLargestN, len(regions) )


import os
mdir, mfile = os.path.split(dmap.data.path)
mname, mext = os.path.splitext ( mfile )


if outputN > 0 :
    print "\n - outputing %d regions by size" % outputN

    import Segger.extract_region_dialog

    for i in range ( outputN ) :
        reg = regions[i]
        print " -- %d, %d voxels" % (i+1, len(reg.points()) )

        newMapName = mname + "_%d.mrc" % (i+1)

        Segger.extract_region_dialog.ExtractNonUI (dmap, dmap, smod, [reg], newMapName, options)
