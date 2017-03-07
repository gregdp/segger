import os.path
import numpy
#from regions import Segmentation, Region, Contact
from sys import stderr



def ReadRegionsFile ( regions_file_path, dmap, smod = None) :

    print "Reading regions ---"

    try :
        e = numpy.fromfile ( regions_file_path, numpy.double )
    except :
        print "could not read:", regions_file_path
        return None

    print " - read ", len(e)

    import regions
    reload (regions)

    if smod == None :
        regions_file = os.path.basename ( regions_file_path )
        smod = regions.Segmentation(regions_file, dmap)
    else :
        print " - found", smod.name
        smod.remove_all_regions()

    smod.path = os.path.dirname ( regions_file_path ) + os.path.sep
    smod.adj_graph = None
    print " - " + smod.path + smod.name

    regions, groups, at = ParseRegions ( e, smod )

    smod.regions = set(groups)
    smod.id_to_region = regions

    smod.rcons = ParseContacts(e, at, regions)

    return smod



def ParseRegions ( e, smod ) :
    
    nregions = int ( e[0] )
    print " - reading %d regions..." % nregions
    at = 1

    regs = {}
    all_regions = {}            # Includes groups.

    import regions
    reload (regions)

    for i in range ( nregions ) :

        try : nvoxels = int ( e[at] )
        except :
            print " - reached end of file before reading all regions"
            break

        at += 1
        rvs = e [ at : at + (nvoxels*3) ]
        at += nvoxels*3
        
        print "Region %d - %d voxels" % (i, nvoxels)

        rpoints = numpy.reshape ( rvs, (nvoxels, 3) ).astype ( numpy.int32 )
        
        #print rpoints

        nparents = int ( e[at] )
        at += 1
        parents = e [ at : at + nparents ].astype ( numpy.int )
        at += nparents

        rid = i+1
        reg = regions.Region ( smod, rid, rpoints[0] )

        smod.mask[rpoints[:,2],rpoints[:,1],rpoints[:,0]] = rid  # set mask at points

        all_regions [ reg.rid ] = reg
        regs [ reg.rid ] = reg

        last_reg = reg
        reg.preg = None

        for pi in parents :

            if pi in all_regions:
                preg = all_regions[pi]
            else :
                preg = regions.Region ( smod, pi )
                preg.max_point = rpoints[0]
                all_regions[pi] = preg

            last_reg.preg = preg

            if preg.cregs.count ( last_reg ) == 0 :
                preg.cregs.append ( last_reg )

            last_reg = preg



    # Regions table only includes top level groups.
    groups = [ reg for reg in all_regions.values() if reg.preg is None ]

    return all_regions, groups, at



def ParseContacts ( e, at, regs ):
    
    import regions
    reload (regions)

    try : ncon = int ( e[at] )
    except :
        print " - reached end of file before reading contacts"
        ncon = 0

    at += 1
    rcons = {}

    if ncon > 0 :

        print " - reading %d contacts..." % ( ncon )
    
        am = e [ at : at + (ncon*4) ]
        cm = numpy.reshape ( am, (ncon, 4) )

        for i in range ( ncon ) :
            c = cm[i]
            rid1, rid2 = int(c[0]), int(c[1])
            #from regions import Contact
            o = regions.Contact(c[2])
            o.D = c[3]

            try : r1 = regs[rid1]
            except : print "File error: contact region id", rid1; continue
            try : r2 = regs[rid2]
            except : print "File error: contact region id", rid2; continue

            if r1 == r2 : print "File error: self contact id", rid1

            if not r1 in rcons : rcons[r1] = {}
            if not r2 in rcons : rcons[r2] = {}

            rcons[r1][r2] = o
            rcons[r2][r1] = o

    print ""

    return rcons



def WriteRegionsFile ( smod, fname = None ) :

    if fname is None:
        # Show save-file dialog.
        def save ( okay, dialog ):
            if okay:
                paths = dialog.getPaths ( )
                if paths:
                    WriteRegionsFile ( smod, paths[0] )

        bname = smod.name [ 0 : smod.name.rfind ('_regions') ]
        prefix = smod.path + bname + "_regions_save_%d"
        uupath = unusedFile ( prefix )
        from OpenSave import SaveModeless
        d = SaveModeless ( title = 'Save Regions',
                           initialdir = smod.path,
                           initialfile = os.path.basename(uupath),
                           command = save )
        return


    tot_write_regions = [0]
    tot_e_size = 1
    for region in smod.regions :
        tot_e_size += RegionSize ( region, tot_write_regions )


    num_cons = 0
    rcons = smod.region_contacts()
    for r1, cr1 in rcons.iteritems () :
        for r2, o in cr1.iteritems () :
            if r1.rid < r2.rid :
                num_cons = num_cons + 1

    tot_e_size = tot_e_size + 1 + 4 * (num_cons)


    print "Writing %d regions, %d grouped" % ( tot_write_regions[0], len(smod.regions) )
    if fname : print " - to", fname

    e = numpy.zeros ( [tot_e_size], numpy.float32 )

    e[0] = float ( tot_write_regions[0] )
    e_at = 1

    rlist = renumberRegions(smod)
    for region in rlist :
        e_at = AddRegion ( smod, region, e, e_at )

    print " - writing %d contacts" % ( num_cons )

    e[e_at] = float ( num_cons )
    e_at = e_at + 1

    #consa = []
    for r1, cr1 in rcons.iteritems () :
        for r2, o in cr1.iteritems () :
            if r1.rid < r2.rid :
                #consa = consa + [r1.rid, r2.rid, o.N, o.D]
                e[e_at+0] = float ( r1.rid )
                e[e_at+1] = float ( r2.rid )
                e[e_at+2] = o.N
                e[e_at+3] = o.D
                e_at = e_at + 4

    #consa = [len(consa)/4] + consa
    #e = numpy.concatenate ( [ e, numpy.array ( consa, numpy.float32 ) ] )

    e.tofile ( fname )
    print "Wrote %s" % os.path.basename(fname)


def AddRegion ( smod, region, e, e_at ) :

    if len(region.cregs) == 0 :

        e[e_at] = float ( region.point_count() )
        e_at = e_at + 1

        for rp in region.points() :
            e[e_at] = rp[0]; e_at = e_at + 1
            e[e_at] = rp[1]; e_at = e_at + 1
            e[e_at] = rp[2]; e_at = e_at + 1

        rp = region
        parents = []
        while rp.preg != None :
            rp = rp.preg
            parents.append ( rp.rid )

        e[e_at] = float ( len(parents) )
        e_at = e_at + 1

        for pi in parents :
            e[e_at] = float ( pi )
            e_at = e_at + 1

    else :
        for creg in region.cregs :
            e_at = AddRegion ( smod, creg, e, e_at )

    return e_at



def RegionSize ( region, tot_write_regions ) :

    if len(region.cregs) == 0 :            

        rp = region
        nparents = 0
        while rp.preg != None :
            rp = rp.preg
            nparents = nparents + 1

        tot_write_regions[0] += 1

        return 1 + region.point_count()*3 + 1 + nparents

    else :
        size = 0
        for creg in region.cregs :
            size = size + RegionSize ( creg, tot_write_regions )
        return size


#
# Number regions having no children in depth first order since that is the
# order they will be written to the file.  Renumber nodes with children
# using higher numbers to keep all region numbers distinct.
#
def renumberRegions(smod):

    newrid = {}
    parents = []
    rlist = list(smod.regions)
    rlist.sort(lambda r1, r2: cmp(r1.rid, r2.rid))
    for r in rlist:
        renumberRegion(r, newrid, parents)
    for r in parents:
        next_id = len(newrid) + 1
        newrid[r] = next_id

    for r, rid in newrid.items():
        r.rid = rid

    smod.id_to_region = dict([(r.id,r) for r in smod.id_to_region.values()])

    return rlist

def renumberRegion(r, rid, parents):

    if len(r.cregs) == 0:
        next_id = len(rid) + 1
        rid[r] = next_id
    else:
        parents.append(r)
        for c in r.cregs:
            renumberRegion(c, rid, parents)

def unusedFile ( path_format ):

    i = 1
    exists = True
    while exists:
        path = path_format % (i,)
        exists = os.path.exists ( path )
        i += 1
    return path
