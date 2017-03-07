
def create_graph ( smod, links ) :

    print "\nCreating graph for %s - %s" % (smod.name, links)

    if hasattr(smod, 'adj_graph') and smod.adj_graph:
        smod.adj_graph.close()
        smod.adj_graph = None

    from VolumePath import Marker_Set, Marker, Link

    gname = smod.name + " graph(%s)" % links
    g = Marker_Set ( gname )
    smod.adj_graph = g
    smod.graph_links = links
    aMap = dict()

    regions = smod.regions
    marker_radius = 0.1 * sum([r.enclosed_volume() ** (1./3) for r in regions]) / len(regions)
#    ijk_to_xyz_transform = smod.point_transform()
    from Matrix import apply_matrix
    for reg in regions :
#        xyz = apply_matrix(ijk_to_xyz_transform, reg.max_point)
        c = reg.center_of_points()
        m = Marker(g, reg.rid, c, reg.color, marker_radius)
        aMap[reg] = m
        m.region = reg
        m.extra_attributes = { 'region_id' : str(reg.rid),
                               'region_size': str(reg.point_count()) }

    link_color = ( .5, .5, .5, 1 )
    link_radius = 0.5 * marker_radius

    from regions import group_contacts
    cons = group_contacts(smod.region_contacts())

    Ns, min_N, max_N = [], None, None
    avgds, min_avgd, max_avgd = [], None, None
    maxds, min_maxd, max_maxd = [], None, None

    # first run through contacts to list contacts and properties
    if links == "avgd" or links == "maxd" or links == "N" :
        for r1 in cons.keys() :
            for r2 in cons[r1].keys() :
                if r2 > r1 :

                    con = cons[r1][r2]

                    #print "link %d -> %d -- N:%.1f, " % (r1.rid, r2.rid, con.N),

                    if con.N < 0.1 :
                        #print "*hmm*"
                        continue
                    else :
                        Ns.append ( con.N )

                    if con.D :
                        #print "D:%.3f, " % con.D,
                        avgd = float(con.D) / (2.0 * float(con.N))
                        avgds.append ( avgd )
                    #else : print "D:*",

                    if con.maximum_density :
                        #print "MaxD:%.3f, " % con.maximum_density
                        maxds.append ( con.maximum_density )
                    #else : print "MaxD:*"


        min_N, max_N = min(Ns), max(Ns)
        min_avgd, max_avgd = min(avgds), max(avgds)
        min_maxd, max_maxd = min(maxds), max(maxds)

        print "Avg densities: %.5f -> %.5f" % (min_avgd, max_avgd)
        print "N: %.1f -> %.1f" % (min_N, max_N)
        print "Maximum densities %.5f -> %.5f" % (min_maxd, max_maxd)


    min_rad = marker_radius * 0.1
    max_rad = marker_radius * 0.75 - min_rad

    for r1 in cons.keys() :
        for r2 in cons[r1].keys() :
            if r2 > r1 :

                con = cons[r1][r2]

                if links == "maxd" and con.maximum_density :
                    # radius proportional to max density at boundary
                    if con.N > 0.1 :
                        maxd = con.maximum_density
                        link_radius_var = min_rad + max_rad * (maxd - min_maxd)/(max_maxd-min_maxd)
                        Link ( aMap[r1], aMap[r2], link_color, link_radius_var )

                elif links == "N" and con.N:
                    # radius of link proportional to area of contact
                    #  - where area of contact ~ #voxels between regions
                    con = cons[r1][r2]
                    link_radius_var = min_rad + max_rad * (con.N - min_N)/(max_N-min_N)
                    Link ( aMap[r1], aMap[r2], link_color, link_radius_var )

                elif links == "avgd" and con.D:
                    # radius of link proportional to average density
                    #  - at boundary
                    if con.N > 0.1 :
                        avgd = float(con.D) / (2.0 * float(con.N))
                        link_radius_var = min_rad + max_rad * (avgd - min_avgd)/(max_avgd-min_avgd)
                        Link ( aMap[r1], aMap[r2], link_color, link_radius_var )
                else :
                    # same link radius for all links
                    Link ( aMap[r1], aMap[r2], link_color, link_radius )


    g.show_model ( True )
    smod.display = True
    smod.regions_scale = 0.5
    smod.display_regions ('Voxel_Surfaces', None, None, True)


def break_selected_links():

    from VolumePath import markerset
    for l in markerset.selected_links():
        l.delete()

def link_selected():

    from VolumePath import selected_markers, Link
    msel = selected_markers()
    if len(msel) == 2:
        m0, m1 = msel
        if m0.marker_set == m1.marker_set:
            link_color = ( .5, .5, .5, 1)
            link_radius = 0.5 * m0.radius
            Link(m0, m1, link_color, link_radius)

def open_skeleton ( smod ) :

    if smod.adj_graph and hasattr(smod.adj_graph, 'path'):
        import os.path
        initdir, initfile = os.path.split(smod.adj_graph.path)
    else:
        initdir, initfile = smod.path, smod.name + ' skeleton'
    def open(o, d, smod = smod):
        open_skeleton_file(o, d, smod)
    import OpenSave
    d = OpenSave.OpenModeless(title = 'Open skeleton',
                              initialdir = initdir,
                              initialfile = initfile,
                              filters = [('Chimera Markers', '*.cmm', '')],
                              multiple = False,
                              command = open)

def open_skeleton_file ( open, dialog, smod ) :

    if not open:
        return

    paths = dialog.getPaths()
    path = paths[0]
    import VolumePath
    g = VolumePath.open_marker_set(path)

    mlist = g.markers()
    for m in mlist:
        m.region = None
        rid = int(m.extra_attributes['region_id'])
        if rid in smod.id_to_region:
            r = smod.id_to_region[rid]
            rsize = int(m.extra_attributes['region_size'])
            if r.point_count() == rsize:
                m.region = r

#TODO: The regions file format renumbers the regions consecutively.  So
# matching based on region id won't work.  Should change file format to
# hdf5 and include region ids, colors, and any other useful per-region info.

    nomatch = [m for m in mlist if m.region is None]
    if nomatch:
        umsg('%d of %d skeleton nodes did not match a region' %
             (len(nomatch), len(mlist)))

    if smod.adj_graph:
        smod.adj_graph.close()
    smod.adj_graph = g


def save_skeleton ( smod ) :

    if smod.adj_graph is None:
        umsg('No skeleton for %s' % smod.name)
        return

    if hasattr(smod.adj_graph, 'path'):
        import os.path
        initdir, initfile = os.path.split(smod.adj_graph.path)
    else:
        initdir, initfile = smod.path, smod.name + ' skeleton'
    import OpenSave
    d = OpenSave.SaveModal(title = 'Save skeleton',
                           initialdir = initdir,
                           initialfile = initfile,
                           filters = [('Chimera Markers', '*.cmm', '')])
    from chimera.tkgui import app
    paths_and_types = d.run(app)
    if paths_and_types:
        path = paths_and_types[0][0]
        out = open(path, 'w')
        from VolumePath import markerset
        markerset.save_marker_sets([smod.adj_graph], out)
        out.close()


def close ( smod ) :

    g = smod.adj_graph
    if g:
        g.close()
        smod.adj_graph = None
        smod.regions_scale = 1.0
        smod.display_regions ('Voxel_Surfaces', None, None, True)


def group_by_skeleton ( smod ) :

    g = smod.adj_graph
    if g is None:
        return

    # Find connected groups of markers.
    msets = connected_markers(g)

    # Find connected regions.
    rgroups = [[m.region for m in mset] for mset in msets]

    # Exclude region groups that are already grouped correctly.
    rgroups = [rgroup for rgroup in rgroups if not is_region_group(rgroup)]

    # Find current most common color for each group.
    colors = [most_common_region_color(rgroup) for rgroup in rgroups]

    # Make split regions have different colors.
    csize = {}
    for rg, c in zip(rgroups, colors):
        if c in csize:
            csize[c] = max(csize[c], len(rg))
        else:
            csize[c] = len(rg)
    from regions import random_color
    for i, c in enumerate(colors):
        if len(rgroups[i]) < csize[c]:
            colors[i] = random_color()

    # Ungroup regions that need regrouping
    remove_parents(concatenate(rgroups), smod)

    # Group connected associated regions
    for rgroup, c in zip(rgroups, colors):
        r = smod.join_regions ( rgroup )
        r.color = c


def remove_parents(regions, smod):

    while True:
        from regions import TopParentRegions
        parents = TopParentRegions ( [r for r in regions if r.preg] )
        if parents:
            smod.ungroup_regions ( parents )
        else:
            break

def is_region_group(regions):

    r0 = regions[0]
    if r0.preg is None:
        return False
    p = r0.preg
    if p.preg:
        return False
    return same_list_elements(p.cregs, regions)

def most_common_region_color(regions):

    ct = {}
    for r in regions:
        c = tuple(r.top_parent().color)
        if c in ct:
            ct[c] += 1
        else:
            ct[c] = 1
    count, color = max([(count, c) for c, count in ct.items()])
    return color

def connected_markers(mset):

    cm = {}
    for l in mset.links():
        m1, m2 = l.marker1, l.marker2
        if m1.region is None or m2.region is None:
            continue
        cm1 = cm.setdefault(m1, set([m1]))
        cm2 = cm.setdefault(m2, set([m2]))
        if not cm2 is cm1:
            cm1.update(cm2)
            for m in cm2:
                cm[m] = cm1
    msets = dict([(id(ms), ms) for ms in cm.values()]).values()

    # Add lone markers.
    for m in mset.markers():
        if m.region and len(m.links()) == 0:
            msets.append(set([m]))

    return msets

def concatenate(lists):

    c = []
    for e in lists:
        c.extend(e)
    return c

def same_list_elements(list1, list2):

    if len(list1) != len(list2):
        return False
    set2 = set(list2)
    for e in list1:
        if not e in set2:
            return False
    return True
