import numpy
from _surface import SurfaceModel

from sys import stderr
from Segger import timing
from time import clock

class Segmentation ( SurfaceModel ):

    def __init__(self, name, volume = None, open = True):

        SurfaceModel.__init__(self)

        self.name = name
        if volume is None:
            self.mask = None
            print " - no mask?"
        else:
            from numpy import zeros, uint32
            self.mask = zeros(volume.data.size[::-1], uint32)
            print " - made mask from", volume.name
        self.regions = set()            # Top level regions
        self.id_to_region = {}
        self.max_region_id = 0
        self.smoothing_level = 0
        self.rcons = None               # Leaf region contacts.
        self.seg_map = volume           # Map being segmented.
        self.map_level = None           # Good contouring level.
        tf = None if volume is None else volume.data.ijk_to_xyz_transform
        print " - ijk_to_xyz tr: "
        print tf
        self.ijk_to_xyz_transform = tf  # mask to physical coords
        self.surface_resolution = 1     # voxels

        self.adj_graph = None           # graph with regions as nodes
        self.graph_links = "uniform"    # how to compute radii of links in graph
        self.regions_scale = 1.0        # for shrinking regions in the graph

        from chimera import openModels as om
        if open:
            om.add([self])
            if not volume is None:
                self.openState.xform = volume.openState.xform

    def volume_data(self):

        v = self.seg_map
        if v and v.__destroyed__:
            v = self.seg_map = None
        return v

    def set_volume_data(self, volume):

        self.seg_map = volume
        if not volume is None:
            self.ijk_to_xyz_transform = volume.data.ijk_to_xyz_transform

    def point_transform(self):

        v = self.volume_data()
        if v:
            return v.data.ijk_to_xyz_transform
        return self.ijk_to_xyz_transform

    def grid_size(self):

        return tuple(self.mask.shape[::-1])

    def grid_origin(self):

        tf = self.point_transform()
        return tuple([tf[a][3] for a in (0,1,2)])

    def grid_step(self):

        tf = self.point_transform()
        from Matrix import length
        step = [length([tf[r][c] for r in (0,1,2)]) for c in (0,1,2)]
        return step

    def voxel_volume(self):

        t = self.point_transform()
        if t is None:
            return 1

        from numpy.linalg import det
        v = det([r[:3] for r in t])
        return v

    def calculate_region_bounds(self):

        import _segment
        b = _segment.region_bounds(self.mask)
        for r in self.childless_regions():
            i = r.rid
            npts = b[i,6]
            if npts > 0:
                r.rbounds = (b[i,:3],b[i,3:6])
                r.npoints = npts
            else:
                self.mask_id = None             # Not found in mask
                r.rbounds = ((1,1,1),(0,0,0))
                r.npoints = 0

    def region_contacts(self, task = None):

        if self.rcons is None:
            import _segment
            if timing: t0 = clock()
            rcon = _segment.region_contacts(self.mask)
            if timing: t1 = clock()
            rcons = {}
            id2r = self.id_to_region
            for i, (r1id, r2id, ncon) in enumerate(rcon):
                r1 = id2r[r1id]
                r2 = id2r[r2id]
                c = Contact(ncon)
                if not r1 in rcons:
                    rcons[r1] = {}
                if not r2 in rcons:
                    rcons[r2] = {}
                rcons[r1][r2] = rcons[r2][r1] = c
                if task and i % 1000 == 0:
                    task.updateStatus('Calculating region contacts %d of %d' %
                                      (i, len(rcon)))
            if timing: t2 = clock()
            v = self.volume_data()
            if v:
                import _segment
                map = v.full_matrix()
                ci, cf = _segment.interface_values(self.mask, map)
                id2r = self.id_to_region
                for (r1id,r2id,nv), (dmax,dsum) in zip(ci,cf):
                    c = rcons[id2r[r1id]][id2r[r2id]]
                    c.maximum_density = dmax
                    c.D = dsum

            if timing: t3 = clock()

            self.rcons = rcons
            cc = sum([len(t) for t in rcons.values()])/2
            print 'Computed %d contacting pairs' % cc
            if timing:
                print 'Time %.2f: contact calc %.2f sec, contact objects %.2f sec, maxima %.2f sec' % (t3-t0,t1-t0,t2-t1,t3-t2)


        return self.rcons

    def contacts_changed(self):

        self.rcons = None

    def open_map(self, open = True):

        print self.name, ":",

        v = self.volume_data()
        if v :
            print "already has map set"
            return None

        if not hasattr(self, 'map_path') and not hasattr(self,'map_name'):
            print "doesn't have map_path or name"
            return None

        #p = self.map_path
        #if p == None :
        #    print "map_path is None"
        #    return None

        dmap = None

        from VolumeViewer import volume_list, open_volume_file
        for v in volume_list() :
            if hasattr(self, 'map_path') and v.data.path == self.map_path :
                dmap = v
                print "found open map by path"
                break
            elif hasattr(self, 'map_name') and self.map_name in v.name :
                dmap = v
                print "found open map by name"
                break

        if dmap == None :
            from os.path import isfile
            if not isfile(self.map_path) :
                print "Map file is not on disk"
                return None

            if open:
                kw = {'model_id': self.id}
            else:
                kw = {'open_models': open}

            maps = open_volume_file(self.map_path, **kw)
            if len(maps) == 0 :
                print "error opening map file"
                return None

            print "loaded map"
            dmap = maps[0]

        if dmap == None :
            print "could not find/open map"
            return None

        self.set_volume_data(dmap)
        if not self.map_level is None:
            dmap.set_parameters(surface_levels = [self.map_level])
            #dmap.show()

        return dmap


    def remove_all_regions(self):

        for r in self.regions:
            r.remove_surface()
        self.regions.clear()
        self.id_to_region.clear()

    def remove_region (self, region,
                       remove_children = False,
                       remove_childless_parents = True,
                       update_surfaces = True,
                       task = None):

        self.remove_regions([region], remove_children,
                            remove_childless_parents, update_surfaces, task)

    def remove_regions (self, regions,
                       remove_children = False,
                       remove_childless_parents = True,
                       update_surfaces = False,
                       task = None):

        rset = set(regions)

        if remove_children:
            for r in regions:
                rset.update(r.all_children())

        if remove_childless_parents:
            check = rset
            while check:
                childless = set()
                for r in check:
                    p = r.preg
                    if p and not p in rset:
                        if len([c for c in p.cregs if not c in rset]) == 0:
                            childless.add(p)
                rset.update(childless)
                check = childless

        parents = set()
        for i, r in enumerate(rset):

            if task and i % 100 == 0:
                task.updateStatus('remove region %d of %d' % (i, len(rset)))

            if r in self.regions:
                self.regions.remove(r)
            del self.id_to_region[r.rid]
            #print " - removed r %d" % r.rid

            if update_surfaces :
                r.remove_surface()

            if len(r.cregs) == 0:
                p = r.points()
                if not p is None:
                    self.mask[p[:,2],p[:,1],p[:,0]] = 0   # zero mask at points
                    self.contacts_changed()

            for c in r.cregs:
                if not c in rset:
                    c.preg = None
                    self.regions.add(c)
                    if update_surfaces:
                        c.make_surface()

            p = r.preg
            if p and not p in rset:
                parents.add(p)

            i += 1

        for p in parents:
            p.remove_children([c for c in p.cregs if c in rset],
                              update_surfaces = update_surfaces)


    def remove_small_regions ( self, minRegionSize, task = None ) :

        #rcons = self.region_contacts(task)
        #regions = SmallConnectedRegions ( self.childless_regions(), rcons, minRegionSize, task )
        #self.remove_regions(regions, task = task)

        rregs = []
        for r in self.regions :
            np = len ( r.points() )
            if np < minRegionSize :
                rregs.append ( r )

        self.remove_regions(rregs, remove_children=True, update_surfaces=True, task = task)
        print 'Removed %d regions smaller than %d voxels' % (len(rregs), minRegionSize)


    def remove_contact_regions ( self, minContactSize, task = None ) :

        rcons = self.region_contacts(task)

        rregs = []
        for r in self.all_regions() :
            csize = 0
            if r in rcons:
                for rc, con in rcons[r].iteritems() :
                    csize += con.N
            if csize < minContactSize :
                rregs.append ( r )

        self.remove_regions(rregs, remove_children=True, update_surfaces=True, task = task)

        print 'Removed %d regions with <%d contact size' % (len(rregs), minContactSize)


    def join_regions ( self, regs, max_point = None, color = None ) :

        if len(regs) == 0 :
            return None

        rid = self.max_region_id + 1
        self.max_region_id += 1

        if max_point is None:
            from numpy import sum, float32, int32
            try :
                ave = sum([r.max_point for r in regs], axis=0).astype(float32) / len(regs)
                max_point = ave.astype(int32)
            except :
                max_point = None

        nreg = Region ( self, rid, children = regs, max_point = max_point )

        if color is None:
            lr = max([(r.point_count(),r) for r in regs])[1]
            color = lr.color
        nreg.color = color

        return nreg


    def find_sym_regions ( self, csyms, task=None ) :

        import Matrix
        import _contour
        import numpy

        centers, syms = csyms
        print "Finding %d-symmetry region groups" % len(syms)

        com = centers[0]
        t_0_com = ( (1.0,0.0,0.0,-com[0]),
                    (0.0,1.0,0.0,-com[1]),
                    (0.0,0.0,1.0,-com[2]) )
        t_to_com = ( (1.0,0.0,0.0,com[0]),
                    (0.0,1.0,0.0,com[1]),
                    (0.0,0.0,1.0,com[2]) )

        rmask = numpy.ones(self.seg_map.data.size[::-1], numpy.uint32) * -1

        for reg in self.regions :
            points = reg.points()
            #i,j,k = points[:,0],points[:,1],points[:,2]
            #nmat[k-lk,j-lj,i-li] = dmat[k,j,i]
            for p in points :
                i,j,k = p
                rmask[k,j,i] = reg.rid

        # a map from region id to a symmetry id
        self.rid_sid = {}
        sid_at = 0

        # sorting the regions from largest to smallest
        # should avoid small regions (prob. noise) from claiming
        # larger regions as their symmetric counterparts
        # since we only process each region once (for efficiency
        # and to avoid claim wars) this would prevent
        # larger regions from finding corresponding large regions
        # that are true symmetric counterparts

        print " - sorting regions by size..."
        task.updateStatus ( 'Sorting regions by size' )
        size_regs = []
        for reg in self.regions : size_regs.append ( [ len(reg.points()), reg] )
        size_regs.sort ()
        size_regs.reverse ()

        for ri, npoints_r in enumerate ( size_regs ) :

            npoints, r = npoints_r

            if ri % 100 == 0 :
                task.updateStatus('Finding symmetric regions %.1f%%' % (
                    100.0 * float(ri) / float(len(self.regions))) )

            try :
                sid = self.rid_sid [ r.rid ]
                continue
            except : pass

            sid_at += 1
            self.rid_sid [ r.rid ] = sid_at

            for si, smat in enumerate ( syms [1 : ] ) :

                # transform the region points using symmetry matrix
                rpoints = r.points ().astype ( numpy.float32 )
                tf = Matrix.multiply_matrices( t_to_com, smat, t_0_com )
                _contour.affine_transform_vertices ( rpoints, tf )

                # and get region ids from the map at the transformed points
                ipoints = numpy.round (rpoints).astype(numpy.int)
                sym_rids = rmask [ ipoints[:,2],ipoints[:,1],ipoints[:,0] ]

                # make a map from sym_rid, a region id that shows up
                # at the trasnformed map indices for the current region
                # to the number of time it appears at transformed indices
                rm = {}
                for sym_rid in sym_rids :
                    try : rm [ sym_rid ] = rm [ sym_rid ] + 1
                    except : rm [ sym_rid ] = 1

                # ignore rids of -1, which are points
                # in the grid where there is no region id
                if -1 in rm : rm.pop ( -1 )

                if 0 :
                    # ignore regions that already have a sid
                    # not sure if this is a good idea yet
                    rm2 = {}
                    for sym_rid, count in rm.iteritems() :
                        if sym_rid in self.rid_sid : continue
                        else : rm2[sym_rid] = count
                    rm = rm2


                unique_sym_rids = rm.keys ()
                if len ( unique_sym_rids ) == 0 : continue

                # take the rid that appears the most as the symmetric counterpart
                sym_rid = unique_sym_rids [ numpy.argmax ( rm.values() ) ]

                self.rid_sid [ sym_rid ] = sid_at
                sreg = self.id_to_region [ sym_rid ]
                sreg.color = r.color

            #if r.rid > 30 : break

        print "%d regions, sid_at: %d - # sids: %d" % (
            len(self.regions), sid_at, len(set(self.rid_sid.values())) )




    def find_sym_regions2 ( self, csyms, regs = None, task=None ) :

        import Matrix
        import _contour
        import numpy

        centers, syms = csyms
        print "Finding %d-symmetry region groups" % len(syms)

        com = centers[0]
        t_0_com = ( (1.0,0.0,0.0,-com[0]),
                    (0.0,1.0,0.0,-com[1]),
                    (0.0,0.0,1.0,-com[2]) )
        t_to_com = ( (1.0,0.0,0.0,com[0]),
                    (0.0,1.0,0.0,com[1]),
                    (0.0,0.0,1.0,com[2]) )


        if task : task.updateStatus('Coloring symmetries - making map with region ids...')

        from time import clock
        t0 = clock()
        print " - making map with region ids"
        rmask = numpy.ones(self.seg_map.data.size[::-1], numpy.uint32) * -1
        t1 = clock()
        print "   - time %.2f sec" % (t1-t0)


        if task : task.updateStatus('Coloring symmetries - setting region ids...')

        t0 = clock()
        print " - setting region ids"
        for reg in self.regions :
            points = reg.points()
            #i,j,k = points[:,0],points[:,1],points[:,2]
            #nmat[k-lk,j-lj,i-li] = dmat[k,j,i]
            for p in points :
                i,j,k = p
                rmask[k,j,i] = reg.rid

        same_count = 0
        t1 = clock()
        print "   - time %.2f sec" % (t1-t0)

        #from CGLutil.AdaptiveTree import AdaptiveTree
        #tree = AdaptiveTree ( scpoints.tolist(), scpoints.tolist(), 4.0)

        log = 1

        if regs == None :
            regs = self.regions
            log = 0

        for ri, reg in enumerate ( regs ) :

            if ri % 10 == 0 and task :
                task.updateStatus('Coloring symmetries - finding symmetric regions: %.1f%%' % (
                    100.0 * float(ri) / float(len(self.regions))) )

            reg.sym_regs = {}

            for si, smat in enumerate ( syms ) :

                if si == 0 : continue


                # transform the region points using symmetry matrix
                rpoints = reg.points ().astype ( numpy.float32 )
                tf = Matrix.multiply_matrices( t_to_com, smat, t_0_com )
                _contour.affine_transform_vertices ( rpoints, tf )

                # and get region ids from the map at the transformed points
                ipoints = numpy.round (rpoints).astype(numpy.int)
                sym_rids = rmask [ ipoints[:,2],ipoints[:,1],ipoints[:,0] ]

                # make a map from sym_rid, a region id that shows up
                # at the transformed map indices for the current region
                # to the number of times it appears at transformed indices

                #print " - reg %d" % reg.rid
                #print sym_rids

                rm = {}
                totN = 0
                for sym_rid in sym_rids :
                    if sym_rid == -1 :
                        continue
                    if sym_rid == reg.rid :
                        same_count += 1
                        continue
                    try :
                        rm [ sym_rid ] = rm [ sym_rid ] + 1
                    except :
                        rm [ sym_rid ] = 1
                    totN += 1

                if totN > 0 :
                    reg.sym_regs[si] = []
                    for srid, num in rm.iteritems() :
                        if float(num) > float(totN) * 0.1 :
                            sreg = self.id_to_region [ srid ]
                            reg.sym_regs[si].append ( sreg )


            # update color from surface, could have been changed by user...
            if reg.surface_piece != None :
                if reg.surface_piece.vertexColors is not None :
                    #clr = r.surface_piece.vertexColors[0]
                    reg.surface_piece.vertexColors = None
                clr = reg.surface_piece.color
                reg.set_color ( clr )


            if log : print " - reg %d: %d sym regs of %d syms" % ( reg.rid, len(reg.sym_regs.keys()), len(syms) )
            for si, sregs in reg.sym_regs.iteritems() :
                if log : print "  -- sym %d - %d regs:" % (si, len(sregs)),
                for sreg in sregs :
                    if log : print sreg.rid,
                    if sreg.surface_piece == None :
                        continue
                    elif sreg.surface_piece.vertexColors is not None :
                        #print " - v color 0 : ", r.surface_piece.vertexColors[0]
                        #clr = r.surface_piece.vertexColors[0]
                        sreg.surface_piece.vertexColors = None

                    sreg.set_color ( reg.color )
                if log : print ""


            #if r.rid > 30 : break
        print " - done, %d same" % same_count
        t2 = clock()
        print "   - time %.2f sec" % (t2-t1)




    def calculate_watershed_regions ( self, mm, thrD, csyms=None, task = None ) :

        self.remove_all_regions()

        self.adj_graph = None

        m = mm.data.full_matrix()
        import _segment
        if timing: t0 = clock()
        if task:
            task.updateStatus('Computing watershed regions for %s' % mm.name)
        from numpy import zeros, uint32
        self.mask = zeros(m.shape, uint32)
        _segment.watershed_regions(m, thrD, self.mask)
        self.map_level = thrD
        if timing: t1 = clock()
        max_points, max_values = _segment.region_maxima(self.mask, m)
        if timing: t2 = clock()
        n = len(max_points)
        regions = []
        for r in range(n):
            if task and r % 1000 == 0:
                task.updateStatus('Created regions %d of %d' % (r, n))
            reg = Region ( self, r+1, max_points[r] )
            regions.append(reg)

        self.regions = set(regions)
        self.id_to_region = dict([(r.rid, r) for r in regions])
        if timing: t3 = clock()

        np = sum([r.point_count() for r in regions])
        print 'Calculated %d watershed regions covering %d grid points' % (len(regions), np)
        if timing:
            print 'Time %.2f: watershed mask %.2f, maxima %.2f, region objects %.2f' % (t3-t0, t1-t0, t2-t1, t3-t2)

        if csyms :
            self.find_sym_regions ( csyms, task )

        return regions


    def smooth_and_group ( self, steps, sdev, min_reg = 1, csyms=None, task=None ) :

        if timing: t0 = clock()

        dmap = self.volume_data()

        from numpy import single as floatc
        sm_mat = dmap.data.full_matrix().astype(floatc)

        rlist = []
        for iti in range ( steps ) :

            if task:
                task.updateStatus('Smoothing step %d of %d' % (iti + 1, steps))
            slev = (iti+1)*sdev

            if timing: t1 = clock()
            ijk_sdev = (slev, slev, slev)
            from VolumeFilter import gaussian
            sm = gaussian.gaussian_convolution (sm_mat, ijk_sdev, task = task)
            if timing: t2 = clock()

            rlist = None
            if csyms :
                rlist = self.group_by_tracking_maxima_sym ( sm, csyms, task )
            else :
                rlist = self.group_by_tracking_maxima ( sm, task )

            for r in rlist:
                r.smoothing_level = slev

            if timing: t3 = clock()
            num_regs = len(self.regions)

            print "Smoothing width %d voxels, %d regions" % (slev, num_regs)
            if timing:
                print 'Time %.2f: smooth %.2f, group %.2f' % (t3-t1, t2-t1, t3-t2)

            if num_regs <= min_reg : break

        self.smoothing_level = steps  * sdev

        print "--------------- sdev: %.1f ---------------" % self.smoothing_level

        return rlist


    def group_by_tracking_maxima ( self, m, task = None ) :

        rlist = list(self.regions)
        pos = numpy.array([r.max_point for r in rlist], numpy.intc)

        import _segment
        _segment.find_local_maxima(m, pos)

        rm = {}
        for reg, pt in zip(rlist, pos):
            pt = tuple(pt)
            if pt in rm:
                rm[pt].append(reg)
            else:
                rm[pt] = [reg]
        groups = [(pt,regs) for pt,regs in rm.items() if len(regs) >= 2]

        rlist = []
        for i, (pt,regs) in enumerate(groups):
            if task and i % 100 == 0:
                s = "Grouping %d of %d" % ( i, len( groups ) )
                task.updateStatus ( s )
            newReg = self.join_regions ( regs, pt )
            rlist.append(newReg)

        return rlist



    def group_by_tracking_maxima_sym ( self, m, csyms, task = None ) :

        rlist = list(self.regions)
        pos = numpy.array([r.max_point for r in rlist], numpy.intc)

        import _segment
        _segment.find_local_maxima(m, pos)

        rm = {}
        for reg, pt in zip(rlist, pos):
            pt = tuple(pt)
            if pt in rm:
                rm[pt].append(reg)
            else:
                rm[pt] = [reg]
        groups = [(pt,regs) for pt,regs in rm.items() if len(regs) >= 2]

        rlist = []
        for i, (pt,regs) in enumerate(groups):
            if task and i % 100 == 0:
                s = "Grouping %d of %d" % ( i, len( groups ) )
                task.updateStatus ( s )

            doJoin = True

            # if any two (or more) regions have the same sid, they should
            # not be joined since that would violate the symmetry condition

            sid_rids = {}
            for r in regs :
                try : sid = self.rid_sid [ r.rid ]
                except : continue
                if sid in sid_rids :
                    # same sid seen more than once
                    doJoin = False
                try : sid_rids [ sid ].append ( r.rid )
                except : sid_rids [ sid ] = [ r.rid ]

            if doJoin :
                newReg = self.join_regions ( regs, pt )
                rlist.append(newReg)

            else :
                print "Trying to join: ", sid_rids
                clr = random_color()
                for sid, rids in sid_rids.iteritems () :
                    for rid in rids :
                        reg = self.id_to_region [ rid ]
                        reg.color = clr

        self.find_sym_regions ( csyms, task )
        return rlist


    def group_connected ( self, regions, min_contact ) :

        # Limit connection map to specified regions.
        cons = group_contacts(self.region_contacts(), regions)

        # Dictionary mapping region to set of connected regions.
        conr = ConnectedRegions ( regions, cons, min_contact )

        # List of sets of connected regions.
        rsets = dict([(id(cr), cr) for r, cr in conr.items()]).values()

        nregs = []
        for rset in rsets :
            if len(rset) > 1 :
                r = self.join_regions ( tuple(rset) )
                nregs.append(r)

        print "Created %d connected regions" % len(nregs)


    def group_connected_n ( self, nsteps, stopAt = 1, regions=None, csyms = None, task = None ) :

        nregs0 = len(self.regions)
        print "Grouping connected - %d regions" % nregs0

        newRegs, delRegs = [], []

        for si in range (nsteps) :

            cons = group_contacts(self.region_contacts(), limit_regions=regions)

            if task :
                task.updateStatus( 'Making contacts list' )

            minN, maxN = 1e5, 0
            clist = []
            for r1, r1cons in cons.iteritems() :
                for r2, con in r1cons.iteritems () :
                    clist.append ( [r1, r2, con] )
                    if con.N < minN : minN = con.N
                    if con.N > maxN : maxN = con.N

            if len(clist) == 0 :
                print " - no more connections, stopping"
                break

            print " - %d cons, N %d - %d, sorting..." % (len(clist), minN, maxN)

            if task :
                task.updateStatus( 'Sorting contacts list' )

            clist.sort ( reverse=True, key=lambda x: x[2].N )


            #rgrouped = {}
            nregs = []
            for r1, r2, con in clist :
                #if not r1 in rgrouped and not r2 in rgrouped : #r1.preg
                if not r1.preg and not r2.preg :
                    r = self.join_regions ( (r1, r2) )
                    #rgrouped[r1], rgrouped[r2] = 1, 1
                    nregs.append(r)
                    if si == 0 : delRegs.extend ( [r1,r2] )

            newRegs = nregs[:]

            print " - connected %d regions, now at %d" % ( len(nregs), len(self.regions) )

            if len(nregs) <= stopAt :
                print " - stopping for %d regions" % stopAt
                break

        return newRegs, delRegs


    def ungroup_regions ( self, regs, task = None ):

        rlist = []
        newRegs = []
        removeRegs = []
        for r in regs :
            if len(r.cregs) == 0 :
                rlist.append(r)
            else :
                rlist.extend(r.cregs)
                newRegs.extend(r.cregs)
                removeRegs.append(r)

        self.remove_regions(removeRegs, update_surfaces = False, task = task)

        print 'Ungrouped %d regions into %d regions' % ( len(regs), len(rlist) )
        return [newRegs, removeRegs]


    def display_regions(self, style = 'Voxel_Surfaces',
                        max_reg = None, task = None,
                        bForce=False):

        dmap = None
        init_mat = None

        print " -- showing ?", len(self.regions), "? regions"

        rlist = list(self.regions)
        rlist.sort(lambda r1,r2: cmp(r2.point_count(),r1.point_count()))

        self.style = style

        for i, reg in enumerate(rlist):

            if not max_reg is None and i >= max_reg :
                reg.remove_surface()
                continue

            if (not bForce) and reg.surface_piece:
                reg.surface_piece.display = True
                continue

            vt = None

            if style == 'Voxel_Surfaces' :

                reg.make_surface (None,None,self.regions_scale,bForce)

            elif style == 'Density_Maxima' :
                rpts = numpy.array ( [r.max_point for r in r.childless_regions()], numpy.float32 )
                import MultiScale
                vt = MultiScale.surface.surface_points ( rpts, 1.0, 0.1, .25, 5 )
                import _contour

                _contour.affine_transform_vertices ( vt[0], self.point_transform() )
                reg.make_surface ( vt[0], vt[1] )

            elif style == 'Iso_Surfaces' :

                if dmap == None:
                    v = self.volume_data()
                    dmap = v.writable_copy ( require_copy=True )
                    init_mat = dmap.data.full_matrix().copy()

                tpoints = reg.map_points()
                vt = None
                print " - reg %d, %d voxels" % (reg.rid, len(tpoints))
                import VolumeData
                sg = VolumeData.zone_masked_grid_data( dmap.data, tpoints, dmap.data.step[0] )
                m = sg.full_matrix()

                cmatrix = dmap.data.full_matrix()
                cmatrix[:,:,:] = m[:,:,:]
                dmap.surface_levels = [ dmap.surface_levels[0] ]
                dmap.region = ( dmap.region[0], dmap.region[1], [1,1,1] )
                from VolumeViewer.volume import Rendering_Options
                ro = Rendering_Options()
                #ro.surface_smoothing = True
                #ro.smoothing_factor = .2
                #ro.smoothing_iterations = 5
                dmap.update_surface ( False, ro )

                surf_sp = None
                surf_sp0 = None
                for sp in dmap.surfacePieces :
                    v, t = sp.geometry
                    #print "- %d vertices, %d tris" % (len(v), len(t))
                    if len(v) == 8 and len(t) == 12 :
                        continue
                    if len(v) == 0 and len(t) == 0 :
                        surf_sp0 = sp
                    else :
                        surf_sp = sp

                try : vt = surf_sp.geometry
                except : vt = surf_sp0.geometry

                cmatrix[:,:,:] = init_mat[:,:,:]
                dmap.data.values_changed()
                reg.make_surface ( vt[0], vt[1] )


            if task and i % 20 == 0:
                task.updateStatus('Making surface for region %d of %d'
                                  % (i, len(rlist)))


        if dmap : dmap.close()


    def color_density(self, regions = None):

        d = self.volume_data()
        if d is None:
            from segment_dialog import umsg
            #segdlg = volume_segmentation_dialog ()
            #if segdlg :
            umsg ("No map - select segmentation & map, then File -> Associate")
            return

        if not regions is None and len(regions) == 0:
            for p in d.surface_piece_list:
                p.vertexColors = None
            return

        m = self.mask
        if m is None:
            return

        ijk_to_xyz = self.point_transform()
        if ijk_to_xyz is None:
            return

        cmap = self.region_colors(regions)
        offset = max(d.data.step)

        color_surface_pieces(d.surface_piece_list, m, ijk_to_xyz, offset, cmap)

        color_solid(d, m, cmap)




    def region_colors(self, regions = None):

        if regions is None:
            regions = self.regions

        nid = self.max_region_id + 1
        from numpy import ones, float32
        colors = ones((nid,4), float32)
        for r in regions:
            for c in r.all_regions():
                clr = ''
                if r.surface_piece == None :
                    continue
                elif r.surface_piece.vertexColors is not None :
                    #print " - v color 0 : ", r.surface_piece.vertexColors[0]
                    clr = r.surface_piece.vertexColors[0]
                else :
                    #print " - s color : ", r.surface_piece.color
                    clr = r.surface_piece.color
                colors[c.rid] = clr
                r.color = clr
        return colors



    def all_regions ( self ):

        return self.id_to_region.values()

    def childless_regions ( self ):

        rlist = [r for r in self.id_to_region.values() if len(r.cregs) == 0]
        return rlist

    def selected_regions ( self ) :

        import Surface
        sregs = self.regions
        regions = [p.region for p in Surface.selected_surface_pieces()
                   if hasattr(p, 'region') and p.region.segmentation is self]
        return regions


    def visible_regions ( self ) :

        import Surface
        sregs = self.regions
        vis = []
        for r in sregs :
            if r.visible() :
                vis.append ( r )
        return vis


    def grouped_regions ( self ) :

        import Surface
        #sregs = self.regions
        regions = [r for r in self.regions if len(r.cregs) > 0]
        return regions


    def change_surface_resolution(self, res, task = None):

        if res == self.surface_resolution:
            return
        self.surface_resolution = res
        for i, r in enumerate(self.regions):
            if task and i % 100 == 0:
                task.updateStatus('Redisplaying surface %d of %d' %
                                  (i, len(self.regions)))
            sp = r.surface_piece
            if sp:
                display = sp.display
                color = sp.color
                r.remove_surface()
                sp = r.make_surface()
                sp.display = display
                sp.color = color

    def close(self):

        from chimera import openModels as om
        om.close ( [self] )

        if self.adj_graph :
            self.adj_graph.close()

class Region:

    def __init__( self, segmentation, rid, max_point = None, children = None ) :

        self.segmentation = segmentation
        self.rid = rid
        segmentation.id_to_region[rid] = self
        segmentation.regions.add(self)
        segmentation.max_region_id = max(rid, segmentation.max_region_id)
        if children is None:
            children = []  # Don't use [] as default since appends change it.
        self.mask_id = None if children else rid
        self.cregs = children   # Child regions
        self.preg = None        # Parent region
        self.smoothing_level = 0

        self.rbounds = None      # (imin, jmin, kmin), (imax, jmax, kmax)
        self.npoints = None
        self.max_point = None     # Position of local maximum
        self.surface_piece = None       # Displayed surface.

        self.color = random_color()         # Surface color
        self.placed = False

        if children:
            for reg in children :
                reg.preg = self
                segmentation.regions.remove(reg)
        self.max_point = max_point

    def bounds ( self ):

        if self.rbounds is None:
            if self.cregs:
                self.rbounds = union_bounds([r.bounds() for r in self.cregs])
            elif self.mask_id is None:
                self.rbounds = ((1,1,1),(0,0,0))        # Empty group
            else:
                self.segmentation.calculate_region_bounds()
        return self.rbounds

    def edge_distance ( self ):

        bmin, bmax = self.bounds()
        s = self.segmentation.grid_size()
        d = min(min(bmin), min([s[a]-bmax[a]-1 for a in (0,1,2)]))
        return d

    def points ( self ):

        if self.cregs :
            cpts = [r.points() for r in self.childless_regions()]
            p = numpy.concatenate(cpts) if cpts else None
        else :
            # would be useful to know what this is doing,
            # it's a little hard to understand from the code...
            (imin, jmin, kmin), (imax, jmax, kmax) = self.bounds()
            m = self.segmentation.mask[kmin:kmax+1,jmin:jmax+1,imin:imax+1]
            import _segment
            p = _segment.region_points(m, self.rid)
            p[:,0] += imin
            p[:,1] += jmin
            p[:,2] += kmin
        return p

    def map_points ( self ) :

        tpoints = numpy.array(self.points(), numpy.float32)
        tf = self.segmentation.point_transform()
        import _contour
        _contour.affine_transform_vertices ( tpoints, tf )
        return tpoints


    def center_of_points ( self, transform = True ) :

        plists = [r.points() for r in self.childless_regions()]
        s = numpy.sum([numpy.sum(plist, axis=0) for plist in plists], axis=0)
        com = s.astype(numpy.float32)
        com /= self.point_count()
        if transform:
            tf = self.segmentation.point_transform()
            import _contour
            _contour.affine_transform_vertices ( com.reshape((1,3)), tf )

        return com


    def region_radius ( self ) :

        from numpy import sqrt, max, sum, square
        rmax = 0
        c = self.center_of_points()
        plists = [r.map_points() for r in self.childless_regions()]
        for p in plists:
            r = sqrt ( max ( sum ( square (p - c), 1 ) ) )
            rmax = max(rmax, r)
        return rmax


    def enclosed_volume ( self ) :

        return self.point_count() * self.segmentation.voxel_volume()


    def set_color ( self, color ):

        self.color = color
        if self.surface_piece:
            self.surface_piece.color = color

    def surface ( self ) :

        return self.surface_piece

    def has_surface ( self ) :

        sp = self.surface_piece
        return sp and not sp.__destroyed__

    # Does this region or any parent have a displayed surface?
    def visible ( self ) :

        sp = self.surface_piece
        shown = sp and not sp.__destroyed__ and sp.display
        if not shown and self.preg:
            shown = self.preg.visible()
        return shown

    def show_surface ( self ):

        sp = self.surface_piece
        if sp:
            sp.display = True

    def hide_surface ( self ):

        sp = self.surface_piece
        if sp:
            sp.display = False

    def make_surface ( self, vertices = None, triangles = None, scale=1.0, bForce=False ):

        if (not bForce) and self.surface_piece:
            return self.surface_piece

        self.remove_surface(including_children = True)

        if vertices is None:
            seg = self.segmentation
            from MultiScale.surface import surface_points
            vertices, triangles, normals = \
                surface_points ( self.points(),
                                 resolution = seg.surface_resolution,
                                 density_threshold = 0.1,
                                 smoothing_factor = .25,
                                 smoothing_iterations = 5 )
            tf = seg.point_transform()

            import numpy
            if numpy.fabs(scale-1.0) > 0.01 :
                com = self.center_of_points ()
                t_0_com = ( (1.0,0.0,0.0,-com[0]),
                            (0.0,1.0,0.0,-com[1]),
                            (0.0,0.0,1.0,-com[2]) )
                t_to_com = ( (1.0,0.0,0.0,com[0]),
                            (0.0,1.0,0.0,com[1]),
                            (0.0,0.0,1.0,com[2]) )
                t_scale = ( (scale,0.0,0.0,0.0),
                            (0.0,scale,0.0,0.0),
                            (0.0,0.0,scale,0.0) )
                import Matrix
                tf = Matrix.multiply_matrices( t_to_com, t_scale, t_0_com, tf )

            import _contour
            _contour.affine_transform_vertices ( vertices, tf )

        rgba = self.top_parent().color
        nsp = self.segmentation.addPiece ( vertices, triangles, rgba )
        # print " - added piece with color", nsp.color

        try : sidstr = ", sym # %d" % self.segmentation.rid_sid [ self.rid ]
        except : sidstr = ""

        nsp.oslName = str(self.rid) + sidstr
        nsp.region = self
        self.surface_piece = nsp
        return nsp

    def remove_surface ( self, including_children = False ):

        p = self.surface_piece
        if p:
            p.model.removePiece ( p )
        self.surface_piece = None

        if including_children:
            for c in self.cregs:
                c.remove_surface(including_children)

    # Contacting top regions.
    def contacting_regions ( self ):

        rset = set()
        rcons = self.segmentation.region_contacts()
        for r in self.childless_regions():
            if r in rcons:
                rset.update(rcons[r].keys())
        cr = set()
        for r in rset:
            cr.add(r.top_parent())
        return tuple(cr)

    def parents ( self ) :

        rlist = []
        r = self
        while r.preg :
            r = r.preg
            rlist.append(r)
        return rlist

    def top_parent ( self ) :

        r = self
        while r.preg :
            r = r.preg
        return r

    def all_regions ( self ):

        rlist = [self]
        for c in self.cregs:
            rlist.extend(c.all_regions())
        return rlist

    def in_group ( self ):

        return (not self.preg is None) or self.has_children()

    def has_children ( self ):

        return len(self.cregs) > 0

    def is_group ( self ):

        return len(self.cregs) > 0

    def children ( self ):

        return self.cregs

    def all_children ( self ):

        rlist = list(self.cregs)
        for c in self.cregs:
            rlist.extend(c.all_children())
        return rlist

    def childless_regions ( self ):

        if self.cregs:
            rlist = []
            for c in self.cregs:
                rlist.extend(c.childless_regions())
        else:
            rlist = [self]
        return rlist

    def point_count ( self ):

        if self.npoints is None:
            if self.cregs:
                self.npoints = sum([r.point_count() for r in self.cregs])
            elif self.mask_id is None:
                self.npoints = 0        # Empty group region
            else:
                self.segmentation.calculate_region_bounds()
                if self.npoints is None:
                    print 'empty region', self.rid
                    self.npoints = 0
        return self.npoints

    def remove_children ( self, cregs, update_surfaces = True):

        for c in cregs:
            self.cregs.remove(c)
        self.children_changed(update_surfaces)

    def children_changed(self, update_surfaces = True):

        # Clear cached point count.
        self.npoints = None
        self.rbounds = None

        # Update displayed surface.
        if self.surface_piece:
            self.remove_surface()
            if update_surfaces and len(self.cregs) > 0:
                self.make_surface()

        # Update parent
        p = self.preg
        if p:
            p.children_changed(update_surfaces)

    def has_attribute(self, name):

        return hasattr(self, 'attrib') and name in self.attrib

    def get_attribute(self, name, default = None):

        return self.attrib[name] if hasattr(self, 'attrib') and name in self.attrib else default

    def set_attribute(self, name, value):

        if self.is_reserved_name(name):
            return False                 # Name clash
        if not hasattr(self, 'attrib'):
            self.attrib = {}
        self.attrib[name] = value
        setattr(self, name, value)
        return True

    def is_reserved_name(self, name):

        return hasattr(self, name) and not (hasattr(self, 'attrib') and name in self.attrib)

    def remove_attribute(self, name):

        if hasattr(self, 'attrib') and name in self.attrib:
          del self.attrib[name]
          delattr(self, name)

    def attributes(self):

        return self.attrib if hasattr(self, 'attrib') else {}


class Contact:

    def __init__ (self, ncontact):

        self.N = ncontact               # number of voxel pairs
        self.D = 0.0                    # divide by 2*N to get average
        self.maximum_density = None


def group_contacts ( rcons, limit_regions = None, task = None ) :

    if not limit_regions is None:
        limit_regions = set(limit_regions)

    # Compute top parents of leaf regions.
    top = dict([(r,r.top_parent()) for r in rcons])

    # Compute contacts of top region groups from leaf contacts.
    cons = {}
    for i, (r1, r1cons) in enumerate(rcons.items()) :
        r1_parent = top[r1]
        if limit_regions is None or r1_parent in limit_regions:
            if task and i % 1000 == 0:
                task.updateStatus('group contacts %d of %d' % (i, len(rcons)))
            r1pcons = cons.setdefault(r1_parent, {})
            for r2, o in r1cons.iteritems() :
                r2_parent = top[r2]
                if r2_parent != r1_parent:
                    if limit_regions is None or r2_parent in limit_regions:
                        if r2_parent in r1pcons:
                            no = r1pcons[r2_parent]
                            no.D += o.D
                            no.N += o.N
                            if (no.maximum_density == None or
                                no.maximum_density < o.maximum_density) :
                                no.maximum_density = o.maximum_density
                        else :
                            r1pcons[r2_parent] = Contact(o.N)

    return cons



def GroupedRegions ( regions, rset = None ):

    if rset is None:
        rset = set ( regions )
    else:
        rset.update( regions )

    for r in regions:
        GroupedRegions ( r.cregs, rset )

    return rset


def SmallRegions ( regions, minRegionSize ) :

    sreg = [r for r in regions if r.point_count() < minRegionSize ]
    return sreg


def SmallConnectedRegions ( regions, rcons, minRegionSize, task = None ) :

    conr = ConnectedRegions ( regions, rcons, task=task )
    csize = {}
    for rset in conr.values():
        s = id(rset)
        if not s in csize:
            csize[s] = sum([c.point_count() for c in rset])
    sireg = [r for r in regions if csize[id(conr[r])] < minRegionSize]
    return sireg


# Calculate dictionary mapping region to set of connected regions.
def ConnectedRegions ( regions, rcons, min_contact = None, task = None ) :

    conr = {}
    cc = 0
    for r in regions:
        if task and cc % 100 == 0:
            task.updateStatus('%d connected region sets' % cc)
        if not r in conr:
            conr[r] = s = set([r])
            bndry = [r]
            cc += 1
            while bndry :
                reg = bndry.pop()
                if reg in rcons:
                    for cr in rcons[reg]:
                        if (min_contact is None or
                            rcons[reg][cr].N >= min_contact):
                            if not cr in conr:
                                s.add(cr)
                                conr[cr] = s
                                bndry.append(cr)
    return conr


def TopParentRegions ( regions ) :

    parents = set()
    for r in regions:
        parents.add ( r.top_parent () )
    return parents

def all_regions ( regions ) :

    rlist = []
    for r in regions:
        rlist.extend(r.all_regions())
    return rlist


def region_bounds(regions):

    rlist = []
    for r in regions:
        rlist.extend(r.childless_regions())
    bounds = [r.bounds() for r in rlist]
    ijk_min = numpy.min([b[0] for b in bounds], axis=0)
    ijk_max = numpy.max([b[1] for b in bounds], axis=0)
    return ijk_min, ijk_max


def random_color(avoid_rgba = None, minimum_rgba_distance = 0.2):

    from random import random as rand
    while True:
        c = ( 0.5*(1+rand()), 0.5*(1+rand()), 0.5*(1+rand()), 1.0 )
        if avoid_rgba is None:
            break
        d = sum([(x0-x1)*(x0-x1) for x0,x1 in zip(c, avoid_rgba)])
        if d > minimum_rgba_distance * minimum_rgba_distance:
            break
    return c


def segmentations():

    from chimera import openModels as om
    slist = om.list(modelTypes = [Segmentation])
    return slist


#
# Group first regions with the smallest drop from one region maximum to the
# contact interface maximum.
#
def group_by_contacts(smod, task = None):

    map = smod.volume_data().full_matrix()

    if task:
        task.updateStatus('Computing interface maxima')
    import _segment
    ci, cf = _segment.interface_values(smod.mask, map)

    if task:
        task.updateStatus('Computing drops')

    for r in smod.all_regions():
        r.color_level = len(smod.regions)

    # Cache top parent, optimization for high grouping depth.
    tp = {}
    for t in smod.regions:
        for c in t.childless_regions():
            tp[c] = t

    id2r = smod.id_to_region
    rc = [(id2r[id1],id2r[id2],dmax) for (id1,id2,nv),(dmax,dsum) in zip(ci,cf)]

    rmax = region_maxima(smod.regions, map)
    drc = [(min(rmax[tp[r1]],rmax[tp[r2]])-cd,r1,r2,cd) for r1,r2,cd in rc]
    drc.sort()
    drc.reverse()

    cc = 0
    rcnt = 0
    drcr = []
    dr = None
    while drc:
        d, r1, r2, cd = drc.pop()
        dl = min(rmax[tp[r1]],rmax[tp[r2]])-cd
        if dl != d:
            # Drop changed due to region merging.
            dr = dl if dr is None else min(dl, dr)
            drcr.append((dl,r1,r2,cd))
            continue
        elif not dr is None and (d > dr or len(drc) == 0):
            # Resort drops.
            drcr.append((dl,r1,r2,cd))
            drc = [(min(rmax[tp[r1]],rmax[tp[r2]])-cd,r1,r2,cd)
                   for d,r1,r2,cd in drc + drcr]
            drc.sort()
            drc.reverse()
            drcr = []
            dr = None
            rcnt += 1
            continue
        cc += 1
        if task and cc % 100 == 0:
            task.updateStatus('Contact %d of %d' % (cc, len(ci)))
        p1 = collapse_chain(tp, r1)
        p2 = collapse_chain(tp, r2)
        if p1 != p2:
            rj = smod.join_regions((p1,p2))
            rj.color_level = len(smod.regions)
            tp[p1] = tp[p2] = rj
            rmax[rj] = max((rmax[p1], rmax[p2]))
    print 'max group depth', maximum_group_depth(smod.regions)
    print 'resorted', rcnt

def connected_subsets(connections):

    subset = {}
    pairs = []
    for e1, e2 in connections:
        s1 = subset.get(e1)
        if s1 is None:
            s1 = subset[e1] = [e1]
        s2 = subset.get(e2)
        if s2 is None:
            s2 = subset[e2] = [e2]
        if not s1 is s2:
            if len(s2) > len(s1):
                t = e1,s1; e1,s1 = e2,s2; e2,s2 = t     # Swap
            s1.extend(s2)
            for e in s2:
                subset[e] = s1
            pairs.append((e1, e2))
    subsets = dict([(id(s),s) for s in subset.values()]).values()
    return subsets, pairs

def region_maxima(regions, map, rmax = None):

    if rmax is None:
        rmax = {}
    for r in regions:
        if r.cregs:
            region_maxima([c for c in r.cregs if not c in rmax], map, rmax)
            rmax[r] = max([rmax[c] for c in r.cregs])
        else:
            rmax[r] = map[r.max_point[2],r.max_point[1],r.max_point[0]]
    return rmax

def maximum_group_depth(regions):

    dmax = 0
    rd = [(r,0) for r in regions]
    while rd:
        r,d = rd.pop()
        if r.cregs:
            rd.extend([(r,d+1) for r in r.cregs])
        elif d > dmax:
            dmax = d
    return dmax

def collapse_chain(d, e):
    v = d[e]
    if v in d and v != e:
        while v in d and d[v] != v:
            v = d[v]
        while e != v:
            n = d[e]
            d[e] = v
            e = n
    return v

def regions_radius ( regions ) :

    from numpy import sqrt, max, sum, square, zeros, float32
    s = zeros((3,), float32)
    n = 0
    for reg in regions:
        for r in reg.childless_regions():
            p = r.map_points()
            s += sum(p, axis=0)
            n += len(p)
    c = s / n
    rmax = 0
    for reg in regions:
        for r in reg.childless_regions():
            p = r.map_points()
            radius = sqrt ( max ( sum ( square (p - c), axis = 1 ) ) )
            if radius > rmax:
                rmax = radius
    return rmax


def mean_and_sd(regions, volume):

    if len(regions) == 0:
        return None

    smod = regions[0].segmentation
    mask_size = smod.grid_size()
    map_size = volume.data.size
    binned_mask = bin_size(mask_size, map_size)
    if binned_mask is None:
        print 'Incompatible mask (%d,%d,%d) and map (%d,%d,%d) sizes.' % (tuple(mask_size) + tuple(map_size))
        return None

    vmat = volume.full_matrix()

    moff = []
    b0, b1, b2 = binned_mask
    for i in range(b0):
        for j in range(b1):
            for k in range(b2):
                moff.append(vmat[k::b2,j::b1,i::b0])

    means = []
    sdevs = []
    from numpy import concatenate
    for reg in regions :
        pts = reg.points()
        i,j,k = pts[:,0], pts[:,1], pts[:,2]
        mval = concatenate([m[k,j,i] for m in moff])
        mean, sd = mval.mean(), mval.std()
        means.append(mean)
        sdevs.append(sd)
        reg.set_attribute('map mean', mean)
        reg.set_attribute('map sd', sd)

    return means, sdevs


def mask_volume(regions, volume) :

    mgrid = mask_data(regions, volume)
    if mgrid == None :
        return None

    import VolumeViewer
    nv = VolumeViewer.volume_from_grid_data ( mgrid, show_data = False,
                                              show_dialog = False )
    nv.copy_settings_from(volume)
    nv.show()
    return nv


def mask_data(regions, volume) :

    if len(regions) == 0:
        return None

    smod = regions[0].segmentation
    mask_size = smod.grid_size()
    map_size = volume.data.size
    binned_mask = bin_size(mask_size, map_size)
    if binned_mask is None:
        return None

    vmat = volume.full_matrix()
    import numpy
    mmat = numpy.zeros_like ( vmat )

    mask_matrix(regions, vmat, mmat, binned_mask)

    import os.path
    name = os.path.splitext ( volume.name )[0] + "_masked"
    d = volume.data
    from VolumeData import Array_Grid_Data
    mgrid = Array_Grid_Data ( mmat, d.origin, d.step, d.cell_angles, name=name)

    #import VolumeViewer
    #nv = VolumeViewer.volume_from_grid_data ( mgrid, show_data = False,
    #                                          show_dialog = False )
    #nv.copy_settings_from(volume)
    #nv.show()
    return mgrid


def masked_matrix (regions, volume) :

    if len(regions) == 0:
        return None

    smod = regions[0].segmentation
    mask_size = smod.grid_size()
    map_size = volume.data.size
    binned_mask = bin_size(mask_size, map_size)
    if binned_mask is None:
        return None

    vmat = volume.full_matrix()
    import numpy
    mmat = numpy.zeros_like ( vmat )

    mask_matrix(regions, vmat, mmat, binned_mask)

    return mmat




def remove_mask_volume(regions, volume) :

    if len(regions) == 0:
        return None

    smod = regions[0].segmentation
    mask_size = smod.grid_size()
    map_size = volume.data.size
    binned_mask = bin_size(mask_size, map_size)
    if binned_mask is None:
        return None

    vmat = volume.full_matrix().copy()
    import numpy
    mmat = numpy.zeros_like ( vmat )

    mask_matrix(regions, mmat, vmat, binned_mask)

    import os.path
    name = os.path.splitext ( volume.name )[0] + "_imasked"
    d = volume.data
    from VolumeData import Array_Grid_Data
    mgrid = Array_Grid_Data ( vmat, d.origin, d.step, d.cell_angles, name=name)
    import VolumeViewer
    nv = VolumeViewer.volume_from_grid_data ( mgrid, show_data = False,
                                              show_dialog = False )
    nv.copy_settings_from(volume)
    nv.show()
    return nv




def mask_matrix(regions, from_matrix, to_matrix, binned_mask = (1,1,1)):

    if tuple(binned_mask) == (1,1,1):
        for reg in regions :
            pts = reg.points()
            i,j,k = pts[:,0], pts[:,1], pts[:,2]
            to_matrix[k,j,i] = from_matrix[k,j,i]
    else:
        moff = []
        b0, b1, b2 = binned_mask
        for i in range(b0):
            for j in range(b1):
                for k in range(b2):
                    moff.append((from_matrix[k::b2,j::b1,i::b0],
                                 to_matrix[k::b2,j::b1,i::b0]))
        for reg in regions :
            pts = reg.points()
            i,j,k = pts[:,0], pts[:,1], pts[:,2]
            for fm, tm in moff:
                tm[k,j,i] = fm[k,j,i]

def bin_size(mask_size, map_size):

    br = [(t/s, t%s) for s,t in zip(mask_size, map_size)]
    for b,r in br:
        if b < 1 or r >= b:
            return None
    bsize = [b for b,r in br]
    return bsize

def select_regions(regions, create_surfaces = True):

    surfs = [r.surface() for r in regions if r.has_surface()]
    if len(surfs) < len(regions) and create_surfaces:
        surfs += make_surfaces([r for r in regions if not r.has_surface()])

    from chimera import selection
    selection.setCurrent(surfs)

def SelectedRegions():

    from chimera import selection
    smods = [s for s in selection.currentGraphs() if isinstance(s,Segmentation)]
    sel = []
    for smod in smods:
        sel.extend(smod.selected_regions())
    return sel

def make_surfaces(regions, surfaces = None, task = None):

    if surfaces is None:
        surfaces = []
    if task:
        for i,r in enumerate(regions):
            if i%20 == 0:
                task.updateStatus('%d of %d' % (i+1,len(regions)))
            surfaces.append(r.make_surface())
    else:
        from chimera import tasks, CancelOperation
        task = tasks.Task('Making surfaces', modal = True)
        try:
            make_surfaces(regions, surfaces, task)
        except CancelOperation:
            pass
        finally:
            task.finished()

    return surfaces

def show_only_regions(rlist):

    segs = [r.segmentation for r in rlist]
    for s in segs:
        for p in s.surfacePieces:
            p.display = False
        s.display = True
    for r in rlist:
        r.show_surface()

def color_surface_pieces(plist, mask, ijk_to_xyz, offset, colormap):

    import Matrix
    xyz_to_ijk = Matrix.invert_matrix(ijk_to_xyz)

    for p in plist:
        v, t = p.geometry
        n = p.normals
        xyz = v - offset*n
        import _interpolate as i
        rid, outside = i.interpolate_volume_data(xyz, xyz_to_ijk, mask,
                                                 'nearest')
        rid = rid.astype(numpy.uint32)    # Interpolation gives float32
        colors = colormap[rid,:]
        p.vertexColors = colors

def color_solid(v, mask, colormap):

    icmap = numpy.empty(colormap.shape, numpy.uint8)
    import _volume
    _volume.colors_float_to_uint(colormap, icmap)

    # Swap red and blue.  colormap is RGBA but need BGRA
    r = icmap[:,0].copy()
    icmap[:,0] = icmap[:,2]
    icmap[:,2] = r

    def cm(colors, slice = None, v=v, mask=mask, cmap=icmap):
        s = list(v.matrix_slice())
        s.reverse()
        m = mask[s]
        ms = m[slice] if slice else m
        # Numpy is 20x slower than _volume C++ routine.
        # colors[...] *= cmap[mask]
        import _volume
        _volume.indices_to_colors(ms, cmap, colors, modulate = True)

    v.mask_colors = cm

    # Update display.
    if v.shown() and v.representation == 'solid':
        v.close_solid()
        v.show()

def union_bounds(bounds):

    b = numpy.array(bounds)
    ijk_min = [b[:,0,a].min() for a in (0,1,2)]
    ijk_max = [b[:,1,a].max() for a in (0,1,2)]
    return ijk_min, ijk_max

# -----------------------------------------------------------------------------
#
def boundary_regions(rset, rcons):

    bset = set()
    for r in rset:
        if r in rcons:
            for b in rcons[r]:
                if not b in rset:
                    bset.add(b)
    return bset

# -----------------------------------------------------------------------------
#
def boundary_groups(rset, rcons):

    bset = boundary_regions(rset, rcons)
    gset = set()
    for r in bset:
        gset.update(r.top_parent().all_regions())
    return gset

# -----------------------------------------------------------------------------
#
def childless_regions(regions):

    c = set()
    for r in regions:
        c.update(r.childless_regions())
    return c
