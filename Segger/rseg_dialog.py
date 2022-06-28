
# Copyright (c) 2020 Greg Pintilie - pintilie@mit.edu

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import chimera
import os
import os.path
import Tkinter
from CGLtk import Hybrid
import VolumeData
import _multiscale
import MultiScale.surface
import _surface
import numpy
import _contour
import Matrix
import VolumeViewer
from sys import stderr
from time import clock
import sets
import FitMap

from axes import prAxes
import regions
import graph
from Segger import showDevTools, timing, seggerVersion

OML = chimera.openModels.list

REG_OPACITY = 0.45


# http://geomalgorithms.com/a06-_intersect-2.html



from segment_dialog import current_segmentation, segmentation_map


class RSeg_Dialog ( chimera.baseDialog.ModelessDialog ):

    title = "rSeg - Radial Segmentation (Segger v" + seggerVersion + ")"
    name = "segger_rseg"
    buttons = ( "Close" )
    help = 'https://cryoem.slac.stanford.edu/ncmi/resources/software/segger'

    def fillInUI(self, parent):

        self.group_mouse_mode = None

        tw = parent.winfo_toplevel()
        self.toplevel_widget = tw
        tw.withdraw()

        parent.columnconfigure(0, weight = 1)

        row = 0

        menubar = Tkinter.Menu(parent, type = 'menubar', tearoff = False)
        tw.config(menu = menubar)

        f = Tkinter.Frame(parent)
        f.grid(column=0, row=row, sticky='ew')
        l = Tkinter.Label(f, text='  ')
        l.grid(column=0, row=row, sticky='w')



        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = "1. Select map in Segment Map dialog, press Segment button.", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)





        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = "2. Optional - for icosahedral (not round) shells: ", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = "     A. Use Tools -> Higher-Order Structure -> Icosahedron Surface.", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)





        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = "     B. Match Icosahedron to current map & segmentation.", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = "     C. From Icosahedron Surface ", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)

            b = Tkinter.Button(ff, text="Find Axes", command=self.Icos)
            b.grid (column=1, row=0, sticky='w', padx=5, pady=1)

        if showDevTools :

            b = Tkinter.Button(ff, text="Line CC", command=self.LineCC)
            b.grid (column=2, row=0, sticky='w', padx=5, pady=1)



        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = "3. Make histogram of distances from center of map to center of each region,", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = "    using", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)

            self.numBins = Tkinter.StringVar(ff)
            self.numBins.set ( "600" )
            e = Tkinter.Entry(ff, width=10, textvariable=self.numBins)
            e.grid(column=1, row=0, sticky='w', padx=5, pady=1)

            l = Tkinter.Label(ff, text = "bins", anchor = 'w')
            l.grid(column=2, row=0, sticky='ew', padx=5, pady=1)

            b = Tkinter.Button(ff, text="Make Histogram", command=self.MakeHist)
            b.grid (column=3, row=0, sticky='w', padx=5, pady=1)


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = "4. Plot histogram (e.g. using plot.ly), find distances with low values.", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = "5. Enter distances at which to separate regions, separated by commas:", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = "   ", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)

            self.segRads = Tkinter.StringVar(ff)

            if 0 or showDevTools :
                self.segRads.set ( "1006" )

            e = Tkinter.Entry(ff, width=40, textvariable=self.segRads)
            e.grid(column=1, row=0, sticky='w', padx=5, pady=1)


            b = Tkinter.Button(ff, text="Group", command=self.Segment)
            b.grid (column=2, row=0, sticky='ew', padx=5, pady=1)


        row += 1
        f = Tkinter.Frame(parent)
        f.grid(column=0, row=row, sticky='ew')
        l = Tkinter.Label(f, text='  ')
        l.grid(column=0, row=row, sticky='w')


        row += 1
        dummyFrame = Tkinter.Frame(parent, relief='groove', borderwidth=1)
        Tkinter.Frame(dummyFrame).pack()
        dummyFrame.grid(row=row,column=0,columnspan=7, pady=7, sticky='we')


        row = row + 1
        self.msg = Tkinter.Label(parent, width = 60, anchor = 'w', justify = 'left', fg="red")
        self.msg.grid(column=0, row=row, sticky='ew', padx=5, pady=1)
        row += 1


    def umsg ( self, txt ) :
        print txt
        self.status ( txt )

    def status ( self, txt ) :
        txt = txt.rstrip('\n')
        self.msg.configure(text = txt)
        self.msg.update_idletasks()




    def Icos ( self ) :

        imod = None
        axmod = None
        for m in chimera.openModels.list() :
            if m.name == "Icosahedron" :
                imod = m
            if m.name == "Icosahedron_Axes" :
                axmod = m


        if axmod == None :
            pass
        else :
            chimera.openModels.close ( [axmod] )


        if imod == None :
            self.umsg ( "No Icosahedron model found - please follow step 2." )
            return


        if len(imod.surfacePieces) <> 1 :
            self.umsg ( "Please set 'Subdivision factor' to 1" )
            return


        print len(imod.surfacePieces[0].geometry[1]), " tris"
        print len(imod.surfacePieces[0].geometry[0]), " verts"

        if len(imod.surfacePieces[0].geometry[1]) <> 20 :
            self.umsg ( "Please set 'Subdivision factor' to 1" )
            return


        self.umsg ( "Building axes..." )


        import _surface
        surf_mod = _surface.SurfaceModel()
        chimera.openModels.add([surf_mod], sameAs = imod)

        import axes; reload (axes)

        self.icos_vecs = []
        from numpy import arccos, pi



        for p in imod.surfacePieces :
            v, t = p.geometry[0], p.geometry[1]
            #print len(v), len(t)

            #for pt in v :
            #    print " - pt: ", pt

            for tri in t :
                #print " - tri: ", tri,
                p1 = v [ tri[0] ]
                p2 = v [ tri[1] ]
                p3 = v [ tri[2] ]
                mp = (p1 + p2 + p3) / 3.0
                pv = chimera.Vector ( mp[0], mp[1], mp[2] )
                r = pv.length
                pv.normalize()
                #print mp
                self.icos_vecs.append ( pv )

                cyl = axes.AddCylinderSolid ( chimera.Vector(0,0,0), pv, r, (.6,.4,.4,1), 10.0, surf_mod )
                cyl.name = "Icosahedron_Axes"

                p1v = chimera.Vector ( p1[0], p1[1], p1[2] ); p1v.normalize ()
                p2v = chimera.Vector ( p2[0], p2[1], p2[2] ); p2v.normalize ()
                p3v = chimera.Vector ( p3[0], p3[1], p3[2] ); p3v.normalize ()

                a1 = arccos ( p1v * pv ) * 180.0 / pi
                a2 = arccos ( p2v * pv ) * 180.0 / pi
                a3 = arccos ( p3v * pv ) * 180.0 / pi

                a12 = arccos ( p1v * p2v ) * 180.0 / pi

                # print a1, a2, a3, a12


        minAng = 1e9
        pv1 = self.icos_vecs[0]
        for pv2 in self.icos_vecs[1:] :
            dp = pv1 * pv2
            ang = arccos ( dp )
            #print ang * 180.0 / pi

        self.umsg ( "Axes built." )




    def MakeHist ( self ) :

        segMap = segmentation_map()
        if segMap == None :
            self.umsg ( "Please select a map in the Segment Map Dialog" )
            return

        import axes
        reload(axes)
        pts, weights = axes.map_points ( segMap )
        print len(pts)

        COM, U, S, V = prAxes ( pts )

        print " - COM : ", COM


        smod = current_segmentation ()
        if smod == None :
            self.umsg ( "Please select a Current Segmentation in the Segment Map dialog" )
            return

        print "Seg has %d regions" % (len(smod.regions))


        if hasattr(self, 'icos_vecs') :
            self.umsg ( "Making (icosahedrally corrected) histogram..." )
        else :
            self.umsg ( "Making histogram..." )

        nregs, last = len(smod.regions), 0
        regs = list(smod.regions)
        distByReg = {}
        for ri, r in enumerate ( regs ) :

            if 0 and r.surface_piece != None :
                if r.surface_piece.display == False :
                    print "i" + ri,
                    continue
            try :
                p = r.center_of_points ()
            except :
                print "+"
                continue

            rvec = chimera.Vector ( p[0], p[1], p[2] ) - chimera.Vector (COM[0], COM[1], COM[2])
            rad = 0.0

            if hasattr(self, 'icos_vecs') :
                for ivec in self.icos_vecs :
                    irad = ivec * rvec
                    if irad > rad :
                        rad = irad
            else :
                rad = rvec.length


            distByReg[r] = rad
            at = int(numpy.floor( 10.0 * (ri+1) / nregs ))
            if at > last :
                #print at,
                if hasattr(self, 'icos_vecs') :
                    self.status ( "Making (icosahedrally corrected) histogram %d regions, at %d" % (len(regs), ri+1) )
                else :
                    self.status ( "Making histogram %d regions, at %d" % (len(regs), ri+1) )
                last = at; at += 1

        print ""

        dists = distByReg.values ()
        maxDist = max (dists) + 0.01
        minDist = min (dists)
        nbins = int ( self.numBins.get() )
        dr = (maxDist - minDist) / float(nbins)
        print "%d dists - max %.2f, min %.2f, nb %d, dr %.2f" % (len(dists), maxDist, minDist, nbins, dr)

        bins = []
        for i in range (nbins) :
            bins.append ( [] )

        print "bad bins: ",
        for regm, rad in distByReg.iteritems() :
            bini = int ( numpy.floor ( (rad - minDist) / dr ) )
            if bini >= len(bins) :
                print bini,
                bini = len(bins)-1
            bins[bini].append (regm)

        print ""



        if 0 :
            f = open ( "rads.txt", "w" )
            for k,regs in enumerate ( bins ) :
                v = len(regs)
                vmin = minDist + k * dr
                vmax = minDist + (k+1) * dr
                rm = .5 * (vmin + vmax)
                vn = v / (4 * 3.14 * rm * rm)
                f.write ( "%d\t%.2f\t%.2f\t%d\t%f\n" % (k, vmin, vmax, v, vn) )
            f.close()

        self.distByReg = distByReg
        #print self.distByReg



        def save ( okay, dialog ):
            if okay:
                paths = dialog.getPaths ( )
                if paths:
                    path = paths[0]
                    self.umsg ( "Saved plot to: " + path )
                    f = open ( path, "w" )
                    for k,regs in enumerate ( bins ) :
                        v = len(regs)
                        vmin = minDist + k * dr
                        vmax = minDist + (k+1) * dr
                        rm = .5 * (vmin + vmax)
                        vn = v / (4 * 3.14 * rm * rm)
                        f.write ( "%.2f,%d\n" % (vmin, v) )
                    f.close()

        from OpenSave import SaveModeless
        SaveModeless ( title = 'Save Histogram',
                       filters = [('TXT', '*.txt', '.txt')],
                       initialfile = "dist_hist.txt", command = save )




    def Segment ( self ) :

        segMap = segmentation_map()
        if segMap == None :
            self.umsg ( "Please select a map in the Segment Map Dialog" )
            return

        smod = current_segmentation ()
        if smod == None :
            self.umsg ( "Please select a Current Segmentation in the Segment Map dialog" )
            return

        print "Seg has %d regions" % (len(smod.regions))


        print "Seg rads:", self.segRads.get()


        if hasattr(self, 'distByReg') :
            print "Found distByReg"
        else :
            self.umsg ( "Make Histogram first." )
            return



        sepRs = []
        for rstr in self.segRads.get().split(",") :
            try :
                radv = float(rstr)
            except :
                self.umsg ( "Error parsing distances; enter only numbers and commas" )
                return

            sepRs.append ( radv )
        sepRs.append ( 1e99 )


        self.umsg ( "Segmenting..." )


        print "Sep rads:", sepRs
        sregs = []
        for r in sepRs :
            sregs.append ( [] )

        for reg, rad in self.distByReg.iteritems() :
            #if reg.surface_piece != None :
            #    if reg.surface_piece.display == False :
            #        continue

            minRad = 0.0
            for i, maxRad in enumerate ( sepRs ) :
                if rad > minRad and rad <= maxRad :
                    sregs[i].append ( reg )
                    break

        for i, regs in enumerate (sregs) :
            print "%d - %d regs" % (i, len(regs))
            if len(regs) > 1 :
                try :
                    smod.join_regions ( regs )
                except :
                    self.umsg ( "An error occurred - regions may have changed - please start again." )
                    smod.display_regions()
                    return

        smod.display_regions()

        self.umsg ( "Done, created %d groups based on radial distances" % len(sregs)  )

        from segment_dialog import volume_segmentation_dialog
        volume_segmentation_dialog().ReportRegionCount ( smod )




    def GetMod ( self, name ) :

        for m in chimera.openModels.list() :
            if m.name == name :
                return m
        return None




    def LineCC ( self ) :

        dmap = segmentation_map()
        if dmap == None :
            umsg ( "No map selected" )
            return


        from chimera import Molecule
        mlist = OML(modelTypes = [Molecule])
        if len(mlist) == 0 :
            umsg ( "No molecule found" )
            return

        mol = mlist[0]

        print "Doing line CC in " + dmap.name + " using mol " + mol.name

        print dmap.openState.xform
        print mol.openState.xform


        rccs = []
        rmap = None
        rmap_pos = None
        rpoints, rpoint_weights = None, None
        xf = None

        resolution = 10.0

        for ri, res in enumerate ( mol.residues ) :
            try :
                cat = res.atomsMap["CA"][0]
            except :
                continue

            if rmap == None :
                rmap = makeMap ( "#%d:%d@CA" % (mol.id, res.id.position)
                                 , resolution, 1, (.5, .5, .5, 1.0), "resmap" )
                rmap_pos = cat.coord().toVector()
                print " - sphere map pos ", rmap_pos
                #rpoints, rpoint_weights = fit_points (rmap)
                rpoints, rpoint_weights = fit_points_old (rmap)
                xf = rmap.openState.xform

                break


        for radi in range ( 0, 1300, 1 ) :

            #d = cat.coord() - rmap_pos
            d = chimera.Vector(0,0,radi) - rmap_pos
            #print chimera.Vector(0,0,radi)
            trx = chimera.Xform.translation ( d )
            #xf = dmap.openState.xform.inverse
            xf2 = xf.__copy__()
            xf2.multiply ( trx )

            rmap.openState.xform = xf2
            break

            if 1 :
                rmap_values = dmap.interpolated_values ( rpoints, xf2 )
                olap, corr = overlap_and_correlation ( rpoint_weights, rmap_values )

                if radi % 100 == 0 :
                    print " %d - overlap: %f, cross-correlation: %f" % (radi, olap, corr)

                rccs.append ( [radi,corr] )
            #print corr,

        #chimera.openModels.close ( rmap )


        fp = open ( "lineCC.txt", "w" )
        for rad, cc in rccs :
            fp.write ( "%d\t%f\n" % (rad, cc) )

        fp.close ()



def overlap_and_correlation ( v1, v2 ):
    import FitMap
    olap, cor = FitMap.overlap_and_correlation ( v1, v2 )[:2]
    return olap, cor





def fit_points_old ( fmap, threshold = None ) :

    f_m = fmap.data.full_matrix();
    size = list(f_m.shape);
    size.reverse()

    points = VolumeData.grid_indices(size, numpy.single)  # i,j,k indices
    _contour.affine_transform_vertices( points, fmap.data.ijk_to_xyz_transform )
    weights = numpy.ravel(f_m).astype(numpy.single)

    threshold = fmap.surface_levels[0]
    #threshold = .3 * max ( numpy.ravel(f_m).astype(numpy.single) )

    ge = numpy.greater_equal(weights, threshold)
    points = numpy.compress(ge, points, 0)
    weights = numpy.compress(ge, weights)
    nz = numpy.nonzero( weights )[0]

    if len(nz) < len (weights) :
        points = numpy.take( points, nz, axis=0 )
        weights = numpy.take(weights, nz, axis=0)

    #mass = numpy.sum(weights, dtype=numpy.single)
    #fmap.rotation_center = numpy.dot(weights,points) / mass

    if 1 : print "FitPoints from %s with threshold %.4f, %d nonzero" % (
        fmap.name, threshold, len(nz) )

    return points, weights



def makeMap ( sel_str, res, gridSpacing, clr, map_name ) :

    cmd = "molmap %s %.3f sigmaFactor 0.187 gridSpacing %.3f replace false" % (
        sel_str, res, gridSpacing )
    #print ">>>", cmd
    chimera.runCommand ( cmd )

    mv = None
    for mod in chimera.openModels.list() :
        ts = mod.name.split()
        if len(ts) > 1 and mod.name.find("map") >=0 and mod.name.find("res") >=0 :
            #print " - found", mod.name
            mv = mod
            mv.name = map_name
            if 0 :
                #print " - saving to:", map_name
                mv.write_file ( map_name, "mrc" )
                xf = mv.openState.xform
                #print " - closing:", map_name
                chimera.openModels.close ( mv )
                mv = VolumeViewer.open_volume_file ( map_name )[0]
                #print " - opened:", mv.name
                mv.openState.xform = xf
            break

    if mv == None :
        umsg ("Map not generated.")
        return

    mv.surface_levels[0] = 0.001

    ro = VolumeViewer.volume.Rendering_Options()
    mv.update_surface ( False, ro )
    for sp in mv.surfacePieces :
        v, t = sp.geometry
        if len(v) == 8 and len(t) == 12 : sp.display = False
        sp.color = ( clr[0], clr[1], clr[2], clr[3] )

    return mv




def show_dialog (closeOld = True):

    from chimera import dialogs

    d = dialogs.find ( RSeg_Dialog.name, create=False )
    if d :
        if closeOld :
            d.toplevel_widget.update_idletasks ()
            d.Close()
            d.toplevel_widget.update_idletasks ()
        else :
            return d

    dialogs.register ( RSeg_Dialog.name, RSeg_Dialog, replace = True)

    d = dialogs.find ( RSeg_Dialog.name, create=True )
    # Avoid transient dialog resizing when created and mapped for first time.
    d.toplevel_widget.update_idletasks ()
    d.enter()

    return d



# -----------------------------------------------------------------------------
#
