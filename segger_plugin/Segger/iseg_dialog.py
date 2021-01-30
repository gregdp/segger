
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
from Segger import dev_menus, timing, seggerVersion

OML = chimera.openModels.list

REG_OPACITY = 0.45


# http://geomalgorithms.com/a06-_intersect-2.html



from segment_dialog import current_segmentation, segmentation_map


class ISeg_Dialog ( chimera.baseDialog.ModelessDialog ):

    title = "iSeg - Icosahedral Segmentation (Segger v" + seggerVersion + ")"
    name = "segger_iseg"
    buttons = ( "Close" )
    help = 'http://ncmi.bcm.edu/ncmi/software/segger/docs'

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
            l = Tkinter.Label(ff, text = "  1. Tools -> Higher-Order Structure -> Icosahedron Surface.", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)




        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = "    - show & match icosahedron to current map (change Orientation if necesary)", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :

            l = Tkinter.Label(ff, text = "  2. Make icosahedral surface mesh", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :

            l = Tkinter.Label(ff, text = "        ", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)

            b = Tkinter.Button(ff, text="Make", command=self.Icos2)
            b.grid (column=1, row=0, sticky='w', padx=5, pady=1)

            b = Tkinter.Button(ff, text="Toggle Display - Mesh/Solid", command=self.ToggleDisp)
            b.grid (column=3, row=0, sticky='w', padx=5, pady=1)


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :

            l = Tkinter.Label(ff, text = "  3. Push outward", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = "        # iterations: ", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)

            self.numIt = Tkinter.StringVar(ff)
            self.numIt.set ( "100" )
            e = Tkinter.Entry(ff, width=7, textvariable=self.numIt)
            e.grid(column=1, row=0, sticky='w', padx=5, pady=1)

            l = Tkinter.Label(ff, text = ", stiffness: ", anchor = 'w')
            l.grid(column=2, row=0, sticky='ew', padx=5, pady=1)

            self.springF = Tkinter.StringVar(ff)
            self.springF.set ( "0.2" )
            e = Tkinter.Entry(ff, width=7, textvariable=self.springF)
            e.grid(column=3, row=0, sticky='w', padx=5, pady=1)


            b = Tkinter.Button(ff, text="Push", command=self.Icos2Push)
            b.grid (column=4, row=0, sticky='w', padx=5, pady=1)



        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = "  -  Set radius:", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)



        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :

            l = Tkinter.Label(ff, text = "        ", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)

            sv = Tkinter.StringVar(ff)
            sv.trace("w", lambda name, index, mode, sv=sv: self.set_rad_changed_cb(sv.get()) )
            self.setRad = sv

            e = Tkinter.Entry(ff, width=7, textvariable=sv )
            e.grid(column=1, row=0, sticky='w', padx=5, pady=1)


            # Radius
            #rs = Hybrid.Scale(ff, '', 1, 1500, 0.01, 1150, length=200)
            #rs.frame.grid(row = row, column = 1, sticky = 'ew', padx=5, pady=1, columnspan=10)
            #rs.entry.config ( width=100 )

            #rs.callback(self.radius_changed_cb)
            #rs.entry.bind('<KeyPress-Return>', self.radius_changed_cb)
            #self.radius = rs


            self.rad = Tkinter.DoubleVar(ff)
            self.rad.set ( 100 )

            smod = self.GetMod ( "Icosahedron Faces"  )
            if smod != None :
                print "Found faces..."
                verts, tris = smod.icosVerts0, smod.icosTris
                p1 = smod.icosVerts [ tris[0][0] ]
                r = numpy.sqrt ( numpy.sum(p1*p1) )
                p1 = smod.icosVerts0 [ tris[0][0] ]
                r0 = numpy.sqrt ( numpy.sum(p1*p1) )
                print " - rad %.4f, orig: %.4f" % (r, r0)
                self.rad.set ( r )


            self.radius = Tkinter.Scale(ff, from_=0, to=1500, variable=self.rad, orient=Tkinter.HORIZONTAL, length=350, command=self.radius_changed_cb)
            self.radius.grid(column=2, row=0, sticky='w', padx=5, pady=1, columnspan=10)



            row = row + 1


            #ff = Tkinter.Frame(f)
            #ff.grid(column=0, row=row, sticky='w')
            #w = Scale(from_=0, to=100, resolution=0.1)




        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :

            l = Tkinter.Label(ff, text = "  5. Cross-correlation / Mask densities between", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = "        start radius: ", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)

            self.startRad = Tkinter.StringVar(ff)
            e = Tkinter.Entry(ff, width=7, textvariable=self.startRad)
            e.grid(column=1, row=0, sticky='w', padx=5, pady=1)

            l = Tkinter.Label(ff, text = ", end radius: ", anchor = 'w')
            l.grid(column=2, row=0, sticky='ew', padx=5, pady=1)

            self.endRad = Tkinter.StringVar(ff)
            e = Tkinter.Entry(ff, width=7, textvariable=self.endRad)
            e.grid(column=3, row=0, sticky='w', padx=5, pady=1)


            b = Tkinter.Button(ff, text="CC", command=self.Icos2CC)
            b.grid (column=4, row=0, sticky='w', padx=5, pady=1)

            #b = Tkinter.Button(ff, text="+CC", command=self.Icos2PushCC)
            #b.grid (column=5, row=0, sticky='w', padx=5, pady=1)





        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = "  6. Radii separated by commas:", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = "   ", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)

            self.segRads = Tkinter.StringVar(ff)

            if 0 or dev_menus :
                self.segRads.set ( "" )

            e = Tkinter.Entry(ff, width=40, textvariable=self.segRads)
            e.grid(column=1, row=0, sticky='w', padx=5, pady=1)


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = "   ", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)

            b = Tkinter.Button(ff, text="Mask Map", command=self.Icos2Map0)
            b.grid (column=1, row=0, sticky='w', padx=5, pady=1)

            b = Tkinter.Button(ff, text="Group Regions", command=self.Segment2)
            b.grid (column=2, row=0, sticky='ew', padx=5, pady=1)


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




    def Icos2 ( self ) :

        imod = self.GetMod ("Icosahedron")

        axmods = []
        for m in chimera.openModels.list() :
            if m.name == "Icosahedron Faces" :
                axmods.append ( m )

        if len(axmods) > 0 :
            chimera.openModels.close ( axmods )

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


        self.umsg ( "Building Icos2" )


        import _surface
        surf_mod = _surface.SurfaceModel()
        surf_mod.name = "Icosahedron Faces"
        chimera.openModels.add([surf_mod], sameAs = imod)

        import axes; reload (axes)

        self.icos_vecs = []
        from numpy import arccos, pi


        for p in imod.surfacePieces :
            v, t = p.geometry[0], p.geometry[1]
            #print len(v), len(t)

            #for pt in v :
            #    print " - pt: ", pt

            surf_mod.icosVerts0 = numpy.copy ( v )
            surf_mod.icosVerts = numpy.copy ( v )
            surf_mod.icosTris = numpy.copy ( t )
            surf_mod.nvecs = numpy.zeros ( (len(t), 3) )
            surf_mod.sps = []


            for ti, tri in enumerate ( t ) :
                #print " - tri: ", tri,
                p1 = v [ tri[0] ]
                p2 = v [ tri[1] ]
                p3 = v [ tri[2] ]

                mp = (p1 + p2 + p3) / 3.0
                pv = chimera.Vector ( mp[0], mp[1], mp[2] )
                r = pv.length
                pv.normalize()
                #print mp
                #self.icos_vecs.append ( pv )
                mp = mp / r

                #cyl = axes.AddCylinderSolid ( chimera.Vector(0,0,0), pv, r, (.6,.4,.4,1), 10.0, surf_mod )
                #cyl.name = "Icosahedron_Axes"

                sp = axes.TriangleMeshDiv ( p1, p2, p3, 50.0, None, None, surf_mod )
                #sp = surf_mod.surfacePieces [ len(surf_mod.surfacePieces)-1 ]
                sp.N = numpy.array ( pv, numpy.float32 )
                #surf_mod.nvecs.append ( mp )
                surf_mod.nvecs[ti] = mp
                surf_mod.sps.append ( sp )
                sp.ind = ti


                #p1v = chimera.Vector ( p1[0], p1[1], p1[2] ); p1v.normalize ()
                #p2v = chimera.Vector ( p2[0], p2[1], p2[2] ); p2v.normalize ()
                #p3v = chimera.Vector ( p3[0], p3[1], p3[2] ); p3v.normalize ()

                #a1 = arccos ( p1v * pv ) * 180.0 / pi
                #a2 = arccos ( p2v * pv ) * 180.0 / pi
                #a3 = arccos ( p3v * pv ) * 180.0 / pi

                #a12 = arccos ( p1v * p2v ) * 180.0 / pi

                #print a1, a2, a3, a12

                #if ti >= 0 :
                #    break

        p1 = surf_mod.icosVerts0 [ surf_mod.icosTris[0][0] ]
        r0 = numpy.sqrt ( numpy.sum(p1*p1) )

        self.umsg ( "Made Icos2 from %d sps in %s -> %d sps, rad %.1f" % (len(imod.surfacePieces), imod.name, len(surf_mod.surfacePieces), r0 ) )

        self.rad.set ( r0 )




    def ToggleDisp ( self ) :

        smod = self.GetMod ( "Icosahedron Faces" )
        if smod == None :
            self.status ( "Did not find Icos2" )
            return

        import _surface
        nmod = _surface.SurfaceModel()
        nmod.name = smod.name

        nmod.icosVerts0 = smod.icosVerts0
        nmod.icosVerts = smod.icosVerts
        nmod.icosTris = smod.icosTris
        nmod.nvecs = smod.nvecs
        nmod.sps = []

        for spi, sp in enumerate ( smod.sps ) :

            v, t = sp.geometry
            #print " sp %d - %d verts, %d tris" % (spi, len(v), len(t) )

            if len(v) > 0 and len(t) > 0 :
                ns = nmod.addPiece ( v, t, sp.color )
                nmod.sps.append ( ns )
                ns.N = sp.N
                ns.ind = spi
                if hasattr ( sp, 'verts0' ) :
                    ns.verts0 = sp.verts0

                if sp.displayStyle == sp.Mesh :
                    ns.displayStyle = sp.Solid
                else :
                    ns.displayStyle = sp.Mesh

        chimera.openModels.close ( [smod] )
        #chimera.openModels.add([nmod], sameAs = smod)
        chimera.openModels.add ( [nmod] )
        smod = nmod


        self.status ( "Toggle Display %s - %d surfaces" % ( smod.name, len(smod.surfacePieces)  ) )



    def NearMaps ( self, sp ) :

        #print " - making near maps"

        verts, tris = sp.geometry

        nmaps = [ None ] * len(verts)
        for vi in range (len(verts)) :
            #nsets[vi] = sets.Set()
            nmaps[vi] = {}

        def setn ( vi1, vi2 ) :
            #nsets[ t[0] ].add ( t[1] )
            s = nmaps[vi1]
            if vi2 not in s :
                v = verts[vi1] - verts[vi2]
                s[vi2] = numpy.sqrt ( numpy.sum(v*v) )

        for t in tris :
            setn ( t[0], t[1] )
            setn ( t[0], t[2] )
            setn ( t[1], t[0] )
            setn ( t[1], t[2] )
            setn ( t[2], t[0] )
            setn ( t[2], t[1] )

        return nmaps




    def Icos2Push ( self ) :

        smod = self.GetMod ( "Icosahedron Faces"  )

        if smod == None :
            self.status ( "Did not find Icos2" )
            return

        N, f = 0, 0.0
        try :
            N = int ( self.numIt.get() )
        except :
            self.umsg ( "Invalid # iterations: " + self.numIt.get() )
            return

        try :
            f = float ( self.springF.get() )
        except :
            self.umsg ( "Invalid stiffness: " + self.springF.get() )
            return



        self.Icos2PushN ( smod, N, f )
        #self.Icos2PushNSym ( smod, 50 )
        self.fi = 2

        self.status ( "Pushing done - %d sps pushed" % len(smod.surfacePieces)  )



    # 700, .2 -- 875,921,964,1005,1025,1039,1150


    def Icos2PushN ( self, smod, N, springf ) :

        print " - pushing %s, %d surfaces - %d iter " % ( smod.name, len(smod.surfacePieces), N )

        for spi, sp in enumerate ( smod.surfacePieces ) :
            verts, tris = sp.geometry

            #print " - surface piece %d points %d tris, " % (len(verts), len(tris)), sp.N

            if not hasattr ( sp, 'nmaps' ) :
                sp.nmaps = self.NearMaps (sp)

            for iter in range ( N ) :                 # SGIV: 600
                for vi in range ( len(verts) ) :

                    nmap = sp.nmaps[vi]

                    f = 0.0
                    if len(nmap) >= 6 :
                        f = 1.0               # SGIV: 1

                    #vv = verts[vi]
                    #vvl = numpy.sqrt ( numpy.sum(vv*vv) )
                    #vv = vv / vvl
                    #fN = numpy.sum(vv*sp.N)

                    fv = 0.1 * sp.N

                    if 1 :
		                for vj, eqd in nmap.iteritems() :
		                    v = verts[vj] - verts[vi]
		                    vl = numpy.sqrt ( numpy.sum(v*v) )
		                    vn = v / vl
		                    ff = vl - eqd
		                    fv = fv + springf * ff * vn              # SGIV: 0.2

                    verts[vi] = verts[vi] + f * fv

                if iter % 10 == 0 :
                    self.status ( "Pushing %d/%d - iter %d/%d - springf %.1f" % (spi+1,len(smod.surfacePieces), iter, N, springf ) )

            sp.geometry = (verts,tris)
            sp.verts0 = numpy.copy ( verts )


    def Icos2PushNSym ( self, smod, N ) :

        print " - pushing - sym - %s, %d surfaces - %d iter " % ( smod.name, len(smod.surfacePieces), N )

        sp = smod.sps[0]
        verts,tris = sp.geometry

        if not hasattr ( sp, 'nmaps' ) :
            sp.nmaps = self.NearMaps (sp)

        for iter in range ( N ) :                 # SGIV: 600
            for vi in range ( len(verts) ) :

                nmap = sp.nmaps[vi]

                f = 0.0
                if len(nmap) >= 6 :
                    f = 1.0               # SGIV: 1

                #vv = verts[vi]
                #vvl = numpy.sqrt ( numpy.sum(vv*vv) )
                #vv = vv / vvl
                #fN = numpy.sum(vv*sp.N)

                fv = 0.1 * sp.N

                if 1 :
	                for vj, eqd in nmap.iteritems() :
	                    v = verts[vj] - verts[vi]
	                    vl = numpy.sqrt ( numpy.sum(v*v) )
	                    vn = v / vl
	                    ff = vl - eqd
	                    fv = fv + 0.2 * ff * vn              # SGIV: 0.2

                verts[vi] = verts[vi] + f * fv

            if iter % 10 == 0 :
                self.status ( "Pushing - iter %d/%d" % ( iter, N ) )

        sp.geometry = (verts,tris)
        sp.verts0 = verts


        verts0, tris0 = smod.icosVerts0, smod.icosTris
        p1 = verts0 [ tris0[0][0] ]
        p2 = verts0 [ tris0[0][1] ]
        p3 = verts0 [ tris0[0][2] ]
        #mp = (p1 + p2 + p3) / 3.0
        a0 = numpy.array ( [p1,p2,p3] )
        #print a0

        import chimera.match

        for ti, tri in enumerate ( smod.icosTris[1:] ) :

            q1 = verts [ tri[0] ]
            q2 = verts [ tri[1] ]
            q3 = verts [ tri[2] ]
            a1 = numpy.array ( [q1,q2,q3] )
            #print a2

            xf = chimera.match.matchPositions ( numpy.array(a1,numpy.float), numpy.array(a0,numpy.float) )

            sp1 = smod.sps[ti]
            verts1, tris1 = sp1.geometry

            newv = numpy.zeros ( (len(verts),3) )

            for vi, v in enumerate ( verts ) :
                tp = xf[0].apply ( chimera.Point( v[0], v[1], v[2] ) )
                #print v, "->", tp
                newv[vi] = numpy.array ( tp )

            sp1.geometry = (newv,tris1)
            sp1.verts0 = newv







    def Icos2PushCC ( self ) :

        smod = self.GetMod ( "Icosahedron Faces"  )

        if smod == None :
            self.status ( "Did not find Icos2" )
            return

        print "Push/CC..."

        self.Icos2PushN ( smod, 100 )

        for i in range ( 20 ) :
            self.Icos2PushN ( smod, 100 )
            self.fi = 200 + i*100
            self.Icos2CC ()
            self.updateIcos2 ( 1110 )


        delattr ( self, 'fi' )



    def Icos2CC ( self ) :

        smod = self.GetMod ( "Icosahedron Faces" )

        if smod == None :
            self.umsg ( "No Icos2 found" )
            return

        dmap = segmentation_map()
        if dmap == None :
            self.umsg ( "No map selected" )
            return


        start_rad, end_rad = 0, 0
        try :
            start_rad = int ( self.startRad.get() )
        except :
            self.umsg ( "Invalid start radius: " + self.startRad.get() )
            return

        try :
            end_rad = int ( self.endRad.get() )
        except :
            self.umsg ( "Invalid end radius: " + self.endRad.get() )
            return


        if end_rad <= start_rad :
            self.umsg ( "End rad should be larger than start rad :) " )
            return


        self.umsg ( "CC in %s" % dmap.name )

        fname = "IcosCC.txt"
        if hasattr ( self, 'fi' ) :
            fname = "IcosCC_%d.txt" % self.fi


        p1 = smod.icosVerts [ smod.icosTris[0][0] ]
        rS = numpy.sqrt ( numpy.sum(p1*p1) )
        print " - rad before: ", rS

        ccs = []
        #fp = open ( fname, "w" )
        for rad in range ( start_rad, end_rad+1 ) :
            self.updateIcos2 ( rad )
            cc = self.IcosCC ( smod, dmap )
            self.status ( "Rad: %d, CC: %.4f" % (rad, cc) )
            #fp.write ( "%d\t%f\n" % (rad, cc) )
            ccs.append ( [rad, cc] )

        #fp.close ()

        self.updateIcos2 ( rS )


        def save ( okay, dialog ):
            if okay:
                paths = dialog.getPaths ( )
                if paths:
                    path = paths[0]
                    self.umsg ( "Saved CCs to: " + path )
                    f = open ( path, "w" )
                    for rad,cc in ccs :
                        f.write ( "%d\t%f\n" % (rad, cc) )
                    f.close()

        from OpenSave import SaveModeless
        SaveModeless ( title = 'Save Cross Correlations',
                       filters = [('TXT', '*.txt', '.txt')],
                       initialfile = "rad_cc.txt", command = save )







    def IcosCC ( self, smod, dmap ) :
        #newv = numpy.zeros_like ( verts )
        numv = len(smod.surfacePieces[0].geometry[0]) * len(smod.surfacePieces)
        #print "%d verts, %d sps, %d points" % ( len(smod.surfacePieces[0].geometry[0]), len(smod.surfacePieces), numv )
        newv = numpy.zeros ( (numv,3) )

        for spi, sp in enumerate ( smod.sps ) :
            verts, tris = sp.geometry
            v0 = spi * len(smod.surfacePieces[0].geometry[0])
            v1 = v0 + len(smod.surfacePieces[0].geometry[0])
            newv[v0:v1] = verts

        #print newv
        map_values = dmap.interpolated_values ( newv, dmap.openState.xform )
        #print map_values

        olap, cor = FitMap.overlap_and_correlation ( numpy.ones_like(map_values), map_values )[:2]
        #print olap, cor
        return cor


    def set_rad_changed_cb ( self, newRad ) :

        #print newRad
        try :
            nrad = int ( newRad )
            self.radius.set ( nrad )
        except :
            pass


    def radius_changed_cb(self, newRad) :

        #radius = self.radius.value(1000)
        #print "Radius: ", newRad
        #self.setRad.set ( newRad )


        radius = int ( newRad )

        self.updateIcos2 ( radius )



    def updateIcos2 ( self, rad ) :

        smod = self.GetMod ( "Icosahedron Faces" )
        if smod == None :
            #self.umsg ( "No Icosahedron2 model found" )
            return

        verts, tris = smod.icosVerts0, smod.icosTris
        p1 = verts [ tris[0][0] ]
        p2 = verts [ tris[0][1] ]
        p3 = verts [ tris[0][2] ]
        mp = (p1 + p2 + p3) / 3.0
        rad0 = numpy.sqrt ( numpy.sum(p1*p1) )
        rad1 = numpy.sqrt ( numpy.sum(mp*mp) )

        fscale = rad / rad0
        sphf = 1.0 - min ( rad, rad0 ) / rad0
        #self.status ( "Rad: %.3f -- rad: %.3f, midRad: %.3f, f: %.3f" % (rad, rad0, rad1, sphf) )

        for spi, sp in enumerate ( smod.surfacePieces ) :
            #sp0 = imod.surfacePieces[spi]
            verts, tris = sp.geometry

            if not hasattr ( sp, 'verts0' ) :
                sp.verts0 = verts
                #print "set init verts"

            #print " - surface piece %d points %d tris, " % (len(verts), len(tris)), sp.N
            newv = numpy.zeros_like ( verts )

            for vi, v in enumerate ( verts ) :
                iv = fscale * sp.verts0[vi]
                newv[vi] = iv
                #vv = v / numpy.sqrt ( numpy.sum (v*v) )
                #sv = vv * min ( rad, rad0 )
                #newv[vi] = sphf * sv + (1.0-sphf) * iv

            sp.geometry = (newv,tris)


        for vi, v in enumerate ( smod.icosVerts0 ) :
            smod.icosVerts[vi] = fscale * smod.icosVerts0[vi]

        #p1 = smod.icosVerts [ tris[0][0] ]
        #r = numpy.sqrt ( numpy.sum(p1*p1) )
        #p1 = smod.icosVerts0 [ tris[0][0] ]
        #r0 = numpy.sqrt ( numpy.sum(p1*p1) )
        #print "Icos - rad %.4f, orig: %.4f" % (r, r0)




    def GetMod ( self, name ) :

        for m in chimera.openModels.list() :
            if m.name == name :
                return m
        return None


    def MakeTNorms ( self, smod ) :

        self.umsg ( "Making triangle norms for %d" % len(smod.sps) )

        for spi, sp in enumerate ( smod.sps ) :

            verts2, tris2 = sp.geometry

            #sp.tdirs = [None] * len(tris2)
            sp.tdirs = numpy.zeros ( ( len(tris2), 3 ) )
            sp.tnorms = [None] * len(tris2)


            for ti, tri in enumerate ( tris2 ) :
                p1 = verts2 [ tri[0] ]
                p2 = verts2 [ tri[1] ]
                p3 = verts2 [ tri[2] ]
                mp = (p1 + p2 + p3) / 3.0
                l = numpy.sqrt ( numpy.sum(mp*mp) )
                sp.tdirs[ti] = mp / l

                v1 = p2 - p1
                v2 = p3 - p1
                N = numpy.cross ( v1, v2 )
                l = numpy.sqrt ( numpy.sum(N*N) )
                sp.tnorms [ti] = N / l



    def MinRad2 ( self, smod ) :
        minr = 1e9
        for sp in smod.surfacePieces :
            verts2, tris2 = sp.geometry
            for v in verts2 :
                r = numpy.sum ( v * v )
                if r < minr :
                    minr = r
        #return numpy.sqrt ( minr )
        return minr

    def MaxRad2 ( self, smod ) :
        maxr = 0
        for sp in smod.surfacePieces :
            verts2, tris2 = sp.geometry
            for v in verts2 :
                r = numpy.sum ( v * v )
                if r > maxr :
                    maxr = r
        #return numpy.sqrt ( maxr )
        return maxr


    def PIsOutside ( self, p, smod ) :

        #print "pt - %d surfps" % len(surfm.surfacePieces)

        #min_i = 0
        #max_d = -1e7
        #max_n = None
        #for nvi, nv in enumerate ( smod.nvecs ) :
        #    d = numpy.dot ( p, nv )
        #    if d > max_d :
        #        min_i = nvi
        #        max_d = d
        #        max_n = nv

        max_i = numpy.argmax ( numpy.sum ( smod.nvecs * p, axis = 1 ) )
        max_n = smod.nvecs [ max_i ]

        tri = smod.icosTris [ max_i ]

        p1 = smod.icosVerts [ tri[0] ]
        #p2 = smod.icosVerts [ tri[1] ]
        #p3 = smod.icosVerts [ tri[2] ]

        #v1 = p2 - p1
        #v2 = p3 - p1
        #N = numpy.cross ( v1, v2 )

        pv = p - p1
        d = numpy.dot ( pv, max_n )
        if d <= 0.0 :
            #print " - inside the tri ", min_i
            return False

        #return True


        sp = smod.sps[max_i]

        #if sp.ind != min_i and not hasattr (sp, 'flagged') :
        #    print sp.ind, "?"
        #    sp.flagged = True


        verts2, tris2 = sp.geometry

        #if not hasattr ( sp, 'tdirs' ) :
        #sp.tdirs = [None] * len(tris2)
        #sp.tnorms = [None] * len(tris2)

        #min_i = 0
        #max_d = -1e7

        #for ti, tri in enumerate ( tris2 ) :
        #    d = numpy.dot ( p, sp.tdirs[ti] )
        #    if d > max_d :
        #        max_d = d
        #        min_i = ti


        max_i = numpy.argmax ( numpy.sum ( sp.tdirs * p, axis = 1 ) )

        tri = tris2[max_i]

        p1 = verts2 [ tri[0] ]
        pv = p - p1
        d = numpy.dot ( pv, sp.tnorms [max_i] )
        if d <= 0.0 :
            #print " - inside the tri ", min_i
            return False


        return True




    def Icos2Map0 ( self ) :

        smod = self.GetMod ( "Icosahedron Faces" )
        if smod == None :
            self.umsg ( "No Icosahedron2 model found" )
            return


        dmap = segmentation_map()
        if dmap == None :
            self.umsg ( "Select a map in Segment Map dialog" )
            return


        sepRs = self.segRads.get().split(",")
        print "Sep rads:", sepRs

        if len(sepRs) != 2 :
            self.umsg ( "Enter two radii separated by a comma" )
            return



        try :
            start_rad = int ( sepRs[0] )
        except :
            self.umsg ( "Invalid start radius: " + sepRs[0] )
            return

        try :
            end_rad = int ( sepRs[1] )
        except :
            self.umsg ( "Invalid end radius: " + sepRs[1] )
            return


        if end_rad <= start_rad :
            self.umsg ( "End rad should be larger than start rad :) " )
            return


        self.umsg ( "Mask %s, %d -> %d" % (dmap.name,start_rad,end_rad) )

        self.MakeTNorms ( smod )


        import time
        start = time.time()

        mm = dmap.full_matrix ()
        #m1 = numpy.zeros_like ( mm )

        # transform to index reference frame of ref_map
        f1 = dmap.data.ijk_to_xyz_transform

        from _contour import affine_transform_vertices as transform_vertices
        #f2 = xform_matrix ( mask_map.openState.xform )
        #f3 = xform_matrix ( ref_map.openState.xform.inverse() )
        #f4 = ref_map.data.xyz_to_ijk_transform
        #tf = multiply_matrices( f2, f1 )
        #tf = multiply_matrices( f3, tf )
        #tf = multiply_matrices( f4, tf )

        nm = numpy.zeros_like ( mm )

        self.updateIcos2 ( start_rad )
        minr, maxr = self.MinRad2 ( smod ), self.MaxRad2 ( smod )
        print " - start rad %d -- min rad %.1f, max rad %.1f" % ( start_rad, numpy.sqrt(minr), numpy.sqrt(maxr))

        done = time.time()
        elapsed = done - start
        print "Took: ", elapsed


        pt = numpy.array ( [[0,0,0]], numpy.float32 )
        p = pt[0]

        for i in range ( dmap.data.size[0] ) :
            self.status ( "Masking %s, outside radius %d, %d/%d" % (dmap.name, start_rad, i+1, dmap.data.size[0]) )
            p[0] = i * f1[0][0] + f1[0][3]
            for j in range ( dmap.data.size[1] ) :
                p[1] = j * f1[1][1] + f1[1][3]
                for k in range ( dmap.data.size[2] ) :
                	#p[2] = k * f1[2][2] + f1[2][3]
                    #pt = numpy.array ( [[i,j,k]], numpy.float32 )
                    #p[0],p[1],p[2] = ti,tj,tk
                    #transform_vertices ( pt, f1 )
					p[2] = k * f1[2][2] + f1[2][3]
					ptr = numpy.sum ( p*p )
					if ptr < minr :
					    pass
					elif ptr > maxr :
					    nm[k,j,i] = mm[k,j,i]
					elif self.PIsOutside ( pt[0], smod ) :
					    nm[k,j,i] = mm[k,j,i]


        self.updateIcos2 ( end_rad )
        minr, maxr = self.MinRad2 ( smod ), self.MaxRad2 ( smod )
        print " - end rad %d -- min rad %.1f, max rad %.1f" % (start_rad, numpy.sqrt(minr), numpy.sqrt(maxr))

        for i in range ( dmap.data.size[0] ) :
            self.status ( "Masking %s, inside radius %d, %d/%d" % (dmap.name, end_rad, i+1, dmap.data.size[0]) )
            p[0] = i * f1[0][0] + f1[0][3]
            for j in range ( dmap.data.size[1] ) :
                p[1] = j * f1[1][1] + f1[1][3]
                for k in range ( dmap.data.size[2] ) :
                    #pt = numpy.array ( [[i,j,k]], numpy.float32 )
                    #p[0],p[1],p[2] = ti,tj,tk
                    #transform_vertices ( pt, f1 )
					p[2] = k * f1[2][2] + f1[2][3]
					ptr = numpy.sum ( p*p )
					if ptr < minr :
					    continue
					elif ptr > maxr :
					    nm[k,j,i] = 0.0
					elif self.PIsOutside ( p, smod ) :
					    nm[k,j,i] = 0.0




        ndata = VolumeData.Array_Grid_Data ( nm, dmap.data.origin, dmap.data.step, dmap.data.cell_angles )
        try : nvg = VolumeViewer.volume.add_data_set ( ndata, None )
        except : nvg = VolumeViewer.volume.volume_from_grid_data ( ndata )
        nvg.name = dmap.name + "__%d--to--%d_fast" % (start_rad, end_rad)

        done = time.time()
        elapsed = done - start
        print "Took: ", elapsed


    def Icos2Map0 ( self ) :

        dmap = segmentation_map()
        if dmap == None :
            self.umsg ( "Select a map in Segment Map dialog" )
            return


        mm = dmap.full_matrix ()
        #m1 = numpy.zeros_like ( mm )

        # transform to index reference frame of ref_map
        f1 = dmap.data.ijk_to_xyz_transform

        nm = numpy.zeros_like ( mm )

        minr, maxr = 300, 400
        pt = numpy.array ( [[0,0,0]], numpy.float32 )
        p = pt[0]

        im, jm, km = dmap.data.size[0]/2, dmap.data.size[1]/2, dmap.data.size[2]/2

        for i in range ( dmap.data.size[0] ) :
            self.status ( "Masking %s %.1f->%.1f, %d/%d" % (dmap.name, minr, maxr, i+1, dmap.data.size[0]) )
            di = abs(i-im) * dmap.data.step[0]

            for j in range ( dmap.data.size[1] ) :
                dj = abs(j-jm) * dmap.data.step[1]

                for k in range ( dmap.data.size[2] ) :
                    dk = abs(k-km) * dmap.data.step[2]
                    r = numpy.sqrt ( di*di + dj*dj + dk*dk )
                    if dk >= minr and dk < maxr :
                        nm[k,j,i] = mm[k,j,i]


        ndata = VolumeData.Array_Grid_Data ( nm, dmap.data.origin, dmap.data.step, dmap.data.cell_angles )
        try : nvg = VolumeViewer.volume.add_data_set ( ndata, None )
        except : nvg = VolumeViewer.volume.volume_from_grid_data ( ndata )
        nvg.name = dmap.name + "__%.0f--to--%.0f" % (minr, maxr)





    def Segment2 ( self ) :

        dmap = segmentation_map()
        if dmap == None :
            self.umsg ( "Please select a map in the Segment Map Dialog" )
            return

        smod = current_segmentation ()
        if smod == None :
            self.umsg ( "Please select a Current Segmentation in the Segment Map dialog" )
            return

        print "Seg has %d regions" % (len(smod.regions))


        imod2 = self.GetMod ( "Icosahedron Faces" )
        if imod2 == None :
            self.umsg ( "No Icosahedron2 model found" )
            return


        sepRs = []
        for rstr in self.segRads.get().split(",") :
            try :
                radv = float(rstr)
            except :
                self.umsg ( "Error parsing distances; enter only numbers and commas" )
                return

            sepRs.append ( radv )

        print "Sep rads:", sepRs
        regs = list(smod.regions)
        sregs = []

        f1 = dmap.data.ijk_to_xyz_transform
        from _contour import affine_transform_vertices as transform_vertices


        self.MakeTNorms ( imod2 )


        for i, srad in enumerate ( sepRs ) :

            self.umsg ( "Segmenting using %s - rad %.1f - %d regs" % ( imod2.name, srad, len(regs) ) )
            self.updateIcos2 ( srad )

            gregs, left_regs = [], []

            for ri, r in enumerate ( regs ) :

                p = r.max_point
                #pt = numpy.array ( [ [ p[2],p[1],p[0] ] ], numpy.float32 )
                pt = numpy.array ( [ [ p[2],p[1],p[0] ] ], numpy.float32 )
                transform_vertices ( pt, f1 )

                c = r.center_of_points()
                ptc = numpy.array ( [ c ], numpy.float32 )

                #print ri, p, c, pt[0]
                #return

                if self.PIsOutside ( ptc[0], imod2 ) :
                    #print " - outside"
                    left_regs.append ( r )
                else :
                    #print " - inside"
                    gregs.append ( r )

                if ri % 1000 == 0 :
                    self.status ( "Segmenting using %s - rad %.1f - %s/%s regs" % ( imod2.name, srad, "{:,}".format(ri), "{:,}".format(len(regs)) ) )


            sregs.append ( gregs )
            regs = left_regs
            print " - rad %.1f - %d regions inside" % ( srad, len(gregs) )

        print " - remaining %d regions" % ( len(regs) )
        sregs.append ( regs )


        for i, regs in enumerate (sregs) :
            self.status ( "Segmenting, layer %d - %d regs" % (i, len(regs)) )
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

    d = dialogs.find ( ISeg_Dialog.name, create=False )
    if d :
        if closeOld :
            d.toplevel_widget.update_idletasks ()
            d.Close()
            d.toplevel_widget.update_idletasks ()
        else :
            return d

    dialogs.register ( ISeg_Dialog.name, ISeg_Dialog, replace = True)

    d = dialogs.find ( ISeg_Dialog.name, create=True )
    # Avoid transient dialog resizing when created and mapped for first time.
    d.toplevel_widget.update_idletasks ()
    d.enter()

    return d



# -----------------------------------------------------------------------------
#
