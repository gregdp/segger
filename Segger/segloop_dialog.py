
# Copyright (c) 2009 Greg Pintilie - pintilie@mit.edu

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
import Surface
import VolumeViewer
from sys import stderr
from time import clock

from axes import prAxes
import regions
import graph
from Segger import dev_menus, timing, seggerVersion
from CGLutil.AdaptiveTree import AdaptiveTree
import random
from VolumePath import Marker_Set, Marker, Link
from _contour import affine_transform_vertices as transform_vertices
from Matrix import xform_matrix, multiply_matrices, chimera_xform, identity_matrix, invert_matrix, shift_and_angle



OML = chimera.openModels.list

REG_OPACITY = 0.45


from segment_dialog import current_segmentation, segmentation_map



aaCodes = {}
aaCodes["G"] = "GLY"
aaCodes["P"] = "PRO"
aaCodes["A"] = "ALA"
aaCodes["V"] = "VAL"
aaCodes["L"] = "LEU"
aaCodes["I"] = "ILE"
aaCodes["M"] = "MET"
aaCodes["C"] = "CYS"
aaCodes["F"] = "PHE"
aaCodes["Y"] = "TYR"
aaCodes["W"] = "TRP"
aaCodes["H"] = "HIS"
aaCodes["K"] = "LYS"
aaCodes["R"] = "ARG"
aaCodes["Q"] = "GLN"
aaCodes["N"] = "ASN"
aaCodes["E"] = "GLU"
aaCodes["D"] = "ASP"
aaCodes["S"] = "SER"
aaCodes["T"] = "THR"


def umsg ( txt ) :
    print txt
    status ( txt )

def status ( txt ) :
    txt = txt.rstrip('\n')
    msg.configure(text = txt)
    msg.update_idletasks()


class Segloop_Dialog ( chimera.baseDialog.ModelessDialog ):

    title = "SegLoop (Segger v" + seggerVersion + ")"
    name = "segger_segloop"
    buttons = ( "Close" )
    help = 'https://github.com/gregdp/segger/wiki'



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
        
        #l = Tkinter.Label(f, text='  ')
        #l.grid(column=0, row=row, sticky='w')

        

        # ---------------------------------------------------------------------------------



        if 0 :
            ff = Tkinter.Frame(f)
            ff.grid(column=0, row=row, sticky='w', pady=5, padx=10)
    
            l = Tkinter.Label(ff, text='  Map:', width=15, anchor=Tkinter.E)
            l.grid(column=0, row=0, sticky='w')
    
            self.dmap = Tkinter.StringVar(parent)
    
            self.mb  = Tkinter.Menubutton ( ff, textvariable=self.dmap, relief=Tkinter.RAISED )
            self.mb.grid (column=1, row=0, sticky='we', padx=5)
            self.mb.menu  =  Tkinter.Menu ( self.mb, tearoff=0, postcommand=self.MapMenu )
            self.mb["menu"]  =  self.mb.menu

            # set first visible map by default
            from VolumeViewer import Volume
            mlist = OML(modelTypes = [Volume])
            for m in mlist :
                if m.display :
                  self.dmap.set ( m.name )
                  self.cur_dmap = m
                  break



        if 1 :
            #row += 1
            ff = Tkinter.Frame(f)
            ff.grid(column=0, row=row, sticky='w', pady=5, padx=10)
    
            l = Tkinter.Label(ff, text='  Structure:', width=15, anchor=Tkinter.E)
            l.grid(column=0, row=0, sticky='w')
     
            self.struc = Tkinter.StringVar(parent)
            self.strucMB  = Tkinter.Menubutton ( ff, textvariable=self.struc, relief=Tkinter.RAISED )
            self.strucMB.grid (column=1, row=0, sticky='we', padx=5)
            self.strucMB.menu  =  Tkinter.Menu ( self.strucMB, tearoff=0, postcommand=self.StrucMenu )
            self.strucMB["menu"]  =  self.strucMB.menu

            from chimera import Molecule
            mlist = OML(modelTypes = [Molecule])
            for m in mlist :
                if m.display :
                  self.struc.set ( m.name )
                  self.cur_mol = m
                  break


            l = Tkinter.Label(ff, text="   Chain:" )
            l.grid(column=2, row=0, sticky='w')

            self.molChain = Tkinter.StringVar(f)
            e = Tkinter.Entry(ff, width=5, textvariable=self.molChain)
            e.grid(column=3, row=0, sticky='w', padx=5)


        if 1 :
            row += 1
            ff = Tkinter.Frame(f)
            ff.grid(column=0, row=row, sticky='w', pady=5, padx=10)
    
            l = Tkinter.Label(ff, text="Output:", width=15, anchor=Tkinter.E)
            l.grid(column=0, row=0, sticky='w')

            vmol = None
            for m in chimera.openModels.list() :
                if type(m) != chimera.Molecule or m.display == False : continue
                vmol = m
                break

            self.scMolName = Tkinter.StringVar(f)
            if vmol != None and hasattr ( vmol, 'openedAs' ) :
                path, molname = os.path.split ( vmol.openedAs[0] )
                fname, fext = os.path.splitext ( molname )
                outName = fname + "_SC" + fext
                self.scMolName.set ( outName )
            else :
                self.scMolName.set ( "mol_sc.pdb" )

            e = Tkinter.Entry(ff, width=30, textvariable=self.scMolName)
            e.grid(column=1, row=0, sticky='w', padx=5)



        row += 1
        dummyFrame = Tkinter.Frame(f, relief='groove', borderwidth=1)
        Tkinter.Frame(dummyFrame).pack()
        dummyFrame.grid(row=row,column=0,columnspan=1, pady=5, sticky='we')






        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w', pady=5, padx=10)
        if 1 :
            l = Tkinter.Label(ff, text="Start res:")
            l.grid(column=0, row=0, sticky='w')
    
            self.startRes = Tkinter.StringVar(f)


            e = Tkinter.Entry(ff, width=6, textvariable=self.startRes)
            e.grid(column=1, row=0, sticky='w', padx=5)


            l = Tkinter.Label(ff, text="End res:")
            l.grid(column=2, row=0, sticky='w')
    
            self.endRes = Tkinter.StringVar(f)


            if 0 or dev_menus :
                self.startRes.set ( "464" ); 
                self.endRes.set ( "492" )
            #self.startRes.set ( "110" ); self.endRes.set ( "112" )


            e = Tkinter.Entry(ff, width=6, textvariable=self.endRes)
            e.grid(column=3, row=0, sticky='w', padx=5)
            


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w', pady=5, padx=10)
        if 1 :
            l = Tkinter.Label(ff, text="Sequence:")
            l.grid(column=0, row=0, sticky='w')
    
            self.sequence = Tkinter.StringVar(f)

            #gp1
            self.sequence.set ( "MADNENRLESILSRFDADWTASDEARREAKNDLFFSRVSQWDDWLSQYTTLQYRGQFDVVRPVVRKLVSEMRQNPIDVLYRPKDGARPDAADVLMGMYRTDMRHNTAKIAVNIAVREQIEAGVGAWRLVTDYEDQSPTSNNQVIRREPIHSACSHVIWDSNSKLMDKSDARHCTVIHSMSQNGWEDFAEKYDLDADDIPSFQNPNDWVFPWLTQDTIQIAEFYEVVEKKETAFIYQDPVTGEPVSYFKRDIKDVIDDLADSGFIKIAERQIKRRRVYKSIITCTAVLKDKQLIAGEHIPIVPVFGEWGFVEDKEVYEGVVRLTKDGQRLRNMIMSFNADIVARTPKKKPFFWPEQIAGFEHMYDGNDDYPYYLLNRTDENSGDLPTQPLAYYENPEVPQANAYMLEAATSAVKEVATLGVDTEAVNGGQVAFDTVNQLNMRADLETYVFQDNLATAMRRDGEIYQSIVNDIYDVPRNVTITLEDGSEKDVQLMAEVVDLATGEKQVLNDIRGRYECYTDVGPSFQSMKQQNRAEILELLGKTPQGTPEYQLLLLQYFTLLDGKGVEMMRDYANKQLIQMGVKKPETPEEQQWLVEAQQAKQGQQDPAMVQAQGVLLQGQAELAKAQNQTLSLQIDAAKVEAQNQLNAARIAEIFNNMDLSKQSEFREFLKTVASFQQDRSEDARANAELLLKGDEQTHKQRMDIANILQSQRQNQPSGSVAETPQ" )
            self.sequence.set ( "YQSIVNDIYDVPRNVTITLEDGSEKDVQLMAEVVDLATGEKQVLNDIRGRY" )
            self.sequence.set ( "YQSIVNDIYDVPRNVTITLEDGSEKDVQLM" )
                        
            #gp9
            #self.sequence.set (  "AAAAAAAAAA" )

            self.seqText = Tkinter.Text(ff, width=50, height=3)
            self.seqText.grid(column=1, row=0, sticky='w', padx=5)
            
            if 0 or dev_menus :
                self.seqText.insert ( Tkinter.END, self.sequence.get() )
            


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w', pady=5, padx=10)
        if 1 :

            #b = Tkinter.Button(ff, text="ListR", command=self.ListRes)
            #b.grid (column=0, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Place Points", command=self.MakeCA)
            b.grid (column=1, row=0, sticky='w', padx=5)
            
            b = Tkinter.Button(ff, text="Make Path", command=self.FindPath)
            b.grid (column=2, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Minimize", command=self.MinGraph)
            b.grid (column=3, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Add Side Chains", command=self.AddSC)
            b.grid (column=4, row=0, sticky='w', padx=5)


        if 1 :

            self.whichGrads = Tkinter.StringVar(f)
            self.whichGrads.set ( "bonds,angles,self,other,othersc,density" )

            if 0 or dev_menus :
                row += 1
                ff = Tkinter.Frame(f)
                ff.grid(column=0, row=row, sticky='w', pady=5, padx=10)
        
                b = Tkinter.Button(ff, text="Gradients", command=self.Gradients)
                b.grid (column=0, row=0, sticky='w', padx=5)
    
                e = Tkinter.Entry(ff, width=30, textvariable=self.whichGrads)
                e.grid(column=1, row=0, sticky='w', padx=5)






        if 0 :
            row += 1
            dummyFrame = Tkinter.Frame(parent, relief='groove', borderwidth=1)
            Tkinter.Frame(dummyFrame).pack()
            dummyFrame.grid(row=row,column=0,columnspan=7, pady=3, sticky='we')
    
            row += 1
            f = Tkinter.Frame(parent)
            f.grid(column=0, row=row, sticky='ew', pady=5, padx=10)
            row += 1
            
            l = Tkinter.Label(f, text='To cite Segger or learn more about it press the Help button', fg="blue")
            l.grid(column=0, row=0, sticky='w')

            
        dummyFrame = Tkinter.Frame(parent, relief='groove', borderwidth=1)
        Tkinter.Frame(dummyFrame).pack()
        dummyFrame.grid(row=row,column=0,columnspan=7, pady=3, sticky='we')
        row += 1
        

        global msg
        msg = Tkinter.Label(parent, width = 60, anchor = 'w', justify = 'left', fg="red", pady=5, padx=10)
        msg.grid(column=0, row=row, sticky='ew')
        self.msg = msg
        
        umsg ( 'Select one or more segmented regions then press "Place Points" to start' )



    def MapMenu ( self ) :

        self.mb.menu.delete ( 0, 'end' )        # Clear menu
        from VolumeViewer import Volume
        mlist = OML(modelTypes = [Volume])
        for m in mlist :
            self.mb.menu.add_radiobutton ( label=m.name, variable=self.dmap,
                                           command=lambda m=m: self.MapSelected(m) )

    def MapSelected ( self, dmap ) :

        self.cur_dmap = dmap
        if dmap:
            dmap.display = True


    def StrucSelected ( self, mol ) :

        self.cur_mol = mol    
        print "Selected " + mol.name
        if mol :
            mol.display = True



    def StrucMenu ( self ) :

        self.strucMB.menu.delete ( 0, 'end' )   # Clear menu

        from chimera import Molecule
        mlist = OML(modelTypes = [Molecule])
        for m in mlist :
            self.strucMB.menu.add_radiobutton ( label=m.name, variable=self.struc,
                                           command=lambda m=m: self.StrucSelected(m) )




    def MakeHelix ( self ) :
        
        segMap = segmentation_map()
        if segMap == None :
            umsg ( "Please select a map in the Segment Map Dialog" )
            return

        smod = current_segmentation ()
        if smod == None :
            umsg ( "Please select a segmentation file in the Segment Map dialog" )
            return

        sregs = smod.selected_regions()
        if len(sregs)==0 : 
            umsg ( "no selected regions found" ); 
            return

        print "Helix for %d regions" % len(sregs)

        import axes
        reload (axes)


        if 0 :
            tpoints = numpy.concatenate ( [r.map_points() for r in sregs], axis=0 )

        else :
            for reg in sregs :
                tpoints = reg.map_points()
                h = Helix ()
                h.COM, h.U, h.S, h.V = axes.prAxes ( tpoints )

                print "Region %d - %f %f %f" % (reg.rid, h.S[0], h.S[1], h.S[2])

                com = numpy.sum(tpoints, axis=0) / len(tpoints)
                comv = numpy.ones_like ( tpoints ) * com
                points = tpoints - comv

                ppoints = points * h.U
                h.Extents = numpy.asarray ( numpy.max ( numpy.abs ( ppoints ), 0 ) )[0]

                print "  - extents: %f %f %f" % (h.Extents[0], h.Extents[1], h.Extents[2])

                h.heightAdj = float ( self.helixLength.get() )
                h.widthF = float ( self.helixWidthF.get() )

                print h.Extents
                h.MakeMod ( segMap )
                



    def ListRes ( self ) :

        mol = None
        for m in chimera.openModels.list() :
            if type(m) != chimera.Molecule or m.display == False : continue
            if m.name == "path" : continue
            mol = m
            break
        
        print "-----------------------------------------------------------------------"
        print "Mol: ", mol.name
        print "-----------------------------------------------------------------------"


        lastRes = None
        lastCat = None
        caDists = []
        for ri, res in enumerate ( mol.residues ) :
            cat = None
            try :
                cat = res.atomsMap["CA"][0]
            except :
                continue
                
            #print "%d - %d %s" % (ri, res.id.position, res.type),

            if lastRes != None :
                if lastRes.id.position != res.id.position - 1 :
                    print "%d - %d %s" % (ri, res.id.position, res.type)
                else :
                    v = cat.coord() - lastCat.coord()
                    print v.length,
                    caDists.append ( v.length )
                    
            
            print ""
            
            lastRes = res
            lastCat = cat

        avg = numpy.average ( caDists )
        stdev = numpy.std ( caDists )

        print "-----------------------------------------------------------------------"
        print " - Average CA dist: %.2f, stdev: %.2f" % (avg, stdev)




    def MakeGraph ( self, points, movef, mol ) :

        m = GetMod ( "graph" )
        if m != None :
            chimera.openModels.close ( [m] )


        from VolumePath import Marker_Set, Marker, Link
        g = Marker_Set ( "graph" )

        print " - graph for %d points" % len(points)


        markers = []
        #cpoints = numpy.zeros ( (len(points), 3) )
        apoints = numpy.zeros ( (len(points), 3) )
        amove = numpy.zeros ( (len(points), 3) )

        for i in range ( len(points) ) :

            p = points[i]
            clr = (.7, .7, .7, 1)
            if movef[i] == 0 : clr = (.7, .2, .2, 1)
            m = Marker(g, i, p, clr, .6)
            markers.append ( m )
            m.index = i
            #cpoints.append ( chimera.Point(p[0], p[1], p[2]) )
            apoints[i] = p
            
            if i >= 2 :
                amove[i] = [1,1,1]



        searchTreeAll = AdaptiveTree ( points, markers, 4.0)

        links = []
        mcons = []

        for i in range( len(points) ) :
            m = markers[i]
            markersNear = searchTreeAll.searchTree ( points[i], 20.0 )
            mcons.append ( markersNear )
            if len(markersNear) > 0 :
                for mnear in markersNear :
                    links.append ( [m, mnear] )

        print " -- %d links" % len(links)

        #import _surface
        #self.surf_mod = _surface.SurfaceModel()
        #chimera.openModels.add([surf_mod], sameAs = nmol)
        #surf_mod.name = "cons"

        #import axes
        #reload (axes)

        for link in links :
            m1, m2 = link[0], link[1]
            if m1.index > m2.index :
                Link ( m1, m2, (.5,.5,.5,1), .3 )
                
                #v = at2.coord() - at1.coord()
                #l = v.length
                #v.normalize()
                #cyl = axes.AddCylinderSolid ( at1.coord().data(), v, l, (1,.2,.2,1), .2, surf_mod )

        g.show_model ( True )

        gmod = GetMod ( "graph" )
        gmod.markers = markers
        gmod.links = links
        gmod.mcons = mcons
        gmod.points = points
        gmod.apoints = apoints
        #gmod.cpoints = cpoints
        gmod.movef = movef
        gmod.amove = amove
        #gmod.opoints = opoints
        #gmod.rmap = rmap
        
        return gmod
        



    def getCAPoints ( self ) :

        segMap = segmentation_map()
        if segMap == None :
            umsg ( "Please select a map in the Segment Map Dialog" )
            return


        mol = None
        for m in chimera.openModels.list() :
            if type(m) != chimera.Molecule or m.display == False : continue
            if m.name == "path" : continue
            mol = m
            break
            
        mol = None
        if hasattr ( self, "cur_mol" ) :
            mol = self.cur_mol

        print " - getting CA points for mol: ", mol.name


        startRi = int ( self.startRes.get() )
        endRi = int ( self.endRes.get() )
        startRes = None
        for r in mol.residues :
            if r.id.chainId == self.molChain.get() and r.id.position == startRi-1 :
                startRes = r
                break
        endRes = None
        for r in mol.residues :
            if r.id.chainId == self.molChain.get() and r.id.position == endRi+1 :
                endRes = r
                break


        g = GetMod ( "graph" )
        g.withoutRes = [startRes, endRes]


        print " - getting other CA from ", mol.name, " - chain ", self.molChain.get()

        numAt = 0
        for r in mol.residues :
            if r.id.chainId == self.molChain.get() :
                if r in g.withoutRes :
                    #print "skipping", r.id.position
                    continue
                else :
                    numAt += 1


        print " - %d CA points" % numAt

        opoints = numpy.zeros ( (numAt, 3) )
    
        rmap = {}
        numAt = 0
        for r in mol.residues :
            if r.id.chainId == self.molChain.get() :
                if r in g.withoutRes :
                    print "skipping", r.id.position
                    rmap[r.id.position] = r
                    continue
                else :
                    try :
                        tp = segMap.openState.xform.inverse().apply ( r.atomsMap["CA"][0].xformCoord() )
                        opoints[numAt] = tp.data()
                    except :
                        print " - res has no CA: ", r.id.position, r.type
                        continue
                    rmap[r.id.position] = r
                    numAt += 1


        CAPointsTree = AdaptiveTree ( opoints.tolist(), opoints.tolist(), 4.0 )

        g.CAPointsTree = CAPointsTree
        g.rmap = rmap





    def getSCPoints ( self ) :

        mol = None
        for m in chimera.openModels.list() :
            if type(m) != chimera.Molecule or m.display == False : continue
            if m.name == "path" : continue
            mol = m
            break

        g = GetMod ( "graph" )

        print " - getting side chain atoms - ", mol.name, " - ", g.name    
        numAt = 0
        for r in mol.residues :
            for at in r.atoms :
                if r.id.chainId == self.molChain.get() :
                    if "CA" in at.name or "H" in at.name :
                        pass
                    else :
                        numAt += 1

        print " - %d atoms" % numAt    

        scpoints = numpy.zeros ( (numAt, 3) )
        numAt = 0
        for r in mol.residues :
            if r.id.chainId == self.molChain.get() :
                for at in r.atoms :
                    if "CA" in at.name or "H" in at.name :
                        pass
                    else :
                        scpoints[numAt] = at.coord().data()
                        numAt += 1


        g.SCPointsTree = AdaptiveTree ( scpoints.tolist(), scpoints.tolist(), 4.0)




    def MakePath ( self, path ) :

        g = GetMod ( "graph" )
        if g == None :
            umsg ( "Did not find points - use 'Place Points' first" )
            return
        
        p = GetMod ( "path" )
        if p != None :
            chimera.openModels.close ( [p] ) 


        P = Marker_Set ( "path" )

        # print " - path for %d points, path length %d" % ( len(g.points), len(path) )


        markers = []

        for i in range ( len(path) ) :

            p = g.apoints [ path[i] ]
            m = Marker(P, path[i], p, (.4, .4, .7, 1), .7)
            markers.append ( m )
            m.index = path[i]
            
            if i > 0 :
                m1 = markers[i-1]
                m2 = markers[i]
                Link ( m1, m2, (.4,.4,.7,1), .4 )

        g.display = False
        P.show_model ( True )



    def ShowGradients ( self, graph ) :

        gra = GetMod ( "gradients" )
        if gra != None :
            chimera.openModels.close ( [gra] )

        import _surface
        amod = _surface.SurfaceModel()
        
        from axes import AddArrow2
    
        #amod = AddAxes ( rad, Extents[0]*f, Extents[1]*f, Extents[2]*f, 1.0, amod )
        amod.name = "gradients"
        
        print " - max grad: ", graph.maxGrad
        
        mul = 5.0 / graph.maxGrad

        for i in range ( len(graph.grads) ) :
            g = graph.grads[i]
            p = graph.apoints[i]
            
            gl = numpy.sqrt ( numpy.sum (g*g) )
            if gl > 1e-3 :
                #sc = numpy.sum ( g )
                g = g / gl
                amod = AddArrow2 ( p, chimera.Vector(g[0], g[1], g[2]), (gl*mul), (.4,.4,.4,1), 0.2, amod )


        chimera.openModels.add([amod], sameAs = graph)




    def PathLength ( self, path ) :

        g = GetMod ( "graph" )
        if g == None :
            umsg ( "Did not find points - use 'Place Points' first" )
            return

        plen = 0.0

        for i in range ( 1, len(path) ) :

            p1 = g.apoints [ path[i-1] ]
            p2 = g.apoints [ path[i] ]
            v = p2 - p1
            vlen = numpy.sqrt ( numpy.sum ( v * v ) )
            plen += vlen

        # print " - path for %d points, length %d, dist %.2f" % ( len(g.points), len(path), plen )

        return plen            



    def MinStep ( self, g, vol ) :

        self.GetGrads ( g, vol, ["bonds", "angles", "self", "other", "othersc", "density"] )
        g.apoints = g.apoints + ( g.grads * g.amove )

        #print " - min step max g: %.3f" % maxg


    def GradMagnAtPosFromTree ( self, p1, tree, cutoff=10, F=1.0 ) :
        gm = 0
        opointsNear = tree.searchTree ( p1.tolist(), cutoff )
        if len(opointsNear) > 0 :
            for p2 in opointsNear :
                v = p2 - p1
                vlen = numpy.sum ( v * v )
                if vlen < 1e-3 :
                    print " - atoms too close!"
                    gm = gm + F * 1000
                else :
                    gm = gm + F / vlen
        return gm


    def GradAtPosFromTree ( self, p1, tree, cutoff=10, F=0.1 ) :
        opointsNear = tree.searchTree ( p1.tolist(), cutoff )
        g = numpy.array ( [0,0,0] )
        if len(opointsNear) > 0 :
            for p2 in opointsNear :
                v = p2 - p1
                vlen = numpy.sqrt ( numpy.sum ( v * v ) )
                if vlen < 1e-5 :
                    print " - atoms too close!"
                    g = g + [10,10,10]
                else :
                    v = v / vlen
                    v1 = v * (-F/vlen)
                    g = g + v1
        return g


    def GetGrads ( self, g, vol, l ) :

        g.grads = numpy.zeros ( (len(g.points), 3) )
        
        #print " - grads: ", l

        if "bonds" in l :        

            deq = 3.82

            for i in range ( 1, len(g.path) ) :
                i1, i2 = g.path[i-1], g.path[i]
                p1, p2 = g.apoints [i1], g.apoints [i2]
                v = p2-p1
                vlen = numpy.sqrt ( numpy.sum ( v * v ) )
                d = vlen - deq
                v = v / vlen
                v1 = v * (d * 0.5)
                v2 = v * (d * -0.5)
                g.grads [i1] = g.grads [i1] + v1
                g.grads [i2] = g.grads [i2] + v2


        if "angles" in l :

            for i in range ( len(g.path) ) :

                i1 = g.path[i]
                cons = g.cons [ i1 ]
                if len(cons) == 1 :
                    i2 = cons[0]
                    p1, p2 = g.apoints [ i1 ], g.apoints [ i2 ]
                    ri1, ri2 = g.resids [ i1 ], g.resids [ i2 ]
                    ri0 = ri1 + (ri1 - ri2)
                    #print "ang ap:%d-%d rid:%d-%d-%d" % (i1,i2, ri0, ri1, ri2),
                    p0 = g.rmap [ ri0 ].atomsMap [ "CA" ][0].coord()
                    v1, v2 = p1 - p0, p2 - p1
                    l1, l2 = numpy.sqrt ( numpy.sum (v1*v1) ), numpy.sqrt ( numpy.sum (v2*v2) )
                    v1, v2 = v1/l1, v2/l2
                    cx = numpy.cross ( v1, v2 )
                    #lcx = numpy.sqrt ( numpy.sum (cx*cx) )
                    vp0, vp2 = numpy.cross ( v1, cx ), numpy.cross ( v2, cx )
                    #print " %.3f|%.3f|%.3f" % (lcx, numpy.sqrt ( numpy.sum (vp0*vp0) ), numpy.sqrt ( numpy.sum (vp2*vp2) ) )
                    g.grads [i2] = g.grads [i2] + (0.2*vp2)
                else :
                    i0, i2 = cons[0], cons[1]
                    p0, p1, p2 = g.apoints [ i0 ], g.apoints [ i1 ], g.apoints [ i2 ]
                    v1, v2 = p1 - p0, p2 - p1
                    l1, l2 = numpy.sqrt(numpy.sum(v1*v1)), numpy.sqrt(numpy.sum(v2*v2))
                    v1, v2 = v1/l1, v2/l2
                    cx = numpy.cross ( v1, v2 )
                    #lcx = numpy.sqrt ( numpy.sum (cx*cx) )
                    vp0, vp2 = numpy.cross ( v1, cx ), numpy.cross ( v2, cx )
                    #print " %.3f|%.3f|%.3f" % (lcx, numpy.sqrt ( numpy.sum (vp0*vp0) ), numpy.sqrt ( numpy.sum (vp2*vp2) ) )
                    g.grads [i0] = g.grads [i0] + (0.05*vp0)
                    g.grads [i2] = g.grads [i2] + (0.05*vp2)


        if "self" in l :

            stree = AdaptiveTree ( g.apoints.tolist(), g.markers, 3.0)

            for i1 in range( len(g.points) ) :
                p1 = g.apoints[i1]
    
                markersNear = stree.searchTree ( p1.tolist(), 10.0 )
                if len(markersNear) > 0 :
                    for mnear in markersNear :
    
                        i2 = mnear.index
                        
                        if ( i2 <= i1 ) :
                            continue
                        
                        if i1 in g.cons[i2] :
                            continue
                        
                        if i2 in g.cons[i1] :
                            continue
                        
                        p2 = g.apoints[i2]
    
                        v = p2 - p1
                        vlen = numpy.sqrt ( numpy.sum ( v * v ) )
                        if ( vlen < 10.0 ) :
                            v = v / vlen
                            v1 = v * (-0.2/vlen)
                            v2 = v * (0.2/vlen)
                                
                            g.grads [i1] = g.grads [i1] + v1
                            g.grads [i2] = g.grads [i2] + v2


        if "other" in l :        
            for i1 in range( len(g.points) ) :
                p1 = g.apoints[i1]
                gvec = self.GradAtPosFromTree ( p1, g.CAPointsTree, cutoff=10, F=0.1 )
                g.grads [i1] = g.grads [i1] + gvec


        if "othersc" in l :        
            for i1 in range( len(g.points) ) :
                p1 = g.apoints[i1]
                gvec = self.GradAtPosFromTree ( p1, g.SCPointsTree, cutoff=10, F=0.1 )
                g.grads [i1] = g.grads [i1] + gvec


        if "density" in l :
            xf = vol.openState.xform
            cpoints = numpy.copy ( g.apoints )

            if 0 :
                for i in range ( len(g.points) ) :
                    p = g.apoints[i]
                    px = xf.apply ( chimera.Point ( p[0], p[1], p[2] ) )
                    cpoints[i] = px.data()

            #transform_vertices ( cpoints.astype (numpy.float32),  )

            data_array, xyz_to_ijk_transform = vol.matrix_and_transform ( chimera.Xform(), subregion = None, step = 1 )
            #xfa = numpy.array ( [[1,0,0,0], [0,1,0,0], [0,0,1,0]] )

            #f4 = ref_map.data.xyz_to_ijk_transform

            tf = multiply_matrices( xyz_to_ijk_transform, xform_matrix ( xf ) )
            
            gradients, outside = VolumeData.interpolate_volume_gradient ( cpoints, tf, data_array )
            gradients = gradients * 0.2
            #print gradients
            
            g.grads = g.grads + gradients

        
        g.maxGrad = numpy.sqrt ( numpy.max ( numpy.sum ( g.grads * g.grads, 1 ) ) )
        #print " - max g: %.3f" %  g.maxGrad




    def MakeCA ( self ) :

        segMap = segmentation_map()
        if segMap == None :
            umsg ( "Please select a map in the Segment Map Dialog" )
            return

        smod = current_segmentation ()
        if smod == None :
            umsg ( "Please select a segmentation file in the Segment Map dialog" )
            return

        sregs = smod.selected_regions()
        if len(sregs)==0 : umsg ( "no selected regions found" ); return


        mol = None
        for m in chimera.openModels.list() :
            if type(m) != chimera.Molecule or m.display == False : continue
            if m.name == "path" : continue
            mol = m
            break
            

        mol = None
        if hasattr ( self, "cur_mol" ) :
            mol = self.cur_mol
        
        if mol == None :
            umsg ( "Please select a Molecule model in the Structure: field" )
            return


        print "Mol: ", mol.name


        close = []
        p = GetMod ( "graph" )
        if p != None : 
            close.append ( p )
        p = GetMod ( "path" )
        if p != None : 
            close.append ( p )
        p = GetMod ( "gradients" )
        if p != None : 
            close.append ( p )
        if len(close) != 0 : 
            chimera.openModels.close ( close ) 




        print " - %d selected regions" % len(sregs)
        tpoints = numpy.concatenate ( [r.map_points() for r in sregs], axis=0 )
        
        print " - %d grid points" %  len(tpoints)
        #print tpoints
        
        startRi = int ( self.startRes.get() )
        endRi = int ( self.endRes.get() )
        
        numRes = endRi - startRi + 1

        dpts = int ( numpy.floor ( float (len(tpoints)) / float(numRes) ) )

        print " - building from %d to %d - %d residues - every %d pt" % (startRi, endRi, numRes, dpts)

        startRes = None
        for r in mol.residues :
            if r.id.chainId == self.molChain.get() :
                if r.id.position == startRi-1 :
                    startRes = r
                    break

        endRes = None
        for r in mol.residues :
            if r.id.chainId == self.molChain.get() :
                if r.id.position == endRi+1 :
                    endRes = r
                    break

        if startRes == None or endRes == None :
            umsg ( "Did not find start and/or end residues" )
            return

        print " - start res: %d %s" % (startRes.id.position, startRes.type)
        print " -   end res: %d %s" % (endRes.id.position, endRes.type)



        points = []
        movef = []

        p = segMap.openState.xform.inverse().apply ( startRes.atomsMap["CA"][0].xformCoord() )
        points.append ( p )
        movef.append ( 0 )

        p = segMap.openState.xform.inverse().apply ( endRes.atomsMap["CA"][0].xformCoord() )
        points.append ( p )
        movef.append ( 0 )

        for i in range (numRes) :
            p = tpoints [i * dpts]
            points.append ( p )
            movef.append ( 1 )


        print " - %d nodes" % len(points)
        g = self.MakeGraph ( points, movef, mol )
        g.withoutRes = [startRes, endRes]
        g.startRes = startRes
        g.endRes = endRes
        g.mol = mol

        umsg ( 'Placed %d points, now press Make Path' % len(points) )
        
    
    
    def FindPathR ( self, g, toIndex, path, visited ) :

        visited [ toIndex ] = 1
        path = numpy.append ( path, toIndex )

        if toIndex == 1 :
            if len(path) == len(visited) :
                g.path = path
                print " - reached 1 - complete"
                return True
            else :
                try :
                    self.pathTry += 1
                except :
                    self.pathTry = 1
                #print " - reached 1 - not complete - ", self.pathTry
                return False

        mAt = g.markers [ toIndex ]
        cmarkers = g.mcons [ toIndex ]

        if len (cmarkers) == 0 :
            print " - nowhere to go at marker %d" % toIndex
            return False

        found = False
        random.shuffle(cmarkers)
        for m in cmarkers :

            if visited[m.index] == 1 :
                continue
            
            #print path, " --> %d" % m.index
            #print " > %d" % m.index,
            if self.FindPathR ( g, m.index, numpy.copy(path), numpy.copy(visited) ) :
                found = True
                break
            else :
                if self.pathTry > 5000 :
                    print " - giving up?"
                    found = False
                    break

        return found



    def FindPath ( self ) :

        g = GetMod ( "graph" )
        if g == None :
            umsg ( "Did not find points - use 'Place Points' first" )
            return


        print "FindPath - using graph with %d points" % len(g.points)

        path = numpy.array ( [], dtype=numpy.int )
        visited = numpy.zeros ( len(g.points), dtype=numpy.int )
        
        if self.FindPathR ( g, 0, path, visited ) == True :
            
            print " - found path: ",
            umsg ( 'Found a path' )
            print g.path
            self.PathLength ( g.path )

            g.resids = numpy.zeros ( len(g.path), numpy.int )
            g.resids[0] = g.startRes.id.position
            
            m1 = g.markers [ 0 ]
            m1.extra_attributes = { 'res' : g.resids[0] }
            
            print " [%d - %d] " % (g.resids[0], g.resids[1])
            
            g.cons = {}
            for i in range ( 1, len(g.path) ) :

                i1, i2 = g.path[i-1], g.path[i]

                g.resids[i2] = g.resids[i1] + 1

                m2 = g.markers [ i2 ]
                m2.extra_attributes = { 'res' : g.resids[i2] }

                try :
                    g.cons[i1].append ( i2 )
                except :
                    g.cons[i1] = [i2]
                try :
                    g.cons[i2].append ( i1 )
                except :
                    g.cons[i2] = [i1]

            self.MakePath ( g.path )


        else :
            umsg ( 'Did not find a path - points may be too far away - try a smaller segment' )

            
            
    def MinGraph ( self ) :

        g = GetMod ( "graph" )
        if g == None :
            umsg ( "Did not find points - use 'Place Points' first" )
            return



        vol = None
        for m in chimera.openModels.list() :
            if type(m) != VolumeViewer.volume.Volume or m.display == False : continue
            vol = m
            break


        vol = segmentation_map()
        if vol == None :
            umsg ( "Please select a map in the Segment Map Dialog" )
            return

        

        print "Map: ", vol.name


        self.getCAPoints ()
        self.getSCPoints ()

        #print " - using graph with %d points" % len(g.points)

        for i in range ( 100 ) :
            status ( 'Minimizing - step %d/%d' % (i+1,100) )
            self.MinStep ( g, vol )

        print " - max g: %.3f, path length: %.3f" % ( g.maxGrad, self.PathLength ( g.path ) )
        umsg ( 'Done minimizing - you can minimize again or press Add Side Chains' )

        self.MakePath ( g.path )
        
        self.Gradients ()



    def Gradients (self) :

        g = GetMod ( "graph" )
        if g == None :
            umsg ( "Did not find points - use 'Place Points' first" )
            return


        vol = segmentation_map()
        if vol == None :
            umsg ( "Please select a map in the Segment Map Dialog" )
            return

        print "Map: ", vol.name


        which = self.whichGrads.get()
        
        self.getCAPoints ()
        self.getSCPoints ()

        self.GetGrads ( g, vol, which.split(",") )
        self.ShowGradients ( g )



    def TreeForAtomsInMol ( self, mol ) :
        print " - making tree for all atoms..."
        numAt = 0
        for r in mol.residues :
            for at in r.atoms :
                if "H" in at.name :
                    pass
                else :
                    numAt += 1

        points = numpy.zeros ( (numAt, 3) )
        numAt = 0
        for r in mol.residues :
            for at in r.atoms :
                if "H" in at.name :
                    pass
                else :
                    points[numAt] = at.coord().data()
                    numAt += 1

        return AdaptiveTree ( points.tolist(), points.tolist(), 4.0)


    def AddSC ( self ) :

        g = GetMod ( "graph" )
        if g == None :
            print " - did not find graph"
            return


        segMap = segmentation_map()
        if segMap == None :
            self.umsg ( "Please select a map in the Segment Map Dialog" )
            return

        mol = None
        for m in chimera.openModels.list() :
            if type(m) != chimera.Molecule or m.display == False : continue
            if m.name == "path" : continue
            mol = m
            break



        mol = None
        if hasattr ( self, "cur_mol" ) :
            mol = self.cur_mol
        
        if mol == None :
            umsg ( "Please select a Molecule model in the Structure: field" )
            return




        print "Add side chains - %s" % mol.name

        resMap = {}
        maxResPosI = 0
        for res in mol.residues :
            if res.id.chainId == self.molChain.get() :
                if res.id.position > maxResPosI :
                    maxResPosI = res.id.position
            try :
                resMap[ res.type ].append (res)
            except :
                resMap[ res.type ] = [ res ]

                
        seq = self.sequence.get()
        seq = self.seqText.get(0.0, Tkinter.END).strip()
        
        print " - sequence: [" + seq + "]"
        plen = len(g.path) - 2
        print " - seq len: ", len(seq), " path len: ", plen
        
        if len(seq) <> plen :
            umsg ( "Sequence length (%d) should match path length (%d)" %(len(seq),plen) )
            return
        
        
        
        scMol = None
        g.newResAtomsTree = None
        
        missingRes = []

        mapxfi = segMap.openState.xform.inverse()


        firstResI = g.resids[ g.path[0] ]
        print " - First res: ", firstResI
        for i in range (1, firstResI) :
            rPosition = i
            
            r = None
            try :
                r = g.rmap [ rPosition ]
            except :
                pass
                
            if r == None :
                missingRes.append ( i ) 
                umsg ("Skipping res %d" % i)
                continue
            
            status ("Copying res %d - %s" % (i, r.type) )
            #print "%d/%d - res:" % (i, r.id.position), r.type
            
            scMol = self.AddResToMol ( r, rPosition, mapxfi, scMol )

            if len ( scMol.residues ) > 1 :
                prevRes = scMol.residues [ -2 ]
                newRes = scMol.residues [ -1 ]
                nb = scMol.newBond ( prevRes.atomsMap["C"][0], newRes.atomsMap["N"][0] )
                nb.display = nb.Stick
                nb.radius = 0.01


        for i in range ( 0, len(g.path) ) :

            i1 = g.path[i]
            rPosition = g.resids[i1]
            #rCode = seq[rPosition-1]
            try :
                rCode = seq[i-1]
            except :
                rCode = 'A'
                print "code??"

            print "%d|%d - res %d - %s" % (i, i1, rPosition, rCode),
            
            try :
                r = g.rmap [ rPosition ]
            except :
                r = None

            if r :
                print " - res:", r.id.position, r.type
                status ("Copying res %d (point %d/%d) - %s" % (rPosition, i, i1, rCode) )
                scMol = self.AddResToMol ( r, rPosition, mapxfi, scMol )

            else :
                print " - no res, adding..."
                umsg ("Adding res %d (point %d/%d) - %s" % (rPosition, i, i1, rCode) )
                if self.AddResSC ( g, scMol, i1, rCode, rPosition, resMap ) == None :
                    break
                #break

            # add N-C connection to last res
            if len ( scMol.residues ) > 1 :

                prevRes = scMol.residues [ -2 ]
                #at1 = prevRes.atomsMap["C"][0]
                #r1 = at1.residue

                newRes = scMol.residues [ -1 ]
                #at2 = newRes.atomsMap["N"][0]
                #r2 = at2.residue

                #print " at1 %s-%s-%s-%d, at2 %s-%s-%s-%d" % (
                #            at1.element.name, at1.name, r1.type, r1.id.position,
                #            at2.element.name, at2.name, r2.type, r2.id.position)

                nb = scMol.newBond ( prevRes.atomsMap["C"][0], newRes.atomsMap["N"][0] )
                nb.display = nb.Stick
                nb.radius = 0.01

            g.newResPointsTree = self.TreeForAtomsInMol ( scMol )


        lastResI = g.resids[ g.path[len(g.path)-1] ]
        print " - continuing with ", (lastResI+1)
        print " - last residue in mol: ", maxResPosI
        
        for i in range (lastResI+1, maxResPosI+1) :
            rPosition = i
            r = None
            try :
                r = g.rmap [ rPosition ]
            except :
                pass

            if r == None :
                umsg ("Skipping res %d" % i)
                missingRes.append ( i ) 
                continue

            #print "%d/%d - res:" % (i, r.id.position), r.type
            status ("Copying res %d - %s" % (i, r.type) )

            scMol = self.AddResToMol ( r, rPosition, mapxfi, scMol )

            if len ( scMol.residues ) > 1 :
                prevRes = scMol.residues [ -2 ]
                newRes = scMol.residues [ -1 ]
                nb = scMol.newBond ( prevRes.atomsMap["C"][0], newRes.atomsMap["N"][0] )
                nb.display = nb.Stick
                nb.radius = 0.01



        #chimera.openModels.add( [scMol] )

        path, molname = os.path.split ( mol.openedAs[0] )
        fname, fext = os.path.splitext ( molname )

        #scMol.name = fname + "_SC" + fext
        oname = self.scMolName.get()
        scMol.name = oname

        fpath = path + os.path.sep + scMol.name

        print " - saving %s %d.%d" % ( fpath, scMol.id, scMol.subid )
        chimera.PDBio().writePDBfile ( [scMol], fpath )
        #chimera.openModels.add ( [scMol] )

        scm = chimera.PDBio().readPDBfile (fpath)
        chimera.openModels.add ( scm )
        scm[0].name = oname
        scm[0].openedAs = [fpath]

        msg = "Saved to: " + oname
        

        if dev_menus :        
            if len(missingRes) > 0 :
                msg += ", missing"
                for mri in missingRes :
                    msg += " %d" % mri


        umsg (  msg )



    def AddResSC ( self, g, scMol, pathI, rCode, rPosition, resMap ) :


        segMap = segmentation_map()
        if segMap == None :
            self.umsg ( "Please select a map in the Segment Map Dialog" )
            return


        try :
            ress = resMap [ aaCodes[rCode] ]
        except :
            print " * did not find res model"
            return None


        resc = ress[0]
        print " - taking res from ", resc.id.position, resc.type

        prevRes = scMol.residues [ -1 ]


        p1 = g.apoints [ pathI ]

        atCA = resc.atomsMap [ "CA" ][0]
        tp = segMap.openState.xform.inverse().apply ( atCA.xformCoord() )
        p0 = numpy.array ( tp.data() )
        #T0 = chimera.Xform.translation ( chimera.Vector (p0[0], p0[1], p0[2]) )
        #T1 = chimera.Xform.translation ( chimera.Vector (p1[0], p1[1], p1[2]) )
        #Txf = T0.multiply ( T1 )
        
        T0 = numpy.matrix ( [
            [ 1, 0, 0, -p0[0] ],
            [ 0, 1, 0, -p0[1] ],
            [ 0, 0, 1, -p0[2] ],
            [ 0, 0, 0, 1 ]  ] )

        T1 = numpy.matrix ( [
            [ 1, 0, 0, p1[0] ],
            [ 0, 1, 0, p1[1] ],
            [ 0, 0, 1, p1[2] ],
            [ 0, 0, 0, 1 ]  ] )
            
        # direction from CA to N in residue being added
        atN = resc.atomsMap [ "N" ][0]
        tp = segMap.openState.xform.inverse().apply ( atN.xformCoord() )
        v1 = numpy.array ( tp.data() ) - p0
        v1l = numpy.sqrt ( numpy.sum(v1*v1) )
        v1 = v1 / v1l
        
        # direction from CA in path to previous residue C
        atC = prevRes.atomsMap [ "C" ][0]
        tp = segMap.openState.xform.inverse().apply ( atC.xformCoord() )
        v2 = numpy.array ( tp.data() ) - p1
        v2l = numpy.sqrt ( numpy.sum(v2*v2) )
        v2 = v2 / v2l
        
        vAx = numpy.cross ( v1, v2 )
        vAxl = numpy.sqrt ( numpy.sum(vAx*vAx) )
        #print " - vax %.3f" % vAxl

        R_NtoC = None
        if vAxl > 1e-6 :
            vAx = vAx / vAxl
            dotl = numpy.dot ( v1, v2 )
            #print " - dot %f" % dotl
            vAng = numpy.arccos ( dotl )
            xfRot = chimera.Xform.rotation ( chimera.Vector(vAx[0],vAx[1],vAx[2]), vAng * 180.0 / numpy.pi )

            X = ( numpy.matrix (xfRot.getOpenGLMatrix()) ).reshape([4,4]).transpose()
            R_NtoC = numpy.matrix ( [
                [ X[0,0], X[0,1], X[0,2], X[0,3] ],
                [ X[1,0], X[1,1], X[1,2], X[1,3] ],
                [ X[2,0], X[2,1], X[2,2], X[2,3] ],
                [      0,      0,      0, 1 ]  ] )

        else :
            print " - CA-N CA-C already aligned"
            R_NtoC = numpy.matrix ( [
                [ 1, 0, 0, 0 ],
                [ 0, 1, 0, 0 ],
                [ 0, 0, 1, 0 ],
                [ 0, 0, 0, 1 ]  ] )



        minSumG = 1e8
        minTxf = None
        N = 10
        for angi in range (N) :

            # rotate around CA - previous C axis to find best side chain placement
            xfRot2 = chimera.Xform.rotation ( chimera.Vector(v2[0],v2[1],v2[2]), angi*360/N )
    
            X = ( numpy.matrix (xfRot2.getOpenGLMatrix()) ).reshape([4,4]).transpose()
            R_ax = numpy.matrix ( [
                [ X[0,0], X[0,1], X[0,2], X[0,3] ],
                [ X[1,0], X[1,1], X[1,2], X[1,3] ],
                [ X[2,0], X[2,1], X[2,2], X[2,3] ],
                [      0,      0,      0, 1 ]  ] )
    
            T = T1 * R_ax * R_NtoC * T0
    
            Txf = chimera.Xform.xform ( T[0,0], T[0,1], T[0,2], T[0,3], T[1,0], T[1,1], T[1,2], T[1,3], T[2,0], T[2,1], T[2,2], T[2,3], True )
            Txf.multiply ( segMap.openState.xform.inverse() )

            sumg = 0.0
            for i, at in enumerate ( resc.atoms ) :
                #tp = segMap.openState.xform.inverse().apply ( at.xformCoord() )
                xpt = numpy.array ( Txf.apply ( at.xformCoord() ).data() )
                sumg += self.GradMagnAtPosFromTree ( xpt, g.CAPointsTree, cutoff=10, F=1 )
                sumg += self.GradMagnAtPosFromTree ( xpt, g.SCPointsTree, cutoff=10, F=1 )
                sumg += self.GradMagnAtPosFromTree ( xpt, g.newResPointsTree, cutoff=10, F=1 )
            
            #print " - %d (%d) - %f" % (angi, angi*360/N, sumg)
            
            if sumg < minSumG :
                minSumG = sumg
                minTxf = Txf


        scMol = self.AddResToMol ( resc, rPosition, minTxf, scMol, False )
        
        return scMol.residues [ -1 ]





    def AddResToMol ( self, res, posi, xf, nmol, keepSS = True ) :

        if nmol == None :
            nmol = chimera.Molecule()
            nmol.name = "SCMOL"

        aMap = dict()
        clr = ( .7, .7, .7, 1 )

        nres = nmol.newResidue (res.type, chimera.MolResId(res.id.chainId, posi))
        # print "New res: %s %d" % (nres.id.chainId, nres.id.position)
        for at in res.atoms :
            atname = at.name
            if at.name == "OT1" :
                atname = "O"
            elif at.name == "OT2" :
                continue
            elif "HT1" in at.name or "HT2" in at.name or "HT3" in at.name :
                continue

            nat = nmol.newAtom (atname, chimera.Element(at.element.number))
            aMap[at] = nat 
            nres.addAtom( nat )

            atpos = at.xformCoord()
            if xf != None :
                atpos = xf.apply ( atpos )

            nat.setCoord ( atpos )
            nat.drawMode = nat.Ball
            nat.radius = 1.5
            nat.color = chimera.MaterialColor( clr[0], clr[1], clr[2], 1.0 )
            nat.display = True

        if keepSS :
            nres.isHelix = res.isHelix
            nres.isHet = res.isHet
            nres.isSheet = res.isSheet
            nres.isStrand = res.isStrand
        else :
            nres.isHelix = False
            nres.isHet = False
            nres.isSheet = False
            nres.isStrand = False
        
        nres.ribbonDisplay = True
        nres.ribbonDrawMode = 2
        nres.ribbonColor = chimera.MaterialColor( clr[0], clr[1], clr[2], 1.0 );
    
        for bond in res.molecule.bonds :
            try :
                nb = nmol.newBond ( aMap[bond.atoms[0]], aMap[bond.atoms[1]] )
            except :
                continue
            nb.display = nb.Stick
            nb.radius = 0.2
    
        
        return nmol



    def ShowRegionAxesSelected ( self ) :
 
        smod = self.CurrentSegmentation()
        if smod == None : return

        sregs = smod.selected_regions()        
        if len(sregs)==0 : print "no selected regions found"; return

        self.ShowRegionsAxes ( sregs )




    def ShowRegionsAxes ( self, regs ) :

        smod = self.CurrentSegmentation()
        if smod is None: return
        
        adjExt = float ( self.helixLength.get() )
        widthF = float ( self.helixWidthF.get() )

        for r in regs :

            sp = r.surface_piece
            try :
                sp.axes.display = True
                chimera.openModels.close ( sp.axes )
            except :
                pass

            tpoints = r.map_points()
            sp.COM, sp.U, sp.S, sp.V = prAxes ( tpoints )

            com = numpy.sum(tpoints, axis=0) / len(tpoints)
            comv = numpy.ones_like ( tpoints ) * com
            points = tpoints - comv

            ppoints = points * sp.U
            sp.Extents = numpy.asarray ( numpy.max ( numpy.abs ( ppoints ), 0 ) )[0] + adjExt

            sp.Extents[0] += 5.0
            sp.Extents[1] += 5.0
            sp.Extents[2] += 5.0

            import axes
            reload (axes)

            if 0 :
            # for ribosome direction
                sp.Extents[1] = sp.Extents[1] * float(self.axesFactor.get())

                sp.axes = axes.AxesMod ( sp.COM, sp.U, sp.Extents, 6, 1.0,alignTo = sp.model )
            else :
                sp.axes = axes.AxesMod ( sp.COM, sp.U, sp.Extents, 1.0, 1.1, alignTo = sp.model )

            sp.axes.name = "region_%d_axes" % r.rid
            
            

    def RegionSizes ( self ) :

        segMap = segmentation_map()
        if segMap == None :
            self.umsg ( "Please select a map in the Segment Map Dialog" )
            return
        
        smod = current_segmentation ()
        if smod == None :
            self.umsg ( "Please select a Current Segmentation in the Segment Map dialog" )
            return
        
        print "Seg has %d regions" % (len(smod.regions))


        nregs, last = len(smod.regions), 0
        regs = list(smod.regions)
        #distByReg = {}
        sizeByReg = {}
        
        for ri, r in enumerate ( regs ) :

            sizeByReg[r] = r.point_count()
            
            at = int(numpy.floor( 10.0 * (ri+1) / nregs ))
            if at > last :
                print at,
                last = at; at += 1

        print ""

        dists = sizeByReg.values ()
        maxDist = max (dists) + 0.01
        minDist = min (dists)
        nbins = int ( self.numBins.get() )
        dr = (maxDist - minDist) / float(nbins)
        print "%d dists - max %.2f, min %.2f, nb %d, dr %.2f" % (len(dists), maxDist, minDist, nbins, dr)

        bins = []
        for i in range (nbins) :
            bins.append ( [] )

        print "bad bins: ",
        for regm, rad in sizeByReg.iteritems() :
            bini = int ( numpy.floor ( (rad - minDist) / dr ) )
            if bini >= len(bins) :
                print bini,
                bini = len(bins)-1
            bins[bini].append (regm)

        print ""

        f = open ( "sizes.txt", "w" )
        for k,regs in enumerate ( bins ) :
            v = len(regs)
            vmin = minDist + k * dr
            vmax = minDist + (k+1) * dr
            rm = .5 * (vmin + vmax)
            vn = v / (4 * 3.14 * rm * rm)
            f.write ( "%d\t%.2f\t%.2f\t%d\t%f\n" % (k, vmin, vmax, v, vn) )
        f.close()
        
		


    def SelectRegions ( self ) :
        
        segMap = segmentation_map()
        if segMap == None :
            umsg ( "Please select a map in the Segment Map Dialog" )
            return

        smod = current_segmentation ()
        if smod == None :
            umsg ( "Please select a segmentation file in the Segment Map dialog" )
            return


        regs = list(smod.regions)
        oregs = []
        mins = int ( self.minSize.get() )
        maxs = int ( self.maxSize.get() )
        
        for ri, r in enumerate ( regs ) :
            i = r.point_count()
            if i >= mins and i <= maxs :
                oregs.append ( r )
        
        if len(oregs) > 10000 :
            umsg ( "Too many regions to select... %d" % ( len(oregs) ) )
            return

        regions.select_regions ( oregs )
        umsg ( "Selected %d regions" % ( len(oregs) ) )


    def FitsSort ( self ) :

        segMap = segmentation_map()
        if segMap == None :
            umsg ( "Please select a map in the Segment Map Dialog" )
            return

        import fit_dialog
        dlg = fit_dialog.fit_segments_dialog ()
        
        
        f = []
        ccs = []
        for fmap, dmap, mat, corr, aI, bI, bC, bO, regions in dlg.list_fits :
            f.append ( [corr, fmap] )
            ccs.append ( corr )

        f.sort ( reverse=True, key=lambda x: x[0] )

        min_cc = min ( ccs )
        max_cc = max ( ccs )

        fp = open ( "score_fits.txt", "w" );

        for corr, fmap in f :
            #print corr, fmap.name
            fp.write ( "%f\t%s\n" % (corr, fmap.name) )
            
            cf = (corr - min_cc) / (max_cc - min_cc)
            rf = 1.0 - cf
            
            for sp in fmap.surfacePieces :
                v, t = sp.geometry
                if len(v) == 8 and len(t) == 12 :
                    sp.display = False
                else :
                    c = sp.color
                    sp.color = ( .7*rf+.3, .3, .7*cf+.3, 1 )
                    
                    if hasattr ( sp, "vertexColors" ) and sp.vertexColors != None :
                        vcolors = []
                        for vc in sp.vertexColors :
                            vcolors.append ( sp.color )
    
                        sp.vertexColors = vcolors
            

        fp.close ()
        

    def FitsToOrigin ( self ) :

        segMap = segmentation_map()
        if segMap == None :
            umsg ( "Please select a map in the Segment Map Dialog" )
            return

        import fit_dialog
        dlg = fit_dialog.fit_segments_dialog ()

        for fmap, dmap, mat, corr, aI, bI, bC, bO, regions in dlg.list_fits :
            fmap.openState.xform = segMap.openState.xform


    def FitsToFit ( self ) :

        segMap = segmentation_map()
        if segMap == None :
            umsg ( "Please select a map in the Segment Map Dialog" )
            return

        import fit_dialog
        dlg = fit_dialog.fit_segments_dialog ()

        for fmap, dmap, mat, corr, aI, bI, bC, bO, regions in dlg.list_fits :
            self.place_map ( fmap, fmap.M, segMap )


    def FitsHide ( self ) :

        segMap = segmentation_map()
        if segMap == None :
            umsg ( "Please select a map in the Segment Map Dialog" )
            return

        import fit_dialog
        dlg = fit_dialog.fit_segments_dialog ()
        
#         from VolumeViewer import Volume
#         mlist = OML(modelTypes = [Volume])
# 
#         fmap = None
#         avgMat = None
#         N = 0.0
# 
#         for m in mlist :
#             if m.display == True :
#                 print m.name

        import Surface

        if self.what.get() == "selected" :
            for fmap, dmap, mat, corr, aI, bI, bC, bO, regions in dlg.list_fits :
                for sp in fmap.surfacePieces :
                    if sp in Surface.selected_surface_pieces() :
                        fmap.display = False

        if self.what.get() == "all" :
            for fmap, dmap, mat, corr, aI, bI, bC, bO, regions in dlg.list_fits :
                fmap.display = False
                    
        if self.what.get() == "visible" :
            for fmap, dmap, mat, corr, aI, bI, bC, bO, regions in dlg.list_fits :
                if fmap.display == True :
                    fmap.display = False
                
                
    def FitsShow ( self ) :

        segMap = segmentation_map()
        if segMap == None :
            umsg ( "Please select a map in the Segment Map Dialog" )
            return

        import fit_dialog
        dlg = fit_dialog.fit_segments_dialog ()
        
        if self.what.get() == "selected" :
            for fmap, dmap, mat, corr, aI, bI, bC, bO, regions in dlg.list_fits :
                for sp in fmap.surfacePieces :
                    if sp in Surface.selected_surface_pieces() :
                        fmap.display = True
                
        if self.what.get() == "all" :
            for fmap, dmap, mat, corr, aI, bI, bC, bO, regions in dlg.list_fits :
                fmap.display = True
                    
        if self.what.get() == "visible" :
            for fmap, dmap, mat, corr, aI, bI, bC, bO, regions in dlg.list_fits :
                if fmap.display == True :
                    fmap.display = True


    def SetColor ( self ) :

        c = self.setColorVar.get().split(",")
        if len(c) != 4 :
            print "bad color"
            return
            
        clr = ( float(c[0]), float(c[1]), float(c[2]), float(c[3]) )


        segMap = segmentation_map()
        if segMap == None :
            umsg ( "Please select a map in the Segment Map Dialog" )
            return

        import fit_dialog
        dlg = fit_dialog.fit_segments_dialog ()


        for fmap, dmap, mat, corr, aI, bI, bC, bO, regions in dlg.list_fits :

            if self.what.get() == "selected" and self.isSelected (fmap) == False :
                continue

            if self.what.get() == "visible" and fmap.display == False :
                continue

            for sp in fmap.surfacePieces :
                v, t = sp.geometry
                if len(v) == 8 and len(t) == 12 :
                    sp.display = False
                else :
                    sp.color = clr
                    if hasattr ( sp, "vertexColors" ) and sp.vertexColors != None :
                        vcolors = []
                        for vc in sp.vertexColors :
                            vcolors.append ( clr )
                        sp.vertexColors = vcolors


    def isSelected ( self, fmap ) :
        for sp in fmap.surfacePieces :
            if sp in Surface.selected_surface_pieces() :
                return True
        return False


    def place_map(self, fmap, mat, dmap):

        import numpy
        tf = numpy.array(mat[:3,:])
        try :
            xf = dmap.openState.xform
        except :
            print "Reference map no longer open"
            return

        from Matrix import xform_matrix, multiply_matrices, chimera_xform, identity_matrix, invert_matrix, shift_and_angle
        xf.multiply(chimera_xform(tf))

        try :
            fmap.openState.xform = xf
        except :
            print "Fitted map no longer open"
            return

        fmap.M = mat

        for mol in fmap.mols :
            mol.openState.xform = xf
        




class Res :
    def __init__ (self) :
        self.pos = None





class Helix :
    def __init__ (self) :
        self.mod = None


    def MakeMod ( self, alignTo = None ) :
        
        import _surface
        surf_mod = _surface.SurfaceModel()
        chimera.openModels.add([surf_mod], sameAs = alignTo)

        import axes
        reload (axes)

        cyl = axes.AddCylinderSolid ( chimera.Vector(0,0,0), chimera.Vector(0,0,1), 1.5*self.Extents[2]+self.heightAdj, (1,.2,.2,1), self.widthF, surf_mod )
        cyl.name = "Helix"
        
        U = self.U
        COM = self.COM
        if U != None :
            R = numpy.array([
                             [  U[0,0], U[0,1], U[0,2], 0.0    ],
                             [  U[1,0], U[1,1], U[1,2], 0.0    ],
                             [  U[2,0], U[2,1], U[2,2], 0.0    ]  ] )
            
            T = numpy.array([
                             [  1.0, 0.0, 0.0, COM[0]   ],
                             [  0.0, 1.0, 0.0, COM[1]   ],
                             [  0.0, 0.0, 1.0, COM[2]   ]  ] )
            
            Ti = numpy.array([
                              [  1.0, 0.0, 0.0, 0  ],
                              [  0.0, 1.0, 0.0, 0   ],
                              [  0.0, 0.0, 1.0, -0.75*self.Extents[2]-self.heightAdj/2.0   ]  ] )
            
            import Matrix
            M = Matrix.multiply_matrices ( R, Ti )
            M = Matrix.multiply_matrices ( T, M )
            

            ps = []
            for p in surf_mod.surfacePieces :
                v, t = numpy.copy(p.geometry[0]), numpy.copy(p.geometry[1])
                ps.append ( [v,t,p.color] )
                surf_mod.removePiece ( p )

            import _contour
            for p in ps :
                _contour.affine_transform_vertices( p[0], M )
                surf_mod.addPiece ( p[0], p[1], p[2] )

            from random import random as rand
            clr = ( rand()*.7, rand()*.7, rand()*.7, 1.0 )
            # for p in axes.surfacePieces : p.color = clr

        #for g in axes.surfacePieces :
        #    g.initial_v = numpy.copy ( g.geometry[0] )
        




def segloop_dialog ( create=False ) :

  from chimera import dialogs


def close_dialog ():

    from chimera import dialogs



def show_dialog ():

    from chimera import dialogs

    d = dialogs.find ( "segger_segloop", create=False )
    if d :
        print " - found old diag"
        d.toplevel_widget.update_idletasks ()
        d.Close()
        d.toplevel_widget.update_idletasks ()

    dialogs.register (Segloop_Dialog.name, Segloop_Dialog, replace = True)

    d = dialogs.find ( "segger_segloop", create=True )
    # Avoid transient dialog resizing when created and mapped for first time.
    d.toplevel_widget.update_idletasks ()
    d.enter()

    return d



def GetMod ( name ) :
    for m in chimera.openModels.list() :
        if m.name == name :
            return m
    return None



# -----------------------------------------------------------------------------
#

