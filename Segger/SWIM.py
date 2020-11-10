
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
import ttk
import tkFont
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

from axes import prAxes
import regions
import graph
from Segger import dev_menus, timing, seggerVersion
from CGLutil.AdaptiveTree import AdaptiveTree

from chimera.resCode import nucleic3to1
from chimera.resCode import protein3to1


devMenus = False


import qscores
reload (qscores)


chargedIons = { "MG":2, "NA":1, "CL":-1, "CA":2, "ZN":2, "MN":2, "FE":3, "CO":2, "NI":2 }

atomColors = {'C' : chimera.MaterialColor (0.565,0.565,0.565),
            'Cbb' : chimera.MaterialColor (0.2,0.6,0.2),
            'S' : chimera.MaterialColor (1.000,1.000,0.188),
            'O' : chimera.MaterialColor (1.000,0.051,0.051),
            'N' : chimera.MaterialColor (0.188,0.314,0.973),
            'P' : chimera.MaterialColor (1.0, 0.502, 0.0),
            'H' : chimera.MaterialColor (0.9,.9,.9),
            ' ' : chimera.MaterialColor (0.2,1,.2),
            "MG" : chimera.MaterialColor (0,1,0),
            "NA" : chimera.MaterialColor (.7,.4,.9),
            "CL" : chimera.MaterialColor (.95,.59,.21), # orange
            "CA" : chimera.MaterialColor (0,1,0),
            "ZN" : chimera.MaterialColor (.52,.60,.25), # dark green
            "MN" : chimera.MaterialColor (0,1,0),
            "FE" : chimera.MaterialColor (.42,.48,.27), # turquise
            "CO" : chimera.MaterialColor (0,1,0),
            "NI" : chimera.MaterialColor (0,1,0)
}



from segment_dialog import current_segmentation, segmentation_map


def umsg ( txt ) :
    print txt
    status ( txt )

def status ( txt ) :
    txt = txt.rstrip('\n')
    msg.configure(text = txt)
    msg.update_idletasks()


# https://android.googlesource.com/toolchain/python/+/243b47fbef58ab866ee77567f2f52affd8ec8d0f/Python-2.7.3/Demo/tkinter/ttk/treeview_multicolumn.py


class SWIM_Dialog ( chimera.baseDialog.ModelessDialog ):

    title = "SWIM: Segment-guided Water and Ion Modeling"
    name = "swim"

    buttons = ("Thresholds", "Go", "Stats", "Options", "Log", "Close")


    help = 'https://cryoem.slac.stanford.edu/ncmi/resources/software'

    def fillInUI(self, parent):

        self.group_mouse_mode = None

        tw = parent.winfo_toplevel()
        self.toplevel_widget = tw
        tw.withdraw()

        parent.columnconfigure(0, weight = 1)
        #parent.columnconfigure(1, weight = 1)

        row = 0

        menubar = Tkinter.Menu(parent, type = 'menubar', tearoff = False)
        tw.config(menu = menubar)




        if 1 :
            ff = Tkinter.Frame(parent)
            ff.grid(column=0, row=row, sticky='w')

            #Tkinter.Grid.columnconfigure(parent, 0, weight=1)
            #Tkinter.Grid.columnconfigure(ff, 0, weight=1)


            l = Tkinter.Label(ff, text='Map: ')
            l.grid(column=0, row=0, sticky='w')

            self.cur_dmap = None
            self.dmap = Tkinter.StringVar(parent)

            self.mb  = Tkinter.Menubutton ( ff, textvariable=self.dmap, relief=Tkinter.RAISED )
            self.mb.grid (column=1, row=0, sticky='we', padx=5)
            self.mb.menu  =  Tkinter.Menu ( self.mb, tearoff=0, postcommand=self.MapMenu )
            self.mb["menu"]  =  self.mb.menu

            ff.columnconfigure(1, weight=1)

            self.cur_dmap = None
            self.SetVisMap ()


        if 1 :
            row += 1

            ff = Tkinter.Frame(parent)
            ff.grid(column=0, row=row, sticky='w') # put we to stretch

            l = Tkinter.Label(ff, text='Model:', anchor=Tkinter.W)
            l.grid(column=0, row=0, sticky='w')

            self.struc = Tkinter.StringVar(parent)
            self.strucMB  = Tkinter.Menubutton ( ff, textvariable=self.struc, relief=Tkinter.RAISED )
            self.strucMB.grid (column=1, row=0, sticky='we', padx=1)
            self.strucMB.menu  =  Tkinter.Menu ( self.strucMB, tearoff=0, postcommand=self.StrucMenu )
            self.strucMB["menu"]  =  self.strucMB.menu

            ff.columnconfigure(1, weight=1)

            self.cur_mol = None
            self.cur_chains = []


        if 1 :

            row += 1

            ff = Tkinter.Frame(parent)
            ff.grid(column=0, row=row, sticky='news')

            self.id_mg = {}

            self.tree = ttk.Treeview(ff)

            #self.tree["columns"]=("one","two","three")
            self.tree.column("#0", width=400, minwidth=200, stretch=Tkinter.YES)
            #self.tree.column("one", width=150, minwidth=150, stretch=Tkinter.NO)
            #self.tree.column("two", width=400, minwidth=200)
            #self.tree.column("three", width=80, minwidth=50, stretch=Tkinter.NO)

            self.tree.heading("#0",text="Chain,Residue,Atom",anchor=Tkinter.W)
            #self.tree.heading("one", text="Date modified",anchor=Tkinter.W)
            #self.tree.heading("two", text="Type",anchor=Tkinter.W)
            #self.tree.heading("three", text="Size",anchor=Tkinter.W)

            #self.tree.pack(side=Tkinter.TOP,fill=Tkinter.X)
            #self.tree.grid(column=0, row=0, sticky='nsew')
            #self.tree.pack(fill=Tkinter.BOTH, expand=1)
            #tree.place(x=0, y=0, relwidth=1, relheight=1)

            self.tree.grid(row = 0, column = 0, sticky='news')
            parent.columnconfigure(0, weight=1)
            parent.rowconfigure(row, weight = 1)
            ff.rowconfigure(0, weight = 1)
            ff.columnconfigure(0, weight=1)

            self.tree.bind('<<TreeviewSelect>>', self.select_mg_cb)
            self.tree.bind('<<TreeviewOpen>>', self.open_mg_cb)
            self.tree.bind('<<TreeviewClose>>', self.close_mg_cb)


        if 1 :
            row += 1
            ff = Tkinter.Frame(parent)
            ff.grid(column=0, row=row, sticky='w')

            b = Tkinter.Button(ff, text="Refresh", command=self.RefreshTree)
            b.grid (column=0, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="Select", command=self.SelectSel)
            b.grid (column=1, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="Select All", command=self.SelectAll)
            b.grid (column=2, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="Show", command=self.ShowSel)
            b.grid (column=3, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="Show Only", command=self.ShowSelOnly)
            b.grid (column=4, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="Show All", command=self.ShowAll)
            b.grid (column=5, row=0, sticky='w', padx=0, pady=1)

            #b = Tkinter.Button(ff, text="Avg", command=self.Average)
            #b.grid (column=1, row=0, sticky='w', padx=0, pady=1)

            #b = Tkinter.Button(ff, text="Open", command=self.Open)
            #b.grid (column=2, row=0, sticky='w', padx=0, pady=1)





        if 1 :

            row += 1
            dummyFrame = Tkinter.Frame(parent, relief='groove', borderwidth=1)
            Tkinter.Frame(dummyFrame).pack()
            dummyFrame.grid(row=row,column=0,columnspan=7, pady=2, sticky='we')

            row += 1
            ff = Tkinter.Frame(parent)
            ff.grid(column=0, row=row, sticky='w')

            um = Hybrid.Checkbutton(ff, 'Place with Mouse (Ctrl+Right Click on Blob)', False)
            um.button.grid(column = 1, row=0, sticky = 'w', padx=5)
            self.use_mouse = um.variable
            um.callback(self.bind_placement_button_cb)


            b = Tkinter.Label(ff, text="Add To Chain:")
            b.grid (column=8, row=0, sticky='w', padx=0, pady=1)

            self.addToChain = Tkinter.StringVar(ff)
            #self.addRess.set ( "vsgtngtkrf" )
            self.addToChain.set ( "" )
            e = Tkinter.Entry(ff, width=2, textvariable=self.addToChain)
            e.grid(column=9, row=0, sticky='w', padx=5, pady=1)



            um = Hybrid.Checkbutton(ff, 'Find Peak', False)
            um.button.grid(column = 10, row=0, sticky = 'w', padx=5)
            self.use_mouse_max = um.variable
            self.use_mouse_max.set(True)
            #um.callback(self.bind_placement_button_cb)


        if 1 :

            row += 1
            ff = Tkinter.Frame(parent)
            ff.grid(column=0, row=row, sticky='w')


            self.guessOpt = Tkinter.StringVar()
            self.guessOpt.set ( 'guess' )

            l = Tkinter.Label(ff, text=' ' )
            l.grid(column=0, row=0, sticky='w')

            c = Tkinter.Radiobutton(ff, text="Guess", variable=self.guessOpt, value = 'guess')
            c.grid (column=1, row=0, sticky='w')


            c = Tkinter.Radiobutton(ff, text="Add:", variable=self.guessOpt, value = 'opt')
            c.grid (column=2, row=0, sticky='w')


            self.addStr = Tkinter.StringVar(ff)
            #self.addRess.set ( "vsgtngtkrf" )
            self.addStr.set ( "W" )
            e = Tkinter.Entry(ff, width=5, textvariable=self.addStr)
            e.grid(column=3, row=0, sticky='w', padx=5, pady=1)

            b = Tkinter.Label(ff, text=" e.g. Mg, Ca, Na, Zn, Mn, Fe, Co, Ni (Ion), W (Water)")
            b.grid (column=4, row=0, sticky='w', padx=0, pady=1)



            #um = Hybrid.Checkbutton(ff, 'Guess', False)
            #um.button.grid(column = 9, row=0, sticky = 'w', padx=5)
            #self.use_mouse_guess = um.variable
            #self.use_mouse_guess.set(True)
            #um.callback(self.bind_placement_button_cb)


            #b = Tkinter.Button(ff, text="SWIM", command=self.Place)
            #b.grid (column=5, row=0, sticky='w', padx=5)






        if devMenus :

            row += 1
            ff = Tkinter.Frame(parent)
            ff.grid(column=0, row=row, sticky='w')


            #b = Tkinter.Button(ff, text="Stats", command=self.Hoh)
            #b.grid (column=2, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(ff, text="ShoW", command=self.HohShow)
            #b.grid (column=46, row=0, sticky='w', padx=5)


            #b = Tkinter.Button(ff, text="Thr", command=self.Thr)
            #b.grid (column=4, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(ff, text="WEx", command=self.HohE)
            #b.grid (column=5, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(ff, text="Asn", command=self.Asn)
            #b.grid (column=51, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Dw", command=self.HohD)
            b.grid (column=6, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Di", command=self.IonD)
            b.grid (column=7, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Da", command=self.AllD)
            b.grid (column=8, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Dwn", command=self.HohDn)
            b.grid (column=9, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Comb", command=self.Combine)
            b.grid (column=10, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Dup", command=self.Duplicates)
            b.grid (column=11, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(ff, text="RMSD", command=self.RMSD)
            #b.grid (column=10, row=0, sticky='w', padx=5)


        if devMenus :
            row += 1
            ff = Tkinter.Frame(parent)
            ff.grid(column=0, row=row, sticky='w')

            b = Tkinter.Label(ff, text="Map Resolution:")
            b.grid (column=1, row=0, sticky='w', padx=0, pady=1)

            self.mapRes = Tkinter.StringVar(ff)
            #self.addRess.set ( "vsgtngtkrf" )
            self.mapRes.set ( "" )
            e = Tkinter.Entry(ff, width=5, textvariable=self.mapRes)
            e.grid(column=2, row=0, sticky='w', padx=5, pady=1)

            b = Tkinter.Label(ff, text="Angstroms")
            b.grid (column=3, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="Estimate", command=self.GuessRes)
            b.grid (column=4, row=0, sticky='w', padx=5)

            b = Tkinter.Label(ff, text="using")
            b.grid (column=5, row=0, sticky='w', padx=0, pady=1)

            self.mapResN = Tkinter.StringVar(ff)
            #self.addRess.set ( "vsgtngtkrf" )
            self.mapResN.set ( "20" )
            e = Tkinter.Entry(ff, width=5, textvariable=self.mapResN)
            e.grid(column=6, row=0, sticky='w', padx=5, pady=1)

            b = Tkinter.Label(ff, text="atoms")
            b.grid (column=7, row=0, sticky='w', padx=0, pady=1)


        if 1 :

            row += 1
            op = Hybrid.Popup_Panel(parent)
            df = op.frame
            df.grid(column = 0, row = row, sticky = 'w')
            #df.grid_remove()
            #ff.columnconfigure(0, weight=1)
            self.optionsPanel = op.panel_shown_variable

            orow = 0
            dummyFrame = Tkinter.Frame(df, relief='groove', borderwidth=1)
            Tkinter.Frame(dummyFrame).pack()
            dummyFrame.grid(row=0,column=orow,columnspan=7, pady=2, sticky='we')

            if 1 :
                orow += 1
                ff = Tkinter.Frame(df)
                ff.grid(column=0, row=orow, sticky='w')

                b = Tkinter.Label(ff, text="Distance ranges (in Angstroms):")
                b.grid (column=0, row=0, sticky='w', padx=0, pady=1)

                orow += 1
                ff = Tkinter.Frame(df)
                ff.grid(column=0, row=orow, sticky='w')

                b = Tkinter.Label(ff, text="   Ions - Min:")
                b.grid (column=1, row=0, sticky='w', padx=0, pady=1)

                self.ionMinD = Tkinter.StringVar(ff)
                #self.addRess.set ( "vsgtngtkrf" )
                self.ionMinD.set ( "1.8" )
                e = Tkinter.Entry(ff, width=5, textvariable=self.ionMinD)
                e.grid(column=2, row=0, sticky='w', padx=5, pady=1)


                b = Tkinter.Label(ff, text=" Max:")
                b.grid (column=3, row=0, sticky='w', padx=0, pady=1)

                self.ionMaxD = Tkinter.StringVar(ff)
                #self.addRess.set ( "vsgtngtkrf" )
                self.ionMaxD.set ( "2.5" )
                e = Tkinter.Entry(ff, width=5, textvariable=self.ionMaxD)
                e.grid(column=4, row=0, sticky='w', padx=5, pady=1)



                #b = Tkinter.Label(ff, text="Distances (in Angstroms):")
                #b.grid (column=0, row=0, sticky='w', padx=0, pady=1)

                b = Tkinter.Label(ff, text="   Water - Min:")
                b.grid (column=1, row=1, sticky='w', padx=0, pady=1)

                self.waterMinD = Tkinter.StringVar(ff)
                #self.addRess.set ( "vsgtngtkrf" )
                self.waterMinD.set ( "2.5" )
                e = Tkinter.Entry(ff, width=5, textvariable=self.waterMinD)
                e.grid(column=2, row=1, sticky='w', padx=5, pady=1)


                b = Tkinter.Label(ff, text=" Max:")
                b.grid (column=3, row=1, sticky='w', padx=0, pady=1)

                self.waterMaxD = Tkinter.StringVar(ff)
                #self.addRess.set ( "vsgtngtkrf" )
                self.waterMaxD.set ( "3.3" )
                e = Tkinter.Entry(ff, width=5, textvariable=self.waterMaxD)
                e.grid(column=4, row=1, sticky='w', padx=5, pady=1)



            if 1 :
                orow += 1
                ff = Tkinter.Frame(df)
                ff.grid(column=0, row=orow, sticky='w')

                um = Hybrid.Checkbutton(ff, 'Put water/ion only when Q-score >', False)
                um.button.grid(column = 1, row=0, sticky = 'w', padx=5)
                self.useQScore = um.variable
                #um.callback(self.bind_placement_button_cb)
                self.useQScore.set(False)

                #b = Tkinter.Label(ff, text="Put water/ion only when Q-score >")
                #b.grid (column=1, row=0, sticky='w', padx=0, pady=1)

                self.placeQ = Tkinter.StringVar(ff)
                #self.addRess.set ( "vsgtngtkrf" )
                self.placeQ.set ( "0.9" )
                e = Tkinter.Entry(ff, width=5, textvariable=self.placeQ)
                e.grid(column=2, row=0, sticky='w', padx=5, pady=1)

                self.qsigma = Tkinter.StringVar(ff)
                self.qsigma.set ( "0.4" )

                if 1 :
                    b = Tkinter.Label(ff, text=" sigma:")
                    b.grid (column=3, row=0, sticky='w', padx=0, pady=1)

                    e = Tkinter.Entry(ff, width=5, textvariable=self.qsigma)
                    e.grid(column=4, row=0, sticky='w', padx=5, pady=1)



        self.optionsPanel.set(False)
        #self.selPanel.set(False)



        row += 1
        dummyFrame = Tkinter.Frame(parent, relief='groove', borderwidth=1)
        Tkinter.Frame(dummyFrame).pack()
        dummyFrame.grid(row=row,column=0,columnspan=7, pady=2, sticky='we')


        row += 1
        global msg
        msg = Tkinter.Label(parent, width = 60, anchor = 'w', justify = 'left', fg="red")
        msg.grid(column=0, row=row, sticky='ew')
        self.msg = msg
        msg.configure(text = "Let's Go SWIMming")


        self.SelectedMgId = None
        self.SetVisMol ()



        callbacks = (self.mouse_down_cb, self.mouse_drag_cb, self.mouse_up_cb)
        #callbacks = (self.mouse_down_cb)
        from chimera import mousemodes
        mousemodes.addFunction('mark swim', callbacks, self.mouse_mode_icon())

        if 1 :
            # bind, unbind in case it was left bound before...
            from chimera import mousemodes
            print " - unbinding mouse..."
            button, modifiers = ('3', ['Ctrl'])
            def_mode = mousemodes.getDefault(button, modifiers)
            mousemodes.setButtonFunction(button, modifiers, def_mode)
            self.bound_button = None





    def Options ( self ) :
        self.optionsPanel.set (not self.optionsPanel.get())


    def Log ( self ) :
        import Idle
        Idle.start_shell()

    def bind_placement_button_cb(self) :

        if self.use_mouse.get() :
            print " - binding mouse..."
            button, modifiers = ('3', ['Ctrl'])
            from chimera import mousemodes
            mousemodes.setButtonFunction(button, modifiers, 'mark swim')
            self.bound_button = (button, modifiers)
        elif self.bound_button:
            print " - unbinding mouse..."
            button, modifiers = self.bound_button
            from chimera import mousemodes
            def_mode = mousemodes.getDefault(button, modifiers)
            mousemodes.setButtonFunction(button, modifiers, def_mode)
            self.bound_button = None


    def mouse_mode_icon(self) :

        import os.path
        icon_path = os.path.join(os.path.dirname(__file__), 'marker.gif')
        from PIL import Image
        image = Image.open(icon_path)
        from chimera import chimage
        from chimera import tkgui
        icon = chimage.get(image, tkgui.app)
        return icon

    def mouse_down_cb(self, viewer, event) :

        print " mouse - "

        #print event.x, event.y
        if 0 :
            print dir(event)
            print event.char
            print event.keycode
            print event.keysym
            print event.keysym_num
            print event.num
            print event.state

        hits = []
        import VolumePath.tracer as tracer

        if 1 :
            from VolumeViewer import volume_list
            hits.extend(tracer.volume_maxima(event.x, event.y, volume_list()))
            print "vol"

        if 0 :
            from VolumeViewer import volume_list
            hits.extend(VolumePath.tracer.volume_plane_intercepts(event.x, event.y, volume_list()))

        if 0 :
            from Surface import surface_models
            hits.extend(tracer.surface_intercepts(event.x, event.y, surface_models()))
            print "surf"

        for C, vol in hits :
            print " --> ", vol.name, " --> %.1f, %.1f, %.1f" % (C[0], C[1], C[2])
            self.PlaceAt ( C, vol )





        #grabbed = (self.move_markers.get() and self.grab_marker(event.x, event.y))
        #if not grabbed:
        #    self.add_marker_at_screen_xy(event.x, event.y)



    def mouse_drag_cb(self, viewer, event):
        shift_mask = 1
        shift = (event.state & shift_mask)
        capslock_mask = 2
        capslock = (event.state & capslock_mask)
        #self.move_or_resize_marker(event.x, event.y, shift, capslock):


    def mouse_up_cb(self, viewer, event):
        #self.ungrab_marker()
        #self.pause_marker_placement = False
        #print "mouse up"
        pass


# To distinguish between water and ions, we use the criteria in the Undowser paper. The paper does not describe exact distances, but the criteria is as follows:
#   1. A placed water that clashes with two or more atoms of the same polarity, and with no nonpolars (C) or opposite polars (O and N), is almost certainly an ion.
#   2. If the 'placed water' clashes (is too close) to negative atoms it is a positive ion
#   3. If the 'placed water' clashes with positive atoms, it is a negative ion
#   4. A doubly charged ion (e.g. Mg++) almost always interacts with at least one fully charged atom (e.g. phosphate or carboxyl O)
#   5. A singly charged ion (e.g. Na+) often interacts with just partial charges (e.g. OH, backbone CO)
# We add to these criteria our previous observations in a high-resolution X-ray structure (pdb:3ajo), the following distances are observed:
#   1. Water atom to nearby polar atoms: 2.8A +/- ~0.4
#   2. Ion to nearby charged/polar atoms: 2.2A +/- ~0.2
# Thus we:
#   - Place ion when distance is 2.2A +/- 0.2 and nearby atoms include
#     * Charged atoms (place double charged ion such as Mg++)
#     * Single or opposite polars (O and N) (place single charged atom such as Na+)
#   - Place water when distance is 2.8 +/- 0.4 and nearby atoms include
#     * Charged atom or polar atom



    def GuessAtom ( self, mol, P, atTree = None, nearAtMap = None, doMsg=True, checkNewAtoms=None ) :

        # mol - molecule to add new ions/waters to
        # P - point on which to consider adding new ion/water
        # atTree - tree of atoms to consider for collisions or nearby
        # nearAtMap - only place ions/waters near these atoms
        # doMsg - create message explaining criteria used to place ion or water

        nearAts = None
        if atTree :
            nearAts = self.AtsWithinPt ( P, 6.0, atTree )
        else :
            #nearAts = [None] * len(mol.atoms)
            nearAts = []
            P = chimera.Point ( P[0], P[1], P[2] )
            for i, at in enumerate(mol.atoms) :
                V = P - at.coord()
                if V.length < 6.0 :
                    nearAts.append ( [V.length, at] )

        newAtsMap = {}
        if checkNewAtoms :
            for at in checkNewAtoms :
                V = P - at.coord()
                if V.length < 6.0 :
                    nearAts.append ( [V.length, at] )
                newAtsMap[at] = 1




        #minDistW, maxDistW = 2.5, 3.3
        #minDistI, maxDistI = 1.9, 2.5

        minDistW, maxDistW = float(self.waterMinD.get()), float(self.waterMaxD.get())
        minDistI, maxDistI = float(self.ionMinD.get()), float(self.ionMaxD.get())


        #minDistW, maxDistW = 2.4, 3.2
        #minDistI, maxDistI = 2.0, 2.4

        # these are nearby protein atoms with ion distances that are typically charged - put Mg++
        chargedAtomsIon = []

        # these are already placed ions within ion/water distances - should place water, not anotheer ion
        ionAtomsIon, ionAtomsWater = [], []

        # these are nearby protein atoms that are polar positive (e.g. N) within ion distances - put Cl-
        posPolarAtomsIon = []

        # these are nearby protein atoms that are polar negative (e.g. O) within ion distances - put Na+
        negPolarAtomsIon = []

        # these are nearby charged or polar atoms (e.g. O, N, S in Cys) within water distances - put water
        chargedAtomsWater, polarAtomsWater = [], []

        # these are non-polar, non-charged atoms that are too close; don't put anything
        collidingAtoms = []

        # is near at least atom in nearAtMap
        isNearAtMap = False

        # is near newAtsMap - put waters but not ions
        isNearNewAtMap = False

        # closest to this chain
        closestChainId, closestChainD = None, 1e9

        # iterate over nearby atoms, adding them to the above lists as appropriate
        for dist, at in nearAts :

            if at.element.name == "H" :
                continue

            #if round(dist*10.0)/10.0 < minDistI :
            if dist < minDistI :
                collidingAtoms.append ( [dist, at] )
                break

            if at.element.name == "C" or (at.element.name == "S" and at.residue.type == "MET") :
                if dist < 2.6 :
                    collidingAtoms.append ( [dist, at] )
                    break

            #if hasattr ( at, 'Q' ) and at.Q < 0.5 :
            #    continue

            #if at.altLoc != '' :
            #    continue

            # only consider points near atoms in nearAtMap, or near new atoms
            # note each atom should only be added to one list, use continue to avoid many if/else
            if nearAtMap != None and at in nearAtMap :
                isNearAtMap = True


            if at in newAtsMap :
                isNearNewAtMap = True

            else :
                # see which chain this is closest to
                if dist < closestChainD :
                    closestChainId = at.residue.id.chainId
                    closestChainD = dist

            if at.residue.type.upper() in chargedIons :
                if dist <= maxDistI :
                    ionAtomsIon.append ( [dist, at] )
                    continue
                if dist < maxDistW :
                    ionAtomsWater.append ( [dist, at] )
                    continue

            chargedAt = False
            polarAt = False
            if at.residue.type == 'HIS' and (at.name == "ND1" or at.name == "NE2") :
                chargedAt = True
            if at.residue.type == "ASP" and (at.name == "OD1" or at.name == "OD2") :
                chargedAt = True
            if at.residue.type == "GLU" and (at.name == "OE1" or at.name == "OE2") :
                chargedAt = True
            if at.residue.type == "LYS" and (at.name == "NZ") :
                chargedAt = True
            if at.residue.type == "ARG" and (at.name == "NH1" or at.name == "NH2") :
                chargedAt = True

            if chargedAt :
                if dist <= maxDistI :
                    chargedAtomsIon.append ( [dist, at] )
                    continue
                if dist <= maxDistW :
                    chargedAtomsWater.append ( [dist, at] )
                    continue

            #elif at.residue.type == "HOH" :
            #    if dist <= maxDistI :
            #        waterAtomsIon.append ( [dist, at] )
            #        continue
            #    if dist <= maxDistW :
            #        waterAtomsWater.append ( [dist, at] )
            #        continue

            else :
                # if not charged, check if polar
                if at.element.name == "N" :
                    polarAt = True
                    if dist <= maxDistI :
                        posPolarAtomsIon.append ( [dist, at] )
                        continue
                if at.element.name == "O" or (at.element.name == "S" and at.residue.type == "CYS") :
                    polarAt = True
                    if dist <= maxDistI :
                        negPolarAtomsIon.append ( [dist, at] )
                        continue

                if polarAt :
                    if dist >= minDistW and dist <= maxDistW :
                        polarAtomsWater.append ( [dist, at] )

        msg = ""
        if doMsg :
            if len(collidingAtoms) > 0 :
                msg = "Clash with:"
                for d, at in collidingAtoms :
                    msg += " " + self.At (at, d)

            else :
                msg += "Near: "
                if len(chargedAtomsIon) > 0 :
                    msg += "\n\nCharged Atoms (at Ion distance):"
                    for d, at in chargedAtomsIon :
                        msg += " \n" + self.At (at, d)

                if len(negPolarAtomsIon) > 0 :
                    msg += "\n\nNegative Polar Atoms (at Ion distance):"
                    for d, at in negPolarAtomsIon :
                        msg += " \n" + self.At (at, d)

                if len(posPolarAtomsIon) > 0 :
                    msg += "\n\nPositive Polar Atoms (at Ion distance):"
                    for d, at in posPolarAtomsIon :
                        msg += " \n" + self.At (at, d)

                if len(chargedAtomsWater) > 0 :
                    msg += "\n\nCharged Atoms (at Water distance):"
                    for d, at in chargedAtomsWater :
                        msg += " \n" + self.At (at, d)

                if len(polarAtomsWater) > 0 :
                    msg += "\n\nPolar Atoms (at Water distance):"
                    for d, at in polarAtomsWater :
                        msg += " \n" + self.At (at, d)

                if len(ionAtomsIon) > 0 :
                    msg += "\n\nIon (at Ion distance):"
                    for d, at in ionAtomsIon :
                        msg += " \n" + self.At (at, d)

                if len(ionAtomsWater) > 0 :
                    msg += "\n\nIon (at Water distance):"
                    for d, at in ionAtomsWater :
                        msg += " \n" + self.At (at, d)


        # use string set by user for type, if in list...
        ionType = "ZN"
        adds = self.addStr.get()
        if adds.upper() in chargedIons :
            ionType = adds.upper()


        atName, atRes = None, None
        clr = None
        placedType = ""

        if isNearAtMap == False and nearAtMap != None and not isNearNewAtMap :
            # a nearAtMap was given so only take points close to those atoms
            # in this case it was not close to any of them...
            pass
        elif len(collidingAtoms) == 0 :
            if len(ionAtomsIon) > 0 or ( len(ionAtomsWater) > 0 and len(negPolarAtomsIon)+len(posPolarAtomsIon) == 0 ) :
                # an ion at ion/water-distance away, likely water
                atName, atRes = "O", "HOH"
                placedType = ""
                clr = (1,0,0)
            elif len(chargedAtomsIon) > 0 :
                # charged atoms at ion distances, likely 2+ ion
                atName, atRes = ionType, ionType
                placedType = "2+ ion"
                clr = (.4,.4,.6)
            elif len(chargedAtomsWater) > 1 :
                # at least 2 charged atoms at water distances, likely 2+ ion
                atName, atRes = ionType, ionType
                placedType = "2+ ion"
                clr = (.4,.4,.6)
            elif len(negPolarAtomsIon)+len(posPolarAtomsIon) > 1 and len(negPolarAtomsIon) > 0 :
                # multiple polar atoms (at least 1 negative), likely 2+ ion
                atName, atRes = ionType, ionType
                placedType = "2+ ion"
                clr = (.4,.4,.6)
            elif len (negPolarAtomsIon) > 0 :
                # negative polar atom and no other ion around, likely 1+ ion
                atName, atRes = "NA", "NA"
                placedType = "1+ ion"
                clr = (.7,.4,.9)
            elif len (posPolarAtomsIon) > 0 :
                # positive polar atom and no other ion around, likely 1+ ion
                atName, atRes = "CL", "CL"
                placedType = "1- ion"
                clr = (0,1,0)
            elif len(polarAtomsWater) > 0 or len(chargedAtomsWater) > 0 :
                atName, atRes = "O", "HOH"
                placedType = ""
                clr = (1,0,0)

        # don't put ions if they are near new atoms but not near any other atoms
        # i.e. only put waters there...
        if isNearNewAtMap and not isNearAtMap :
            if atName != "O" :
                atName, atRes = None, None
                clr = None
                placedType = ""


        msgFull = msg
        msg = ""

        if doMsg :
            if atName != None :
                msg = "Placed %s %s/%s" % (placedType, atName, atRes)
            else :
                if len(collidingAtoms) > 0 :
                    msg = msgFull
                else :
                    msg = "Not placed - Not near any atoms (check distances in Options)"

        return msg, msgFull, atName, atRes, closestChainId, clr


    def At ( self, at, d ) :
        rt = at.residue.type
        #if rt in protein3to1 :
        #    rt = protein3to1[rt]
        #elif rt in nucleic3to1 :
        #    rt = nucleic3to1[rt]

        return " %.1fA to atom %s (element %s) in residue %s  %d, chain %s" % (d, at.name, at.element.name, rt, at.residue.id.position, at.residue.id.chainId)



    def AtsWithin (self, ats, R, atTree) :

        nearAts = []
        R2 = R * R
        for at in ats :
            pt = at.coord()
            vPt = numpy.array ( pt.data() )
            opointsNear = atTree.searchTree ( [pt[0], pt[1], pt[2]], R )
            if len(opointsNear) > 0 :
                for p in opointsNear :
                    try :
                        v = vPt - p.coord().data()
                    except :
                        continue
                    sqSum = numpy.sum ( v * v )
                    if sqSum < R2 :
                        nearAts.append (p)

        return nearAts


    def AtsWithinXf (self, ats, R, atTree) :

        nearAts = []
        R2 = R * R
        for at in ats :
            pt = at.xformCoord()
            vPt = numpy.array ( pt.data() )
            opointsNear = atTree.searchTree ( [pt[0], pt[1], pt[2]], R )
            if len(opointsNear) > 0 :
                for p in opointsNear :
                    try :
                        v = vPt - p.xformCoord().data()
                    except :
                        continue
                    sqSum = numpy.sum ( v * v )
                    if sqSum < R2 :
                        nearAts.append (p)

        return nearAts


    def AtsWithinPt (self, pt, R, atTree) :

        nearAts = []
        R2 = R * R

        vPt = numpy.array ( [pt[0], pt[1], pt[2]] )
        opointsNear = atTree.searchTree ( [pt[0], pt[1], pt[2]], R )
        if len(opointsNear) > 0 :
            for p in opointsNear :
                try :
                    v = vPt - p.coord().data()
                except :
                    continue
                sqSum = numpy.sum ( v * v )
                if sqSum < R2 :
                    nearAts.append ( [numpy.sqrt(sqSum), p] )

        return nearAts


    def AtsWithinPtXf (self, pt, R, atTree) :

        nearAts = []
        R2 = R * R

        vPt = numpy.array ( [pt[0], pt[1], pt[2]] )
        opointsNear = atTree.searchTree ( [pt[0], pt[1], pt[2]], R )
        if len(opointsNear) > 0 :
            for p in opointsNear :
                try :
                    v = vPt - p.xformCoord().data()
                except :
                    continue
                sqSum = numpy.sum ( v * v )
                if sqSum < R2 :
                    nearAts.append ( [numpy.sqrt(sqSum), p] )

        return nearAts



    def PlaceAt ( self, pt, dmap ) :

        mol = self.cur_mol
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return []

        #chainId = self.chain.get()

        #aname, chainId = self.addRess.get().split(".")
        aname = self.addStr.get()
        chainId = self.addToChain.get()
        if len(chainId) > 1 :
            umsg ( "Enter a single character in 'Add To Chain' field" )
            return



        P = chimera.Point(pt[0], pt[1], pt[2])
        P = dmap.openState.xform.inverse().apply(P)

        if self.use_mouse_max.get() :
            pts, avgMapV = PtsToMax ( [ [P[0], P[1], P[2]] ], dmap )
            maxPt = pts[0]
            V = maxPt - P
            #print " - diff to max: %.3f" % V.length
            P = maxPt

        P = dmap.openState.xform.apply(P)
        P = mol.openState.xform.inverse().apply(P)


        if self.guessOpt.get() == 'guess' :

            #print " - guessing..."

            msg, msgFull, atName, resName, closestChainId, clr = self.GuessAtom (mol, [P[0],P[1],P[2]] )

            if chainId == None or len(chainId) == 0 :
                chainId = closestChainId

            atRi = 0
            for r in mol.residues :
                if r.id.chainId == chainId and r.id.position > atRi :
                    atRi = r.id.position

            atRi += 1

            print ""
            umsg ( "Placing %s in chain %s, position %d, for map %s" % (aname, chainId, atRi, dmap.name) )


            status ( msg )
            print ""
            print msg
            print ""
            print msgFull
            print ""

            #print msg

            if atName != None :

                nres = mol.newResidue (resName, chimera.MolResId(chainId, atRi))
                nat = mol.newAtom (atName, chimera.Element(atName))

                nres.addAtom( nat )
                nat.setCoord ( P )
                #nat.drawMode = nat.Ball
                nat.color = chimera.MaterialColor( clr[0], clr[1], clr[2], 1.0 )
                nat.display = True

                nat.radius = 1.46
                if atName.lower() == "o" :
                    nat.drawMode = nat.EndCap
                else :
                    #print "ball"
                    nat.drawMode = nat.Ball

                #nat.drawMode = nat.EndCap


        else :

            if len(chainId) == 0 :
                chainId = "_"

            atRi = 0
            for r in mol.residues :
                if r.id.chainId == chainId and r.id.position > atRi :
                    atRi = r.id.position

            atRi += 1

            print ""
            umsg ( "Placing %s in chain %s, position %d, for map %s" % (aname, chainId, atRi, dmap.name) )

            nres, nat = None, None
            if aname.lower() == "w" :
                nres = mol.newResidue ("HOH", chimera.MolResId(chainId, atRi))
                nat = mol.newAtom ("O", chimera.Element('O'))
            else :
                nres = mol.newResidue (aname, chimera.MolResId(chainId, atRi))
                nat = mol.newAtom (aname, chimera.Element(aname))


            nres.addAtom( nat )
            nat.setCoord ( P )
            nat.display = True

            #nat.drawMode = nat.Ball

            if aname.lower() == "w" :
                nat.color = chimera.MaterialColor( 1.0, 0.0, 0.0, 1.0 )
                nat.drawMode = 2 # nat.EndCap
            else :
                nat.color = chimera.MaterialColor( 0.0, 1.0, 0.0, 1.0 )
                nat.drawMode = 3 # nat.Ball

            nat.drawMode = nat.EndCap


        self.RefreshTree ()


    def SetVisMap ( self ) :
        dmap = None
        mlist = chimera.openModels.list(modelTypes = [VolumeViewer.volume.Volume])
        for m in mlist :
            if m.display and not "sel_masked" in m.name :
                dmap = m
                break

        if dmap == None :
            if len(mlist) > 0 :
                dmap = mlist[0]

        if dmap != None :
            self.dmap.set ( dmap.name + " (%d)" % dmap.id )
            self.cur_dmap = dmap


    def MapMenu ( self ) :

        self.mb.menu.delete ( 0, 'end' )        # Clear menu
        from VolumeViewer import Volume
        mlist = chimera.openModels.list(modelTypes = [Volume])
        for m in mlist :
            self.mb.menu.add_radiobutton ( label=m.name + " (%d)"%m.id, variable=self.dmap,
                                           command=lambda m=m: self.MapSelected(m) )

    def SetMapMenu (self, dmap):

        mname = dmap.name if dmap else ''
        self.dmap.set(mname)
        self.cur_dmap = dmap
        #print "Set map menu to ", dmap.name

    def MapSelected ( self, dmap ) :

        self.cur_dmap = dmap
        if dmap:
            dmap.display = True
        print "Map: %s" % dmap.name



    def SetVisMol ( self ) :
        mol = None
        mlist = chimera.openModels.list(modelTypes = [chimera.Molecule])
        for m in mlist :
            if m.display :
                mol = m
                break

        if mol == None :
            if len(mlist) > 0 :
                mol = mlist[0]

        if mol != None :
            self.struc.set ( mol.name + " (%d)" % mol.id )
            self.cur_mol = mol
            SetBBAts ( mol )
            print "Mol: %s" % self.cur_mol.name

            self.RefreshTree()


    def StrucMenu ( self ) :
        self.strucMB.menu.delete ( 0, 'end' )   # Clear menu
        mlist = chimera.openModels.list(modelTypes = [chimera.Molecule])
        for m in mlist :
            self.strucMB.menu.add_radiobutton ( label=m.name+" (%d)"%m.id, variable=self.struc,
                                           command=lambda m=m: self.StrucSelected(m) )


    def StrucSelected ( self, mol ) :

        self.cur_mol = mol
        print "Selected ", mol.name, " - ", mol.id
        if mol :

            mlist = chimera.openModels.list(modelTypes = [chimera.Molecule])
            for m in mlist :
                m.display = False

            mol.display = True
            SetBBAts ( mol )

            #print "Mol: %s" % self.cur_mol.name

            self.RefreshTree()


    def select_mg_cb (self, event):


        print "Sel:", self.tree.selection()
        print "Focus:", self.tree.focus()


        to = self.tree.focus()

        if to in self.toChain :
            print " -- Chain:", self.toChain[to]
        elif to in self.toRes :
            res = self.toRes[to]
            try :
                print " -- Res: %d.%s.%s" % (res.id.position, res.type, res.id.chainId)
            except :
                pass

        elif to in self.toRess :
            ress = self.toRess[to]
            print " -- %d res" % len(ress)


        return

        if self.SelectedMgId != self.tree.focus() :
            self.SelectedMgId = self.tree.focus()
            import Segger.ar_mg_dialog
            reload ( Segger.ar_mg_dialog )
            Segger.ar_mg_dialog.show_dialog().Refresh ()


    def GetSelAtoms ( self ) :

        atoms = []
        for to in self.tree.selection () :

            if to in self.toChain :
                #print " -- Chain:", self.toChain[to]
                for res in self.cur_mol.residues :
                    if res.id.chainId == self.toChain[to] :
                        atoms.extend ( res.atoms )
            elif to in self.toRes :
                res = self.toRes[to]
                #print " -- Res: %d.%s.%s" % (res.id.position, res.type, res.id.chainId)
                atoms.extend ( res.atoms )
            elif to in self.toRess :
                ress = self.toRess[to]
                #print " -- %d res" % len(ress)
                for res in ress :
                    try :
                        atoms.extend ( res.atoms )
                    except :
                        status ( "Atoms not found, molecule may have changed" )
                        pass

        return atoms



    def SelectSel ( self ) :
        print " - selecting..."

        if self.cur_mol == None :
            umsg ( "Select a molecule first" )
            return

        ats = self.GetSelAtoms ()
        umsg ( "Selecting %d atoms" % len(ats) )

        chimera.selection.clearCurrent ()
        chimera.selection.addCurrent ( ats )



    def SelectAll ( self ) :
        print " - selecting all"

        if self.cur_mol == None :
            umsg ( "Select a molecule first" )
            return

        ats = self.GetSelAtoms ()
        umsg ( "Selecting %d atoms" % len(self.cur_mol.atoms) )

        chimera.selection.clearCurrent ()
        chimera.selection.addCurrent ( self.cur_mol.atoms )



    def ShowSel ( self ) :
        print " - showing sel..."

        if self.cur_mol == None :
            umsg ( "Select a molecule first" )
            return

        SetBBAts ( self.cur_mol )

        for m in chimera.openModels.list () :
            if type(m) == chimera.Molecule and m != self.cur_mol :
                m.display = False

        ats = self.GetSelAtoms ()
        umsg ( "Showing %d atoms" % len(self.cur_mol.atoms) )

        for at in ats :
            r = at.residue
            if r.isProt or r.isNA :
                r.ribbonDisplay = True
                for at in r.atoms :
                    at.display = False

            at.display = True
            self.ColorAt ( at )

        for bond in self.cur_mol.bonds :
            bond.display = bond.Smart



    def ShowSelOnly ( self ) :

        print " - showing sel..."

        if self.cur_mol == None :
            umsg ( "Select a molecule first" )
            return

        SetBBAts ( self.cur_mol )

        for m in chimera.openModels.list () :
            if type(m) == chimera.Molecule and m != self.cur_mol :
                m.display = False

        ats = self.GetSelAtoms ()
        umsg ( "Showing %d atoms" % len(self.cur_mol.atoms) )

        atm = {}
        resm = {}
        for at in ats :
            atm[at] = 1
            resm[at.residue] = 1

        for res in self.cur_mol.residues :

            if res in resm :
                if res.isProt or res.isNA :
                    res.ribbonDisplay = True
                    for at in res.atoms :
                        at.display = False
                else :
                    for at in res.atoms :
                        if at in atm :
                            at.display = True
                            self.ColorAt ( at )
                        else :
                            at.display = False

            else :
                if res.isProt or res.isNA :
                    res.ribbonDisplay = False
                for at in res.atoms :
                    at.display = False

        for bond in self.cur_mol.bonds :
            bond.display = bond.Smart


    def ColorAt ( self, at ) :
        if at.element.name.upper() in chargedIons :
            at.drawMode = at.Ball
            at.radius = 1.46
        else :
            at.drawMode = at.EndCap
        try :
            at.color = atomColors[at.element.name.upper()]
        except :
            at.color = atomColors[' ']

    def ShowAll ( self ) :

        print " - showing all..."

        if self.cur_mol == None :
            umsg ( "Select a molecule first" )
            return

        SetBBAts ( self.cur_mol )


        for res in self.cur_mol.residues :
            if res.isProt or res.isNA :
                res.ribbonDisplay = True
                for at in res.atoms :
                    at.display = False
            else :
                for at in res.atoms :
                    at.display = True
                    self.ColorAt ( at )


        for bond in self.cur_mol.bonds :
            bond.display = bond.Smart




    def open_mg_cb (self, event):
        #print "open"
        #print self.tree.selection()
        #print self.tree.focus()
        pass

    def close_mg_cb (self, event):
        #print "close"
        #print self.tree.selection()
        #print self.tree.focus()
        pass

    def mg_b1 (self, event):
        #print "b1"
        #print self.tree.selection()
        pass

    def mg_b1_up (self, event):
        #print "b1 up"
        #print self.tree.selection()
        pass



    def RefreshTree ( self ) :

        """ Updates treeview with water/ion atoms in self.cur_mol """

        if 0 :
            # Level 1
            at1=self.tree.insert("", 1, "", text="1" )
            self.tree.insert("", 2, "", text="2")

            # Level 2
            self.tree.insert(at1, "end", "", text="1.1", values=("t1.1"))
            self.tree.insert(at1, "end", "", text="1.2", values=("t1.2"))
            self.tree.insert(at1, "end", "", text="1.3", values=("t1.3"))

        #from os import listdir
        #from os.path import isfile, join

        self.tree.delete(*self.tree.get_children())

        #print "Refresh with: %s" % self.cur_mol.name

        SetBBAts ( self.cur_mol )

        cress = {}
        ress = []

        for res in self.cur_mol.residues :

            if res.id.chainId in cress :
                cress[res.id.chainId].append ( res )
            else :
                cress[res.id.chainId] = [ res ]

            if res.isProt or res.isNA :
                continue

            ress.append ( [res.id.position, res] )

        if 0 :
            ress.sort ( reverse=False, key=lambda x: int(x[0]) )
            for ri, res in ress :
                tRes = self.tree.insert("", "end", "", text="%d.%s - %d atoms" % (ri, res.type, len(res.atoms)) )

        else :

            self.toRes = {}
            self.toRess = {}
            self.toChain = {}

            cids = cress.keys()
            cids.sort()
            for ci in cids :
                ress = cress[ci]
                protRes, naRes, molRes = [], [], []
                for r in ress :
                    if r.isProt :
                        protRes.append ( [r.id.position, r] )
                    elif r.isNA :
                        naRes.append ( [r.id.position,  r] )
                    else :
                        molRes.append ( [r.id.position, r] )

                label = "Chain %s" % ci
                if 0 :
                    if nMol != 0 and nProt == 0 and nNA == 0 :
                        # only Ligands
                        label = "Chain %s: %d Molecules" % (ci, nMol)
                    elif nProt != 0 and nMol == 0 and nNA == 0 :
                        label = "Chain %s: %d Protein residues" % (ci, nProt)
                    elif nNA != 0 and nMol == 0 and nProt == 0 :
                        label = "Chain %s: %d Nucleotides" % (ci, nNA)
                    else :
                        label = "Chain %s: %d Residues, %d Nucleotides, %d Other" % (ci, nProt, nNA, nMol)

                chainTO = self.tree.insert("", "end", "", text=label )
                self.toChain[chainTO] = ci

                if len(molRes) > 0 :
                    molTO = self.tree.insert(chainTO, "end", "", text="%d molecules" % len(molRes) )
                    self.toRess[molTO] = [r for ri, r in molRes]

                    molRes.sort ()
                    for ri, res in molRes :
                        resTO = self.tree.insert(molTO, "end", "", text="%d.%s - %d atoms" % (ri, res.type, len(res.atoms)) )
                        self.toRes[resTO] = res

                if len(protRes) > 0 :
                    protTO = self.tree.insert(chainTO, "end", "", text="%d residues" % len(protRes) )
                    self.toRess[protTO] = [r for ri, r in protRes]

                    protRes.sort ()
                    for ri, res in protRes :
                        resTO = self.tree.insert(protTO, "end", "", text="%d.%s - %d atoms" % (ri, res.type, len(res.atoms)) )
                        self.toRes[resTO] = res

                if len(naRes) > 0 :
                    naTO = self.tree.insert(chainTO, "end", "", text="%d nucleotides" % len(naRes) )
                    self.toRess[naTO] = [r for ri, r in naRes]

                    naRes.sort ()
                    for ri, res in naRes :
                        resTO = self.tree.insert(naTO, "end", "", text="%d.%s - %d atoms" % (ri, res.type, len(res.atoms)) )
                        self.toRes[resTO] = res

                self.tree.item(chainTO, open=False)


    def Average ( self ) :

        for eid in self.tree.selection() :
            print eid
            e = self.id_mg [eid]  # {'fname':fname, 'fpath':fpath, 'vdata':vdata }
            #print e['fpath']
            path, fname = os.path.split ( e['fpath'] )

            #from ar_mg_proc import avg_frames

            avg_frames ( e['vdata'] )





    def HohE ( self ) :

        print "hoh - figure in Q-scores paper"

        mol = self.cur_mol
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return []

        chainId = self.chain.get()
        chimera.selection.clearCurrent()

        s = {184:1,280:1,278:1,183:1,236:1,281:1,357:1,282:1,399:1}
        s = {184:1,280:1,278:1,183:1,236:1,281:1,357:1}

        for res in mol.residues :
            if res.type == "HOH" or res.type == "MG" :
                if res.id.position in s :
                    for at in res.atoms :
                        at.display = True
                        if res.id.position == 183 or res.id.position == 184 :
                            at.drawMode = at.Sphere
                        chimera.selection.addCurrent ( at )
                else :
                    for at in res.atoms :
                        at.display = False



    def HohD_ ( self ) :

        print "hoh-D - distances between HOH atoms using same residue numbers"

        mols = []
        for m in chimera.openModels.list() :
            if m.display == True and type(m) == chimera.Molecule :
                mols.append ( m )

        m1, m2 = mols

        print "M1: %s" % m1.name
        print "M2: %s" % m2.name

        atm1 = {}
        for at in m1.atoms :
            if at.residue.type == "HOH" :
                aid = "%d.%s.%s" % (at.residue.id.position, at.residue.id.chainId, at.name)
                atm1[aid] = at

        ds = []
        rm, N = 0.0, 0.0
        for at2 in m2.atoms :
            if at2.residue.type == "HOH" :
                aid = "%d.%s.%s" % (at2.residue.id.position, at2.residue.id.chainId, at2.name)
                at1 = atm1[aid]

                p1 = at1.xformCoord() # m2.openState.xform.inverse().apply (at1.xformCoord())
                p2 = at2.xformCoord()
                v = p1 - p2
                ds.append ( v.length )
                rm += v.length * v.length
                N += 1.0

        rmsd = numpy.sqrt (rm/N)
        print "%.0f atoms, min: %2f, max: %.2f, avg: %.2f, rmsd: %.2f" % (N, min(ds), max(ds), numpy.average(ds), rmsd)


        ds = []
        rm, N = 0.0, 0.0
        nsame = 0
        for at2 in m2.atoms :
            if at2.residue.type == "HOH" and at2.Q >= 0.7 :
                aid = "%d.%s.%s" % (at2.residue.id.position, at2.residue.id.chainId, at2.name)
                at1 = atm1[aid]

                p1 = at1.xformCoord() # m2.openState.xform.inverse().apply (at1.xformCoord())
                p2 = at2.xformCoord()
                v = p1 - p2
                ds.append ( v.length )
                rm += v.length * v.length
                N += 1.0
                if v.length < 0.25 :
                    nsame += 1

        rmsd = numpy.sqrt (rm/N)
        print "%.0f atoms, min: %2f, max: %.2f, avg: %.2f, rmsd: %.2f -- %d same" % (N, min(ds), max(ds), numpy.average(ds), rmsd, nsame)





    def HohD ( self ) :

        print "hoh-D - distances between HOH atoms - nearest search"

        mols = []
        for m in chimera.openModels.list() :
            if m.display == True and type(m) == chimera.Molecule :
                mols.append ( m )

        m1, m2 = mols

        print "M1: %s" % m1.name
        print "M2: %s" % m2.name


        ats = [at for at in m1.atoms if at.residue.type == "HOH"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d HOH atoms / %d ats" % ( len(ats), len(m1.atoms) )
        atTree = AdaptiveTree ( points.tolist(), ats, 2.0)


        Ds = {}

        def addD ( t, d ) :
            if not t in Ds :
                Ds[t] = numpy.zeros ( 16 )
            i = int ( numpy.round(d*5.0) )
            if i < 31 :
                Ds[t][i] += 1

        num = 0
        sum, N = 0.0, 0.0
        for at2 in m2.atoms :
            if at2.residue.type == "HOH" :

                nearAts = self.AtsWithin ( [at2], 3.0, atTree )
                for nat in nearAts :
                    d = (nat.coord() - at2.coord()).length
                    addD ( "-", d )
                    if d <= 1 :
                        num += 1
                        sum += d*d
                        N += 1.0


        print ""
        print "Distances:"

        s = ""
        for i in range ( 16 ) :
            s = s + "\t%.2f" % (i/5.0)
        print s

        for t, dists in Ds.iteritems () :
            s = t
            for n in dists :
                s = s + "\t%d" % n
            print s

        print ""
        print "%d within 1A, rmsd: %.2f" % (num,numpy.sqrt(sum/N))




    def Duplicates ( self ) :

        for m in chimera.openModels.list() :
            if m.display == True and type(m) == chimera.Molecule :

                print m.name
                SetBBAts ( m )

                ats = []
                for at in m.atoms :
                    if at.residue.isProt or at.residue.isNA :
                        continue
                    ats.append ( at )

                print " - %d ats" % len(ats)
                points = _multiscale.get_atom_coordinates ( ats, transformed = True )
                print " - search tree: %d ion atoms / %d ats" % ( len(ats), len(m.atoms) )
                atTree = AdaptiveTree ( points.tolist(), ats, 2.0)

                for at in m.atoms :
                    nats = self.AtsWithinXf ( [at], 1.0, atTree )
                    for nat in nats :
                        if nat != at :
                            v = nat.xformCoord() - at.xformCoord()
                            print " - %s.%s.%d -- %s.%s.%d -- %.2f" % (at.name, at.residue.id.chainId, at.residue.id.position, nat.name, nat.residue.id.chainId, nat.residue.id.position, v.length)




    def Combine ( self ) :

        dmap = self.cur_dmap
        xf = dmap.openState.xform.inverse()

        print "hoh-D - distances between HOH atoms - nearest search -- "

        cids = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

        nmol = None

        mols = []
        mi = 0
        for m in chimera.openModels.list() :
            if m.display == True and type(m) == chimera.Molecule :
                mols.append ( m )

                atTree = None
                if nmol != None :
                    SetBBAts ( nmol )
                    ats = []
                    for at in nmol.atoms :
                        if at.residue.isProt or at.residue.isNA :
                            continue
                        ats.append ( at )

                    points = _multiscale.get_atom_coordinates ( ats, transformed = True )
                    print " - search tree: %d ion atoms / %d ats" % ( len(ats), len(nmol.atoms) )
                    atTree = AdaptiveTree ( points.tolist(), ats, 2.0)

                print "%s %d -> %s" % (m.name, mi, cids[mi])
                nmol = self.AddChain ( nmol, m, "A", cids[mi], xf )
                nmol = self.AddChain ( nmol, m, "a", cids[mi].lower(), xf, atTree )
                mi += 1




    def AddChain ( self, toMol, fromMol, cid, ncid, xf, atTree=None ) :

        if toMol == None :
            toMol = chimera.Molecule()
            toMol.name = "combine"
            chimera.openModels.add ( [toMol] )

        aMap = dict()
        from random import random as rand
        clr = ( rand(), rand(), rand() )

        for res in fromMol.residues :
            if res.id.chainId == cid :

                clash = False
                for at in res.atoms :
                    atc = at.xformCoord()
                    if atTree != None :
                        nats = self.AtsWithinPtXf ( atc, 1.0, atTree )
                        if len(nats) > 0 :
                            clash = True
                            break

                if clash :
                    print " - not adding clashing res %s.%d.%s" % (res.type, res.id.position, res.id.chainId)
                    continue

                nres = toMol.newResidue (res.type, chimera.MolResId(ncid, res.id.position))
                # print "New res: %s %d" % (nres.id.chainId, nres.id.position)
                for at in res.atoms :
                    nat = toMol.newAtom (at.name, chimera.Element(at.element.number))
                    # todo: handle alt
                    aMap[at] = nat
                    nres.addAtom( nat )
                    nat.setCoord ( xf.apply (at.xformCoord()) )
                    nat.drawMode = nat.Sphere
                    nat.color = chimera.MaterialColor( clr[0], clr[1], clr[2], 1.0 )
                    nat.display = False
                    nat.altLoc = at.altLoc
                    nat.occupancy = at.occupancy
                    nat.bfactor = at.bfactor

                nres.isHelix = res.isHelix
                nres.isHet = res.isHet
                nres.isSheet = res.isSheet
                nres.isStrand = res.isStrand
                nres.ribbonDisplay = True
                nres.ribbonDrawMode = 2
                nres.ribbonColor = chimera.MaterialColor( clr[0], clr[1], clr[2], 1.0 );

        for bond in fromMol.bonds :
            try :
                nb = toMol.newBond ( aMap[bond.atoms[0]], aMap[bond.atoms[1]] )
                nb.display = nb.Smart
            except :
                pass

        return toMol


    def IonD ( self ) :

        print "hoh-D - distances between HOH atoms - nearest search -- "

        mols = []
        for m in chimera.openModels.list() :
            if m.display == True and type(m) == chimera.Molecule :
                mols.append ( m )

        m1, m2 = mols

        print "M1: %s" % m1.name
        print "M2: %s" % m2.name

        ats = []
        for at in m1.atoms :
            if at.residue.type.upper() in chargedIons :
                ats.append ( at )

        points = _multiscale.get_atom_coordinates ( ats, transformed = True )
        print " - search tree: %d ion atoms / %d ats" % ( len(ats), len(m1.atoms) )
        atTree = AdaptiveTree ( points.tolist(), ats, 2.0)


        Ds = {}

        def addD ( t, d ) :
            if not t in Ds :
                Ds[t] = numpy.zeros ( 16 )
            i = int ( numpy.round(d*5.0) )
            if i < 31 :
                Ds[t][i] += 1


        nums, numd = {}, {}

        selAts = {}

        for at2 in m2.atoms :
            if at2.residue.type.upper() in chargedIons :

                #print at2.element.name, at2.residue.type

                nearAts = self.AtsWithinXf ( [at2], 3.0, atTree )
                for nat in nearAts :
                    d = (nat.xformCoord() - at2.xformCoord()).length

                    if nat.element.name == at2.element.name :
                        addD ( nat.element.name, d )

                    if d <= 1 :
                        selAts[at2] = 1
                        selAts[nat] = 1
                        if nat.element.name in nums :
                            nums[nat.element.name] += 1
                            numd[nat.element.name] += d
                        else :
                            nums[nat.element.name] = 1
                            numd[nat.element.name] = d


        sats = []
        for m in mols :
            SetBBAts (m)
            for at in m.atoms :
                if at.residue.isProt or at.residue.isNA :
                    continue
                if at.residue.type == "HOH" :
                    continue
                if at not in selAts :
                    sats.append ( at )
        chimera.selection.clearCurrent ()
        chimera.selection.addCurrent ( sats )


        print ""
        print "Distances:"

        s = ""
        for i in range ( 16 ) :
            s = s + "\t%.2f" % (i/5.0)
        print s

        for t, dists in Ds.iteritems () :
            s = t
            for n in dists :
                s = s + "\t%d" % n
            print s

        print ""
        for aname, num in nums.iteritems() :
            print " - %s - %d within 1A, avgd %0.2f" % (aname, num, numd[aname]/float(num) )


    def AllD ( self ) :

        print "hoh-D - distances between HOH atoms - nearest search"

        mols = []
        for m in chimera.openModels.list() :
            if m.display == True and type(m) == chimera.Molecule :
                mols.append ( m )

        m1, m2 = mols

        print "M1: %s" % m1.name
        print "M2: %s" % m2.name

        ionOrW = { "MG":2, "NA":1, "CA":2, "ZN":2, "MN":2, "FE":3, "CO":2, "NI":2, "HOH":0 }

        ats = []
        for at in m1.atoms :
            if at.residue.type.upper() in ionOrW :
                ats.append ( at )

        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d ion atoms / %d ats" % ( len(ats), len(m1.atoms) )
        atTree = AdaptiveTree ( points.tolist(), ats, 2.0)


        Ds = {}

        def addD ( t, d ) :
            if not t in Ds :
                Ds[t] = numpy.zeros ( 16 )
            i = int ( numpy.round(d*5.0) )
            if i < 31 :
                Ds[t][i] += 1

        num = 0
        for at2 in m2.atoms :
            if at2.residue.type.upper() in ionOrW :

                nearAts = self.AtsWithin ( [at2], 3.0, atTree )
                for nat in nearAts :
                    d = (nat.coord() - at2.coord()).length
                    addD ( "-", d )
                    if d <= 1 :
                        num += 1


        print ""
        print "Distances:"

        s = ""
        for i in range ( 16 ) :
            s = s + "\t%.2f" % (i/5.0)
        print s

        for t, dists in Ds.iteritems () :
            s = t
            for n in dists :
                s = s + "\t%d" % n
            print s

        print ""
        print "%d within 1A" % num



    def HohDn ( self ) :

        print "hoh-D - distances between HOH atoms - nearest search"

        m1 = None
        mols = []
        for m in chimera.openModels.list() :
            if m.display == True and type(m) == chimera.Molecule :
                mols.append ( m )

        #if len(mols) < 2 :
        #    umsg ( "Make at least two molecules visible" )
        #    return

        #m1, mols = mols[0], mols[1:]

        print "\nUsing %s as ref - %d HOH" % (m1.name, len([r for r in m1.residues if r.type == "HOH"]))

        for at in m1.atoms :
            if at.residue.type == "HOH" :
                at.nclose = 0
                at.aclose = []

        for m2 in mols :

            print "  %s - %d HOH" % (m2.name, len([r for r in m2.residues if r.type == "HOH"]))

            ats = [at for at in m2.atoms if at.residue.type == "HOH"]
            points = _multiscale.get_atom_coordinates ( ats, transformed = False )
            print "   - search tree: %d HOH atoms / %d ats" % ( len(ats), len(m2.atoms) )
            atTree = AdaptiveTree ( points.tolist(), ats, 2.0)

            for at in m2.atoms :
                if at.residue.type == "HOH" :
                    at.display = False
                    at.drawMode = 2

            for at1 in m1.atoms :
                if at1.residue.type == "HOH" :

                    at1xc = m2.openState.xform.inverse().apply(at1.xformCoord())
                    atsNear = atTree.searchTree ( at1xc.data(), 3.0 )
                    num = 0
                    for nat in atsNear :
                        d = (nat.coord() - at1xc).length
                        if d <= 1.0 :
                            num += 1
                            at1.aclose.append (nat)
                    if num > 0 :
                        at1.nclose += 1

        ns = [0] * (len(mols)+1)
        for at in m1.atoms :
            if at.residue.type == "HOH" :
                ns[at.nclose] += 1
                at.drawMode = 2
                if at.nclose == len(mols) :
                    at.display = True
                else :
                    at.display = False

        print ""
        print "Numbers:"
        print ns, " / ", numpy.sum(ns)







    def Hoh_ ( self ) :

        print ""
        print "Test solvent atoms for Q-scores (make distributions) and"
        print "distances to other atoms"
        print ""

        for m in chimera.openModels.list() :
            if type(m) == chimera.Molecule :
                SetBBAts ( m )
                num={}
                for r in m.residues :
                    if r.isProt or r.isNA :
                        continue
                    if r.type in num :
                        num[r.type] += 1
                    else :
                        num[r.type] = 1
                print m.name
                for t, n in num.iteritems() :
                    print " - ", t, n

        mol = self.cur_mol
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return []

        print " - in mol: %s" % mol.name

        chainId = self.chain.get()

        dmap = self.cur_dmap
        print " - scale map: %s" % dmap.name

        points = _multiscale.get_atom_coordinates ( mol.atoms, transformed = False )
        print " - search tree: %d ats" % ( len(mol.atoms) )
        atTree = AdaptiveTree ( points.tolist(), mol.atoms, 2.0)

        Ds = {}
        Qs = {}

        def addD ( t, d ) :
            if not t in Ds :
                Ds[t] = numpy.zeros ( 31 )
            i = int ( numpy.round(d*5.0) )
            if i < 31 :
                Ds[t][i] += 1

        def addQ ( t, q ) :
            if not t in Qs :
                Qs[t] = numpy.zeros ( 11 )
            i = int ( max (numpy.floor(q*10.0), 0) )
            if i > 10 :
                i = 10
            Qs[t][i] += 1


        for r in self.cur_mol.residues :

            if r.id.chainId != chainId :
                continue

            rid = "%d.%s" % (r.id.position, r.id.chainId)

            #if not r.isProt and not r.isNA :
            if r.type == "HOH" :

                #at = r.atoms[0]

                at = None
                for a in r.atoms :
                    if a.element.name == "O" :
                        at = a

                if at == None :
                    print " - O not found in HOH %d.%s" % (r.id.position, r.id.chainId)
                    continue

                #totAt += 1

                if 1 :
                    addQ ( 'HOH', at.Q )
                    if at.Q < 0.7 :
                        deletAts[at] = 1
                        #continue
                        pass

                nearAts = self.AtsWithin ( [at], 6.0, atTree )
                for nat in nearAts :

                    if nat == at or nat.element.name == "H" :
                        continue

                    d = (nat.coord() - at.coord()).length

                    if d < 2.0 and nat.residue.isProt :
                        print " - Hoh res %d.%s may overlap %s.%s.%d.%s - d: %.2f" % (at.residue.id.position, at.residue.id.chainId, nat.name, nat.residue.type, nat.residue.id.position, nat.residue.id.chainId, d)
                        deletAts[at] = 1
                        pass

                    #if d < 2.0 and nat.residue.id.chainId != at.residue.id.chainId :
                    #    print " - hoh res %d may overlap at %s.%s.%d.%s" % (at.residue.id, nat.name, nat.residue.type, nat.residue.id.position, nat.residue.id.chainId)

                    #if d < 2.0 :
                    #    print " - Hoh res %d may overlap at %s.%s.%d.%s - d: %.2f" % (at.residue.id.position, nat.name, nat.residue.type, nat.residue.id.position, nat.residue.id.chainId, d)

                    if d < 2.0 and nat.residue.type == "HOH" and at != nat :
                        print " - Hoh res %d.%s may overlap %s.%s.%d.%s - d: %.2f - " % (at.residue.id.position, at.residue.id.chainId, nat.name, nat.residue.type, nat.residue.id.position, nat.residue.id.chainId, d)
                        deletAts[at] = 1
                        pass

                    if nat.element.name == "O" :
                        if nat.residue.type == "HOH" :
                            addD ( "HOH-HOH", d )


        print ""
        print "Distances:"

        s = ""
        for i in range ( 31 ) :
            s = s + "\t%.2f" % (i/5.0)
        print s

        for t, dists in Ds.iteritems () :
            s = t
            for n in dists :
                s = s + "\t%d" % n
            print s





    def Stats ( self ) :

        print ""
        print "Test solvent atoms for Q-scores (make distributions) and"
        print "distances to other atoms"
        print ""

        if 0 :
            for m in chimera.openModels.list() :
                if type(m) == chimera.Molecule :
                    SetBBAts ( m )
                    num={}
                    for r in m.residues :
                        if r.isProt or r.isNA :
                            continue
                        if r.type in num :
                            num[r.type] += 1
                        else :
                            num[r.type] = 1
                    print m.name
                    for t, n in num.iteritems() :
                        print " - ", t, n

        mol = self.cur_mol
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return []

        print " - in mol: %s" % mol.name


        dmap = self.cur_dmap
        if dmap == None :
            umsg ("Select a map first")
            return []


        print " - in map: %s" % dmap.name
        minD, maxD = qscores.MinMaxD ( dmap )

        self.Log()

        umsg ( "Making statistics on ions and waters..." )

        ats = [at for at in mol.atoms if not at.element.name == "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(mol.atoms) )
        allAtTree = AdaptiveTree ( points.tolist(), ats, 2.0 )

        #points = _multiscale.get_atom_coordinates ( mol.atoms, transformed = False )
        #print " - search tree: %d ats" % ( len(mol.atoms) )
        #atTree = AdaptiveTree ( points.tolist(), mol.atoms, 2.0)

        Ds = {}
        Qs = {}
        Hoh_Hoh, Hoh_O = [], []
        Mg_Hoh, Mg_O, Mg_N = [], [], []

        def addD ( t, d ) :
            if not t in Ds :
                Ds[t] = numpy.zeros ( 21 )
            i = int ( numpy.round(d*5.0) )
            if i < 21 :
                Ds[t][i] += 1

        def addQ ( t, q ) :
            if not t in Qs :
                Qs[t] = numpy.zeros ( 11 )
            i = int ( max (numpy.floor(q*10.0), 0) )
            if i > 10 :
                i = 10
            Qs[t][i] += 1


        SetBBAts ( mol )
        doRes = []
        doAts = {}
        for res in mol.residues :
            # only looking for ions or water molecules which should have just one heavy atom
            if res.isProt or res.isNA :
                continue
            rats = [at for at in res.atoms if not at.element.name == "H"]
            if len(rats) == 1 :
                at = rats[0]
                doRes.append ( [res, at] )
                doAts[at] = 1

        if len(doRes) == 0 :
            umsg ( "No water or ions found in %s?" % mol.name )
            return

        numSame = 0
        deletAts = {}
        atI = 0
        for res, atom in doRes:

            if 1 or not hasattr ( atom, 'Q' ) :
                #at.Q = 0.0
                atom.Q = qscores.Qscore ( [atom], dmap, 0.4, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )


            rtype = "H2O" if res.type.upper() == "HOH" else res.type.upper()

            addQ ( rtype, atom.Q )

            #if at.residue.id.position == 200 and at.residue.id.chainId == "K" :
            #    print " - Q: %.3f" % atom.Q

            if atom.Q < 0.9 :
                deletAts[atom] = 1
                #continue
                pass

            nearAts = self.AtsWithin ( [atom], 6.0, allAtTree )
            for nat in nearAts :

                if nat == atom :
                    continue

                d = (nat.coord() - atom.coord()).length
                if d < 0.2 :
                    # duplicates from applying symmetry - keep just one
                    # keep atom with lowest chainId
                    at1, at2 = nat, atom
                    keepAt = at1 if at1.residue.id.chainId < at2.residue.id.chainId else at2
                    delAt = at1 if keepAt == at2 else at2
                    deletAts[delAt] = 1
                    #if not hasattr ( nat, 'keep' ) :
                    #    deletAts[nat] = 1
                    #at.keep = True
                    numSame += 1
                    continue

                #if d < 2.0 and nat.residue.isProt :
                #    print " - Hoh res %d.%s may overlap %s.%s.%d.%s - d: %.2f" % (at.residue.id.position, at.residue.id.chainId, nat.name, nat.residue.type, nat.residue.id.position, nat.residue.id.chainId, d)
                #    deletAts[at] = 1

                #if d < 2.0 and nat.residue.id.chainId != at.residue.id.chainId :
                #    print " - hoh res %d may overlap at %s.%s.%d.%s" % (at.residue.id, nat.name, nat.residue.type, nat.residue.id.position, nat.residue.id.chainId)

                #if d < 2.0 :
                #    print " - Hoh res %d may overlap at %s.%s.%d.%s - d: %.2f" % (at.residue.id.position, nat.name, nat.residue.type, nat.residue.id.position, nat.residue.id.chainId, d)

                #if d < 2.0 and nat.residue.type == "HOH" and nat.residue.id.chainId != at.residue.id.chainId :
                #    print " - Hoh res %d.%s may overlap %s.%s.%d.%s - d: %.2f - " % (at.residue.id.position, at.residue.id.chainId, nat.name, nat.residue.type, nat.residue.id.position, nat.residue.id.chainId, d)
                #    deletAts[at] = 1

                if nat.element.name == "O" :
                    if nat.residue.type == "HOH" :
                        addD ( "%s-H2O" % rtype, d )
                        addD ( "H2O-%s" % rtype, d )
                        #if d < 3.5 :
                        #    Hoh_Hoh.append ( d )
                    else :
                        addD ( "%s-O" % rtype, d )
                        #nr = nat.residue
                        #if d > 2.0 and d < 3.5 :
                        #    Hoh_O.append ( d )
                        #    #print " Hoh-O res %d.%s %s - %d.%s %s - d %.2f" % (r.id.position, r.id.chainId, r.type, nr.id.position, nr.id.chainId, nr.type, d)

                elif nat in doAts :
                    nrtype = "H2O" if nat.residue.type.upper() == "HOH" else nat.residue.type.upper()
                    addD ( "%s-%s" % (rtype, nrtype), d )
                    addD ( "%s-%s" % (nrtype, rtype), d )

                else :
                    #if nat.element.name == "N" :
                    addD ( "%s-%s" % (rtype, nat.element.name), d )
                    addD ( "%s-%s" % (nat.element.name, rtype), d )

                #if nat.element.name == "C" :
                #    addD ( "%s-C" % rtype, d )


            atI += 1
            if atI % 10 == 0 :
                status ( "Making statistics on ions and waters... at %d/%d" % (atI, len(doRes)) )
                print ".",


        status ( "Statistics done - open Log (IDLE) for results" )

        print ""
        print " - # same: %d" % numSame


        if 1 :
            print ""
            print " - deleting %d ats" % len(deletAts.keys())
            for at in deletAts.keys() :
                #print " - %s in res %s %d chain %s" % (at.name, at.residue.type, at.residue.id.position, at.residue.id.chainId)
                if len(at.residue.atoms) == 1 :
                    mol.deleteResidue ( at.residue )
                else :
                    mol.deleteAtom ( at )


        #print ""
        #print "Type\tAvg\tStd"

        #for l, ds in [ ["HOH-HOH", Hoh_Hoh], ["HOH-O", Hoh_O] ] :
        #    print "%s\t%f\t%f" % (l, numpy.average(ds), numpy.std(ds))

        #for l, ds in [ ["Mg-HOH", Mg_Hoh], ["Mg-O", Mg_O], ["Mg-N", Mg_N] ] :
        #    print "%s\t%f\t%f" % (l, numpy.average(ds), numpy.std(ds))

        #print " - tot HOH/MG: %d" % (totAt)

        print ""
        print ""
        print "Distances:"

        s = ""
        for i in range ( 21 ) :
            s = s + "\t%.2f" % (i/5.0)
        print s

        done = {}
        types = Ds.keys()
        types.sort()
        for typ in types :
            t1, t2 = typ.split("-")
            if "%s-%s" % (t2, t1) in done :
                continue
            done[typ] = 1
            f = 2 if t1 == t2 else 1 # counted twice if t1 == t2
            s = typ
            for n in Ds[typ] :
                s = s + "\t%d" % (n/f)
            print s

        print ""
        print "Q-scores:"

        s = ""
        for i in range ( 11 ) :
            s = s + "\t%.1f" % ( i/10.0 )
        print s

        for t, qs in Qs.iteritems () :
            s = t
            for n in qs :
                s = s + "\t%d" % n
            print s






    def HohShow ( self ) :

        print "hoh - show"

        mol = self.cur_mol
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return []

        chainId = self.chain.get()

        dmap = self.cur_dmap
        print " - scale map: %s" % dmap.name


        totAt, showAt = 0, 0

        tot = {}


        for r in self.cur_mol.residues :

            #if r.id.chainId != chainId :
            #    continue

            rid = "%d.%s" % (r.id.position, r.id.chainId)

            #if not r.isProt and not r.isNA :
            if not r.isProt and not r.isNA :

                for at in r.atoms :

                    totAt += 1

                    if at.Q < 0.6 :
                        at.display = False
                        try :
                            tot[at.element.name] += 1
                        except :
                            tot[at.element.name] = 1
                    else :
                        at.display = True
                        showAt += 1

        umsg ( "Showing %d/%d solvent atoms" % (showAt, totAt) )
        for tp, n in tot.iteritems() :
            print tp, n



    def GuessRes ( self ) :

        print ""

        dmap = self.cur_dmap
        if self.cur_dmap == None :
            umsg ("Select a map first")
            return []

        mol = self.cur_mol
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return []

        #chainId = self.chain.get()


        nats = int(self.mapResN.get())
        status ( "Estimating resolution of %s using %d atoms" % (dmap.name, nats) )


        ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(self.cur_mol.atoms) )
        allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)

        minD, maxD = qscores.MinMaxD ( self.cur_dmap )

        qscores.SetBBAts (self.cur_mol)
        bbAts = [at for at in self.cur_mol.atoms if at.isBB == True]
        scAts = [at for at in self.cur_mol.atoms if at.isBB == False]

        nats = nats/2

        import random
        atoms = []
        atoms = atoms + random.sample ( bbAts, min(nats,len(bbAts)) )
        atoms = atoms + random.sample ( scAts, min(nats,len(scAts)) )

        avgQ, N = 0.0, 0.0
        for ati, at in enumerate(atoms) :


            if 0 and hasattr ( at, 'Q' ) :
                avgQ += at.Q
                continue

            rr = qscores.RadCC ( [at], dmap, sigma=0.5, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
            #CC, CCm, yds, err = rr
            CC, CCm = rr

            rr = qscores.PtCC ( at.coord().data(), mol.openState.xform, dmap, sigma=0.5, allAtTree=allAtTree, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
            #CC, CCm, yds, err = rr
            CC2, CCm2 = rr

            print " - %d - %.4f,%.4f - %.4f,%.4f" % (ati, CC, CCm, CC2, CCm2)

            at.Q = CCm
            avgQ += at.Q
            N += 1.0


            status ( "Estimating resolution of %s using %d atoms - at %d" % (dmap.name, nats, ati) )

        avgQ = avgQ / N
        avgR = (avgQ-1.1244)/-0.1794

        umsg ( "Average Q=%0.2f -> res %.2f" % (avgQ, avgR) )

        self.mapRes.set ( "%.2f" % avgR )



    def Go ( self ) :

        from chimera import tasks, CancelOperation
        import traceback
        task = tasks.Task("Placing water/ions", modal = True)

        try :
            self.Go_ ( task=task )

        except Exception, err:
            umsg ( "Canceled" )
            print Exception, err
            traceback.print_exc()
            return

        finally :
            task.finished()


    def Go_ ( self, task = None ) :

        segMap = segmentation_map()
        if segMap == None :
            umsg ( "Please select a map in the Segment Map Dialog" )
            return

        smod = current_segmentation ()
        if smod == None :
            umsg ( "Please select a segmentation file in the Segment Map dialog" )
            return

        mol = self.cur_mol
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return []

        #chainId = self.chain.get()

        dmap = self.cur_dmap
        if self.cur_dmap == None :
            umsg ("Select a map first")
            return []


        #print " -- chain %s" % chainId
        toChain = self.addToChain.get()
        if len(toChain) > 1 :
            umsg ( "Enter a single character in 'Add To Chain' field" )
            return

        nearAtMap = {}
        for at in chimera.selection.currentAtoms() :
            nearAtMap[at] = 1

        if len(nearAtMap) == 0 :
            umsg ( "Select atoms near which ions/waters should be placed" )
            return []




        umsg ( "Placing water/ions in map: %s, model: %s" % (segMap.name, mol.name) )
        print ".",
        if task : task.updateStatus( "Placing water/ions in map: %s, model: %s" % (segMap.name, mol.name) )


        atTree = None
        if 0 :
            points = _multiscale.get_atom_coordinates ( mol.atoms, transformed = False )
            print " - search tree: %d ats" % ( len(mol.atoms) )
            atTree = AdaptiveTree ( points.tolist(), mol.atoms, 2.0)

        regs = list(smod.regions)

        minD, maxD = qscores.MinMaxD ( segMap )
        print " - mind: %.3f, maxd: %.3f" % (minD, maxD)


        useQ = self.useQScore.get()
        try :
            minQ = float(self.placeQ.get())
            sigQ = float(self.qsigma.get())
        except :
            umsg ( "Check Q-score and sigma, should be numbers..." )
            return


        msg = ""
        if useQ :
            print " - using min Q-score: %.2f sigma %.2f" % (minQ, sigQ)
            #msg = "minQ %.2f sigma %.2f" % (minQ, sigQ)

        umsg ( "Placing water in map: %s, model: %s, %d regions %s" % (segMap.name, mol.name, len(regs), msg) )

        import time
        startt = time.time()

        n_regs = []
        for reg in regs :
            npts = len(reg.points())
            if reg.surface_piece:
                reg.surface_piece.display = False
            if npts > 3 :
                n_regs.append ( [npts, reg] )

        # give larger regions more priority...
        n_regs.sort ( reverse=True, key=lambda x: x[0] )

        addPts = []
        addW = []
        addI = []
        skipped, skippedQ = 0, 0
        numW, numI = 0, 0
        xfI = segMap.openState.xform


        # a temporary molecule to add new waters/ions so they can be considered
        # when adding new waters/ions
        nmol = chimera.Molecule()

        regi = 0
        for numRegPts, reg in n_regs :

            if regi % 10 == 0 :
                ts = qscores.TimeLeftStr (regi, len(n_regs), time.time() - startt )
                s1 = '{:,}'.format(regi) + "/" + '{:,}'.format(len(n_regs))
                s2 = '{:,}'.format(skipped) + "/" + '{:,}'.format(skippedQ)
                s3 = '{:,}'.format(numW)
                s4 = '{:,}'.format(numI)
                status ( "At region %s (%s) - %s waters, %s ions so far, eta: %s" % (s1, s2, s3, s4, ts)  )
                if task : task.updateStatus( "At region %s (%s) - %s waters, %s ions so far, eta: %s" % (s1, s2, s3, s4, ts) )
                #print ".",

            regi += 1


            P, ctr = None, None
            if 0 :

                # take the center of all points in the region
                # - may not be close to highest value or peak
                ctr = reg.center_of_points()
                P = chimera.Point(ctr[0],ctr[1],ctr[2])

            elif 1 :
                # take the highest value in the region
                # - may not be close to the center, but it is the peak

                ctr, maxD = None, -1e9
                rpts = reg.map_points()
                map_values = segMap.interpolated_values ( rpts, segMap.openState.xform )
                #print map_values
                #break
                maxValPt = None
                for pt, val in zip(rpts, map_values) :
                    if val > maxD :
                        maxValPt, maxD = pt, val

                if 0 :
                    P = chimera.Point(maxValPt[0], maxValPt[1], maxValPt[2])
                    ctr = [maxValPt[0], maxValPt[1], maxValPt[2]]

                else :
                    #print pt
                    # go to interpolated maximum...
                    # - interestingly this can be a bit different than the voxel with
                    # - the highest value
                    pts, avgMapV = PtsToMax ( [maxValPt], segMap )
                    maxPt = pts[0]

                    #maxPt_ = [maxPt[0], maxPt[1], maxPt[2]]
                    #map_values = segMap.interpolated_values ( [maxPt_], segMap.openState.xform )
                    #print map_values
                    #print "|%.5f -> %.5f|" % (maxD, avgMapV)
                    #break

                    # if movement is too large to maximum, likely this
                    # is not a well-separated blob, so ignore it
                    V = maxPt - chimera.Point(maxValPt[0], maxValPt[1], maxValPt[2])
                    #if maxD > avgMapV :
                    #    print "|%.5f -> %.5f|" % (maxD, avgMapV),
                    #    skipped += 1
                    #    continue
                    if V.length > segMap.data.step[0] :
                        #print "|%.1f|" % V.length,
                        skipped += 1
                        continue

                    P = maxPt
                    ctr = [P[0], P[1], P[2]]


                    if useQ :
                        # check Q-score of pt
                        qs = qscores.QscorePt ( ctr, xfI, segMap, sigQ, allAtTree=None, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
                        if qs < minQ :
                            skippedQ += 1
                            continue


                #print pt
                #break


            else :
                # take the center of mass
                rpts = reg.map_points()
                #rpts = reg.points()
                map_values = segMap.interpolated_values ( rpts, segMap.openState.xform )
                #print map_values
                #break
                ctr, sum = numpy.array ( [0,0,0] ), 0.0
                for pt, val in zip(rpts, map_values) :
                    ctr += pt * val
                    sum += val
                ctr = ctr / sum
                P = chimera.Point(ctr[0],ctr[1],ctr[2])



            # check already added waters (which are not in atTree)
            # - switched to adding on the fly, so this is not needed
            if 0 :
                clash = False
                for atName, resName, clr, reg, P0 in addPts :
                    d = P - P0
                    if d.length < 2.2 :
                        clash = True
                        break

                if clash :
                    continue

            msg, msgFull, atName, resName, closestChainId, clr = self.GuessAtom ( mol, P, atTree=atTree, nearAtMap=nearAtMap, doMsg=False, checkNewAtoms=nmol.atoms )

            if atName != None :

                # add atom to new molecule, to be used in checkNewAtoms list
                nres = nmol.newResidue (resName, chimera.MolResId("A", len(nmol.residues)+1))
                nat = nmol.newAtom (atName, chimera.Element(atName))
                nres.addAtom( nat )
                nat.setCoord ( P )

                addPts.append ( [atName, resName, clr, reg, P, closestChainId] )
                if atName == 'O' :
                    numW += 1
                    addW.append ( [atName, resName, clr, reg, P, closestChainId] )
                else :
                    numI += 1
                    addI.append ( [atName, resName, clr, reg, P, closestChainId] )



        # add ions first

        #toChain = chainId.lower()
        print " - adding %d ions to %s, skipped %d/%d regions (move/Q)" % (len(addI), toChain, skipped, skippedQ)

        largestResIdForChain = {}
        for r in mol.residues :
            if not r.id.chainId in largestResIdForChain :
                largestResIdForChain[r.id.chainId] = r.id.position
            else :
                largestResIdForChain[r.id.chainId] = max(r.id.position, largestResIdForChain[r.id.chainId])


        for atName, resName, clr, reg, P, closestChainId in addI :

            cid = "_"
            if toChain == None or len(toChain) == 0 :
                cid = closestChainId
            else :
                cid = toChain

            if not cid in largestResIdForChain :
                # new chain...
                largestResIdForChain[cid] = 0

            i = largestResIdForChain[cid] + 1
            largestResIdForChain[cid] = i


            nres = mol.newResidue (resName, chimera.MolResId(cid, i))
            nat = mol.newAtom (atName, chimera.Element(atName))

            nres.addAtom( nat )
            nat.setCoord ( P )
            #nat.drawMode = nat.Ball
            nat.drawMode = 2
            nat.color = chimera.MaterialColor( clr[0], clr[1], clr[2], 1.0 )
            nat.display = True

            nat.radius = 1.46
            nat.drawMode = nat.EndCap if atName.lower() == "o" else nat.Ball
            nat.color = atomColors[atName.upper()] if atName.upper() in atomColors else atomColors[' ']


        # then add waters

        #toChain = "w"
        print " - adding %d waters to %s, skipped %d regions" % (len(addW), toChain, skipped)


        for atName, resName, clr, reg, P, closestChainId in addW :

            cid = "_"
            if toChain == None or len(toChain) == 0 :
                cid = closestChainId
            else :
                cid = toChain

            if not cid in largestResIdForChain :
                # new chain...
                largestResIdForChain[cid] = 0

            i = largestResIdForChain[cid] + 1
            largestResIdForChain[cid] = i

            nres = mol.newResidue (resName, chimera.MolResId(cid, i))
            nat = mol.newAtom (atName, chimera.Element(atName))

            nres.addAtom( nat )
            nat.setCoord ( P )
            nat.drawMode = nat.EndCap
            nat.color = chimera.MaterialColor( clr[0], clr[1], clr[2], 1.0 )
            nat.display = True

            nat.radius = 1.46
            nat.drawMode = nat.EndCap if atName.lower() == "o" else nat.Ball
            nat.color = atomColors[atName.upper()] if atName.upper() in atomColors else atomColors[' ']

            i += 1




        status ( "Added %d waters, %d ions - done" % (len(addW), len(addI)) )

        if 1 :

            qs = ""
            if useQ :
                qs = "_Q%.2f_%.2f" % (minQ, sigQ)

            thr = 0.0
            if hasattr ( segMap, "segmentThreshold" ) :
                thr = segMap.segmentThreshold
            else :
                thr = segMap.surface_levels[0]

            molPath = os.path.splitext(mol.openedAs[0])[0]
            nname = molPath + "_thr%.4f%s_[%dw]_[%di].pdb" % (thr, qs, len(addW), len(addI))

            print ""
            print "Saving pdb waters ->", nname
            chimera.PDBio().writePDBfile ( [mol], nname )

            print ""


        self.RefreshTree ()



    def Thresholds ( self ) :


        #mol = self.cur_mol
        #if self.cur_mol == None :
        #    umsg ("Select a molecule first")
        #    return []

        #chainId = self.chain.get()

        dmap = self.cur_dmap
        print " - scale map: %s" % dmap.name

        if dmap == None :
            umsg ( "Open & select a map first..." )
            return


        M = dmap.data.full_matrix()
        sdev = numpy.std(M)
        avg = numpy.average(M)
        thr = dmap.surface_levels[0]

        print ""
        umsg ( "Calculating sigma in %s..." % dmap.name )
        print "Avg: %.3f, sdev: %.3f, thr: %.4f [%.4f sdev above mean]" % (avg, sdev, thr, (thr-avg)/sdev)
        print " - 0.5 sdev above avg: %.4f" % (avg + 0.5 * sdev)
        print " - 1 sdev above avg: %.4f" % (avg + 1.0 * sdev)
        print " - 2 sdev above avg: %.4f" % (avg + 2.0 * sdev)
        print " - 3 sdev above avg: %.4f" % (avg + 3.0 * sdev)

        sig1 = avg + 1.0 * sdev
        sig2 = avg + 2.0 * sdev
        sig3 = avg + 3.0 * sdev

        umsg ( "1-sigma=[%.4f], 2-sigma:[%.4f], 3-sigma:[%.4f] -- in %s" % (sig1, sig2, sig3, dmap.name) )







def PtsToMax ( pts, dmap ) :

    from numpy import array, ones

    #import _multiscale
    #fpoints = _multiscale.get_atom_coordinates ( mol.atoms, transformed = False )
    apts = array ( pts, dtype=numpy.float32 )

    wts = ones ( len(pts), numpy.float32 )

    darray = dmap.data.matrix()

    from FitMap import locate_maximum, overlap_and_correlation

    xyz_to_ijk_tf = dmap.data.xyz_to_ijk_transform
    #map_values, outside = VolumeData.interpolate_volume_data(pts, xyz_to_ijk_tf, darray)
    #olap0, cc0, other = overlap_and_correlation ( wts, map_values )

    #print ( " ToMax -- CC %.3f" % (cc0) )

    move_tf, stats = locate_maximum(apts, wts,
                                    darray, xyz_to_ijk_tf,
                                    max_steps = 1000,
                                    ijk_step_size_min = 0.01,
                                    ijk_step_size_max = 0.5,
                                    optimize_translation = True,
                                    optimize_rotation = True,
                                    metric = 'sum product',
                                    request_stop_cb = None)

    from Matrix import chimera_xform
    xf = chimera_xform ( move_tf )
    #xT, xR = xf_2_M ( chimera_xform ( move_tf ) )
    #M = M * xT * xR

    #corr = stats['correlation']
    #print ( " -- fit CC %.3f -> %.3f, shift [%.3f], rot [%.3f], %d steps" % (cc0, corr, stats['shift'], stats['angle'], stats['steps']) )
    #print stats

    mpts = [None] * len(pts)
    for i, pt in enumerate(pts) :
        pt = chimera.Point ( apts[i][0], apts[i][1], apts[i][2] )
        xpt = xf.apply (pt)
        #mpts[i] = [xpt[0], xpt[1], xpt[2]]
        mpts[i] = xpt

    return mpts, stats['average map value']







def SetBBAts ( mol ) :

    #if hasattr ( mol, "bbats" ) :
    #    return
    #mol.bbats = True

    #print " - setting bbAts in %s" % mol.name
    for r in mol.residues :

        #r.isProt = "C" in r.atomsMap and "CA" in r.atomsMap and "N" in r.atomsMap
        #r.isProt = "CA" in r.atomsMap
        #r.isNA = "O3'" in r.atomsMap and "O5'" in r.atomsMap

        from chimera.resCode import nucleic3to1
        from chimera.resCode import protein3to1
        protein3to1['HSD'] = protein3to1['HIS']

        r.isProt = r.type in protein3to1
        r.isNA = r.type in nucleic3to1

        if r.isProt :
            r.rtype = "prot"
        elif r.isNA :
            r.rtype = "na"
        else :
            r.rtype = "?"


        if r.isNA :
            try :
                if nucleic3to1[r.type] == "G" :
                    r.baseAt = r.atomsMap["N9"][0]
                elif nucleic3to1[r.type] == "C" :
                    r.baseAt = r.atomsMap["N1"][0]
                elif nucleic3to1[r.type] == "A" :
                    r.baseAt = r.atomsMap["N9"][0]
                elif nucleic3to1[r.type] == "U" :
                    r.baseAt = r.atomsMap["N1"][0]
            except :
                #print " - baseAt not found - "
                pass


        r.bbAtoms = []
        r.scAtoms = []

        if r.isProt :
            for a in r.atoms :
                if a.element.name == "H" :
                    a.isBB, a.isSC = False, False
                    continue
                n = a.name
                a.isBB = n=="C" or n=="CA" or n=="O" or n=="N" or n=="OT1" or n=="OT2"
                a.isSC = not a.isBB
                if a.isBB :
                    r.bbAtoms.append ( a )
                else :
                    r.scAtoms.append ( a )

                a.isSugar, a.isBase = False, False

        elif r.isNA :
            for a in r.atoms :
                if a.element.name == "H" :
                    a.isBB, a.isSC = False, False
                    continue
                n = a.name

                a.isBB = n=="P" or n=="O1P" or n=="O2P" or n=="OP1" or n=="OP2" or n=="O5'" or n=="C5'" or n=="O3'"
                a.isSugar = n=="C1'" or n=="C2'" or n=="O4'" or n=="O2'" or n=="C3'" or n=="C4'"
                a.isBB = a.isBB or a.isSugar

                a.isBase = False

                if nucleic3to1[r.type] == "G" :
                    a.isBase = n=="N9" or n=="C8" or n=="N7" or n=="C5" or n=="C4" or n=="C6" or n=="O6" or n=="N1" or n=="C2" or n=="N2" or n=="N3"

                elif nucleic3to1[r.type] == "C" :
                    a.isBase = n=="N1" or n=="C2" or n=="O2" or n=="N3" or n=="C4" or n=="N4" or n=="C5" or n=="C6"

                elif nucleic3to1[r.type] == "A" :
                    a.isBase = n=="N9" or n=="C8" or n=="N7" or n=="C5" or n=="C4" or n=="N3" or n=="C2" or n=="N1" or n=="C6" or n=="N6"

                elif nucleic3to1[r.type] == "U" :
                    a.isBase = n=="N1" or n=="C2" or n=="O2" or n=="N3" or n=="C4" or n=="O4" or n=="C5" or n=="C6"

                else :
                    #print " -x- NA res %d.%s is ?" % (r.id.position, r.type)
                    break

                a.isSC = a.isBase

                #if nucleic3to1[r.type] == "G" :
                #    r.isBase = n=="" or n=="" or n=="" or n=="" or n=="" or n=="" or n=="" or n=="" or n="" or n="" or n=""
                #    r.baseAt = r.atomsMap["N9"][0]

                if a.isBB :
                    r.bbAtoms.append ( a )
                else :
                    r.scAtoms.append ( a )

        else :
            for a in r.atoms :
                a.isBB, a.isSC, a.isSugar, a.isBase = False, False, False, False





def dialog ( ) :
	from chimera import dialogs
	return dialogs.find ( SWIM_Dialog.name, create=False )


def show_dialog ( closeOld = True ):

	from chimera import dialogs

	d = dialogs.find ( SWIM_Dialog.name, create=False )
	if d :
		if closeOld :
			d.toplevel_widget.update_idletasks ()
			d.Close()
			d.toplevel_widget.update_idletasks ()
		else :
			# is there a way to bring it to front?
			return d

	dialogs.register (SWIM_Dialog.name, SWIM_Dialog, replace = True)

	d = dialogs.find ( SWIM_Dialog.name, create=True )
	d.toplevel_widget.update_idletasks ()
	d.enter()

	return d
