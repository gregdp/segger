
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
from Segger import showDevTools, timing, seggerVersion
from CGLutil.AdaptiveTree import AdaptiveTree

from chimera.resCode import nucleic3to1
from chimera.resCode import protein3to1

#devMenus = False
#showDevTools = False

import qscores
reload (qscores)

chargedIons = { "MG":2, "NA":1, "CL":-1, "CA":2, "ZN":2, "MN":2, "FE":3, "CO":2, "NI":2, "K":1 }

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

    buttons = ("Thr", "Go", "Q", "Stats", "Options", "Log")


    help = 'https://github.com/gregdp/segger/blob/master/tutorials/Segger%20Tutorial%209%20-%20SWIM.pdf'

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
            self.mb.grid (column=1, row=0, sticky='we', padx=2)
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
            self.tree.column("#0", width=50, minwidth=50, stretch=Tkinter.YES)
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

            #b = Tkinter.Button(ff, text="R", command=self.RefreshTree)
            #b.grid (column=0, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="Sel", command=self.SelectSel)
            b.grid (column=1, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="All", command=self.SelectAll)
            b.grid (column=2, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="W&I", command=self.SelW)
            b.grid (column=3, row=0, sticky='w', padx=1, pady=1)

            b = Tkinter.Button(ff, text="Show", command=self.ShowSel)
            b.grid (column=5, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="Only", command=self.ShowSelOnly)
            b.grid (column=6, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="All", command=self.ShowAll)
            b.grid (column=7, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="Z", command=self.Zone)
            b.grid (column=8, row=0, sticky='w', padx=1, pady=1)

            self.zoneRad = Tkinter.StringVar(ff)
            #self.addRess.set ( "vsgtngtkrf" )
            self.zoneRad.set ( "5" )
            e = Tkinter.Entry(ff, width=2, textvariable=self.zoneRad)
            e.grid(column=9, row=0, sticky='w', padx=1, pady=1)

            #b = Tkinter.Button(ff, text="Avg", command=self.Average)
            #b.grid (column=1, row=0, sticky='w', padx=0, pady=1)

            #b = Tkinter.Button(ff, text="Open", command=self.Open)
            #b.grid (column=2, row=0, sticky='w', padx=0, pady=1)





        if 1 :

            row += 1
            dummyFrame = Tkinter.Frame(parent, relief='groove', borderwidth=1)
            Tkinter.Frame(dummyFrame).pack()
            dummyFrame.grid(row=row,column=0,columnspan=1, pady=2, sticky='we')

            row += 1
            ff = Tkinter.Frame(parent)
            ff.grid(column=0, row=row, sticky='w')

            um = Hybrid.Checkbutton(ff, 'Place with Mouse (Ctrl+Click)', False)
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
            self.addStr.set ( "MG" )
            e = Tkinter.Entry(ff, width=5, textvariable=self.addStr)
            e.grid(column=3, row=0, sticky='w', padx=5, pady=1)

            b = Tkinter.Label(ff, text=" e.g. Mg, Ca, Na, Zn, Fe, W (Water)")
            b.grid (column=4, row=0, sticky='w', padx=0, pady=1)



            #um = Hybrid.Checkbutton(ff, 'Guess', False)
            #um.button.grid(column = 9, row=0, sticky = 'w', padx=5)
            #self.use_mouse_guess = um.variable
            #self.use_mouse_guess.set(True)
            #um.callback(self.bind_placement_button_cb)


            #b = Tkinter.Button(ff, text="SWIM", command=self.Place)
            #b.grid (column=5, row=0, sticky='w', padx=5)



        if showDevTools :

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

            b = Tkinter.Button(ff, text="Di", command=self.HohD2)
            b.grid (column=7, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Da", command=self.AllD)
            b.grid (column=8, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="wDn", command=self.HohDn)
            b.grid (column=9, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="SNi", command=self.SNi)
            b.grid (column=10, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="SN", command=self.SN)
            b.grid (column=11, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(ff, text="Mg", command=self.Mg)
            #b.grid (column=12, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="set", command=self.SetSel)
            b.grid (column=12, row=0, sticky='w', padx=5)


        if showDevTools :

            row += 1
            ff = Tkinter.Frame(parent)
            ff.grid(column=0, row=row, sticky='w')

            b = Tkinter.Button(ff, text="Comb", command=self.Combine)
            b.grid (column=7, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Flip", command=self.Flip)
            b.grid (column=8, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="vis W&I", command=self.SelVisW)
            b.grid (column=10, row=0, sticky='w', padx=1, pady=1)

            b = Tkinter.Button(ff, text="Dup", command=self.Duplicates)
            b.grid (column=11, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="RxWI", command=self.RelaxWaterIons)
            b.grid (column=12, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(ff, text="Aw", command=self.AlignWater)
            #b.grid (column=13, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(ff, text="<", command=self.AlignWaterBack)
            #b.grid (column=14, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(ff, text="RMSD", command=self.RMSD)
            #b.grid (column=10, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Ds", command=self.SwimD)
            b.grid (column=14, row=0, sticky='w', padx=5)


        if 0 and showDevTools :
            row += 1
            ff = Tkinter.Frame(parent)
            ff.grid(column=0, row=row, sticky='w')

            b = Tkinter.Label(ff, text="Map Res:")
            b.grid (column=1, row=0, sticky='w', padx=0, pady=1)

            self.mapRes = Tkinter.StringVar(ff)
            #self.addRess.set ( "vsgtngtkrf" )
            self.mapRes.set ( "" )
            e = Tkinter.Entry(ff, width=5, textvariable=self.mapRes)
            e.grid(column=2, row=0, sticky='w', padx=5, pady=1)

            b = Tkinter.Label(ff, text="A")
            b.grid (column=3, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="Est.", command=self.GuessRes)
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



            l = Tkinter.Label(ff, text='Map: ')
            l.grid(column=0, row=0, sticky='w')

            self.cur_dmap = None
            self.dmap = Tkinter.StringVar(parent)

            self.mb  = Tkinter.Menubutton ( ff, textvariable=self.dmap, relief=Tkinter.RAISED )
            self.mb.grid (column=1, row=0, sticky='we', padx=2)
            self.mb.menu  =  Tkinter.Menu ( self.mb, tearoff=0, postcommand=self.MapMenu )
            self.mb["menu"]  =  self.mb.menu


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
            dummyFrame.grid(row=0,column=orow,columnspan=1, pady=2, sticky='we')

            if 1 :
                orow += 1
                ff = Tkinter.Frame(df)
                ff.grid(column=0, row=orow, sticky='w')

                b = Tkinter.Label(ff, text="Distance ranges (in Angstroms):")
                b.grid (column=0, row=0, sticky='w', padx=0, pady=1)

                orow += 1
                ff = Tkinter.Frame(df)
                ff.grid(column=0, row=orow, sticky='w')

                b = Tkinter.Label(ff, text="   Ion distances: from ")
                b.grid (column=1, row=0, sticky='w', padx=0, pady=1)

                self.ionMinD = Tkinter.StringVar(ff)
                #self.addRess.set ( "vsgtngtkrf" )
                self.ionMinD.set ( "1.8" )
                e = Tkinter.Entry(ff, width=5, textvariable=self.ionMinD)
                e.grid(column=2, row=0, sticky='w', padx=5, pady=1)

                b = Tkinter.Label(ff, text="A to ")
                b.grid (column=3, row=0, sticky='w', padx=0, pady=1)

                self.ionMaxD = Tkinter.StringVar(ff)
                #self.addRess.set ( "vsgtngtkrf" )
                self.ionMaxD.set ( "2.5" )
                e = Tkinter.Entry(ff, width=5, textvariable=self.ionMaxD)
                e.grid(column=4, row=0, sticky='w', padx=5, pady=1)

                b = Tkinter.Label(ff, text="A")
                b.grid (column=5, row=0, sticky='w', padx=0, pady=1)

                #b = Tkinter.Label(ff, text="Distances (in Angstroms):")
                #b.grid (column=0, row=0, sticky='w', padx=0, pady=1)

                b = Tkinter.Label(ff, text="   Water distances: from ")
                b.grid (column=1, row=1, sticky='w', padx=0, pady=1)

                self.waterMinD = Tkinter.StringVar(ff)
                #self.addRess.set ( "vsgtngtkrf" )
                self.waterMinD.set ( "2.5" )
                e = Tkinter.Entry(ff, width=5, textvariable=self.waterMinD)
                e.grid(column=2, row=1, sticky='w', padx=5, pady=1)


                b = Tkinter.Label(ff, text="A to")
                b.grid (column=3, row=1, sticky='w', padx=0, pady=1)

                self.waterMaxD = Tkinter.StringVar(ff)
                #self.addRess.set ( "vsgtngtkrf" )
                self.waterMaxD.set ( "3.4" )
                e = Tkinter.Entry(ff, width=5, textvariable=self.waterMaxD)
                e.grid(column=4, row=1, sticky='w', padx=5, pady=1)

                b = Tkinter.Label(ff, text="A")
                b.grid (column=5, row=1, sticky='w', padx=0, pady=1)


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
                self.placeQ.set ( "0.7" )
                e = Tkinter.Entry(ff, width=5, textvariable=self.placeQ)
                e.grid(column=2, row=0, sticky='w', padx=5, pady=1)

                self.qsigma = Tkinter.StringVar(ff)
                self.qsigma.set ( "0.6" )

                if 1 :
                    b = Tkinter.Label(ff, text=" sigma:")
                    b.grid (column=3, row=0, sticky='w', padx=0, pady=1)

                    e = Tkinter.Entry(ff, width=5, textvariable=self.qsigma)
                    e.grid(column=4, row=0, sticky='w', padx=5, pady=1)

            if 1 :

                orow += 1
                ff = Tkinter.Frame(df)
                ff.grid(column=0, row=orow, sticky='w')

                l = Tkinter.Label(ff, text='Half Map A: ')
                l.grid(column=0, row=0, sticky='w')

                self.cur_dmap_h1 = None
                self.dmap_h1 = Tkinter.StringVar(parent)
                self.dmap_h1.set ( "" )

                self.mb_h1  = Tkinter.Menubutton ( ff, textvariable=self.dmap_h1, relief=Tkinter.RAISED )
                self.mb_h1.grid (column=1, row=0, sticky='we', padx=2)
                self.mb_h1.menu  =  Tkinter.Menu ( self.mb_h1, tearoff=0, postcommand=self.MapMenu_H1 )
                self.mb_h1["menu"]  =  self.mb_h1.menu


                orow += 1
                ff = Tkinter.Frame(df)
                ff.grid(column=0, row=orow, sticky='w')

                l = Tkinter.Label(ff, text='Half Map B: ')
                l.grid(column=0, row=0, sticky='w')

                self.cur_dmap_h2 = None
                self.dmap_h2 = Tkinter.StringVar(parent)
                self.dmap_h2.set ( "" )

                self.mb_h2  = Tkinter.Menubutton ( ff, textvariable=self.dmap_h2, relief=Tkinter.RAISED )
                self.mb_h2.grid (column=1, row=0, sticky='we', padx=2)
                self.mb_h2.menu  =  Tkinter.Menu ( self.mb_h2, tearoff=0, postcommand=self.MapMenu_H2 )
                self.mb_h2["menu"]  =  self.mb_h2.menu


            if showDevTools :
                orow += 1
                ff = Tkinter.Frame(df)
                ff.grid(column=0, row=orow, sticky='w')

                b = Tkinter.Label(ff, text="Select ")
                b.grid (column=3, row=0, sticky='w', padx=0, pady=1)

                self.selAtsType = Tkinter.StringVar(ff)
                #self.addRess.set ( "vsgtngtkrf" )
                self.selAtsType.set ( "MG,HOH" )
                e = Tkinter.Entry(ff, width=10, textvariable=self.selAtsType)
                e.grid(column=4, row=0, sticky='w', padx=5, pady=1)

                b = Tkinter.Label(ff, text="  near  ")
                b.grid (column=5, row=0, sticky='w', padx=0, pady=1)

                self.selNearAts = Tkinter.StringVar(ff)
                #self.addRess.set ( "vsgtngtkrf" )
                self.selNearAts.set ( "O,N7" )
                e = Tkinter.Entry(ff, width=10, textvariable=self.selNearAts)
                e.grid(column=6, row=0, sticky='w', padx=5, pady=1)

                b = Tkinter.Button(ff, text="Sel", command=self.SelNear)
                b.grid (column=7, row=0, sticky='w', padx=1, pady=1)




        self.optionsPanel.set(False)
        #self.selPanel.set(False)



        row += 1
        dummyFrame = Tkinter.Frame(parent, relief='groove', borderwidth=1)
        Tkinter.Frame(dummyFrame).pack()
        dummyFrame.grid(row=row,column=0,columnspan=1, pady=2, sticky='we')


        row += 1
        global msg
        msg = Tkinter.Label(parent, width = 10, anchor = 'w', justify = 'left', fg="red")
        msg.grid(column=0, row=row, sticky='ew')
        self.msg = msg
        msg.configure(text = "Press Help below for more information")


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
            button, modifiers = ('1', ['Ctrl'])
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

            # min/max distance for water (W) and ion (I) from GUI
            # by default, ion is min:1.8 max:2.5, water is min:2.5 max:3.5
            minDistW, maxDistW = float(self.waterMinD.get()), float(self.waterMaxD.get())
            minDistI, maxDistI = float(self.ionMinD.get()), float(self.ionMaxD.get())

            # use string set by user for type of ion
            # by default it's actually MG
            ionType = "MG"
            adds = self.addStr.get()
            if adds.upper() in chargedIons :
                ionType = adds.upper()

            msg, msgFull, atName, resName, closestChainId, clr = GuessAtom (mol, [P[0],P[1],P[2]], atGrid=None, nearAtMap=None, doMsg=True, minDistI=minDistI, maxDistI=maxDistI, minDistW=minDistW, maxDistW=maxDistW, ionType=ionType )

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
            nat.radius = 1.46
            if aname.lower() == "w" :
                nat.drawMode = 2 # nat.EndCap
                nat.color = chimera.MaterialColor( 1.0, 0.0, 0.0, 1.0 )
            else :
                if aname.upper() in atomColors :
                    nat.color = atomColors[aname.upper()]
                else :
                    nat.color = chimera.MaterialColor( 0.0, 1.0, 0.0, 1.0 )
                nat.drawMode = 3 # nat.Ball

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
        # if dmap: dmap.display = True
        print "Map: %s" % dmap.name



    def MapMenu_H1 ( self ) :

        self.mb_h1.menu.delete ( 0, 'end' )        # Clear menu
        from VolumeViewer import Volume
        mlist = chimera.openModels.list(modelTypes = [Volume])
        self.mb_h1.menu.add_radiobutton ( label="", variable=self.dmap_h1,
                                       command=lambda : self.MapSelected_H1(None) )
        for m in mlist :
            self.mb_h1.menu.add_radiobutton ( label=m.name + " (%d)"%m.id, variable=self.dmap_h1,
                                           command=lambda m=m: self.MapSelected_H1(m) )

    def MapSelected_H1 ( self, dmap ) :

        self.cur_dmap_h1 = dmap
        if dmap == None :
            print "Half Map 1: none"
        else :
            print "Half Map 1: %s" % dmap.name


    def MapMenu_H2 ( self ) :

        self.mb_h2.menu.delete ( 0, 'end' )        # Clear menu
        from VolumeViewer import Volume
        mlist = chimera.openModels.list(modelTypes = [Volume])
        self.mb_h2.menu.add_radiobutton ( label="", variable=self.dmap_h2,
                                       command=lambda : self.MapSelected_H2(None) )
        for m in mlist :
            self.mb_h2.menu.add_radiobutton ( label=m.name + " (%d)"%m.id, variable=self.dmap_h2,
                                           command=lambda m=m: self.MapSelected_H2(m) )

    def MapSelected_H2 ( self, dmap ) :

        self.cur_dmap_h2 = dmap
        if dmap == None :
            print "Half Map 2: none"
        else :
            print "Half Map 2: %s" % dmap.name


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


        #print "Sel:", self.tree.selection()
        #print "Focus:", self.tree.focus()


        to = self.tree.focus()

        if to in self.toChain :
            #print " -- Chain:", self.toChain[to]
            pass
        elif to in self.toRes :
            res = self.toRes[to]
            try :
                pass
                #print " -- Res: %d.%s.%s" % (res.id.position, res.type, res.id.chainId)
            except :
                pass

        elif to in self.toRess :
            ress = self.toRess[to]
            #print " -- %d res" % len(ress)


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


    def SelW ( self ) :
        print " - selecting w&i"

        if self.cur_mol == None :
            umsg ( "Select a molecule first" )
            return

        #ats = self.GetSelAtoms ()
        #umsg ( "Selecting %d atoms" % len(self.cur_mol.atoms) )

        from chimera.resCode import protein3to1
        from chimera.resCode import nucleic3to1

        chimera.selection.clearCurrent ()

        for at in self.cur_mol.atoms :
            #if at.residue.type in protein3to1 :
            #    continue
            #elif at.residue.type in nucleic3to1 :
            #    continue
            #else :
            #    chimera.selection.addCurrent ( at )

            if at.residue.type.upper() in chargedIons :
                chimera.selection.addCurrent ( at )
            elif at.residue.type.upper() == "HOH" :
                chimera.selection.addCurrent ( at )



    def SelNear ( self ) :

        if self.cur_mol == None :
            umsg ( "Select a molecule first" )
            return

        selNames = self.selAtsType.get().split(',')
        print " - selecting: ", selNames

        nearNs = self.selNearAts.get().split(',')
        if len(self.selNearAts.get()) == 0 :
            nearNs = []
        print " - near: ", nearNs

        #ats = self.GetSelAtoms ()
        #umsg ( "Selecting %d atoms" % len(self.cur_mol.atoms) )

        from chimera.resCode import protein3to1
        from chimera.resCode import nucleic3to1

        ats = [at for at in self.cur_mol.atoms if at.element.name != "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d atoms / %d ats" % ( len(ats), len(self.cur_mol.atoms) )
        atTree = AdaptiveTree ( points.tolist(), ats, 2.0)


        chimera.selection.clearCurrent ()

        for at in self.cur_mol.atoms :
            if at.residue.type.upper() in selNames :
                nearAts = self.AtsWithin ( [at], 3.0, atTree )
                for nat in nearAts :
                    if len(nearNs) == 0 or nat.name in nearNs :
                        chimera.selection.addCurrent ( at )
                        chimera.selection.addCurrent ( nat )




    def SelVisW ( self ) :
        print " - selecting w&i visible"

        if self.cur_mol == None :
            umsg ( "Select a molecule first" )
            return

        #ats = self.GetSelAtoms ()
        #umsg ( "Selecting %d atoms" % len(self.cur_mol.atoms) )

        from chimera.resCode import protein3to1
        from chimera.resCode import nucleic3to1

        chimera.selection.clearCurrent ()

        for at in self.cur_mol.atoms :
            #if at.residue.type in protein3to1 :
            #    continue
            #elif at.residue.type in nucleic3to1 :
            #    continue
            #else :
            #    chimera.selection.addCurrent ( at )

            if at.display == True :
                if 0 and at.residue.type.upper() in chargedIons :
                    chimera.selection.addCurrent ( at )
                elif at.residue.type.upper() == "HOH" :
                    chimera.selection.addCurrent ( at )



    def ShowSel ( self ) :
        print " - showing sel..."

        if self.cur_mol == None :
            umsg ( "Select a molecule first" )
            return

        #SetBBAts ( self.cur_mol )

        for m in chimera.openModels.list () :
            if type(m) == chimera.Molecule and m != self.cur_mol :
                m.display = False

        ats = self.GetSelAtoms ()
        umsg ( "Showing %d atoms" % len(self.cur_mol.atoms) )

        for at in ats :
            r = at.residue
            if r.type in protein3to1 or r.type in nucleic3to1 :
                r.ribbonDisplay = True
                for at in r.atoms :
                    at.display = False
            else :
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


    def Zone ( self ) :

        print "Zone:", self.zoneRad.get()

        try :
            rad = float ( self.zoneRad.get() )
        except :
            umsg ( "Enter a number for zone radius" )
            return

        atoms = chimera.selection.currentAtoms()
        if len(atoms) == 0 :
            umsg ( "Nothing selected" )
            return

        if self.cur_dmap == None :
            umsg ( "Select a Map" )
            return

        dmap = self.cur_dmap
        m = atoms[0].molecule

        from _multiscale import get_atom_coordinates
        points = get_atom_coordinates ( atoms, transformed = True )

        nname = os.path.splitext(dmap.name)[0] + "_Z%.0f_" % rad + ".mrc"
        cmap = PtsToMap ( points, dmap, rad, nname, clr=(.7,.7,.7,.2) )

        umsg ( "Made zone map: " + nname )
        dmap.display = False

        M = dmap.data.full_matrix()
        sdev = numpy.std(M)
        avg = numpy.average(M)

        cmap.surface_levels = [avg + 2.0 * sdev]
        chimera.runCommand ( "vol #%d style surface region all step 1" % cmap.id )


        #chimera.openModels.add ( [cmap] )

        #dpath = os.path.splitext(m.openedAs[0])[0] + "_chain_" + cid + ".mrc"
        #print " -> ", dpath
        #cmap.write_file ( dpath, "mrc" )



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
                protRes, naRes, molRes, hohRes, ionRes = [], [], [], [], []
                for r in ress :
                    if r.isProt :
                        protRes.append ( [r.id.position, r] )
                    elif r.isNA :
                        naRes.append ( [r.id.position,  r] )
                    elif r.type == "HOH" :
                        hohRes.append ( [r.id.position,  r] )
                    elif len(r.atoms) == 1 :
                        ionRes.append ( [r.id.position,  r] )
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

                if len(molRes) > 0 :
                    molTO = self.tree.insert(chainTO, "end", "", text="%d molecules" % len(molRes) )
                    self.toRess[molTO] = [r for ri, r in molRes]

                    molRes.sort ()
                    for ri, res in molRes :
                        resTO = self.tree.insert(molTO, "end", "", text="%d.%s - %d atoms" % (ri, res.type, len(res.atoms)) )
                        self.toRes[resTO] = res

                if len(ionRes) > 0 :
                    molTO = self.tree.insert(chainTO, "end", "", text="%d ions" % len(ionRes) )
                    self.toRess[molTO] = [r for ri, r in ionRes]

                    ionRes.sort ()
                    for ri, res in ionRes :
                        resTO = self.tree.insert(molTO, "end", "", text="%d.%s - %d atoms" % (ri, res.type, len(res.atoms)) )
                        self.toRes[resTO] = res

                if len(hohRes) > 0 :
                    molTO = self.tree.insert(chainTO, "end", "", text="%d water" % len(hohRes) )
                    self.toRess[molTO] = [r for ri, r in hohRes]

                    hohRes.sort ()
                    for ri, res in hohRes :
                        resTO = self.tree.insert(molTO, "end", "", text="%d.%s - %d atoms" % (ri, res.type, len(res.atoms)) )
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
        print "%d within 1A, rmsd: %.6f" % (num,numpy.sqrt(sum/N))



    def HohD2 ( self ) :

        print "hoh-D - distances between HOH atoms - nearest search"

        mols = []
        for m in chimera.openModels.list() :
            if m.display == True and type(m) == chimera.Molecule :
                mols.append ( m )

        m1, m2 = mols

        num, pp, rmsd = self.wiDistNum ("HOH", m1, m2)




    def wiAtoms ( self, m1, tp ) :

        ats1 = []
        for at in m1.atoms :
            itype = None
            rtype = at.residue.type

            if rtype.upper() in chargedIons :
                itype = "%d" % chargedIons[rtype.upper()]
            else :
                itype = rtype

            if itype.upper() == tp.upper() or rtype.upper() == tp.upper() :
                ats1.append ( at )

            elif tp == "ion" and rtype.upper() in chargedIons :
                ats1.append ( at )

        return ats1


    def wiDistNum ( self, tp, m1, m2 ) :

        print "M1: %s" % m1.name,
        print "M2: %s" % m2.name

        ats1 = self.wiAtoms ( m1, tp )
        if len(ats1) == 0 :
            print " - no atoms/res of type %s in %s" % (tp, m1.name)
            return 0.0, 0.0, 0.0

        ats2 = self.wiAtoms ( m2, tp )
        if len(ats2) == 0 :
            print " - no atoms/res of type %s in %s" % (tp, m2.name)
            return 0.0, 0.0, 0.0

        points = _multiscale.get_atom_coordinates ( ats1, transformed = True )
        #print " - search tree: %d %s atoms / %d ats" % ( len(ats1), tp, len(m1.atoms) )
        atTree = AdaptiveTree ( points.tolist(), ats1, 2.0)

        Ds = {}

        def addD ( t, d ) :
            if not t in Ds :
                Ds[t] = numpy.zeros ( 16 )
            i = int ( numpy.round(d*5.0) )
            if i < 31 :
                Ds[t][i] += 1

        num = 0
        sum, N = 0.0, 0.0
        for at2 in ats2 :
            nearAts = self.AtsWithinXf ( [at2], 1.0, atTree )
            minD = 1e9
            for nat in nearAts :
                d = (nat.xformCoord() - at2.xformCoord()).length
                addD ( "-", d )
                if d < 1.0 :
                    if d < minD :
                        minD = d
            if minD < 1e8 :
                num += 1
                sum += minD*minD
                N += 1.0

        #print ""
        #print "Distances:"

        s = ""
        for i in range ( 16 ) :
            s = s + "\t%.2f" % (i/5.0)
        #print s

        for t, dists in Ds.iteritems () :
            s = t
            for n in dists :
                s = s + "\t%d" % n
            #print s

        #pp = 100.0 * float(num) / float ( min(len(ats1),len(ats2)) )
        pp = 100.0 * float(num) / float ( len(ats1) )
        rmsd = numpy.sqrt(sum/N) if N > 0 else 0.0

        #print ""
        print "%d/%d|%d %.0f%% within 1A, rmsd: %.6f" % (num, len(ats1), len(ats2), pp, rmsd)

        return num, pp, rmsd






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




    def RelaxWaterIons ( self ) :

        dmap, mol = self.cur_dmap, self.cur_mol
        print " - map: %s" % dmap.name

        M = dmap.data.full_matrix()
        sdev = numpy.std(M)
        avg = numpy.average(M)
        thr = dmap.surface_levels[0]
        thr3 = avg + 3.0 * sdev
        print "    - 3 sdev above avg: %.4f" % thr3


        print " - model: %s" % mol.name

        wiRes = [r for r in mol.residues if (r.type == "HOH" or r.type in chargedIons)]
        print " - %d water/ion" % len(wiRes)

        ats = [at for at in mol.atoms if not at.element.name == "H"]

        import gridm; reload(gridm)
        g1 = gridm.Grid ()
        g1.FromAtomsLocal ( ats, 5.0 )


        from time import time

        for res in wiRes :
            at = res.atoms[0]
            if not hasattr ( at, 'coord0' ) :
                at.coord0 = at.coord()

        for i in range ( 10 ) :

            startt = time()
            totD, totN = 0.0, 0
            totD2, totN2 = 0.0, 0

            for res in wiRes :
                at = res.atoms[0]

                map_vals = dmap.interpolated_values ( [at.coord().data()], dmap.openState.xform )

                nearAts = g1.AtsNearPtLocal ( at.coord() )
                mv = chimera.Vector(0,0,0)

                if res.type == "HOH" : # and res.id.position == 491 :
                    #print "%d - %s - %.3f - %d near" % (res.id.position, res.type, map_vals[0], len(nearAts))
                    for nat, v in nearAts :
                        if nat == at :
                            continue
                        vn = v / -v.length
                        if nat.element.name == "C" :
                            if v.length < 3.4 :
                                #print " - C %.2f" % v.length
                                mv += vn * 0.1
                        elif nat.residue.type == "HOH" :
                            if abs(nat.occupancy - 1.0) < 0.01 :
                                if v.length < 2.6 :
                                    mv += vn * 0.01
                        elif nat.residue.type in chargedIons :
                            if v.length < 2.3 :
                                mv += vn * 0.01
                        elif nat.element.name == "N" or nat.element.name == "O" :
                            #hasH = False
                            #for natn in nat.neighbors :
                            #    if natn.element.name == "H" :
                            #        hasH = True
                            #        break
                            #if hasH :
                            #    pass
                            #else :
                            if v.length < 2.6 :
                                mv += vn * 0.01

                if res.type == "MG" : # and res.id.position == 491 :
                    for nat, v in nearAts :
                        if nat == at :
                            continue
                        vn = v / -v.length
                        if nat.element.name == "C" :
                            if v.length < 3.4 :
                                #print " - C %.2f" % v.length
                                mv += vn * 0.1
                        elif nat.residue.type == "MG" :
                            if v.length < 3.4 :
                                print " MG %d -- MG %d %.2f" % (res.id.position, nat.residue.id.position, v.length)
                                mv += vn * 0.01
                                #return
                        elif nat.element.name == "N" or nat.element.name == "O" :
                            if v.length < 2.2 :
                                mv += vn * 0.01

                #print mv
                if mv.length > 1e-6 :
                    npos = at.coord() + mv
                    vals = dmap.interpolated_values ( [npos.data()], dmap.openState.xform )
                    if vals[0] > thr3 :
                        g1.MoveAtomLocal ( at, npos )
                        if res.type == "HOH" :
                            totD += mv.length
                            totN += 1
                        else :
                            totD2 += mv.length
                            totN2 += 1

            print " - moved %d/%.2f -- %d/%.2f -- %.2f" % (totN, totD, totN2, totD2, time()-startt)

        sumMv, maxMv, maxMoveAt = 0.0, 0.0, None
        for res in wiRes :
            at = res.atoms[0]
            v = at.coord() - at.coord0
            sumMv += v.length
            if v.length > maxMv :
                maxMv = v.length
                maxMoveAt = at

        sumMv /= float(len(wiRes))
        print " - avg move %.3f, max %.3f by %d" % (sumMv, maxMv, maxMoveAt.residue.id.position)


    def AlignWater ( self ) :

        res = chimera.selection.currentResidues()[0]
        print res.type, res.id.position, len(res.atoms), len(res.okNearAtoms)

        atO = res.atomsMap["O"][0]
        atH1 = res.atomsMap["H1"][0]
        atH2 = res.atomsMap["H2"][0]
        v1 = atH1.coord() - atO.coord()
        v2 = atH2.coord() - atO.coord()
        mp = (atH1.coord().toVector() + atH2.coord().toVector())/2.0
        mp = chimera.Point ( mp.x, mp.y, mp.z )
        ax1 = mp - atO.coord(); ax1.normalize()
        ax2 = chimera.cross ( v1, v2 ); ax2.normalize()
        ax3 = chimera.cross ( ax1, ax2 ); ax3.normalize()

        xf = chimera.Xform.translation ( atO.coord().toVector() * -1.0 )
        xf.premultiply ( chimera.Xform.rotation (ax1, 90.0) )
        xf.premultiply ( chimera.Xform.rotation(ax2, 180.0) )
        xf.premultiply ( chimera.Xform.translation ( atO.coord().toVector() ) )

        e1pos = xf.apply ( atH1.coord() )
        e2pos = xf.apply ( atH2.coord() )

        nmol = chimera.Molecule()
        nmol.name = "electrons"
        nres = nmol.newResidue ("E", chimera.MolResId('E', 1))
        e1 = nmol.newAtom ("H", chimera.Element('H')); nres.addAtom( e1 ); e1.setCoord ( e1pos )
        e2 = nmol.newAtom ("H", chimera.Element('H')); nres.addAtom( e2 ); e2.setCoord ( e2pos )
        #eO = nmol.newAtom ("O", chimera.Element('O')); nres.addAtom( eO ); e2.setCoord ( e2pos )
        chimera.openModels.add ( [nmol] )
        e1.color = chimera.MaterialColor( .3, .3, .3, .5 )
        e2.color = chimera.MaterialColor( .3, .3, .3, .5 )

        atO.coord0 = atO.coord()
        atH1.coord0 = atH1.coord()
        atH2.coord0 = atH2.coord()

        ptsN, ptsW = [], []
        ptsN.append ( atO.coord().data() )
        ptsW.append ( atO.coord().data() )

        res.okNearAtoms.sort ( reverse=False, key=lambda x: x[0] )
        atHi = 0; hatoms = [atH1, atH2]
        atEi = 0; eatoms = [e1, e2]
        for d, nat in res.okNearAtoms :
            print "%s.%d - %.2f" % (nat.name, nat.residue.id.position, d)
            if nat.element.name == "N" :
                hAts = [a for a in nat.bondsMap.keys() if a.element.name == "H"]
                if len(hAts) == 0 :
                    # acceptor, align H
                    v = nat.coord() - atO.coord(); v.normalize()
                    ptsN.append ( atO.coord() + v )
                    ptsW.append ( hatoms[atHi].coord() )
                    atHi += 1
                elif len(hAts) == 1 :
                    print " - N has 1 H?"
                elif len(hAts) == 2 :
                    # donor, align electron
                    v = nat.coord() - atO.coord(); v.normalize()
                    ptsN.append ( atO.coord() + v )
                    ptsW.append ( eatoms[atEi].coord() )
                    atEi += 1

        print "%d match pos" % len(ptsN)
        ptsN = numpy.array ( ptsN )
        ptsW = numpy.array ( ptsW )

        #xf, rmsd = chimera.match.matchPositions(fixedPts, movePts)
        xf, rmsd = chimera.match.matchPositions(ptsN, ptsW)
        print rmsd
        for at in res.atoms + [e1, e2] :
            at.setCoord ( xf.apply ( at.coord() ) )


    def AlignWaterBack ( self ) :

        res = chimera.selection.currentResidues()[0]
        for at in res.atoms :
            at.setCoord ( at.coord0 )



    def Combine_ ( self ) :

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
                #nmol = self.AddChain ( nmol, m, "A", cids[mi], xf )
                #nmol = self.AddChain ( nmol, m, "A", cids[mi].lower(), xf, atTree )
                nmol = self.AddChain ( nmol, m, "A", cids[mi], xf, atTree )
                mi += 1


    def Combine2 ( self ) :

        dmap = self.cur_dmap
        xf = dmap.openState.xform.inverse()

        ms = []
        for m in chimera.openModels.list() :
            if type(m) == chimera.Molecule :
                ms.append ( m )
                for r in m.residues :
                    r.ribbonDisplay = False
                    for at in r.atoms :
                        at.display = False

        print ms[0].name
        print ms[1].name
        print ms[2].name

        ats1 = [at for at in ms[1].atoms if not at.element.name == "H"]
        ats2 = [at for at in ms[2].atoms if not at.element.name == "H"]

        import gridm; reload(gridm)
        g1 = gridm.Grid (); g1.FromAtomsLocal ( ats1, 1.0 )
        g2 = gridm.Grid (); g2.FromAtomsLocal ( ats2, 1.0 )

        for at in ats1 + ats2 + ms[0].atoms :
            at.display = False

        dAts = []
        mAts = []
        for at in ms[0].atoms :

            if at.residue.type.upper() == "HOH" or at.residue.type.upper() in chargedIons :

                found1, found2 = None, None
                found1_, found2_ = None, None

                for nat, v in g1.AtsNearPtLocal ( at.coord() ) :
                    if nat.residue.type.upper() == at.residue.type.upper() :
                        found1 = nat
                    else :
                        found1_ = nat

                for nat, v in g2.AtsNearPtLocal ( at.coord() ) :
                    if nat.residue.type.upper() == at.residue.type.upper() :
                        found2 = nat
                    else :
                        found2_ = nat

                if found1 != None or found2 != None :
                    take = True
                else :
                    dAts.append ( at )

                if found1_ != None or found2_ != None :
                    take = True
                    at.display = True
                    if found1_ : found1_.display = True
                    if found2_ : found2_.display = True
                    if found1 : found1.display = True
                    if found2 : found2.display = True
                else :
                    mAts.append ( at )

        chimera.selection.clearCurrent ()
        chimera.selection.addCurrent ( mAts )

        # https://www.youtube.com/watch?v=nB5GDITfvZw


    def Combine ( self ) :

        dmap = self.cur_dmap
        xf = dmap.openState.xform.inverse()

        ms = []
        for m in chimera.openModels.list() :
            if type(m) == chimera.Molecule and m.display == True :
                ms.append ( m )
                for r in m.residues :
                    r.ribbonDisplay = False
                    for at in r.atoms :
                        at.display = False

        print "m1:", ms[0].name
        print "m2:", ms[1].name

        M2, M1 = ms

        ats1 = [at for at in M1.atoms if not at.element.name == "H"]

        for at in M2.atoms + M1.atoms :
            if at.residue.type.upper() == "HOH" or at.residue.type.upper() in chargedIons :
                at.display = False

        import gridm; reload(gridm)
        g1 = gridm.Grid ()
        g1.FromAtoms ( ats1, 1.0 )

        delAts = []
        mAts, wAts, iAts = [], [], []
        sumDiff, sumN, maxD, minD = 0.0, 0.0, 0.0, 1e9
        for at in M2.atoms :

            if at.residue.type.upper() == "HOH" or at.residue.type.upper() in chargedIons :

                found1 = None
                found1_ = None

                for nat, v in g1.AtsNearPt ( at.xformCoord() ) :
                    if nat.residue.type.upper() == at.residue.type.upper() :
                        if v.length <= 1.0 :
                            found1 = nat
                        #print ".",
                    else :
                        if nat.residue.type.upper() == "HOH" or nat.residue.type.upper() in chargedIons :
                            found1_ = nat

                if found1 != None :
                    take = True
                    if 0 :
                        at.display = True
                        found1.display = True
                    #if found1_ : found1_.display = True
                    #if found1 : found1.display = True
                    v = found1.xformCoord() - at.xformCoord()
                    sumDiff += v.length
                    sumN += 1.0
                    maxD = max (maxD, v.length)
                    minD = min ( minD, v.length)
                    if at.residue.type.upper() == "HOH" :
                        wAts.append ( found1 )
                    else :
                        iAts.append ( found1 )
                else :
                    delAts.append ( at )
                    if found1_ != None :
                        #nat = found1_
                        if 1 :
                            at.display = True
                            found1_.display = True
                        mAts.append ( at )
                        print "%s.%d(%d) - %.2f - %s.%d(%d)" % (at.residue.type, at.residue.id.position, at.molecule.id, v.length, found1_.residue.type, found1_.residue.id.position, found1_.molecule.id)

        chimera.selection.clearCurrent ()
        chimera.selection.addCurrent ( mAts )

        print " - water/ion move avg %.2f, min %.2f, max %.2f, --%d hoh, %d ion, %d mix--" % (sumDiff / sumN, minD, maxD, len(wAts), len(iAts), len(mAts))


    def Flip ( self ) :

        dmap = self.cur_dmap
        xf = dmap.openState.xform.inverse()

        ms = []
        for m in chimera.openModels.list() :
            if type(m) == chimera.Molecule and m.display == True :
                ms.append ( m )
                for r in m.residues :
                    r.ribbonDisplay = False
                    for at in r.atoms :
                        at.display = False

        print "m1:", ms[0].name
        print "m2:", ms[1].name

        M2, M1 = ms

        ats1 = [at for at in M1.atoms if not at.element.name == "H"]

        for at in M2.atoms + M1.atoms :
            if at.residue.type.upper() == "HOH" or at.residue.type.upper() in chargedIons :
                at.display = False

        import gridm; reload(gridm)
        g1 = gridm.Grid ()
        g1.FromAtoms ( ats1, 1.0 )

        delAts = []
        mAts1, mAts2, wAts, iAts = [], [], [], []
        sumDiff, sumN, maxD, minD = 0.0, 0.0, 0.0, 1e9
        for at in M2.atoms :

            if at.residue.type.upper() == "HOH" or at.residue.type.upper() in chargedIons :

                found1 = None
                found1_ = None

                for nat, v in g1.AtsNearPt ( at.xformCoord() ) :
                    if nat.residue.type.upper() == at.residue.type.upper() :
                        if v.length <= 1.0 :
                            found1 = nat
                        #print ".",
                    else :
                        if nat.residue.type.upper() == "HOH" or nat.residue.type.upper() in chargedIons :
                            found1_ = nat

                if found1 != None :
                    take = True
                    at.display = True
                    found1.display = True
                    #if found1_ : found1_.display = True
                    #if found1 : found1.display = True
                    v = found1.xformCoord() - at.xformCoord()
                    sumDiff += v.length
                    sumN += 1.0
                    maxD = max (maxD, v.length)
                    minD = min ( minD, v.length)
                    if at.residue.type.upper() == "HOH" :
                        wAts.append ( found1 )
                    else :
                        iAts.append ( found1 )
                else :
                    delAts.append ( at )
                    #at.display = True
                    if found1_ != None :
                        nat = found1_
                        mAts1.append ( at )
                        mAts2.append ( found1_ )
                        print "%s.%d(%d) - %.2f - %s.%d(%d)" % (at.residue.type, at.residue.id.position, at.molecule.id, v.length, found1_.residue.type, found1_.residue.id.position, found1_.molecule.id)

        chimera.selection.clearCurrent ()
        chimera.selection.addCurrent ( mAts1 )

        print " - water/ion move avg %.2f, min %.2f, max %.2f, --%d hoh, %d ion, %d mixed--" % (sumDiff / sumN, minD, maxD, len(wAts), len(iAts), len(mAts1))

        print "from %s.%d" % (M2.name, M2.id)
        for at in mAts1 :
            if at.residue.type == "MG" :
                print "%s.%d(%d) -> H2O" % (at.residue.type, at.residue.id.position, at.molecule.id)
                self.FlipToWater ( at )

        print "from %s.%d" % (M1.name, M1.id)
        for at in mAts2 :
            if at.residue.type == "MG" :
                print "%s.%d(%d) -> H2O" % (at.residue.type, at.residue.id.position, at.molecule.id)
                self.FlipToWater ( at )


    def FlipToWater ( self, at ) :
        mol = at.molecule
        C, rid, cid, r = at.coord(), at.residue.id.position, at.residue.id.chainId, at.radius
        mol.deleteResidue ( at.residue )
        nres = mol.newResidue ( "HOH", chimera.MolResId(cid, rid))
        nat = mol.newAtom ( "O", chimera.Element("O") )
        nres.addAtom( nat )
        nat.setCoord ( C )
        nat.drawMode = nat.EndCap
        nat.radius = r
        nat.color = atomColors["O"] if "O" in atomColors else atomColors[' ']
        nat.display = True



    # try to match distances between waters/ions to distances between nearest NT atoms
    def Combine_local_RMS ( self ) :

        print "__"

        dmap = self.cur_dmap
        xf = dmap.openState.xform.inverse()

        ms = []
        rmsMol = None
        for m in chimera.openModels.list() :
            if type(m) == chimera.Molecule and m.display == True :
                if "rmsf" in m.name :
                    print "rmsMol: %s" % m.name
                    rmsMol = m
                else :
                    ms.append ( m )

        print "m1:", ms[0].name
        print "m2:", ms[1].name

        M1, M2 = ms

        ats1 = [at for at in M1.atoms if not at.element.name == "H"]

        for at in M2.atoms + M1.atoms :
            if at.residue.type.upper() == "HOH" or at.residue.type.upper() in chargedIons :
                at.display = True

        import gridm; reload(gridm)
        g1 = gridm.Grid ()
        g1.FromAtoms ( ats1, 4.0 )

        from gridm import Grid
        m2Grid = Grid()
        m2Ats = []
        for at in M2.atoms :
            if not at.element.name == "H" :
                if not at.residue.type.upper() == "HOH" :
                    if not at.residue.type.upper() in chargedIons :
                        m2Ats.append ( at )
        m2Grid.FromAtoms ( m2Ats, 4.0 )

        m1AtMap = {}
        for at in M1.atoms :
            if not at.element.name == "H" :
                if not at.residue.type.upper() == "HOH" :
                    if not at.residue.type.upper() in chargedIons :
                        aid = "%s%d%s%s" % (at.residue.id.chainId, at.residue.id.position, at.name, at.altLoc)
                        m1AtMap[aid] = at

        m2AtMap = {}
        for at in M2.atoms :
            if not at.element.name == "H" :
                if not at.residue.type.upper() == "HOH" :
                    if not at.residue.type.upper() in chargedIons :
                        aid = "%s%d%s%s" % (at.residue.id.chainId, at.residue.id.position, at.name, at.altLoc)
                        m2AtMap[aid] = at

        rmsAtMap = {}
        for at in rmsMol.atoms :
            if not at.element.name == "H" :
                if not at.residue.type.upper() == "HOH" :
                    if not at.residue.type.upper() in chargedIons :
                        aid = "%s%d%s%s" % (at.residue.id.chainId, at.residue.id.position, at.name, at.altLoc)
                        rmsAtMap[aid] = at

        #m1Grid = Grid()
        #m1Grid.FromAtoms ( [at for at in M1.atoms if (not at.element.name == "H"] )

        delAts = []
        mAts = []
        sumDiff, sumN, maxD, minD = 0.0, 0.0, 0.0, 1e9

        fList, nfList = [], []
        fList_, nfList_ = [], []

        for at in M2.atoms :

            if at.residue.type.upper() == "HOH" or at.residue.type.upper() in chargedIons :

                found1 = None
                found1_ = None

                for nat, v in g1.AtsNearPt ( at.xformCoord(), 1.0 ) :
                    if nat.residue.type.upper() == at.residue.type.upper() :
                        if v.length <= 1.0 :
                            found1 = nat
                        #print ".",
                    #else :
                    #    found1_ = nat

                m1At = None
                nearAts = m2Grid.AtsNearPt ( at.xformCoord() )
                nearestAt = None
                if len(nearAts) > 0 :
                    nearAts.sort ( reverse=False, key=lambda x: x[1].length )
                    nat, nearestD = nearAts[0]
                    naid = "%s%d%s%s" % (nat.residue.id.chainId, nat.residue.id.position, nat.name, nat.altLoc)
                    #m1At = rmsAtMap [naid]
                    m1At = m1AtMap [naid]
                    nearestAt = nat

                if found1 != None :
                    take = True
                    at.display = False
                    found1.display = False
                    #if found1_ : found1_.display = True
                    #if found1 : found1.display = True
                    v = found1.xformCoord() - at.xformCoord()
                    sumDiff += v.length
                    sumN += 1.0
                    maxD = max (maxD, v.length)
                    minD = min ( minD, v.length)

                    #vRes = m1At.xformCoord() - nat.xformCoord()
                    #print "%f\t%f" % (v.length, vRes.length)
                    #print "%f\t%f" % (v.length, m1At.bfactor)
                    if nearestAt and m1At :
                        #fList.append ( m1At.bfactor )
                        fList.append ( (nearestAt.xformCoord() - m1At.xformCoord()).length )

                else :
                    delAts.append ( at )

                    if nearestAt and m1At :
                        #nfList.append ( m1At.bfactor )
                        nfList.append ( (nearestAt.xformCoord() - m1At.xformCoord()).length )
                    #at.display = True

                #if found1_ != None :
                #    take = True
                #else :
                #    mAts.append ( at )

        chimera.selection.clearCurrent ()
        chimera.selection.addCurrent ( delAts )

        print " - water/ion move avg %.2f, min %.2f, max %.2f" % (sumDiff / sumN, minD, maxD)

        print "match\t%f\t%f" % (numpy.mean(fList), numpy.std(fList))
        print "no match\t%f\t%f" % (numpy.mean(nfList), numpy.std(nfList))


        fList, nfList = [], []
        fList_, nfList_ = [], []

        C1 = numpy.array ( [.33,.56,.88] )
        #C2 = numpy.array ( [.92,.20,.15] )
        C2 = numpy.array ( [.99,.99,.3] )

        m1Res, m2Res = None, None

        for res in rmsMol.residues :
            hasWI = False
            sumB, numB = 0.0, 0.0
            sumD, numD = 0.0, 0.0

            if res.type == "MG" :
                continue

            for at in res.atoms :
                nats = g1.AtsNearPt ( at.xformCoord(), 3.4 )
                for nat, v in nats :
                    if nat.residue.type.upper() == "HOH" or nat.residue.type.upper() in chargedIons :
                        hasWI = True
                sumB += at.bfactor
                numB += 1.0

                aid = "%s%d%s%s" % (at.residue.id.chainId, at.residue.id.position, at.name, at.altLoc)
                m1at = m1AtMap[aid]
                m2at = m2AtMap[aid]
                sumD += (m1at.xformCoord() - m2at.xformCoord()).length

                m1Res = m1at.residue
                m2Res = m2at.residue

            avgB = sumB/numB
            f = min(10.0, avgB) / 10.0
            C = (1.0-f)*C1 + f * C2
            m1Res.ribbonColor = chimera.MaterialColor ( C[0], C[1], C[2], 1.0  )
            m2Res.ribbonColor = chimera.MaterialColor ( C[0], C[1], C[2], 1.0  )
            res.ribbonColor = chimera.MaterialColor ( C[0], C[1], C[2], 1.0  )

            if hasWI :
                #fList.append ( sumB/numB )
                fList.append ( sumD/numB )
            else :
                #nfList.append ( sumB/numB )
                nfList.append ( sumD/numB )

            #print "%f\t%f" % (sumD/numB, sumB/numB)

        print "wi\t%d\t%f\t%f" % (len(fList), numpy.mean(fList), numpy.std(fList))
        print "no\t%d\t%f\t%f" % (len(nfList), numpy.mean(nfList), numpy.std(nfList))



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

        print "hoh-D - distances between ion atoms - nearest search -- "

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

        if len(mols) < 2 :
            umsg ( "Make at least two molecules visible" )
            return

        m1, mols = mols[0], mols[1:]

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







    def Stats0 ( self ) :

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

        #chainId = self.chain.get()

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

            #if r.id.chainId != chainId :
            #    continue

            #rid = "%d.%s" % (r.id.position, r.id.chainId)

            #if not r.isProt and not r.isNA :
            if r.type == "HOH" :

                #at = r.atoms[0]

                at = None
                for a in r.atoms :
                    if a.element.name == "O" :
                        at = a

                if at == None :
                    #print " - O not found in HOH %d.%s" % (r.id.position, r.id.chainId)
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



    def Mg ( self ) :

        mol = self.cur_mol
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return []

        print " - in mol: %s" % mol.name

        ats = [at for at in mol.atoms if not at.element.name == "H"]

        import grid; reload(grid)
        agrid = grid.Grid ()
        agrid.FromAtomsLocal ( ats, 2.8 )

        mgByNumHoh = []
        for at in ats :

            if at.residue.type == "MG" :
                nats = agrid.AtsNearPtLocal ( at.coord() )
                numHoh = 0
                hohAts = []
                nAts = {}
                for nat, v in nats :
                    if nat.residue.type == "HOH" :
                        numHoh += 1
                        hohAts.append ( nat )
                    else :
                        atType = nat.element.name
                        if atType in nAts :
                            nAts[atType].append ( nat )
                        else :
                            nAts[atType] = [nat]

                #if not "N" in nAts : continue
                if not "O" in nAts : continue

                mgByNumHoh.append ( [numHoh, at, hohAts, nAts] )

        mgByNumHoh.sort ( reverse=True, key=lambda x: x[0] )
        for numHoh, at, hohAts, nAts in mgByNumHoh :
            if numHoh > 0 :
                print "%s - %d.%s - %d" % (at.residue.type, at.residue.id.position, at.residue.id.chainId, numHoh ),
                if "O" in nAts :
                    nat = nAts["O"][0]
                    print " -- %s.%d.%s" % (nat.residue.type, nat.residue.id.position, nat.residue.id.chainId),
                print ""



    def SetSel ( self ) :

        toType = self.addStr.get()
        print toType

        for at in chimera.selection.currentAtoms() :
            print " %d %s -> %s" % ( at.residue.id.position, at.residue.type, toType )



    def SwimD ( self ) :

        D = {}
        D_ = {}
        for m in chimera.openModels.list() :
            if type(m) == chimera.Molecule :
                if m.display == True :
                    self.AddDs ( m, D, D_ )

        print ""
        print ""
        print "Near atoms, using full atom name"
        print ""
        print ""
        print "Ion/H2O\tAtom\tMean\tStDev\tMin\tMax\t# counted"
        for rtype, nnames in D.iteritems () :
            rt = "H2O" if rtype == "HOH" else rtype
            rt = self.IonLabel ( rt )
            print ""
            print "%s" % rt,
            for nname, ds in nnames.iteritems() :
                mean, stdev, min, max = numpy.mean (ds), numpy.std (ds), numpy.min (ds), numpy.max(ds)
                print "\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%d" % (nname, mean, stdev, min, max, len(ds))

        print ""
        print ""
        print "Near atoms, using atom's element name"
        print ""
        print ""
        print "Ion/H2O\tElement\tMean\tStDev\tMin\tMax\t# counted"
        for rtype, nnames in D_.iteritems () :
            rt = "H2O" if rtype == "HOH" else rtype
            rt = self.IonLabel ( rt )
            print ""
            print "%s" % rt,
            for nname, ds in nnames.iteritems() :
                mean, stdev, min, max = numpy.mean (ds), numpy.std (ds), numpy.min (ds), numpy.max(ds)
                print "\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%d" % (nname, mean, stdev, min, max, len(ds))


    def IonLabel ( self, rt ) :
        if rt.upper() in chargedIons :
            if chargedIons[rt.upper()] > 0 :
                rt += "(+%d)" % chargedIons[rt.upper()]
            else :
                rt += "(%d)" % chargedIons[rt.upper()]
        return rt

    def AddDs ( self, mol, D, D_, log=False ) :

        print ""
        print "%s" % mol.name

        ats = [at for at in mol.atoms if not at.element.name == "H"]
        import gridm; reload(gridm)
        atGrid = gridm.Grid ()
        atGrid.FromAtomsLocal ( ats, 3.5 )

        for res in mol.residues :
            if res.type == "HOH" or res.type.upper() in chargedIons :
                for at in res.atoms :
                    if at.element.name == "H" : continue
                    nats = atGrid.AtsNearPtLocal ( at.coord() )
                    for nat, v in nats :
                        if nat == at : continue

                        #if v.length > 2.5 and at.residue.type.upper() in chargedIons :
                        #    continue

                        if nat.element.name == "C" : continue

                        nname = nat.name
                        if nname == "OP1" or nname == "OP2" : nname = "OP"
                        if nname == "OE1" or nname == "OE2" : nname = "OE"
                        if nname == "OD1" or nname == "OD2" : nname = "OD"
                        if nname == "ND1" or nname == "ND2" : nname = "ND"
                        if nname == "NH1" or nname == "NH2" : nname = "NH"
                        if nname == "NE1" or nname == "NE2" : nname = "NE"
                        if nname == "BR1" or nname == "BR2" or nname == "BR3" : nname = "BR"

                        if nat.residue.type == "HOH" : nname = "O (H2O)"
                        nname = self.IonLabel ( nname )

                        if not at.residue.type in D : D[at.residue.type] = {}
                        if not nname in D[at.residue.type] : D[at.residue.type][nname] = []
                        D[at.residue.type][nname].append ( v.length )

                        ename = nat.element.name
                        if nat.residue.type == "HOH" : ename = "O (H2O)"
                        ename = self.IonLabel ( ename )

                        if not at.residue.type in D_ : D_[at.residue.type] = {}
                        if not ename in D_[at.residue.type] : D_[at.residue.type][ename] = []
                        D_[at.residue.type][ename].append ( v.length )

        if log :
            print "Ion/H2O\tAtom\tMean\tStDev\tMin\tMax\t# counted"
            for rtype, nnames in D.iteritems () :
                print "%s" % rtype,
                for nname, ds in nnames.iteritems() :
                    mean, stdev, min, max = numpy.mean (ds), numpy.std (ds), numpy.min (ds), numpy.max(ds)
                    print "%s\t%.2f\t%.2f\t%.2f\t%.2f\t%d" % (nname, mean, stdev, min, max, len(ds))








    def Stats ( self ) :

        print ""
        print "Test solvent atoms for Q-scores (make distributions) and"
        print "distances to other atoms"
        print ""

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
        import qscores
        minD, maxD = qscores.MinMaxD ( dmap )

        self.Log()

        umsg ( "Making statistics on ions and waters..." )

        ats = [at for at in mol.atoms if not at.element.name == "H"]
        #points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        #print " - search tree: %d/%d ats" % ( len(ats), len(mol.atoms) )
        #allAtTree = AdaptiveTree ( points.tolist(), ats, 2.0 )

        import gridm; reload(gridm)
        g1 = gridm.Grid ()
        g1.FromAtomsLocal ( ats, 5.0 )
        print " - %d atom grid" % len(ats)

        #points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        #import gridm
        #reload(gridm)
        #ptGrid = gridm.Grid()
        #ptGrid.FromPoints ( points, 4.0 )
        #print " - %d pts grid" % len(points)

        #import gridm; reload(gridm)
        #g1 = gridm.Grid ()
        #g1.FromAtoms ( ats, 4.0 )

        SetBBAts ( mol )
        doRes = []
        doAts = {}
        for res in mol.residues :
            # only looking for ions or water molecules which should have just one heavy atom
            if res.isProt or res.isNA :
                continue

            #if res.type == "HOH" or res.type.upper() in chargedIons :
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
        avgQs = {}
        ncMap = {}
        ncAts, ncsAts = [], []

        ncsW = {}
        ncsW["base"] = { "base":[], "sugar":[], "bb":[] };
        ncsW["sugar"] = { "base":[], "sugar":[], "bb":[] };
        ncsW["bb"] = {  "base":[], "sugar":[], "bb":[] }

        ncsMg = {}
        ncsMg["base"] = { "base":[], "sugar":[], "bb":[] };
        ncsMg["sugar"] = { "base":[], "sugar":[], "bb":[] };
        ncsMg["bb"] = {  "base":[], "sugar":[], "bb":[] }

        ncsW_1 = { "base":0, "sugar":0, "bb":0 }
        ncsMg_1 = { "base":0, "sugar":0, "bb":0 }

        NT = {}

        Qscores, QRes, VRes = {}, {}, {}
        Dists = {}
        DistVs = {}

        ats1, ats2, ats3 = [], [], []
        qs1, qs2, qs3 = [], [], []
        vs1, vs2, vs3 = [], [], []

        def addD ( tp, dist ) :
            if not tp in Dists :
                Dists[tp] = []
            Dists[tp].append ( dist )

        for res, atom in doRes:

            if 1 or not hasattr ( atom, 'Q' ) :
                #at.Q = 0.0
                #atom.Q = qscores.Qscore ( [atom], dmap, 0.6, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )

                #atPt = atom.coord()
                #atPt = [atPt.x, atPt.y, atPt.z]
                #xfI = atom.molecule.openState.xform
                #atom.Q = qscores.QscorePt3 ( atPt, xfI, dmap, 0.6, ptGrid=ptGrid, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )

                atom.Q = qscores.QscoreG ( [atom], dmap, 0.6, agrid=g1, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )

            if 1 :
                pt = [atom.coord()[0], atom.coord()[1], atom.coord()[2]]
                atom.V = dmap.interpolated_values ( [pt], dmap.openState.xform )
                #print map_values
                #print "|%.5f -> %.5f|" % (maxD, avgMapV)
                #break

            if atom.Q < 0.9 :
                #deletAts[atom] = 1
                #continue
                pass

            rtype = "H2O" if res.type.upper() == "HOH" else res.type.upper()

            if not rtype in Qscores : Qscores[rtype] = []
            Qscores[rtype].append ( atom.Q )


            itype = None
            if rtype.upper() in chargedIons :
                itype = "%d" % chargedIons[rtype.upper()]
            else :
                itype = rtype

            if not itype in avgQs :
                avgQs[itype] = [ atom.Q, 1.0 ]
            else :
                avgQs[itype][0] += atom.Q
                avgQs[itype][1] += 1.0

            #nearAts = self.AtsWithin ( [atom], 6.0, allAtTree )
            #nearAts = allAtTree.searchTree ( atom.coord().data(), 6.0 )

            closestD, closestAt = 7.0, None
            numClose, numCloseSolvent = 0, 0
            numCloseBBW, numCloseSugarW, numCloseBaseW = 0, 0, 0
            numCloseBBI, numCloseSugarI, numCloseBaseI = 0, 0, 0

            closeBBW, closeSugarW, closeBaseW = {}, {}, {}
            closeBBI, closeSugarI, closeBaseI = {}, {}, {}

            res.okNearAtoms = []
            #for nat in nearAts :
            nearAts = g1.AtsNearPtLocal ( atom.coord(), 5.0 )
            #if agrid.NumAtsNearAtLocal(at,D=outRad) < 1 :
            for nat, v in nearAts :

                if nat == atom :
                    continue

                #v = nat.coord() - atom.coord()
                d = v.length

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

                if d < closestD :
                    closestD, closestAt = d, nat

                if nat.element.name == "C" :
                    continue
                if nat.element.name == "H" :
                    continue

                nname = nat.name
                #if nname == "OP1" or nname == "OP2" : nname = "OP"

                isMG = res.type == "MG" and d < 2.5
                isHOH = res.type == "HOH" and d <= 3.5

                if isMG or isHOH :
                    if not nat.residue.type in NT :
                        NT[nat.residue.type] = {}
                    if not nname in NT[nat.residue.type] :
                        NT[nat.residue.type][nname] = {}
                    if not res.type in NT[nat.residue.type][nname] :
                        NT[nat.residue.type][nname][res.type] = [d, 1]
                    else :
                        NT[nat.residue.type][nname][res.type][0] += d
                        NT[nat.residue.type][nname][res.type][1] += 1

                if d < 3.2 :
                    numClose += 1
                    if 0 and nat.residue.type in chargedIons or nat.residue.type == "HOH" :
                        numCloseSolvent += 1
                    elif 0 and nat.isSC and nat.element.name == "O" :
                        numCloseSolvent += 1

                if res.type == "MG" and d < 2.5 :
                    if nat.isBB :
                        numCloseBBI += 1
                    elif nat.isSugar :
                        numCloseSugarI += 1
                    elif nat.isBase :
                        numCloseBaseI += 1
                    res.okNearAtoms += [ [d, nat] ]

                if res.type == "HOH" and d <= 3.5 :
                    if nat.isBB :
                        numCloseBBW += 1
                        closeBBW[nat.residue] = 1
                    elif nat.isSugar :
                        numCloseSugarW += 1
                    elif nat.isBase :
                        numCloseBaseW += 1
                    res.okNearAtoms += [ [d, nat] ]

            if res.type == "HOH" :

                if numCloseBaseW + numCloseSugarW + numCloseBBW <= 1 :
                    ncsW_1["base"] += numCloseBaseW
                    ncsW_1["sugar"] += numCloseSugarW
                    ncsW_1["bb"] += numCloseBBW

                if numCloseBaseW > 1 :
                    ncsW["base"]["base"] += [ [atom.Q, res] ]
                if numCloseBaseW > 0 and numCloseSugarW > 0 :
                    ncsW["base"]["sugar"] += [ [atom.Q, res] ]
                    ncsW["sugar"]["base"] += [ [atom.Q, res] ]
                if numCloseBaseW > 0 and numCloseBBW > 0 :
                    ncsW["base"]["bb"] += [ [atom.Q, res] ]
                    ncsW["bb"]["base"] += [ [atom.Q, res] ]

                if numCloseSugarW > 1 :
                    ncsW["sugar"]["sugar"] += [ [atom.Q, res] ]
                if numCloseSugarW > 0 and numCloseBBW > 0 :
                    ncsW["sugar"]["bb"] += [ [atom.Q, res] ]
                    ncsW["bb"]["sugar"] += [ [atom.Q, res] ]

                #if numCloseBBW > 1 :
                if len(closeBBW.keys()) > 1 :
                    ncsW["bb"]["bb"] += [ [atom.Q, res] ]
                    #print " HOH bb-bb %d.%s" % (res.id.position, res.id.chainId)

            elif res.type == "MG"  :

                if numCloseBaseI + numCloseSugarI + numCloseBBI <= 1 :
                    ncsMg_1["base"] += numCloseBaseI
                    ncsMg_1["sugar"] += numCloseSugarI
                    ncsMg_1["bb"] += numCloseBBI
                    if 0 and numCloseBaseI > 0 :
                        print " Mg-base: %d.%s" % (res.id.position, res.id.chainId)

                if numCloseBaseI > 1 :
                    ncsMg["base"]["base"] += [ [atom.Q, res] ]
                if numCloseBaseI > 0 and numCloseSugarI > 0 :
                    ncsMg["base"]["sugar"] += [ [atom.Q, res] ]
                    ncsMg["sugar"]["base"] += [ [atom.Q, res] ]
                if numCloseBaseI > 0 and numCloseBBI > 0 :
                    ncsMg["base"]["bb"] += [ [atom.Q, res] ]
                    ncsMg["bb"]["base"] += [ [atom.Q, res] ]

                if numCloseSugarI > 1 :
                    ncsMg["sugar"]["sugar"] += [ [atom.Q, res] ]
                if numCloseSugarI > 0 and numCloseBBI > 0 :
                    ncsMg["sugar"]["bb"] += [ [atom.Q, res] ]
                    ncsMg["bb"]["sugar"] += [ [atom.Q, res] ]

                if numCloseBBI > 1 :
                    ncsMg["bb"]["bb"] += [ [atom.Q, res] ]

            if numClose in ncMap :
                ncMap[numClose] += 1
            else :
                ncMap[numClose] = 1
            ncAts.append ( [numClose, atom] )
            ncsAts.append ( [numCloseSolvent, atom] )

            #if res.type == "MG" :
            #    print "MG %d.%s - %.2f" % (res.id.position, res.id.chainId, closestD)


            if closestAt != None :
                nat = closestAt
                d = closestD

                pt = [nat.coord()[0], nat.coord()[1], nat.coord()[2]]
                closestAt.V = dmap.interpolated_values ( [pt], dmap.openState.xform )[0]
                #pts = [ at.coord() for at in nat.residue.atoms if not at.element.name == "H" ]
                #Vs = dmap.interpolated_values ( pts, dmap.openState.xform )
                #closestAt.V = numpy.sum(Vs) / float ( len(Vs) )

                if d < 2.2 :
                    qs1.append ( atom.Q )
                    vs1.append ( atom.V )
                elif d < 2.8 :
                    qs2.append ( atom.Q )
                    vs2.append ( atom.V )
                else :
                    qs3.append ( atom.Q )
                    vs3.append ( atom.V )

                isMG = res.type == "MG" and d <= 2.2
                isHOH = res.type == "HOH" and d >= 2.8 and d <= 3.5

                isProt = False
                #if closestAt.element.name == "N" :
                #    for bat in closestAt.neighbors :
                #        if bat.element.name == "H" :
                #            isProt = True

                if isMG :
                    if not rtype in QRes : QRes[rtype] = []
                    QRes[rtype].append ( [atom.Q, res.id.position] )
                    if not rtype in VRes : VRes[rtype] = []
                    VRes[rtype].append ( [atom.V-closestAt.V, res.id.position] )
                if isHOH :
                    if 1 or isProt :
                        if not rtype in QRes : QRes[rtype] = []
                        QRes[rtype].append ( [atom.Q, res.id.position] )
                        if not rtype in VRes : VRes[rtype] = []
                        VRes[rtype].append ( [atom.V-closestAt.V, res.id.position] )

                if 1 or isMG or (isHOH and isProt) :
                    dround = "%.1f" % d
                    if not dround in DistVs :
                        DistVs[dround] = []
                    DistVs[dround].append ( atom.V-closestAt.V )

                if nat.element.name == "O" :
                    if 1 and nat.residue.type == "HOH" :
                        addD ( "%s-H2O" % rtype, d )
                        addD ( "H2O-%s" % rtype, d )
                        if d < 2.0 :
                            r = res
                            nr = nat.residue
                            print " -- res %d.%s %s - %d.%s %s - d %.2f - occ %.2f, %.2f" % (r.id.position, r.id.chainId, r.type, nr.id.position, nr.id.chainId, nr.type, d, atom.occupancy, nat.occupancy)
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

                    #if rtype == "MG" and nat.element.name == "N" and d < 3.0 :
                    #    print "%.2f - %s - %s.%d.%s.%s" % (d, rtype, nat.name, nat.residue.id.position, nat.residue.type, nat.residue.id.chainId)

                #if nat.element.name == "C" :
                #    addD ( "%s-C" % rtype, d )


            atI += 1
            if atI % 30 == 0 :
                status ( "Making statistics on ions and waters... at %d/%d" % (atI, len(doRes)) )
                #print ".",


        status ( "Statistics done - open Log (IDLE) for results" )

        print ""
        print " - # same: %d" % numSame


        print ""
        print " - %d ats to delete" % len(deletAts.keys())
        if 0 :
            print " - doing delete" % len(deletAts.keys())
            for at in deletAts.keys() :
                #print " - %s in res %s %d chain %s" % (at.name, at.residue.type, at.residue.id.position, at.residue.id.chainId)
                if len(at.residue.atoms) == 1 :
                    mol.deleteResidue ( at.residue )
                else :
                    mol.deleteAtom ( at )


        print ""
        print ""
        print "Distances:"

        edges = numpy.array ( range ( 0, 21 ) ) / 5.0 + 0.1

        s = ""
        for ei, e in enumerate ( edges ) :
            s = s + "\t%.1f-%.1f" % (e, e+0.2)
        print s

        s = ""
        for e in edges :
            s = s + "\t%.2f" % e
        print s

        done = {}
        types = Dists.keys()
        types.sort()
        for typ in types :
            t1, t2 = typ.split("-")
            if "%s-%s" % (t2, t1) in done :
                continue
            done[typ] = 1
            f = 2 if t1 == t2 else 1 # counted twice if t1 == t2
            s = typ
            hist = numpy.histogram ( Dists[typ], bins=edges )[0]
            #print hist
            for n in hist :
                s = s + "\t%d" % (n/f)
            print s


        edges = numpy.array ( range ( 0, 11 ) ) / 10.0

        print ""
        print "Q-scores:"
        s = ""
        for e in edges :
            s = s + "\t%.1f-%.1f" % (e, e+0.1)
        print s

        for tp, qscores in Qscores.iteritems () :
            hist = numpy.histogram ( qscores, bins=edges )[0]
            s = tp
            for n in hist :
                s = s + "\t%d" % n
            print s

        print ""
        print "Type\tAvg.Q\tNum"
        print mol.name.replace ( ".maxit.pdb", "" )
        for tp, qn in avgQs.iteritems() :
            print "%s\t%.2f\t%.0f" % ( tp, qn[0]/qn[1], qn[1] )

        print ""

        print "#close atoms\t# ions/H2O"
        for nc, num in ncMap.iteritems () :
            print "%d\t%d" % ( nc, num )

        print ""

        if 0 :
            print "\nBy close atoms..."
            ncAts.sort ( reverse=True, key=lambda x: x[0] )
            for nc, at in ncAts [1:10] :
                print "%d - %s(%s),%d.%s" % (nc, at.name, at.residue.type, at.residue.id.position, at.residue.id.chainId)

        if 0 :
            print "\nBy close atoms (solvent)..."
            ncsAts.sort ( reverse=True, key=lambda x: x[0] )
            for nc, at in ncsAts [1:10] :
                print "%d - %s(%s),%d.%s" % (nc, at.name, at.residue.type, at.residue.id.position, at.residue.id.chainId)

        if 0 :
            print "\nW"
            for t1 in ["base", "sugar", "bb"] :
                print "%s\t%d" % (t1, ncsW_1[t1])

            print "Base\tSugar\tBB"
            for t1 in ["base", "sugar", "bb"] :
                for t2 in ["base", "sugar", "bb"] :
                    print "%d\t" % len(ncsW[t1][t2]),
                print ""
            print ""

            print "\nMg"
            for t1 in ["base", "sugar", "bb"] :
                print "%s\t%d" % (t1, ncsMg_1[t1])
            print "Base\tSugar\tBB"
            for t1 in ["base", "sugar", "bb"] :
                for t2 in ["base", "sugar", "bb"] :
                    print "%d\t" % len(ncsMg[t1][t2]),
                print ""
            print ""

        if 0 :
            print ""
            for i1, t1 in enumerate ( ["base", "sugar", "bb"] ) :
                for i2, t2 in enumerate ( ["base", "sugar", "bb"] ) :
                    if i2 >= i1 :
                        qrs = ncsW[t1][t2]
                        qrs.sort ( reverse=True, key=lambda x: x[0] )
                        print "%s - H2O - %s [%d]:" % (t1, t2, len(qrs))
                        for q, res in qrs [0:20] :
                            nats = res.okNearAtoms
                            print "%.2f/%d.%s/%d" % (q, res.id.position, res.id.chainId, len(nats)),
                            #for nat in nats :
                            #    print ".%s-%d" % (nat.name, nat.residue.id.position),
                        print ""
                        print ""

            print ""
            for i1, t1 in enumerate ( ["base", "sugar", "bb"] ) :
                for i2, t2 in enumerate ( ["base", "sugar", "bb"] ) :
                    if i2 >= i1 :
                        qrs = ncsMg[t1][t2]
                        qrs.sort ( reverse=True, key=lambda x: x[0] )
                        print "%s - ion - %s [%d]:" % (t1, t2, len(qrs)),
                        for q, res in qrs [0:20] :
                            nats = res.okNearAtoms
                            print "%.2f/%d.%s/%d" % (q, res.id.position, res.id.chainId, len(nats)),
                        print ""


        if 0 :
            for ntype, nats in NT.iteritems () :
                if ntype == "A" or ntype == "G" or ntype == "C" or ntype == "U" :
                    print ""
                    print ntype
                    print "At.Name\t# MG\t# H2O"
                    for nat in nats.keys () :
                        dMG, numMG = nats[nat]["MG"] if "MG" in nats[nat] else [0, 0]
                        dHOH, numHOH = nats[nat]["HOH"] if "HOH" in nats[nat] else [0, 0]
                        print "%s\t%d\t%d" % (nat, numMG, numHOH)




        if 0 :
            if not os.path.isfile ("/Users/greg/Desktop/txt.txt") :
                fp = open ( "/Users/greg/Desktop/txt.txt", "a" )
                fp.write ( "Mol Name\t_AvgQ(H2O)_\t_#(H2O)_\t_AvgQ(+2)_\t_#(+2)_\t_AvgQ(+1)_\t_#(+1)_\t_AvgQ(-1)_\t_#(-1)_\t_AvgQ(+3)_\t_#(+3)_\n" )
                fp.close()
            fp = open ( "/Users/greg/Desktop/txt.txt", "a" )
            fp.write ( "%s" % mol.name.replace(".maxit.pdb", "") )
            for tp in ["H2O", "2", "1", "-1", "3"] :
                if tp in avgQs :
                    qn = avgQs[tp]
                    fp.write ( "\t%.2f\t%.0f" % ( qn[0]/qn[1], qn[1] ) )
                else :
                    fp.write (  "\t\t" )
            fp.write ( "\n" )
            fp.close()


        print ""
        print "Top 10 by Q-score / type"
        for tp, qres in QRes.iteritems () :
            print "%s" % tp,
            qres.sort ( reverse=True )
            for q, ri in qres[0:10] :
                print "\t%.2f\t%d" % (q, ri)

        print ""
        print "Q-score stats by type"
        for tp, qres in QRes.iteritems () :
            qs = [item[0] for item in qres]
            mean, std, num = numpy.mean(qs), numpy.std(qs), len(qs)
            print "%s\t%.4f\t%.4f\t%d" % (tp, mean, std, num)

        print ""
        print "Density value stats by type"
        for tp, vres in VRes.iteritems () :
            vs = [item[0] for item in vres]
            mean, std, num = numpy.mean(vs), numpy.std(vs), len(vs)
            print "%s\t%.4f\t%.4f\t%d" % (tp, mean, std, num)

        print ""
        print "Density values by distance..."
        for dr, vs in DistVs.iteritems() :
            mean, std, num = numpy.mean(vs), numpy.std(vs), len(vs)
            print "%s\t%.4f\t%.4f\t%d" % (dr, mean, std, num)

        print ""
        print "1.8 -- 2.2 %d %d" % ( len(qs1), len(vs1) )
        print "%.3f\t%.3f\t%.3f\t%.3f" % ( numpy.mean(qs1), numpy.std(qs1), numpy.mean(vs1), numpy.std(vs1) )

        print "2.2 -- 2.8  %d %d" % ( len(qs2), len(vs2) )
        print "%.3f\t%.3f\t%.3f\t%.3f" % ( numpy.mean(qs2), numpy.std(qs2), numpy.mean(vs2), numpy.std(vs2) )

        print "2.8 -- 3.6  %d %d" % ( len(qs3), len(vs3) )
        print "%.3f\t%.3f\t%.3f\t%.3f" % ( numpy.mean(qs3), numpy.std(qs3), numpy.mean(vs3), numpy.std(vs3) )

    def S1 ( self ) :

        print ""

        for m in chimera.openModels.list() :

            if type(m) != chimera.Molecule :
                continue

            print "\n\n-------------- %s -------- " % m.name

            self.cur_mol = m
            self.Stats()



    def SN ( self ) :

        print ""

        mols = []
        for m in chimera.openModels.list() :
            if type(m) == chimera.Molecule :
                mols.append ( m )

        tp = "HOH"

        fp = open ( "/Users/greg/Desktop/%s.txt" % tp, "w" )
        fp.write ( "\t" )
        for m in mols :
            if m.display == False : continue
            mname = m.name.replace(".maxit.pdb", "").replace("T020", "").replace("EM0", "_")
            fp.write ( "\t%s" % mname )
        fp.write ( "\n\t" )
        for m in mols :
            if m.display == False : continue
            ats = self.wiAtoms ( m, tp )
            fp.write ( "\t%d" % len(ats) )
        fp.write ( "\n" )

        for m1 in mols :

            if m1.display == False :
                continue

            ats1 = self.wiAtoms ( m1, tp )
            mname = m1.name.replace(".maxit.pdb", "").replace("T020", "").replace("EM0", "_")
            fp.write ( "%s\t%d" % (mname, len(ats1)) )

            for m2 in mols :

                if m2.display == False :
                    continue

                num, pp, rmsd = self.wiDistNum ( tp, m1, m2)
                #fp.write ( "\t%d (%.0f%%) %.2f" % (num, pp, rmsd) )
                #fp.write ( "\t%d (%.0f%%)" % (num, pp) )
                #fp.write ( "\t%.0f" % pp )
                fp.write ( "\t%d" % num )
                #break

            fp.write ( "\n" )

        fp.close()



    def SNi ( self ) :

        print ""

        mols = []
        for m in chimera.openModels.list() :
            if type(m) == chimera.Molecule :
                mols.append ( m )

        tp = "MG"

        fp = open ( "/Users/greg/Desktop/%s.txt" % tp, "w" )
        fp.write ( "\t" )
        for m in mols :
            if m.display == False : continue
            mname = m.name.replace(".mrc", "").replace("Con2-", "").replace("Con1_", "_")
            fp.write ( "\t%s" % mname )
        fp.write ( "\n\t" )
        for m in mols :
            if m.display == False : continue
            ats = self.wiAtoms ( m, tp )
            fp.write ( "\t%d" % len(ats) )
        fp.write ( "\n" )

        for m1 in mols :

            if m1.display == False :
                continue

            ats1 = self.wiAtoms ( m1, tp )
            mname = m1.name.replace(".mrc", "").replace("Con2-", "").replace("Con1_", "_")
            fp.write ( "%s\t%d" % (mname, len(ats1)) )

            for m2 in mols :

                if m2.display == False :
                    continue

                num, pp, rmsd = self.wiDistNum ( tp, m1, m2)
                #fp.write ( "\t%d (%.0f%%) %.2f" % (num, pp, rmsd) )
                #fp.write ( "\t%d (%.0f%%)" % (num, pp) )
                #fp.write ( "\t%.0f" % pp )
                fp.write ( "\t%d" % num )
                #break

            fp.write ( "\n" )

        fp.close()


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
        SetBBAts ( mol )

        dmap = self.cur_dmap
        if self.cur_dmap == None :
            umsg ("Select a map first")
            return []

        #print " -- chain %s" % chainId
        toChain = self.addToChain.get()
        if len(toChain) > 1 :
            umsg ( "Enter a single character in 'Add To Chain' field" )
            return

        nearAtoms = chimera.selection.currentAtoms()
        if len(nearAtoms) == 0 :
            umsg ( "Select atoms near which ions/waters should be placed" )
            return []

        # min/max distance for water (W) and ion (I) from GUI
        # by default, ion is min:1.8 max:2.5, water is min:2.5 max:3.5
        minDistW, maxDistW = float(self.waterMinD.get()), float(self.waterMaxD.get())
        minDistI, maxDistI = float(self.ionMinD.get()), float(self.ionMaxD.get())

        umsg ( "Placing water/ions in map: %s, model: %s ... " % (segMap.name, mol.name) )
        if task :
            print " - got task..."
            task.updateStatus( "Placing water/ions in map: %s, model: %s ... " % (segMap.name, mol.name) )

        #mapBase = os.path.splitext ( segMap.name )[0]
        hMapA = self.cur_dmap_h1 # GetMod_ ( mapBase + "_half_A.mrc" )
        hMapB = self.cur_dmap_h2 # GetMod_ ( mapBase + "_half_B.mrc" )
        if hMapA and hMapB :
            print " - using half map A:", hMapA.name
            print " - using half map B:", hMapB.name

        useQ = self.useQScore.get()
        try :
            minQ = float(self.placeQ.get())
            sigQ = float(self.qsigma.get())
        except :
            umsg ( "Check Q-score and sigma, should be numbers..." )
            return


        umsg ( "Placing water/ions in map: %s, model: %s, %d regions ... " % (segMap.name, mol.name, len(smod.regions)) )

        if useQ :
            print " - using min Q-score: %.2f sigma %.2f" % (minQ, sigQ)
            if hMapA :
                print "   - in half map A: %s" % hMapA.name
            if hMapB :
                print "   - in half map B: %s" % hMapB.name

        addW, addI = goSWIM ( segMap, smod, mol, nearAtoms, toChain, minDistI, maxDistI, minDistW, maxDistW, hMapA, hMapB, useQ, minQ, sigQ, task )

        status ( "Added %d waters, %d ions - done" % (len(addW), len(addI)) )

        if 1 :

            qs = ""
            if useQ : qs += "_Q%.2f_%.2f" % (minQ, sigQ)
            if hMapA : qs += "__hA"
            if hMapB : qs += "__hB"

            thr = 0.0
            if hasattr ( segMap, "segmentThreshold" ) :
                thr = segMap.segmentThreshold
            else :
                thr = segMap.surface_levels[0]

            mapName = os.path.splitext ( segMap.name )[0]
            print " - map name: %s" % mapName

            molPath = os.path.splitext(mol.openedAs[0])[0]
            nname = molPath + "__%s__thr%.3f%s__%d-water__%d-ion.pdb" % (mapName, thr, qs, len(addW), len(addI))

            print ""
            print "Saving pdb waters ->", nname
            chimera.PDBio().writePDBfile ( [mol], nname )

            print ""


        self.RefreshTree ()


    def Thr ( self ) :


        #mol = self.cur_mol
        #if self.cur_mol == None :
        #    umsg ("Select a molecule first")
        #    return []

        #chainId = self.chain.get()

        dmap = self.cur_dmap
        print " - scale map: %s" % dmap.name

        if dmap == None :
            umsg ( "Select a map first..." )
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


        #dmap.surface_levels[0] = sig2
        #chimera.runCommand ( "vol #%d style surface region all step 1" % dmap.id )



    def Q ( self ) :

        dmap = self.cur_dmap
        if dmap == None :
            umsg ( "Select a map first..." )
            return

        atoms = chimera.selection.currentAtoms ()
        if len(atoms) == 0 :
            umsg ( "Select some atoms to calcualte Q-scores for" )
            return

        mol = atoms[0].molecule
        #selAtom = selAts[0]
        #r = selAtom.residue
        #print "Res: %s - %d.%s - %s - Atom: %s" % (r.type, r.id.position, r.id.chainId, r.molecule.name, selAtom.name)

        sigma = float(self.qsigma.get())

        from time import time

        #sigma = 0.4
        #print " - in map: %s" % self.cur_dmap.name
        #print " - mol: %s" % mol.name
        #print " - sigma: %.2f" % sigma

        ats = [at for at in mol.atoms if not at.element.name == "H"]

        import gridm; reload(gridm)
        g1 = gridm.Grid ()
        gstart = time()
        g1.FromAtomsLocal ( ats, 3.0 )
        gend = time()
        print " - %d ats in %.3f sec" % (len(ats), gend-gstart)

        if 0 :

            start = time()
            allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)

            print "Tree %.3f sec:" % (time()-start)
            opointsNear = allAtTree.searchTree ( atoms[0].coord(), 3.0 )
            foundNearPt = False
            for at in opointsNear :
                v = atoms[0].coord() - at.coord()
                print " - %s.%d.%s -- %.3f" % (at.name, at.residue.id.position, at.residue.id.chainId, v.length)

            print "Grid %.3f sec:" % (gend-gstart)
            nearAts = g1.AtsNearPtLocal ( atoms[0].coord(), 3.0 )
            #if agrid.NumAtsNearAtLocal(at,D=outRad) < 1 :
            for at, v in nearAts :
                print " - %s.%d.%s -- %.3f" % (at.name, at.residue.id.position, at.residue.id.chainId, v.length)

            print "?"


        minD, maxD = qscores.MinMaxD ( dmap )
        #print " - minD %.3f, maxD %.3f" % (minD, maxD)


        if 1 :

            points = _multiscale.get_atom_coordinates ( ats, transformed = False )
            allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)

            start = time()
            msg = ""
            qatoms = [at for at in atoms if not at.element.name == "H"]
            for ai, at in enumerate ( qatoms ) :
                at.Q = qscores.Qscore ( [at], dmap, sigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
                msg += "%s.%d.%s:%.3f" % (at.name, at.residue.id.position, at.residue.id.chainId, at.Q)
                if self.cur_dmap_h1 != None :
                    q1 = qscores.Qscore ( [at], self.cur_dmap_h1, sigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
                    msg += "/%.3f" % (q1)
                if self.cur_dmap_h2 != None :
                    q2 = qscores.Qscore ( [at], self.cur_dmap_h2, sigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
                    msg += "/%.3f" % (q2)
                msg += " "

            #umsg ( msg )
            dur = time()-start
            print msg, "%.3fsec" % dur

        if 1 :
            start = time()
            msg = ""
            qatoms = [at for at in atoms if not at.element.name == "H"]
            for ai, at in enumerate ( qatoms ) :
                at.Q = qscores.QscoreG ( [at], dmap, sigma, agrid=g1, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
                msg += "Q-score of %s.%d.%s : %.3f" % (at.name, at.residue.id.position, at.residue.id.chainId, at.Q)
                if self.cur_dmap_h1 != None :
                    q1 = qscores.QscoreG ( [at], self.cur_dmap_h1, sigma, agrid=g1, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
                    msg += ", in Map A: %.3f" % (q1)
                if self.cur_dmap_h2 != None :
                    q2 = qscores.QscoreG ( [at], self.cur_dmap_h2, sigma, agrid=g1, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
                    msg += ", in Map B:%.3f" % (q2)
                msg += " "

            dur = time()-start
            print msg, "%.3fsec" % dur
            umsg ( msg )




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





def PtsToMap ( points, dmap, atomRad, nname, showMesh = False, clr = (0.7, 0.7, 0.7, 0.2) ) :

    #_contour.affine_transform_vertices ( points, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
    #mdata = VolumeData.zone_masked_grid_data ( dmap.data, points, atomRad )

    import _contour
    points1 = numpy.copy ( points )
    _contour.affine_transform_vertices ( points1, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
    points0 = numpy.copy ( points1 )
    _contour.affine_transform_vertices ( points1, dmap.data.xyz_to_ijk_transform )

    bound = 5
    li,lj,lk = numpy.min ( points1, axis=0 ) - (bound, bound, bound)
    hi,hj,hk = numpy.max ( points1, axis=0 ) + (bound, bound, bound)

    n1 = hi - li + 1
    n2 = hj - lj + 1
    n3 = hk - lk + 1

    #print " - bounds - %d %d %d --> %d %d %d --> %d %d %d" % ( li,lj,lk, hi,hj,hk, n1,n2,n3 )

    #nmat = numpy.zeros ( (n1,n2,n3), numpy.float32 )
    #dmat = dmap.full_matrix()

    nstep = (dmap.data.step[0], dmap.data.step[1], dmap.data.step[2] )
    #nstep = (fmap.data.step[0]/2.0, fmap.data.step[1]/2.0, fmap.data.step[2]/2.0 )

    nn1 = int ( round (dmap.data.step[0] * float(n1) / nstep[0]) )
    nn2 = int ( round (dmap.data.step[1] * float(n2) / nstep[1]) )
    nn3 = int ( round (dmap.data.step[2] * float(n3) / nstep[2]) )

    O = dmap.data.origin
    #print " - %s origin:" % dmap.name, O
    nO = ( O[0] + float(li) * dmap.data.step[0],
           O[1] + float(lj) * dmap.data.step[1],
           O[2] + float(lk) * dmap.data.step[2] )

    #print " - new map origin:", nO

    nmat = numpy.zeros ( (nn1,nn2,nn3), numpy.float32 )
    ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles )

    #print " - fmap grid dim: ", numpy.shape ( fmap.full_matrix() )
    #print " - new map grid dim: ", numpy.shape ( nmat )

    npoints = VolumeData.grid_indices ( (nn1, nn2, nn3), numpy.single)  # i,j,k indices
    _contour.affine_transform_vertices ( npoints, ndata.ijk_to_xyz_transform )

    dvals = dmap.interpolated_values ( npoints, dmap.openState.xform )
    #dvals = dmap.interpolated_values ( npoints, chimera.Xform.identity() )
    #dvals = dmap.interpolated_values ( npoints, dmap.openState.xform.inverse() )
    #dvals = numpy.where ( dvals > threshold, dvals, numpy.zeros_like(dvals) )
    #nze = numpy.nonzero ( dvals )

    nmat = dvals.reshape( (nn3,nn2,nn1) )
    #f_mat = fmap.data.full_matrix()
    #f_mask = numpy.where ( f_mat > fmap.surface_levels[0], numpy.ones_like(f_mat), numpy.zeros_like(f_mat) )
    #df_mat = df_mat * f_mask

    ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles )
    #try : nv = VolumeViewer.volume.add_data_set ( ndata, None )
    #except : nv = VolumeViewer.volume.volume_from_grid_data ( ndata )

    #nv.openState.xform = dmap.openState.xform

    mdata = VolumeData.zone_masked_grid_data ( ndata, points0, atomRad )
    gdata = VolumeData.Array_Grid_Data ( mdata.full_matrix(), nO, nstep, dmap.data.cell_angles, name = nname )
    nv = VolumeViewer.volume.volume_from_grid_data ( gdata )
    nv.openState.xform = dmap.openState.xform

    nv.name = nname
    #dmap.display = False
    nv.region = ( nv.region[0], nv.region[1], [1,1,1] )

    nv.surface_levels[0] = dmap.surface_levels[0]

    ro = VolumeViewer.volume.Rendering_Options()
    if 1 :
        ro.smoothing_factor = .3
        ro.smoothing_iterations = 2
        ro.surface_smoothing = True
    #ro.square_mesh = True
    #ro.line_thickness = 2
    nv.update_surface ( False, ro )
    #setro (ro)
    for sp in nv.surfacePieces :
        v, t = sp.geometry
        if len(v) == 8 and len(t) == 12 :
            sp.display = False
        else :
            if showMesh :
                sp.color = (.5, .5, .5, 1.0)
                sp.displayStyle = sp.Mesh
            else :
                sp.color = clr

    return nv





def GuessAtom ( mol, P, atGrid=None, nearAtMap=None, doMsg=True, minDistI=1.8, maxDistI=2.5, minDistW=2.5, maxDistW=3.4, ionType="MG" ) :

    nearAts = None

    # find the nearest atoms (in mol)
    if atGrid != None :
        # if a grid is given, use it as it will be quick
        nearAts = atGrid.AtsNearPtLocal ( P )
        #print "%d" % len(nearAts),

    else :
        # otherwise do slower all-atom search
        #nearAts = [None] * len(mol.atoms)
        nearAts = []
        P = chimera.Point ( P[0], P[1], P[2] )
        for i, at in enumerate(mol.atoms) :
            if not at.element.name == "H" :
                V = P - at.coord()
                if V.length < 6.0 :
                    nearAts.append ( [at, V] )


    #R = lambda : None
    isNearAtMap = False
    collidingAtoms = []
    closestChainId, closestChainD = None, 1e9

    # number of nearby atoms within ion distance (positive or negative)
    posAtomsIonD, negAtomsIonD = [], []
    # number of nearby atoms within water distance (positive or negative)
    posAtomsWaterD, negAtomsWaterD = [], []
    # number of nearby ions within water distance (positive or negative)
    ionAtomsIonD, ionAtomsWaterD = [], []

    #R.hbAcceptors, R.hbDonors = [], []

    for at, v in nearAts :

        dist = v.length
        if at.element.name == "H" : continue
        #if hasattr ( at, 'Q' ) and at.Q < 0.1 : continue
        #if at.altLoc != '' : continue

        # if carbon atom too close, mark as clash/collision
        if at.element.name == "C" and dist < 2.6 :
            collidingAtoms.append ( [dist, at] )
            #print "c",

        # for other atoms, if close than minDistI, mark as collision
        if dist < minDistI :
            collidingAtoms.append ( [dist, at] )

        # check if close to a selected atom, only placing water/ions next to these
        if nearAtMap != None and at in nearAtMap :
            isNearAtMap = True

        # keep track of which chain is closest, the placed/water will be in this chain
        # unless otherwise specified
        if dist < closestChainD :
            closestChainId, closestChainD = at.residue.id.chainId, dist

        # add to count depending on nearby atom type and distance
        if at.residue.type.upper() in chargedIons :
            if dist < maxDistI : ionAtomsIonD.append ( [dist, at] )
            elif dist < maxDistW : ionAtomsWaterD.append ( [dist, at] )
        elif at.element.name == "N" :
            hAts = [a for a in at.bondsMap.keys() if a.element.name == "H"]
            if len(hAts) > 0 :
                if dist < maxDistI : posAtomsIonD.append ( [dist, at] )
                elif dist < maxDistW : posAtomsWaterD.append ( [dist, at] )
            else :
                if dist < maxDistI : negAtomsIonD.append ( [dist, at] )
                elif dist < maxDistW : negAtomsWaterD.append ( [dist, at] )
        elif at.element.name == "O" and at.residue.type.upper() == "HOH" :
            # todo... look at coordination?
            pass
        elif at.element.name == "O" or (at.element.name == "S" and at.residue.type == "CYS") :
            if dist < maxDistI : negAtomsIonD.append ( [dist, at] )
            elif dist < maxDistW : negAtomsWaterD.append ( [dist, at] )

    # generate a message listing nearby atoms
    msg = ""
    if doMsg :
        if len(collidingAtoms) > 0 :
            msg = "Clash:"
            for d, at in collidingAtoms :
                msg += " " + At (at, d)
                at.display = True

        else :
            msg += "Near: "

            if len(negAtomsIonD) > 0 :
                msg += "\n\n (-) Atoms (at Ion distance):"
                for d, at in negAtomsIonD :
                    msg += " \n" + At (at, d)

            if len(posAtomsIonD) > 0 :
                msg += "\n\n (+) Atoms (at Ion distance):"
                for d, at in posAtomsIonD :
                    msg += " \n" + At (at, d)

            if len(ionAtomsIonD) > 0 :
                msg += "\n\nIon (at Ion distance):"
                for d, at in ionAtomsIonD :
                    msg += " \n" + At (at, d)

            if len(negAtomsWaterD) > 0 :
                msg += "\n\n (-) Atoms (at Water distance):"
                for d, at in negAtomsWaterD :
                    msg += " \n" + At (at, d)

            if len(posAtomsWaterD) > 0 :
                msg += "\n\n (+) Atoms (at Water distance):"
                for d, at in posAtomsWaterD :
                    msg += " \n" + At (at, d)

            if len(ionAtomsWaterD) > 0 :
                msg += "\n\nIon (at Water distance):"
                for d, at in ionAtomsWaterD :
                    msg += " \n" + At (at, d)


    atName, atRes = None, None
    clr = None
    placedType = ""

    if isNearAtMap == False and nearAtMap != None :
        # a nearAtMap was given so only take points close to those atoms
        # in this case it was not close to any of them...
        pass

    elif len(collidingAtoms) == 0 :

        if len(negAtomsIonD) > 0 and len(posAtomsIonD) == 0 and len(posAtomsWaterD) == 0 :
            # next to atom, ion distance away, no positive atoms nearby (H atoms)
            atName, atRes = ionType, ionType
            placedType = "2+ ion"
            clr = (0,1,0)
        elif len(negAtomsIonD) > 0 :
            pass
        elif 0 and len (posAtomsIonD) > 0 :
            # next to positive atom
            atName, atRes = "CL", "CL"
            placedType = "1- ion"
            clr = (0,1,0)
        elif len(ionAtomsIonD) > 0 or len(ionAtomsWaterD) > 0 :
            # next to ion, put water
            atName, atRes = "O", "HOH"
            placedType = ""
            clr = (1,0,0)
        #elif len(negAtomsWaterD) >= 4 :
        #    # next to at least 4 atoms water distance away - can't be water?
        #    atName, atRes = ionType, ionType
        #    placedType = "2+ ion"
        #    clr = (0,1,0)
        elif len(negAtomsWaterD) > 0 or len(posAtomsWaterD) > 0 :
            # next to at least 1 atom, water distance away
            atName, atRes = "O", "HOH"
            placedType = ""
            clr = (1,0,0)

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



def RegPt ( reg, segMap, mode="tomax" ) :

    P, ctr, val = None, None, None
    # what point in the region to use...
    if mode == "ctr" :

        # use the center of all points in the region
        # - may not be close to highest value or peak
        ctr = reg.center_of_points()
        P = chimera.Point(ctr[0],ctr[1],ctr[2])

    elif mode == "hval" :
        # use the highest value in the region
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

        P = chimera.Point(maxValPt[0], maxValPt[1], maxValPt[2])
        ctr = [maxValPt[0], maxValPt[1], maxValPt[2]]
        val = maxD

    elif mode == "tomax" :

        #ctr = reg.center_of_points()
        #ctrP = chimera.Point(ctr[0], ctr[1], ctr[2])

        ctr, maxD = None, -1e9
        rpts = reg.map_points()
        map_values = segMap.interpolated_values ( rpts, segMap.openState.xform )
        #print map_values
        #break
        maxValPt = None
        for pt, val in zip(rpts, map_values) :
            if val > maxD :
                maxValPt, maxD = pt, val

        ctr = maxValPt
        ctrP = chimera.Point(ctr[0], ctr[1], ctr[2])

        #print pt
        # go to interpolated maximum...
        # - interestingly this can be different than the voxel with
        # - the highest value, due to interpolation used
        pts, avgMapV = PtsToMax ( [ctr], segMap )
        maxPt = pts[0]
        maxValP = chimera.Point ( maxPt[0], maxPt[1], maxPt[2] )

        #maxPt_ = [maxPt[0], maxPt[1], maxPt[2]]
        #map_values = segMap.interpolated_values ( [maxPt_], segMap.openState.xform )
        #print map_values
        #print "|%.5f -> %.5f|" % (maxD, avgMapV)
        #break

        # if movement is too large to maximum, likely this
        # is not a well-separated blob, so ignore it
        V = maxValP - ctrP
        #if maxD > avgMapV :
        #    print "|%.5f -> %.5f|" % (maxD, avgMapV),
        #    skipped += 1
        #    continue
        if V.length <= segMap.data.step[0] :
            #print "|%.1f|" % V.length,
            #skipped += 1
            #continue
            P = maxValP
            val = avgMapV
        else :
            P = None
            val = V.length

    elif mode == "com" :
        # use the center of mass
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

    return [P, val]



def goSWIM ( segMap, smod, mol, nearAtoms, toChain='', \
                minDistI=1.8, maxDistI=2.5, minDistW=2.5, maxDistW=3.5, \
                hMapA=None, hMapB=None, \
                useQ=True, minQ=0.7, sigQ=0.6, \
                task=None ) :

    ats = [at for at in mol.atoms if not at.element.name == "H"]
    import gridm; reload(gridm)
    atGrid = gridm.Grid ()
    atGrid.FromAtomsLocal ( ats, maxDistW )
    print " - made atoms grid with %d atoms" % len(ats)

    nearAtMap = {}
    for at in nearAtoms :
        nearAtMap[at] = 1

    minD, maxD = qscores.MinMaxD ( segMap )
    print " - map mind: %.3f, maxd: %.3f" % (minD, maxD)

    print " - processing regions..."
    regs = list(smod.regions)
    n_regs = []
    for reg in regs :
        npts = len(reg.points())
        if reg.surface_piece:
            reg.hide_surface()
        if npts > 3 :
            P, val = RegPt ( reg, segMap )
            if P != None :
                n_regs.append ( [val, P, reg] )

    # sort regions by value
    n_regs.sort ( reverse=True, key=lambda x: x[0] )

    addPts = []
    addW = []
    addI = []
    skipped, skippedQ, skippedQA, skippedQB = 0, 0, 0, 0
    numW, numI = 0, 0
    xfI = segMap.openState.xform

    # a temporary molecule - needed? new waters and ions are put to the grid
    nmol = chimera.Molecule()

    from time import time
    startt = time()

    regi = 0
    for val, P, reg in n_regs :

        modi = 10 if task else 1000
        if regi % modi == 0 :
            ts = qscores.TimeLeftStr (regi, len(n_regs), time() - startt )
            s1 = '{:,}'.format(regi) + "/" + '{:,}'.format(len(n_regs))
            s2 = '{:,}'.format(skipped) + "/" + '{:,}'.format(skippedQ)
            s3 = '{:,}'.format(numW)
            s4 = '{:,}'.format(numI)
            msg = "At region %s (%s) - %s waters, %s ions so far, eta: %s" % (s1, s2, s3, s4, ts)
            if task :
                task.updateStatus( msg )
                status ( msg )
            else :
                print " -", msg

        regi += 1

        if useQ :
            # check Q-score of pt
            ctr = [P[0], P[1], P[2]]
            qs = qscores.QscorePt ( ctr, xfI, segMap, sigQ, allAtTree=None, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
            if qs < minQ :
                skippedQ += 1
                continue
            if hMapA :
                qsA = qscores.QscorePt ( ctr, xfI, hMapA, sigQ, allAtTree=None, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
                if qsA < minQ :
                    skippedQA += 1
                    continue
            if hMapB :
                qsB = qscores.QscorePt ( ctr, xfI, hMapB, sigQ, allAtTree=None, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
                if qsB < minQ :
                    skippedQB += 1
                    continue


        msg, msgFull, atName, resName, closestChainId, clr = GuessAtom ( mol, P, atGrid=atGrid, nearAtMap=nearAtMap, doMsg=False, minDistI=1.8, maxDistI=2.5, minDistW=2.5, maxDistW=3.4, ionType="MG" )

        if atName != None :

            if 1 and atName == 'CL' :
                continue
            if 1 and atName == "NA" :
                continue

            if closestChainId == None :
                print " - new %s atom doesn't have nearest chain" % atName
                continue

            # add atom to new molecule, to be used in checkNewAtoms list
            nres = nmol.newResidue (resName, chimera.MolResId(closestChainId, len(nmol.residues)+1))
            nat = nmol.newAtom (atName, chimera.Element(atName))
            nres.addAtom( nat )
            nat.setCoord ( P )

            atGrid.AddAtomsLocal ( [nat] )
            nearAtMap[nat] = 1

            addPts.append ( [atName, resName, clr, reg, P, closestChainId] )
            if atName == 'O' :
                numW += 1
                addW.append ( [atName, resName, clr, reg, P, closestChainId] )
            else :
                numI += 1
                addI.append ( [atName, resName, clr, reg, P, closestChainId] )


    # add ions first

    natGrid = None
    #import gridm; reload(gridm)
    natGrid = gridm.Grid ()
    natGrid.FromAtomsLocal ( [at for at in mol.atoms if not at.element.name == "H"], minDistW )
    print " - done new atoms grid"

    #toChain = chainId.lower()
    print " - adding %d ions to chain %s, skipped %d/%d regions (move/Q)" % (len(addI), toChain, skipped, skippedQ)
    print " - skipped A %d, skipped B %d" % (skippedQA, skippedQB)

    largestResIdForChain = {}
    for r in mol.residues :
        if not r.id.chainId in largestResIdForChain :
            largestResIdForChain[r.id.chainId] = r.id.position
        else :
            largestResIdForChain[r.id.chainId] = max(r.id.position, largestResIdForChain[r.id.chainId])

    newIAts = []
    for atName, resName, clr, reg, P, closestChainId in addI :

        reg.show_surface()
        if reg.surface_piece :
            if atName.upper() == "NA" :
                reg.surface_piece.color = (1,0,1,1)
            elif atName.upper() == "CL" :
                reg.surface_piece.color = (1,0,1,1)
            else :
                reg.surface_piece.color = (0,1,0,1)

        cid = "_"
        if toChain == None or len(toChain) == 0 :
            cid = closestChainId
        else :
            cid = toChain

        if not cid in largestResIdForChain : largestResIdForChain[cid] = 0
        i = largestResIdForChain[cid] + 1
        largestResIdForChain[cid] = i

        if cid == None :
            print " - at %s doesn't have closest chain" % atName
            continue

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

        natGrid.AddAtomsLocal ( [nat] )
        newIAts.append ( nat )


    # then add waters
    print " - adding %d waters to chain %s, skipped %d regions" % (len(addW), toChain, skipped)

    newWAts = []
    for atName, resName, clr, reg, P, closestChainId in addW :

        reg.show_surface()
        if reg.surface_piece :
            reg.surface_piece.color = (1,0,0,1)

        cid = "_"
        if toChain == None or len(toChain) == 0 :
            cid = closestChainId
        else :
            cid = toChain

        if not cid in largestResIdForChain : largestResIdForChain[cid] = 0
        i = largestResIdForChain[cid] + 1
        largestResIdForChain[cid] = i

        if cid == None :
            print " - at %s doesn't have closest chain" % atName
            continue

        # test if already added a water within
        nearAts = natGrid.AtsNearPtLocal ( P )
        for nearAt, v in nearAts :
            if nearAt.residue.type == "HOH" and v.length < minDistW :

                print " - skipping water, too close to alrady placed water..."
                continue

                print " %d.%s -- %.2f -- %d.%s " %  (nres.id.position, nres.type, v.length, nearAt.residue.id.position, nearAt.residue.type)
                nat.occupancy = nat.occupancy / 2.0
                nearAt.occupancy = nearAt.occupancy / 2.0


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

        natGrid.AddAtomsLocal ( [nat] )
        newWAts.append ( nat )


    return newWAts, newIAts






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
                    a.isBB, a.isSC, a.isSugar, a.isBase = False, False, False, False
                    continue

                n = a.name

                a.isBB = n=="P" or n=="O1P" or n=="O2P" or n=="OP1" or n=="OP2" or n=="O5'" or n=="C5'" or n=="O3'"
                a.isSugar = n=="C1'" or n=="C2'" or n=="O4'" or n=="O2'" or n=="C3'" or n=="C4'"
                #a.isBB = a.isBB or a.isSugar
                a.isBase = not a.isBB and not a.isSugar
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



def GetMod ( name ) :
    for m in chimera.openModels.list() :
        if m.name == name :
            return m
    return None


def GetMod_ ( name ) :
    for m in chimera.openModels.list() :
        if os.path.splitext(m.name)[0] == os.path.splitext(name)[0] :
            return m
    return None


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
