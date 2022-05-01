
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
import time

from axes import prAxes
import regions
import graph
from Segger import dev_menus, timing, seggerVersion
from CGLutil.AdaptiveTree import AdaptiveTree

from chimera.resCode import nucleic3to1
from chimera.resCode import protein3to1

devMenus = False

import qscores
#reload (qscores)

import mmcif
reload (mmcif)

import molref
reload (molref)

import molbuild
reload ( molbuild )


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


phPath = "/Users/greg/_mol/phenix-1.18.2-3874/build/bin/"


# https://android.googlesource.com/toolchain/python/+/243b47fbef58ab866ee77567f2f52affd8ec8d0f/Python-2.7.3/Demo/tkinter/ttk/treeview_multicolumn.py


class SegMod_Dialog ( chimera.baseDialog.ModelessDialog ):

    title = "SegMod: Segment-guided Modeling"
    name = "segmod"

    buttons = ( 'Pro', 'NA', 'Lig', 'Thr', "Tree", "Log", "Close")
    buttons = ( 'Pro', 'NA', 'Lig', 'Thr', "Log")
    buttons = ( "M", "S", 'Mod', 'Thr', "Log", "Close")

    help = 'https://github.com/gregdp/segger'

    def fillInUI(self, parent):

        self.group_mouse_mode = None

        tw = parent.winfo_toplevel()
        self.toplevel_widget = tw
        tw.withdraw()

        self.parent = parent


        menubar = Tkinter.Menu(parent, type = 'menubar', tearoff = False)
        tw.config(menu = menubar)

        file_menu_entries = (
            ('Open ...', self.OpenModel),
            ('Save ...', self.SaveModel)
            )
        fmenu = Hybrid.cascade_menu(menubar, 'File', file_menu_entries)

        from chimera.tkgui import aquaMenuBar
        aquaMenuBar(menubar, parent, row = 0, columnspan=3)


        parent.columnconfigure(0, weight = 1)
        #parent.columnconfigure(1, weight = 1)

        row = 1

        #menubar = Tkinter.Menu(parent, type = 'menubar', tearoff = False)
        #tw.config(menu = menubar)

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
            self.mb.grid (column=1, row=0, sticky='we', padx=1)
            self.mb.menu  =  Tkinter.Menu ( self.mb, tearoff=0, postcommand=self.MapMenu )
            self.mb["menu"]  =  self.mb.menu

            ff.columnconfigure(1, weight=1)

            self.cur_dmap = None


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


            #b = Tkinter.Button(ff, text="Ca", command=self.CaBlam)
            #b.grid (column=4, row=0, sticky='w', padx=1)


        if 1 :

            row += 1

            cp = Hybrid.Popup_Panel(parent)
            cpf = cp.frame
            cpf.grid(row = row, column = 0, sticky = 'news')
            cpf.grid_remove()
            #cpf.columnconfigure(0, weight=1)
            cpf.rowconfigure(0, weight=1)
            cpf.columnconfigure(0, weight=1)
            self.treePanel = cp.panel_shown_variable
            self.treePanel.set(True)

            orow = 0

            ff = Tkinter.Frame(cpf)
            ff.grid(column=0, row=0, sticky='news')

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


            ff = Tkinter.Frame(cpf)
            ff.grid(column=0, row=1, sticky='w')

            b = Tkinter.Button(ff, text="Refresh", command=self.RefreshTree)
            b.grid (column=0, row=1, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="Select", command=self.SelectSel)
            b.grid (column=1, row=1, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="All", command=self.SelectAll)
            b.grid (column=2, row=1, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="Show", command=self.ShowSel)
            b.grid (column=3, row=1, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="Hide", command=self.HideSel)
            b.grid (column=4, row=1, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="Only", command=self.ShowSelOnly)
            b.grid (column=5, row=1, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="All", command=self.ShowAll)
            b.grid (column=6, row=1, sticky='w', padx=0, pady=1)


            #b = Tkinter.Button(ff, text="Avg", command=self.Average)
            #b.grid (column=1, row=0, sticky='w', padx=0, pady=1)

            #b = Tkinter.Button(ff, text="Open", command=self.Open)
            #b.grid (column=2, row=0, sticky='w', padx=0, pady=1)


        if 1 :
            row += 1
            dummyFrame = Tkinter.Frame(parent, relief='groove', borderwidth=1)
            Tkinter.Frame(dummyFrame).pack()
            dummyFrame.grid(row=row,column=0,columnspan=3, pady=2, sticky='we')


        row += 1
        cp = Hybrid.Popup_Panel(parent)
        cpf = cp.frame
        cpf.grid(row = row, column = 0, sticky = 'w')
        cpf.grid_remove()
        #cpf.columnconfigure(0, weight=1)
        #cpf.rowconfigure(0, weight=1)
        #cpf.columnconfigure(0, weight=1)
        self.modPanel = cp.panel_shown_variable
        self.modPanel.set(True)

        orow = 0

        if 1 :

            orow += 1
            ff = Tkinter.Frame(cpf)
            ff.grid(column=0, row=orow, sticky='w')

            b = Tkinter.Button(ff, text="Zone", command=self.Zone)
            b.grid (column=6, row=0, sticky='w', padx=1, pady=1)

            self.zoneRad = Tkinter.StringVar(ff)
            self.zoneRad.set ( "2" )
            e = Tkinter.Entry(ff, width=2, textvariable=self.zoneRad)
            e.grid(column=7, row=0, sticky='w', padx=1, pady=1)

            b = Tkinter.Label(ff, text="A ->")
            b.grid (column=8, row=0, sticky='w', padx=0, pady=1)

            self.zoneMapName = Tkinter.StringVar(ff)
            self.zoneMapName.set ( "" )
            e = Tkinter.Entry(ff, width=25, textvariable=self.zoneMapName)
            e.grid(column=9, row=0, sticky='w', padx=1, pady=1)


        if 1 :

            orow += 1
            ff = Tkinter.Frame(cpf)
            ff.grid(column=0, row=orow, sticky='w')

            b = Tkinter.Label(ff, text=" Selected")
            b.grid (column=1, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="To Map", command=self.SelRegsToMap)
            b.grid (column=2, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="", command=self.HideASel)
            b.grid (column=3, row=0, sticky='w', padx=1)

            b = Tkinter.Label(ff, text=" Atoms")
            b.grid (column=4, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="Near", command=self.ShowAtsNearSelRegs)
            b.grid (column=5, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="In", command=self.ShowAtsInSelRegs)
            b.grid (column=6, row=0, sticky='w', padx=1)

            b = Tkinter.Label(ff, text=" Regions")
            b.grid (column=7, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="Near", command=self.ShowRegsNrAts)
            b.grid (column=8, row=0, sticky='w', padx=1)


        if 1 :

            orow += 1
            ff = Tkinter.Frame(cpf)
            ff.grid(column=0, row=orow, sticky='w')

            b = Tkinter.Label(ff, text=" Unmodeled Regions")
            b.grid (column=1, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="Find", command=self.ShowDiffRegs)
            b.grid (column=5, row=0, sticky='w', padx=1)

            if devMenus :
                b = Tkinter.Button(ff, text="<", command=self.ShowLastDiffRegs)
                b.grid (column=6, row=0, sticky='w', padx=0)

            b = Tkinter.Button(ff, text="Size >", command=self.RegsSize)
            b.grid (column=7, row=0, sticky='w', padx=1)

            self.regSize = Tkinter.StringVar(ff)
            self.regSize.set ( "20" )
            e = Tkinter.Entry(ff, width=2, textvariable=self.regSize)
            e.grid(column=8, row=0, sticky='w', padx=1 )


        if 1 :

            orow += 1
            ff = Tkinter.Frame(cpf)
            ff.grid(column=0, row=orow, sticky='w')

            b = Tkinter.Label(ff, text=" To Chain")
            b.grid (column=1, row=0, sticky='w', padx=0, pady=1)

            self.addToChain = Tkinter.StringVar(ff)
            self.addToChain.set ( "" )
            e = Tkinter.Entry(ff, width=5, textvariable=self.addToChain)
            e.grid(column=2, row=0, sticky='w', padx=1, pady=1)

            b = Tkinter.Button(ff, text="sel", command=self.AddSelRes)
            b.grid (column=8, row=0, sticky='w', padx=1)

            um = Hybrid.Checkbutton(ff, 'rename', False)
            um.button.grid(column = 9, row=0, sticky = 'w', padx=1)
            self.renameAdd = um.variable

            um = Hybrid.Checkbutton(ff, 'at end', True)
            um.button.grid(column = 10, row=0, sticky = 'w', padx=1)
            self.addAtEnd = um.variable

            #b = Tkinter.Button(ff, text="dif", command=self.DiffSelRes)
            #b.grid (column=11, row=0, sticky='w', padx=1)

            if devMenus :
                b = Tkinter.Button(ff, text="RC", command=self.RenameChains)
                b.grid (column=11, row=0, sticky='w', padx=1)

            if 0 :
                b = Tkinter.Button(ff, text="M", command=self.ModSel)
                b.grid (column=12, row=0, sticky='w', padx=1)

                b = Tkinter.Button(ff, text="S", command=self.SegToggle)
                b.grid (column=13, row=0, sticky='w', padx=1)

        if devMenus :  # protein

            orow += 1
            ff = Tkinter.Frame(cpf)
            ff.grid(column=0, row=orow, sticky='w')

            l = Tkinter.Label(ff, text=' Protein:' )
            l.grid(column=0, row=0, sticky='w')

            self.addResName = Tkinter.StringVar(ff)
            self.addResName.set ( "C" )
            e = Tkinter.Entry(ff, width=5, textvariable=self.addResName)
            e.grid(column=1, row=0, sticky='w', padx=1, pady=1)

            #b = Tkinter.Button(ff, text="+", command=self.AddRes)
            #b.grid (column=2, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="set", command=self.SetRes)
            b.grid (column=2, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="-", command=self.AddLoop)
            b.grid (column=3, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="~", command=self.AddHelix)
            b.grid (column=4, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="=", command=self.AddSheet)
            b.grid (column=5, row=0, sticky='w', padx=1)

            self.numStrands = Tkinter.StringVar(ff)
            self.numStrands.set ( "1" )
            e = Tkinter.Entry(ff, width=2, textvariable=self.numStrands)
            e.grid(column=6, row=0, sticky='w', padx=1, pady=1)

            b = Tkinter.Button(ff, text="c-n", command=self.ConnectResCN)
            b.grid (column=10, row=0, sticky='w', padx=1)

            #b = Tkinter.Button(ff, text="#", command=self.RenumberRes)
            #b.grid (column=11, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="R", command=self.PutResRota)
            b.grid (column=12, row=0, sticky='w', padx=1)

            #b = Tkinter.Button(ff, text="S", command=self.AddSheet)
            #b.grid (column=5, row=0, sticky='w', padx=1)

            #b = Tkinter.Button(ff, text="Ca", command=self.CaBlam)
            #b.grid (column=10, row=0, sticky='w', padx=1)



        if devMenus :
            orow += 1

            ff = Tkinter.Frame(cpf)
            ff.grid(column=0, row=orow, sticky='w')

            #b = Tkinter.Button(ff, text="c", command=self.ConnectRes)
            #b.grid (column=13, row=0, sticky='w', padx=2)

            l = Tkinter.Label(ff, text=' Nucleic:' )
            l.grid(column=20, row=0, sticky='w')

            self.addNucName = Tkinter.StringVar(ff)
            self.addNucName.set ( "g" )
            e = Tkinter.Entry(ff, width=5, textvariable=self.addNucName)
            e.grid(column=21, row=0, sticky='w', padx=1, pady=1)

            b = Tkinter.Button(ff, text="+", command=self.AddNA)
            b.grid (column=22, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="c", command=self.ConnectRes)
            b.grid (column=23, row=0, sticky='w', padx=2)

            b = Tkinter.Button(ff, text="Guess", command=self.NaGuess)
            b.grid (column=24, row=0, sticky='w', padx=2)


        if 1 :
            orow += 1
            ff = Tkinter.Frame(cpf)
            ff.grid(column=0, row=orow, sticky='w')

            l = Tkinter.Label(ff, text=' Ligand:' )
            l.grid(column=0, row=0, sticky='w')

            self.addMolName = Tkinter.StringVar(ff)
            #self.addRess.set ( "vsgtngtkrf" )
            self.addMolName.set ( "PTQ" )
            e = Tkinter.Entry(ff, width=5, textvariable=self.addMolName)
            e.grid(column=1, row=0, sticky='w', padx=1, pady=1)

            b = Tkinter.Button(ff, text="+", command=self.AddLigand)
            b.grid (column=2, row=0, sticky='w', padx=1)



        if devMenus :
            orow += 1
            ff = Tkinter.Frame(cpf)
            ff.grid(column=0, row=orow, sticky='w')

            #l = Tkinter.Label(f, text=' Fit by:')
            #l.grid(column=0, row=0, sticky='w')

            l = Tkinter.Label(ff, text=' SegFit:')
            l.grid(column=1, row=0, sticky='w')


            self.rotaSearch = Tkinter.IntVar()
            self.rotaSearch.set ( "rota" )

            #l = Tkinter.Label(f, text=' ', width=5)
            #l.grid(column=0, row=0, sticky='w')

            c = Tkinter.Radiobutton(ff, text="PCA", variable=self.rotaSearch, value = "pca")
            c.grid (column=2, row = 0, sticky='w')

            #l = Tkinter.Label(f, text=' ', width=5)
            #l.grid(column=0, row=0, sticky='w')

            c = Tkinter.Radiobutton(ff, text="Ctr+", variable=self.rotaSearch, value = "rota")
            c.grid (column=3, row = 0, sticky='w')

            self.rotaSearchNum = Tkinter.StringVar(ff, "100")
            e = Tkinter.Entry(ff, width=3, textvariable=self.rotaSearchNum)
            e.grid(column=4, row=0, sticky='w', padx=1)

            l = Tkinter.Label(ff, text='rotations')
            l.grid(column=5, row=0, sticky='w')

            b = Tkinter.Button(ff, text="Res", command=self.SegFitSel)
            b.grid (column=6, row=0, sticky='w', padx=1)




        if 1 :
            orow += 1
            ff = Tkinter.Frame(cpf)
            ff.grid(column=0, row=orow, sticky='w')

            l = Tkinter.Label(ff, text=' Torsion - random step size')
            l.grid(column=1, row=0, sticky='w')

            self.randSearchSize = Tkinter.StringVar(ff, "10")
            e = Tkinter.Entry(ff, width=5, textvariable=self.randSearchSize)
            e.grid(column=2, row=0, sticky='w', padx=1)

            #l = Tkinter.Label(ff, text='step')
            #l.grid(column=3, row=0, sticky='w')

            b = Tkinter.Button(ff, text="Fit", command=self.TorFitRes)
            b.grid (column=8, row=0, sticky='w', padx=1)

            #b = Tkinter.Button(ff, text="Bonds", command=self.TorFitSel)
            #b.grid (column=9, row=0, sticky='w', padx=1)

            #b = Tkinter.Button(ff, text="BB", command=self.TorFitBB)
            #b.grid (column=10, row=0, sticky='w', padx=1)


        if 1 :

            orow += 1
            ff = Tkinter.Frame(cpf)
            ff.grid(column=0, row=orow, sticky='w')

            l = Tkinter.Label(ff, text=' Torsion - gradient')
            l.grid(column=1, row=0, sticky='w')

            b = Tkinter.Button(ff, text="Fit", command=self.TorFitGrads)
            b.grid (column=2, row=0, sticky='w', padx=1)

            #b = Tkinter.Button(ff, text="Ats", command=self.TorFitGradsSel)
            #b.grid (column=11, row=0, sticky='w', padx=1)

            #b = Tkinter.Button(ff, text="BB", command=self.TorFitGradsBB)
            #b.grid (column=12, row=0, sticky='w', padx=1)

            l = Tkinter.Label(ff, text='  bond rotate')
            l.grid(column=9, row=0, sticky='w')

            b = Tkinter.Button(ff, text="-", command=self.RotBondL)
            b.grid (column=10, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="+", command=self.RotBondR)
            b.grid (column=11, row=0, sticky='w', padx=1)


        if 0 :

            row += 1
            ff = Tkinter.Frame(parent)
            ff.grid(column=0, row=row, sticky='w')

            l = Tkinter.Label(ff, text='Other:')
            l.grid(column=1, row=0, sticky='w')

            b = Tkinter.Button(ff, text="Ca", command=self.CaBlam)
            b.grid (column=10, row=0, sticky='w', padx=1)


        if devMenus :

            row += 1
            ff = Tkinter.Frame(parent)
            ff.grid(column=0, row=row, sticky='w')

            b = Tkinter.Label(ff, text=" Ref:")
            b.grid (column=20, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="<", command=self.RefBack)
            b.grid (column=21, row=0, sticky='w', padx=2)

            b = Tkinter.Button(ff, text="||", command=self.StopRef)
            b.grid (column=22, row=0, sticky='w', padx=2)

            b = Tkinter.Button(ff, text="1", command=self.RefStep1)
            b.grid (column=23, row=0, sticky='w', padx=2)

            b = Tkinter.Button(ff, text=">", command=self.StartRef)
            b.grid (column=24, row=0, sticky='w', padx=2)

            b = Tkinter.Label(ff, text="MapF:")
            b.grid (column=25, row=0, sticky='w', padx=0, pady=1)

            self.mapF = Tkinter.StringVar(ff)
            #self.addRess.set ( "vsgtngtkrf" )
            self.mapF.set ( "0.1" )
            e = Tkinter.Entry(ff, width=3, textvariable=self.mapF)
            e.grid(column=26, row=0, sticky='w', padx=2, pady=1)


        row += 1
        dummyFrame = Tkinter.Frame(parent, relief='groove', borderwidth=1)
        Tkinter.Frame(dummyFrame).pack()
        dummyFrame.grid(row=row,column=0,columnspan=7, pady=2, sticky='we')


        row += 1
        global msg
        msg = Tkinter.Label(parent, width = 30, anchor = 'w', justify = 'left', fg="red")
        msg.grid(column=0, row=row, sticky='ew')
        self.msg = msg
        msg.configure(text = "Welcome to segmod... press Help button below for info")


        self.SelectedMgId = None
        self.SetVisMap ()
        self.SetVisMol ()



        #callbacks = (self.mouse_down_cb, self.mouse_drag_cb, self.mouse_up_cb)
        ##callbacks = (self.mouse_down_cb)
        #from chimera import mousemodes
        #mousemodes.addFunction('mark swim', callbacks, self.mouse_mode_icon())

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

    def Tree ( self ) :
        self.treePanel.set ( not self.treePanel.get() )

    def Mod ( self ) :
        self.modPanel.set( not self.modPanel.get() )
        pass

    def Thr ( self ) :

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



    def bind_placement_button_cb_segmod (self) :

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



    def At ( self, at, d ) :
        rt = at.residue.type
        #if rt in protein3to1 :
        #    rt = protein3to1[rt]
        #elif rt in nucleic3to1 :
        #    rt = nucleic3to1[rt]

        return " %.1fA to atom %s (element %s) in residue %s  %d, chain %s" % (d, at.name, at.element.name, rt, at.residue.id.position, at.residue.id.chainId)



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
            self.zoneMapName.set ( os.path.splitext ( dmap.name )[0] + "_z2.mrc" )


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

        self.zoneMapName.set ( os.path.splitext ( dmap.name )[0] + "_z2.mrc" )



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
            #SetBBAts ( mol )
            print "Mol: %s" % self.cur_mol.name

            self.RefreshTree()


    def StrucMenu ( self ) :
        self.strucMB.menu.delete ( 0, 'end' )   # Clear menu
        mlist = chimera.openModels.list(modelTypes = [chimera.Molecule])
        for m in mlist :
            self.strucMB.menu.add_radiobutton ( label=m.name+" (%d)"%m.id, variable=self.struc,
                                           command=lambda m=m: self.StrucSelected(m) )

        self.strucMB.menu.add_radiobutton ( label="<new>", variable=self.struc,
                                       command=lambda m=None: self.StrucSelected(None) )

    def StrucSelected ( self, mol ) :

        self.cur_mol = mol
        if mol :
            print "Selected ", mol.name, " - ", mol.id
            mlist = chimera.openModels.list(modelTypes = [chimera.Molecule])

            for m in mlist :
                m.display = False

            mol.display = True
            #SetBBAts ( mol )
            #print "Mol: %s" % self.cur_mol.name

        self.RefreshTree()


    def select_mg_cb (self, event):

        #print "Sel:", self.tree.selection()
        #print "Focus:", self.tree.focus()

        to = self.tree.focus()

        if to in self.toChain :
            print " -- Chain:", self.toChain[to]
            self.addToChain.set( self.toChain[to] )

        elif to in self.toRes :
            res = self.toRes[to]
            #try :
            #    print " -- Res: %d.%s.%s" % (res.id.position, res.type, res.id.chainId)
            #except :
            #    pass

            self.addToChain.set( res.id.chainId )

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
            else :
                at.display = True
                self.ColorAt ( at )

        for bond in self.cur_mol.bonds :
            bond.display = bond.Smart


    def HideASel ( self ) :

        for a in chimera.selection.currentAtoms () :
            res = a.residue
            res.ribbonDisplay = False
            for at in res.atoms :
                at.display = False



    def HideSel ( self ) :
        print " - showing sel..."

        if self.cur_mol == None :
            umsg ( "Select a molecule first" )
            return

        SetBBAts ( self.cur_mol )

        for m in chimera.openModels.list () :
            if type(m) == chimera.Molecule and m != self.cur_mol :
                m.display = False

        ats = self.GetSelAtoms ()
        umsg ( "Hiding %d atoms" % len(self.cur_mol.atoms) )

        for at in ats :
            r = at.residue
            if r.isProt or r.isNA :
                r.ribbonDisplay = False
                for at in r.atoms :
                    at.display = False
            else :
                at.display = False

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

        if self.cur_mol == None :
            return

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

                if hasattr ( self.cur_mol, 'chainDescr' ) and ci in self.cur_mol.chainDescr :
                    label += " -- " + ", ".join ( self.cur_mol.chainDescr[ci] )

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



    def GetSelRegsAt ( self ) :

        selAt, selRegs = None, []
        import _surface
        import _molecule
        for c in chimera.selection.currentContents()[0] :
            if type(c) == _surface.SurfacePiece :
                #print " - sp",
                selSp = c
                if hasattr ( selSp, 'region' ) :
                    selRegs.append ( selSp.region )
                    #selReg = selSp.region
                    #print " - reg: %d" % selSp.region.rid
                else :
                    print "?"
            elif type(c) == _molecule.Atom :
                selAt = c
                #print " - atom: %s" % selAt.name

        return selRegs, selAt


    def GetSelRegsAts ( self ) :

        selAts, selRegs = [], []
        import _surface
        import _molecule
        for c in chimera.selection.currentContents()[0] :
            if type(c) == _surface.SurfacePiece :
                print " - sp",
                selSp = c
                if hasattr ( selSp, 'region' ) :
                    selRegs.append ( selSp.region )
                    #selReg = selSp.region
                    print " - reg: %d" % selSp.region.rid
                else :
                    print "?"
            elif type(c) == _molecule.Atom :
                selAts.append ( c )
                #print " - atom: %s" % selAt.name

        return selRegs, selAts




    def ShowAtsNearSelRegs ( self ) :

        selRegs, selAt = self.GetSelRegsAt()

        ats = None
        if len(selRegs) > 0 :
            umsg ( "Finding atoms near selected regions..." )
            ats = molref.AtomsNearRegs ( selRegs, 4.0 )

        elif len(chimera.selection.currentAtoms()) > 0 :
            ats = chimera.selection.currentAtoms()
            umsg ( "Finding atoms near %d selected atoms..." % len(ats) )
            ats = molref.AtomsNearAtoms ( ats, 4.0 )

        if ats == None or len(ats) == 0 :
            umsg ( "Nothing selected or no atoms close" )
            return

        mols = {}
        ress = {}

        for at in ats :
            mols[at.molecule] = 1
            ress[at.residue] = 1

        umsg ( "Sel regions near %d atoms, %d residue, %d models" % (len(ats), len(ress), len(mols)) )

        for mol in mols :
            SetBBAts ( mol )

        for res in ress :
            res.ribbonDisplay = True
            for at in res.atoms :
                if at.isSC or at.name == "CA" :
                    at.display = True
                    try :
                        at.color = atomColors[at.element.name.upper()]
                    except :
                        at.color = atomColors[' ']
                #else :
                #    at.display = False


    def ShowAtsInSelRegs ( self ) :

        regs, selAt = self.GetSelRegsAt()

        if len(regs) == 0 :
            umsg ( "No regions selected?" )
            return

        dmap = regs[0].segmentation.seg_map
        rdata = molbuild.RegsData  ( regs )
        rmat = rdata.matrix()

        #zoneR = numpy.min(dmap.data.step)
        #rpoints = numpy.concatenate ( [r.map_points() for r in regs], axis=0 ).astype ( numpy.float32 )
        #rdata = VolumeData.zone_masked_grid_data ( dmap.data, rpoints, zoneR )
        #rmat = rdata.matrix()
        #gdata = VolumeData.Array_Grid_Data ( ndata.full_matrix(), dmap.data.origin, dmap.data.step, dmap.data.cell_angles, name = "atom masked" )
        #nv = VolumeViewer.volume.volume_from_grid_data ( rdata )
        #nv.name = "regs"
        apoints = numpy.zeros ( [1, 3], numpy.float32 )

        mols = []
        #print "atoms near in:"
        for m in chimera.openModels.list() :
            if type(m) == chimera.Molecule and m.display == True :

                print m.name

                for r in m.residues :

                    show = False
                    for at in r.atoms :

                        apoints[0] = dmap.openState.xform.inverse().apply ( at.xformCoord() )
                        #apoints[0] = at.coord()

                        values, outside = VolumeData.interpolate_volume_data ( apoints, rdata.xyz_to_ijk_transform, rmat )

                        if values[0] > 0 :
                            #print values
                            #break
                            show = True

                    for at in r.atoms :
                        at.display = show
                    r.ribbonDisplay = show





    def RegsNearAts ( self, ats, seg, segMap, task=None ) :

        #points = _multiscale.get_atom_coordinates ( ats, transformed = True )
        #print " - search tree: %d atoms" % ( len(ats) )
        #plist = points.tolist()
        #atsTree = AdaptiveTree ( plist, plist, 2.0)

        import gridm; reload ( gridm )
        agridm = gridm.Grid ()
        agridm.FromAtoms ( ats, segMap.data.step[0] )

        regs = []
        import _contour
        ri = 0
        for reg in seg.regions :

            points = reg.points().astype ( numpy.float32 )
            _contour.affine_transform_vertices ( points, segMap.data.ijk_to_xyz_transform )
            _contour.affine_transform_vertices ( points, Matrix.xform_matrix( segMap.openState.xform ) )

            nearAt = False
            for pt in points :
                nats = agridm.AtsNearPt ( chimera.Point(pt[0], pt[1], pt[2]) )
                if len(nats) > 0 :
                    nearAt = True
                    break

            if nearAt :
                regs.append ( reg )

            if ri % 1000 == 0 :
                if task != None :
                    task.updateStatus ( "Finding regions in %s, at %d/%d" % ( seg.name, ri+1, len(seg.regions) ) )
                else :
                    status ( "Seg: %s, %d/%d regions" % ( seg.name, ri+1, len(seg.regions) ) )
                    print ".",
            ri += 1

        return regs


    def ShowRegsNrAts ( self ) :
        print "show regs near ats"

        ats = chimera.selection.currentAtoms ()
        ats = [at for at in ats if not at.element.name == "H"]

        if len(ats) == 0 :
            umsg ( "No atoms selected?" )
            return

        seg = current_segmentation()
        segMap = segmentation_map ()

        if seg == None or segMap == None :
            umsg ( "No segmentation in Segger dialog?" )
            return

        umsg ( "Seg: %s, %d regions" % ( seg.name, len(seg.regions) ) )

        chimera.selection.clearCurrent ()
        for reg in seg.regions :
            reg.hide_surface()

        regs = self.RegsNearAts ( ats, seg, segMap, task=None )
        for reg in regs :
            reg.show_surface()
            if hasattr ( reg, 'surface_piece' ) and reg.surface_piece :
                chimera.selection.addCurrent ( reg.surface_piece )

        umsg ( "Selected %d/%d regions near %d atoms" % (len(regs), len(seg.regions), len(ats)) )


    def ShowDiffRegs ( self ) :

        print "diff regs"

        seg = current_segmentation()
        segMap = segmentation_map ()

        if seg == None or segMap == None :
            umsg ( "Segment first in Segger dialog" )
            return

        if self.cur_mol == None :
            umsg ( "Select a molecule" )
            return

        ats = chimera.selection.currentAtoms ()

        if 0 :
            ats = []
            for m in chimera.openModels.list() :
                if type(m) == chimera.Molecule and m.display == True :
                    ats.extend ( m.atoms )

        ats = self.cur_mol.atoms
        ats = [at for at in ats if not at.element.name == "H"]

        if len(ats) == 0 :
            umsg ( "No atoms selected?" )
            return


        umsg ( "Seg: %s, %d regions" % ( seg.name, len(seg.regions) ) )

        from chimera import tasks, CancelOperation
        task = tasks.Task('Finding unmodeled regions', modal = True)

        try :
            regs = self.RegsNearAts ( ats, seg, segMap, task=task )
        except CancelOperation:
            print "canceled"
        finally :
            task.finished()
            print "done"


        regmap = {}
        for reg in regs :
            regmap[reg] = 1

        self.diffRegs = []
        for reg in seg.regions :
            if not reg in regmap :
                self.diffRegs.append ( reg )

        self.ShowLastDiffRegs ()

        umsg ( "Selected %d/%d regions not near %d atoms" % (len(self.diffRegs), len(seg.regions), len(ats)) )


    def ShowLastDiffRegs ( self ) :

        seg = current_segmentation()
        segMap = segmentation_map ()

        if seg == None or segMap == None :
            umsg ( "No segmentation in Segger dialog?" )
            return

        if not hasattr ( self, 'diffRegs' ) :
            umsg ( "No difference regions found, press 'Diff' first" )
            return

        for reg in seg.regions :
            reg.hide_surface()

        chimera.selection.clearCurrent ()
        for reg in self.diffRegs :
            if hasattr ( reg, 'surface_piece' ) and reg.surface_piece :
                if len(reg.points()) > int ( self.regSize.get() ) :
                    reg.show_surface()
                    reg.surface_piece.color = ( 1.0, 0.0, 0.0, 1.0 )
                #else :
                #    reg.surface_piece.color = ( .7, .7, .7, 1.0 )


        #chimera.selection.addCurrent ( reg.surface_piece )

    def RegsSize ( self ) :

        self.ShowLastDiffRegs ()
        return

        selRegs, selAt = self.GetSelRegsAt()

        msg = "sizes:"
        for reg in selRegs[0:10] :
            msg += " %d" % len( reg.points() )

        umsg ( msg )

        seg = current_segmentation()
        segMap = segmentation_map ()

        for reg in seg.regions :

            if hasattr ( reg, 'surface_piece' ) and reg.surface_piece :
                if reg.surface_piece.display :
                    if len(reg.points()) > int ( self.regSize.get() ) :
                        reg.surface_piece.color = ( 1.0, 0.0, 0.0, 1.0 )
                    else :
                        reg.surface_piece.color = ( .7, .7, .7, 1.0 )





    def AddSelRes ( self ) :

        if self.cur_mol == None :
            self.cur_mol = chimera.Molecule()
            self.cur_mol.name = "new"
            #umsg ("Select a molecule first")
            #return []
            chimera.openModels.add ( [self.cur_mol] )
            self.struc.set ( self.cur_mol.name + " (%d)" % self.cur_mol.id )
            self.RefreshTree()

        toMol = self.cur_mol
        #chainId = self.chain.get()
        toChain = self.addToChain.get().strip().replace(" ", "")

        if len(toChain) == 0 :
            #umsg ( "enter a chain Id" )
            self.addToChain.set("A")
            toChain = self.addToChain.get()

        asType = None
        if self.renameAdd.get() :
            molToAdd = self.addMolName.get().upper().strip().replace(" ", "")
            if len(molToAdd) > 0 :
                asType = molToAdd

        selRes = chimera.selection.currentResidues()
        if len(selRes) == 0 :
            umsg ( "Select residue(s) to add to selected mol and chain" )
            return

        print " - sorting %d res" % len(selRes)
        resi = sorted(selRes, key=lambda r: r.id.position, reverse=False)

        print " - making rmap..."
        rmap = {}
        for r in self.cur_mol.residues :
            if r.id.chainId == toChain :
                rmap[r.id.position] = r

        print " - adding ress..."
        for res in resi :

            rid = None
            if not self.addAtEnd.get() :
                if res.id.position in rmap :
                    print " - res %d found, skipping" % res.id.position
                    rres = rmap[res.id.position]
                    if rres.type != res.type :
                        print " - sel res %s,%d != %s,%d" % (res.type, res.id.position, rres.type, rres.id.position)
                    continue
                else :
                    rid = res.id.position

            xf = res.molecule.openState.xform
            xf.premultiply ( toMol.openState.xform.inverse() )
            nres = molref.AddResToMol ( res, toMol, toChain, xf, withoutAtoms=[], rid=rid, asType=asType )


            status ( "Added res %s.%d.%s to %s chain %s position %d" % (res.type, res.id.position, res.id.chainId, toMol.name, toChain, nres.id.position) )

            rmap[nres.id.position] = nres
            if nres.type in protein3to1 :
                nres.ribbonColor = res.ribbonColor
                if nres.id.position-1 in rmap :
                    pres = rmap[nres.id.position-1]
                    if pres.type in protein3to1 :
                        nb = self.cur_mol.newBond ( pres.atomsMap['C'][0], nres.atomsMap['N'][0] )
                        nres.ribbonDisplay = True
                        #print " - connected to prev prot res %d" % pres.id.position
                        for at in nres.atoms :
                            at.display = False
                if nres.id.position+1 in rmap :
                    fres = rmap[nres.id.position+1]
                    if fres.type in protein3to1 :
                        nb = self.cur_mol.newBond ( nres.atomsMap['C'][0], fres.atomsMap['N'][0] )
                        nres.ribbonDisplay = True
                        #print " - connected to next prot res %d" % fres.id.position
                        for at in nres.atoms :
                            at.display = False
            elif nres.type in nucleic3to1 :
                nres.ribbonColor = res.ribbonColor
                if nres.id.position-1 in rmap :
                    pres = rmap[nres.id.position-1]
                    if pres.type in nucleic3to1 :
                        nb = self.cur_mol.newBond ( pres.atomsMap["O3'"][0], nres.atomsMap['P'][0] )
                        nres.ribbonDisplay = True
                        #print " - connected to prev prot res %d" % pres.id.position
                        for at in nres.atoms :
                            at.display = False
                if nres.id.position+1 in rmap :
                    fres = rmap[nres.id.position+1]
                    if fres.type in nucleic3to1 :
                        nb = self.cur_mol.newBond ( nres.atomsMap['P'][0], fres.atomsMap["O3'"][0] )
                        nres.ribbonDisplay = True
                        #print " - connected to next prot res %d" % fres.id.position
                        for at in nres.atoms :
                            at.display = False


        self.RefreshTree ()

    def RenameChains ( self ) :

        print "renaming chains for %s" % self.cur_mol.name
        from _multiscale import get_atom_coordinates

        fromMol = None
        for m in chimera.openModels.list() :
            if type(m) == chimera.Molecule and m != self.cur_mol and m.display == True :
                fromMol = m
                break

        print " - from %s" % fromMol.name
        chAts = {}
        for res in fromMol.residues :
            if res.id.chainId in chAts :
                chAts[res.id.chainId].extend ( res.atoms )
            else :
                chAts[res.id.chainId] = res.atoms [:]

        print " - chains: ",
        chCtr = []
        for ch, atoms in chAts.iteritems() :
            print ch,
            points = get_atom_coordinates ( atoms, transformed = True )
            com = numpy.sum(points, axis=0) / len(points)
            chCtr.append ( [ch, chimera.Vector ( com[0], com[1], com[2] )] )

        print ""

        from random import random
        chAts = {}
        chRes = {}
        chCol = {}
        for res in self.cur_mol.residues :
            if res.id.chainId in chAts :
                chAts[res.id.chainId].extend ( res.atoms )
                chRes[res.id.chainId].append ( res )
            else :
                chAts[res.id.chainId] = res.atoms [:]
                chRes[res.id.chainId] = [res]
                chCol[res.id.chainId] = chimera.MaterialColor ( random(), random(), random(), 1.0 )

        nmol = chimera.Molecule()
        nmol.name = self.cur_mol.name + "_new_chain_names"

        atMap = {}
        for ch, atoms in chAts.iteritems() :
            points = get_atom_coordinates ( atoms, transformed = True )
            com = numpy.sum(points, axis=0) / len(points)
            comv = chimera.Vector ( com[0], com[1], com[2] )

            minD, toCh = 1e9, None
            for chi, ctr in chCtr :
                v = ctr - comv
                if v.length < minD :
                    minD, toCh = v.length, chi

            ress = chRes[ch]
            print "%s -> %s, %d residues" % (ch, toCh, len(ress))
            for res in ress :
                nres = nmol.newResidue ( res.type, chimera.MolResId(toCh, res.id.position))
                res.isProt = res.type in protein3to1
                res.isNA = res.type in nucleic3to1
                if res.isProt or res.isNA :
                    nres.ribbonDisplay = True
                    nres.ribbonColor = chCol[ch]
                for at in res.atoms :
                    nat = nmol.newAtom (at.name, chimera.Element(at.element.number))
                    atMap[at] = nat
                    nres.addAtom( nat )
                    nat.drawMode = nat.EndCap
                    nat.setCoord ( at.coord() )
                    if res.isProt or res.isNA :
                        nat.display = False
                    else :
                        nat.display = True
                        if nat.element.name.upper() in atomColors :
                            nat.color = atomColors[nat.element.name.upper()]
                        else :
                            nat.color = chimera.MaterialColor ( random(), random(), random(), 1.0 )


        for bond in self.cur_mol.bonds :
            nb = nmol.newBond ( atMap[bond.atoms[0]], atMap[bond.atoms[1]] )
            nb.display = nb.Smart
            nb.drawMode = nb.Stick

        chimera.openModels.add ( [nmol] )




    def S ( self ) :

        from Segger import regions
        for m in chimera.openModels.list() :
            if type (m) == regions.Segmentation :
                m.display = not m.display


    def M ( self ) :

        print ""
        if hasattr ( self, 'showMods' ) :
            for m in self.showMods :
                m.display = True
            del self.showMods
        else :
            self.showMods = []
            for m in chimera.openModels.list() :
                if type(m) != chimera.Molecule :
                    if m.display :
                        self.showMods.append ( m )
                        m.display = False





    def DiffSelRes0 ( self ) :

        toMol = self.cur_mol
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return []

        #chainId = self.chain.get()
        toChain = self.addToChain.get().strip().replace(" ", "")
        if len(toChain) != 1 :
            umsg ( "enter a chain Id" )
            return

        asType = None
        if self.renameAdd.get() :
            molToAdd = self.addMolName.get().upper().strip().replace(" ", "")
            if len(molToAdd) > 0 :
                asType = molToAdd

        selRes = chimera.selection.currentResidues()
        if len(selRes) == 0 :
            umsg ( "Select residue(s) to add to selected mol and chain" )
            return

        resi = sorted(selRes, key=lambda r: r.id.position, reverse=False)

        rmap = {}
        for r in self.cur_mol.residues :
            if r.id.chainId == toChain :
                rmap[r.id.position] = r

        delRes = []

        selStr = ":"

        for res in resi :
            if not res.id.position in rmap :
                delRes.append ( res )
                selStr += "%s," % res.id.position

        print " - not present:", selStr[0:-1]

        print " - different:"
        selStr = ":"
        for res in resi :
            if res.id.position in rmap :
                rres = rmap[res.id.position]
                if rres.type != res.type :
                    selStr += "%s," % res.id.position
                    print "[%s,%d - %s,%d]" % (res.type,res.id.position,rres.type,rres.id.position),

        print ""
        print selStr[0:-1]


        #for r in delRes :
        #    self.cur_mol.deleteResidue ( r )


        self.RefreshTree ()



    def DiffSelRes ( self ) :

        mols = []
        for m in chimera.openModels.list() :
            if m.display == True :
                mols.append ( m )

        m0 = mols[0]
        print "Mol: %s" % m0.name
        crmap = {}
        for r in m0.residues :
            if not r.id.chainId in crmap :
                crmap[r.id.chainId] = {}
            crmap[r.id.chainId][r.id.position] = r

        from chimera.resCode import nucleic3to1
        from chimera.resCode import protein3to1

        maxRmsd = 0.0
        rmsds = []
        for m in mols[1:] :
            print " --  %s" % m.name
            for r in m.residues :
                if r.type.upper() in nucleic3to1 or r.type.upper() in protein3to1 :
                    if r.id.chainId in crmap :
                        rmap = crmap[r.id.chainId]
                        if r.id.position in rmap :
                            r0 = rmap[r.id.position]
                            if r.type != r0.type :
                                print "[%s,%d.%s - %s,%d.%s]" % (r.type,r.id.position,r.id.chainId,r0.type,r0.id.position,r0.id.chainId)
                                r.rmsd = None
                            else :
                                r.rmsd = 0.0
                                for at in r.atoms :
                                    at0 = r0.atomsMap[at.name][0]
                                    v = at0.xformCoord() - at.xformCoord()
                                    r.rmsd += v.length
                                r.rmsd = r.rmsd / float(len(r.atoms))
                                rmsds.append ( r.rmsd )
                                if r.rmsd > maxRmsd :
                                    maxRmsd = r.rmsd
                    else :
                        r.rmsd = None
                else :
                    #print r.type
                    pass

        print " - max rmsd: %.3f" % maxRmsd
        std = 2.0 * numpy.std ( rmsds )
        std = maxRmsd
        print " - std: %.3f" % std

        for m in mols[1:] :
            for r in m.residues :
                if r.type.upper() in nucleic3to1 or r.type.upper() in protein3to1 :
                    if not hasattr ( r, 'rmsd' ) :
                        r.ribbonColor = chimera.MaterialColor ( .7, .7, .7, 1.0 )
                    elif r.rmsd == None :
                        r.ribbonColor = chimera.MaterialColor ( 1.0, 1.0, .2, 1.0 )
                    else :
                        C = min ( r.rmsd / std, 1.0 )
                        r.ribbonColor = chimera.MaterialColor ( C, 0.0, 0.0, 1.0 )



    def SelRegsToMap ( self ) :

        selRegs, selAt = self.GetSelRegsAt()
        if len(selRegs) > 0 :
            molbuild.RegsToMap ( selRegs )



    def AddLigand ( self ) :

        toMol = self.cur_mol
        if self.cur_mol == None :
            umsg ("Select a model to add to")
            return []

        dmap = self.cur_dmap

        #chainId = self.chain.get()
        toChain = self.addToChain.get().strip().replace(" ", "")
        molToAdd = self.addMolName.get().upper().strip().replace(" ", "")

        selRegs, selAt = self.GetSelRegsAt()

        if selAt :
            toMol = selAt.molecule
            if len(toChain) == 0 :
                toChain = selAt.residue.id.chainId

        msg = "Adding ligand mol: %s" % molToAdd

        if len(selRegs) > 0 :
            msg += " - in %d region(s)" % len(selRegs)

            dmap = molbuild.RegsToMap ( selRegs )
            dmap.delAfter = True

            nearAts = molref.AtomsNearRegs ( selRegs, 4.0 )
            if len(nearAts) > 0 :
                toMol = nearAts[0].molecule
                if len(toChain) == 0 :
                    toChain = nearAts[0].residue.id.chainId

        else :
            umsg ( "Segment and select some regions" )
            return

        if toMol :
            print " - to %s, chain %s" % (toMol.name , toChain)
            msg += " - to %s, chain %s" % (toMol.name , toChain)

        if self.cur_dmap :
            print " - in map %s" % dmap.name
            msg += " - in map %s" % dmap.name

        status ( msg  )

        molref.AddMol ( molToAdd, selAt, dmap, selRegs, toMol, toChain )

        if hasattr ( dmap, 'delAfter' ) :
            chimera.openModels.close ( dmap )

        self.RefreshTree ()


    def AddLoop ( self ) :
        print "loop..."

        dmap = self.cur_dmap
        if dmap == None :
            umsg ( "Select a map?" )
            return

        ress = chimera.selection.currentResidues ()
        selAt = chimera.selection.currentAtoms()

        if len(selAt) == 1 and selAt[0].name == "N" :
            seq = self.addResName.get()
            umsg ( "Adding connected [loop] ress - %s" % seq )
            molbuild.AddRess_N ( seq, selAt[0].residue )

        elif len(ress) > 0 :
            for r in ress :
                r.isHelix = False
                r.isSheet = False




    def AddHelix ( self ) :
        print "helix..."

        dmap = self.cur_dmap
        if dmap == None :
            umsg ( "Select a map?" )
            return

        ress = chimera.selection.currentResidues ()
        selAt = chimera.selection.currentAtoms()

        if len(selAt) == 1 and selAt[0].name == "N" :
            seq = self.addResName.get()
            umsg ( "Adding connected [helix] ress - %s" % seq )
            molbuild.AddRess_N ( seq, selAt[0].residue, helix=True )

        elif len(ress) > 0 :
            for r in ress :
                atN, atC, atCA, atO = r.atomsMap["N"][0], r.atomsMap["C"][0], r.atomsMap["CA"][0], r.atomsMap["O"][0]
                r.isHelix = True
                r.isSheet = False

                bPhi = atN.bondsMap[atCA]
                atC_ = None
                for at, b in atN.bondsMap.iteritems() :
                    if at.name == "C" and at.residue.id.position == atN.residue.id.position-1 :
                        atC_ = at
                if atC_ != None :
                    d = molref.diha ( atC_, atN, atCA, atC )
                    print " - res %d phi %.2f" % (r.id.position, d )
                    atsF, atsB = molref.BondAts ( atN, atCA )
                    if len(atsF) < len(atsB) :
                        #print "<"
                        molref.RotBond ( atN, atCA, atsF, -57 - d )
                    else :
                        #print ">"
                        molref.RotBond ( atN, atCA, atsB, - (-57 - d) )

                d = molref.diha ( atN, atCA, atC, atO )
                print " - res %d psi %.2f" % (r.id.position, d )
                atsF, atsB = molref.BondAts ( atCA, atC )
                if len(atsF) < len(atsB) :
                    print "<"
                    molref.RotBond ( atCA, atC, atsF, 133 - d ) # -47 + 180 = 133
                else :
                    print ">"
                    molref.RotBond ( atCA, atC, atsB, - (133 - d) )

        else :
            toMol = self.cur_mol
            if self.cur_mol == None :
                #umsg ("Select a molecule first")
                #return []
                toMol = chimera.Molecule()
                toMol.name = "new"
                chimera.openModels.add ( [toMol] )

            #chainId = self.chain.get()
            toChain = self.addToChain.get().strip().replace(" ", "")

            if len(toChain) == 0 :
                # find a new chain id that's not used...
                cids = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
                for r in toMol.residues :
                    if r.id.chainId in cids :
                        cids = cids.replace ( r.id.chainId, "" )
                toChain = cids[0]
                self.addToChain.set ( toChain )

            selRegs, selAt = self.GetSelRegsAt()

            if len(selRegs) == 0 :
                umsg ( "Select some regions..." )
                return

            surfMod = GetSegMod ()

            hx = molbuild.Helix ()
            hx.Make ( selRegs, toMol, toChain, dmap, surfMod )




    def AddSheet ( self ) :
        print "sheet..."

        dmap = self.cur_dmap
        if dmap == None :
            umsg ( "Select a map?" )
            return

        ress = chimera.selection.currentResidues ()
        selAt = chimera.selection.currentAtoms()

        for r in ress :
            r.isHelix = False
            r.isSheet = True


    def SetRes ( self ) :

        ress = chimera.selection.currentResidues()
        seq = self.addResName.get()

        toResI = None
        try :
            toResI = int(seq)
        except :
            pass

        if toResI != None :
            if len(ress) == 1 :
                umsg ( "Renumbering residues at %d -> %d" % (ress[0].id.position, toResI) )
                molbuild.SetResI ( ress[0], toResI )
                return
            else :
                umsg ( "Select (one) residue to start numbering at" )
                return

        umsg ( "Setting %d res - seq: %s" % (len(ress), seq) )

    	from SwapRes import swap, SwapResError
        from chimera.resCode import protein1to3

        for i, r in enumerate (ress) :
            if seq[i].upper() in protein1to3 :

                swap(r, protein1to3[seq[i]], preserve=False, bfactor=False)

                for at in r.atoms :
                    at.drawMode = at.EndCap
                    at.display = True # not showRibbon
                    at.color = atomColors[at.element.name if at.element.name in atomColors else " "]


    def ConnectResCN ( self ) :

        from chimera.resCode import protein3to1, protein1to3

        atN, atC = None, None
        for at in chimera.selection.currentAtoms () :
            if at.name == "C" and at.residue.type in protein3to1 :
                atC = at
            elif at.name == "N" and at.residue.type in protein3to1 :
                atN = at

        if atN == None or atC == None :
            umsg ( "select N and C atoms" )
            return

        molbuild.ConnectRess ( atN, atC )


    def RenumberRes ( self ) :

        print " - renumber"


    def PutResRota ( self ) :

        print "rota..."




    def AddNA ( self ) :

        toMol = self.cur_mol
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return []

        dmap = self.cur_dmap

        #chainId = self.chain.get()
        toChain = self.addToChain.get().strip().replace(" ", "")
        resToAdd = self.addNucName.get().upper().strip().replace(" ", "")

        selRegs, selAt = self.GetSelRegsAt()

        if selAt :
            toMol = selAt.molecule
            if len(toChain) == 0 :
                toChain = selAt.residue.id.chainId
        else :
            umsg ( "Select an atom..." )
            return


        #msg = "Adding %s" % molToAdd
        print "Adding", resToAdd

        if len(selRegs) > 0 :
            print " - in %d region(s)" % len(selRegs)
            dmap = RegsToMap ( selRegs )
            dmap.delAfter = True

        else :
            if self.cur_dmap :
                print " - in map %s" % dmap.name

        #if toMol :
        #    msg += " to %s chain %s" % (toMol.name , toChain)

        #umsg ( msg  )

        molref.AddNuc ( resToAdd, selAt, dmap, selRegs )

        if hasattr ( dmap, 'delAfter' ) :
            chimera.openModels.close ( dmap )

        self.RefreshTree ()






    def ConnectRes ( self ) :

        print ""
        print "Connect"

        selAts = chimera.selection.currentAtoms()

        if len(selAts) != 2 :
            umsg ( "Select two atoms" )
            return

        at1, at2 = selAts

        mol1 = at1.molecule
        mol2 = at2.molecule

        if mol1 != mol2 :
            umsg ( "Not same molecule" )
            return

        nb = mol1.newBond ( at1, at2 )
        nb.display = nb.Smart
        nb.drawMode = nb.Stick


    def DelSel ( self ) :

        print ""
        print "Del"

        selAts = chimera.selection.currentAtoms()
        if len(selAts) > 0 :
            pass

        selBonds = chimera.selection.currentBonds()
        for b in selBonds :
            b.atoms[0].molecule.deleteBond(b)




    def NaGuess ( self ) :

        print ""
        print "NaGuess"

        selAts = chimera.selection.currentAtoms()

        if len(selAts) == 0 :
            umsg ( "Select something" )
            return

        at = selAts[0]
        res = at.residue
        mol = at.molecule

        print "Guess - at %s.%d.%s" % (res.type, res.id.position, res.id.chainId)



    def RotBondL ( self ) :

        b = chimera.selection.currentBonds()
        if len(b) != 1 :
            umsg ( "select one bond" )
            return

        bond = b[0]

        atsF, atsB = molref.BondAts ( bond.atoms[0], bond.atoms[1] )
        if len(atsF) < len(atsB) :
            #print "<"
            molref.RotBond ( bond.atoms[0], bond.atoms[1], atsF, 5.0 )
        else :
            #print ">"
            molref.RotBond ( bond.atoms[0], bond.atoms[1], atsB, -5.0 )



    def RotBondR ( self ) :

        b = chimera.selection.currentBonds()
        if len(b) != 1 :
            umsg ( "select one bond" )
            return

        bond = b[0]

        atsF, atsB = molref.BondAts ( bond.atoms[0], bond.atoms[1] )
        if len(atsF) < len(atsB) :
            #print "<"
            molref.RotBond ( bond.atoms[0], bond.atoms[1], atsF, -5.0 )
        else :
            #print ">"
            molref.RotBond ( bond.atoms[0], bond.atoms[1], atsB, +5.0 )


    def TorFitRes ( self ) :

        print "TorFit - Random"

        selRegs, selAt = self.GetSelRegsAt()

        if selAt == None :
            umsg ( "Select an atom" )
            return

        dmap = self.cur_dmap
        if len(selRegs) > 0 :
            dmap = molbuild.RegsToMap ( selRegs )
            dmap.delAfter = True
            print " - in regs map: %s" % dmap.name

        stepSize = float(self.randSearchSize.get())
        print " - step size: %.2f" % stepSize

        #res = selAt.residue
        #print " - residue %s %d.%s" % (res.type, res.id.position, res.id.chainId)
        ress = chimera.selection.currentResidues()

        for r in ress :
            if r.type in protein3to1 or r.type in nucleic3to1 :
                umsg ( "Only for ligands..." )
                return

        molref.TorFitRand ( ress, dmap, stepSize, parent=self.parent )

        if 0 :
            from chimera import tasks, CancelOperation
            task = tasks.Task('Torsion Fit', modal = True)

            try :
                molref.TorFitRand ( ress, dmap, stepSize, parent=None, task=task )
            except CancelOperation:
                print "canceled"
            finally :
                task.finished()
                print "done"



        if hasattr ( dmap, 'delAfter' ) :
            chimera.openModels.close ( dmap )


    def TorFitSel ( self ) :

        selRegs, selAt = self.GetSelRegsAt()

        sbonds = chimera.selection.currentBonds()
        if len(sbonds) == 0 :
            umsg ( "Select one or more bonds" )
            return

        dmap = self.cur_dmap
        if len(selRegs) > 0 :
            dmap = RegsToMap ( selRegs )
            dmap.delAfter = True
            dmap.display = False

        stepSize = float(self.randSearchSize.get())

        print " - %d bonds sel" % len(sbonds)

        molref.TorFitRSel ( sbonds, dmap, stepSize )

        if hasattr ( dmap, 'delAfter' ) :
            chimera.openModels.close ( dmap )


    def TorFitBB ( self ) :

        selRegs, selAt = self.GetSelRegsAt()

        ress = chimera.selection.currentResidues()
        if len(ress) == 0 :
            umsg ( "Select one or more ress" )
            return

        SetBBAts ( ress[0].molecule )

        amap = {}
        for r in ress :
            for at in r.atoms :
                if at.isBB :
                    amap[at] = 1

        sbonds = []
        for bond in ress[0].molecule.bonds :
            if bond.atoms[0] in amap or bond.atoms[1] in amap :
                sbonds.append ( bond )


        dmap = self.cur_dmap
        if len(selRegs) > 0 :
            dmap = RegsToMap ( selRegs )
            dmap.delAfter = True
            dmap.display = False

        stepSize = float(self.randSearchSize.get())

        print " - %d bonds sel" % len(sbonds)

        molref.TorFitRSel ( sbonds, dmap, stepSize )

        if hasattr ( dmap, 'delAfter' ) :
            chimera.openModels.close ( dmap )





    def TorFitGrads ( self ) :

        print "TorFit - Grads"
        selRegs, selAt = self.GetSelRegsAt()

        dmap = self.cur_dmap
        delAfter = False
        if len(selRegs) > 0 :
            dmap = molbuild.RegsToMap ( selRegs )
            delAfter = True
            dmap.display = False

        ress = chimera.selection.currentResidues()

        if len(ress) == 0 :
            umsg ( "Select some residues..." )
            return

        for r in ress :
            if r.type in protein3to1 or r.type in nucleic3to1 :
                umsg ( "Only for ligands..." )
                return

        molref.TorFitGrads ( ress, dmap )

        if delAfter :
            chimera.openModels.close ( dmap )



    def TorFitGradsSel ( self ) :

        selRegs, selAts = self.GetSelRegsAts()

        dmap = self.cur_dmap
        if len(selRegs) > 0 :
            dmap = RegsToMap ( selRegs )
            dmap.delAfter = True
            dmap.display = False

        selAts = chimera.selection.currentAtoms()
        if len(selAts) == 0 :
            umsg ( "Select an atom" )
            return

        molref.TorFitGrads ( selAts[0].residue, dmap, useAtoms=selAts )

        if hasattr ( dmap, 'delAfter' ) :
            chimera.openModels.close ( dmap )



    def TorFitGradsBB ( self ) :

        selRegs, selAt = self.GetSelRegsAt()

        dmap = self.cur_dmap
        if len(selRegs) > 0 :
            dmap = RegsToMap ( selRegs )
            dmap.delAfter = True
            dmap.display = False

        ress = chimera.selection.currentResidues()
        if len(ress) == 0 :
            umsg ( "Select one or more ress" )
            return

        SetBBAts ( ress[0].molecule )

        amap = {}
        for r in ress :
            for at in r.atoms :
                if at.isBB :
                    amap[at] = 1

        sbonds = []
        for bond in ress[0].molecule.bonds :
            if bond.atoms[0] in amap or bond.atoms[1] in amap :
                sbonds.append ( bond )


        molref.TorFitGradsBonds ( sbonds, dmap, useAtoms=None )

        if hasattr ( dmap, 'delAfter' ) :
            chimera.openModels.close ( dmap )



    def SegFitSel ( self ) :

        print "fitting sel..."

        selRegs, selAts = self.GetSelRegsAts()

        dmap = self.cur_dmap
        if len(selRegs) > 0 :
            dmap = RegsToMap ( selRegs )
            dmap.delAfter = True

        if len(selAts) == None :
            umsg ( "Select some atoms to align" )
            return

        molref.SegFitRes ( selAts[0].residue, dmap, selRegs, useAts=selAts )

        if hasattr ( dmap, 'delAfter' ) :
            chimera.openModels.close ( dmap )

        self.RefreshTree ()



    def RFitRes ( self, res ) :

        print ""




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


        nname = self.zoneMapName.get()

        saveFile = True
        if len(nname) == 0 :
            saveFile = False
            mods = {}
            for m in chimera.openModels.list() :
                mods[m.name] = m
            for i in range ( 10000 ) :
                nname = os.path.splitext(self.cur_dmap.name)[0] + "_Z%.0f_%d" % (rad,i+1) + ".mrc"
                if not nname in mods :
                    break

        from _multiscale import get_atom_coordinates
        points = get_atom_coordinates ( atoms, transformed = True )
        cmap = PtsToMap ( points, self.cur_dmap, rad, nname, clr=(.7,.7,.7,.2) )

        umsg ( "Made zone map: " + nname )
        self.cur_dmap.display = False

        chimera.runCommand ( "vol #%d style surface region all step 1" % cmap.id )

        if saveFile > 0 :
            mdir, mfile = os.path.split(self.cur_dmap.data.path)
            dpath = os.path.join ( mdir, nname )
            print " -> ", dpath
            cmap.write_file ( dpath, "mrc" )




    def CaBlam ( self ) :

        print ""

        if self.cur_mol == None :
            umsg ( "Select an open model" )
            return


        args = [phPath+'phenix.cablam', "pdb_infile=%s" % self.cur_mol.openedAs[0] ]
        print " - running:"
        for arg in args :
            print arg,
        print ""

        mpath = None
        try :
            mpath, mname = os.path.split(self.cur_mol.openedAs[0])
        except :
            print "-"


        rmap = {}
        for r in self.cur_mol.residues :
            rmap["%s%d"%(r.id.chainId, r.id.position)] = r
            if hasattr ( r, 'ribbonColor' ) :
                r.ribbonColor = chimera.MaterialColor ( .1, .1, .1, 1.0 )

        import subprocess
        #fout = open ( dpath + '/' + '_0_adp.log', "w" )
        p = subprocess.Popen(args, stdout=subprocess.PIPE, cwd=mpath)
        #p.wait()
        #fout.close()
        #out = p.stdout.read()
        i = 0
        for l in p.stdout.readlines() :
            #print i, l,

            if "SUMMARY: " in l :
                ll = l.replace ( "SUMMARY: ", "" )
                print ll,
                continue

            ts = l.split()
            if len(ts) > 2 :
                cid = ts[0]

                try :
                    ri = int ( ts[1] )
                except :
                    continue

                if len(cid) > 3 :
                    continue

                ris = "%s%d"%(cid, ri)
                if not ris in rmap :
                    print "res %d.%s - not found" % (ri, cid)
                else :
                    r = rmap[ris]
                    if "Disfavored" in l :
                        r.ribbonColor = chimera.MaterialColor ( .9, .9, .2, 1.0 )
                    elif "Outlier" in l :
                        r.ribbonColor = chimera.MaterialColor ( .9, .2, .2, 1.0 )
                    elif 0 and "alpha" in l :
                        r.ribbonColor = chimera.MaterialColor ( .9, .2, .9, 1.0 )
                    elif 0 and "beta" in l :
                        r.ribbonColor = chimera.MaterialColor ( .2, .2, .9, 1.0 )
                    else :
                        r.ribbonColor = chimera.MaterialColor ( .6, .6, .6, 1.0 )


            i += 1



    def load ( self, okay, dialog ):
        if okay:
            paths = dialog.getPaths ( )
            if paths:
                path = paths[0]
                umsg ( "Load: " + path )

                if os.path.splitext ( path )[1] == ".cif" :
                    mmcif.LoadMol ( path, log=True )




    def OpenModel ( self ) :

        init = None
        mol = None
        for m in chimera.openModels.list() :
            if type(m) == chimera.Molecule and m.display == True and hasattr ( m, 'openedAs' ) :
                init = os.path.split ( m.openedAs[0] ) [0]
                break
            if type(m) == VolumeViewer.volume.Volume :
                init = os.path.split( m.data.path ) [0]

        if init == None :
            init = "/Users/greg/Box Sync/_data"

        print "init: %s" % init

        if 1 :

            from OpenSave import OpenModeless
            OpenModeless ( title = 'Open Model',
                           #filters = [('TXT', '*.txt', '.txt')],
                           filters = [],
                           initialfile = init, command = self.load )

        else :

            fpath = "/Users/greg/Box Sync/_data/problems/emd_30342/7cec.cif"




    def save ( self, okay, dialog ):
        if okay:
            paths = dialog.getPaths ( )
            if paths:
                path = paths[0]
                umsg ( "Save: " + path )


    def SaveModel ( self ) :
        print "save"

        mol = None
        for m in chimera.openModels.list() :
            if type(m) == chimera.Molecule and m.display == True :
                mol = m

        if 1 :
            init = ""
            if mol :
                init = mol.openedAs[0]
            else :
                init = "/Users/greg/Box Sync/_data"

            from OpenSave import SaveModeless
            SaveModeless ( title = 'Save Model',
                           #filters = [('TXT', '*.txt', '.txt')],
                           filters = [],
                           initialfile = init, command = self.save )

        else :
            mol = None
            for m in chimera.openModels.list() :
                if type(m) == chimera.Molecule :
                    mol = m

            fpath = "/Users/greg/Box Sync/_data/problems/emd_30342/7cec_Q.cif"
            print "Writing %s -> %s" % (mol.name, fpath)

            mmcif.WriteMol ( mol, fpath )





    def RefBack ( self ) :

        molref.RefBack ()


    def RefStep1 ( self ) :

        ress = chimera.selection.currentResidues()
        if len(ress) == 0 :
            umsg ( "Select some residues..." )
            return

        molref.RefStart ( ress, self.cur_dmap )

        startt = time.time()
        molref.RefStep (self.cur_dmap, self.mapF)
        dur = time.time() - startt

        molref.RefPut()

        s = molref.RefE (self.cur_dmap)

        status ( "%.3fs / " % dur + s )



    def StopRef ( self ) :

        #print ""
        self.doRef = False
        self.queueTo.put ("stop")


    def StartRef ( self ) :

        ress = chimera.selection.currentResidues()
        if len(ress) == 0 :
            umsg ( "Select some residues..." )
            return

        molref.RefStart ( ress, self.cur_dmap )

        self.doRef = True
        #self.parent.after(30, self.RefStep)

        import Queue
        self.queue = Queue.Queue()
        self.queueTo = Queue.Queue()
        self.queueTo.put ( "go" )

        molref.RefThread(self.queue, self.queueTo, self.cur_dmap, float(self.mapF.get()) ).start()
        self.parent.after(10, self.process_queue)


    def process_queue(self):
        import Queue

        msg, gotMsg = "", False
        try:
            msg = self.queue.get(0)
            gotMsg = True

        except Queue.Empty:
            pass

        if gotMsg :
            status ( msg )
            molref.RefPut()
            self.queueTo.put ("go")

        if self.doRef :
            #molref.RefThread(self.queue, self.cur_dmap, float(self.mapF.get()) ).start()
            self.parent.after(10, self.process_queue)
        else :
            print "Stop"




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



def GetMod ( name ) :
    for m in chimera.openModels.list() :
        if m.name == name :
            return m
    return None


def GetSegMod () :

    segMod = GetMod ( "SegMod" )

    if segMod != None :
        return segMod

    import _surface
    segMod = _surface.SurfaceModel()
    chimera.openModels.add ( [segMod], sameAs = None )
    segMod.name = "SegMod"
    segMod.mods = []

    print "Created SegMod"

    return segMod




def dialog ( ) :
	from chimera import dialogs
	return dialogs.find ( SegMod_Dialog.name, create=False )



def show_dialog ( closeOld = True ):

	from chimera import dialogs

	d = dialogs.find ( SegMod_Dialog.name, create=False )
	if d :
		if closeOld :
			d.toplevel_widget.update_idletasks ()
			d.Close()
			d.toplevel_widget.update_idletasks ()
		else :
			# is there a way to bring it to front?
			return d

	dialogs.register (SegMod_Dialog.name, SegMod_Dialog, replace = True)

	d = dialogs.find ( SegMod_Dialog.name, create=True )
	d.toplevel_widget.update_idletasks ()
	d.enter()

	return d
