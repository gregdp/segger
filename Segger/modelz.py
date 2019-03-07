
# Copyright (c) 2018 Greg Pintilie - gregp@slac.stanford.edu

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
import tkFont
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
import FitMap
from sys import stderr
from time import clock
import _contour
import chimera.match

from axes import prAxes
import _multiscale
from CGLutil.AdaptiveTree import AdaptiveTree
import random
from VolumePath import Marker_Set, Marker, Link
from _contour import affine_transform_vertices as transform_vertices
from Matrix import xform_matrix, multiply_matrices, chimera_xform, identity_matrix, invert_matrix, shift_and_angle
import struct


from Rotamers import getRotamers
from chimera.resCode import protein1to3


OML = chimera.openModels.list

devMenu = True


atomColors = {'C' : chimera.MaterialColor (0.565,0.565,0.565),
            'Cbb' : chimera.MaterialColor (0.2,0.6,0.2),
            'S' : chimera.MaterialColor (1.000,1.000,0.188),
            'O' : chimera.MaterialColor (1.000,0.051,0.051),
            'N' : chimera.MaterialColor (0.188,0.314,0.973),
            'P' : chimera.MaterialColor (1.0, 0.502, 0.0),
            'H' : chimera.MaterialColor (0.9,.9,.9) 
            }




def umsg ( txt ) :
    print txt
    status ( txt )

def status ( txt ) :
    txt = txt.rstrip('\n')
    msg.configure(text = txt)
    msg.update_idletasks()


class ModelZ_Dialog ( chimera.baseDialog.ModelessDialog ):

    title = "ModelZ v1.1; press 'Help' at the bottom right corner for more details and to cite this tool."
    name = "modelz"
    buttons = ( "Close" )
    help = 'https://cryoem.slac.stanford.edu/ncmi/resources/software/modelz'


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

        self.InitVars()


        if 1 :
            #row += 1
            ff = Tkinter.Frame(f)
            ff.grid(column=0, row=row, sticky='nsew', pady=0, padx=0)

            Tkinter.Grid.columnconfigure(f, 0, weight=1)
            Tkinter.Grid.columnconfigure(ff, 0, weight=1)

            Tkinter.Grid.rowconfigure(f, row, weight=1)
            Tkinter.Grid.rowconfigure(ff, 0, weight=1)


            self.Canvas = Tkinter.Canvas(ff, height=80)
            self.Canvas.grid(column=0, row=0, sticky='nsew')

            self.modX = 10; self.modY = 10; self.modH = 30
            self.seqX = 10; self.seqY = 45; self.seqH = 30

            self.Canvas.bind("<ButtonPress-1>", lambda event : self.B1_Down ( event ) )
            self.Canvas.bind("<Control-ButtonPress-1>", lambda event : self.B1_Down_Ctrl ( event ) )
            self.Canvas.bind("<Shift-ButtonPress-1>", lambda event : self.B1_Down_Shift ( event ) )
            self.Canvas.bind("<Option-ButtonPress-1>", lambda event : self.B1_Down_Alt ( event ) )
            self.Canvas.bind("<Alt-ButtonPress-1>", lambda event : self.B1_Down_Alt ( event ) )
            self.Canvas.bind("<ButtonPress-2>", lambda event : self.B2_Down (event) )
            self.Canvas.bind("<ButtonPress-3>", lambda event : self.B3_Down (event) )
            self.Canvas.bind("<ButtonRelease-1>", lambda event : self.B1_Up ( event ) )
            self.Canvas.bind("<Control-ButtonRelease-1>", lambda event : self.B1_Up_Ctrl ( event ) )
            self.Canvas.bind("<Shift-ButtonRelease-1>", lambda event : self.B1_Up_Shift ( event ) )
            self.Canvas.bind("<Alt-ButtonRelease-1>", lambda event : self.B1_Up_Alt ( event ) )
            self.Canvas.bind("<Option-ButtonRelease-1>", lambda event : self.B1_Up_Alt ( event ) )

            self.Canvas.bind("<ButtonRelease-2>", lambda event : self.B2_Up (event) )
            self.Canvas.bind("<Option-ButtonRelease-2>", lambda event : self.B2_Up_Alt (event) )
            self.Canvas.bind("<Alt-ButtonRelease-2>", lambda event : self.B2_Up_Alt (event) )
            self.Canvas.bind("<Control-ButtonRelease-2>", lambda event : self.B2_Up_Ctrl (event) )
            self.Canvas.bind("<Command-ButtonRelease-2>", lambda event : self.B2_Up_Comm (event) )
            self.Canvas.bind("<Shift-ButtonRelease-2>", lambda event : self.B2_Up_Shift (event) )

            self.Canvas.bind("<ButtonRelease-3>", lambda event : self.B2_Up (event) )
            self.Canvas.bind("<Option-ButtonRelease-3>", lambda event : self.B2_Up_Alt (event) )
            self.Canvas.bind("<Alt-ButtonRelease-3>", lambda event : self.B2_Up_Alt (event) )
            self.Canvas.bind("<Control-ButtonRelease-3>", lambda event : self.B2_Up_Ctrl (event) )
            self.Canvas.bind("<Command-ButtonRelease-3>", lambda event : self.B2_Up_Comm (event) )
            self.Canvas.bind("<Shift-ButtonRelease-3>", lambda event : self.B2_Up_Shift (event) )

            self.Canvas.bind("<B1-Motion>", lambda event : self.B1_Drag ( event ) )
            self.Canvas.bind("<B2-Motion>", lambda event : self.B2_Drag ( event ) )
            self.Canvas.bind("<B3-Motion>", lambda event : self.B3_Drag ( event ) )
            self.Canvas.bind("<Motion>", lambda event : self.Mouse_Move ( event ) )
            self.Canvas.bind("<Configure>", lambda event : self.Canvas_Config (event) )
            self.Canvas.bind("<Leave>", lambda event : self.Canvas_Leave (event) )
            self.Canvas.bind("<MouseWheel>", lambda event : self.Canvas_Wheel (event) )


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w', pady=0, padx=5)

        if 1 :
            ff = Tkinter.Frame(f)
            ff.grid(column=0, row=row, sticky='w', pady=5, padx=10)

            l = Tkinter.Label(ff, text='Map:', anchor=Tkinter.W)
            l.grid(column=0, row=0, sticky='w')

            self.dmap = Tkinter.StringVar(parent)
            self.dmapMB  = Tkinter.Menubutton ( ff, textvariable=self.dmap, relief=Tkinter.RAISED, width=20 )
            self.dmapMB.grid (column=1, row=0, sticky='we', padx=5)
            self.dmapMB.menu  =  Tkinter.Menu ( self.dmapMB, tearoff=0, postcommand=self.MapMenu )
            self.dmapMB["menu"]  =  self.dmapMB.menu

            self.cur_dmap = None
            self.SetVisMap ()


            l = Tkinter.Label(ff, text='Model:', anchor=Tkinter.W)
            l.grid(column=2, row=0, sticky='w')

            self.struc = Tkinter.StringVar(parent)
            self.strucMB  = Tkinter.Menubutton ( ff, textvariable=self.struc, relief=Tkinter.RAISED, width=20 )
            self.strucMB.grid (column=3, row=0, sticky='we', padx=5)
            self.strucMB.menu  =  Tkinter.Menu ( self.strucMB, tearoff=0, postcommand=self.StrucMenu )
            self.strucMB["menu"]  =  self.strucMB.menu

            self.cur_mol = None
            self.cur_chains = []
            self.SetVisMol ()


            l = Tkinter.Label(ff, text=" Chain:" )
            l.grid(column=4, row=0, sticky='w')

            self.chain = Tkinter.StringVar(parent)
            self.chainMB  = Tkinter.Menubutton ( ff, textvariable=self.chain, relief=Tkinter.RAISED, width=4 )
            self.chainMB.grid (column=5, row=0, sticky='we', padx=5)
            self.chainMB.menu  =  Tkinter.Menu ( self.chainMB, tearoff=0, postcommand=self.ChainMenu )
            self.chainMB["menu"]  =  self.chainMB.menu

            if len ( self.cur_chains ) > 0 :
                self.chain.set ( self.cur_chains[0] )
                #self.ShowCh ( self.cur_chains[0] )
                self.GetSeq ()
                

            b = Tkinter.Button(ff, text="Show Chain", command=self.AllChain)
            b.grid (column=6, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Show All", command=self.AllChains)
            b.grid (column=7, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(ff, text="RandColor", command=self.RandColorChains )
            #b.grid (column=7, row=0, sticky='w', padx=5)

            #l = Tkinter.Label(ff, text=' Z-Scores:', fg="#777")
            #l.grid(column=8, row=0, sticky='e')

            #b = Tkinter.Button(ff, text="SSE", command=self.SSE)
            #b.grid (column=9, row=0, sticky='w', padx=5)

            if 1 :
                oft = Hybrid.Checkbutton(ff, 'Ribbon', False)
                #oft.button.grid(column = 12, row = 0, sticky = 'w')
                self.showRibbon = oft.variable
                #self.showRibbon.set ( 1 )


            #oft = Hybrid.Checkbutton(ff_color, 'SSE', False)
            #oft.button.grid(column = 21, row = 0, sticky = 'w')
            #self.colorSSE = oft.variable
            #self.colorSSE.set ( 0 )

            #oft = Hybrid.Checkbutton(ff_color, 'SC', False)
            #oft.button.grid(column = 22, row = 0, sticky = 'w')
            #self.colorSC = oft.variable
            #self.colorSC.set ( 0 )

            #oft = Hybrid.Checkbutton(ff_color, 'Rand', False)
            #oft.button.grid(column = 23, row = 0, sticky = 'w')
            #self.colorRand = oft.variable
            #self.colorRand.set ( 0 )

            #oft = Hybrid.Checkbutton(ff_color, 'Map', False)
            #oft.button.grid(column = 24, row = 0, sticky = 'w')
            #self.colorMap = oft.variable
            #self.colorMap.set ( 1 )

            #b = Tkinter.Button(ff_color, text="Update", command=self.DoColor)
            #b.grid (column=25, row=0, sticky='w', padx=5)



            #b = Tkinter.Button(ff, text="Res-B", command=self.ResB)
            #b.grid (column=11, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(ff, text="Mask", command=self.Mask)
            #b.grid (column=12, row=0, sticky='w', padx=5)


        if 1 :

            l = Tkinter.Label(ff, text='        Zoom:', fg="#777")
            l.grid(column=35, row=0, sticky='e')

            b = Tkinter.Button(ff, text="-", command=self.ZoomMinus)
            b.grid (column=36, row=0, sticky='w', padx=0)

            b = Tkinter.Button(ff, text="+", command=self.ZoomPlus)
            b.grid (column=37, row=0, sticky='w', padx=0)

            b = Tkinter.Button(ff, text="<", command=self.ZoomBegin)
            b.grid (column=38, row=0, sticky='w', padx=0)

            b = Tkinter.Button(ff, text=">", command=self.ZoomEnd)
            b.grid (column=39, row=0, sticky='w', padx=0)


        self.whichScore = Tkinter.StringVar()
        self.whichScore.set ( 'Z' )

        if 0 :

            row += 1
            ff = Tkinter.Frame(f)
            ff.grid(column=0, row=row, sticky='w', pady=0, padx=5)

            b = Tkinter.Button(ff, text="Score:", command=self.CalcZScores )
            b.grid (column=10, row=0, sticky='w', padx=5)


            c = Tkinter.Radiobutton(ff, text="Z", variable=self.whichScore, value = 'Z')
            c.grid (column=11, row=0, sticky='w')
            
            c = Tkinter.Radiobutton(ff, text="CC", variable=self.whichScore, value = 'CC')
            c.grid (column=12, row=0, sticky='w')

            c = Tkinter.Radiobutton(ff, text="CC (Mean)", variable=self.whichScore, value = 'CCm')
            c.grid (column=13, row=0, sticky='w')

            #c = Tkinter.Radiobutton(ff, text="Q", variable=self.whichScore, value = 'Q')
            #c.grid (column=14, row=0, sticky='w')


        if 1 :

            row += 1

            ff = Tkinter.Frame(f)
            ff.grid(column=0, row=row, sticky='w', pady=0, padx=5)


            fff = Tkinter.Frame(ff, borderwidth=1, padx=2, pady=2, relief=Tkinter.GROOVE)
            fff.grid(column=10, row=0, sticky='e', pady=0, padx=5)

            l = Tkinter.Label(fff, text='Calculate:', fg="#777")
            l.grid(column=10, row=0, sticky='e')

            b = Tkinter.Button(fff, text="Z-scores", command=self.CalcZScores )
            b.grid (column=11, row=0, sticky='w', padx=5)

            b = Tkinter.Button(fff, text="CC", command=self.CalcCCScores )
            b.grid (column=12, row=0, sticky='w', padx=5)

            b = Tkinter.Button(fff, text="CC(mean)", command=self.CalcCCmScores )
            b.grid (column=13, row=0, sticky='w', padx=5)


        if 1 :

            self.colorMod = Tkinter.StringVar()
            self.colorMod.set ( 'sc' )
            
            fff = Tkinter.Frame(ff, borderwidth=1, padx=2, pady=2, relief=Tkinter.GROOVE)
            fff.grid(column=20, row=0, sticky='e', pady=0, padx=5)

            b = Tkinter.Button(fff, text="Color:", command=self.DoColor)
            b.grid (column=20, row=0, sticky='w', padx=5)
            
            c = Tkinter.Radiobutton(fff, text="Backbone", variable=self.colorMod, value = 'bb')
            c.grid (column=21, row=0, sticky='w')
            
            c = Tkinter.Radiobutton(fff, text="Side Chains", variable=self.colorMod, value = 'sc')
            c.grid (column=22, row=0, sticky='w')
            
            c = Tkinter.Radiobutton(fff, text="Random", variable=self.colorMod, value = 'rand')
            c.grid (column=23, row=0, sticky='w')


            l = Tkinter.Label(fff, text='', fg="#000")
            l.grid(column=25, row=0, sticky='ens')

        if 1 :
            #row += 1
            #ff = Tkinter.Frame(f)
            #ff.grid(column=0, row=row, sticky='w', pady=0, padx=5)

            ff = Tkinter.Frame(ff, borderwidth=1, padx=2, pady=2, relief=Tkinter.GROOVE)
            ff.grid(column=30, row=0, sticky='e', pady=0, padx=5)

            #l = Tkinter.Label(ff, text='Select (Ctrl+Click+Drag On Sequence) Show:', fg="#000")
            l = Tkinter.Label(ff, text='Show:', fg="#000")
            l.grid(column=35, row=0, sticky='ens')

            oft = Hybrid.Checkbutton(ff, 'Ribbon', True)
            oft.button.grid(column = 36, row = 0, sticky = 'w')
            self.showRibbon = oft.variable
            #self.showRibbon.set ( 1 )

            oft = Hybrid.Checkbutton(ff, 'Side Chains', True)
            oft.button.grid(column = 37, row = 0, sticky = 'w')
            self.showAtoms = oft.variable
            #self.showRibbon.set ( 1 )

            oft = Hybrid.Checkbutton(ff, 'Mesh', False)
            oft.button.grid(column = 38, row = 0, sticky = 'w')
            self.showMesh = oft.variable
            #self.showRibbon.set ( 1 )

            #oft = Hybrid.Checkbutton(ff, 'Preserve', False, command=self.cb)
            #oft.button.grid(column = 39, row = 0, sticky = 'w')
            #self.preserveSel = oft.variable
            self.preserveSel = Tkinter.IntVar()
            oft = Tkinter.Checkbutton( ff, text="Preserve", variable=self.preserveSel, command=self.preserveSelCb)
            oft.grid(column = 39, row = 0, sticky = 'w')
            #self.showRibbon.set ( 1 )

            #b = Tkinter.Button(ff, text="Clear", command=self.ClearSel)
            #b.grid (column=40, row=0, sticky='w', padx=5)

            #self.keepExMap = Tkinter.IntVar()
            #self.keepExMap.set(0)
            #oft = Tkinter.Checkbutton( ff, text="Keep Extracted Maps", variable=self.keepExMap, command=self.keepExMapCb)
            #oft.grid(column = 40, row = 0, sticky = 'w')

        if 0 and devMenu :

            row += 1
            ff = Tkinter.Frame(f)
            ff.grid(column=0, row=row, sticky='w', pady=0, padx=5)
            
            
            b = Tkinter.Button(ff, text="Asp", command=self.asp )
            b.grid (column=1, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Extract Res", command=self.Extract )
            b.grid (column=2, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Align 1", command=self.AlignRes1 )
            b.grid (column=3, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Align 2", command=self.AlignRes2 )
            b.grid (column=4, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Avg", command=self.Avg )
            b.grid (column=5, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Close", command=self.CloseExtracted )
            b.grid (column=6, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Rad1", command=self.RadScore )
            b.grid (column=7, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="RadBB", command=self.RadScoreBB )
            b.grid (column=8, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Rads", command=self.RadScores )
            b.grid (column=9, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="ExA", command=self.ExCustA )
            b.grid (column=10, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="ExB", command=self.ExCustB )
            b.grid (column=11, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="ExC", command=self.ExCustC )
            b.grid (column=12, row=0, sticky='w', padx=5)


        dummyFrame = Tkinter.Frame(parent, relief='groove', borderwidth=1)
        Tkinter.Frame(dummyFrame).pack()
        dummyFrame.grid(row=row,column=0,columnspan=7, pady=3, sticky='we')
        row += 1


        global msg
        msg = Tkinter.Label(parent, width = 60, anchor = 'w', justify = 'left', fg="red", pady=5, padx=10)
        msg.grid(column=0, row=row, sticky='ew')
        self.msg = msg

        #umsg ( 'Select one or more segmented regions then press "Place Points" to start' )

        
    def InitVars ( self ) :

        self.mag = 13
        self.seqt = []
        self.boldSeqT = None
        self.drag = ''

        #self.sheetBaseClr = numpy.array ( [50.0,205.0,50.0] )
        #self.sheetClr = numpy.array ( [204.0,255.0,204.0] )
        self.sheetBaseClr = numpy.array ( [55.0,55.0,150.0] )
        self.sheetClr = numpy.array ( [150.0,150.0,250.0] )
        self.sheetClrD = self.sheetClr - self.sheetBaseClr

        self.helixBaseClr = numpy.array ( [150.0,50.0,50.0] )
        self.helixClr = numpy.array ( [255.0,150.0,150.0] )
        self.helixClrD = self.helixClr - self.helixBaseClr
        
        c = self.helixBaseClr; self.helix1 = "#" + struct.pack('BBB',c[0],c[1],c[2]).encode('hex')
        c = self.helixClr;     self.helix2 = "#" + struct.pack('BBB',c[0],c[1],c[2]).encode('hex')
        
        self.switch = "#522"
        
        c = self.sheetBaseClr; self.strand1 = "#77F"
        c = self.sheetClr;     self.strand2 = "#77F"

        c = self.sheetBaseClr; self.sheet1 = "#" + struct.pack('BBB',c[0],c[1],c[2]).encode('hex')
        c = self.sheetClr;     self.sheet2 = "#" + struct.pack('BBB',c[0],c[1],c[2]).encode('hex')

        self.loop1 = "#999"
        
        self.selColor = "#7e7"


        self.font = tkFont.Font(family='Courier', size=(self.mag), weight='normal')
        #self.boldFont = tkFont.Font(family='Courier', size=(self.mag+4), weight='bold')
        self.tw = self.font.measure ( "a" )
        
        self.seq = ""

        #self.OrderMods ()


    def SetVisMap ( self ) :
        dmap = None
        mlist = OML(modelTypes = [VolumeViewer.volume.Volume])
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
        self.dmapMB.menu.delete ( 0, 'end' )   # Clear menu
        mlist = OML(modelTypes = [VolumeViewer.volume.Volume])
        for m in mlist :
            self.dmapMB.menu.add_radiobutton ( label=m.name+" (%d)"%m.id, variable=self.dmap,
                                command=lambda m=m: self.MapSelected(m) )


    def MapSelected ( self, dmap ) :

        self.cur_dmap = dmap    
        print "Selected " + dmap.name

        self.GetSeq ()
        self.ZoomBegin ()


    def GetChains ( self, mol ) :
        ct = {}
        for r in mol.residues: 
            ct[r.id.chainId] = 1
        clist = ct.keys()
        clist.sort()
        return clist
        

    def SetVisMol ( self ) :
        mol = None
        mlist = OML(modelTypes = [chimera.Molecule])
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
            self.cur_chains = self.GetChains ( mol )
            SetBBAts ( mol )


    def StrucSelected ( self, mol ) :

        self.cur_mol = mol
        print "Selected ", mol.name, " - ", mol.id
        if mol :

            mlist = OML(modelTypes = [chimera.Molecule])
            for m in mlist :
                m.display

            mol.display = True

            self.cur_chains = self.GetChains ( mol )

            if len(self.cur_chains) == 0 :
                self.chain.set ( "" )
            elif self.chain.get() in self.cur_chains :
                print " - ch " + self.chain.get() + " already sel"
                self.ShowCh ( self.chain.get() )
            else :
                self.chain.set ( self.cur_chains[0] )
                self.ShowCh ( self.chain.get() )
            
            self.GetSeq ()
            self.ZoomBegin ()
            SetBBAts ( mol )


    
    def ChainSelected ( self, ch ) :
        print " - sel chain: ", ch, self.chain.get()
        self.ShowCh ( ch )
        self.GetSeq ()
        self.ZoomBegin ()


    def StrucMenu ( self ) :
        self.strucMB.menu.delete ( 0, 'end' )   # Clear menu
        mlist = OML(modelTypes = [chimera.Molecule])
        for m in mlist :
            self.strucMB.menu.add_radiobutton ( label=m.name+" (%d)"%m.id, variable=self.struc,
                                           command=lambda m=m: self.StrucSelected(m) )

    def ChainMenu ( self ) :
        self.chainMB.menu.delete ( 0, 'end' )   # Clear menu
        print " - chain menu"
        print self.cur_chains
        for ch in self.cur_chains :
            self.chainMB.menu.add_radiobutton ( label=ch, variable=self.chain, 
                                            command=lambda ch=ch: self.ChainSelected(ch) )

        
    
    def DoColor ( self ) :
        
        print "color...", self.colorMod.get()
        
        #colSC = self.colorSC.get()
        #colRand = self.colorRand.get()
        
        if self.colorMod.get() == "rand" :
            self.RandColorChains()
        else :
            self.UpdateModColor ()

        #if self.colorMap.get() :
        #    self.UpdateSurfColor ()



    def UpdateSurfColor ( self ) :

        print " - surf of %s, by %s" % ( self.cur_dmap.name, self.cur_mol.name )

        numAt = 0
        for r in self.cur_mol.residues :
            for at in r.atoms :
                if "H" in at.name :
                    pass
                else :
                    numAt += 1

        allAtPos = numpy.zeros ( (numAt, 3) )
        allAts = [None] * numAt

        numAt = 0
        for r in self.cur_mol.residues :
            for at in r.atoms :
                if "H" in at.name :
                    pass
                else :
                    allAtPos[numAt] = at.coord().data()
                    allAts[numAt] = at
                    at.allPtI = numAt
                    numAt += 1


        print " - tree with %d ats" % numAt
        allAtTree = AdaptiveTree ( allAtPos.tolist(), allAts, 4.0)
        print " - done"
        
        
        


    def UpdateModColor ( self ) :

        ress = []
        try :
            ress = self.seqRes
        except :
            pass

        if len ( ress ) == 0 :
            umsg ( "No molecule/chain selected?" )
            return

        if not hasattr (self, 'scores') :
            umsg ( "No scores - press 'Z-Scores' button first" )
            return

        foundScore = False            
        for sc in self.scores :
            if sc != None :
                foundScore = True
        
        if not foundScore :
            umsg ( "No scores - press 'Z-Scores' button first" )
            return


        ac = { 'O' : chimera.MaterialColor( .9, .2, .2, 1.0 ),
                'C' : chimera.MaterialColor( .7, .7, .7, 1.0 ),
                'N' : chimera.MaterialColor( .2, .2, .9, 1.0 ),
                'H' : chimera.MaterialColor( 1, 1, 1, 1.0 ),
                'S' : chimera.MaterialColor( .9, .9, 0, 1.0 ),
                ' ' : chimera.MaterialColor( .2, .2, .2, 1.0 ),
                 }

        minScore, maxScore = 0,0
        colorSC = self.colorMod.get() == "sc"
        if colorSC : 
            minScore, maxScore = self.minSCscore, self.maxSCscore
        else :
            minScore, maxScore = self.minBBscore, self.maxBBscore

        cH = numpy.array( [0.0,1.0,0.0] )
        cL = numpy.array( [1.0,0.0,0.0] )

        for ri, r in enumerate ( self.seqRes ) :
            sc = None
            #sc = self.scores[ri] if colorSC else self.scores2[ri]
            sc = r.scZ if colorSC else r.bbZ

            if sc == None  :
                r.ribbonColor = chimera.MaterialColor ( .7, .7, .7, 1.0 )
                for at in r.atoms :
                    #at.color = r.ribbonColor
                    try :
                        at.color = ac[at.name[0]]
                    except :
                        at.color = ac[' ']

            else :
                h = (sc - minScore) / (maxScore - minScore)
                if h > 1 : h = 1
                if h < 0 : h = 0
                c = h * cH + (1-h) * cL
                r.ribbonColor = chimera.MaterialColor ( c[0], c[1], c[2], 1.0 )
                for at in r.atoms :
                    #at.color = r.ribbonColor
                    try :
                        at.color = ac[at.name[0]]
                    except :
                        at.color = ac[' ']
            
                #ra = r.scZ, r.bbZ

        
        
    
    def RandColorChains ( self ) :
    
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        m = self.cur_mol

        from random import random as rand
    
        ct = {}
        for r in m.residues: ct[r.id.chainId] = 1
        clist = ct.keys()
        clist.sort()
        chains_clrs = {}
        cnames = ""
    
        for ci, cid in enumerate ( clist ) :
            clr = ( rand()*.8+.1, rand()*.8+.1, rand()*.8+.1 )
            chains_clrs[cid] = chimera.MaterialColor ( clr[0], clr[1], clr[2], 1.0 )
            cnames = cnames + cid
    
        print "%s - color ribbon for %d chains -" % ( m.name, len(cnames) ), cnames
    
        # color atoms
        for r in m.residues :
            clr = chains_clrs[r.id.chainId]
            r.ribbonColor = clr
            for at in r.atoms :
                at.color = clr
    
    
    def AllChain ( self ) :

        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        chainId = self.chain.get()
        if len(chainId) == 0 :
            umsg ("Select a chain first")
            return
            
        umsg ( "Showing mol %s chain %s" % (self.cur_mol.name, chainId) )

        #ct = {}
        #for r in self.cur_mol.residues: ct[r.id.chainId] = 1
        #clist = ct.keys()
        #clist.sort()

        for r in self.cur_mol.residues :
            if r.id.chainId == chainId :
                if ("CA" in r.atomsMap and "N" in r.atomsMap and "C" in r.atomsMap) or ("O3'" in r.atomsMap and "O5'" in r.atomsMap)  :
                    r.ribbonDisplay = True
                    r.ribbonDrawMode = 2
                else :
                    r.ribbonDisplay = False
                    for at in r.atoms :
                        at.drawMode = at.Ball
                        at.display = True
            else :
                if ("CA" in r.atomsMap and "N" in r.atomsMap and "C" in r.atomsMap) or ("O3'" in r.atomsMap and "O5'" in r.atomsMap)  :
                    r.ribbonDisplay = False
                    r.ribbonDrawMode = 2
                else :
                    r.ribbonDisplay = False
                    for at in r.atoms :
                        at.drawMode = at.Ball
                        at.display = False


    def AllChains ( self ) :
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        m = self.cur_mol
        
        #ct = {}
        #for r in m.residues: ct[r.id.chainId] = 1
        #clist = ct.keys()
        #clist.sort()

        for r in m.residues :
            if ("CA" in r.atomsMap and "N" in r.atomsMap and "C" in r.atomsMap) or ("O3'" in r.atomsMap and "O5'" in r.atomsMap)  :
                r.ribbonDisplay = True
                r.ribbonDrawMode = 2
            else :
                r.ribbonDisplay = False
                for at in r.atoms :
                    #at.drawMode = at.Ball
                    at.display = True
        
    
    def ShowCh ( self, ch ) :

        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        print " - showing chain:", ch

        m = self.cur_mol
        print " - cur mol:", m.name

        ct = {}
        for r in m.residues: ct[r.id.chainId] = 1
        clist = ct.keys()
        clist.sort()

        for r in m.residues :
            show = True if r.id.chainId == ch else False
            if ("CA" in r.atomsMap and "N" in r.atomsMap and "C" in r.atomsMap) or ("O3'" in r.atomsMap and "O5'" in r.atomsMap) :
                r.ribbonDisplay = show
                #r.ribbonDrawMode = 2
                for at in r.atoms :
                    at.display = False
            else :
                r.ribbonDisplay = False
                for at in r.atoms :
                    #at.drawMode = at.Ball
                    at.display = show


    def GetMod ( self, name ) :
        for m in chimera.openModels.list() :
            if name != None and len(name) > 0 :
                if m.name == name :
                    return m
            else :
                if m.display == True :
                    return m
        return None



    def GetSeq ( self ) :
        
        if self.cur_mol == None :
            umsg ( "No selected molecule" )
            return
        
        if len ( self.chain.get() ) == 0 :
            umsg ( "No selected chain" )
            return
            
        self.RemoveSeq ()

        try :
            print self.cur_mol.name
        except :
            print " - mol may have been closed"
            return

        self.GetSeqFromStruc ( self.cur_mol, self.chain.get() )
        
        if len(self.seq) > 0 :

            print "-- seq from open mol -- %d res" % len(self.seq)
            print self.seq

            self.seqt = []
            self.seqSheetR = [None] * len(self.seq)
            self.seqHelixR = [None] * len(self.seq)
            self.seqScoreR = [None] * len(self.seq)
            self.seqScoreR2 = [None] * len(self.seq)
            self.scores2 = [None] * len(self.seq)
            self.scores = [None] * len(self.seq)
            
            self.UpdateSeqFont ()
            
            return True
        
        return False



    def RemoveSeq  (self) :
        
        if self.seq == "" :
            return

        for si in range ( len(self.seq) ) :
            res = self.seq[si]
            pred = self.pred[si]
            conf = float ( self.conf[si] ) / 10.0

            if pred == 'E' :
                if self.seqSheetR[si] != None :
                    self.Canvas.delete ( self.seqSheetR[si] )

            elif pred == 'H' :
                if self.seqHelixR[si] != None :
                    self.Canvas.delete ( self.seqHelixR[si] )

            if self.seqScoreR[si] != None :
                self.Canvas.delete ( self.seqScoreR[si] )

            if self.seqScoreR2[si] != None :
                self.Canvas.delete ( self.seqScoreR2[si] )


        # box showing selected Residue
        if hasattr ( self, 'seqMouseR' ) :
            self.Canvas.delete ( self.seqMouseR )
            del self.seqMouseR

        if hasattr ( self, 'seqText' ) :
            self.Canvas.delete ( self.seqText )
            del self.seqText
            
        self.seqSel = None
        self.seq = ""
        self.UpdateSeqSel ()



    def GetSeqFromStruc ( self, mol, chainId ) :

        print "Getting seq from %s, %s" % (mol.name, chainId)

        self.conf = ""
        self.pred = ""
        self.seq = ""
        self.seqRes = []

        from chimera.resCode import protein3to1
        from chimera.resCode import nucleic3to1
        protein3to1['HSD'] = protein3to1['HIS']

        rids = {}
        for r in mol.residues :
            if r.id.chainId == chainId :
                if r.type in protein3to1 or r.type in nucleic3to1 :
                    rids[r.id.position] = r


        ris = rids.keys()
        ris.sort()

        for ri in ris :
            r = rids[ri]
            if r.type in protein3to1 :
                self.seq = self.seq + protein3to1[r.type]
                self.conf = self.conf + "9"
                self.predi = "C"
                if r.isSheet : 
                    self.predi = "E"
                if r.isHelix :
                    self.predi = "H"
                self.pred = self.pred + self.predi
                self.seqRes.append ( r )
            elif r.type in nucleic3to1 :
                self.seq = self.seq + nucleic3to1[r.type]
                self.conf = self.conf + "9"
                self.predi = "C"
                self.pred = self.pred + self.predi
                self.seqRes.append ( r )
    
    


    def SSE ( self ) :

        print "sse"
        #self.GetFromMol ( mod, chainId )


    def CurRes ( self ) :

        #self.GetFromMol ( mod, chainId )
        
        if self.cur_mol == None :
            umsg ( "No selected molecule" )
            return []
        
        if self.cur_dmap == None :
            umsg ( "No selected map" )
            return []
            
        if len ( self.chain.get() ) == 0 :
            umsg ( "No selected chain" )
            return []
        
        from chimera.resCode import protein3to1
        protein3to1['HSD'] = protein3to1['HIS']
        
        rids = {}
        for r in self.cur_mol.residues :
            if r.id.chainId == self.chain.get() :
                if r.type in protein3to1 :
                    rids[r.id.position] = r
        
        print " - %d residues" % len(rids.values())
        return [ rids[6] ]
        #return rids.values ()
        
        
        
    def CalcZScores ( self ) :
        self.CalcScores ( 'Z' )
    
    def CalcCCScores ( self ) :
        self.CalcScores ( 'CC' )
    
    def CalcCCmScores ( self ) :
        self.CalcScores ( 'CCm' )


    def CalcScores ( self, t ) :

        self.whichScore.set ( t )

        from chimera import tasks, CancelOperation
        task = tasks.Task('Calculating %s-scores' % t, modal = True)

        try:
            self.CalcZScores_( task )
        except CancelOperation:
            umsg('Cancelled %s-scores' % t)
        finally:
            task.finished()


    def CalcZScores_ ( self, task ) :

        ress = []
        try :
            ress = self.seqRes
        except :
            pass
            
        if len ( ress ) == 0 :
            umsg ( "No molecule/chain selected?" )
            return

        self.scores2 = [None] * len(self.seqRes)
        scoreI = 0

        status ( "Getting secondary structure elements..." )

        ok = True
        try :
            print self.cur_dmap.name
        except :
            status ( "Selected map not found; please choose another map" )
            self.dmap.set ("")
            ok = False
            
        try :
            print self.cur_mol.name
        except :
            status ( "Selected model not found; please choose another model" )
            self.struc.set ("")
            self.chain.set ("")
            self.RemoveSeq ()
            ok = False
        
        if not ok :
            return


        resolution = 3.0 * self.cur_dmap.data.step[0]
        #resolution = 3.0
        umsg ( "Calculating backbone Z-scores..." )


        #if self.whichScore.get() == "Q" :
        self.bbCC, self.bbCCM = ccBB ( self.cur_mol, self.seqRes, resolution, self.cur_dmap )
        print " - all bb cc: %.3f %.3f" % (self.bbCC, self.bbCCM)

        self.bbAvgD = avgdBB ( self.cur_mol, self.seqRes, self.cur_dmap )
        print " - all bb avgd: %.3f" % self.bbAvgD


        zscores2 = []

        if 0 : # old - use SSE segments
            sses = SSEs ( self.seqRes )
            #print " - ",len(sses),"sse for ", len(ress), "res"

            atI = 1

            for el in sses :
                si, ei, ss, elRess = el

                if atI % 10 == 0 :
                    status ( "BB scores: %d/%d" % (atI,len(sses) ) )
                atI += 1

                #if 1 or (startRes < 129 and endRes > 129) :
                startResI, endResI, sseType, ress = el
                #print " : %d-%d, %s, %d res" % (startResI, endResI, sseType, len(ress))

                score = None
                if self.whichScore.get() == "Z" :
                    score, ccs = zBB ( self.cur_mol, ress, resolution, self.cur_dmap )

                elif self.whichScore.get() == "CC" :
                    score, ccs = zBB ( self.cur_mol, ress, resolution, self.cur_dmap )

                elif self.whichScore.get() == "CCm" :
                    cc, score = ccBB ( self.cur_mol, ress, resolution, self.cur_dmap )

                elif self.whichScore.get() == "Q" :
                    cc, ccm = ccBB ( self.cur_mol, ress, resolution, self.cur_dmap )
                    score = cc / self.bbCC

                #print ss, si, "-", ei, zscore

                if score != None :
                    zscores2.append ( score )

                for r in elRess :
                    r.bbZ = score
                    self.scores2[scoreI] = score
                    scoreI += 1

        else :

            bbs = BBsegs ( self.seqRes )
            W = 3
            atRes = 0

            for bb in bbs :
                print "%d res, %d-%d" % (len(bb),bb[0].id.position,bb[-1].id.position)

                for ri, r in enumerate ( bb ) :
                    firstRi = max ( 0, ri-(W-1)/2 )
                    lastRi = min ( len(bb)-1, ri+(W-1)/2 )
                    ress = bb[firstRi:lastRi+1]

                    #zscore, ccs = zBB ( self.cur_mol, ress, resolution, self.cur_dmap )

                    score = None
                    if self.whichScore.get() == "Z" :
                        score, ccs = zBB ( self.cur_mol, ress, resolution, self.cur_dmap )
    
                    elif self.whichScore.get() == "CC" :
                        score, ccs = zBB ( self.cur_mol, ress, resolution, self.cur_dmap )
    
                    elif self.whichScore.get() == "CCm" :
                        cc, score = ccBB ( self.cur_mol, ress, resolution, self.cur_dmap )
    
                    elif self.whichScore.get() == "Q" :
                        cc, ccm = ccBB ( self.cur_mol, ress, resolution, self.cur_dmap )
                        score = cc / self.bbCC

                    #print "  %d : %d - %d, %.3f" % (ri, firstRi, lastRi, zscore)
                    if atRes % 10 == 0 :
                        #status ( "Backbone - residue %d/%d" % (atRes,len(self.seqRes) ) )
                        #print "%d/%d" % (atRes,len(self.seqRes))
                        task.updateStatus('Calculating backbone score - %d/%d' % (atRes,len(self.seqRes)))

                        #print "."

                    atRes += 1

                    if score != None :
                        zscores2.append ( score )

                    r.bbZ = score
                    #r.CCS = ccs
                    self.scores2[scoreI] = score
                    scoreI += 1


        #print zscores2

        print " - %d res, min %.2f max %.2f, avg %.2f" % (len(ress), min(zscores2), max(zscores2), numpy.average(zscores2) )
        self.avgScore2 = numpy.average ( zscores2 )

        doRes = []

        doAllResInMol = False

        if doAllResInMol :
            for res in self.cur_mol.residues :
                if "CA" in res.atomsMap and "N" in res.atomsMap and "C" in res.atomsMap :
                    doRes.append ( res )
            
            print "++++ added all %d res from %s ++++" % (len(doRes), self.cur_mol.name)

        else :
            for r in self.seqRes :
                try :
                    blah
                    ra = r.scZ
                except :
                    doRes.append ( r )
	


        #doRes = self.seqRes
        #doRes = self.CurRes()
        print " - need score for %d res" % len(doRes)

        umsg ( "Calculating Side Chains / Bases Z-scores..." )

        sczScores = []
        if len(doRes) > 0 :
            sczScores = evalSC ( self.cur_dmap, self.cur_mol, doRes, None, self.whichScore.get(), self.bbCC, self.bbAvgD, task )
            #avgA, stdA = numpy.average ( A ), numpy.std ( A )
            #umsg ( "Avg side chain Z-score: %.3f" % ( avgA ) )

        if not doAllResInMol :
            doRes = self.seqRes 

        self.scores = [None] * len(doRes)
        for ri, r in enumerate ( doRes ) :
            self.scores[ri] = r.scZ
            
        scores = [x for x in self.scores if x is not None]

        self.minScore = min ( scores )
        self.maxScore = max ( scores )
        self.avgScore = numpy.average ( scores )

        print " - %d res, min %.2f max %.2f, avg %.2f" % (len(doRes),self.minScore,self.maxScore, self.avgScore)

        score = None
        if self.whichScore.get() == "Z" or self.whichScore.get() == "Q" :
            self.minSCscore, self.maxSCscore = 0,2
            self.minBBscore, self.maxBBscore = 0,4
        else :
            self.minSCscore, self.maxSCscore = 0,1
            self.minBBscore, self.maxBBscore = 0,1
        

        bbRes = numpy.power ( numpy.e, (self.avgScore2 - 8.0334) / -4.128 ) # y = -4.128ln(x) + 8.0334 
        scRes = numpy.power ( numpy.e, (self.avgScore - 4.8261) / -3.097 ) # y = -3.097ln(x) + 4.8261
        #scRes = (self.avgScore2 - 3.507) / -0.721
        #bbRes = (self.avgScore - 6.1234) / -0.9191

        umsg ( "Average BB Z-score: %.2f (%.1fA), Average Side Chain Z-score: %.2f (%.1fA)" % (self.avgScore2, bbRes, self.avgScore, scRes) )

        self.UpdateSeq ()



        sByType = {}
        rByType = {}
        for r in doRes :
            if r.scZ != None :
                if not r.type in sByType :
                    rByType[r.type] = []
                    sByType[r.type] = []
                rByType[r.type].append ( [r.scZ, r] )
                sByType[r.type].append ( [r.scZ] )

        avgs = []
        for rtype, ra in sByType.iteritems () :
            avgs.append ( [numpy.average (ra), rtype] )

        from chimera.resCode import protein3to1
        from chimera.resCode import nucleic3to1
        avgs.sort ( reverse=True, key=lambda x: x[0] )


        #mpath, mname = os.path.split ( dmap.data.path )
        dname, dext = os.path.splitext ( self.cur_dmap.data.path )
        #mfname = os.path.split ( self.cur_mol.openedAs[0] )[-1]
        #mname, mext = os.path.splitext ( mfname )

        avgm, numt = {}, {}
        for avgScore, rtype in avgs :

            rscores = rByType[rtype]
            rscores.sort ( reverse=True, key=lambda x: x[0] )
            hr = rscores[0]
            R = hr[1]
            highestScore = hr[0]
            numRes = len(rscores)
            
            if R.isProt :
                print "%s\t%s\t%d\t%f\t%d\t.%s\t%f" % (rtype, protein3to1[rtype], numRes, avgScore, R.id.position, R.id.chainId, highestScore)
            else :
                print "%s\t%s\t%d\t%f\t%d\t.%s\t%f" % (rtype, nucleic3to1[rtype], numRes, avgScore, R.id.position, R.id.chainId, highestScore)

            avgm[rtype] = avgScore
            numt[rtype] = numRes
        

        if 0 :
            ofname = "%s__%s__scz_rtype.txt" % (dname, self.cur_mol.name)
            print " -> ", ofname
            fp = open ( ofname, "w" )
    
            for rt in ["PHE", "PRO", "ILE", "LEU", "VAL"] : # , "GLY", , "ALA"
                fp.write ( "%s\t%d\t%f\n" % (rt, numt[rt], avgm[rt]) )
    
            for rt in ["MET"] :
                fp.write ( "%s\t%d\t%f\n" % (rt, numt[rt], avgm[rt]) )
    
            for rt in ["HIS", "ARG", "LYS", "TRP", "CYS"] : # 
                try :
                    fp.write ( "%s\t%d\t%f\n" % (rt, numt[rt], avgm[rt]) )
                except :
                    print " - no %s" % rt
    
            for rt in ["GLN", "ASN", "THR"] :
                fp.write ( "%s\t%d\t%f\n" % (rt, numt[rt], avgm[rt]) )
    
            for rt in ["TYR", "GLU", "ASP", "SER"] :
                fp.write ( "%s\t%d\t%f\n" % (rt, numt[rt], avgm[rt]) )
    
            fp.close()
    
    
    def RtypeOut ( self, avgScore, rtype, rByType, fout ) :
        pass
        
        


    def UpdateSeqFont ( self ) :
        # http://stackoverflow.com/questions/4296249/how-do-i-convert-a-hex-triplet-to-an-rgb-tuple-and-back

        if not hasattr ( self, 'seq' ) :
            print " - update seq font - no seq"
            return

        #print "seq len %d, text w %d" % ( len(self.seq), self.tw )

        # boxes for BBs
        x_at = self.seqX
        y_at = self.seqY + self.seqH/2

        y0 = self.seqY+5
        y1 = self.seqY+self.seqH-5

        for si in range ( len(self.seq) ) :
            res = self.seq[si]
            pred = self.pred[si]
            conf = float ( self.conf[si] ) / 10.0

            if pred == 'E' :
                x0 = self.seqX + si * self.tw
                x1 = x0 + self.tw
                #self.Canvas.coords ( self.seqMouseR, x0, y0, x1, y1 )
                #self.Canvas.itemconfigure ( self.seqMouseR, state=Tkinter.NORMAL )

                if self.seqSheetR[si] == None :
                    c = self.sheetBaseClr + self.sheetClrD * conf
                    clr = "#" + struct.pack('BBB',c[0],c[1],c[2]).encode('hex')
                    self.seqSheetR[si] = self.Canvas.create_rectangle(x0, y0, x1, y1, outline=clr, fill=clr)
                else :
                    self.Canvas.coords ( self.seqSheetR[si], x0, y0, x1, y1 )

            elif pred == 'H' :
                x0 = self.seqX + si * self.tw
                x1 = x0 + self.tw

                if self.seqHelixR[si] == None :
                    c = self.helixBaseClr + self.helixClrD * conf
                    clr = "#" + struct.pack('BBB',c[0],c[1],c[2]).encode('hex')
                    self.seqHelixR[si] = self.Canvas.create_rectangle(x0, y0, x1, y1, outline=clr, fill=clr)
                else :
                    self.Canvas.coords ( self.seqHelixR[si], x0, y0, x1, y1 )



        # box showing selected Residue
        if hasattr ( self, 'seqMouseR' ) :
            self.Canvas.coords ( self.seqMouseR, 0, 0, 0, 0 )
        else :
            self.seqMouseR = self.Canvas.create_rectangle(0, 0, 0, 0, outline="#aab", fill="#bbc", state=Tkinter.HIDDEN)



        x_at = self.seqX
        y_at = self.seqY + self.seqH/2

        if hasattr ( self, 'seqText' ) :
            self.Canvas.coords ( self.seqText, x_at, y_at )
            self.Canvas.itemconfigure ( self.seqText, font=self.font )
        else :
            self.seqText = self.Canvas.create_text( x_at, y_at, text=self.seq, font=self.font, anchor='w')


        #self.UpdateSeqSel ()




    def UpdateSeq ( self ) :
        
        if not hasattr ( self, 'seq' ) :
            print " - update seq - no seq"
            return

        x_at = self.seqX
        y_at = self.seqY + self.seqH/2
        
        if hasattr ( self, 'seqText' ) :
            self.Canvas.coords ( self.seqText, x_at, y_at )
        else :
            self.seqText = self.Canvas.create_text( x_at, y_at, text=self.seq, font=self.font, anchor='w')

        if 1 :
            y0 = self.seqY+5
            y1 = self.seqY+self.seqH-5
            
            cH = numpy.array( [50,200,50] )
            cL = numpy.array( [200,50,50] )

            for si in range ( len(self.seq) ) :
                #if i >= len ( self.seqt ) :
                #    t = self.Canvas.create_text( x_at, y_at, text=self.seq[i], font=self.font)
                #    self.seqt.append ( t )
                #else :
                #    t = self.seqt [ i ]
                #    self.Canvas.coords ( t, x_at, y_at )
                # x_at += self.tw
                    
                pred = self.pred[si]
                if pred == 'E' :
                    if self.seqSheetR[si] != None :
                        x0 = self.seqX + si * self.tw
                        x1 = x0 + self.tw
                        self.Canvas.coords ( self.seqSheetR[si], x0, y0, x1, y1 )

                elif pred == 'H' :
                    if self.seqHelixR[si] != None :
                        x0 = self.seqX + si * self.tw
                        x1 = x0 + self.tw
                        self.Canvas.coords ( self.seqHelixR[si], x0, y0, x1, y1 )

                sc = None
                try :
                    sc = self.scores[si]
                except :
                    #continue
                    pass
                
                if sc == None :
                    if self.seqScoreR[si] != None :
                        self.Canvas.delete ( self.seqScoreR[si] )
                    self.seqScoreR[si] = None
                else :
                    xx0 = self.seqX + si * self.tw + 2
                    xx1 = xx0 + self.tw - 2
                    h = (sc - self.minSCscore) / (self.maxSCscore - self.minSCscore)
                    if h > 1 : h = 1
                    if h < 0 : h = 0
                    Y, H = self.modY, (self.modH/2 - 2)
                    yy0, yy1 = numpy.ceil(Y+H - H*h), numpy.floor(Y+H)
                    if self.seqScoreR[si] != None :
                        self.Canvas.coords ( self.seqScoreR[si], xx0, yy0, xx1, yy1 )
                    else :
                        #c = self.helixBaseClr + self.helixClrD * conf
                        c = h * cH + (1-h) * cL
                        clr = "#" + struct.pack('BBB',c[0],c[1],c[2]).encode('hex')
                        self.seqScoreR[si] = self.Canvas.create_rectangle(xx0, yy0, xx1, yy1, outline=clr, fill=clr)
                        
                bb = None
                try :
                    bb = self.scores2[si]
                except :
                    #continue
                    pass

                if bb == None :
                    if self.seqScoreR2[si] != None :
                        self.Canvas.delete ( self.seqScoreR2[si] )
                    self.seqScoreR2[si] = None
                else :
                    xx0 = self.seqX + si * self.tw + 2
                    xx1 = xx0 + self.tw - 2
                    h = (bb - self.minBBscore) / (self.maxBBscore - self.minBBscore)
                    if h > 1 : h = 1
                    if h < 0 : h = 0
                    Y, H = self.modY, self.modH/2
                    #yy0, yy1 = Y+H, Y+H+H*h #upside down chart
                    yy0, yy1 = numpy.ceil(Y+H+H-H*h), numpy.floor(Y+H+H)
                    if self.seqScoreR2[si] != None :
                        self.Canvas.coords ( self.seqScoreR2[si], xx0, yy0, xx1, yy1 )
                    else :
                        #c = self.helixBaseClr + self.helixClrD * conf
                        c = h * cH + (1-h) * cL
                        clr = "#" + struct.pack('BBB',c[0],c[1],c[2]).encode('hex')
                        self.seqScoreR2[si] = self.Canvas.create_rectangle(xx0, yy0, xx1, yy1, outline=clr, fill=clr)


            self.UpdateSeqSel ()




    def SeqRec ( self, sel ) :
        y0 = self.seqY+5
        y1 = self.seqY+self.seqH-5

        x0 = self.seqX + sel[0] * self.tw
        x1 = self.seqX + (sel[1]+1) * self.tw
        
        return x0, y0, x1, y1


    def UpdateSeqSel ( self ) :
        
        if not hasattr ( self, 'seqSel' ) :
            return
        
        if self.seqSel == None :
            if hasattr(self, 'seqSelRect') :
                self.Canvas.delete ( self.seqSelRect )
                self.seqSelRect = None
            return

        x0, y0, x1, y1 = self.SeqRec ( self.seqSel )

        if hasattr(self, 'seqSelRect') and self.seqSelRect != None :
            self.Canvas.coords ( self.seqSelRect, x0, y0, x1, y1  )
        else :
            #c = self.helixBaseClr + self.helixClrD * conf
            #clr = "#" + struct.pack('BBB',c[0],c[1],c[2]).encode('hex')
            self.seqSelRect = self.Canvas.create_rectangle(x0, y0, x1, y1, outline=self.selColor,  width=3)





    def B1_Down (self, event):
        self.drag = ''

        #print "b1 _", event.x, event.y
        if self.isInSeq ( event.x, event.y ) :
            self.drag = 'seq'
        self.last_x = event.x
        self.last_y = event.y


    def B1_Down_Ctrl ( self, event ) :
        #print "b1 _ <ctrl>", event.x, event.y
        self.drag = ''
        
        if self.isInSeq ( event.x, event.y ) :
            self.drag = 'seqSel'

            if hasattr ( self, 'seqSel' ) and self.seqSel != None :
                self.prevSeqSel = self.seqSel
            else :
                self.prevSeqSel = None

            #print "sel seq..."
            seqI = ( event.x - self.seqX ) / self.tw
            status ( "Start sequence sel at %d" % (seqI+1) )
            self.seqSel = [seqI, seqI]
            self.UpdateSeqSel ()
            
        self.last_x = event.x
        self.last_y = event.y


    def B1_Down_Shift ( self, event ) :
        print "B1 down - shift"

        self.drag = ''

        if self.isInSeq ( event.x, event.y ) :
            if hasattr ( self, 'seqSel' ) and self.seqSel != None :
                seqI = ( event.x - self.seqX ) / self.tw
                if seqI >= self.seqSel[0] and seqI <= self.seqSel[1] :
                    self.drag = "con"
                    if not hasattr ( self, 'conLine' ) or self.conLine == None :
                        self.conLine = self.Canvas.create_line( event.x, event.y, event.x, event.y, fill="red", dash=(1, 1), width=2)
                    status ( "In selected sequence at %d" % seqI )


    def B1_Down_Alt ( self, event ) :
        print "B1 down - alt"

        self.drag = ''

        if self.isInMod ( event.x, event.y ) :
            self.dragMod = self.SelectedMod ( event.x, event.y )
            if self.dragMod != None :
                if self.dragMod.type == "Helix" :
                    self.drag = 'modRot'
                    self.dragStartX = event.x



    def B1_Up_Ctrl ( self, event ) :
        print "b1 up - ctrl - ", event.x, event.y
        self.B1_Up ( event )
        
        
    def B1_Up_Shift ( self, event ) :
        print "b1 up - shift - "
        self.B1_Up ( event )

    def B1_Up_Alt ( self, event ) :
        print "b1 up - alt - "
        self.B1_Up ( event )
        

    def B1_Up (self, event):
        print "b1 up - ", event.x, event.y

        if self.drag == 'seqSel' and hasattr ( self, 'seqSel' ) :
            status ( "Selected: %d-%d" % (self.seqSel[0], self.seqSel[1]) )
        
            if hasattr ( self, 'prevSeqSel' ) and self.prevSeqSel != None :
                if self.seqSel[0] == self.seqSel[1] :
                    self.seqSel = None
                    self.prevSeqSel = None
                    self.UpdateSeqSel ()
                    status ( "Cleared sequence selection" )
                    chimera.selection.clearCurrent ()

            if self.seqSel != None :
                m, cid = self.cur_mol, self.chain.get()
                if m != None :
                    startI = self.seqRes [ max(self.seqSel[0],0) ].id.position
                    endI = self.seqRes [ min(self.seqSel[1],len(self.seqRes)-1) ].id.position
                    selStr = "#%d:%d-%d.%s" % (m.id,startI,endI,cid)
                    
                    self.lastSelStr = "%d-%d.%s" % (startI,endI,cid)
                    
                    if hasattr ( self, 'prevSel' ) and self.preserveSel.get () :
                        for s in self.prevSel :
                            print " -s- adding to sel:", s
                            selStr = selStr + "," + s 
                    else :
                        self.prevSel = []

                    if self.preserveSel.get () :
                        self.prevSel.append ( "%d-%d.%s" % (startI,endI,cid) )
                        print " - added to selection list..."
                    
                    umsg ( "Selected: " + selStr )
                    sel = chimera.selection.OSLSelection ( selStr )
                    chimera.selection.setCurrent ( sel )
                    #chimera.selection.addCurrent ( sel )
                    self.ShowSel ()

                else :
                    status ( "no model visible" )

            #else :
            #    print "cleared past sel"
            #    self.prevSel = []


        elif self.drag == 'modSel' :
            status ( 'Selected %d mods' % len(self.selMods) )

        elif self.drag == 'con' :
            selMod = None
            if hasattr ( self, 'selModPiece' ) and self.selModPiece != None :
                selMod = self.selModPiece
                self.selModPiece = None
            else :
                return

            if hasattr ( self, 'conLine' ) and self.conLine != None :
                self.Canvas.delete ( self.conLine )
                self.conLine = None
    
            status ( "connected to %s" % selMod.type )

            selMod.seq = self.seqSel
            selMod.numRes = (self.seqSel[1] - self.seqSel[0] + 1)
            selMod.MakeMod ()

            self.UpdateMod ()
        
        self.drag = ''
        print "mod: ", self.modX, " seq:", self.seqX


    def preserveSelCb (self) :
        print "Preserve set to ", self.preserveSel.get()      
        if self.preserveSel.get() :
            print " - setting current selection to preserve..."  
            if hasattr ( self, 'lastSelStr' ) :
                self.prevSel = [ self.lastSelStr ]
        else :
            print " - clearing current"  
            self.prevSel = []
            
    #def keepExMapCb (self) :
    #    print "Kep ex map set to ", self.keepExMap.get()      


    def ClearSel ( self ) :
        self.prevSel = []
        self.seqSel = None
        self.prevSeqSel = None
        self.UpdateSeqSel ()
        status ( "Cleared sequence selection" )
        chimera.selection.clearCurrent ()




    def ExCustA ( self ) :

        if self.cur_dmap == None :
            umsg ("Select a map first")
            return

        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return
 
        #selStr = "#%d:80-87.I,171-184.I,227-239.I,80-87.A,171-184.A,227-239.A,80-87.B,171-184.B,227-239.B,80-87.J,171-184.J,227-239.J,80-87.H,171-184.H,227-239.H" % self.cur_mol.id
        selStr = "#%d:80-87.I,171-184.I,227-239.I,80-87.A,171-184.A,227-239.A,80-87.J,171-184.J,227-239.J" % self.cur_mol.id
        
        umsg ( "Selected: " + selStr )
        sel = chimera.selection.OSLSelection ( selStr )
        chimera.selection.setCurrent ( sel )
        self.ShowSel()

 
    def ExCustB ( self ) :

        if self.cur_dmap == None :
            umsg ("Select a map first")
            return

        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return
 
        selStr = "#%d:428-435.F,365-376.F,396-402.F,428-435.I,365-376.I,396-402.I" % self.cur_mol.id
        #selStr = "#%d:428-435.A,365-376.A,396-402.A,428-435.H,365-376.H,396-402.H" % self.cur_mol.id

        
        umsg ( "Selected: " + selStr )
        sel = chimera.selection.OSLSelection ( selStr )
        chimera.selection.setCurrent ( sel )
        self.ShowSel()

    def ExCustC ( self ) :

        if self.cur_dmap == None :
            umsg ("Select a map first")
            return

        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return
 
        #selStr = "#%d:10:548-558.I,520-530.I,548-558.F,520-530.F" % self.cur_mol.id
        selStr = "#%d:428-435.F,365-376.F,396-402.F,428-435.I,365-376.I,396-402.I,548-558.I,520-530.I,548-558.F,520-530.F" % self.cur_mol.id

        
        umsg ( "Selected: " + selStr )
        sel = chimera.selection.OSLSelection ( selStr )
        chimera.selection.setCurrent ( sel )
        self.ShowSel()



 
    def ShowSel ( self ) :
        
        #showRibbon = self.showRibbon.get()
        showRibbon = self.showRibbon.get()
        showSC = self.showAtoms.get()

        atoms = []
        scores = []
        selResM = {}
        for r in chimera.selection.currentResidues () :
            rid = "%d.%s" % (r.id.position, r.id.chainId)
            selResM [rid] = 1

        if self.cur_mol == None :
            return

        if 1 or not hasattr ( self.cur_mol, 'bbats' ) :
            SetBBAts(self.cur_mol)
            self.cur_mol.bbats = True


        for r in self.cur_mol.residues :
            rid = "%d.%s" % (r.id.position, r.id.chainId)
            if rid in selResM :

                if hasattr (r, 'scZ') and r.scZ != None :
                    scores.append(r.scZ)

                r.ribbonDisplay = showRibbon

                for at in r.atoms :
                    if at.element.name == "H" :
                        at.display = False
                    elif at.isSC :
                        if showSC :
                            at.display = True
                            atoms.append ( at )
                        else :
                            at.display = False
                    else :
                        at.display = True
                        atoms.append ( at )
                    if at.element.name in atomColors :
                        if at.isBB :
                            at.color = atomColors[at.element.name]
                            if at.element.name == "C" :
                                at.color = atomColors['Cbb']
                        else :
                            at.color = atomColors[at.element.name]

            else :
                r.ribbonDisplay = False
                for at in r.atoms :
                    at.display = False


        #for bond in self.seqRes[0].molecule.bonds :
        #    bond.display = bond.Smart
            #if bond.atoms[0] in atMap and bond.atoms[1] in atMap :
            #    #bond.display = bond.Smart
            #    bond.display = bond.Smart
            #else :
            #    #bond.display = bond.Never
            #    bond.display = bond.Smart
            

        if len(atoms) > 0 :

            from _multiscale import get_atom_coordinates
            points = get_atom_coordinates ( atoms, transformed = True )
            COM, U, S, V = prAxes ( points )

            if 1 :            
                p0 = numpy.array ( chimera.viewer.camera.center )
                p1 = numpy.array ( [ COM[0], COM[1], COM[2] ] )
                for i in range (10) :
                    f = float(i) / 9.0
                    f1, f2 = 2.0*f*f*f-3.0*f*f+1.0, 3*f*f-2*f*f*f
                    P = p0 * f1 + p1 * f2
                    chimera.viewer.camera.center = (P[0],P[1],P[2])
                    print ".",
                print ""


            atomRad = 2.0 # float ( self.maskWithSelDist.get() )
            print " - %d selected atoms, mask at %.2f" % ( len(atoms), atomRad )
            dmap = self.cur_dmap

            mlist = OML(modelTypes = [VolumeViewer.volume.Volume])

            for m in mlist :
                if "sel_masked" in m.name :
                    chimera.openModels.close ( [m] )

            if len ( atoms ) > 0 and dmap != None :
                #points = get_atom_coordinates ( atoms, transformed = False )
                self.PtsToMap ( points, dmap, atomRad, dmap.name + " sel_masked", False )
                if self.showMesh.get () :
                    self.PtsToMap ( points, dmap, atomRad, dmap.name + " sel_masked_mesh", True )

            if len(scores) > 0 :
                umsg ( "%d residue scores, avg score %.1f" % ( len(scores), numpy.average(scores) ) )

        else :
            umsg ( "no atoms selected, try reselecting the model and chain..." )



 
 

    def PtsToMap0 ( self, points, dmap, atomRad, nname, neg = 1.0 ) :
        import _contour
        _contour.affine_transform_vertices ( points, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
        mdata = VolumeData.zone_masked_grid_data ( dmap.data, points, atomRad )

        gdata = VolumeData.Array_Grid_Data ( mdata.full_matrix()*neg, dmap.data.origin, dmap.data.step, dmap.data.cell_angles, name = nname )
        nv = VolumeViewer.volume.volume_from_grid_data ( gdata )
        nv.name = nname
        dmap.display = False
        nv.region = ( nv.region[0], nv.region[1], [1,1,1] )
        nv.surface_levels[0] = dmap.surface_levels[0]
        ro = VolumeViewer.volume.Rendering_Options()
        #ro.smoothing_factor = .5
        #ro.smoothing_iterations = 20
        #ro.surface_smoothing = True
        nv.update_surface ( False, ro )
        for sp in nv.surfacePieces :
            v, t = sp.geometry
            if len(v) == 8 and len(t) == 12 :
                sp.display = False
            else :
                sp.color = (0.7, 0.7, 0.7, 0.3)


    def PtsToMap ( self, points, dmap, atomRad, nname, showMesh = False ) :

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

        print " - bounds - %d %d %d --> %d %d %d --> %d %d %d" % ( li,lj,lk, hi,hj,hk, n1,n2,n3 )

        #nmat = numpy.zeros ( (n1,n2,n3), numpy.float32 )
        #dmat = dmap.full_matrix()

        nstep = (dmap.data.step[0], dmap.data.step[1], dmap.data.step[2] )
        #nstep = (fmap.data.step[0]/2.0, fmap.data.step[1]/2.0, fmap.data.step[2]/2.0 )

        nn1 = int ( round (dmap.data.step[0] * float(n1) / nstep[0]) )
        nn2 = int ( round (dmap.data.step[1] * float(n2) / nstep[1]) )
        nn3 = int ( round (dmap.data.step[2] * float(n3) / nstep[2]) )

        O = dmap.data.origin
        print " - %s origin:" % dmap.name, O
        nO = ( O[0] + float(li) * dmap.data.step[0],
               O[1] + float(lj) * dmap.data.step[1],
               O[2] + float(lk) * dmap.data.step[2] )

        print " - new map origin:", nO

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
        gdata = VolumeData.Array_Grid_Data ( mdata.full_matrix(), nO, nstep, dmap.data.cell_angles, name = "atom masked" )
        nv = VolumeViewer.volume.volume_from_grid_data ( gdata )
        nv.openState.xform = dmap.openState.xform
        
        nv.name = nname
        dmap.display = False
        nv.region = ( nv.region[0], nv.region[1], [1,1,1] )
        nv.surface_levels[0] = dmap.surface_levels[0]
        ro = VolumeViewer.volume.Rendering_Options()
        ro.smoothing_factor = .3
        ro.smoothing_iterations = 2
        ro.surface_smoothing = False
        ro.square_mesh = True
        ro.line_thickness = 2
        nv.update_surface ( False, ro )
        setro (ro)
        for sp in nv.surfacePieces :
            v, t = sp.geometry
            if len(v) == 8 and len(t) == 12 :
                sp.display = False
            else :
                if showMesh :
                    sp.color = (.5, .5, .5, 1.0)
                    sp.displayStyle = sp.Mesh
                else :
                    sp.color = (0.7, 0.7, 0.7, 0.3)


    def B1_Drag (self, event):
        #print "b1m ", event.x, event.y

        if self.drag == 'seq' :
            d = event.x - self.last_x
            self.seqX += d
            #GetSegMod().seqX = self.seqX
            self.UpdateSeq ()
        elif self.drag == 'mod' :
            d = event.x - self.last_x
            self.modX += d
            #GetSegMod().modX = self.modX
            self.UpdateMod ()
        elif self.drag == 'seqSel' :
            if hasattr ( self, 'seqSel' ):
                seqI = ( event.x - self.seqX ) / self.tw
                if seqI > self.seqSel[0] :
                    self.seqSel[1] = seqI
                elif seqI < self.seqSel[1] :
                    self.seqSel[0] = seqI
                status ( "Sequence selected %d - %d" % (self.seqSel[0]+1, self.seqSel[1]+1) )
                self.UpdateSeqSel ()
        elif self.drag == 'con' :
            x1, y1, x2, y2 = self.Canvas.coords ( self.conLine )
            self.Canvas.coords ( self.conLine, x1, y1, event.x, event.y )
            self.SelectedModClr ( event.x, event.y )
        elif self.drag == "modRot" :
            dx = event.x - self.dragStartX
            self.dragStartX = event.x
            self.dragMod.Rotate ( dx )



        self.last_x = event.x
        self.last_y = event.y


    def B2_Down (self, event):
        print "b2 - down"



    
    def B2_Up (self, event):
        print "b2  - up", event.x, event.y
    
        if hasattr ( self, 'selModPiece' ) and self.selModPiece != None :
            
            if self.selModPiece.type == "Loop" :
                self.selModPiece.MakeMod ()
                
            else :
                self.selModPiece.switch = not self.selModPiece.switch
                self.selModPiece.MakeMod ()
                self.UpdateMod ()
    

    def B2_Up_Ctrl (self, event):
        print "b2  - up - control", event.x, event.y
        if hasattr ( self, 'selModPiece' ) and self.selModPiece != None :
            if self.selModPiece.type == "Loop" :
                MakeLoopMod1 ( self.selModPiece )
                #MakeLoopMod ( self.selModPiece )



    def B2_Up_Alt (self, event):
        print "b2  - up - alt", event.x, event.y
        if hasattr ( self, 'selModPiece' ) and self.selModPiece != None :
            if self.selModPiece.type == "Loop" :
                LoopPathOpt ( self.selModPiece, self.refUseMap.get() )
                

    def B2_Up_Shift (self, event):
        print "b2  - up - alt", event.x, event.y
        if hasattr ( self, 'selModPiece' ) and self.selModPiece != None :
            if self.selModPiece.type == "Loop" :
                LoopPathOpt ( self.selModPiece, self.refUseMap.get() )



    def B2_Up_Comm (self, event):
        print "b2  - up - command", event.x, event.y
    



    def B2_Drag (self, event):
        #print "b2m ", event.x, event.y
        pass



    def B3_Down (self, event):

        print "b3 _", event.x, event.y
        


    
    def B3_Up (self, event):
        print "b3 ^", event.x, event.y
        self.B2_Up ( event )


    def B3_Drag (self, event):
        #print "b3m ", event.x, event.y
        pass

    
    def isInSeq ( self, x, y ) :
        if y >= self.seqY and y <= self.seqY + self.seqH :
            return True
        else :
            return False
    
    def isInMod ( self, x, y ) :
        if y >= self.modY and y <= self.modY + self.modH :
            return True
        else :
            return False


    def Mouse_Move (self, event):
        #print "mod m ", event.x, event.y
        #self.Canvas.coords ( self.seqMouseLine, event.x,self.seqY,event.x,self.seqY+self.seqH )

        if self.isInSeq ( event.x, event.y ) and hasattr ( self, 'seq') and len(self.seq) > 0 :

            if hasattr ( self, 'seqRec' ) and hasattr ( self, 'tw' ) and hasattr ( self, 'seqMouseR' ) :
                self.Canvas.itemconfigure ( self.seqRec, state=Tkinter.NORMAL )

                si = ( event.x - self.seqX ) / self.tw
                if si < 0 :
                    si = 0
                if si < len ( self.seq ) :
                    res = self.seqRes [ si ]
                    resEnd = self.seqRes [ len(self.seqRes) - 1 ]
                    
                    try :
                        status ( "Sequence: %s/%s %d/%d" % ( self.seq[si], res.type, res.id.position, resEnd.id.position ) )
                    except :
                        status ( "model not found" )
                        self.chain.set("")
                        self.struc.set("")
                        self.RemoveSeq ()
                        return

                    y0 = self.seqY+5
                    y1 = self.seqY+self.seqH-5
                    if event.y >= y0 and event.y <= y1 and hasattr ( self, 'seqMouseR' ) :
                        x0 = self.seqX + si * self.tw
                        x1 = x0 + self.tw
                        self.Canvas.coords ( self.seqMouseR, x0, y0, x1, y1 )
                        self.Canvas.itemconfigure ( self.seqMouseR, state=Tkinter.NORMAL )
                    else :
                        self.Canvas.itemconfigure ( self.seqMouseR, state=Tkinter.HIDDEN )

            else :
                self.Canvas.itemconfigure ( self.seqRec, state=Tkinter.HIDDEN )

                if hasattr ( self, 'seqMouseR' ) :
                    self.Canvas.itemconfigure ( self.seqMouseR, state=Tkinter.HIDDEN )


        self.last_x = event.x
        self.last_y = event.y


    def Canvas_Leave ( self, event ) :
        #self.Canvas.coords ( self.seqMouseLine, 0,0,0,0 )
        pass


    def Canvas_Config (self, event) :
        #print "mod cfg ", event.width, event.height
        self.W = event.width
        self.H = event.height

        #self.Canvas.delete("all")
        if 1 :
            if hasattr(self, 'backRec') :
                self.Canvas.coords (self.backRec, 0, 0, self.W, self.H)
            else :
                self.backRec = self.Canvas.create_rectangle(0, 0, self.W, self.H, outline="#eee", fill="#eee")
                #self.seqMouseLine = self.Canvas.create_line(0, 0, 0, 0, fill="#66a")

            if hasattr ( self, 'seqRec' ) :
                self.Canvas.coords ( self.seqRec, 0, self.seqY, self.W, self.seqY+self.seqH )
            else :
                self.seqRec = self.Canvas.create_rectangle(0, self.seqY, self.W, self.seqY+self.seqH, outline="#ddd", fill="#ddd" )

            self.Canvas.tag_lower(self.seqRec)
            self.Canvas.tag_lower(self.backRec)


    def Canvas_Wheel ( self, event ) :

        if self.isInSeq (self.last_x, self.last_y) :
            
            #self.seqX += event.delta * 10

            self.mag = self.mag + event.delta
            if self.mag > 15 : self.mag = 15
            if self.mag < 2 : self.mag = 2

            self.font = tkFont.Font(family='Courier', size=(self.mag), weight='normal')
            #self.boldFont = tkFont.Font(family='Courier', size=(self.mag+4), weight='bold')
            self.tw = self.font.measure ( "a" )

            #GetSegMod().seqX = self.seqX
            self.UpdateSeqFont ()
            self.UpdateSeq ()
            
            # ['__doc__', '__module__', 'char', 'delta', 'height', 'keycode', 'keysym', 'keysym_num', 'num', 'send_event', 'serial', 'state', 'time', 'type', 'widget', 'width', 'x', 'x_root', 'y', 'y_root']
            #print dir(event)
            #print event.delta
            status ( "Mag: %d" % self.mag )
    


    def ZoomMinus ( self ) :
        self.mag = self.mag - 1
        if self.mag > 15 : self.mag = 15
        if self.mag < 2 : self.mag = 2
        #print "w ", event.delta, " mag: ", self.mag

        self.font = tkFont.Font(family='Courier', size=(self.mag), weight='normal')
        #self.boldFont = tkFont.Font(family='Courier', size=(self.mag+4), weight='bold')
        self.tw = self.font.measure ( "a" )

        self.UpdateSeqFont ()
        self.UpdateSeq ()
        status ( "Magnification: %d" % self.mag )



    def ZoomPlus ( self ) :
        self.mag = self.mag + 1
        if self.mag > 15 : self.mag = 15
        if self.mag < 2 : self.mag = 2
        #print "w ", event.delta, " mag: ", self.mag

        self.font = tkFont.Font(family='Courier', size=(self.mag), weight='normal')
        #self.boldFont = tkFont.Font(family='Courier', size=(self.mag+4), weight='bold')
        self.tw = self.font.measure ( "a" )

        self.UpdateSeqFont ()
        #self.UpdateSeq ()
        self.UpdateSeq ()

        status ( "Magnification: %d" % self.mag )


    def ZoomBegin ( self ) :        
        self.seqX = 10
        self.UpdateSeq ()

    def ZoomEnd ( self ) :
        self.seqX = - ( len(self.seq) - 50 ) * self.tw
        self.UpdateSeq ()




    def isSelected ( self, fmap ) :
        for sp in fmap.surfacePieces :
            if sp in Surface.selected_surface_pieces() :
                return True
        return False




    def RadScore (self) :
        
        selRes = chimera.selection.currentResidues()
        if len ( selRes ) == 0 :
            return
            
        dmap = self.cur_dmap
        

        r = selRes[0]
        print "Res: %s - %d.%s - %s" % (r.type, r.id.position, r.id.chainId, r.molecule.name)

        if not hasattr ( r.molecule, 'bbats' ) :
            SetBBAts(r.molecule)
            r.molecule.bbats = True

        ptsMol = GetMod ( "SD points" )
        if ptsMol :
            chimera.openModels.remove ( [ptsMol] )


        points = _multiscale.get_atom_coordinates ( r.molecule.atoms, transformed = False )
        print " - tree with %d ats" % len(r.molecule.atoms)
        allAtTree = AdaptiveTree ( points.tolist(), r.molecule.atoms, 1.0)

        #allAtTree = None

        r.Sd = SdevRes ( selRes, dmap, allAts=1, bbAts=0, scAts=0, allAtTree=allAtTree, show=True, toRAD=5, dRAD=0.2 )     



    def RadScoreBB (self) :
        
        selRes = chimera.selection.currentResidues()
        if len ( selRes ) == 0 :
            return
            
        dmap = self.cur_dmap
        

        r = selRes[0]
        print "Res: %s - %d.%s - %s" % (r.type, r.id.position, r.id.chainId, r.molecule.name)

        if not hasattr ( r.molecule, 'bbats' ) :
            SetBBAts(r.molecule)
            r.molecule.bbats = True
            
        ptsMol = GetMod ( "SD points" )
        if ptsMol :
            chimera.openModels.remove ( [ptsMol] )
            
    
        points = _multiscale.get_atom_coordinates ( r.molecule.atoms, transformed = False )
        print " - tree with %d ats" % len(r.molecule.atoms)
        allAtTree = AdaptiveTree ( points.tolist(), r.molecule.atoms, 1.0)
        
        #allAtTree = None
        
        r.Sd = SdevRes ( selRes, dmap, allAts=0, bbAts=1, scAts=0, allAtTree=allAtTree, show=True, toRAD=5, dRAD=0.2 )     
    


    def RadScores (self) :
        
        ress = []
        try :
            ress = self.seqRes
        except :
            pass
            
        if len ( ress ) == 0 :
            umsg ( "No molecule/chain selected?" )
            return

       
        ok = True
        try :
            print self.cur_dmap.name
        except :
            status ( "Selected map not found; please choose another map" )
            self.dmap.set ("")
            ok = False
            
        try :
            print self.cur_mol.name
        except :
            status ( "Selected model not found; please choose another model" )
            self.struc.set ("")
            self.chain.set ("")
            self.RemoveSeq ()
            ok = False
        
        if not ok :
            return


        doRes = []

        doAllResInMol = False

        if doAllResInMol :
            for res in self.cur_mol.residues :
                doRes.append ( res )
            print "-- added all %d res from %s --" % (len(doRes), self.cur_mol.name)
        
        else :
            for r in self.seqRes :
                try :
                    blah
                    ra = r.scZ
                except :
                    doRes.append ( r )
	
        print " - need score for %d res" % len(doRes)

        umsg ( "Calculating sdev scores" )


        points = _multiscale.get_atom_coordinates ( self.cur_mol.atoms, transformed = False )
        print " - tree with %d ats" % len(self.cur_mol.atoms)
        allAtTree = AdaptiveTree ( points.tolist(), self.cur_mol.atoms, 1.0)
        #allAtTree = None
        
        scores, scoresBB, scoresSC = [], [], []
        for ri, r in enumerate ( ress ) :
            #r.sd = SdevRes ( [r], self.cur_dmap, 1, 0, 0, allAtTree, show=0, toRAD=5, dRAD=0.2 )
            r.scZ = SdevRes ( [r], self.cur_dmap, 0, 0, 1, allAtTree, show=0, toRAD=5, dRAD=0.2 )
            r.bbZ = SdevRes ( [r], self.cur_dmap, 0, 1, 0, allAtTree, show=0, toRAD=5, dRAD=0.2 )

            #r.Sd = SdevRes ( selRes, dmap, allAts=0, bbAts=1, scAts=0, allAtTree=allAtTree, show=True, toRAD=5, dRAD=0.2 )     

            #scores.append ( r.sd )
            scoresBB.append ( r.bbZ )
            scoresSC.append ( r.scZ )
            
            if ri % 10 == 0 :
                status ( "Calculating sdev scores - res %d/%d" % (ri, len(ress)) )
                print ".",

        print ""

        self.scores, self.scores2 = scoresSC, scoresBB

        #sc = [x for x in scores if x is not None]
        scSC = [x for x in scoresSC if x is not None]
        scBB = [x for x in scoresBB if x is not None]

        print " - %d res, min %.2f max %.2f, avg %.2f" % (len(doRes), min(sc), max(sc), numpy.average(sc))
        print " - avg %.1f, side chain %.1f, backbone %.1f" % (numpy.average(scores), numpy.average(scSC), numpy.average(scBB) )


        #bbRes = numpy.power ( numpy.e, (self.avgScore2 - 8.0334) / -4.128 ) # y = -4.128ln(x) + 8.0334 
        #scRes = numpy.power ( numpy.e, (self.avgScore - 4.8261) / -3.097 ) # y = -3.097ln(x) + 4.8261
        #scRes = (self.avgScore2 - 3.507) / -0.721
        #bbRes = (self.avgScore - 6.1234) / -0.9191

        #umsg ( "Average BB Z-score: %.2f (%.1fA), Average Side Chain Z-score: %.2f (%.1fA)" % (self.avgScore2, bbRes, self.avgScore, scRes) )

        self.UpdateSeq ()



        sByType = {}
        rByType = {}
        for r in doRes :
            if r.scZ != None :
                if not r.type in sByType :
                    rByType[r.type] = []
                    sByType[r.type] = []
                rByType[r.type].append ( [r.scZ, r] )
                sByType[r.type].append ( [r.scZ] )

        avgs = []
        for rtype, ra in sByType.iteritems () :
            avgs.append ( [numpy.average (ra), rtype] )

        from chimera.resCode import protein3to1
        from chimera.resCode import nucleic3to1
        avgs.sort ( reverse=True, key=lambda x: x[0] )
        for avgScore, rtype in avgs :

            rscores = rByType[rtype]
            rscores.sort ( reverse=True, key=lambda x: x[0] )
            hr = rscores[0]
            R = hr[1]
            highestScore = hr[0]
            numRes = len(rscores)
            
            if R.isProt :
                print "%s\t%s\t%d\t%f\t%d\t.%s\t%f" % (rtype, protein3to1[rtype], numRes, avgScore, R.id.position, R.id.chainId, highestScore)
            else :
                print "%s\t%s\t%d\t%f\t%d\t.%s\t%f" % (rtype, nucleic3to1[rtype], numRes, avgScore, R.id.position, R.id.chainId, highestScore)
        




    def AlignRes1 ( self ) :
        
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        chainId = self.chain.get()
        if len(chainId) == 0 :
            umsg ("Select a chain first")
            return

        if self.cur_dmap == None :
            umsg ("Select a map first")
            return

            
        #SetBBAts ( self.cur_mol )
        last_x = 0.0
        last_y = 0.0

        r0, exR0, xtR0 = None, None, None
        
        alAts = []
        if self.exType == "ASP" : alAts = ["CG","OD1","OD2"]
        if self.exType == "GLU" : alAts = ["CD","OE1","OE2"]
        if self.exType == "TYR" : alAts = ["OH","CE1","CE2","CD1","CD2","CG","CB"]
        

        for r in self.cur_mol.residues :
            if r.id.chainId == chainId and r.type == self.exType :
                print " - res %s %d" % (r.type, r.id.position)

                if r0 == None :
                    r0 = r

                    r.exMaps[0].display = True
                    r.exMaps[1].display = False

                    #r.xtMaps[0].display = False
                    #r.xtMaps[1].display = False

                    for at in r.atoms :
                        at.display = at.name in alAts

                else :

                    exR0 = r0.exMol.residues[0]
                    exR = r.exMol.residues[0]
                    ats0, ats = [], []
                    for atName in alAts :
                        ats0.append ( exR0.atomsMap[atName][0] ) 
                        ats.append ( exR.atomsMap[atName][0] ) 
                        
                    for at in r.atoms :
                        at.display = at.name in alAts

                    #aCG0, aOD10, aOD20 = exR0.atomsMap['CG'][0], exR0.atomsMap['OD1'][0], exR0.atomsMap['OD2'][0],
                    #aCG, aOD1, aOD2 = exR.atomsMap['CG'][0], exR.atomsMap['OD1'][0], exR.atomsMap['OD2'][0],

                    #xf, rmsd = chimera.match.matchPositions ( pts_o, pts_c )
                    #xf, rmsd = chimera.match.matchAtoms ( [aCG0, aOD10, aOD20], [aCG, aOD1, aOD2] )
                    xf, rmsd = chimera.match.matchAtoms ( ats0, ats )
                    print " - rmsd: ", rmsd

                    #from _multiscale import get_atom_coordinates
                    #points = get_atom_coordinates ( atoms, transformed = True )

                    #exR.xf0 = r.exMol.openState.xform

                    mxf = r0.exMol.openState.xform
                    mxf.multiply ( xf )
                    r.exMol.openState.xform = mxf
                    r.exMaps[0].openState.xform = mxf
                    r.exMaps[1].openState.xform = mxf
                    r.exMaps[0].display = True
                    r.exMaps[1].display = False
                    
                    #r.xtMaps[0].display = False
                    #r.xtMaps[1].display = False


                    #break


    def AlignRes2 ( self ) :
        
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        chainId = self.chain.get()
        if len(chainId) == 0 :
            umsg ("Select a chain first")
            return

        if self.cur_dmap == None :
            umsg ("Select a map first")
            return

            
        #SetBBAts ( self.cur_mol )
        last_x = 0.0
        last_y = 0.0

        r0, exR0, xtR0 = None, None, None
        

        for r in self.cur_mol.residues :
            if r.id.chainId == chainId and r.type == "ASP" :
                print " - res %s %d" % (r.type, r.id.position)
                
                if r0 == None :
                    r0 = r
                
                    r.exMaps[0].display = False
                    r.exMaps[1].display = False

                    r.xtMaps[0].display = True
                    r.xtMaps[1].display = False


                else :
                    
                    r.exMaps[0].display = False
                    r.exMaps[1].display = False

                    exR0 = r0.xtMol.residues[0]
                    aCB0, aCG0, aOD10, aOD20 = exR0.atomsMap['CB'][0], exR0.atomsMap['CG'][0], exR0.atomsMap['OD1'][0], exR0.atomsMap['OD2'][0],

                    exR = r.xtMol.residues[0]
                    aCB, aCG, aOD1, aOD2 = exR.atomsMap['CB'][0], exR.atomsMap['CG'][0], exR.atomsMap['OD1'][0], exR.atomsMap['OD2'][0],

                    #xf, rmsd = chimera.match.matchPositions ( pts_o, pts_c )
                    xf, rmsd = chimera.match.matchAtoms ( [aCB0, aCG0, aOD10, aOD20], [aCB, aCG, aOD1, aOD2] )
                    print " - rmsd: ", rmsd

                    #from _multiscale import get_atom_coordinates
                    #points = get_atom_coordinates ( atoms, transformed = True )
                    
                    #exR.xf0 = r.exMol.openState.xform

                    mxf = r0.xtMol.openState.xform
                    mxf.multiply ( xf )
                    r.xtMol.openState.xform = mxf
                    r.xtMaps[0].openState.xform = mxf
                    r.xtMaps[1].openState.xform = mxf
                    r.xtMaps[0].display = True
                    r.xtMaps[1].display = False


                    #break

                    

    def Avg ( self ) :

        print " -- finding base map --- "
        largestMap = None
        maxD = 0
        for m in OML(modelTypes = [VolumeViewer.volume.Volume]) :
            if m.display == True :
                d = numpy.sum ( m.data.size )
                if d > maxD :
                    maxD = d
                    largestMap = m
        
        print " - largest map: ", largestMap.name
        dmap = largestMap
        dmap.display = False


        fmap = None
        avgMat = dmap.data.full_matrix()
        N = 0.0

        print " ----------- Averaging... ---------------------"

        for m in OML(modelTypes = [VolumeViewer.volume.Volume]) :
            if m.display == True and m != dmap :
                print m.name

                df_mat = self.Map2Map ( m, dmap )
                m.display = False
                N = N + 1.0
                avgMat = avgMat + df_mat


        print " ----------- n=%f ---------------------" % N

        avgMat = avgMat / N
        df_data = VolumeData.Array_Grid_Data ( avgMat, dmap.data.origin, dmap.data.step, dmap.data.cell_angles, name="avg" )

        MapFromData ( df_data, "Avg", dmap, False )
        MapFromData ( df_data, "Avg", dmap, True )


        #df_v = VolumeViewer.volume.volume_from_grid_data ( df_data )
        #df_v.name = "Avg"
        #df_v.openState.xform = dmap.openState.xform

        #nv = self.ShrinkMap ( df_v, 1e-3 )



    def Map2Map ( self, densitiesFromMap, toGridOfMap, mask = False ) :

        fmap = toGridOfMap
        dmap = densitiesFromMap

        import _contour
        n1, n2, n3 = fmap.data.size[0], fmap.data.size[1], fmap.data.size[2]
        f_points = VolumeData.grid_indices( (n1,n2,n3), numpy.single )  # i,j,k indices
        _contour.affine_transform_vertices( f_points, fmap.data.ijk_to_xyz_transform )

        d_vals = dmap.interpolated_values ( f_points, fmap.openState.xform )
        df_mat = d_vals.reshape( (n3,n2,n1) )

        if mask :
            f_mat = fmap.data.full_matrix()
            f_mask = numpy.where ( f_mat > fmap.surface_levels[0], numpy.ones_like(f_mat), numpy.zeros_like(f_mat) )
            df_mat = df_mat * f_mask

        return df_mat




    def CloseExtracted ( self ) :
        
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        chainId = self.chain.get()
        if len(chainId) == 0 :
            umsg ("Select a chain first")
            return

        if self.cur_dmap == None :
            umsg ("Select a map first")
            return

            
        for r in self.cur_mol.residues :
        
            if hasattr ( r, "exMaps" ) :
                chimera.openModels.close ( r.exMaps ); del r.exMaps

            if hasattr ( r, "xtMaps" ) :
                chimera.openModels.close ( r.xtMaps ); del r.xtMaps
                
            if hasattr ( r, "exMol" ) :
                chimera.openModels.close ( [r.exMol] ); del r.exMol

            if hasattr ( r, "xtMol" ) :
                chimera.openModels.close ( [r.xtMol] ); del r.xtMol
        
        for m in chimera.openModels.list() :
            if m.name == "Avg" or m.name == "Avg_mesh" :
                chimera.openModels.close ( [m] )
            
                
    


    def Extract ( self ) :
        
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        chainId = self.chain.get()
        if len(chainId) == 0 :
            umsg ("Select a chain first")
            return

        if self.cur_dmap == None :
            umsg ("Select a map first")
            return

            
        #SetBBAts ( self.cur_mol )
        last_x = 0.0
        last_y = 0.0


        print "Extracting - %s - %s - %s" % (self.cur_dmap.name, self.cur_mol.name, chainId)
        
        #self.exType = "TYR"
        #self.exType = "GLU"
        self.exType = "ASP"
        
        yzAts = { "ASP" : ["CB","CG","OD1"],
                  "GLU" : ["CG","CD","OE1"],
                  "TYR" : ["CB","CZ","CD1"],
        }

        for r in self.cur_mol.residues :
        
            if r.id.chainId == chainId and r.type == self.exType :
                
                print " - res %s %d" % (r.type, r.id.position)

                self.ExtractRes ( r, self.cur_mol, self.cur_dmap, last_x, last_y, yzAts[self.exType] )

                #self.ExtendRes ( r, self.cur_mol, self.cur_dmap, last_x, -8.0, thrF=0.8 )

                last_x += 7.0
                
                #break

                
                
    def ExtractRes ( self, r, mol, dmap, atX, atY, xyAts ) :

        nmol, nres = CopyRess ( [r] )
        nmol.name = mol.name + "_%s_%d" % (r.type, r.id.position)
        chimera.openModels.add ( [nmol] )
        nmol.openState.xform = mol.openState.xform

        for at in nmol.atoms :
            #at.drawMode = 3
            if at.element.name in atomColors : at.color = atomColors[at.element.name]
            #at.radius = at.radius * 0.8

        mname = dmap.name + "_%s_%d" % (r.type, r.id.position)

        #aCB, aCG, aOD1 = r.atomsMap['CB'][0], r.atomsMap['CG'][0], r.atomsMap['OD1'][0]
        aCB, aCG, aOD1 = r.atomsMap[xyAts[0]][0], r.atomsMap[xyAts[1]][0], r.atomsMap[xyAts[2]][0]

        dmap, mmap = ExtractDen ( r.atoms, dmap, mname, boundRad=2.0, showMesh=True )
        r.exMol = nmol
        r.exMaps = [dmap, mmap]

        X = aOD1.coord() - aCB.coord(); X.normalize()
        Y = aCG.coord() - aCB.coord(); Y.normalize()
        Z = chimera.cross ( X, Y ); Z.normalize()
        X = chimera.cross ( Y, Z ); Y.normalize()

        xf = chimera.Xform.coordFrame ( X, Y, Z, aCB.coord(), True ).inverse()
        xf.premultiply ( chimera.Xform.translation(atX, atY, 0) )
        
        nmol.openState.xform = xf
        dmap.openState.xform = xf
        if mmap : mmap.openState.xform = xf


    def ExtendRes ( self, r, mol, dmap, atX, atY, thrF=0.75 ) :

        nmol, nres = CopyRess ( [r] )
        nmol.name = mol.name + "_%s_%d_ext" % (r.type, r.id.position)
        chimera.openModels.add ( [nmol] )
        nmol.openState.xform = mol.openState.xform

        for at in nmol.atoms :
            at.drawMode = 3
            if at.element.name in atomColors : at.color = atomColors[at.element.name]
            at.radius = at.radius * 0.8

        mname = dmap.name + "_%s_%d_ext" % (r.type, r.id.position)


        R = nres[0]
        R.O, R.N, R.C, R.CA = R.atomsMap["O"][0], R.atomsMap["N"][0], R.atomsMap["C"][0], R.atomsMap["CA"][0]
        R.CB, R.CG, R.OD1, R.OD2 = R.atomsMap["CB"][0], R.atomsMap["CG"][0], R.atomsMap["OD1"][0], R.atomsMap["OD2"][0]

        bones = []
        bones.append ( Bone(R.CA, R.N, R.CB) )
        bones.append ( Bone(R.CA, R.C, R.CB) )
        bones.append ( Bone(R.C, R.O, R.CA) )
        
        bones.append ( Bone(R.CA, R.CB, R.N) )
        bones.append ( Bone(R.CG, R.CB, R.OD1) )
        bones.append ( Bone(R.CG, R.OD1, R.OD2) )
        bones.append ( Bone(R.CG, R.OD2, R.OD1) )

        for bi, bo in enumerate ( bones ) :
            if GetMod ( "bone_%d.mrc" % bi ) != None : chimera.openModels.close ( "bone_%d.mrc" % bi )
            if GetMod ( "bone_%d.mrc_mesh" % bi ) != None : chimera.openModels.close ( "bone_%d.mrc_mesh" % bi )
            bo.dmap = BoneMap ( bo, dmap, 1.0, "bone_%d.mrc" % bi, show = False, showMesh=True )

        v1 = R.CB.coord() - R.CA.coord(); v1.normalize()
        v2 = R.CB.coord() - R.CG.coord(); v2.normalize()
        ang = numpy.arccos ( v1*v2 ) * 180.0/numpy.pi
        ax = chimera.cross ( v1, v2 ); ax.normalize()

        print "CB-CG: %.2f" % (-ang + 180)

        T = chimera.Xform.translation ( R.CB.coord().toVector() )
        T.multiply ( chimera.Xform.rotation ( ax, -ang + 180 ) )
        T.multiply ( chimera.Xform.translation ( R.CB.coord().toVector()*-1.0 ) )

        for an in ["CG", "OD1", "OD2"] :
            at = R.atomsMap[an][0]
            at.setCoord ( T.apply (at.coord()) )

        #MoldMap2 ( bones, rmaps[0], rmaps[1] )
        
        d1 = diha ( R.N, R.CB, R.CG, R.OD1 )
        d2 = diha ( R.N, R.CB, R.CG, R.OD2 )
        ang = d1 if numpy.abs(d1) < numpy.abs(d2) else d2
        print "CG dihedral - ", d1, d2, " -> ", ang
        ax = R.CG.coord() - R.CB.coord(); ax.normalize()

        T = chimera.Xform.translation ( R.CG.coord().toVector() )
        T.multiply ( chimera.Xform.rotation ( ax, -ang ) )
        T.multiply ( chimera.Xform.translation ( R.CG.coord().toVector()*-1.0 ) )
    
        for an in ["OD1", "OD2"] :
            at = R.atomsMap[an][0]
            at.setCoord ( T.apply (at.coord()) )
            
        dmap, dmesh = MapForAtoms ( R.atoms, dmap, mname, showMesh=True, thrF=thrF )
        MoldMap2 ( bones, dmap, dmesh )
        r.xtMol = nmol
        r.xtMaps = [dmap, dmesh]

    
        X = R.OD1.coord() - R.CB.coord(); X.normalize()
        Y = R.CG.coord() - R.CB.coord(); Y.normalize()
        Z = chimera.cross ( X, Y ); Z.normalize()
        X = chimera.cross ( Y, Z ); Y.normalize()

        xf = chimera.Xform.coordFrame ( X, Y, Z, R.CB.coord(), True ).inverse()
        xf.premultiply ( chimera.Xform.translation(atX, atY, 0) )
        
        nmol.openState.xform = xf
        dmap.openState.xform = xf
        if dmesh : dmesh.openState.xform = xf
    
        
        
    
    def asp ( self ) :
    
        N = 1

        framei = 0
        mpath = "/Users/greg/Desktop/frames"
        import os
        for f in os.listdir ( mpath ) :
            if f.endswith(".png") :
                os.remove( mpath + "/" + f )

        dmap, mol = VisMapMod()
        resolution = 3.0 * dmap.data.step[0]

        print "Map: %s, mol: %s" % (dmap.name, mol.name)
        res = chimera.selection.currentResidues()[0]
        print " - res: %s %d.%s" % (res.type, res.id.position, res.id.chainId)
        z = None

        nname = "%s_%d" % ( res.type, res.id.position )

        #for na in ["ASP","molded.mrc","skinned.mrc"] :
        #    m = GetMod ( na )
        #    if m != None :
        #        chimera.openModels.close ( [m] )


        nmol = GetMod ( nname + ".pdb" )
        if nmol == None :
            nmol, nres = CopyRess ( [res] )
            nmol.name = nname + ".pdb"
            chimera.openModels.add ( [nmol] )
            nmol.openState.xform = mol.openState.xform

            xf = nmol.openState.xform
            #xf.multiply ( chimera.Xform.translation ( 0,0,5 ) )
            nmol.openState.xform = xf
            
            for at in nmol.atoms:
                at.drawMode = 3
                if at.element.name in atomColors :
                    at.color = atomColors[at.element.name]
                    at.radius = at.radius * 0.8
    
        nres = nmol.residues
        R = nres[0]

        R.O = R.atomsMap["O"][0]
        R.N = R.atomsMap["N"][0]
        R.C = R.atomsMap["C"][0]
        R.CA = R.atomsMap["CA"][0]
        R.CB = R.atomsMap["CB"][0]
        R.CG = R.atomsMap["CG"][0]
        R.OD1 = R.atomsMap["OD1"][0]
        R.OD2 = R.atomsMap["OD2"][0]


        bones = []
        bones.append ( Bone(R.CA, R.N, R.CB) )
        bones.append ( Bone(R.CA, R.C, R.CB) )
        bones.append ( Bone(R.C, R.O, R.CA) )
        
        bones.append ( Bone(R.CA, R.CB, R.N) )
        bones.append ( Bone(R.CG, R.CB, R.OD1) )
        bones.append ( Bone(R.CG, R.OD1, R.OD2) )
        bones.append ( Bone(R.CG, R.OD2, R.OD1) )
    
        for bi, bo in enumerate ( bones ) :
            if GetMod ( "bone_%d.mrc" % bi ) != None : chimera.openModels.close ( "bone_%d.mrc" % bi )
            if GetMod ( "bone_%d.mrc_mesh" % bi ) != None : chimera.openModels.close ( "bone_%d.mrc_mesh" % bi )
            bo.dmap = BoneMap ( bo, dmap, 1.0, "bone_%d.mrc" % bi, show = False, showMesh=True )


        v1 = R.CB.coord() - R.CA.coord(); v1.normalize()
        v2 = R.CB.coord() - R.CG.coord(); v2.normalize()
        ang = numpy.arccos ( v1*v2 ) * 180.0/numpy.pi
        print ang
        ax = chimera.cross ( v1, v2 ); ax.normalize()

        dmap.display = False
        mol.display = False

        NB = 2
        #N = 90
        toAng = -ang + 180
        dAng = toAng / float(N)

        print "CB-CG: %.2f/%.2f deg" % (toAng, dAng)

        rmaps = None

        for i in range ( N ) :

            print i, 

            T = chimera.Xform.translation ( R.CB.coord().toVector() )
            #T.multiply ( chimera.Xform.rotation ( ax, -ang + 180 ) )
            T.multiply ( chimera.Xform.rotation ( ax, dAng ) )
            T.multiply ( chimera.Xform.translation ( R.CB.coord().toVector()*-1.0 ) )

            for an in ["CG", "OD1", "OD2"] :
                at = R.atomsMap[an][0]
                at.setCoord ( T.apply (at.coord()) )
   
            #SkinMap ( R.atoms, bones, NB, dmap, 2.0, "skinned.mrc", True)
            #MoldMap ( R.atoms, bones, dmap, "molded.mrc", showMesh=True )
            
            if rmaps == None :
                rmaps = MapForAtoms ( R.atoms, dmap, nname+".mrc", showMesh=True )
            #    for m in rmaps :
            #        if m != None :
            #            m.openState.xform = nmol.openState.xform
            
            MoldMap2 ( bones, rmaps[0], rmaps[1] )
            
    
            if N > 1 :
                chimera.viewer.postRedisplay()
                self.toplevel_widget.update_idletasks ()
                chimera.printer.saveImage ( mpath + "/%06d.png" % framei )
                framei += 1
        
        print ""
    
        if 1 :
    
            d1 = diha ( R.N, R.CB, R.CG, R.OD1 )
            d2 = diha ( R.N, R.CB, R.CG, R.OD2 )
            ang = d1 if numpy.abs(d1) < numpy.abs(d2) else d2
            print "CG dihedral - ", d1, d2, " -> ", ang
            ax = R.CG.coord() - R.CB.coord(); ax.normalize()
    
            toAng = -ang
            dAng = toAng / float( max(N/2,1) )
            print "CG dihedral -- %.2f/%.2f deg" % (toAng, dAng)
    
            for i in range ( max(N/2,1) ) :
            
                print i, 
    
                T = chimera.Xform.translation ( R.CG.coord().toVector() )
                T.multiply ( chimera.Xform.rotation ( ax, dAng ) )
                T.multiply ( chimera.Xform.translation ( R.CG.coord().toVector()*-1.0 ) )
            
                for an in ["OD1", "OD2"] :
                    at = R.atomsMap[an][0]
                    at.setCoord ( T.apply (at.coord()) )
                    
                #print "%d bones" % len(bones)
                #PtsToMapSkinD ( R.atoms, bones, NB, dmap, 2.0, "skinned.mrc", True)
                #MoldMap ( R.atoms, bones, dmap, "molded.mrc", showMesh=True )
                MoldMap2 ( bones, rmaps[0], rmaps[1] )
    
                if N > 1 :
                    chimera.viewer.postRedisplay()
                    self.toplevel_widget.update_idletasks ()
                    chimera.printer.saveImage ( mpath + "/%06d.png" % framei )
                    framei += 1
    
    
    
        
        if N > 1 :
            args = [ "/Users/greg/_mol/Chimera.app/Contents/Resources/bin/ffmpeg", "-r", "30", 
                "-i", mpath + "/%06d.png", "-y", "-qscale", "1", "-b", "9000", "-vcodec", "mpeg4",  # mpeg4 libx264
                "-f", "mov", mpath+"/__ares.mov" ]
            
            print "- running: "
            for a in args : print a,
            print ""
            
            import subprocess
            subprocess.call ( args )
            print "done!\n"
        
    

    
def AddSpherePts ( pts, clr, rad ) : 
    
    from chimera import elements, Coord, Atom, MolResId

    ptsMol = GetMod ( "SD points" )

    res = None
    if ptsMol == None:
        from chimera import Molecule, openModels
        ptsMol = Molecule()
        ptsMol.name = "SD points"
        ptsMol.isRealMolecule = False
        openModels.add ( [ptsMol], noprefs = True )
        res = ptsMol.newResidue('marker', chimera.MolResId('1', 1) )
    else :
        res = ptsMol.residues[0]

    for pt in pts :
        a = ptsMol.newAtom('', elements.H)
        res.addAtom(a)

        a.setCoord ( chimera.Point(*pt) )  # ( chimera.Point(*xyz) )
        a.radius = rad
        a.drawMode = Atom.Sphere
        a.color = chimera.MaterialColor ( *clr )
        a.surfaceCategory = 'markers'



def SpherePts ( ctr, rad, N ) :

    thetas, phis = [], []
    from math import acos, sin, cos, sqrt, pi
    for k in range ( 1, N+1 ) :
        h = -1.0 + ( 2.0*float(k-1)/float(N-1) )
        phis.append ( acos(h) )
        thetas.append ( 0 if k == 1 or k == N else
                        (thetas[k-2] + 3.6/sqrt(N*(1.0-h**2.0))) % (2*pi) )

    pts = [None] * N
    for i, theta, phi in zip ( range(N), thetas, phis ):
        v = chimera.Vector (sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi))
        pt = ctr + v * rad
        pts[i] = pt
    
    return pts
    




def SdevRes ( ress, dmap, allAts, bbAts, scAts, allAtTree = None, show=0, toRAD=4.0, dRAD=0.2 ) :

    pts = []
    for res in ress :
        for at in res.atoms :
            if allAts or (scAts and at.isSC) or (bbAts and at.isBB) :
                p = at.coord()
                pts.append ( [p[0], p[1], p[2]] )
    
    if len(pts) == 0 :
        return None

    RD_ = []             
    d_vals = dmap.interpolated_values ( pts, ress[0].molecule.openState.xform )
    avg = numpy.average ( d_vals )
    if show : print "0\t%f\t%f" % (avg,avg)
    RD_.append ( [0,avg] )

    if show :
        #for r, clr in [ [0.5], [1.0], [1.5] ] :
        T = 0.5
        #for RAD, clr, arad in [ [0.5, (.8,.2,.2,T), 0.05], [1, (.2,.8,.2,T), 0.075], [1.5, (.2,.2,.8,T), 0.1], [2.0, (.8,.2,.8,T), 0.1] ] :
        for RAD, clr, arad in [ [0.5, (.8,.2,.2,T), 0.05], [1, (.2,.8,.2,T), 0.075], [1.5, (.2,.2,.8,T), 0.1] ] :

            pts = []
            for res in ress :
                for at in res.atoms :
                    if allAts or (scAts and at.isSC) or (bbAts and at.isBB) :
                        outPts = SpherePts ( at.coord(), RAD, 10 )
                        for pt in outPts :
                            if allAtTree != None :
                                opointsNear = allAtTree.searchTree ( [pt[0], pt[1], pt[2]], RAD*0.75 )
                                if len(opointsNear) > 0 :
                                    continue
                                else :
                                    pts.append ( [pt[0], pt[1], pt[2]] )
                            else :
                                pts.append ( [pt[0], pt[1], pt[2]] )

            AddSpherePts ( pts, clr, arad )

            d_vals = dmap.interpolated_values ( pts, res.molecule.openState.xform )
            avg = numpy.average ( d_vals )
            RD_.append ( [RAD,avg] )
            #if show : print "%.1f\t%f" % (RAD, avg)

    #dRAD, toRAD, RAD = 0.2, 1.8, 0.1
    RAD = 0.2

    while RAD < toRAD + 0.1 :
        pts = []
        for res in ress :
            for at in res.atoms :
                if allAts or (scAts and at.isSC) or (bbAts and at.isBB) :
                    outPts = SpherePts ( at.coord(), RAD, 10 )
                    for pt in outPts :
                        if allAtTree != None :
                            opointsNear = allAtTree.searchTree ( [pt[0], pt[1], pt[2]], RAD*0.75 )
                            if len(opointsNear) > 0 :
                                continue
                            else :
                                pts.append ( [pt[0], pt[1], pt[2]] )
                        else :
                            pts.append ( [pt[0], pt[1], pt[2]] )

        d_vals = dmap.interpolated_values ( pts, res.molecule.openState.xform )
        avg = numpy.average ( d_vals )
        gv = RD_[0][1] * numpy.exp ( -0.5 * numpy.power(RAD/0.8,2) )
        if show : print "%.1f\t%f\t%f" % (RAD, avg, gv)
        RD_.append ( [RAD,avg] )
        RAD += dRAD

    #minSd = opt0 ( RD_, 0.0001 )
    #if show : print " SD: %.1f" % minSd
    
    minSd = opt ( RD_, 0.0001 )
    if show : print " SD: %.4f" % minSd
    
    
    
    return minSd
    

def err ( XYz, sd ) :

    y0 = XYz[0][1]
    err = 0
    for x,y in XYz[1:] :
        yd = y - y0 * numpy.exp ( -0.5 * numpy.power(x/sd,2) )
        err += yd * yd
    #err /= float(len(XYz))
    return err


def opt0 ( RD_ ) :

    sd = 0.1
    y0 = RD_[0][1]
    minSd, minErr, N = None, 1e99, float ( len(RD_)-1 )
    while sd < 10.0 :
    
        err = 0
        for x,y in RD_[1:] :
            yd = y - y0 * numpy.exp ( -0.5 * numpy.power(x/sd,2) )
            err += yd * yd
        err /= N
        
        #print err
        
        if err < minErr :
            minErr = err
            minSd = sd
        
        sd += 0.1

def opt ( V, maxErr ) :

    dd = 1.0
    sdAt = 0.1
    lastE = err ( V, sdAt )
    #while True :
    for i in range(100000) :
        sdAt += dd
        e = err ( V, sdAt )
        #print "%d %.2f %f %.4f" % (i, sdAt, numpy.log(e), dd)
        if e > lastE :
            dd *= -0.75
            if abs(dd) < maxErr :
                return sdAt
        lastE = e
    return sdAt
    


def CurMolAndChain () :

    segModDialog = modelz_dialog ()
    if segModDialog != None :
        
        if segModDialog.cur_mol == None :
            segModDialog.cur_mol = chimera.Molecule()
            segModDialog.cur_mol.name = "Model"
            #chimera.openModels.add ( [mol], noprefs = True )
            chimera.openModels.add ( [segModDialog.cur_mol] )
            segModDialog.struc.set ( segModDialog.cur_mol.name )

            try :
                segModDialog.cur_mol.openState.xform = chimera.openModels.list()[0].openState.xform
            except :
                pass
        
        chainId = segModDialog.molChain.get()
        if len(chainId) == 0 :
            chainId = "A"
            segModDialog.molChain.set ( chainId )
        
        return segModDialog.cur_mol, chainId
    
    return None, ""


def VisMapMod () :
    
    mol, map = None, None
    
    for m in OML(modelTypes = [chimera.Molecule]) :
        if m.display :
            mol = m

    for m in OML(modelTypes = [VolumeViewer.volume.Volume]) :
        if m.display :
            map = m
    
    return map, mol
    


def ZScoresVis ( ) :
    
    map, mol = VisMapMod()

    if mol != None and map != None :
        ZScores ( mol, map)
    else :
        print "Did not find visible mol and map"



def ZScores ( mol, map ) :

    resolution = 3.0 * map.data.step[0]
    print "Mol: %s, Map: %s -- res %.1f" % (mol.name, map.name, resolution)
    
    SetBBAts ( mol )


    cmap = {}
    for r in mol.residues :
        
        if r.id.chainId in cmap :
            cmap[r.id.chainId].append ( [r.id.position, r] )
        else :
            cmap[r.id.chainId] = [ [r.id.position, r] ]


    #ress = cmap['0']

    allBB, allSC = [], []

    for cid, ress in cmap.iteritems() :
        print " - chain %s" % cid

        ress.sort ()
        ares = [el[1] for el in ress]

        zscores = []
        if 0 :
            sses = SSEs ( ares )
            for el in sses :
                si, ei, ss, elRess = el
                zscore, ccs = zBB ( mol, elRess, resolution, map )
                #print ss, si, "-", ei, zscore
                if zscore != None :
                    zscores.append ( zscore )
                for r in elRess :
                    r.bbZ = zscore

        else :
            bbs = BBsegs ( self.seqRes )
            W = 3
            print " - %d BB segments" % len(bbs)
            for bb in bbs :
                print "  %d res, %d-%d" % (len(bb),bb[0].id.position,bb[-1].id.position)

                for ri, r in enumerate ( bb ) :
                    firstRi = max ( 0, ri-(W-1)/2 )
                    lastRi = min ( len(bb)-1, ri+(W-1)/2 )
                    ress = bb[firstRi:lastRi+1]
                    zscore, ccs = zBB ( self.cur_mol, ress, resolution, self.cur_dmap )
                    if zscore != None :
                        zscores.append ( zscore )

        avgBB = 0
        if len(zscores) > 0 :
            avgBB = numpy.average(zscores)
            allBB.extend ( zscores )
            #print " - BB - min %.2f max %.2f, avg %.2f" % (min(zscores), max(zscores), avgBB )
        #else :
        #    print " - BB - no zscores?"


        avgSC = 0
        zscores = evalSC ( map, mol, ares, None )
        if len(zscores) > 0 :
            avgSC = numpy.average(zscores)
            #print " - SC - min %.2f max %.2f, avg %.2f" % (min(zscores), max(zscores), numpy.average(zscores) )
            allSC.extend ( zscores )
        #else :
        #    print " - SC - no zscores?"
                
        print "Chain %s - %d res - avgBB %.2f, avgSC %.2f" % ( cid, len(ares), avgBB, avgSC )


    print ""

    avgBB = 0
    if len(avgBB) > 0 :
        avgBB = numpy.average(allBB)
        print "BB All - %d scores - min %.2f max %.2f, avg %.2f" % (len(allBB), min(allBB), max(allBB), avgBB )
    else :
        print "BB - no zscores?"

    avgSC = 0
    if len(allSC) > 0 :
        avgSC = numpy.average(allSC)
        print "SC All - %d scores - min %.2f max %.2f, avg %.2f" % (len(allSC), min(allSC), max(allSC), avgSC )
    else :
        print "SC - no zscores?"
    
    print ""

    



def BBsegs ( ress ) :

    bbs = []
    
    firstRi, atRi = 0, 1
    for r in ress[1:] :
        if ress[atRi].id.position > ress[atRi-1].id.position + 1 or r.rtype == "?" :
            bbs.append ( ress[firstRi:atRi] )
            firstRi = atRi
        atRi += 1

    bbs.append ( ress[firstRi:atRi] )
    
    return bbs




def SSEs ( allRess ) :

    if len(allRess) < 1 :
        return []

    sses, ss = [], ""

    res, rStart = allRess[0], allRess[0]
    #print "  - at first res / pos: %d " % res.id.position
    if res.isHelix :
        ss = "H"
    elif res.isSheet or res.isStrand :
        ss = "E"
    else :
        ss = "_"

    ress = [ res ]
    lastRes = rStart
    for res in allRess [1:] :

        if res.id.position > lastRes.id.position + 1 :
            print " - gap at", res.id.position
            sses.append ( [rStart.id.position, lastRes.id.position, ss, ress] )
            ress = []
            rStart = res
            if res.isHelix :
                ss = "H"
            elif res.isSheet or res.isStrand :
                ss = "E"
            else :
                ss = "_"

        if res.isHelix :
            if ss != "H" :
                #print "%s -> H - at %d rid %d | %d->%d, %d res" % (ss, i, res.id.position, rStart.id.position, lastRes.id.position, len(ress))
                sses.append ( [rStart.id.position, lastRes.id.position, ss, ress] )
                ress = []
                rStart = res
                ss = "H"
        elif res.isSheet or res.isStrand :
            if ss != "E" :
                #print "%s -> E - at %d rid %d | %d->%d, %d res" % (ss, i, res.id.position, rStart.id.position, lastRes.id.position, len(ress))
                sses.append ( [rStart.id.position, lastRes.id.position, ss, ress] )
                ress = []
                rStart = res
                ss = "E"
        else :
            if ss == "H" or ss == "E" :
                #print "%s -> _ at %d rid %d | %d->%d, %d res" % (ss, i, res.id.position, rStart.id.position, lastRes.id.position, len(ress))
                sses.append ( [rStart.id.position, lastRes.id.position, ss, ress] )
                ress = []
                rStart = res
                ss = "_"

        ress.append ( res )
        lastRes = res

    #print "Done at rid %d - %s | %d->%d, %d res" % ( res.id.position, ss, rStart.id.position, res.id.position, len(ress))
    sses.append ( [rStart.id.position, res.id.position, ss, ress] )
    return sses




def evalSC ( dmap, mol, ress, whichRes, whichScore="Z", bbCC=1, bbAvgD=1, task=None ) :

    A = []
    resolution = 3.0 * dmap.data.step[0]

    for ri, res in enumerate ( ress ) :

        if whichScore == "Z" :

            if 1 :
                if res.isProt :
                    res.scZ = zRotSideChain ( mol, res, resolution, dmap, show=False )
                    
                    if res.scZ is not None :
                        avgdSC = avgdAts (mol, res.scAtoms, dmap )
                        avgdBB = avgdAts (mol, res.bbAtoms, dmap )
                        res.scZ = res.scZ * avgdSC * avgdBB / (bbAvgD*bbAvgD)

                elif res.isNA :
                    res.scZ = zRotBase ( mol, res, resolution, dmap, show=False )
                else :
                    print "?_%d.%s_%s" % (res.id.position, res.id.chainId, res.type)
                    res.scZ = 0
            
            else :
                res.scZ = zShakeSC ( mol, res, resolution, dmap, show=False )
            
            
            if res.scZ != None :
                A.append ( res.scZ )
        
        else :
            cc, ccm = ccSC ( mol, res, resolution, dmap )
            if whichScore == "CC" :
                res.scZ = cc
            elif whichScore == "CCm" :
                res.scZ = ccm
            elif whichScore == "Q" :
                res.scZ = cc / bbCC
                
        if task :
            task.updateStatus('Calculating side chain score - %d/%d' % (ri,len(ress)))


    #avgA, stdA = numpy.average ( A ), numpy.std ( A )
    #umsg ( "Avg side chain Z-score: %.3f" % ( avgA ) )
    return A



def MoveSC () :

    map, mol = VisMapMod()
    resolution = 3.0 * map.data.step[0]

    print "Map: %s, mol: %s" % (map.name, mol.name)
    res = chimera.selection.currentResidues()[0]
    print " - res: %s %d.%s" % (res.type, res.id.position, res.id.chainId)
    z = None
    
    if 1 :
        if res.isProt :
            z = zRotSideChain ( mol, res, resolution, map, True )
        elif res.isNA :
            z = zRotBase ( mol, res, resolution, map, True )

    else :
        z = zShakeSC ( mol, res, resolution, map, True )

    print z




def zShakeSC ( mol, res, resolution, dmap, show=False ) :

    atoms = res.scAtoms

    if len(atoms) < 1 :
        #print " - no sc atoms" % len(atoms)
        return None

    score0 = 0
    scores, scorest = [], []
    T = 1
    trange = [-T*1.0, 0.0, T*1.0]
    #trange = [-T*2.0, -T, 0.0, T, T*2.0]
    
    fout = None
    if show :
        fout = open ("/Users/greg/Desktop/sc.txt", "w")

    moved = False

    for xx in trange :
        for yy in trange :
            for zz in trange :

                v = chimera.Vector(xx,yy,zz)
                xfT = chimera.Xform.translation ( chimera.Vector(xx,yy,zz) )

                molg = MyMolMapX ( mol, atoms, resolution, dmap.data.step[0], xfT )

                fpoints, fpoint_weights = fit_points_g ( molg )
                map_values = dmap.interpolated_values ( fpoints, mol.openState.xform )
                olap, corr1, corr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )

                if numpy.fabs(xx) < .01 and numpy.fabs(yy) < .01 and numpy.fabs(zz) < .01 :
                    score0 = corr1
                else :
                    scores.append ( corr1 )
                    if fout :

                        #if not moved :
                        nmol, cress = CopyRess ( [res] )
                        for nr in cress :
                            for nat in nr.atoms :
                                try :
                                    nat.setCoord ( xfT.apply ( nat.coord() ) )
                                except :
                                    pass
                        #chimera.openModels.add ( [nmol] )
                        nmol.name = "S_%.0f_%.0f_%.0f" % (xx,yy,zz)
                        moved = True

                        scorest.append ( [corr1, [xx,yy,zz], nmol] )


    if fout :
        scorest.sort ()
        #scorest.reverse ()
        scorest = scorest[0:len(scorest)/2]
        if fout :
            fout.write ( "%.0f,%.0f,%.0f\t%f\n" % (0,0,0, score0) )
            for sc, t, nmol in scorest:
                fout.write ( "%.0f,%.0f,%.0f\t%f\n" % (t[0],t[1],t[2], sc) )
                chimera.openModels.add ( [nmol] )
                SetBBAts ( nmol )
                for at in nmol.atoms :
                    at.display = at.isSC
                        

        fout.close()
        
    if 1 :
        scores.sort ()
        #scores.reverse ()
        scores = scores[0:len(scores)/2]

    #print ""
    avg = numpy.average ( scores )  #numpy.average ( scores[1:] )
    stdev = numpy.std ( scores ) #numpy.std ( scores[1:] )
    if stdev < 1e-8 :
        #print " - nostdev"
        return None
    zscore = (score0 - avg) / stdev #(scores[0] - avg) / stdev
    #print " - scores: avg %.4f, std %.4f, z-score %.4f" % (avg, stdev, zscore )
    #fout.close()
    
    return zscore




def zRotSideChain ( mol, r, resolution, dmap, show=False ) :

    r.CA, r.CB, r.CG = None, None, None
    try :
        r.CA = r.atomsMap["CA"][0]
        r.CB = r.atomsMap["CB"][0]
    except :
        pass

    if "CG" in r.atomsMap :
        r.CG = r.atomsMap["CG"][0]
    elif "CG1" in r.atomsMap :
        r.CG = r.atomsMap["CG1"][0]
    elif "CG2" in r.atomsMap :
        r.CG = r.atomsMap["CG2"][0]
    elif "OG" in r.atomsMap :
        r.CG = r.atomsMap["OG"][0]
    elif "SG" in r.atomsMap :
        r.CG = r.atomsMap["SG"][0]

    if r.CA == None or r.CB == None or r.CG == None :
        #print r.type, " - no ats"
        return None

    resolution = 3.0 * dmap.data.step[0]

    scores = []

    #molg = MyMolMap ( mol, r.atoms, resolution, dmap.data.step[0] )
    #fpoints, fpoint_weights = fit_points_g ( molg )
    #map_values = dmap.interpolated_values ( fpoints, mol.openState.xform )
    #olap_0, corr1_0, corr2_0 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )

    rats = r.scAtoms
    nrats = []
    for at in rats :
        try :
            at.p0 = at.coord()
            nrats.append ( at )
        except :
            pass

    fout = None
    if show :
        fout = open ("/Users/greg/Desktop/sc.txt", "w")

    #for ri, rmol in enumerate ( rmols[0:10] ) :
    for deg in range (0, 360, 36) :

        RotAts ( nrats, r.CA, r.CB, deg )

        if fout :
            nmol, cress = CopyRess ( [r] )
            chimera.openModels.add ( [nmol] )
            nmol.name = "SC %d %.0f" % (r.id.position, deg)
            nr = nmol.residues[0]
            SetBBAts ( nmol )
            for at in nr.atoms :
                if at.isBB :
                    at.display = False
                else :
                    at.display = True

        corr = ResCC ( mol, nrats, resolution, dmap )
        scores.append ( corr )

        for at in nrats :
            at.setCoord ( at.p0 )
        
    if fout :
        for sci, sc in enumerate ( scores ):
            fout.write ( "%d\t%f\n" % (sci*36, sc) )

        fout.close()

    zscore1 = None
    if len(scores) > 3 :
        avg = numpy.average ( scores[1:] )
        stdev = numpy.std ( scores[1:] )
        zscore1 = ( (scores[0] - avg) / stdev ) if stdev > 1e-5 else 0
        #print " -1- avg %.4f, std %.4f, z-score %.4f" % (avg, stdev, zscore1 )
    

    return zscore1




def zRotBase ( mol, r, resolution, dmap, show=False ) :

    resolution = 3.0 * dmap.data.step[0]

    scores = []

    rats = r.scAtoms
    nrats = []
    for at in rats :
        try :
            if at.element.name == "H" :
                continue
            at.p0 = at.coord()
            nrats.append ( at )
        except :
            pass

    fout = None
    if show :
        fout = open ("/Users/greg/Desktop/sc.txt", "w")

    #for ri, rmol in enumerate ( rmols[0:10] ) :
    for deg in range (0, 360, 36) :

        RotAts ( nrats, r.atomsMap["C1'"][0], r.baseAt, deg )

        if fout :
            nmol, cress = CopyRess ( [r] )
            chimera.openModels.add ( [nmol] )
            nmol.name = "SC %d %.0f" % (r.id.position, deg)
            nr = nmol.residues[0]
            SetBBAts ( nmol )
            for at in nr.atoms :
                if at.isBB :
                    at.display = False
                else :
                    at.display = True

        corr = ResCC ( mol, nrats, resolution, dmap )
        scores.append ( corr )

        for at in nrats :
            at.setCoord ( at.p0 )
        
    if fout :
        for sci, sc in enumerate ( scores ):
            fout.write ( "%d\t%f\n" % (sci*36, sc) )

        fout.close()

    zscore1 = None
    if len(scores) > 3 :
        avg = numpy.average ( scores[1:] )
        stdev = numpy.std ( scores[1:] )
        zscore1 = ( (scores[0] - avg) / stdev ) if stdev > 1e-5 else 0
        #print " -1- avg %.4f, std %.4f, z-score %.4f" % (avg, stdev, zscore1 )
    

    return zscore1







def MoveBB () :

    map, mol = VisMapMod()
    resolution = 3.0 * map.data.step[0]

    print "Map: %s, mol: %s" % (map.name, mol.name)
    z, cc = zBB ( mol, chimera.selection.currentResidues(), resolution, map, True )
    print z



def zBB ( mol, ress, resolution, dmap, show=False ) :

    atoms = []
    for r in ress :
        #if 'C' in r.atomsMap : atoms.append ( r.atomsMap['C'][0] )
        #if 'N' in r.atomsMap : atoms.append ( r.atomsMap['N'][0] )
        #if 'CA' in r.atomsMap : atoms.append ( r.atomsMap['CA'][0] )
        #if 'O' in r.atomsMap : atoms.append ( r.atomsMap['O'][0] )
        atoms.extend ( r.bbAtoms )
        #atoms.extend ( r.scAtoms )

    if len(atoms) < 1 :
        #print " - no atoms" % len(atoms)
        return [0,0]

    score0 = 0
    scores, scorest = [], []
    T = 2
    trange = [-T*1.0, 0.0, T*1.0]
    #trange = [-T*2.0, -T, 0.0, T, T*2.0]

    fout = None
    if show :
        fout = open ("/Users/greg/Desktop/sse.txt", "w")

    moved = False

    for xx in trange :
        for yy in trange :
            for zz in trange :

                v = chimera.Vector(xx,yy,zz)
                xfT = chimera.Xform.translation ( chimera.Vector(xx,yy,zz) )

                molg = MyMolMapX ( mol, atoms, resolution, dmap.data.step[0], xfT )

                fpoints, fpoint_weights = fit_points_g ( molg )
                map_values = dmap.interpolated_values ( fpoints, mol.openState.xform )
                olap, corr1, corr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )

                if numpy.fabs(xx) < .01 and numpy.fabs(yy) < .01 and numpy.fabs(zz) < .01 :
                    score0 = corr2
                else :
                    scores.append ( corr2 )
                    if fout :
                        scorest.append ( [corr2, [xx,yy,zz]] )

                        if not moved :
                            nmol, cress = CopyRess ( ress )
                            for nr in cress :
                                for nat in nr.atoms :
                                    try :
                                        nat.setCoord ( xfT.apply ( nat.coord() ) )
                                    except :
                                        pass
                            chimera.openModels.add ( [nmol] )
                            nmol.name = "T_%.0f_%.0f_%.0f" % (xx,yy,zz)
                            moved = True


    if fout :
        scorest.sort ()
        scorest.reverse ()
        scorest = scorest[len(scorest)/2:]
        if fout :
            fout.write ( "%.0f,%.0f,%.0f\t%f\n" % (0,0,0, score0) )
            for sc, t in scorest:
                fout.write ( "%.0f,%.0f,%.0f\t%f\n" % (t[0],t[1],t[2], sc) )

        fout.close()
        
    if 0 :
        scores.sort ()
        scores.reverse ()
        scores = scores[len(scores)/2:]

    #print ""
    avg = numpy.average ( scores )  #numpy.average ( scores[1:] )
    stdev = numpy.std ( scores ) #numpy.std ( scores[1:] )
    if stdev < 1e-8 :
        #print " - nostdev"
        return [0,0]
    zscore = (score0 - avg) / stdev #(scores[0] - avg) / stdev
    #print " - scores: avg %.4f, std %.4f, z-score %.4f" % (avg, stdev, zscore )
    #fout.close()
    
    return [zscore, score0]



def avgdBB ( mol, ress, dmap, show=False ) :
    
    atoms = []
    for r in ress :
        atoms.extend ( r.bbAtoms )

    if len(atoms) < 1 :
        #print " - no atoms" % len(atoms)
        return 0

    from _multiscale import get_atom_coordinates
    apos = get_atom_coordinates(atoms, transformed = False)

    dvals = dmap.interpolated_values ( apos, mol.openState.xform )
    
    #print dvals
    return numpy.average(dvals)


def avgdAts ( mol, atoms, dmap, show=False ) :
    

    if len(atoms) < 1 :
        #print " - no atoms" % len(atoms)
        return 0

    from _multiscale import get_atom_coordinates
    apos = get_atom_coordinates(atoms, transformed = False)

    dvals = dmap.interpolated_values ( apos, mol.openState.xform )
    
    #print dvals
    return numpy.average(dvals)
    



def ccBB ( mol, ress, resolution, dmap, show=False ) :

    atoms = []
    for r in ress :
        #if 'C' in r.atomsMap : atoms.append ( r.atomsMap['C'][0] )
        #if 'N' in r.atomsMap : atoms.append ( r.atomsMap['N'][0] )
        #if 'CA' in r.atomsMap : atoms.append ( r.atomsMap['CA'][0] )
        #if 'O' in r.atomsMap : atoms.append ( r.atomsMap['O'][0] )
        atoms.extend ( r.bbAtoms )
        #atoms.extend ( r.scAtoms )

    if len(atoms) < 1 :
        #print " - no atoms" % len(atoms)
        return [0,0]

    molg = MyMolMapX ( mol, atoms, resolution, dmap.data.step[0], chimera.Xform.identity() )

    fpoints, fpoint_weights = fit_points_g ( molg )
    map_values = dmap.interpolated_values ( fpoints, mol.openState.xform )
    olap, corr1, corr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )

    
    return [corr1, corr2]




def ccSC ( mol, r, resolution, dmap, show=False ) :

    if len(r.scAtoms) < 1 :
        #print " - no atoms" % len(atoms)
        return [0,0]

    molg = MyMolMapX ( mol, r.scAtoms, resolution, dmap.data.step[0], chimera.Xform.identity() )

    fpoints, fpoint_weights = fit_points_g ( molg )
    map_values = dmap.interpolated_values ( fpoints, mol.openState.xform )
    olap, corr1, corr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )

    
    return [corr1, corr2]







def CopyRess ( res ) :

    nmol = chimera.Molecule()
    ress = [None] * len ( res )

    aMap = dict()
    for ri, r in enumerate ( res ) :
        nres = nmol.newResidue (r.type, chimera.MolResId(r.id.chainId, r.id.position))
        ress[ri] = nres
        for at in r.atoms :
            nat = nmol.newAtom (at.name, chimera.Element(at.element.number))
            aMap[at] = nat
            nres.addAtom( nat )
            p = chimera.Point ( at.coord().x, at.coord().y, at.coord().z )
            nat.setCoord ( p )
            nat.coord0 = chimera.Point ( at.coord().x, at.coord().y, at.coord().z )
            #if at.name == "C" or at.name == 'CA' or at.name == 'O' or at.name == "N" :
            #    at.display = False


            
    for bond in res[0].molecule.bonds :
        try :
            nb = nmol.newBond ( aMap[bond.atoms[0]], aMap[bond.atoms[1]] )
            nb.display = nb.Smart
        except :
            pass
        
    for r in ress :
        r.CA, r.CB, r.CG = None, None, None
        try :
            r.CA = r.atomsMap["CA"][0]
            r.CB = r.atomsMap["CB"][0]
            r.CG = r.atomsMap["CG"][0]
        except :
            pass
    
    return nmol, ress



def RotAts (rats, a1, a2, deg) :
    
    # phi: N -> CA
    p1, p2 = a1.coord(), a2.coord()
    v = p2 - p1; v.normalize()

    xf = chimera.Xform.translation ( p1.toVector() )
    xf.multiply ( chimera.Xform.rotation ( v, deg ) )
    xf.multiply ( chimera.Xform.translation ( p1.toVector() * -1.0 ) )

    #for at in res.atoms :
    #    if at.name != 'C' and at.name != 'CA' and at.name != 'N' and at.name != 'CB' and at.name != 'O' :
    for at in rats :
            at.setCoord ( xf.apply (at.coord()) )
            




def molecule_grid_dataX (m0, atoms, resolution, step, pad, xfT, cutoff_range, sigma_factor, transforms = [], csys = None):

    from _multiscale import get_atom_coordinates
    xyz = get_atom_coordinates(atoms, transformed = True)

    # Transform coordinates to local coordinates of the molecule containing
    # the first atom.  This handles multiple unaligned molecules.
    # Or if on_grid is specified transform to grid coordinates.
    #m0 = atoms[0].molecule
    xf = m0.openState.xform
    xf.multiply ( xfT )
    import Matrix as M
    M.transform_points(xyz, M.xform_matrix(xf.inverse()))
    if csys:
        xf.premultiply(csys.xform.inverse())
    tflist = M.coordinate_transform_list(transforms, M.xform_matrix(xf))

    anum = [a.element.number for a in atoms]

    molecules = set([a.molecule for a in atoms])
    if len(molecules) > 1:
        name = 'molmap res %.3g' % (resolution,)
    else:
        name = 'molmap %s res %.3g' % (m0.name, resolution)

    grid = bounding_grid(xyz, step, pad, tflist)
    grid.name = name

    sdev = resolution * sigma_factor
    add_gaussians(grid, xyz, anum, sdev, cutoff_range, tflist)

    #return grid, molecules
    return grid
        

def MyMolMapX ( m0, atoms, resolution, step, xf ) :

    #from MoleculeMap import molecule_grid_data
    from math import sqrt, pi
    from chimera import openModels as om
    from VolumeViewer import volume_from_grid_data

    atoms = tuple(atoms)

    pad = 3*resolution
    cutoff_range = 5 # in standard deviations
    sigma_factor = 1/(pi*sqrt(2)) # standard deviation / resolution
    transforms,csys = [], None
    display_threshold = 0.95
    
    return molecule_grid_dataX (m0, atoms, resolution, step, pad, xf, cutoff_range, sigma_factor, transforms, csys)



def MyMolMap ( m0, atoms, resolution, step ) :

    #from MoleculeMap import molecule_grid_data
    from math import sqrt, pi
    from chimera import openModels as om
    from VolumeViewer import volume_from_grid_data

    atoms = tuple(atoms)

    pad = 3*resolution
    cutoff_range = 5 # in standard deviations
    sigma_factor = 1/(pi*sqrt(2)) # standard deviation / resolution
    transforms,csys = [], None
    display_threshold = 0.95
    
    return molecule_grid_data(m0, atoms, resolution, step, pad, None, cutoff_range, sigma_factor, transforms, csys)




def molecule_grid_data(m0, atoms, resolution, step, pad, on_grid,
                       cutoff_range, sigma_factor,
                       transforms = [], csys = None):



    from _multiscale import get_atom_coordinates
    xyz = get_atom_coordinates(atoms, transformed = True)

    # Transform coordinates to local coordinates of the molecule containing
    # the first atom.  This handles multiple unaligned molecules.
    # Or if on_grid is specified transform to grid coordinates.
    #m0 = atoms[0].molecule
    xf = on_grid.openState.xform if on_grid else m0.openState.xform
    import Matrix as M
    M.transform_points(xyz, M.xform_matrix(xf.inverse()))
    if csys:
        xf.premultiply(csys.xform.inverse())
    tflist = M.coordinate_transform_list(transforms, M.xform_matrix(xf))

    anum = [a.element.number for a in atoms]

    molecules = set([a.molecule for a in atoms])
    if len(molecules) > 1:
        name = 'molmap res %.3g' % (resolution,)
    else:
        name = 'molmap %s res %.3g' % (m0.name, resolution)

    if on_grid:
        from numpy import float32
        grid = on_grid.region_grid(on_grid.region, float32)
    else:
        grid = bounding_grid(xyz, step, pad, tflist)
    grid.name = name

    sdev = resolution * sigma_factor
    add_gaussians(grid, xyz, anum, sdev, cutoff_range, tflist)

    #return grid, molecules
    return grid




def ResCC ( mol, rats, resolution, dmap ) :

    molg = MyMolMap ( mol, rats, resolution, dmap.data.step[0] )
    
    #if 0 :
    #    fmap = VolumeViewer.volume.volume_from_grid_data ( molg )
    #    fmap.name = "res molmap!"
    #    fpoints, fpoint_weights = fit_points(fmap, False)
    #    map_values = dmap.interpolated_values ( fpoints, fmap.openState.xform )
    #    olap, corr1, corr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    #    scores.append ( corr1 )
    #    chimera.openModels.close ( [fmap] )
    #else :

    fpoints, fpoint_weights = fit_points_g ( molg )
    map_values = dmap.interpolated_values ( fpoints, mol.openState.xform )
    olap, corr1, corr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    return corr1



def fit_points_g (fdata, threshold = 1e-5):

    mat = fdata.full_matrix()

    import _volume
    points = _volume.high_indices(mat, threshold)
    fpoints = points.astype(numpy.single)
    fpoint_weights = mat[points[:,2],points[:,1],points[:,0]]

    nz = numpy.nonzero( fpoint_weights )[0]
    if len(nz) < len (fpoint_weights) :
        fpoints = numpy.take( fpoints, nz, axis=0 )
        fpoint_weights = numpy.take(fpoint_weights, nz, axis=0)

    transform_vertices( fpoints, fdata.ijk_to_xyz_transform )

    if 0 : print "FitPoints from %s with threshold %.4f, %d nonzero" % (
        fmap.name, threshold, len(nz) )

    return fpoints, fpoint_weights






# -----------------------------------------------------------------------------
#
def bounding_grid(xyz, step, pad, transforms):

    xyz_min, xyz_max = point_bounds(xyz, transforms)
    origin = [x-pad for x in xyz_min]
    from math import ceil
    shape = [int(ceil((xyz_max[a] - xyz_min[a] + 2*pad) / step)) for a in (2,1,0)]
    from numpy import zeros, float32
    matrix = zeros(shape, float32)
    from VolumeData import Array_Grid_Data
    grid = Array_Grid_Data(matrix, origin, (step,step,step))
    return grid


# -----------------------------------------------------------------------------
#
def add_gaussians(grid, xyz, weights, sdev, cutoff_range, transforms = []):

    from numpy import zeros, float32, empty
    sdevs = zeros((len(xyz),3), float32)
    for a in (0,1,2):
        sdevs[:,a] = sdev / grid.step[a]

    import Matrix as M
    if len(transforms) == 0:
        transforms = [M.identity_matrix()]
    from _gaussian import sum_of_gaussians
    ijk = empty(xyz.shape, float32)
    matrix = grid.matrix()
    for tf in transforms:
        ijk[:] = xyz
        M.transform_points(ijk, M.multiply_matrices(grid.xyz_to_ijk_transform, tf))
        sum_of_gaussians(ijk, weights, sdevs, cutoff_range, matrix)

    from math import pow, pi
    normalization = pow(2*pi,-1.5)*pow(sdev,-3)
    matrix *= normalization



# -----------------------------------------------------------------------------
#
def point_bounds(xyz, transforms = []):

    from _multiscale import bounding_box
    if transforms :
        from numpy import empty, float32
        xyz0 = empty((len(transforms),3), float32)
        xyz1 = empty((len(transforms),3), float32)
        txyz = empty(xyz.shape, float32)
        import Matrix as M
        for i, tf in enumerate(transforms) :
            txyz[:] = xyz
            M.transform_points(txyz, tf)
            xyz0[i,:], xyz1[i,:] = bounding_box(txyz)
        xyz_min, xyz_max = xyz0.min(axis = 0), xyz1.max(axis = 0)
    else:
        xyz_min, xyz_max = bounding_box(xyz)

    return xyz_min, xyz_max
    




# ---------------------------------------------------------------------------------




def SkinMap ( atoms, bones, N, dmap, atomRad, nname, showMesh = False ) :

    from _multiscale import get_atom_coordinates
    points = get_atom_coordinates ( atoms, transformed = True )

    import _contour
    points0 = numpy.copy ( points )
    _contour.affine_transform_vertices ( points0, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )

    
    nn3, nn2, nn1 = dmap.data.size

    npoints = VolumeData.grid_indices ( (int(nn1), int(nn2), int(nn3) ), numpy.single)  # i,j,k indices
    _contour.affine_transform_vertices ( npoints, dmap.data.ijk_to_xyz_transform )
    _contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( dmap.openState.xform ) )
    _contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( atoms[0].molecule.openState.xform.inverse() ) )
    
    for bo in bones :
        bo.MakeFrame ()

    if N == 1 :
        for pi, p in enumerate ( npoints ) :

            cbone, minDist = None, 1e9
            for bo in bones :
                d = bo.DistToPoint ( p )

                if d < minDist :
                    minDist = d
                    cbone = bo

            pt = cbone.SkinPoint ( p )
            npoints[pi] = pt

    else :

        for pi, p in enumerate ( npoints ) :

            dbos = []
            for bo in bones :
                dbos.append ( [bo.DistToPoint ( p ), bo] )
                
            dbos.sort()
            
            totD = 0.0
            sp = numpy.array ( [0,0,0] )
            for i in range ( N ) :
                d, bo = dbos[i]
                sp = sp + numpy.array ( bo.SkinPoint ( p ) ) * d
                totD += d

            npoints[pi] = sp / totD


    _contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( atoms[0].molecule.openState.xform ) )
    _contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )

    dvals = dmap.interpolated_values ( npoints, dmap.openState.xform )
    nmat = dvals.reshape( (nn3,nn2,nn1) )
    ndata = VolumeData.Array_Grid_Data ( nmat, dmap.data.origin, dmap.data.step, dmap.data.cell_angles )
    
    mdata = VolumeData.zone_masked_grid_data ( ndata, points0, atomRad )
    
    MapFromData ( mdata, nname, dmap, False )
    if showMesh :
        MapFromData ( mdata, nname, dmap, True )



def ExtractDen ( atoms, dmap, nname, boundRad = 2.0, showMesh = False) :

    from _multiscale import get_atom_coordinates
    points1 = get_atom_coordinates ( atoms, transformed = False )
    #COM, U, S, V = prAxes ( points )

    bound = 4.0
    li,lj,lk = numpy.min ( points1, axis=0 ) - (bound, bound, bound)
    hi,hj,hk = numpy.max ( points1, axis=0 ) + (bound, bound, bound)

    nstep = (dmap.data.step[0], dmap.data.step[1], dmap.data.step[2] )

    n1 = int ( numpy.ceil ( (hi - li + 1) / nstep[0] ) )
    n2 = int ( numpy.ceil ( (hj - lj + 1) / nstep[1] ) )
    n3 = int ( numpy.ceil ( (hk - lk + 1) / nstep[2] ) )

    O = chimera.Point ( li, lj, lk  )
    #O = atoms[0].molecule.openState.xform.apply ( O )

    #print " - new map origin:", nO

    npoints = VolumeData.grid_indices ( (n1, n2, n3), numpy.single)  # i,j,k indices
    S = dmap.data.step

    _contour.affine_transform_vertices ( npoints, ((S[0], 0.0, 0.0, O[0]), (0.0, S[1], 0.0, O[1]), (0.0, 0.0, S[1], O[2])) )
    #_contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )

    dvals = dmap.interpolated_values ( npoints, atoms[0].molecule.openState.xform )
    nmat = dvals.reshape( (n3,n2,n1) )

    ndata = VolumeData.Array_Grid_Data ( nmat, O, nstep, dmap.data.cell_angles, name = nname )
  
    #_contour.affine_transform_vertices ( points1, Matrix.xform_matrix( atoms[0].molecule.openState.xform ) )
    #_contour.affine_transform_vertices ( points, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
    mdata = VolumeData.zone_masked_grid_data ( ndata, points1, boundRad )

    dmap = MapFromData ( mdata, nname, dmap, False )
    dmap.openState.xform = atoms[0].molecule.openState.xform
    dmesh = None
    
    if showMesh :
        dmesh = MapFromData ( mdata, nname, dmap, True )
        dmesh.openState.xform = atoms[0].molecule.openState.xform
    
    return [dmap, dmesh]






def BoneMap ( bone, dmap, atomRad, nname, show = False, showMesh = False ) :

    #_contour.affine_transform_vertices ( points, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
    #mdata = VolumeData.zone_masked_grid_data ( dmap.data, points, atomRad )
    
    from _multiscale import get_atom_coordinates
    atoms = [bone.a1, bone.a2]
    points = get_atom_coordinates ( atoms, transformed = True )

    import _contour
    points1 = numpy.copy ( points )
    _contour.affine_transform_vertices ( points1, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
    points0 = numpy.copy ( points1 )
    _contour.affine_transform_vertices ( points1, dmap.data.xyz_to_ijk_transform )

    bound = int ( numpy.ceil( atomRad / dmap.data.step[0] ) ) + 1
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

    wmat = numpy.zeros ( (nn3,nn2,nn1), numpy.float32 )
    ndata = VolumeData.Array_Grid_Data ( wmat, nO, nstep, dmap.data.cell_angles )

    npoints = VolumeData.grid_indices ( (nn1, nn2, nn3), numpy.single)  # i,j,k indices
    npointsi = numpy.copy ( npoints )
    _contour.affine_transform_vertices ( npoints, ndata.ijk_to_xyz_transform )
    _contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( dmap.openState.xform ) )
    _contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( atoms[0].molecule.openState.xform.inverse() ) )

    for pi, p in enumerate ( npoints ) :
        
        i,j,k = npointsi[pi]
        d = bone.DistToPoint ( p )
        if d < atomRad :
            wmat[k,j,i] = 1.0
        else :
            wmat[k,j,i] = 1.0 / numpy.power (1+d-atomRad,8)

    _contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( atoms[0].molecule.openState.xform ) )
    _contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )

    dvals = dmap.interpolated_values ( npoints, dmap.openState.xform )
    nmat = dvals.reshape( (nn3,nn2,nn1) )

    bone.ndata = VolumeData.Array_Grid_Data ( nmat*wmat, nO, nstep, dmap.data.cell_angles, name = nname )
    bone.xfmod = dmap
    
    if show :

        from random import random as rand
        clr = ( rand()*.5+.1, rand()*.5+.1, rand()*.5+.1 )

        bone.dmap = MapFromData ( bone.ndata, nname, dmap, showMesh, color = clr )
        bone.dmap.openState.xform = dmap.openState.xform





def MoldMap ( atoms, bones, dmap, nname, showMesh = False ) :

    
    ndata = dmap.data
    nn3, nn2, nn1 = dmap.data.size
    nO = dmap.data.origin
    nmat = numpy.zeros ( (nn3,nn2,nn1), numpy.float32 )
    nstep = (dmap.data.step[0], dmap.data.step[1], dmap.data.step[2] )

    if 1 :
        ndata = DataForAtoms ( atoms, dmap ) 

    npoints = VolumeData.grid_indices ( (nn1, nn2, nn3), numpy.single)  # i,j,k indices
    _contour.affine_transform_vertices ( npoints, ndata.ijk_to_xyz_transform )

    #_contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( dmap.openState.xform ) )
    #_contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( atoms[0].molecule.openState.xform.inverse() ) )

    for bone in bones :

        npointsc = numpy.copy ( npoints )

        _contour.affine_transform_vertices ( npointsc, Matrix.xform_matrix( bone.Xf().inverse() ) )
        _contour.affine_transform_vertices ( npointsc, Matrix.xform_matrix( bone.Xf0() ) )

        #_contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( atoms[0].molecule.openState.xform ) )
        #_contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )

        _contour.affine_transform_vertices ( npointsc, bone.ndata.xyz_to_ijk_transform )

        p2mt = Matrix.xform_matrix ( chimera.Xform.identity() )
        #dvals, outvals = VolumeData.interpolate_volume_data ( npointsc, p2mt, bone.dmap.data.matrix(), method='linear' )
        dvals, outvals = VolumeData.interpolate_volume_data ( npointsc, p2mt, bone.ndata.matrix(), method='linear' )

        bmat = dvals.reshape( (nn3,nn2,nn1) )
        #nmat = nmat + bmat
        nmat = numpy.maximum ( nmat, bmat )

    #nmat = nmat / float ( len(bones) )

    ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles, name = nname )

    MapFromData ( ndata, nname, dmap, False )
    if showMesh :
        MapFromData ( ndata, nname, dmap, True )




def MoldMap2 ( bones, dmap, dmesh ) :

    
    ndata = dmap.data
    nn1, nn2, nn3 = dmap.data.size
    nO = dmap.data.origin
    nmat = numpy.zeros ( (nn3,nn2,nn1), numpy.float32 )
    nstep = (dmap.data.step[0], dmap.data.step[1], dmap.data.step[2] )

    npoints = VolumeData.grid_indices ( (nn1, nn2, nn3), numpy.single)  # i,j,k indices
    _contour.affine_transform_vertices ( npoints, ndata.ijk_to_xyz_transform )

    _contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( dmap.openState.xform ) )
    _contour.affine_transform_vertices ( npoints, Matrix.xform_matrix (bones[0].a1.molecule.openState.xform.inverse()) )

    for bone in bones :

        npointsc = numpy.copy ( npoints )

        _contour.affine_transform_vertices ( npointsc, Matrix.xform_matrix( bone.Xf().inverse() ) )
        _contour.affine_transform_vertices ( npointsc, Matrix.xform_matrix( bone.Xf0() ) )

        _contour.affine_transform_vertices ( npointsc, Matrix.xform_matrix (bone.a1.molecule.openState.xform) )
        _contour.affine_transform_vertices ( npointsc, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )

        _contour.affine_transform_vertices ( npointsc, bone.ndata.xyz_to_ijk_transform )

        p2mt = Matrix.xform_matrix ( chimera.Xform.identity() )
        #dvals, outvals = VolumeData.interpolate_volume_data ( npointsc, p2mt, bone.dmap.data.matrix(), method='linear' )
        dvals, outvals = VolumeData.interpolate_volume_data ( npointsc, p2mt, bone.ndata.matrix(), method='linear' )

        bmat = dvals.reshape( (nn3,nn2,nn1) )
        #nmat = nmat + bmat
        nmat = numpy.maximum ( nmat, bmat )

    #nmat = nmat / float ( len(bones) )

    #ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles, name = nname )
    dmap.data.full_matrix()[:,:,:] = nmat[:,:,:]
    dmap.data.values_changed()
    MapUp ( dmap, False )

    if dmesh != None :
        dmesh.data.full_matrix()[:,:,:] = nmat[:,:,:]
        dmesh.data.values_changed()
        MapUp ( dmesh, True )




def DataForAtoms ( atoms, dmap, nname = "data for atoms" ) :

    from _multiscale import get_atom_coordinates
    points = get_atom_coordinates ( atoms, transformed = True )

    points1 = numpy.copy ( points )
    _contour.affine_transform_vertices ( points1, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
    #points0 = numpy.copy ( points1 )
    _contour.affine_transform_vertices ( points1, dmap.data.xyz_to_ijk_transform )

    bound = 5
    li,lj,lk = numpy.min ( points1, axis=0 ) - (bound, bound, bound)
    hi,hj,hk = numpy.max ( points1, axis=0 ) + (bound, bound, bound)

    n1 = hi - li + 1
    n2 = hj - lj + 1
    n3 = hk - lk + 1

    nstep = (dmap.data.step[0], dmap.data.step[1], dmap.data.step[2] )

    nn1 = int ( round (dmap.data.step[0] * float(n1) / nstep[0]) )
    nn2 = int ( round (dmap.data.step[1] * float(n2) / nstep[1]) )
    nn3 = int ( round (dmap.data.step[2] * float(n3) / nstep[2]) )

    O = dmap.data.origin
    nO = ( O[0] + float(li) * dmap.data.step[0],
           O[1] + float(lj) * dmap.data.step[1],
           O[2] + float(lk) * dmap.data.step[2] )

    nmat = numpy.zeros ( (nn3,nn2,nn1), numpy.float32 )
    ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles, name = nname )
    return ndata



def MapForAtoms ( atoms, dmap, nname, showMesh=False, thrF = 1.0 ) :

    ndata = DataForAtoms ( atoms, dmap, nname ) 

    m1 = MapFromData ( ndata, nname, dmap, False, thrF=thrF )
    m2 = None

    if showMesh :
        m2 = MapFromData ( ndata, nname, dmap, True, thrF=thrF )

    return [m1,m2]


def MapUp (dmap, showMesh = False, color=(.7,.7,.7,1)) :

    ro = VolumeViewer.volume.Rendering_Options()
    ro.smoothing_factor = .3
    ro.smoothing_iterations = 2
    ro.surface_smoothing = False
    ro.square_mesh = True
    ro.line_thickness = 1

    dmap.update_surface ( False, ro )
    for sp in dmap.surfacePieces :
        v, t = sp.geometry
        if len(v) == 8 and len(t) == 12 :
            sp.display = False
        else :
            if showMesh :
                sp.color = (color[0]/2.0, color[1]/2.0, color[2]/2.0, 1.0)
                sp.displayStyle = sp.Mesh
            else :
                sp.color = (color[0], color[1], color[2], 0.1)


def MapFromData ( ndata, nname, dmap, showMesh, thrF=1.0, color=(.7,.7,.7,1) ) :

    if showMesh :
        m = GetMod ( nname + "_mesh" )
        if m != None :
            chimera.openModels.close ( [m] )
    else :
        m = GetMod ( nname )
        if m != None :
            chimera.openModels.close ( [m] )


    nv = VolumeViewer.volume.volume_from_grid_data ( ndata )
    nv.openState.xform = dmap.openState.xform
    nv.name = nname
    if showMesh :
        nv.name = nname + "_mesh"
    nv.region = ( nv.region[0], nv.region[1], [1,1,1] )
    nv.surface_levels[0] = dmap.surface_levels[0] * thrF

    MapUp(nv, showMesh, color)
    return nv



def diha ( a1, a2, a3, a4 ) :
    #n1 = vnorm ( a1.coord(), a2.coord(), a3.coord() )
    #n2 = vnorm ( a2.coord(), a3.coord(), a4.coord() )
    #return numpy.arccos ( n2 * n1 * -1.0 ) * 180.0 / numpy.pi

    # http://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
    b1 = a2.coord() - a1.coord()
    b2 = a3.coord() - a2.coord()
    b3 = a4.coord() - a3.coord()
    
    n1 = chimera.cross ( b1, b2 ); n1.normalize()
    n2 = chimera.cross ( b2, b3 ); n2.normalize()
    m1 = chimera.cross ( n1, b2 ); m1.normalize()
    
    x = n1 * n2
    y = m1 * n2
    
    return -1.0 * numpy.arctan2 ( y, x) * 180.0 / numpy.pi
    
    
def angle ( a1, a2, a3 ) :
    n1 = a1.coord() - a2.coord()
    n2 = a3.coord() - a2.coord()
    return numpy.arccos ( (n2/n1.length) * (n1/n2.length) )  * 180.0 / numpy.pi


class Bone (object) :

    def __init__ (self, a1, a2, a3) :
        BoneInit ( self, a1, a2, a3 )
        
    def CS ( self ) :
        return CS ( a1.coord(), a2.coord(), a3.coord() )
    
    def CS0 ( self ) :
        return CS ( a1.coord0, a2.coord0, a3.coord0 )
        
    def Xf ( self ) :
        X,Y,Z = CS ( self.a1.coord(), self.a2.coord(), self.a3.coord() )
        return chimera.Xform.coordFrame ( X, Y, Z, self.a1.coord(), True )

    def Xf0 ( self ) :
        X,Y,Z = CS ( self.a1.coord0, self.a2.coord0, self.a3.coord0 )
        return chimera.Xform.coordFrame ( X, Y, Z, self.a1.coord0, True )
        
    def MakeFrame ( self ) :
        BoneMakeFrame ( self )
    
    def DistToPoint ( self, pt ) :
        return BoneDistToPoint ( self, pt )
        
    def SkinPoint ( self, pt ) :
        return BoneSkinPoint ( self, pt )
    

def BoneInit (bo, a1, a2, a3) :
    bo.a1, bo.a2, bo.a3 = a1, a2, a3
    bo.X0, bo.Y0, bo.Z0 = CS ( a1.coord0, a2.coord0, a3.coord0 )
    bo.F0 = chimera.Xform.coordFrame ( bo.X0, bo.Y0, bo.Z0, bo.a1.coord0, True )

def BoneMakeFrame ( bo ) :
    bo.X, bo.Y, bo.Z = CS ( bo.a1.coord(), bo.a2.coord(), bo.a3.coord() )
    bo.F = chimera.Xform.coordFrame ( bo.X, bo.Y, bo.Z, bo.a1.coord(), True )
    bo.F = bo.F.inverse()


def CS ( p1, p2, p3 ) :
    X = p2 - p1; X.normalize()
    Y = p3 - p1; Y.normalize()
    Z = chimera.cross ( X, Y ); Z.normalize()
    Y = chimera.cross ( Z, X ); Y.normalize()
    return X,Y,Z


def BoneDistToPoint ( bo, pt ) :
    
    pt = chimera.Point(pt[0], pt[1], pt[2])
    V = bo.a2.coord() - bo.a1.coord()
    v = pt - bo.a1.coord()
    t = V * v
    if t < 0.0 :
        return v.length
    elif t > 1.0 :
        return (pt-bo.a2.coord()).length
    else :
        lp = bo.a1.coord() + (V*t)
        return (pt-lp).length


def BoneSkinPoint ( bo, pt ) :
    
    #bo.X, bo.Y, bo.Z = CS ( bo.a1.coord(), bo.a2.coord(), bo.a3.coord() )
    #x = chimera.Xform.coordFrame ( bo.X, bo.Y, bo.Z, bo.a1.coord(), True )
    #x = x.inverse()
    #y = chimera.Xform.coordFrame ( bo.X0, bo.Y0, bo.Z0, bo.a1.coord0, True )
    
    pt = chimera.Point ( pt[0], pt[1], pt[2] )
    pt = bo.F.apply ( pt )
    pt = bo.F0.apply ( pt )
    return [pt[0], pt[1], pt[2]]






# ---------------------------------------------------

def modelz_dialog ( create=False ) :

    from chimera import dialogs
    d = dialogs.find ( "modelz", create=False )
    return d



def close_dialog () :
    from chimera import dialogs


def setro (ro) :
    from chimera import dialogs
    d = dialogs.find ( "volume viewer", create=False )
    if d :
        d.surface_options_panel.set_gui_from_rendering_options (ro)
        #d.redisplay_needed_cb()


def vold () :
    from chimera import dialogs
    d = dialogs.find ( "volume viewer", create=False )
    d.surface_options_panel.line_thickness.set(2)
    d.redisplay_needed_cb()
    set_gui_from_rendering_options
    


def show_dialog () :

    from chimera import dialogs

    d = dialogs.find ( "modelz", create=False )
    if d :
        print " - found old diag"
        d.toplevel_widget.update_idletasks ()
        d.Close()
        d.toplevel_widget.update_idletasks ()

    dialogs.register (ModelZ_Dialog.name, ModelZ_Dialog, replace = True)

    d = dialogs.find ( "modelz", create=True )
    # Avoid transient dialog resizing when created and mapped for first time.
    d.toplevel_widget.update_idletasks ()
    d.enter()

    return d



def GetMod ( name ) :
    for m in chimera.openModels.list() :
        if m.name == name :
            return m
    return None



def SetBBAts ( mol ) :
    
    #if hasattr ( mol, "bbats" ) :
    #    return
    #mol.bbats = True

    print " - setting bbAts in %s" % mol.name
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
            if nucleic3to1[r.type] == "G" :
                r.baseAt = r.atomsMap["N9"][0]
            elif nucleic3to1[r.type] == "C" :
                r.baseAt = r.atomsMap["N1"][0]
            elif nucleic3to1[r.type] == "A" :
                r.baseAt = r.atomsMap["N9"][0]
            elif nucleic3to1[r.type] == "U" :
                r.baseAt = r.atomsMap["N1"][0]


        r.bbAtoms = []
        r.scAtoms = []

        if r.isProt :
            for a in r.atoms :
                n = a.name
                a.isBB = n=="C" or n=="CA" or n=="O" or n=="N"
                a.isSC = not a.isBB
                if a.isBB :
                    r.bbAtoms.append ( a )
                else :
                    r.scAtoms.append ( a )

        elif r.isNA :
            for a in r.atoms :
                n = a.name

                a.isBB = n=="P" or n=="O1P" or n=="O2P" or n=="O5'" or n=="C5'" or n=="O3'"
                a.isSugar = n=="C1'" or n=="C2'" or n=="C3'" or n=="C4'" or n=="O4'" or n=="O2'"
                
                if nucleic3to1[r.type] == "G" :
                    a.isBase = n=="N9" or n=="C8" or n=="N7" or n=="C5" or n=="C4" or n=="C6" or n=="O6" or n=="N1" or n=="C2" or n=="N2" or n=="N3"
                    
                elif nucleic3to1[r.type] == "C" :
                    a.isBase = n=="N1" or n=="C2" or n=="O2" or n=="N3" or n=="C4" or n=="N4" or n=="C5" or n=="C6"
                    
                elif nucleic3to1[r.type] == "A" :
                    a.isBase = n=="N9" or n=="C8" or n=="N7" or n=="C5" or n=="C4" or n=="N3" or n=="C2" or n=="N1" or n=="C6" or n=="N6"

                elif nucleic3to1[r.type] == "U" :
                    a.isBase = n=="N1" or n=="C2" or n=="O2" or n=="N3" or n=="C4" or n=="O4" or n=="C5" or n=="C6"
            
                #if nucleic3to1[r.type] == "G" :
                #    r.isBase = n=="" or n=="" or n=="" or n=="" or n=="" or n=="" or n=="" or n=="" or n="" or n="" or n=""
                #    r.baseAt = r.atomsMap["N9"][0]
                
                a.isSC = not a.isBB and not a.isSugar
                if a.isBB :
                    r.bbAtoms.append ( a )
                else :
                    r.scAtoms.append ( a )

    

#def GetVisibleMol () :
#    for m in chimera.openModels.list() :
#        if m.display == True and type(m) == chimera.Molecule :
#            return m
#    return None

NA = {
    "A" : { "baseAtoms" : ["","",""] }
}


class NA ( object ):

    type















