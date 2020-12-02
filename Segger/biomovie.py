
# Copyright (c) 2020 Greg Pintilie - gregdp@gmail.com
# LICENCE - please see: https://opensource.org/licenses/MIT


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
import VolumeViewer.volume
from sys import stderr
from time import clock
import VolumeViewer
import json
import ttk

import axes
reload(axes)

OML = chimera.openModels.list


REG_OPACITY = 0.45


def umsg ( txt ) :
    print txt
    status ( txt )

def status ( txt ) :
    txt = txt.rstrip('\n')
    msg.configure(text = txt)
    msg.update_idletasks()


devMenus = False


class BioMovie ( chimera.baseDialog.ModelessDialog ) :

    title = "BioMovie 0.9.2"
    name = "BioMovie"
    buttons = ( "Close" )
    help = 'https://github.com/gregdp/biomovie'


    def fillInUI(self, parent) :

        # these vars shoul be overwritted from the movie script
        self.framesPath = "/Users/greg/dev/mol/chimera/frames/"
        self.ffmpegPath = None
        self.scriptFun = None
        self.movieFormat = ".mov"



        tw = parent.winfo_toplevel()
        self.toplevel_widget = tw
        tw.withdraw()


        row = 0

        menubar = Tkinter.Menu(parent, type = 'menubar', tearoff = False)
        tw.config(menu = menubar)

        #f = Tkinter.Frame(parent)
        #f.grid(column=0, row=row, sticky='ew')

        #l = Tkinter.Label(f, text='  ')
        #l.grid(column=0, row=row, sticky='w')

        parent.columnconfigure(0, weight = 1)


        if 1 :
            row += 1
            ff = Tkinter.Frame(parent)
            ff.grid(column=0, row=row, sticky='w')

            l = Tkinter.Label(ff, text=" ")
            l.grid(column=0, row=0, sticky='w')

            c = Hybrid.Checkbutton(ff, 'Save', False )
            c.button.grid (column=1, row=0, sticky='w')
            self.makeMovie = c.variable

            c = Hybrid.Checkbutton(ff, 'Stop', False )
            c.button.grid (column=2, row=0, sticky='w')
            self.stopMovie = c.variable

            l = Tkinter.Label(ff, text="  Movie Name: ")
            l.grid(column=3, row=0, sticky='w')

            self.movieName = Tkinter.StringVar(parent)
            self.movieName.set ( "_movie_name" )
            e = Tkinter.Entry(ff, width=24, textvariable=self.movieName)
            e.grid(column=4, row=0, sticky='w', padx=5, pady=5)

            b = Tkinter.Button(ff, text="Go", command=self.Go)
            b.grid (column=5, row=0, sticky='w', padx=5)


        if 0 :
            row += 1
            ff = Tkinter.Frame(parent)
            ff.grid(column=0, row=row, sticky='w')

            l = Tkinter.Label(ff, text=" ")
            l.grid(column=0, row=0, sticky='w')

            b = Tkinter.Button(ff, text="Cycle", command=self.Cycle)
            b.grid (column=1, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(ff, text="FromInitPos", command=self.FromInitPos)
            #b.grid (column=2, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Threshold", command=self.CycleThr)
            b.grid (column=3, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Rot-X", command=self.RotateX)
            b.grid (column=4, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Rot-Y", command=self.RotateY)
            b.grid (column=5, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Rot-Z", command=self.RotateZ)
            b.grid (column=6, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Rock", command=self.Rock)
            b.grid (column=8, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(ff, text="Conf", command=self.Conf)
            #b.grid (column=9, row=0, sticky='w', padx=5)




        if 1 :
            row += 1
            ff = Tkinter.Frame(parent)
            ff.grid(column=0, row=row, sticky='news')

            self.id_keyName = {}

            self.tree = ttk.Treeview(ff)

            #self.tree["columns"]=("one","two","three")
            self.tree.column("#0", width=300, minwidth=100, stretch=Tkinter.YES)
            #self.tree.column("one", width=150, minwidth=150, stretch=Tkinter.NO)
            #self.tree.column("two", width=400, minwidth=200)
            #self.tree.column("three", width=80, minwidth=50, stretch=Tkinter.NO)

            self.tree.heading("#0",text="Views",anchor=Tkinter.W)
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

            l = Tkinter.Label(ff, text="View Name: ")
            l.grid(column=0, row=0, sticky='w')

            self.keyName = Tkinter.StringVar(parent)
            self.keyName.set ( "Side" )
            e = Tkinter.Entry(ff, width=15, textvariable=self.keyName)
            e.grid(column=1, row=0, sticky='w', padx=5, pady=5)

            b = Tkinter.Button(ff, text="Save", command=self.SetKey)
            b.grid (column=2, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(ff, text="Apply", command=self.ApplyKey)
            #b.grid (column=4, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Delete", command=self.DeleteKey)
            b.grid (column=4, row=0, sticky='w', padx=5)

            #l = Tkinter.Label(ff, text="   ")
            #l.grid(column=5, row=0, sticky='w')

            b = Tkinter.Button(ff, text="Load Views", command=self.GetKeys)
            b.grid (column=6, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(ff, text="Open", command=self.OpenFiles)
            #b.grid (column=5, row=0, sticky='w', padx=5)

            # keys = views



        if 1 :
            row += 1
            ff = Tkinter.Frame(parent)
            ff.grid(column=0, row=row, sticky='w')

            l = Tkinter.Label(ff, text="Activate: ")
            l.grid(column=0, row=0, sticky='w')

            b = Tkinter.Button(ff, text="Sel", command=self.ActivateSel)
            b.grid (column=1, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="All", command=self.ActivateAll)
            b.grid (column=2, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Inv", command=self.InvertSel)
            b.grid (column=3, row=0, sticky='w', padx=5)

            l = Tkinter.Label(ff, text="Move: ")
            l.grid(column=10, row=0, sticky='w')

            self.pushA = Tkinter.StringVar(parent)
            self.pushA.set ( "10" )
            e = Tkinter.Entry(ff, width=4, textvariable=self.pushA)
            e.grid(column=11, row=0, sticky='w', padx=5, pady=5)

            b = Tkinter.Button(ff, text="+", command=self.Back10)
            b.grid (column=12, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="-", command=self.Forward10)
            b.grid (column=13, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Ctr", command=self.ComSel)
            b.grid (column=20, row=0, sticky='w', padx=5)


        if 0 :
            row += 1
            ff = Tkinter.Frame(parent)
            ff.grid(column=0, row=row, sticky='w')

            #l = Tkinter.Label(ff, text="Activate: ")
            #l.grid(column=3, row=0, sticky='w')

            b = Tkinter.Button(ff, text="Activate Sel", command=self.ActivateSel)
            b.grid (column=8, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Activate All", command=self.ActivateAll)
            b.grid (column=9, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Invert", command=self.InvertSel)
            b.grid (column=10, row=0, sticky='w', padx=5)

            if 0 :
                l = Tkinter.Label(ff, text="Hide: ")
                l.grid(column=3, row=0, sticky='w')

                b = Tkinter.Button(ff, text="H-Sel", command=self.HideSel)
                b.grid (column=10, row=0, sticky='w', padx=5)

                b = Tkinter.Button(ff, text="H-All", command=self.HideAll)
                b.grid (column=11, row=0, sticky='w', padx=5)


            b = Tkinter.Button(ff, text="Center on Sel", command=self.ComSel)
            b.grid (column=12, row=0, sticky='w', padx=5)



        if devMenus :
            row += 1
            ff = Tkinter.Frame(parent)
            ff.grid(column=0, row=row, sticky='w')

            #l = Tkinter.Label(ff, text="Activate: ")
            #l.grid(column=0, row=0, sticky='w')

            b = Tkinter.Button(ff, text="Sph", command=self.MembraneSphere)
            b.grid (column=1, row=0, sticky='w', padx=5)



        #row += 1
        #f = Tkinter.Frame(parent)
        #f.grid(column=0, row=row, sticky='ew')
        #l = Tkinter.Label(f, text='  ')
        #l.grid(column=0, row=row, sticky='w')


        row += 1
        dummyFrame = Tkinter.Frame(parent, relief='groove', borderwidth=1)
        Tkinter.Frame(dummyFrame).pack()
        dummyFrame.grid(row=row,column=0,columnspan=7, pady=7, sticky='we')


        row = row + 1
        global msg
        msg = Tkinter.Label(parent, width = 50, anchor = 'w', justify = 'left', fg="red")
        msg.grid(column=0, row=row, sticky='ew')
        self.msg = msg
        row += 1


        self.GetKeys()




    def select_mg_cb (self, event):
        print "Sel:", self.tree.selection()
        print "Focus:", self.tree.focus()

        to = self.tree.focus()


        kname = self.id_keyName[to]
        print kname

        self.keyName.set (kname)
        self.ApplyKey ()



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



    def Go ( self ) :

        print "Go - "

        if self.scriptFun == None :
            umsg ("No script set - run 'execfile ([path to script])' in IDLE first")

        else :
            self.scriptFun()


    def GetModels (self) :

        self.oms = {}
        self.visOms = []
        self.allOms = []
        for mod in chimera.openModels.list () :
            om = AnimatableModel ()
            om.FromMod ( mod )
            self.oms[om.mod.name] = om
            self.allOms.append ( om )
            if mod.display :
                self.visOms.append ( om )

        print len(self.oms), "amods,", len(self.visOms), "visible"


    def GetMap (self, mname) :

        am = None
        for mod in chimera.openModels.list () :
            if mod.name == mname :
                if am != None :
                    print "WARNING: two models with same name found."
                    return None
                om = AnimatableModel ()
                om.FromMod ( mod )
                am = om.FromMap()
        if am == None :
            print "ERROR: asked for model with name %s, which was not found (i.e. is not an open model); returning a None object... This will likely cause an exception!" % mname
        return am

    def GetMolecule (self, mname) :

        am = None
        for mod in chimera.openModels.list () :
            if mod.name == mname :
                if am != None :
                    print "WARNING: two models with same name %s found; will confuse movie script..." % mname
                    return None
                om = AnimatableModel ()
                om.FromMod ( mod )
                am = om.FromMol ()
        if am == None :
            print "WARNING: asked for model with name %s, which was not found; returning a None object..." % mname
        return am


    def GetMod (self, mname) :

        am = None
        for mod in chimera.openModels.list () :
            if mod.name == mname :
                if am != None :
                    print "WARNING: two models with same name found."
                    return None
                om = AnimatableModel ()
                om.FromMod ( mod )

                if type(mod) == VolumeViewer.volume.Volume :
                    am = om.FromMap()
                elif type(mod) == chimera.Molecule :
                    am = om.FromMol()

        if am == None :
            print "ERROR: asked for model with name %s, which was not found (i.e. is not an open model); returning a None object... This will likely cause an exception!" % mname
        return am




    def MembraneSphere ( self ) :

        print "sphere?"



    def GetVisMods (self) :

        self.visMods = []
        for mod in chimera.openModels.list () :
            if mod.display :
                om = AnimatableModel ()
                om.FromMod ( mod )
                self.visMods.append ( om )

        print len(self.visMods), " visible"



    def CycleThr ( self ) :

        self.ClearFrames ()
        self.GetAllAMods()

        dmap = self.oms["groel_r16.mrc"].FromMap ()
        #mod = self.oms["c1_p.pdb"].FromMol ()

        MOVIE = Movie ( self, "groel_r16_thr" )

        t = 0; d = 90;
        #t += d; d = 90

        #MOVIE.add ( Show (t, [dmap] ) )
        #MOVIE.add ( RotateMove (t, t+d, [dmap, mod], dmap.comv, [0,1,0], 360.0, [0,0,0], itype="linear" ) )

        MOVIE.add ( VaryThr (t, t+d, dmap, 0.398, .184) )
        t+=d


        MOVIE.make ()



    def RotateZ ( self ) :

        self.ClearFrames ()
        self.GetAllAMods()


        #m1 = self.oms['groel_e16_5143.mrc'].FromMap ()
        #MOVIE = Movie ( self, "groel_e16_rotate" )


        mods = []
        dmap = None
        for m in chimera.openModels.list() :
            if m.display == True and type(m) == VolumeViewer.volume.Volume :
                mp = self.oms[m.name].FromMap ()
                mods.append ( mp )
                if dmap == None :
                    dmap = mp
                print " --v ", m.name

        for m in chimera.openModels.list() :
            if m.display == True and type(m) == chimera.Molecule :
                mods.append ( self.oms[m.name].FromMol () )
                print " --m ", m.name

        #m1 = self.oms['groel_e4_6422.mrc'].FromMap ()
        #mm = self.oms['1xck_B.pdb'].FromMol ()

        MOVIE = Movie ( self, self.movieName.get() )


        t = 0; d = 360
        #MOVIE.add ( RotateM (t, t+d, mods, dmap.comv, [0,-1,0], 360.0, itype="linear") ) # chimera.viewer.camera.center
        MOVIE.add ( RotateM (t, t+d, mods, chimera.viewer.camera.center, [0,0,1], 360.0, itype="linear") ) # chimera.viewer.camera.center

        if 0 :
            t += d; d = 90
            MOVIE.add ( RotateM (t, t+d, mods, dmap.comv, [-1,0,0], 90.0, itype="linear") ) # chimera.viewer.camera.center

            t += d; d = 360
            MOVIE.add ( RotateM (t, t+d, mods, dmap.comv, [0,0,1], 360.0, itype="linear") ) # chimera.viewer.camera.center

            t += d; d = 90
            MOVIE.add ( RotateM (t, t+d, mods, dmap.comv, [1,0,0], 90.0, itype="linear") ) # chimera.viewer.camera.center


        MOVIE.make ()




    def RotateY ( self ) :

        self.ClearFrames ()
        self.GetAllAMods()


        #m1 = self.oms['groel_e16_5143.mrc'].FromMap ()
        #MOVIE = Movie ( self, "groel_e16_rotate" )


        mods = []
        dmap = None
        for m in chimera.openModels.list() :
            if m.display == True and type(m) == VolumeViewer.volume.Volume :
                mp = self.oms[m.name].FromMap ()
                mods.append ( mp )
                if dmap == None :
                    dmap = mp
                print " --v ", m.name

        for m in chimera.openModels.list() :
            if m.display == True and type(m) == chimera.Molecule :
                mods.append ( self.oms[m.name].FromMol () )
                print " --m ", m.name

        #m1 = self.oms['groel_e4_6422.mrc'].FromMap ()
        #mm = self.oms['1xck_B.pdb'].FromMol ()

        MOVIE = Movie ( self, self.movieName.get() )

        #dm = self.oms["mt_threed_07symsf_08122019.hdf"]
        #print dm.mod.name

        frameMul = 30

        t = 0; d = 10 * frameMul
        #MOVIE.add ( RotateM (t, t+d, mods, dm.comv, [0,-1,0], 360.0, itype="linear") ) # chimera.viewer.camera.center
        MOVIE.add ( RotateM (t, t+d, mods, chimera.viewer.camera.center, [0,-1,0], 360.0, itype="linear") ) # chimera.viewer.camera.center
        #MOVIE.add ( RotateMove (t, t+d, mods, dm, dm.comp, [0,-1,0], 360.0, [0,0,0], itype="linear" ) )

        if 0 :
            t += d; d = 90
            MOVIE.add ( RotateM (t, t+d, mods, dmap.comv, [-1,0,0], 90.0, itype="linear") ) # chimera.viewer.camera.center

            t += d; d = 360
            MOVIE.add ( RotateM (t, t+d, mods, dmap.comv, [0,0,1], 360.0, itype="linear") ) # chimera.viewer.camera.center

            t += d; d = 90
            MOVIE.add ( RotateM (t, t+d, mods, dmap.comv, [1,0,0], 90.0, itype="linear") ) # chimera.viewer.camera.center


        MOVIE.make ()





    def RotateX ( self ) :

        self.ClearFrames ()
        self.GetAllAMods()


        #m1 = self.oms['groel_e16_5143.mrc'].FromMap ()
        #MOVIE = Movie ( self, "groel_e16_rotate" )


        mods = []
        dmap = None
        for m in chimera.openModels.list() :
            if m.display == True and type(m) == VolumeViewer.volume.Volume :
                print " --v ", m.name
                mp = self.oms[m.name].FromMap ()
                mods.append ( mp )
                if dmap == None :
                    dmap = mp

            elif m.display == True and type(m) == chimera.Molecule :
                print " --m ", m.name
                mods.append ( self.oms[m.name].FromMol () )

            elif m.display == True and type(m) == _surface.SurfaceModel :
                print " --s ", m.name
                mods.append ( self.oms[m.name].FromSurf () )

            else :
                print " - ", m.name, type(m)


        #m1 = self.oms['groel_e4_6422.mrc'].FromMap ()
        #mm = self.oms['1xck_B.pdb'].FromMol ()

        MOVIE = Movie ( self, self.movieName.get() )

        frameMul = 30

        t = 0; d = 10 * frameMul
        #MOVIE.add ( RotateM (t, t+d, mods, dmap.comv, [0,-1,0], 360.0, itype="linear") ) # chimera.viewer.camera.center
        MOVIE.add ( RotateM (t, t+d, mods, chimera.viewer.camera.center, [-1,0,0], 360.0, itype="linear") ) # chimera.viewer.camera.center

        if 0 :
            t += d; d = 90
            MOVIE.add ( RotateM (t, t+d, mods, dmap.comv, [-1,0,0], 90.0, itype="linear") ) # chimera.viewer.camera.center

            t += d; d = 360
            MOVIE.add ( RotateM (t, t+d, mods, dmap.comv, [0,0,1], 360.0, itype="linear") ) # chimera.viewer.camera.center

            t += d; d = 90
            MOVIE.add ( RotateM (t, t+d, mods, dmap.comv, [1,0,0], 90.0, itype="linear") ) # chimera.viewer.camera.center


        MOVIE.make ()




    def RotateWithMods ( self ) :

        self.ClearFrames ()
        self.GetAllAMods()


        #m1 = self.oms['groel_e8_2221.mrc'].FromMap ()
        m1 = self.oms['groel_r8.mrc'].FromMap ()

        mols = []
        mols.append ( self.oms['1xck_A_f10.pdb'].FromMol () )
        mols.append ( self.oms['1xck_A_f6.pdb'].FromMol () )
        mols.append ( self.oms['1xck_A_f2.pdb'].FromMol () )
        mols.append ( self.oms['1xck_A_f5.pdb'].FromMol () )
        mols.append ( self.oms['1xck_A_f7.pdb'].FromMol () )
        mols.append ( self.oms['1xck_A_f8.pdb'].FromMol () )
        mols.append ( self.oms['1xck_A_f12.pdb'].FromMol () )
        mols.append ( self.oms['1xck_A_f14.pdb'].FromMol () )
        mols.append ( self.oms['1xck_A_f4.pdb'].FromMol () )
        mols.append ( self.oms['1xck_A_f9.pdb'].FromMol () )
        mols.append ( self.oms['1xck_A_f3.pdb'].FromMol () )
        mols.append ( self.oms['1xck_A_f1.pdb'].FromMol () )
        mols.append ( self.oms['1xck_A_f11.pdb'].FromMol () )
        mols.append ( self.oms['1xck_A_f13.pdb'].FromMol () )

        MOVIE = Movie ( self, "groel_r8_rotate_mods" )


        t = 0; d = 360
        MOVIE.add ( Hide(t, mols) )
        MOVIE.add ( RotateM (t, t+d, [m1]+mols, m1.comv, [0,1,0], 360.0, itype="linear") )
        #t += d

        for ti in range ( 14 ) :
            MOVIE.add ( Show (t,  [ mols[ti] ]  ) )
            t += 25


        t += 60; d = 90
        MOVIE.add ( RotateM (t, t+d, [m1]+mols, m1.comv, [1,0,0], 90.0, itype="linear") )

        t += d+60; d = 90
        MOVIE.add ( RotateM (t, t+d, [m1]+mols, m1.comv, [1,0,0], -90.0, itype="linear") )
        t += d


        MOVIE.make ()



    def Rock ( self ) :

        self.ClearFrames ()
        self.GetAllAMods()

        mods = []
        dmap = []
        surfs = []

        for m in chimera.openModels.list() :
            if m.display == True and type(m) == VolumeViewer.volume.Volume :
                mods.append ( self.oms[m.name].FromMap () )
                print " --v ", m.name

        for m in chimera.openModels.list() :
            if m.display == True and type(m) == chimera.Molecule :
                mods.append ( self.oms[m.name].FromMol () )
                print " --m ", m.name

        for m in chimera.openModels.list() :
            if m.display == True and type(m) == _surface.SurfaceModel :
                mods.append ( self.oms[m.name].FromSurf () )
                print " --s ", m.name


        #m1 = self.oms['groel_e4_6422.mrc'].FromMap ()
        #mm = self.oms['1xck_B.pdb'].FromMol ()

        MOVIE = Movie ( self, self.movieName.get() )


        frameMul = 30

        t = 0; d = 10*frameMul

        c = chimera.viewer.camera.center
        p = mods[0].mod.openState.xform.inverse().apply ( chimera.Point(c[0],c[1],c[2]) )

        for i in range (1) :
            MOVIE.add ( Rock (t, t+d, mods, mods[0], p, [0,-1,0], 45.0, 1.0, itype="linear") )
            t +=d


        MOVIE.make ()



    def SetKey ( self ) :

        K = self.keyName.get()
        umsg ("Setting: " + self.keyName.get())

        if not hasattr ( self, "akeys" ) :
            self.akeys = {}

        self.akeys[K] = {}

        self.akeys[K]['camCenter'] = chimera.viewer.camera.center
        self.akeys[K]['camExtent'] = chimera.viewer.camera.extent

        #print " center: ", self.akeys[K]['camCenter']
        #print " extent: ", self.akeys[K]['camExtent']

        for mod in chimera.openModels.list() :
            xf = mod.openState.xform
            self.akeys[K][mod.name] = Matrix.xform_matrix ( xf )
            #mod.kxf[K] = Matrix.xform_matrix ( xf )

        self.UpdateModKeys()
        self.WriteKeys()


    def KeysPath ( self ) :

        kpath = None

        # look for open maps first, get their folder
        for m in chimera.openModels.list() :
            if type(m) == VolumeViewer.volume.Volume :
                mdir, mpfile = os.path.split(m.data.path)
                kpath = os.path.join ( mdir, "_views.txt" )
                break

        if kpath != None :
            return kpath

        # otherwise try open models
        for m in chimera.openModels.list() :
            if type(m) == chimera.Molecule :
                if hasattr ( m, 'openedAs' ) :
                    path, molname = os.path.split ( m.openedAs[0] )
                    kpath = os.path.join ( path, "_views.txt" )
                    break

        return kpath



    def WriteKeys ( self ) :

        kpath = self.KeysPath()
        if kpath == None :
            return

        try :
            with open(kpath, 'w') as outfile:
                json.dump(self.akeys, outfile)
        except :
            #umsg ( "did not save _keys.txt file" )
            print " - did not save _views.txt file"
            return
        print ( "_views.txt saved: " + kpath )



    def GetKeys ( self ) :

        kpath = self.KeysPath()
        if kpath == None :
            return

        self.akeys = {}
        try :
            fin = open(kpath, 'r')
            self.akeys = json.load(fin)
        except :
            print "no keys"
            #print data

        self.UpdateModKeys ()

        print ( "_views.txt loaded: " + kpath )



    def UpdateModKeys ( self ) :

        self.tree.delete(*self.tree.get_children())
        self.id_keyName = {}

        if not hasattr ( self, "akeys" ) :
            self.akeys = {}
            return

        knames = self.akeys.keys()
        knames.sort()

        for kname in knames :
            kid = self.tree.insert("", "end", "", text=kname)
            self.id_keyName[kid] = kname

        for mod in chimera.openModels.list() :
            mod.kxf = {}
            for kname in self.akeys.keys() :
                if mod.name in self.akeys[kname] :
                    mod.kxf[kname] = self.akeys[kname][mod.name]
                else :
                    if len(self.akeys[kname]) > 0 :
                        mn = self.akeys[kname].keys()[0]
                        mod.kxf[kname] = self.akeys[kname][mn]



    def DeleteKey ( self ) :

        fout = self.KeysPath()
        if fout == None :
            umsg ( "Open a map/model before setting a key..." )
            return


        to = self.tree.focus()

        if len(to) == 0 :
            umsg ( 'No view selected' )
            return

        K = self.id_keyName[to]
        umsg ( "Deleting: " + K )

        if not hasattr ( self, "akeys" ) :
            self.akeys = {}
            return

        if K in self.akeys :
            del self.akeys[K]

        self.UpdateModKeys ()
        self.WriteKeys()





    def ApplyKey ( self ) :

        K = self.keyName.get()
        print " - applying key: ", self.keyName.get()

        xf0 = None
        if K in self.akeys :

            for mod in chimera.openModels.list() :
                if not hasattr ( mod, 'kxf' ) :
                    print " - ", mod.name, "has no keys?"
                    mod.kxf = {}

                if K in mod.kxf :
                    mod.xf0 = chimera.Xform(mod.openState.xform)
                    mod.xf1 = Matrix.chimera_xform ( mod.kxf[K] )
                    #mod.openState.xform = Matrix.chimera_xform ( mod.kxf[K] )
                    if xf0 == None :
                        xf0 = mod.kxf[K]
                else :
                    mod.xf0 = chimera.Xform(mod.openState.xform)
                    mod.kxf[K] = xf0
                    mod.xf1 = Matrix.chimera_xform ( mod.kxf[K] )


            self.startCenter = chimera.Point ( *chimera.viewer.camera.center )
            self.startExtent = chimera.viewer.camera.extent

            #print " from center: ", self.startCenter
            #print " from extent: ", self.startExtent

            if 'camCenter' in self.akeys[K] :
                self.toCenter = chimera.Point ( *self.akeys[K]['camCenter'] )
                self.toExtent = self.akeys[K]['camExtent']
                #print " to center: ", self.akeys[K]['camCenter']
                #print " to extent: ", self.akeys[K]['camExtent']
            else :
                print " to center: ?"
                print " to extent: ?"
                self.toCenter = self.startCenter
                self.toExtent = self.startExtent

        else :
            umsg ( "view " + K + " not found" )
            return


        self.InterpToKey ()



    def InterpToKey ( self ) :


        for mod in chimera.openModels.list() :

            if " -surface.for.chain- " in mod.name :
                continue

            if not hasattr ( mod, "COM" ) :
                if type(mod) == chimera.Molecule :
                    sel = chimera.selection.OSLSelection ( "#%d" % mod.id )
                    atoms = sel.atoms()
                    from _multiscale import get_atom_coordinates
                    points = get_atom_coordinates ( atoms, transformed = False )
                    mod.COM, mod.U, mod.S, mod.V = prAxes ( points )
                    mod.comp = chimera.Point ( mod.COM[0], mod.COM[1], mod.COM[2] )
                    print " %s (map) -- " % mod.name, mod.comp
                elif type(mod) == VolumeViewer.volume.Volume :
                    pts, weights = map_points ( mod )
                    if len(pts) == 0 :
                        pts, weights = map_points ( mod, False )
                    mod.COM, mod.U, mod.S, mod.V = prAxes ( pts )
                    mod.comp = chimera.Point ( mod.COM[0], mod.COM[1], mod.COM[2] )
                    print " %s (mol) -- " % mod.name, mod.comp
                elif type(mod) == _surface.SurfaceModel :
                    mod.comp, brad = axes.SurfCtrRad (mod)
                    print " %s (mol) -- " % mod.name, mod.comp

            if not hasattr ( mod, 'comp' ) :
                continue

            mod.fromPos = mod.xf0.apply ( mod.comp )
            mod.toPos = mod.xf1.apply ( mod.comp )
            mod.tvec = mod.toPos - mod.fromPos


        from quaternion import Quaternion, slerp

        #from chimera import tasks, CancelOperation
        #task = tasks.Task("Going to '%s' View 1/20" % self.keyName.get(), modal = True)

        #try :
        N = 20
        for i in range ( N ) :
            f = i / float(N-1)
            #f1 = 1.0 - f0
            f1, f2 = 2.0*f*f*f-3.0*f*f+1.0, 3*f*f-2*f*f*f

            ctr = self.startCenter + (self.toCenter - self.startCenter)*f2
            chimera.viewer.camera.center = (ctr[0], ctr[1], ctr[2])
            chimera.viewer.camera.extent = self.startExtent + (self.toExtent - self.startExtent)*f2
            #print chimera.viewer.camera.extent

            for mod in chimera.openModels.list() :

                if not hasattr ( mod, 'fromPos' ) :
                    continue

                if " -surface.for.chain- " in mod.name :
                    continue

                t0 = mod.xf0.getTranslation ()
                q0 = Quaternion ()
                q0.fromXform ( mod.xf0 )
                #q0i = q0.inverse ()

                t1 = mod.xf1.getTranslation ()
                q1 = Quaternion ()
                q1.fromXform ( mod.xf1 )

                Q = slerp ( q0, q1, f2 )
                Q.normalize()

                xf = chimera.Xform.translation ( mod.fromPos.toVector() + mod.tvec * f2 )
                xf.multiply ( Q.Xform () )
                xf.multiply ( chimera.Xform.translation ( mod.comp.toVector() * -1.0 ) )

                mod.openState.xform = xf

                if hasattr (mod, 'surfMods') :
                    for cid, m in mod.surfMods.iteritems() :
                        m.openState.xform = xf

            print ".",
            #task.updateStatus( "Going to '%s' View %d/%d" % (self.keyName.get(),i+1,N) )
            chimera.viewer.postRedisplay()
            self.toplevel_widget.update_idletasks ()

        #except CancelOperation:
        #    print "Going to view canceled"

        #finally:
        #    print "Going to view done"
        #    task.finished()

        print ""



    def OpenFiles ( self ) :

        dm = chimera.openModels.list()[0]
        #print dm.name, dm.data.path

        mdir, mpfile = os.path.split(dm.data.path)
        #print " - ", mdir
        #print " - ", mpfile

        fout = mdir + "/anim.txt"
        print " -> ", fout

        self.akeys = {}
        try :
            fin = open(fout, 'r')
            self.akeys = json.load(fin)
        except :
            print "no keys"
            #print data

        print "Views:"
        print self.akeys.keys()

        k0 = self.akeys['0']
        files = k0.keys()
        files.sort()
        for f in files :
            if f == dm.name :
                print " - skipping", f
            else :
                print f
                chimera.openModels.open ( mdir + "/" + f )



        else :
            umsg ( "key " + K + " not found" )




    def Back10 ( self ) :

        a = float ( self.pushA.get() )

        for mod in chimera.openModels.list() :
            if mod.openState.active == True :
                xf = mod.openState.xform
                xf.premultiply ( chimera.Xform.translation(0,0,-a) )
                mod.openState.xform = xf


    def Forward10 ( self ) :

        a = float ( self.pushA.get() )

        for mod in chimera.openModels.list() :
            if mod.openState.active == True :
                xf = mod.openState.xform
                xf.premultiply ( chimera.Xform.translation(0,0,+a) )
                mod.openState.xform = xf


    def ActivateSel ( self ) :

        print "a-sel"

        amods = {}

        for mod in chimera.openModels.list() :
            #print mod.name, chimera.selection.containedInCurrent ( mod )
            if chimera.selection.containedInCurrent ( mod ) :
                mod.openState.active = True
                #print mod.name
            else :
                mod.openState.active = False

    def InvertSel ( self ) :

        print "inv-sel"

        amods = {}
        sel = []

        for mod in chimera.openModels.list() :
            #print mod.name, chimera.selection.containedInCurrent ( mod )
            if chimera.selection.containedInCurrent ( mod ) :
                mod.openState.active = False
                #print mod.name
            else :
                mod.openState.active = True
                sel.append ( mod )

        chimera.selection.clearCurrent ()
        for m in sel :
            #print " - selecting: %s" % m.name
            chimera.selection.addCurrent ( m )


    def ActivateAll ( self ) :

        print "a-all"

        amods = {}

        for mod in chimera.openModels.list() :
            mod.openState.active = True

    def HideSel ( self ) :

        print "h-sel"

        amods = {}

        for mod in chimera.openModels.list() :
            #print mod.name, chimera.selection.containedInCurrent ( mod )
            if chimera.selection.containedInCurrent ( mod ) :
                mod.display = False


    def HideAll ( self ) :

        print "h-all"

        amods = {}

        for mod in chimera.openModels.list() :
            mod.display = False


    def ComSel ( self ) :

        selMod = None
        for mod in chimera.openModels.list() :
            #print mod.name, chimera.selection.containedInCurrent ( mod )
            if chimera.selection.containedInCurrent ( mod ) :
                selMod = mod

        print "Sel:", selMod.name

        om = AnimatableModel ()
        om.FromMod ( selMod )
        #print om.COM

        if type(selMod) == VolumeViewer.volume.Volume :
            amod = om.FromMap()
            p0 = amod.mod.openState.xform.apply ( amod.comPt )
            print p0
            #print amod.comp_wc

            chimera.openModels.cofr = amod.comPt_wc

        elif type(selMod) == chimera.Molecule :

            selAts = chimera.selection.currentAtoms()
            #com, N =
            #for at in selAts :

            chimera.openModels.cofr = selAts[0].xformCoord()

            #amod = om.FromMol()




    def Reset ( self ) :

        print "Resetting"
        self.fri = 0



    def SaveFrame (self, framei, ncopies = 1 ) :

        chimera.viewer.postRedisplay()
        self.toplevel_widget.update_idletasks ()
        #chimera.printer.saveImage ( self.framesPath + "%06d.png" % framei )

        import os
        if not os.path.isdir (self.framesPath) :
            try :
                print "Making folder for frames:"
                print " - ", self.framesPath
                os.mkdir ( self.framesPath )
            except :
                umsg ("Could not make folder for frames... stopping. Please specify folder in your movie script using the framesPath variable.")
                self.stopMovie.set(True)
                return

        fname = os.path.join ( self.framesPath, "%06d.png" % framei )
        chimera.printer.saveImage ( fname )

        #import shutil
        #for i in range (ncopies-1) :
        #    framei += 1
        #    shutil.copy ( self.framesPath + "%06d.png" % (framei-1), self.framesPath + "%06d.png" % framei )


    def SaveKeyFrame (self, framei, text ) :

        chimera.viewer.postRedisplay()
        self.toplevel_widget.update_idletasks ()
        chimera.printer.saveImage ( self.framesPath + "_k_%06d_%s.png" % (framei,text) )



    def ClearFrames (self) :

        import os

        if not os.path.isdir (self.framesPath) :
            return

        for f in os.listdir ( self.framesPath ) :
            if f.endswith(".png") :
                os.remove( os.path.join(self.framesPath,f) )


    def MakeMovie (self, name = "movie") :

        if self.ffmpegPath == None :
            self.ffmpegPath = FindFFmpeg()
            if self.ffmpegPath == None :
                print " - did not find ffmpeg path"
                return

        from os.path import join, split
        framesPath = join( self.framesPath, "%06d.png" )
        moviePath = join ( split(self.framesPath)[0], name + self.movieFormat )

        print " - frames from:", framesPath
        print " - movie file:", moviePath

        if self.movieFormat == ".mov" :
            args = [ self.ffmpegPath, "-r", "30", "-i", framesPath, "-y",
                "-qscale", "1", "-b", "9000", "-vcodec", "mpeg4",
                "-f", "mov", moviePath ]
        else :
            args = [ self.ffmpegPath, "-r", "30", "-i", framesPath, "-y",
                "-qscale", "1", "-b", "9000", "-vcodec", "mpeg4",
                "-f", "mp4", moviePath ]

        #    "-f", "mp4", self.framesPath + "_"+name+".mp4" ]

        print "- running: "
        for a in args : print a,
        print ""

        import subprocess
        subprocess.call ( args )
        print "done!\n"




def FindFFmpeg () :

    # backtrack through path to chimera script to find executable
    print "finding ffmpeg exec:"
    import sys
    atPath = sys.argv[0]
    chiPath = None

    #print " - start ", atPath

    for i in range (100) :

        # go back along the path, look for the ffmpeg binary
        atPath = os.path.split ( atPath )[0]
        #print " --- at %s" % atPath

        # mac
        from os.path import join as J
        tryPath = J(J(J(J(atPath,'Contents'),'Resources'),'bin'), 'ffmpeg')
        if os.path.isfile(tryPath) :
            print " - found ffmpeg path: %s" % tryPath
            chiPath = tryPath
            break

        # unix
        tryPath = J(J(atPath,'bin'),'ffmpeg')
        if os.path.isfile(tryPath) :
            print " - found ffmpeg path: %s" % tryPath
            chiPath = tryPath
            break

        # Windows
        tryPath = J(J(atPath,'bin'),'ffmpeg.exe')
        if os.path.isfile(tryPath) :
            print " - found ffmpeg path: %s" % tryPath
            chiPath = tryPath
            break


        if len(atPath) == 0 :
            break

    return chiPath




# ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------


# This is the base object for all other classes
# - basically sets the interpolation factor f based on start and end frame

class Frames(object) :

    def __init__ (self, startStep, endStep) :
        self.start = startStep
        self.end = endStep
        #print "   - frames - %d|%d" % (self.start, self.end)
        AddAction (self)

    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        if self.end <= self.start :
            self.f = 1
        else :
            self.f = float (stepAt - self.start) / float (self.end - self.start)

        #print " - (%d|%d|%d) f:%.2f" % (self.start, stepAt, self.end, self.f),



class VaryThr (Frames) :

    def __init__ (self, startStep, endStep, animMod, startThr, endThr) :
        #print " - vary thr - %s - %.3f|%.3f" % (animMod.mod.name, self.startThr, self.endThr)
        super(VaryThr, self).__init__(startStep, endStep)
        self.animMod = animMod
        self.startThr = startThr
        self.endThr = endThr


    def step (self, stepAt) :

        if stepAt < self.start or stepAt > self.end :
            return

        #print " - VT step",
        super(VaryThr, self).step(stepAt)

        thr = self.startThr + self.f * (self.endThr - self.startThr)
        #print " - thr:%.3f" % thr

        if type (self.animMod) == list :
          for om in self.animMod :
            om.SetSurfThr ( thr )
        else :
          self.animMod.SetSurfThr ( thr )



class VaryAlpha (Frames) :

    def __init__ (self, startStep, endStep, animMod, startA, endA) :
        #print " - vary A - %s - %.3f|%.3f" % (animMod.mod.name, self.startA, self.endA)
        super(VaryAlpha, self).__init__(startStep, endStep)
        self.animMod = animMod
        self.startA = startA
        self.endA = endA


    def step (self, stepAt) :

        if stepAt < self.start or stepAt > self.end :
            return

        #print " - A step",
        super(VaryAlpha, self).step(stepAt)

        self.alphaAt = self.startA + self.f * (self.endA - self.startA)
        #print " - a:%.3f" % a

        if type (self.animMod) == list :
            for om in self.animMod :
                om.SetAlpha ( self.alphaAt )
        else :
            self.animMod.SetAlpha ( self.alphaAt )



class SetColor (Frames) :

    def __init__ (self, atStep, animMod, clr) :
        super(SetColor, self).__init__(atStep, atStep)
        self.start = self.end = atStep
        self.animMod = animMod
        self.toColor = clr
        self.triggered = False


    def step (self, stepAt) :

        if stepAt < self.start or self.triggered :
            return

        self.triggered = True

        #print " - set color:", self.animMod.mod.name, self.toColor

        if type (self.animMod) == list :

            for om in self.animMod :

                if type(om.mod) == VolumeViewer.volume.Volume :
                    om.SetSurfColor ( self.toColor[0], self.toColor[1], self.toColor[2], self.toColor[3] )
                    om.colorAt = self.toColor

                elif type(om.mod) == chimera.Molecule :

                    TODO

                    self.animMod.dispMode = self.toDispMode
                    self.animMod.color = None
                    self.animMod.alphaAt = 1

                    self.animMod.colors = {}
                    for res in self.animMod.mod.residues :
                        if res.ribbonColor != None :
                            self.animMod.colors[res.id.chainId] = res.ribbonColor.rgba()
                            self.animMod.alphaAt = res.ribbonColor.rgba()[3]


                    if type (self.animMod) == list :
                        for om in self.animMod :
                            om.UpdateDisp ()

                    else :
                        self.animMod.UpdateDisp ()

        else :

            if type(self.animMod.mod) == VolumeViewer.volume.Volume :
                self.animMod.SetSurfColor ( self.toColor[0], self.toColor[1], self.toColor[2], self.toColor[3] )
                self.animMod.mod.colorAt = self.toColor

            else :
                TODO




class VaryColor (Frames) :

    def __init__ (self, startStep, endStep, animMod, startC, endC) :
        #print " - vary A - %s - %.3f|%.3f" % (animMod.mod.name, self.startA, self.endA)
        super(VaryColor, self).__init__(startStep, endStep)
        self.animMod = animMod
        self.startC = numpy.array ( startC )
        self.endC = numpy.array ( endC )


    def step (self, stepAt) :

        if stepAt < self.start or stepAt > self.end :
            return

        #print " - A step",
        super(VaryColor, self).step(stepAt)

        C = self.startC + self.f * (self.endC - self.startC)
        self.colorAt = [C[0], C[1], C[2], C[3]]

        if type (self.animMod) == list :
            for om in self.animMod :

                if type(om.mod) == VolumeViewer.volume.Volume :
                    om.SetSurfColor ( C[0], C[1], C[2], C[3] )

                elif type(om.mod) == chimera.Molecule :
                    for r in om.mod.residues :
                        r.ribbonColor = chimera.MaterialColor ( C[0], C[1], C[2], C[3] )

                else :
                    TODO

                om.colorAt = self.colorAt

        else :
            if type(self.animMod.mod) == VolumeViewer.volume.Volume :
                self.animMod.SetSurfColor ( C[0], C[1], C[2], C[3] )
                self.animMod.colorAt = self.colorAt

            elif type(self.animMod.mod) == chimera.Molecule :

                #for r in self.animMod.mod.residues :
                #    r.ribbonColor = chimera.MaterialColor ( C[0], C[1], C[2], C[3] )

                if hasattr (amod.mod, 'surfMods') :
                    for cid, mod in amod.mod.surfMods.iteritems() :
                        #mod.openState.xform = xf_to_pos
                        TODO
                        self.alphaAt = color[3]
                        self.SetModSurfColor ( mod, (color[0], color[1], color[2], self.alphaAt) )


            else :
                TODO




class SetColorRes (Frames) :

    def __init__ (self, atStep, selStr, cmap) :
        super(SetColorRes, self).__init__(atStep, atStep)
        self.start = self.end = atStep
        self.selStr = selStr
        self.cmap = cmap
        self.triggered = False


    def step (self, stepAt) :

        if stepAt < self.start or self.triggered :
            return

        self.triggered = True

        #print " - set color:", self.animMod.mod.name, self.toColor

        sel = chimera.selection.OSLSelection ( self.selStr )
        for r in sel.residues() :
            C = self.cmap[r.id.chainId]
            r.ribbonColor = chimera.MaterialColor ( C[0], C[1], C[2], C[3] )






class SetDisp (Frames) :

    def __init__ (self, atStep, animMod, dmode) :
        super(SetDisp, self).__init__(atStep, atStep)
        self.start = self.end = atStep
        self.animMod = animMod
        self.toDispMode = dmode
        self.triggered = False


    def step (self, stepAt) :

        if stepAt < self.start or self.triggered :
            return

        self.triggered = True

        print " - set disp:", self.animMod.mod.name, self.toDispMode


        if type(self.animMod.mod) == VolumeViewer.volume.Volume :

            if type (self.animMod) == list :
                for om in self.animMod :
                    om.SetMapDisplay ( self.toDispMode )
                    om.dispAt = self.toDispMode

            else :
                self.animMod.SetMapDisplay ( self.toDispMode )
                self.animMod.toDisp = self.toDispMode


        elif type(self.animMod.mod) == chimera.Molecule :

            self.animMod.dispMode = self.toDispMode
            self.animMod.color = None
            self.animMod.alphaAt = 1

            self.animMod.colors = {}
            for res in self.animMod.mod.residues :
                if res.ribbonColor != None :
                    self.animMod.colors[res.id.chainId] = res.ribbonColor.rgba()
                    self.animMod.alphaAt = res.ribbonColor.rgba()[3]


            if type (self.animMod) == list :
                for om in self.animMod :
                    om.UpdateDisp ()

            else :
                self.animMod.UpdateDisp ()






class SetAlpha (Frames) :

    def __init__ (self, atStep, animMod, alpha) :
        super(SetAlpha, self).__init__(atStep, atStep)
        self.start = self.end = atStep
        self.animMod = animMod
        self.toAlpha = alpha
        self.triggered = False


    def step (self, stepAt) :

        if stepAt < self.start or self.triggered :
            return

        self.triggered = True


        if type (self.animMod) == list :
            for om in self.animMod :
                om.SetAlpha ( self.toAlpha )
                om.alphaAt = self.toAlpha

        else :
            print " - set A:%.3f %s" % (self.toAlpha, self.animMod.mod.name)
            self.animMod.SetAlpha ( self.toAlpha )
            self.animMod.alphaAt = self.toAlpha




class SetThr (Frames) :

    def __init__ (self, atStep, animMod, thr) :
        super(SetThr, self).__init__(atStep, atStep)
        self.start = self.end = atStep
        self.animMod = animMod
        self.toThr = thr
        self.triggered = False


    def step (self, stepAt) :

        if stepAt < self.start or self.triggered :
            return

        self.triggered = True

        print " - set thr:%.3f" % self.toThr

        if type (self.animMod) == list :
          for om in self.animMod :
            om.SetSurfThr ( self.toThr )
        else :
          self.animMod.SetSurfThr ( self.toThr )





class Hide (Frames) :

    def __init__ (self, atStep, animMod) :
        super(Hide, self).__init__(atStep, atStep)
        self.animMod = animMod
        self.start = self.end = atStep
        self.triggered = False
        #if type (animMod) == list :
        #    print " - hide - LIST - %.3f" % (self.start)
        #else :
        #    print " - hide - %s - %.3f" % (animMod.mod.name, self.start)



    def step (self, stepAt) :

        if stepAt < self.start or self.triggered :
            return

        self.triggered = True

        if type (self.animMod) == list :
            for om in self.animMod :
                om.Hide ()

                if hasattr ( om.mod, 'morphMod' ) :
                    print "closing morph mod...", om.mod.morphMod.name
                    chimera.openModels.close ( [om.mod.morphMod] )
                    del om.mod.morphMod

        else :
            self.animMod.Hide()

            om = self.animMod
            if hasattr ( om.mod, 'morphMod' ) :
                print "closing morph mod...", om.mod.morphMod.name
                chimera.openModels.close ( [om.mod.morphMod] )
                del om.mod.morphMod


class Show (Frames) :

    def __init__ (self, atStep, animMod) :
        super(Show, self).__init__(atStep, atStep)
        self.animMod = animMod
        self.start = self.end = atStep
        self.triggered = False
        #if type (animMod) == list :
        #    print " - show - LIST - %.3f" % (self.start)
        #else :
        #    print " - show - %s - %.3f" % (animMod.mod.name, self.start)


    def step (self, stepAt) :

        if stepAt < self.start or self.triggered :
            return

        self.triggered = True

        if type (self.animMod) == list :
          for om in self.animMod :
            om.Show ()
        else :
          self.animMod.Show()


class Select :

    def __init__ (self, atStep, animMod) :
        self.animMod = animMod
        self.start = self.end = atStep
        self.triggered = False
        if type (animMod) == list :
            print " - select - LIST - %.3f" % (self.start)
        else :
            print " - select - %s - %.3f" % (animMod.mod.name, self.start)


    def step (self, stepAt) :

        if stepAt < self.start or self.triggered :
            return

        self.triggered = True

        chimera.selection.clearCurrent ()

        if type (self.animMod) == list :
            mods = []
            for om in self.animMod :
                mods.append ( om.mod )
            chimera.selection.addCurrent ( mods )
        else :
            chimera.selection.addCurrent ( [self.animMod.mod] )



class ColorContacts :

    def __init__ (self, atStep, animMods) :
        self.animMods = animMods
        self.start = self.end = atStep
        self.triggered = False


    def step (self, stepAt) :

        if stepAt < self.start or self.triggered :
            return

        self.triggered = True

        for i in range ( len(self.animMods) ) :
            for j in range ( len(self.animMods) ) :
                if i != j :
                    ColorMapsByContact ( self.animMods[i].mod, self.animMods[j].mod )





# Interpolate models from one xform to another
# - xforms are interpolated by separating rotation and center movement
# - each model moves about its own center of mass, so they can appear disjoint
# - this is good for 'exploding' views

class ToView (Frames) :

    def __init__ ( self, startStep, endStep, animMods, toKey, atype="cubic" ) :

        #print " - xf interp - %s" % (animMod.mod.name)
        super(ToView, self).__init__(startStep, endStep)

        self.toKey = toKey
        self.amods = animMods
        self.atype = atype


    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        if stepAt == self.start :

            for amod in self.amods :

                xf0 = amod.mod.openState.xform
                xf1 = Matrix.chimera_xform ( amod.mod.kxf[self.toKey] )

                endCOM_LC = chimera.Point ( amod.COM[0], amod.COM[1], amod.COM[2] )
                endCOM_WC = xf1.apply ( endCOM_LC )

                startCOM_LC = chimera.Point ( amod.COM[0], amod.COM[1], amod.COM[2] )
                startCOM_WC = xf0.apply ( startCOM_LC )

                amod.t_vec = endCOM_WC - startCOM_WC
                amod.to_pos = endCOM_WC
                amod.comlc = startCOM_LC
                amod.comwc = startCOM_WC
                amod.xf0 = xf0
                amod.xf1 = xf1


        #super(Frames, self).step(stepAt)
        super(ToView, self).step(stepAt)

        f = self.f
        # linear interpolation
        f1, f2 = (1.0-f), f

        if self.atype == "cubic" :
            # cubic interpolation
            f1, f2 = 2.0*f*f*f-3.0*f*f+1.0, 3*f*f-2*f*f*f

        for amod in self.amods :

            from quaternion import Quaternion, slerp
            t0 = amod.xf0.getTranslation ()
            q0 = Quaternion ()
            q0.fromXform ( amod.xf0 )
            #q0i = q0.inverse ()

            t1 = amod.xf1.getTranslation ()
            q1 = Quaternion ()
            q1.fromXform ( amod.xf1 )

            Q = slerp ( q0, q1, f2 )
            Q.normalize()

            tv = amod.t_vec * f2
            xf_to_pos = chimera.Xform.translation ( amod.comwc.toVector() + tv )
            xf_to_pos.multiply ( Q.Xform () )
            xf_to_pos.multiply ( chimera.Xform.translation ( amod.comlc.toVector() * -1.0 ) )

            amod.mod.openState.xform = xf_to_pos

            if hasattr (amod.mod, 'surfMods') :
                for cid, mod in amod.mod.surfMods.iteritems() :
                    mod.openState.xform = xf_to_pos


# Interpolate models from one xform to another
# - xforms are interpolated by separating rotation and center movement
# - each model moves about the center of mass of the ctrMod
# - this is good for keeping models together in a complex


class XfInterpKs (Frames) :

    # same as above but keeps mods together, rotates around ctrMod's COM

    def __init__ ( self, startStep, endStep, ctrMod, animMods, fromKey, toKey, atype="cubic" ) :

        #print " - xf interp - %s" % (animMod.mod.name)
        super(XfInterpKs, self).__init__(startStep, endStep)

        #self.startMod = start
        #self.endMod = end
        self.animMod = ctrMod
        self.fromKey = fromKey
        self.toKey = toKey
        self.amods = animMods
        self.atype = atype


    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        if stepAt == self.start :

            xf0 = self.animMod.mod.openState.xform
            if self.fromKey != None :
                xf0 = Matrix.chimera_xform ( self.animMod.mod.kxf[self.fromKey] )

            xf1 = Matrix.chimera_xform ( self.animMod.mod.kxf[self.toKey] )

            endCOM_LC = chimera.Point ( self.animMod.COM[0], self.animMod.COM[1], self.animMod.COM[2] )
            endCOM_WC = xf1.apply ( endCOM_LC )

            startCOM_LC = chimera.Point ( self.animMod.COM[0], self.animMod.COM[1], self.animMod.COM[2] )
            startCOM_WC = xf0.apply ( startCOM_LC )

            self.t_vec = endCOM_WC - startCOM_WC
            self.to_pos = endCOM_WC
            self.comlc = startCOM_LC
            self.comwc = startCOM_WC
            self.xf0 = xf0
            self.xf1 = xf1


        #super(Frames, self).step(stepAt)

        xf0 = self.xf0
        xf1 = self.xf1

        from quaternion import Quaternion, slerp
        t0 = xf0.getTranslation ()
        q0 = Quaternion ()
        q0.fromXform ( xf0 )
        #q0i = q0.inverse ()

        t1 = xf1.getTranslation ()
        q1 = Quaternion ()
        q1.fromXform ( xf1 )


        super(XfInterpKs, self).step(stepAt)
        f = self.f
        # linear interpolation
        f1, f2 = (1.0-f), f

        if self.atype == "cubic" :
            # cubic interpolation
            f1, f2 = 2.0*f*f*f-3.0*f*f+1.0, 3*f*f-2*f*f*f

        #pos = t0 * f1  + t1 * f2

        if 0 :
            s = q0.s * f1  + q1.s * f2
            v = q0.v * f1  + q1.v * f2
            Q = Quaternion ( s, v )
        else :
            Q = slerp ( q0, q1, f2 )

        Q.normalize()

        tv = self.t_vec * f2
        xf_to_pos = chimera.Xform.translation ( self.comwc.toVector() + tv )
        xf_to_pos.multiply ( Q.Xform () )
        xf_to_pos.multiply ( chimera.Xform.translation ( self.comlc.toVector() * -1.0 ) )

        #self.animMod.mod.openState.xform = xf_to_pos

        for amod in self.amods :
            amod.mod.openState.xform = xf_to_pos

            if hasattr (amod.mod, 'surfMods') :
                for cid, mod in amod.mod.surfMods.iteritems() :
                    mod.openState.xform = xf_to_pos






class SetView (Frames) :

    def __init__ ( self, startStep, animMods, toKey ) :
        #print " - xf interp - %s" % (animMod.mod.name)
        super(SetView, self).__init__(startStep, startStep)
        self.toKey = toKey
        self.amods = animMods


    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        if stepAt == self.start :

            for amod in self.amods :

                if not hasattr ( amod.mod, 'kxf' ) :
                    print "SetView: %s doesn't have views set" % amod.mod.name
                elif not self.toKey in amod.mod.kxf :
                    print "SetView: %s doesn't have view %s" % (amod.mod.name, self.toKey)
                else :
                    xf1 = Matrix.chimera_xform ( amod.mod.kxf[self.toKey] )
                    amod.mod.openState.xform = xf1




class SetXf (Frames) :

    def __init__ ( self, startStep, amods, xfMod ) :

        print " - xf set - %s" % (xfMod.mod.name)
        super(SetXf, self).__init__(startStep, startStep)

        self.xfMod = xfMod
        self.amods = amods


    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        if stepAt == self.start :
            for amod in self.amods :
                amod.mod.openState.xform = self.xfMod.mod.openState.xform




class ModInterp (Frames) :

    def __init__ ( self, startStep, endStep, mod, mod0, mod1 ) :
        #print " - xf interp - %s" % (animMod.mod.name)
        super(ModInterp, self).__init__(startStep, endStep)
        self.mod = mod
        self.mod0 = mod0
        self.mod1 = mod1




    def getResMap ( self, mol ) :
        rmap = {}
        for res in mol.residues:
            if res.id.chainId in rmap :
                rmap[res.id.chainId][res.id.position] = res
            else :
                rmap[res.id.chainId] = {}
                rmap[res.id.chainId][res.id.position] = res
        return rmap


    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        if stepAt == self.start :

            #print " - mod interp start: "
            #print "    - %s, %d atoms" % (self.mod.mod.name, len(self.mod.mod.atoms))
            #print "    - %s, %d atoms" % (self.mod0.mod.name, len(self.mod0.mod.atoms))
            #print "    - %s, %d atoms" % (self.mod1.mod.name , len(self.mod1.mod.atoms))

            self.R = self.getResMap ( self.mod.mod )
            self.R0 = self.getResMap ( self.mod0.mod )
            self.R1 = self.getResMap ( self.mod1.mod )

            for cid, rm in self.R.iteritems() :
                #print cid,
                rm0 = self.R0[cid]
                for ri, res in rm.iteritems() :
                    res0 = rm0[ri]
                    for at in res.atoms :
                        try :
                            at0 = res0.atomsMap[at.name][0]
                        except :
                            print " - could not find res %d %s, at %s -- res %d %s" % (res.id.position, res.type, at.name, res0.id.position, res0.type)
                            blah
                        at.setCoord( at0.coord() )
            #print ""


        super(ModInterp, self).step(stepAt)
        f = self.f
        # linear interpolation
        f1, f2 = (1.0-f), f

        # cubic interpolation
        f1, f2 = 2.0*f*f*f-3.0*f*f+1.0, 3*f*f-2*f*f*f


        #for at, at0, at1 in zip ( self.mod.mod.atoms, self.mod0.mod.atoms, self.mod1.mod.atoms ) :
        #    v = at0.coord().toVector() * f1 + at1.coord().toVector() * f2
            #at.setCoord ( chimera.Point( v[0], v[1], v[2] ) )
            #at.setCoord ( chimera.Point( *v ) )

        for cid, rm in self.R.iteritems() :
            #print "cid"
            rm0 = self.R0[cid]
            rm1 = self.R1[cid]
            for ri, res in rm.iteritems() :
                res0 = rm0[ri]
                res1 = rm1[ri]

                if hasattr(self.mod0, 'colors') :
                    C = self.mod0.colors[cid]
                    res.ribbonColor = chimera.MaterialColor ( C[0], C[1], C[2], 1.0 )

                #if res0.type != res1.type :
                #    print " - res %s[%d.%s.%s] ~~ %s[%d.%s.%s]" % (self.mod0.mod.name, res0.id.position, res0.type, res0.id.chainId, self.mod1.mod.name, res1.id.position, res1.type, res1.id.chainId)
                #    haha
                for at in res.atoms :
                    try :
                        at0 = res0.atomsMap[at.name][0]
                    except :
                        print " - did not find atom %s,%d(%s).%s in %s, res %d(%s)" % (at.name, res.id.position, res.type, res.id.chainId, res0.molecule.name, res0.id.position, res0.type )
                        return

                    try :
                        at1 = res1.atomsMap[at.name][0]
                    except :
                        print " - did not find atom %s,%d(%s).%s in %s, res %d(%s)" % (at.name, res.id.position, res.type, res.id.chainId, res1.molecule.name, res1.id.position, res1.type )
                        return

                    v = at0.coord().toVector() * f1 + at1.coord().toVector() * f2
                    at.setCoord ( chimera.Point( *v ) )

                #for atn in ['C', 'N', 'CA'] :
                #    at0 = res0.atomsMap[atn][0]
                #    at1 = res1.atomsMap[atn][0]
                #    #at = res.atomsMap[atn][0]



        self.mod.UpdateDisp ()



class UpdateShape (Frames) :

    def __init__ ( self, startStep, endStep, shapeMod ) :

        #print " - xf interp - %s" % (animMod.mod.name)
        super(ModInterp, self).__init__(startStep, endStep)

        #self.startMod = start
        #self.endMod = end
        self.mod = shapeMod

    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        if stepAt == self.start :

            print " - update shape start"
            print self.mod.shape


        super(ModInterp, self).step(stepAt)
        f = self.f
        # linear interpolation
        f1, f2 = (1.0-f), f

        # cubic interpolation
        f1, f2 = 2.0*f*f*f-3.0*f*f+1.0, 3*f*f-2*f*f*f





class HideRess (Frames) :

    def __init__ ( self, startStep, endStep, amod, ress=None ) :

        #print " - xf interp - %s" % (animMod.mod.name)
        super(HideRess, self).__init__(startStep, endStep)

        #self.startMod = start
        #self.endMod = end
        self.amod = amod
        if ress == None :
            self.ress = range ( 0, len(self.amod.mod.residues) )
        else :
            self.ress = ress


    def getResMap ( self, mol ) :
        rmap = {}
        for res in mol.residues:
            if res.id.chainId in rmap :
                rmap[res.id.chainId][res.id.position] = res
            else :
                rmap[res.id.chainId] = {}
                rmap[res.id.chainId][res.id.position] = res
        return rmap


    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        if stepAt == self.start :

            #print " - setting ribbon off for %d res" % len(self.amod.mod.residues)
            for res in self.amod.mod.residues :
                #nres.ribbonDisplay, nres.ribbonDrawMode = False, 2
                res.ribbonDisplay = False

                for at in res.atoms :
                    at.display = False

            self.lastResShown = 0


        super(HideRess, self).step(stepAt)
        f = self.f
        # linear interpolation
        f1, f2 = (1.0-f), f

        if 0 :
            nRes = float ( len( self.ress ) )
            toResI = int( numpy.round(nRes * f) )

            for ri in self.ress [self.lastResShown : toResI] :
                self.amod.mod.residues[ri].ribbonDisplay = True
                #res.ribbonDisplay = True

            self.lastResShown = toResI
            #self.mod.UpdateDisp ()



class ShowRess (Frames) :

    def __init__ ( self, startStep, endStep, amod, ress=None ) :

        #print " - xf interp - %s" % (animMod.mod.name)
        super(ShowRess, self).__init__(startStep, endStep)

        #self.startMod = start
        #self.endMod = end
        self.amod = amod
        print " ->> ", amod.mod.name
        if ress == None :
            self.ress = range ( 0, len(self.amod.mod.residues) )
        else :
            self.ress = ress

        self.lastResShown = 0


    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        super(ShowRess, self).step(stepAt)
        f = self.f
        # linear interpolation
        f1, f2 = (1.0-f), f

        # cubic interpolation
        #f1, f2 = 2.0*f*f*f-3.0*f*f+1.0, 3*f*f-2*f*f*f

        nRes = float ( len( self.ress ) )
        toResI = int( numpy.round(nRes * f) )

        for ri in self.ress [self.lastResShown : toResI] :
            res = self.amod.mod.residues[ri]
            res.ribbonDisplay = True

        self.lastResShown = toResI
        #self.mod.UpdateDisp ()



class ShowRessAts (Frames) :

    def __init__ ( self, startStep, endStep, amod, ress=None ) :
        super(ShowRessAts, self).__init__(startStep, endStep)

        self.amod = amod
        if ress == None :
            self.ress = range ( 0, len(self.amod.mod.residues) )
        else :
            self.ress = ress

        self.lastResShown = 0


    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        if 0 and stepAt == self.start :
            for res in self.amod.mod.residues :
                res.ribbonDisplay = False
                for at in res.atoms :
                    at.display = False

        super(ShowRessAts, self).step(stepAt)
        f = self.f
        # linear interpolation
        f1, f2 = (1.0-f), f

        # cubic interpolation
        #f1, f2 = 2.0*f*f*f-3.0*f*f+1.0, 3*f*f-2*f*f*f

        nRes = float ( len( self.ress ) )
        toResI = int( numpy.round(nRes * f) )

        for ri in self.ress [self.lastResShown : toResI] :
            res = self.amod.mod.residues[ri]
            res.ribbonDisplay = True
            for at in res.atoms :
                at.display = True

        self.lastResShown = toResI



class ShowAts (Frames) :

    def __init__ ( self, startStep, amod, ress ) :

        #print " - xf interp - %s" % (animMod.mod.name)
        super(ShowAts, self).__init__(startStep, startStep)

        #self.startMod = start
        #self.endMod = end
        self.amod = amod
        self.ress = ress



    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        if stepAt == self.start :

            ac = { 'O' : chimera.MaterialColor( .9, .2, .2, 1.0 ),
                    'C' : chimera.MaterialColor( .7, .7, .7, 1.0 ),
                    'N' : chimera.MaterialColor( .2, .2, .9, 1.0 ),
                    'H' : chimera.MaterialColor( 1, 1, 1, 1.0 ),
                    ' ' : chimera.MaterialColor( .2, .2, .2, 1.0 ),
                     }

            rmap = {}
            for res in self.amod.mod.residues :
                rmap[res.id.position] = res

            for ri in self.ress :
                res = rmap[ri]
                for at in res.atoms :
                    at.drawMode = at.EndCap
                    at.display = True
                    try :
                        at.color = ac[at.name[0]]
                    except :
                        at.color = ac[" "]



class ShowSideChains (Frames) :

    def __init__ ( self, startStep, selStr ) :

        super(ShowSideChains, self).__init__(startStep, startStep)

        self.selStr = selStr




    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        if stepAt == self.start :

            ac = { 'O' : chimera.MaterialColor( .9, .2, .2, 1.0 ),
                    'C' : chimera.MaterialColor( .7, .7, .7, 1.0 ),
                    'N' : chimera.MaterialColor( .2, .2, .9, 1.0 ),
                    'H' : chimera.MaterialColor( 1, 1, 1, 1.0 ),
                    ' ' : chimera.MaterialColor( .2, .2, .2, 1.0 ),
                     }

            sel = chimera.selection.OSLSelection ( self.selStr )
            for r in sel.residues() :
                for at in r.atoms :
                    at.drawMode = at.EndCap
                    at.display = True
                    try :
                        at.color = ac[at.name[0]]
                    except :
                        at.color = ac[" "]


class HideSideChains (Frames) :

    def __init__ ( self, startStep, selStr ) :

        #print " - xf interp - %s" % (animMod.mod.name)
        super(HideSideChains, self).__init__(startStep, startStep)

        self.selStr = selStr




    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        if stepAt == self.start :

            ac = { 'O' : chimera.MaterialColor( .9, .2, .2, 1.0 ),
                    'C' : chimera.MaterialColor( .7, .7, .7, 1.0 ),
                    'N' : chimera.MaterialColor( .2, .2, .9, 1.0 ),
                    'H' : chimera.MaterialColor( 1, 1, 1, 1.0 ),
                    ' ' : chimera.MaterialColor( .2, .2, .2, 1.0 ),
                     }

            sel = chimera.selection.OSLSelection ( self.selStr )
            for r in sel.residues() :
                for at in r.atoms :
                    at.display = False





class HideAts (Frames) :

    def __init__ ( self, startStep, amod, ress ) :

        #print " - xf interp - %s" % (animMod.mod.name)
        super(HideAts, self).__init__(startStep, startStep)

        #self.startMod = start
        #self.endMod = end
        self.amod = amod
        self.ress = ress



    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        if stepAt == self.start :

            rmap = {}
            for res in self.amod.mod.residues :
                rmap[res.id.position] = res

            for ri in self.ress :
                res = rmap[ri]
                for at in res.atoms :
                    at.display = False






# Rotate models around an axis


class Rotate (Frames) :

    def __init__ ( self, startStep, endStep, animMods, refMod, ctrPt, axis, totDeg, itype="cubic" ) :
        # startStep : starting time step
        # endStep: ending time step
        # animMods: which models to animate
        # refMod: reference model; the transform of this model at the starting
        #         time step will be used to determine center of rotation
        # ctrPt: center of rotation; refMod transform will be applied to this
        # axis: axis of rotation
        # totDeg: how many degrees to rotate in total
        # itype: interpolation type - either cubic or linear

        super(Rotate, self).__init__(startStep, endStep)

        self.animMods = animMods
        self.totDeg = totDeg
        self.axis = chimera.Vector ( axis[0], axis[1], axis[2] )
        self.ctrPt = chimera.Point ( ctrPt[0], ctrPt[1], ctrPt[2] )
        self.itype = itype
        self.refMod = refMod


    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        from quaternion import Quaternion

        if stepAt == self.start :
            self.ctrPtWC = self.refMod.mod.openState.xform.apply ( self.ctrPt )
            for amod in self.animMods :
                amod.xf0 = amod.mod.openState.xform

        #super(Frames, self).step(stepAt)
        super(Rotate, self).step(stepAt)
        f = self.f
        f1, f2 = 2.0*f*f*f-3.0*f*f+1.0, 3*f*f-2*f*f*f  # cubic interpolation
        if not "cubic" == self.itype :
            f1, f2 = (1.0-f), f                        # linear interpolation

        deg = self.totDeg * f2
        #print " - at deg ", deg

        for amod in self.animMods :

            xf = chimera.Xform ( amod.xf0 )
            xf.premultiply ( chimera.Xform.translation(self.ctrPtWC.toVector() * -1.0)  )
            xf.premultiply ( chimera.Xform.rotation ( self.axis, deg )  )
            xf.premultiply ( chimera.Xform.translation(self.ctrPtWC.toVector())  )

            amod.mod.openState.xform = xf

            if hasattr (amod.mod, 'surfMods') :
                for cid, mod in amod.mod.surfMods.iteritems() :
                    mod.openState.xform = xf






class Rock (Frames) :

    def __init__ ( self, startStep, endStep, animMods, comMod, COM, axis, totDeg, numCycles, itype="cubic" ) :
        print " - rock M"
        print "   axis: ", axis
        print "   totDeg: ", totDeg
        print "   numCyc: ", numCycles
        super(Rock, self).__init__(startStep, endStep)

        self.animMods = animMods
        self.totDeg = float ( totDeg )
        self.N = float ( numCycles )
        self.axis = chimera.Vector ( axis[0], axis[1], axis[2] )
        self.comMod = comMod
        self.comlc = chimera.Point ( COM[0], COM[1], COM[2] )
        self.itype = itype


    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        from quaternion import Quaternion

        if stepAt == self.start :
            #print " - rotate M - first step",
            self.comwc = self.comMod.mod.openState.xform.apply ( self.comlc )
            for amod in self.animMods :
                amod.xf0 = amod.mod.openState.xform
                amod.q0 = Quaternion ()
                amod.q0.fromXform ( amod.xf0 )
                amod.comlc = chimera.Point ( amod.COM[0], amod.COM[1], amod.COM[2] )
                amod.comwcv = amod.xf0.apply ( amod.comlc ).toVector() - self.comwc.toVector()

        #super(Frames, self).step(stepAt)

        super(Rock, self).step(stepAt)
        f = self.f
        f1, f2 = 2.0*f*f*f-3.0*f*f+1.0, 3*f*f-2*f*f*f         # cubic interpolation
        if not "cubic" == self.itype :
            f1, f2 = (1.0-f), f                                   # linear interpolation

        #deg = self.totDeg * f2
        deg = self.totDeg * numpy.sin ( f2 * 2.0 * numpy.pi * self.N )
        #print " - at deg ", deg

        for amod in self.animMods :

            xf_to_0 = chimera.Xform.translation ( amod.comlc.toVector() * -1.0 )
            xf_to_P = chimera.Xform.translation ( amod.comwcv )
            xf_rot0 = amod.q0.Xform ()
            xf_rotR = chimera.Xform.rotation ( self.axis, deg )
            xf_to_pos = chimera.Xform.translation ( self.comwc.toVector() )

            #xf_to_pos.multiply ( xf_to_com )
            xf_to_pos.multiply ( xf_rotR )
            xf_to_pos.multiply ( xf_to_P )
            xf_to_pos.multiply ( xf_rot0 )
            xf_to_pos.multiply ( xf_to_0 )

            amod.mod.openState.xform = xf_to_pos

            if hasattr (amod.mod, 'surfMods') :
                for cid, mod in amod.mod.surfMods.iteritems() :
                    mod.openState.xform = xf_to_pos




class RockAts (Frames) :

    def __init__ ( self, startStep, endStep, animMods, comMod, selStr, axis, totDeg, numCycles, itype="cubic" ) :
        super(RockAts, self).__init__(startStep, endStep)
        #print " - rock Ats"
        #print "   axis: ", axis
        #print "   totDeg: ", totDeg
        #print "   numCyc: ", numCycles

        self.animMods = animMods
        self.totDeg = float ( totDeg )
        self.N = float ( numCycles )
        self.axis = chimera.Vector ( axis[0], axis[1], axis[2] )
        self.comMod = comMod
        #self.comlc = chimera.Point ( COM[0], COM[1], COM[2] )
        self.selStr = selStr
        self.itype = itype


    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        from quaternion import Quaternion

        if stepAt == self.start :
            sel = chimera.selection.OSLSelection ( self.selStr )

            from _multiscale import get_atom_coordinates
            points = get_atom_coordinates ( sel.atoms(), transformed = True )
            COM, U, S, V = prAxes ( points )
            comp = chimera.Point ( COM[0], COM[1], COM[2] )

            #print " - rock Ats - first step - %s, " % self.selStr, COM

            self.comwc = comp # self.comMod.mod.openState.xform.apply ( comp )
            for amod in self.animMods :
                amod.xf0 = amod.mod.openState.xform
                amod.q0 = Quaternion ()
                amod.q0.fromXform ( amod.xf0 )
                amod.comlc = chimera.Point ( amod.COM[0], amod.COM[1], amod.COM[2] )
                amod.comwcv = amod.xf0.apply ( amod.comlc ).toVector() - self.comwc.toVector()

        #super(Frames, self).step(stepAt)

        super(RockAts, self).step(stepAt)
        f = self.f
        f1, f2 = 2.0*f*f*f-3.0*f*f+1.0, 3*f*f-2*f*f*f         # cubic interpolation
        if not "cubic" == self.itype :
            f1, f2 = (1.0-f), f                                   # linear interpolation

        #deg = self.totDeg * f2
        deg = self.totDeg * numpy.sin ( f2 * 2.0 * numpy.pi * self.N )
        #print " - at deg ", deg

        for amod in self.animMods :

            xf_to_0 = chimera.Xform.translation ( amod.comlc.toVector() * -1.0 )
            xf_to_P = chimera.Xform.translation ( amod.comwcv )
            xf_rot0 = amod.q0.Xform ()
            xf_rotR = chimera.Xform.rotation ( self.axis, deg )
            xf_to_pos = chimera.Xform.translation ( self.comwc.toVector() )

            #xf_to_pos.multiply ( xf_to_com )
            xf_to_pos.multiply ( xf_rotR )
            xf_to_pos.multiply ( xf_to_P )
            xf_to_pos.multiply ( xf_rot0 )
            xf_to_pos.multiply ( xf_to_0 )

            amod.mod.openState.xform = xf_to_pos

            if hasattr (amod.mod, 'surfMods') :
                for cid, mod in amod.mod.surfMods.iteritems() :
                    mod.openState.xform = xf_to_pos






class Scale (Frames) :

    def __init__ ( self, startStep, endStep, targetScale ) :

        print " - xf scale"
        super(Scale, self).__init__(startStep, endStep)

        self.targetScale = targetScale
        self.scaleSign = 1.0
        self.scaleIncr = targetScale / (endStep - startStep)



    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        print " - scale step %f" % self.targetScale
        #super(Frames, self).step(stepAt)

        cmd = "scale %f" % self.targetScale
        chimera.runCommand ( cmd )




class Cycle (Frames) :

    def __init__ ( self, startStep, endStep, animMods, dstep ) :
        print " - Cycle M"
        print "   dstep: ", dstep
        print "   num models: ", len(animMods)
        super(Cycle, self).__init__(startStep, endStep)

        self.animMods = animMods
        self.dstep = dstep
        self.atMod = 0


    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        super(Cycle, self).step(stepAt)
        #f = self.f
        #rangef = float (self.end - self.start)

        import numpy
        modf = numpy.floor ( float(self.atMod) / float(self.dstep) )
        modi = int ( modf ) % len(self.animMods)

        #print " step %d (%d -> %d) d %d, modi: %d" % (stepAt, self.start, self.end, self.dstep, modi)

        for ai, amod in enumerate (self.animMods) :

            if ai == modi :
                amod.mod.display = True
            else :
                amod.mod.display = False

        self.atMod += 1



# chimera.runCommand ( "clip hither -2" )

class Clip (Frames) :

    def __init__ ( self, startStep, endStep, totalC ) :
        print " - Clip"
        print "   totalC: ", totalC
        super(Clip, self).__init__(startStep, endStep)

        self.totalC = totalC
        self.dc = float(totalC) / float(endStep-startStep+1)
        self.at = 0


    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        super(Clip, self).step(stepAt)
        self.at += self.dc

        #print " clip step %d (%d -> %d) d %.5f, total %d" % (stepAt, self.start, self.end, self.dc, self.totalC)
        print " clip %d - %d/%d" % (stepAt, self.at, self.totalC)
        chimera.runCommand ( "clip hither %.5f" % self.dc )



class ClipOff (Frames) :

    def __init__ ( self, startStep ) :
        print " - Clip Off"
        super(ClipOff, self).__init__(startStep, startStep)

        self.done = False


    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        super(ClipOff, self).step(stepAt)

        if self.done :
            return

        print " clip off step %d (%d -> %d)" % (stepAt, self.start, self.end)
        chimera.runCommand ( "clip off" )
        self.done = True




class VolMorph (Frames) :

    # er_dna

    def __init__ ( self, startStep, endStep, startM, endM ) :
        print " - morph %s, %.1f -> %s, %.1f" % ( startM.mod.name, startM.mod.surface_levels[0], endM.mod.name, endM.mod.surface_levels[0] )
        super(VolMorph, self).__init__(startStep, endStep)
        self.startM = startM
        self.endM = endM
        self.df_v = None


    def step ( self, stepAt ) :

        if stepAt < self.start or stepAt > self.end :
            return

        from quaternion import Quaternion


        if stepAt == self.start :
            print " - morph - first step"

            startM = self.startM.mod
            endM = self.endM.mod

            self.start_mat = startM.data.full_matrix()
            self.startSurfaceLevel = startM.surface_levels[0]
            f_mask = numpy.where ( self.start_mat > startM.surface_levels[0], numpy.ones_like(self.start_mat), numpy.zeros_like(self.start_mat) )
            self.startVol = numpy.sum ( f_mask ) * startM.data.step[0] * startM.data.step[1] * startM.data.step[2]
            print " - start %s thr %.3f vol: %.3f" % (startM.name, startM.surface_levels[0], self.startVol)
            print startM.surface_levels

            self.startColor = numpy.array ( startM.surfacePieces[0].color )
            #print " - start color: ", self.startColor

            self.end_mat = endM.data.full_matrix().copy()
            self.endSurfaceLevel = endM.surface_levels[0]
            f_mask = numpy.where ( self.end_mat > endM.surface_levels[0], numpy.ones_like(self.end_mat), numpy.zeros_like(self.end_mat) )
            self.endVol = numpy.sum ( f_mask ) * endM.data.step[0] * endM.data.step[1] * endM.data.step[2]
            print " - end %s thr %.3f vol: %.3f" % (endM.name, endM.surface_levels[0], self.endVol)
            print endM.surface_levels

            self.endColor = numpy.array ( endM.surfacePieces[0].color )
            #print " - end color: ", self.endColor

            fmap = endM
            dmap = startM
            import _contour
            n1, n2, n3 = fmap.data.size[0], fmap.data.size[1], fmap.data.size[2]
            f_points = VolumeData.grid_indices( (n1,n2,n3), numpy.single )  # i,j,k indices
            _contour.affine_transform_vertices( f_points, fmap.data.ijk_to_xyz_transform )

            d_vals = dmap.interpolated_values ( f_points, fmap.openState.xform )
            self.start_mat = d_vals.reshape( (n3,n2,n1) )

            #f_mask = numpy.where ( self.end_mat > endM.surface_levels[0], numpy.ones_like(self.end_mat), numpy.zeros_like(self.end_mat) )
            #self.endVol = numpy.sum ( f_mask ) * endM.data.step[0] * endM.data.step[1] * endM.data.step[2]
            #print " - end vol after interp: %.3f" % self.endVol

            self.endM.mod.morphMod = endM.writable_copy(require_copy = True, name = endM.name+'__morph')
            startM.display = False
            endM.display = False


        #try :
        #    chimera.openModels.close ( self.df_v )
        #except :
        #    pass

        super(VolMorph, self).step(stepAt)
        f = self.f
        #f1, f2 = 2.0*f*f*f-3.0*f*f+1.0, 3*f*f-2*f*f*f         # cubic interpolation
        #if not "cubic" == self.itype :
        f1, f2 = (1.0-f), f                                   # linear interpolation
        #print f1, f2

        df_mat = self.start_mat * f1 + self.end_mat * f2

        #M = self.df_v.data.full_matrix()
        #M[:,:,:] = df_mat[:,:,:]
        #self.df_v.data.values_changed()

        morphMod = self.endM.mod.morphMod

        morphMod.data.full_matrix()[:,:,:] = df_mat[:,:,:]
        morphMod.data.values_changed()

        #sf = self.df_v.surface_level_for_enclosed_volume ( self.endVol )
        sf = self.startSurfaceLevel * f1 + self.endSurfaceLevel * f2
        #print " - new surf level for end vol: %.2f" % sf
        morphMod.surface_levels = [sf]

        col = self.startColor * f1 + self.endColor * f2

        ro = VolumeViewer.volume.Rendering_Options()
        morphMod.update_surface ( False, ro )
        for sp in morphMod.surfacePieces :
            v, t = sp.geometry
            if len(v) == 8 and len(t) == 12 :
                sp.display = False
            if len(v) == 0 and len(t) == 0 :
                sp.display = False
            else :
                sp.color = (col[0],col[1],col[2],col[3])


        if stepAt == self.end :
            #print " - last step!"
            self.endM.mod.display = True
            chimera.openModels.close ( [self.endM.mod.morphMod] )
            del self.endM.mod.morphMod










# -----------------------------------------------------------------------------


# the active movie as global to avoid passing this around
bioMovieActiveMovie = None

def AddAction (action) :
    global bioMovieActiveMovie
    if bioMovieActiveMovie == None :
        print "No movie created, hence action was not added to any movie"
        print " - create a movie first, e.g. movie = biomovie.Movie(biomovie_dialog) "
    else :
        bioMovieActiveMovie.add ( action )


class Movie :


    def __init__ (self, dlg ) :
        self.anims = []
        self.start = 10000000
        self.end = 0
        self.dlg = dlg
        self.keys = {}

        global bioMovieActiveMovie
        bioMovieActiveMovie = self


    def add ( self, anim ) :
        if anim.start < self.start :
            self.start = anim.start
        if anim.end > self.end :
            self.end = anim.end
        self.anims.append ( anim )


    def addKey ( self, framei, text ) :
        self.keys[framei] = text


    def make (self, saveMovie = True) :

      saveMovie = self.dlg.makeMovie.get()
      stop = self.dlg.stopMovie.get()

      if saveMovie :
        self.dlg.ClearFrames ()

      fri = 0

      for i in range ( self.start, self.end+1 ) :

        saveMovie = self.dlg.makeMovie.get()
        stop = self.dlg.stopMovie.get()
        #print "%d - " % i

        if stop :
            print "Stopped"
            break

        for anim in self.anims :
            anim.step ( i )


        chimera.viewer.postRedisplay()
        self.dlg.toplevel_widget.update_idletasks ()


        if i in self.keys :
            self.dlg.SaveKeyFrame ( i, self.keys[i] )

        if i % 30 == 0 :
            print "%d/%d" % (i+1,self.end+1)
        else :
            print ".",

        if saveMovie :
            self.dlg.SaveFrame ( fri, 1 )
            fri += 1

      if saveMovie and not self.dlg.stopMovie.get() :
          self.dlg.MakeMovie ( self.dlg.movieName.get() )


    def run (self) :
        self.make ( False )






class AnimatableModel :

    def __init__ (self) :
        self.mod = None


    def FromMod ( self, fm ) :

        self.mod = fm
        self.xf0 = fm.openState.xform
        self.xf = fm.openState.xform
        #print "New AM:", fm.name

        if type(fm) == VolumeViewer.volume.Volume :
            self.type = "map"
        elif type(fm) == chimera.Molecule :
            self.type = "mol"

        if 0 and hasattr ( self.mod, "COM" ) :
            self.COM = self.mod.COM
            self.comPt = chimera.Point ( self.COM[0], self.COM[1], self.COM[2] )
            self.comVec = self.comp.toVector ()
            self.comPt_wc = self.mod.openState.xform.apply ( self.comPt )
            print " - ", self.comPt_wc


    def FromMap ( self ) :

        self.pts, self.weights = map_points ( self.mod )
        if len(self.pts) == 0 :
            self.pts, self.weights = map_points ( self.mod, False )
        #print len(self.pts)
        self.COM, self.U, self.S, self.V = prAxes ( self.pts )
        #print " - " + self.mod.name + ", COM : ", self.COM

        self.comPt = chimera.Point ( self.COM[0], self.COM[1], self.COM[2] )
        self.comVec = self.comPt.toVector ()

        self.comPt_wc = self.mod.openState.xform.apply ( self.comPt )
        self.COMWC = [self.comPt_wc[0], self.comPt_wc[1], self.comPt_wc[2]]

        return self




    def FromMol (self) :

        #if hasattr ( self, "COM" ) :
        #    return

        sel = chimera.selection.OSLSelection ( "#%d" % self.mod.id )
        atoms = sel.atoms()

        #print self.mod.name, " - ", len(atoms), "atoms"

        from _multiscale import get_atom_coordinates
        points = get_atom_coordinates ( atoms, transformed = False )
        self.COM, self.U, self.S, self.V = prAxes ( points )

        self.comPt = chimera.Point ( self.COM[0], self.COM[1], self.COM[2] )
        self.comVec = self.comPt.toVector ()
        self.comPt_wc = self.mod.openState.xform.apply ( self.comPt )

        self.dispMode = ["ribbon"]

        # print " - com: ", self.comp

        return self




    def FromSurf (self) :

        print self.mod.name, " - from surf"

        self.COM = numpy.array ( [ 0,0,0 ], numpy.float32 )
        N = 0.0;
        rad = 0.0;

        for sp in self.mod.surfacePieces :
            for p in sp.geometry[0] :
                self.COM = self.COM + p;
                N = N + 1.0;
                r = numpy.sqrt ( (p**2).sum() )
                if r > rad :
                    rad = r

        self.COM = self.COM / N;
        self.comPt = chimera.Point ( self.COM[0], self.COM[1], self.COM[2] )
        self.comVec = self.comPt.toVector ()
        # print " - com: ", self.comp

        return self




    def Rotate00 ( self, deg=5.0, center = None, axis = [0,0,1] ) :

        if ( center == None ) :
            center = self.COM

        rxf = chimera.Xform.rotation ( chimera.Vector(axis[0],axis[1],axis[2]), deg )
        txf0 = chimera.Xform.translation ( chimera.Vector(-center[0],-center[1],-center[2]) )
        txf = chimera.Xform.translation ( chimera.Vector(center[0],center[1],center[2]) )

        self.xf.multiply ( txf )
        self.xf.multiply ( rxf )
        self.xf.multiply ( txf0 )
        self.mod.openState.xform = self.xf


    def Show ( self ) :

        self.mod.display = True

        if hasattr (self.mod, 'surfMods') :
            #self.mod.display = False
            self.mod.display = True
            for cid, mod in self.mod.surfMods.iteritems() :
                try :
                    self.mod.surfMods[cid].display = True
                except :
                    pass


    def Hide ( self ) :

        self.mod.display = False

        if hasattr (self.mod, 'surfMods') :
            for cid, mod in self.mod.surfMods.iteritems() :
                try :
                    self.mod.surfMods[cid].display = False
                except :
                    pass



    def SetSurfThr ( self, thr ) :

        if self.mod.display == False :
            return

        if type ( self.mod ) != VolumeViewer.volume.Volume :
            return

        self.SetModThr ( self.mod, thr )
        if hasattr ( self, 'toAlpha' ) :
            self.SetAlpha ( self.toAlpha )



    def SetSurfColor ( self, r, g, b, a ) :

        self.SetModSurfColor ( self.mod, (r,g,b,a) )



    def SetModThr ( self, mod, thr ) :

        if mod.display == False :
            return

        if type ( mod ) != VolumeViewer.volume.Volume :
            return

        mod.region = ( mod.region[0], mod.region[1], [1,1,1] )
        mod.surface_levels[0] = thr

        ro = VolumeViewer.volume.Rendering_Options()
        ro.smoothing_factor = .2
        ro.smoothing_iterations = 5
        ro.surface_smoothing = True

        mod.update_surface ( False, ro )

        for sp in mod.surfacePieces :
            v, t = sp.geometry
            if len(v) == 8 and len(t) == 12 :
                sp.display = False


    def SetModSurfColor ( self, mod, clr ) :

        for sp in mod.surfacePieces :
            v, t = sp.geometry
            if len(v) == 8 and len(t) == 12 :
                sp.display = False
            else :
                sp.color = clr


    def SetMapDisplay ( disp ) :

        for sp in mod.surfacePieces :
            v, t = sp.geometry
            if len(v) == 8 and len(t) == 12 :
                sp.display = False
            else :
                if disp == "mesh" :
                    sp.displayStyle = sp.Mesh
                    sp.lineThickness = 2.0
                else :
                    sp.displayStyle = sp.Solid




    def SetAlpha ( self, a ) :

        import Segger
        import Segger.regions

        if type(self.mod) == Segger.regions.Segmentation :
            print "seg set a"

            for r in self.mod.regions :
                if r.has_surface():
                    cr,cg,cb = r.surface_piece.color[:3] #r.color[:3]
                    r.surface_piece.color = ( cr, cg, cb, a )

        elif type(self.mod) == VolumeViewer.volume.Volume :

            for sp in self.mod.surfacePieces :
                v, t = sp.geometry
                if len(v) == 8 and len(t) == 12 :
                    sp.display = False
                else :
                    c = sp.color
                    sp.color = ( c[0], c[1], c[2], a )
                    #sp.vertexColors = None

                    if hasattr ( sp, "vertexColors" ) and sp.vertexColors != None :
                        vcolors = []
                        for vc in sp.vertexColors :
                            vcolors.append ( (vc[0], vc[1], vc[2], a) )

                        sp.vertexColors = vcolors

        elif type(self.mod) == chimera.Molecule :
            if hasattr (self.mod, 'surfMods') :
                #print " - set a: ", self.mod.name, a
                for cid, mod in self.mod.surfMods.iteritems() :
                    #chimera.openModels.close ( [mod] )
                    color = self.colors[cid]
                    self.SetModSurfColor ( self.mod.surfMods[cid], (color[0], color[1], color[2], a) )
                    self.alphaAt = a
                    #print cid,
                #print ""

        elif type(self.mod) == _surface.SurfaceModel :
            for sp in self.mod.surfacePieces :
                v, t = sp.geometry
                if len(v) == 8 and len(t) == 12 :
                    sp.display = False
                else :
                    c = sp.color
                    sp.color = ( c[0], c[1], c[2], a )


    def SetXform ( self, f ) :

        from quaternion import Quaternion
        t0 = self.xf0.getTranslation ()
        q0 = Quaternion ()
        q0.fromXform ( self.xf0 )

        t1 = self.xf.getTranslation ()
        q1 = Quaternion ()
        q1.fromXform ( self.xf )


        # linear interpolation
        f1, f2 = (1.0-f), f

        # cubic interpolation
        f1, f2 = 2.0*f*f*f-3.0*f*f+1.0, 3*f*f-2*f*f*f

        pos = t0 * f1  + t1 * f2
        s = q0.s * f1  + q1.s * f2
        v = q0.v * f1  + q1.v * f2

        Q = Quaternion ( s, v )
        Q.normalize()

        xf = Q.Xform ()
        # print "- com pos: ", pos

        #tr0 = chimera.Xform.translation ( -sm.COM )
        tr = chimera.Xform.translation ( pos )

        #xf.multiply ( tr0 )
        xf.premultiply ( tr )

        #if rxf1.ref_mod :
        #    xf.premultiply ( rxf1.ref_mod.openState.xform )

        self.mod.openState.xform = xf


    def UpdateDisp (self) :

        if not hasattr(self, 'dispMode') or self.dispMode == None :
            return

        if type( self.mod ) == chimera.Molecule :

            if self.dispMode[0] == "ribbon" :

                #self.mod.display = True

                if hasattr (self.mod, 'surfMods') :
                    for cid, smod in self.mod.surfMods.iteritems() :
                        chimera.openModels.close ( [smod] )

                self.mod.surfMods = {}

            elif self.dispMode[0] == "surf" :

                #self.mod.display = False
                #self.mod.display = True

                if hasattr (self.mod, 'surfMods') :
                    for cid, smod in self.mod.surfMods.iteritems() :
                        chimera.openModels.close ( [smod] )

                self.mod.surfMods = {}

                if hasattr(self, 'colors') :
                    for cid, color in self.colors.iteritems() :

                        modName = self.mod.name + " -surface.for.chain- " + cid
                        closeMods = []
                        for m in chimera.openModels.list() :
                            if m.name == modName :
                                closeMods.append ( m )
                                break
                        if len(closeMods) > 0 :
                            chimera.openModels.close ( closeMods )

                        res, step, thr = self.dispMode[1], self.dispMode[2], self.dispMode[3]
                        self.mod.surfMods[cid] = self.GenStrucMap ( cid, step, res )
                        self.mod.surfMods[cid].name = modName
                        self.SetModThr ( self.mod.surfMods[cid], thr )
                        self.alphaAt = color[3]
                        self.SetModSurfColor ( self.mod.surfMods[cid], (color[0], color[1], color[2], self.alphaAt) )
                        #print " %d" % len(self.mod.surfMods[cid].surfacePieces),
                    #print "."



    def GenStrucMap ( self, cid, step, res ) :

        #cmd = "molmap #%s %f sigmaFactor 0.187 gridSpacing %f replace false" % ( mol.id, res, step )
        cmd = "molmap #%s:.%s %f gridSpacing %f replace false" % ( self.mod.id, cid, res, step )
        #print " -", cmd
        chimera.runCommand ( cmd )

        mv = None
        for mod in chimera.openModels.list() :
            ts = mod.name.split()
            if len(ts) > 1 and mod.name.find("map") >=0 and mod.name.find("res") >=0 :
                #print " - found", mod.name
                mv = mod
                break

        if mv == None :
            print "- molmap not found..."
            return None

        return mv





# calculates 'principal axes' of a set of points
def prAxes ( points ) :

    com = numpy.sum(points, axis=0) / len(points)
    C = chimera.Vector ( com[0], com[1], com[2] )

    comv = numpy.ones_like ( points ) * com
    points = points - comv

    i = numpy.matrix ( [[1,0,0], [0,1,0], [0,0,1]] )
    ii = i * numpy.sum ( numpy.multiply ( points, points ) )
    p_t = numpy.transpose(points)
    td = numpy.tensordot ( points, p_t, axes=[0,1] )

    I0 = ii - td

    try :
        U, S, V = numpy.linalg.svd( I0 )
    except :
        print "- error computing SVD - prob. singular matrix"
        return []

    #U[0,0] = U[0,0] * -1.0
    #U[1,0] = U[1,0] * -1.0
    #U[2,0] = U[2,0] * -1.0

    #U[0,2] = U[0,2] * -1.0
    #U[1,2] = U[1,2] * -1.0
    #U[2,2] = U[2,2] * -1.0

    return [C, U, S, V]



# returns grid points in the map above a given threshold value
def map_points (fmap, useThreshold = True):

    from _contour import affine_transform_vertices as transform_vertices

    mat = fmap.data.full_matrix()
    threshold = fmap.surface_levels[0]

    if useThreshold == False :
        #threshold = -1e9
        threshold = 1e-5
        #print " - not using threshold"

    import _volume
    points = _volume.high_indices(mat, threshold)
    fpoints = points.astype(numpy.single)
    fpoint_weights = mat[points[:,2],points[:,1],points[:,0]]

    nz = numpy.nonzero( fpoint_weights )[0]
    if len(nz) < len (fpoint_weights) :
        fpoints = numpy.take( fpoints, nz, axis=0 )
        fpoint_weights = numpy.take(fpoint_weights, nz, axis=0)

    transform_vertices( fpoints, fmap.data.ijk_to_xyz_transform )

    if 0 : print "FitPoints from %s with threshold %.4f, %d nonzero" % (
        fmap.name, threshold, len(nz) )

    return fpoints, fpoint_weights





# ---------------------------------------------------------------------------------------------------------


def get_dialog () :

    from chimera import dialogs
    d = dialogs.find ( "BioMovie", create=False )
    return d



def close_dialog ():

	from chimera import dialogs
	d = dialogs.find ( "BioMovie", create=False )

	if d :
		print " - found dialog"
		d.toplevel_widget.update_idletasks ()
		d.Close()
		d.toplevel_widget.update_idletasks ()
	else :
		print " - did not find dialog"



def show_dialog ():

    close_dialog()

    from chimera import dialogs
    dialogs.register ("BioMovie", BioMovie, replace = True)
    d = dialogs.find ( "BioMovie", create=True )

    # Avoid transient dialog resizing when created and mapped for first time.
    d.toplevel_widget.update_idletasks ()
    d.enter()

    return d



def getMod ( name ) :
    for mol in chimera.openModels.list () :
        if mol.name == name :
            return mol
    return None

def getModById ( id ) :
    for mol in chimera.openModels.list () :
        try : mol.id
        except : continue
        if mol.id == id :
            return mol
    return None

def visMods () :
    for mol in chimera.openModels.list () :
        try :
            mol.shown()
        except :
            continue
        if mol.shown() == True :
            print mol.id, " ", mol.name

    return None
