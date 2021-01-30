
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

from axes import prAxes
import regions
import graph
from Segger import dev_menus, timing, seggerVersion

OML = chimera.openModels.list

REG_OPACITY = 0.45


from segment_dialog import current_segmentation, segmentation_map



def umsg ( txt ) :
    print txt
    status ( txt )


def status ( txt ) :
    txt = txt.rstrip('\n')
    msg.configure(text = txt)
    msg.update_idletasks()




class ProMod_Dialog ( chimera.baseDialog.ModelessDialog ):

    title = "ProMod - Probabilistic Models (Segger v" + seggerVersion + ")"
    name = "segger_promod"
    buttons = ( "Close" )
    help = 'https://github.com/gregdp/segger'

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
            l = Tkinter.Label(ff, text = "1. Open all models to be considered, make them visible, hide other models", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)



        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = "2. Find (closest-to) average model", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)


            b = Tkinter.Button(ff, text="Find Average Model", command=self.AvgMod)
            b.grid (column=1, row=0, sticky='w', padx=5, pady=1)


            self.avgModLabel = Tkinter.Label(ff, text = " ", anchor = 'w')
            self.avgModLabel.grid(column=2, row=0, sticky='ew', padx=5, pady=1)



        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = "3. Calculate standard deviations at each residue ", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)


            b = Tkinter.Button(ff, text="Calculate", command=self.Calc)
            b.grid (column=1, row=0, sticky='w', padx=5, pady=1)


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = " - standard deviations are stored for each residue atom as the b-factor", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)



        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = " - use Tools -> Depiction -> Render by Attribute to show deviations using", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        if 1 :
            l = Tkinter.Label(ff, text = "    color and/or ribbon thickness. See tutorial by pressing Help below.", anchor = 'w')
            l.grid(column=0, row=0, sticky='ew', padx=5, pady=1)



        row += 1
        f = Tkinter.Frame(parent)
        f.grid(column=0, row=row, sticky='ew')
        l = Tkinter.Label(f, text='  ')
        l.grid(column=0, row=row, sticky='w')


        row += 1
        dummyFrame = Tkinter.Frame(parent, relief='groove', borderwidth=1)
        Tkinter.Frame(dummyFrame).pack()
        dummyFrame.grid(row=row,column=0,columnspan=7, pady=7, sticky='we')


        global msg
        row = row + 1
        msg = Tkinter.Label(parent, width = 60, anchor = 'w', justify = 'left', fg="red")
        msg.grid(column=0, row=row, sticky='ew', padx=5, pady=1)
        row += 1






    def Calc ( self ) :


        if hasattr ( self, 'avgMod' ) and hasattr ( self, 'mods' ) and len(self.mods) > 0 and self.avgMod != None :
            print "Average model: %s -- %d mods" % ( self.avgMod.name, len(self.mods) )
        else :
            umsg ("Find Average Model first.")
            return


        avgMod = self.avgMod
        mods = self.mods

        umsg ( "Calculating standard deviations..." )

        vars = []

        for ri, avgRes in enumerate ( avgMod.residues ) :


            status ( "Res %d/%d" % (ri+1,len(avgMod.residues)) )


            for avgAt in avgRes.atoms :

                mean = 0.0

                for m in mods :
                    res = m.residues[ri]
                    cat = res.atomsMap[avgAt.name][0]
                    v = cat.coord() - avgAt.coord()
                    d = v.length * v.length
                    mean += d

                mean /= len(mods)
                stdev = numpy.sqrt ( mean )
                vars.append ( stdev )

                for m in mods :
                    res = m.residues[ri]
                    cat = res.atomsMap[avgAt.name][0]
                    cat.bfactor = stdev


        umsg ( "%d models, %d residues - min variance %.2f, max variance %.2f" % (
                    len(mods), len(avgMod.residues), numpy.min(vars), numpy.max(vars) ) )



    def Calc_CA ( self ) :


        if hasattr ( self, 'avgMod' ) and hasattr ( self, 'mods' ) and len(self.mods) > 0 and self.avgMod != None :
            print "Average model: %s -- %d mods" % ( self.avgMod.name, len(self.mods) )
        else :
            umsg ("Find Average Model first.")
            return


        avgMod = self.avgMod
        mods = self.mods

        umsg ( "Calculating standard deviations..." )

        vars = []

        for ri, resAvg in enumerate ( avgMod.residues ) :
            try :
                catAvg = resAvg.atomsMap["CA"][0]
            except :
                continue


            mean = 0.0

            for m in mods :
                res = m.residues[ri]
                cat = res.atomsMap["CA"][0]
                v = cat.coord() - catAvg.coord()
                d = v.length * v.length
                mean += d

            mean /= len(mods)
            stdev = numpy.sqrt ( mean )

            vars.append ( stdev )

            for m in mods :
                res = m.residues[ri]
                for at in res.atoms :
                    at.bfactor = stdev
                    #at.occupancy = stdev


        umsg ( "%d models, %d residues - min variance %.2f, max variance %.2f" % (
                    len(mods), len(avgMod.residues), numpy.min(vars), numpy.max(vars) ) )




    def AvgMod0 ( self ) :

        self.avgMod = None
        self.mods = []
        import numpy

        for m in chimera.openModels.list() :
            if type (m) == chimera.Molecule and m.display == True:
                self.mods.append ( m )

        N = len(self.mods)

        if N < 2 :
            umsg ( "At least 2 models are needed - make sure they are shown" )
            self.avgModLabel.configure ( text = "" )
            return



        mod0 = self.mods[0]
        numRes = len(mod0.residues)

        umsg ( "Finding average of %d mods, %d residues" % ( len(self.mods), len(mod0.residues) ) )

        avgPs = numpy.zeros ( [len(mod0.residues), 3] )

        for mod in self.mods :
            #print " - mod: %s, %d residues" % ( mod.name, len(mod.residues) )

            if numRes <> len(mod.residues) :
                umsg ("All models should have the same number of residues")
                self.avgModLabel.configure ( text = "" )
                return

            for ri, res in enumerate ( mod.residues ) :
                cat = None
                try :
                    cat = res.atomsMap["CA"][0]
                except :
                    #print "carbon alpha not found in res ", ri, res.id.position
                    #return None
                    pass

                if cat :
                    avgPs[ri] += cat.coord().data()


        N = float ( len(self.mods) )
        for ri, res in enumerate ( mod0.residues ) :
            avgPs[ri] /= N

            #if ri == 0 :
            #    print " r0 avg pos: ", avgPs[ri]



        minDist = -1.0
        minMod = None

        for mod in self.mods :

            #print " - mod: %s, %d residues" % ( mod.name, len(mod.residues) ),
            modDist = 0.0

            for ri, res in enumerate ( mod.residues ) :
                try :
                    cat = res.atomsMap["CA"][0]
                except :
                    #print "carbon alpha not found in mod %s res " % mod.name, ri, res.id.position
                    #return None
                    continue

                dv = avgPs[ri] - cat.coord().data()
                modDist += numpy.sum ( dv * dv )

            #print ", dist: ", modDist

            if minMod == None or modDist < minDist :
                minMod = mod
                minDist = modDist

        print "Avg mod: %s, min dist to avg: %.2f" % (minMod.name, minDist)

        self.avgMod = minMod

        self.avgModLabel.configure ( text = " found: %s" % minMod.name )
        umsg ( "Average of %d models is %s" % (len(self.mods), minMod.name) )


        return minMod, avgPs





    def AvgMod ( self ) :

        self.avgMod = None
        self.mods = []
        import numpy

        for m in chimera.openModels.list() :
            if type (m) == chimera.Molecule and m.display == True:
                self.mods.append ( m )

        N = len(self.mods)

        if N < 2 :
            umsg ( "At least 2 models are needed - make sure they are shown" )
            self.avgModLabel.configure ( text = "" )
            return



        mod0 = self.mods[0]
        numRes = len(mod0.residues)

        umsg ( "Finding average of %d mods, %d residues" % ( len(self.mods), len(mod0.residues) ) )
        print "."

        #avgPs = numpy.zeros ( [len(mod0.atoms), 3] )
        avg = {}

        for mod in self.mods :
            #print " - mod: %s, %d residues" % ( mod.name, len(mod.residues) )

            for res in mod.residues :
                for at in res.atoms :
                    if not res.id.chainId in avg :
                        avg[res.id.chainId] = {}
                    if not res.id.position in avg[res.id.chainId] :
                        avg[res.id.chainId][res.id.position] = {}
                    if not at.name in avg[res.id.chainId][res.id.position] :
                        avg[res.id.chainId][res.id.position][at.name] = []

                    avg[res.id.chainId][res.id.position][at.name].append ( numpy.array ( at.coord().data() ) )


        for ci, rmap in avg.iteritems () :
            for ri, amap in rmap.iteritems () :
                for aname, plist in amap.iteritems () :
                    if len(plist) <> len(self.mods) :
                        print " - at %s_%d.%s has only %d/%d pos" % ( aname, ri, ci, len(plist), len(self.mods) )

                    avgp = numpy.array ( [0,0,0] )
                    for p in plist :
                        avgp += p
                    avgp /= float ( len(plist) )



        minDist = -1.0
        minMod = None

        for mod in self.mods :

            #print " - mod: %s, %d residues" % ( mod.name, len(mod.residues) ),
            modDist = 0.0

            for ri, res in enumerate ( mod.residues ) :
                for at in res.atoms :
                    avgPos = avg[res.id.chainId][res.id.position][at.name]
                    dv = numpy.array ( at.coord().data() ) - avgPos
                    modDist += numpy.sum ( dv * dv )

            #print ", dist: ", modDist

            if minMod == None or modDist < minDist :
                minMod = mod
                minDist = modDist

        print "Avg mod: %s, min dist to avg: %.2f" % (minMod.name, minDist)

        self.avgMod = minMod

        self.avgModLabel.configure ( text = " found: %s" % minMod.name )
        umsg ( "Average of %d models is %s" % (len(self.mods), minMod.name) )


        return minMod



def Bring () :

    print "bring..."
    fromm, tom = None, None
    for m in chimera.openModels.list() :
        if type (m) == chimera.Molecule and m.display == True:
            if "promod" in m.name :
                fromm = m
            else :
                tom = m

    print " - from: %s" % fromm.name
    print " -   to: %s" % tom.name

    bfs = []
    rid = {}
    for r in fromm.residues :
        rid[r.id.position] = r
        for at in r.atoms :
            bfs.append ( at.bfactor )

    print "devs mean: %.3f" % numpy.average(bfs)
    print "devs std: %.3f" % numpy.std(bfs)
    print "devs 3sig: %.3f" % (numpy.average(bfs) + 3.0*numpy.std(bfs))

    for r in tom.residues :
        rf = rid[r.id.position]
        for at in r.atoms :
            at.bfactor = rf.atomsMap[at.name][0].bfactor


def show_dialog (closeOld = True):

    from chimera import dialogs

    d = dialogs.find ( ProMod_Dialog.name, create=False )
    if d :
        if closeOld :
            d.toplevel_widget.update_idletasks ()
            d.Close()
            d.toplevel_widget.update_idletasks ()
        else :
            return d

    dialogs.register ( ProMod_Dialog.name, ProMod_Dialog, replace = True)

    d = dialogs.find ( ProMod_Dialog.name, create=True )
    # Avoid transient dialog resizing when created and mapped for first time.
    d.toplevel_widget.update_idletasks ()
    d.enter()

    return d



# -----------------------------------------------------------------------------
#
