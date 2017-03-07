
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


class Extract_Region_Dialog ( chimera.baseDialog.ModelessDialog ):

    title = "Extract Densities (Segger v" + seggerVersion + ")"
    name = "extract region"
    
    if dev_menus :
        buttons = ('EQ', 'Extract', "Close")
    else :
        buttons = ('Extract', "Close")
    
    
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

        #l = Tkinter.Label(f, text='  ')
        #l.grid(column=0, row=row, sticky='w')
        #row += 1


        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')

        l = Tkinter.Label(ff, text='  Extract densities from map:')
        l.grid(column=0, row=0, sticky='w')

        self.dmap = Tkinter.StringVar(parent)

        self.mb  = Tkinter.Menubutton ( ff, textvariable=self.dmap, relief=Tkinter.RAISED )
        self.mb.grid (column=1, row=0, sticky='we', padx=5)
        self.mb.menu  =  Tkinter.Menu ( self.mb, tearoff=0, postcommand=self.MapMenu )
        self.mb["menu"]  =  self.mb.menu


        # set the segmentation map (if any) as the initial default
        self.cur_dmap = segmentation_map()
        st1, st2, st3 = "", "", ""
        d1, d2, d3 = "", "", ""
        bwidth = ""
        if self.cur_dmap :
            # print "Current segmentation map: ", self.cur_dmap.name
            self.dmap.set ( self.cur_dmap.name )
            st = self.cur_dmap.data.step
            
            st1 = "%g" % self.cur_dmap.data.step[0]
            st2 = "%g" % self.cur_dmap.data.step[1]
            st3 = "%g" % self.cur_dmap.data.step[2]

            d1 = "%g" % self.cur_dmap.data.size[0]
            d2 = "%g" % self.cur_dmap.data.size[1]
            d3 = "%g" % self.cur_dmap.data.size[2]
        
            s1, s2, s3 = self.cur_dmap.data.step[0], self.cur_dmap.data.step[1], self.cur_dmap.data.step[2]
            bwidth = "%.0f" % (numpy.sqrt(s1*s1 + s2*s2 + s3*s3)*5.0)


        #row += 1
        #l = Tkinter.Label(f, text='  ')
        #l.grid(column=0, row=row, sticky='w')


        row += 1
        l = Tkinter.Label(f, text='  Dimensions of new map:')
        l.grid(column=0, row=row, sticky='w')


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        
        l = Tkinter.Label(ff, text=' ', width=5)
        l.grid(column=0, row=0, sticky='w')

        self.newMapDimOption = Tkinter.StringVar()
        self.newMapDimOption.set ( 'box' )

        l = Tkinter.Radiobutton(ff, text="Same as map from which densities are extracted: ", variable=self.newMapDimOption, value = 'same')
        l.grid (column=1, row=0, sticky='w')

        self.oMapDim1 = Tkinter.StringVar(ff, d1)
        e = Tkinter.Entry(ff, width=5, textvariable=self.oMapDim1, state="disabled")
        e.grid(column=2, row=0, sticky='w', padx=5)

        self.oMapDim2 = Tkinter.StringVar(ff, d2)
        e = Tkinter.Entry(ff, width=5, textvariable=self.oMapDim2, state="disabled")
        e.grid(column=3, row=0, sticky='w', padx=5)

        self.oMapDim3 = Tkinter.StringVar(ff, d3)
        e = Tkinter.Entry(ff, width=5, textvariable=self.oMapDim3, state="disabled")
        e.grid(column=4, row=0, sticky='w', padx=5)


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        l = Tkinter.Label(ff, text=' ', width=5)
        l.grid(column=0, row=0, sticky='w')

        c = Tkinter.Radiobutton(ff, text="Cube around region, with border of", variable=self.newMapDimOption, value = 'cube')
        c.grid (column=1, row=0, sticky='w')

        self.borderWidth = Tkinter.StringVar(ff, st1)
        self.borderWidth.set ( "8" )
        e = Tkinter.Entry(ff, width=5, textvariable=self.borderWidth)
        e.grid(column=2, row=0, sticky='w', padx=5)

        c = Tkinter.Label(ff, text='voxels')
        c.grid (column=3, row=0, sticky='w')

        c = Tkinter.Radiobutton(ff, text="Box around region, with border of", variable=self.newMapDimOption, value = 'box')
        c.grid (column=1, row=1, sticky='w')

        e = Tkinter.Entry(ff, width=5, textvariable=self.borderWidth)
        e.grid(column=2, row=1, sticky='w', padx=5)

        c = Tkinter.Label(ff, text='voxels')
        c.grid (column=3, row=1, sticky='w')


        #row += 1
        #l = Tkinter.Label(f, text='  ')
        #l.grid(column=0, row=row, sticky='w')


        if 0 :
            row += 1
            l = Tkinter.Label(f, text='  Voxel size in new map:')
            l.grid(column=0, row=row, sticky='w')
    
    
            row += 1
            ff = Tkinter.Frame(f)
            ff.grid(column=0, row=row, sticky='w')
    
            self.newGridStep = Tkinter.StringVar()
            self.newGridStep.set ( 'same' )
    
            l = Tkinter.Label(ff, text=' ', width=5)
            l.grid(column=0, row=0, sticky='w')
    
            c = Tkinter.Radiobutton(ff, text="Same as map from which densities are extracted: ", variable=self.newGridStep, value = 'same')
            c.grid (column=1, row=0, sticky='w')
    
            self.oMapStep1 = Tkinter.StringVar(ff, st1)
            e = Tkinter.Entry(ff, width=5, textvariable=self.oMapStep1, state="disabled")
            e.grid(column=2, row=0, sticky='w', padx=5)
    
            self.oMapStep2 = Tkinter.StringVar(ff, st2)
            e = Tkinter.Entry(ff, width=5, textvariable=self.oMapStep2, state="disabled")
            e.grid(column=3, row=0, sticky='w', padx=5)
    
            self.oMapStep3 = Tkinter.StringVar(ff, st3)
            e = Tkinter.Entry(ff, width=5, textvariable=self.oMapStep3, state="disabled")
            e.grid(column=4, row=0, sticky='w', padx=5)
    
    
            row += 1
            ff = Tkinter.Frame(f)
            ff.grid(column=0, row=row, sticky='w')
    
            l = Tkinter.Label(ff, text=' ', width=5)
            l.grid(column=0, row=0, sticky='w')
    
            c = Tkinter.Radiobutton(ff, text="Custom: ", variable=self.newGridStep, value = 'custom')
            c.grid (column=1, row=0, sticky='w')
    
            self.newMapStep1 = Tkinter.StringVar(ff, ".5")
            e = Tkinter.Entry(ff, width=5, textvariable=self.newMapStep1)
            e.grid(column=2, row=0, sticky='w', padx=5)
    
            self.newMapStep2 = Tkinter.StringVar(ff, ".5")
            e = Tkinter.Entry(ff, width=5, textvariable=self.newMapStep2)
            e.grid(column=3, row=0, sticky='w', padx=5)
    
            self.newMapStep3 = Tkinter.StringVar(ff, ".5")
            e = Tkinter.Entry(ff, width=5, textvariable=self.newMapStep3)
            e.grid(column=4, row=0, sticky='w', padx=5)


        #row += 1
        #l = Tkinter.Label(f, text='  ')
        #l.grid(column=0, row=row, sticky='w')


        row += 1
        l = Tkinter.Label(f, text='  Which densities to extract:')
        l.grid(column=0, row=row, sticky='w')

        self.whichDensities = Tkinter.StringVar()
        self.whichDensities.set ( 'inside' )

        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        l = Tkinter.Label(ff, text=' ', width=5)
        l.grid(column=0, row=0, sticky='w')

        c = Tkinter.Radiobutton(ff, text="Only densities inside the selected region(s)", variable=self.whichDensities, value = 'inside')
        c.grid (column=1, row=0, sticky='w')

        c = Tkinter.Radiobutton(ff, text="All densities inside the bounds of the new map", variable=self.whichDensities, value = 'all')
        c.grid (column=1, row=1, sticky='w')


        row += 1
        l = Tkinter.Label(f, text='  What maps to create:')
        l.grid(column=0, row=row, sticky='w')

        self.whatMaps = Tkinter.StringVar()
        self.whatMaps.set ( 'combined' )

        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        l = Tkinter.Label(ff, text=' ', width=5)
        l.grid(column=0, row=0, sticky='w')

        c = Tkinter.Radiobutton(ff, text="A map spanning all selected regions", variable=self.whatMaps, value = 'combined')
        c.grid (column=1, row=0, sticky='w')

        c = Tkinter.Radiobutton(ff, text="A map for each selected region", variable=self.whatMaps, value = 'each')
        c.grid (column=1, row=1, sticky='w')


        #row += 1
        #l = Tkinter.Label(f, text='  ')
        #l.grid(column=0, row=row, sticky='w')


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        l = Tkinter.Label(ff, text=' ', width=1)
        l.grid(column=0, row=0, sticky='w')

        c = Hybrid.Checkbutton(ff, 'Add fall-off densities outside region boundary with witdth ~', False )
        c.button.grid (column=1, row=0, sticky='w')
        self.addDropOff = c.variable

        self.dropOffWidth = Tkinter.StringVar(ff, bwidth)
        e = Tkinter.Entry(ff, width=5, textvariable=self.dropOffWidth)
        e.grid(column=2, row=0, sticky='w', padx=5)

        l = Tkinter.Label(ff, text='Angstroms')
        l.grid(column=3, row=0, sticky='w')



        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        l = Tkinter.Label(ff, text=' ', width=1)
        l.grid(column=0, row=0, sticky='w')

        c = Hybrid.Checkbutton(ff, 'Low-pass result with gaussian filter width: ', False )
        c.button.grid (column=1, row=0, sticky='w')
        self.gaussLP = c.variable

        self.gaussLPWidth = Tkinter.StringVar(ff, "2.0")
        e = Tkinter.Entry(ff, width=5, textvariable=self.gaussLPWidth)
        e.grid(column=2, row=0, sticky='w', padx=5)

        #l = Tkinter.Label(ff, text='Angstroms')
        #l.grid(column=3, row=0, sticky='w')

        
        if 1 :

            row += 1
            ff = Tkinter.Frame(f)

            ff.grid(column=0, row=row, sticky='w')

            l = Tkinter.Label(ff, text=' ', width=1)
            l.grid(column=0, row=0, sticky='w')
    
            c = Hybrid.Checkbutton(ff, 'Smooth mask - gaussian width: ', False )
            c.button.grid (column=1, row=0, sticky='w')
            self.smoothMask = c.variable
    
            self.smoothMaskWidth = Tkinter.StringVar(ff, "10")
            e = Tkinter.Entry(ff, width=5, textvariable=self.smoothMaskWidth)
            e.grid(column=2, row=0, sticky='w', padx=5)
        
            l = Tkinter.Label(ff, text='Angstroms')
            l.grid(column=3, row=0, sticky='w')
        

        if dev_menus :
            row += 1
            ff = Tkinter.Frame(f)
            ff.grid(column=0, row=row, sticky='w')
            l = Tkinter.Label(ff, text=' ', width=1)
            l.grid(column=0, row=0, sticky='w')
    
            c = Hybrid.Checkbutton(ff, 'Add noise with mean: ', False )
            c.button.grid (column=1, row=0, sticky='w')
            self.addNoise = c.variable
    
            self.noiseMean = Tkinter.StringVar(ff, "0.0")
            e = Tkinter.Entry(ff, width=5, textvariable=self.noiseMean)
            e.grid(column=2, row=0, sticky='w', padx=5)
    
            l = Tkinter.Label(ff, text=', st.dev.:')
            l.grid(column=3, row=0, sticky='w')
    
            self.noiseStDev = Tkinter.StringVar(ff, "1.0")
            e = Tkinter.Entry(ff, width=5, textvariable=self.noiseStDev)
            e.grid(column=4, row=0, sticky='w', padx=5)


        
        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w')
        l = Tkinter.Label(ff, text=' ', width=1)
        l.grid(column=0, row=0, sticky='w')

        c = Hybrid.Checkbutton(ff, 'Resample result on grid of another map: ', False )
        c.button.grid (column=1, row=0, sticky='w')
        self.resample = c.variable

        self.resampleMap = Tkinter.StringVar(parent)
        #self.resampleMap.set ( "Not selected" )

        self.mbr  = Tkinter.Menubutton ( ff, textvariable=self.resampleMap, relief=Tkinter.RAISED )
        self.mbr.grid (column=2, row=0, sticky='we', padx=5)
        self.mbr.menu  =  Tkinter.Menu ( self.mbr, tearoff=0, postcommand=self.ResampleMapMenu )
        self.mbr["menu"]  =  self.mbr.menu



        if 1 :
            row += 1
            ff = Tkinter.Frame(f)
            ff.grid(column=0, row=row, sticky='w')
            l = Tkinter.Label(ff, text=' ', width=1)
            l.grid(column=0, row=0, sticky='w')
    
            c = Hybrid.Checkbutton(ff, 'Save map(s) with name: ', False )
            c.button.grid (column=1, row=0, sticky='w')
            self.saveMaps = c.variable
    
            if self.cur_dmap == None :
                base = ""
            else :
                base = os.path.splitext(self.cur_dmap.name)[0] + "_%s.mrc"
            self.saveMapsBaseName = Tkinter.StringVar(ff, base)
            e = Tkinter.Entry(ff, width=30, textvariable=self.saveMapsBaseName)
            e.grid(column=2, row=0, sticky='w', padx=5)


            row += 1
            ff = Tkinter.Frame(f)
            ff.grid(column=0, row=row, sticky='w')

            l = Tkinter.Label(ff, text='(%s in name -> region id, %d -> incrementing numbers starting with 1)', width=60, padx=20)
            l.grid(column=0, row=0, sticky='w')
    



        
        row += 1
        dummyFrame = Tkinter.Frame(parent, relief='groove', borderwidth=1)
        Tkinter.Frame(dummyFrame).pack()
        dummyFrame.grid(row=row,column=0,columnspan=7, pady=7, sticky='we')

        row = row + 1

        global msg
        msg = Tkinter.Label(parent, width = 60, anchor = 'w', justify = 'left', fg="red")
        msg.grid(column=0, row=row, sticky='ew')
        self.msg = msg
        row += 1



    def MapMenu ( self ) :

        self.mb.menu.delete ( 0, 'end' )        # Clear menu
        from VolumeViewer import Volume
        mlist = OML(modelTypes = [Volume])
        for m in mlist :
            self.mb.menu.add_radiobutton ( label=m.name, variable=self.dmap,
                                           command=lambda m=m: self.MapSelected(m) )

    def ResampleMapMenu ( self ) :

        self.mbr.menu.delete ( 0, 'end' )        # Clear menu
        from VolumeViewer import Volume
        mlist = OML(modelTypes = [Volume])
        for m in mlist :
            self.mbr.menu.add_radiobutton ( label=m.name, variable=self.resampleMap,
                                           command=lambda m=m: self.ResampleMapSelected(m) )


    def SetMapMenu (self, dmap):
        mname = dmap.name if dmap else ''
        self.dmap.set(mname)
        self.cur_dmap = dmap
        #print "Set map menu to ", dmap.name


    def ResampleMapSelected ( self, dmap ) :
        print "selected resample map: ", dmap.name
        self.resampleMapMod = dmap


    def MapSelected ( self, dmap ) :

        self.cur_dmap = dmap
        if dmap:
            #dmap.display = True

            self.oMapDim1.set ( "%d" % dmap.data.size[0] )
            self.oMapDim2.set ( "%d" % dmap.data.size[1] )
            self.oMapDim3.set ( "%d" % dmap.data.size[2] )

            #self.oMapStep1.set ( "%g" % dmap.data.step[0] )
            #self.oMapStep2.set ( "%g" % dmap.data.step[1] )
            #self.oMapStep3.set ( "%g" % dmap.data.step[2] )

        else :
            self.oMapDim1.set ( "" )
            self.oMapDim2.set ( "" )
            self.oMapDim3.set ( "" )
        
            #self.oMapStep1.set ( "" )
            #self.oMapStep2.set ( "" )
            #self.oMapStep3.set ( "" )






    def Extract ( self ) :

        fromMap = self.cur_dmap
        if fromMap == None :
            umsg ( "Please select a map to extract densities from" )
            return

        segMap = segmentation_map()
        segMod = current_segmentation ()

        if segMap == None :
            umsg ( "Please select a map in the Segment Map Dialog" )
            return

        if segMod == None :
            umsg ( "Please select a segmentation in the Segment Map Dialog" )
            return

        umsg ( "Extracting densities from " + fromMap.name + " based on selected regions in " + segMap.name )

        if segMap == fromMap :
            print "Same map!"

        else:
            print "Different maps!"
            

        self.rri = 0


        regs = segMod.selected_regions()
        if len(regs)==0 :
            umsg ( "no selected regions found" ); return
            
        if self.whatMaps.get() == "combined" :
            self.Extract2 ( fromMap, segMap, segMod, regs )
        
        else :
            for reg in regs :
                self.Extract2 ( fromMap, segMap, segMod, [reg] )
        
        print "done"
        

    

    def Extract2 ( self, fromMap, segMap, segMod, regs ) :

        reg_str = ""
        for r in regs : reg_str = reg_str + "_r%d" % r.rid
        print reg_str

        mdata = None
        ndata = None
        ndata_in = None


        if self.newMapDimOption.get() == "same" :

            print " - same dimensions as original map"

            if hasattr(self, 'smoothMask') and self.smoothMask.get () :
                print " - taking inner densities for soft mask"
                ndata_in = self.dataMaskedWithSelectedRegions ( segMod, segMap, fromMap, regs )
                if ndata_in == None : return

            if self.whichDensities.get() == "inside" :
                print " - want inner densities only"
                ndata = self.dataMaskedWithSelectedRegions ( segMod, segMap, fromMap, regs )
                if ndata == None : return

            else :
                print " - want all densities"
                # nothing to do, since they want the same-dimension map, with all the densities
                ndata = fromMap.data


        elif self.newMapDimOption.get() == "cube" or self.newMapDimOption.get() == "box" :

            print " - shrinking to region bounds "

            #dims = self.boundsOfSelectedRegions ( segMod, fromMap )
            #if dims == None : return
            #li,lj,lk, hi,hj,hk, n1,n2,n3 = dims

            mdata = self.dataMaskedWithSelectedRegions ( segMod, segMap, fromMap, regs )
            if mdata == None : return

            regsm = mdata.matrix()
            nze = numpy.nonzero ( regsm )
            # print nze

            li = numpy.min ( nze[0] )
            lj = numpy.min ( nze[1] )
            lk = numpy.min ( nze[2] )

            hi = numpy.max ( nze[0] )
            hj = numpy.max ( nze[1] )
            hk = numpy.max ( nze[2] )

            bound = int ( self.borderWidth.get() )
            li = li - bound; lj = lj - bound; lk = lk - bound
            hi = hi + bound; hj = hj + bound; hk = hk + bound

            n1 = hi - li + 1
            n2 = hj - lj + 1
            n3 = hk - lk + 1

            if n1 % 2 == 1 : n1 += 1
            if n2 % 2 == 1 : n2 += 1
            if n3 % 2 == 1 : n3 += 1

            if self.newMapDimOption.get() == "cube" :
                n = max ( n1, n2, n3 )
                li -= int ( (n-n1) / 2 )
                lj -= int ( (n-n2) / 2 )
                lk -= int ( (n-n3) / 2 )
                n1, n2, n3 = n, n, n


            print " - bounds of selected regions: %d %d %d --> %d %d %d --> %d %d %d" % ( li,lj,lk, hi,hj,hk, n1,n2,n3 )


            nmat = numpy.zeros ( (n1,n2,n3), numpy.float32 )
            nmat_in = numpy.zeros ( (n1,n2,n3), numpy.float32 )
            #dmat = dmap.full_matrix()

            dmat = fromMap.full_matrix()
            print "map grid dim: ", numpy.shape ( dmat )
            print "masked grid dim: ", numpy.shape ( regsm )
            print "new map grid dim: ", numpy.shape ( nmat )
            
            if hasattr(self, 'smoothMask') and self.smoothMask.get () :
                # copy the densities from insize the regions only
                for ii in range ( len(nze[0]) ) :
                    i,j,k = nze[0][ii], nze[1][ii], nze[2][ii]
                    #nmat[k-lk,j-lj,i-li] = regsm[k,j,i]
                    nmat_in[i-li,j-lj,k-lk] = regsm[i,j,k]

            if self.whichDensities.get() == "inside" :
                # copy the densities from insize the regions only
                for ii in range ( len(nze[0]) ) :
                    i,j,k = nze[0][ii], nze[1][ii], nze[2][ii]
                    #nmat[k-lk,j-lj,i-li] = regsm[k,j,i]
                    nmat[i-li,j-lj,k-lk] = regsm[i,j,k]

            else :
                # copy the densities inside the entire block
                for i in range ( li, hi ) :
                    for j in range ( lj, hj ) :
                        for k in range ( lk, hk ) :
                            #val = dmat[k,j,i]
                            #nmat[k-lk,j-lj,i-li] = val
                            try :
                                #nmat[k-lk,j-lj,i-li] = dmat[k,j,i]
                                nmat[i-li,j-lj,k-lk] = dmat[i,j,k]
                            except :
                                #print "iout (%d,%d,%d) (%d,%d,%d)" % (k, j, i, k-lk,j-lj,i-li)
                                pass
    
            O = fromMap.data.origin
            print "origin:", O
            nO = ( O[0] + float(lk) * fromMap.data.step[0],
                   O[1] + float(lj) * fromMap.data.step[1],
                   O[2] + float(li) * fromMap.data.step[2] )
            
            print "new origin:", nO

            ndata = VolumeData.Array_Grid_Data ( nmat, nO, fromMap.data.step, fromMap.data.cell_angles )
            ndata_in = VolumeData.Array_Grid_Data ( nmat_in, nO, fromMap.data.step, fromMap.data.cell_angles )




        if 0 and self.newGridStep.get() == "custom" :

            st1 = float ( self.newMapStep1.get() )
            st2 = float ( self.newMapStep2.get() )
            st3 = float ( self.newMapStep3.get() )

            print " - new step: %d %d %d", st1, st2, st3

            # make a temp map from which the density values will be interpolated

            try : nv = VolumeViewer.volume.add_data_set ( ndata, None )
            except : nv = VolumeViewer.volume.volume_from_grid_data ( ndata )
            nv.name = "TempMap"

            n1 = int ( numpy.ceil ( ndata.size[0] * ( ndata.step[0] / st1 ) ) )
            n2 = int ( numpy.ceil ( ndata.size[1] * ( ndata.step[1] / st2 ) ) )
            n3 = int ( numpy.ceil ( ndata.size[2] * ( ndata.step[2] / st3 ) ) )

            print " - new dimensions: %d %d %d" % (n1, n2, n3)

            # make a new matrix with the desired voxel size
            nmat = numpy.zeros ( (n3,n2,n1), numpy.float32 )
            ndata = VolumeData.Array_Grid_Data ( nmat, ndata.origin, (st1,st2,st3), ndata.cell_angles )

            npoints = VolumeData.grid_indices ( (n1,n2,n3), numpy.single)  # i,j,k indices
            _contour.affine_transform_vertices ( npoints, ndata.ijk_to_xyz_transform )

            dvals = nv.interpolated_values ( npoints, None, method = 'linear' )  # dmap.openState.xform )
            #dvals = numpy.where ( dvals > threshold, dvals, numpy.zeros_like(dvals) )
            #nze = numpy.nonzero ( dvals )

            nmat = dvals.reshape( (n3,n2,n1) )
            #f_mat = fmap.data.full_matrix()
            #f_mask = numpy.where ( f_mat > fmap.surface_levels[0], numpy.ones_like(f_mat), numpy.zeros_like(f_mat) )
            #df_mat = df_mat * f_mask

            ndata = VolumeData.Array_Grid_Data ( nmat, ndata.origin, (st1,st2,st3), ndata.cell_angles )

            nv.close ()


        try : nv = VolumeViewer.volume.add_data_set ( ndata, None )
        except : nv = VolumeViewer.volume.volume_from_grid_data ( ndata )
        nv.name = "Extracted Densities Map"

        if segMap == fromMap :
            nv.name = os.path.splitext(fromMap.name)[0] + reg_str + ".mrc"
        else :
            nv.name = os.path.splitext(fromMap.name)[0] + "_with_" + os.path.splitext(segMap.name)[0] + reg_str + ".mrc"




        if 0 and ndata_in != None :
            try : nvi = VolumeViewer.volume.add_data_set ( ndata_in, None )
            except : nvi = VolumeViewer.volume.volume_from_grid_data ( ndata_in )
            inm = nvi.full_matrix()
            in_mask = numpy.where ( inm > 0, numpy.ones_like(inm), numpy.zeros_like(inm) )
            chimera.openModels.close ( [nvi] )


        nvm = nv.full_matrix()
        f_mask = numpy.where ( nvm > 0, numpy.ones_like(nvm), numpy.zeros_like(nvm) )
        smMaskM = f_mask


        if self.addDropOff.get () :

            print "\n---adding dropoff---\n"

            if  ( 0 ) :
                # use gaussian
                s = nv.data.step
                width = max (nv.data.step) #float ( self.dropOffWidth.get() )
                width = numpy.sqrt ( s[0]*s[0] + s[1]*s[1] + s[2]*s[2] ) #float ( self.dropOffWidth.get() )
                from VolumeFilter import gaussian
                gvol = gaussian.gaussian_convolve (nv, width )
                gvol.name = nv.name # + "_g"
                gvm = gvol.full_matrix()

                ngvm = f_mask * gvm + nvm

            elif ( 1 ) :
                # heat equation - boundary conditions:
                # - fixed densities at region boundary
                # - 0 outside fall-off width

                s = nv.data.step # A/pixel
                diag_l = numpy.sqrt ( s[0]*s[0] + s[1]*s[1] + s[2]*s[2] ) # A/pixel
                desired_width = float ( self.dropOffWidth.get() ) # A
                num_it = desired_width / diag_l # how many iterations will reach the desired width
                print " - diagonal width: %.3f, desired width: %.3f, # iterations: %.3f" % (diag_l, desired_width, num_it)

                numit = int ( numpy.ceil ( num_it ) )
                gvm = nvm.copy();
                for i in range (numit) :
                    nv_1 = numpy.roll(gvm, 1, axis=0)
                    nv_2 = numpy.roll(gvm, -1, axis=0)
                    nv_3 = numpy.roll(gvm, 1, axis=1)
                    nv_4 = numpy.roll(gvm, -1, axis=1)
                    nv_5 = numpy.roll(gvm, 1, axis=2)
                    nv_6 = numpy.roll(gvm, -1, axis=2)
                    gvm = 1.0/6.0 * ( nv_1 + nv_2 + nv_3 + nv_4 + nv_5 + nv_6 )
                    gvm = f_mask * gvm + nvm
                    umsg ("Adding drop-off - iteration %d" % i)
                    
                if 0 :

                  # equilibrate while keep everything outside that has been reached so far (the fall-off width) 0
                  o_mask = numpy.where ( gvm > 0, numpy.ones_like(gvm), numpy.zeros_like(gvm) )
                  f_mask = f_mask * o_mask
  
                  for i in range (1000) :
                      nv_1 = numpy.roll(gvm, 1, axis=0)
                      nv_2 = numpy.roll(gvm, -1, axis=0)
                      nv_3 = numpy.roll(gvm, 1, axis=1)
                      nv_4 = numpy.roll(gvm, -1, axis=1)
                      nv_5 = numpy.roll(gvm, 1, axis=2)
                      nv_6 = numpy.roll(gvm, -1, axis=2)
                      gvm2 = f_mask * (1.0/6.0 * ( nv_1 + nv_2 + nv_3 + nv_4 + nv_5 + nv_6 )) + nvm
                      if numpy.sum ( numpy.abs(gvm2 - gvm) ) < 1e-3 :
                          print " - stopped equilibration after %d iterations" % i
                          gvm = gvm2
                          break
                      gvm = gvm2

                ngvm = f_mask * gvm + nvm


            else :
                from numpy import maximum as MMAX
                # fraction of max
                gvm = nvm.copy();
                for i in range (20) :
                    nv_1 = numpy.roll(gvm, 1, axis=0)
                    nv_2 = numpy.roll(gvm, -1, axis=0)
                    nv_3 = numpy.roll(gvm, 1, axis=1)
                    nv_4 = numpy.roll(gvm, -1, axis=1)
                    nv_5 = numpy.roll(gvm, 1, axis=2)
                    nv_6 = numpy.roll(gvm, -1, axis=2)
                    gvm = 4.0/6.0 * ( MMAX(nv_1, MMAX(nv_2, MMAX(nv_3, MMAX(nv_4, MMAX(nv_5,nv_6))))) )

                ngvm = f_mask * gvm + nvm

            ndata = VolumeData.Array_Grid_Data ( ngvm, nv.data.origin, nv.data.step, nv.data.cell_angles )
            try : nvg = VolumeViewer.volume.add_data_set ( ndata, None )
            except : nvg = VolumeViewer.volume.volume_from_grid_data ( ndata )
            nvg.name = nv.name

            chimera.openModels.close ( [nv] )
            nv = nvg



        if self.gaussLP.get() :
            lpw = float ( self.gaussLPWidth.get() )
            print " - gauss lp width:", lpw
            
            from VolumeFilter import gaussian
            gvol = gaussian.gaussian_convolve (nv, lpw )
            gvm = gvol.full_matrix()

            chimera.openModels.close ( [nv] )
            nv = gvol



        if hasattr(self, 'smoothMask') and self.smoothMask.get () :
            smw = float ( self.smoothMaskWidth.get() )
            print " smooth mask width: ", smw

            mm = ndata_in.full_matrix()
            maskM = numpy.where ( mm > 0, numpy.ones_like(mm), numpy.zeros_like(mm) )
            outM = numpy.where ( mm > 0, numpy.zeros_like(mm), numpy.ones_like(mm) )

            smMaskM = maskM.copy();

            if 0 :
                nmat = nv.data.full_matrix () * smMaskM;
                smdata = VolumeData.Array_Grid_Data ( nmat, ndata_in.origin, ndata_in.step, ndata_in.cell_angles )
                try : smv = VolumeViewer.volume.add_data_set ( smdata, None )
                except : smv = VolumeViewer.volume.volume_from_grid_data ( smdata )
                smv.name = "Hard masked"

            if 0 :
                smn = int ( self.smoothMaskWidth.get() )
                for i in range (smn) :
                    umsg ( "smooth mask iter %d/%d" % (i,smn) )
                    nv_1 = numpy.roll(smMaskM, 1, axis=0)
                    nv_2 = numpy.roll(smMaskM, -1, axis=0)
                    nv_3 = numpy.roll(smMaskM, 1, axis=1)
                    nv_4 = numpy.roll(smMaskM, -1, axis=1)
                    nv_5 = numpy.roll(smMaskM, 1, axis=2)
                    nv_6 = numpy.roll(smMaskM, -1, axis=2)
                    smMaskM = 1.0/6.0 * ( nv_1 + nv_2 + nv_3 + nv_4 + nv_5 + nv_6 )
                    smMaskM = outM * smMaskM + maskM
            else :
                print " - smoothing mask with Gaussian blur..."

                maskData0 = VolumeData.Array_Grid_Data ( smMaskM, ndata_in.origin, ndata_in.step, ndata_in.cell_angles )
                try : maskMap0 = VolumeViewer.volume.add_data_set ( maskData0, None )
                except : maskMap0 = VolumeViewer.volume.volume_from_grid_data ( maskData0 )

                from VolumeFilter import gaussian
                maskMapG = gaussian.gaussian_convolve (maskMap0, smw )
                smMaskM = maskMapG.full_matrix().copy()
                
                maskMap0.close() #.name = "mask map 0"
                maskMapG.close() #.name = "mask map G%.0f" % smw


            if 0 :
                print "saving smooth mask"
                mdata = VolumeData.Array_Grid_Data ( smMaskM, ndata_in.origin, ndata_in.step, ndata_in.cell_angles )
                try : mv = VolumeViewer.volume.add_data_set ( mdata, None )
                except : mv = VolumeViewer.volume.volume_from_grid_data ( mdata )
                mv.name = "Smooth Mask!"

            denMat = nv.data.full_matrix ()
            if 0 :
                print "making noise"
                from numpy.random import standard_normal as srand
                mean = float ( self.noiseMean.get() )
                stdev = float ( self.noiseStDev.get() )
                s=nv.data.size
                denMat = srand ( (s[2],s[1],s[0]) ) * stdev - (numpy.ones_like(nvm) * mean)

            if 0 :
                print "MAKING MASK ALL 1s"
                smMaskM = numpy.where ( smMaskM > 0, numpy.ones_like(smMaskM), numpy.zeros_like(smMaskM) )

            nmat = denMat * smMaskM;
            smdata = VolumeData.Array_Grid_Data ( nmat, ndata_in.origin, ndata_in.step, ndata_in.cell_angles )
            try : smv = VolumeViewer.volume.add_data_set ( smdata, None )
            except : smv = VolumeViewer.volume.volume_from_grid_data ( smdata )
            
            nv.close ()
            nv = smv

        nv.openState.xform = fromMap.openState.xform



        if hasattr(self, 'addNoise') and self.addNoise.get () :

            mean = float ( self.noiseMean.get() )
            stdev = float ( self.noiseStDev.get() )

            print "\n---adding noise mean:",mean, " stdev:", stdev, "---\n"

            nvm = nv.full_matrix()
            #f_mask = numpy.where ( nvm > 0, numpy.zeros_like(nvm), numpy.ones_like(nvm) )

            from numpy.random import standard_normal as srand
            s=nv.data.size

            noisem = srand ( (s[2],s[1],s[0]) ) * stdev - (numpy.ones_like(nvm) * mean)
            ngvm = noisem + nvm
            #ngvm = smMaskM * noisem
            #ngvm = noisem

            ndata = VolumeData.Array_Grid_Data ( ngvm, nv.data.origin, nv.data.step, nv.data.cell_angles, name=nv.name + "_noise" )
            try : nvg = VolumeViewer.volume.add_data_set ( ndata, None )
            except : nvg = VolumeViewer.volume.volume_from_grid_data ( ndata )
            nvg.name = nv.name + "_noise"

            chimera.openModels.close ( [nv] )
            nv = nvg



        if self.resample.get () :
            print "Resampling to", self.resampleMapMod.name
            nvr = place_map_resample ( nv, self.resampleMapMod )
            nv.close ()
            nv = nvr
            nv.openState.xform = self.resampleMapMod.openState.xform


        if "%s" in self.saveMapsBaseName.get() :
            nv.name = self.saveMapsBaseName.get() % reg_str
        elif "%d" in self.saveMapsBaseName.get() :
            rri = 0
            try : 
                rri = self.rri + 1
            except :
                rri = 1
            self.rri = rri
            nv.name = self.saveMapsBaseName.get() % rri
        else :
            nv.name = self.saveMapsBaseName.get()


        ndata = VolumeData.Array_Grid_Data ( nv.full_matrix(), nv.data.origin, nv.data.step, nv.data.cell_angles, name=nv.name )
        try : nvg = VolumeViewer.volume.add_data_set ( ndata, None )
        except : nvg = VolumeViewer.volume.volume_from_grid_data ( ndata )
        nvg.name = nv.name
        nv.close()
        nv = nvg

        umsg ( "Done - created " + nvg.name )
        
        if self.saveMaps.get() :
            mdir, mfile = os.path.split(fromMap.data.path)
            dpath = mdir + "/" + nv.name
            print "Saving extracted map to", dpath
            nv.write_file ( dpath, "mrc" )




    def EQ0 ( self ) :

        print " - equalize 0 "

        dmap = self.cur_dmap
        if dmap == None :
            umsg ( "Please select a map to extract densities from" )
            return

        print " - map: ", dmap.name


        mat = dmap.data.full_matrix()
        data = dmap.data

        print " - smoothing..."
        from VolumeFilter import gaussian
        gvol = gaussian.gaussian_convolve (dmap, 40.0 )
        matg = gvol.full_matrix()

        matg = numpy.where ( mat > 0.0001, matg, numpy.ones_like(matg) )

        mat2 = mat / matg

        data2 = VolumeData.Array_Grid_Data ( mat2, data.origin, data.step, data.cell_angles )
        try : dmap2 = VolumeViewer.volume.add_data_set ( data2, None )
        except : dmap2 = VolumeViewer.volume.volume_from_grid_data ( data2 )
        
        dmap2.openState.xform = dmap.openState.xform
        dmap2.name = dmap.name + "__EQ"


    def EQ ( self ) :

        print " - equalize "

        fromMap = self.cur_dmap
        if fromMap == None :
            umsg ( "Please select a map to extract densities from" )
            return

        segMap = segmentation_map()
        segMod = current_segmentation ()

        if segMap == None :
            umsg ( "Please select a map in the Segment Map Dialog" )
            return

        if segMod == None :
            umsg ( "Please select a segmentation in the Segment Map Dialog" )
            return

        umsg ( "Equalizing " + fromMap.name + " regions in " + segMap.name )


        regs = segMod.selected_regions()
        if len(regs)==0 :
            umsg ( "no selected regions found" ); return
        
        print " %d regions" % len(regs)
        

        sumMat = None
        

        print "Collecting stats |%d regions|" % len(regs)        
        maxD = 0.0
        minD = 10e7
        for ri, reg in enumerate ( regs ) :
            rdata = self.dataMaskedWithSelectedRegions ( segMod, segMap, fromMap, [reg] )
            rmat = rdata.full_matrix ()
            weights = rmat.ravel()
            smin = numpy.min (weights)
            sdev = numpy.std (weights)
            savg = numpy.average(weights)
            smax = numpy.max (weights)
            if smax > maxD :
                maxD = smax
            if smin < minD :
                minD = smin
            
            print " %d" % (ri+1),

        print ""
        print "Max: ", maxD
        print "Min: ", minD
        

        print "Equalizing and adding regions maps... |%d|" % len(regs)        
        for ri, reg in enumerate (regs) :
        
            rdata = self.dataMaskedWithSelectedRegions ( segMod, segMap, fromMap, [reg] )
            
            rmat = rdata.full_matrix ()
            
            weights = rmat.ravel()
            smin = numpy.min (weights)
            sdev = numpy.std (weights)
            savg = numpy.average(weights)
            smax = numpy.max (weights)            
            
            rmat = rmat - ( numpy.ones_like(rmat) * smin )

            #df_mat = numpy.where ( df_mat > thr, df_mat, numpy.zeros_like(df_mat) )

            rmat = rmat * ( maxD / smax )
            rmat = rmat + ( numpy.ones_like(rmat) * minD )

            if 0 :
                rmf = rmat.flatten()
                imhist,bins = numpy.histogram ( rmf, 20, normed=True )
                cdf = imhist.cumsum() #cumulative distribution function
                cdf = 10.0 * cdf / cdf[-1] #normalize
                #use linear interpolation of cdf to find new pixel values
                rmat = numpy.interp ( rmf, bins[:-1], cdf )
                rmat = rmat.reshape(segMap.data.full_matrix().shape)


            if sumMat == None :
                sumMat = rmat
            else :
                sumMat = sumMat + rmat
            

            print " %d" % (ri+1),
            #print " - r %d - (%.6f,%.6f) |%.6f| +/- %.6f" % ( reg.rid, smin, smax, savg, sdev )

        print ""

        gvm = sumMat.copy();
        for i in range ( 5 ) :
            nv_1 = numpy.roll(gvm, 1, axis=0)
            nv_2 = numpy.roll(gvm, -1, axis=0)
            nv_3 = numpy.roll(gvm, 1, axis=1)
            nv_4 = numpy.roll(gvm, -1, axis=1)
            nv_5 = numpy.roll(gvm, 1, axis=2)
            nv_6 = numpy.roll(gvm, -1, axis=2)
            gvm = 1.0/6.0 * ( nv_1 + nv_2 + nv_3 + nv_4 + nv_5 + nv_6 )
            #gvm =  gvm + nvm


        df_data = VolumeData.Array_Grid_Data ( gvm, segMap.data.origin, segMap.data.step, segMap.data.cell_angles )
        try : df_v = VolumeViewer.volume.add_data_set ( df_data, None )
        except : df_v = VolumeViewer.volume.volume_from_grid_data ( df_data )
        df_v.name = "Equalized"
        df_v.openState.xform = segMap.openState.xform


        if "%s" in self.saveMapsBaseName.get() :
            df_v.name = self.saveMapsBaseName.get() % "_EQ"
        elif "%d" in self.saveMapsBaseName.get() :
            df_v.name = self.saveMapsBaseName.get() % 0
        else :
            df_v.name = self.saveMapsBaseName.get()
        
        umsg ( "Done - created " + df_v.name )
        
        if self.saveMaps.get() :
            mdir, mfile = os.path.split(fromMap.data.path)
            dpath = mdir + "/" + df_v.name
            print "Saving extracted map to", dpath
            df_v.write_file ( dpath, "mrc" )




    def dataMaskedWithSelectedRegions ( self, segModel, segMap, fromMap, regs ) :

        points = regs[0].points().astype ( numpy.float32 )
        for r in regs[1:] :
            npoints = r.points().astype ( numpy.float32 )
            points = numpy.concatenate ( [points, npoints], axis=0 )

        _contour.affine_transform_vertices ( points, segMap.data.ijk_to_xyz_transform )
        _contour.affine_transform_vertices ( points, Matrix.xform_matrix( segModel.openState.xform ) )
        _contour.affine_transform_vertices ( points, Matrix.xform_matrix( fromMap.openState.xform.inverse() ) )

        d = min ( segMap.data.step ) * 0.99
        ndata = VolumeData.zone_masked_grid_data ( fromMap.data, points, d )

        return ndata

          

    def MaskAnotherMapWRegions ( self ) :

        if len(self.dmap.get()) == 0 : umsg ("Please select a map from which density will be taken" ); return
        dmap = self.SegmentationMap()
        if dmap == None : umsg ( "%s is not open" % self.dmap.get() ); return

        smod = self.CurrentSegmentation()
        if smod is None : return

        regs = smod.selected_regions()
        if len(regs)==0 :
            umsg ( "Please select one ore more regions" )
            return

        points = regs[0].points().astype ( numpy.float32 )
        for r in regs[1:] :
            npoints = r.points().astype ( numpy.float32 )
            points = numpy.concatenate ( [points, npoints], axis=0 )

        _contour.affine_transform_vertices ( points, dmap.data.ijk_to_xyz_transform )
        _contour.affine_transform_vertices ( points, Matrix.xform_matrix( smod.openState.xform ) )
        _contour.affine_transform_vertices ( points, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )

        sg = VolumeData.zone_masked_grid_data ( dmap.data, points, dmap.data.step[0] )

        try : gv = VolumeViewer.volume.add_data_set ( sg, None )
        except : gv = VolumeViewer.volume.volume_from_grid_data ( sg )
        gv.openState.xform = dmap.openState.xform
        #chimera.openModels.add ( [gv] )
        gv.name = "Masked"





    def MaskMapWRegionsCube ( self ) :

        # thsi is useful for input to EMAN fitting procedures which
        # requires a cube map

        if len(self.dmap.get()) == 0 : umsg ("Please select a map from which density will be taken" ); return
        dmap = self.SegmentationMap()
        if dmap == None : umsg ( "%s is not open" % self.dmap.get() ); return

        smod = self.CurrentSegmentation()
        if smod is None : return

        regs = smod.selected_regions()
        if len(regs)==0 :
            umsg ( "Please select one ore more regions" )
            return

        if 0 :
            points = regs[0].points().astype ( numpy.float32 )
            for r in regs[1:] :
                npoints = r.points().astype ( numpy.float32 )
                points = numpy.concatenate ( [points, npoints], axis=0 )

        for rri, reg in enumerate ( regs ) :

            print " ---- Region %d/%d ---- " % (rri+1, len(regs))

            points = reg.points().astype ( numpy.float32 )

            _contour.affine_transform_vertices ( points, dmap.data.ijk_to_xyz_transform )
            _contour.affine_transform_vertices ( points, Matrix.xform_matrix( smod.openState.xform ) )
            _contour.affine_transform_vertices ( points, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )

            sg = VolumeData.zone_masked_grid_data ( dmap.data, points, dmap.data.step[0] )
            regsm = sg.matrix()

            nze = numpy.nonzero ( regsm )

            # print nze

            li = numpy.min ( nze[0] )
            lj = numpy.min ( nze[1] )
            lk = numpy.min ( nze[2] )

            hi = numpy.max ( nze[0] )
            hj = numpy.max ( nze[1] )
            hk = numpy.max ( nze[2] )

            ci = int ( numpy.ceil ( (hi + li) / 2 ) )
            cj = int ( numpy.ceil ( (hj + lj) / 2 ) )
            ck = int ( numpy.ceil ( (hk + lk) / 2 ) )

            n1 = hi - li + 1
            n2 = hj - lj + 1
            n3 = hk - lk + 1

            n = 120 # max ( n1, n2, n3 ) + 4
            n2 = int ( numpy.ceil ( n / 2 ) )

            li = ci - n2; lj = cj - n2; lk = ck - n2
            hi = ci + n2; hj = cj + n2; hk = ck + n2

            print "Bounds - %d %d %d --> %d %d %d --> %d %d %d (%d)" % ( li, lj, lk, hi, hj, hk, n1, n2, n3, n )

            umsg ( "Saving %d regions to mrc file..." % len(regs) )

            nmat = numpy.zeros ( (n,n,n), numpy.float32 )
            #dmat = dmap.full_matrix()

            print "map grid dim: ", numpy.shape ( dmap.full_matrix() )
            print "masked grid dim: ", numpy.shape ( regsm )
            print "new map grid dim: ", numpy.shape ( nmat )


            #regs_name = ""
            for ii in range ( len(nze[0]) ) :
                i,j,k = nze[0][ii], nze[1][ii], nze[2][ii]
                mapVal = regsm[i,j,k]
                nmat[i-li,j-lj,k-lk] = mapVal

            O = dmap.data.origin
            print "origin:", O
            if 1 :
                nO = ( O[0] + float(lk) * dmap.data.step[0],
                       O[1] + float(lj) * dmap.data.step[1],
                       O[2] + float(li) * dmap.data.step[2] )
                print "new origin:", nO
            else :
                nO = ( -float(n2) * dmap.data.step[0],
                       -float(n2) * dmap.data.step[1],
                       -float(n2) * dmap.data.step[2] )
                print "new origin:", nO
            

            ndata = VolumeData.Array_Grid_Data ( nmat, nO, dmap.data.step, dmap.data.cell_angles )
            try : nv = VolumeViewer.volume.add_data_set ( ndata, None )
            except : nv = VolumeViewer.volume.volume_from_grid_data ( ndata )

            suff = "_CubeRid%d.mrc" % reg.rid

            nv.name = os.path.splitext (dmap.name) [0] + suff
            nv.openState.xform = dmap.openState.xform

            path = os.path.splitext (dmap.data.path) [0] + suff
            nv.write_file ( path, "mrc" )


    def ExtractMapWRegionsCube ( self ) :

        if len(self.dmap.get()) == 0 : umsg ("Please select a map from which density will be taken" ); return
        dmap = self.SegmentationMap()
        if dmap == None : umsg ( "%s is not open" % self.dmap.get() ); return

        smod = self.CurrentSegmentation()
        if smod is None : return

        regs = smod.selected_regions()
        if len(regs)==0 :
            umsg ( "Please select one ore more regions" )
            return

        if 0 :
            points = regs[0].points().astype ( numpy.float32 )
            for r in regs[1:] :
                npoints = r.points().astype ( numpy.float32 )
                points = numpy.concatenate ( [points, npoints], axis=0 )

        for rri, reg in enumerate ( regs ) :

            (li,lj,lk), (hi,hj,hk) = regions.region_bounds( [reg] )

            ci = int ( numpy.ceil ( (hi + li) / 2 ) )
            cj = int ( numpy.ceil ( (hj + lj) / 2 ) )
            ck = int ( numpy.ceil ( (hk + lk) / 2 ) )

            n1 = hi - li + 1
            n2 = hj - lj + 1
            n3 = hk - lk + 1

            n = 62 # max ( n1, n2, n3 ) + 4
            n2 = int ( numpy.ceil ( n / 2 ) )

            li = ci - n2; lj = cj - n2; lk = ck - n2
            hi = ci + n2; hj = cj + n2; hk = ck + n2

            #bound = 2
            #li = li - bound; lj = lj - bound; lk = lk - bound
            #hi = hi + bound; hj = hj + bound; hk = hk + bound

            print "Bounds - %d %d %d --> %d %d %d --> %d %d %d, %d" % ( li, lj, lk, hi, hj, hk, n1,n2,n3, n )

            umsg ( "Saving %d regions to mrc file..." % len(regs) )

            #nmat = numpy.zeros ( (n3,n2,n1), numpy.float32 )
            nmat = numpy.zeros ( (n,n,n), numpy.float32 )
            dmat = dmap.full_matrix()

            #regs_name = ""
            for i in range ( li, hi ) :
                for j in range ( lj, hj ) :
                    for k in range ( lk, hk ) :
                        try :
                            nmat[k-lk,j-lj,i-li] = dmat[k,j,i]
                        except :
                            pass

            O = dmap.data.origin
            print "origin:", O
            nO = O
            if 1 :
                nO = ( O[0] + float(li) * dmap.data.step[0],
                       O[1] + float(lj) * dmap.data.step[1],
                       O[2] + float(lk) * dmap.data.step[2] )
            else :            
                nO = ( -float(n2) * dmap.data.step[0],
                       -float(n2) * dmap.data.step[1],
                       -float(n2) * dmap.data.step[2] )
                print "new origin:", nO

            ndata = VolumeData.Array_Grid_Data ( nmat, nO, dmap.data.step, dmap.data.cell_angles )
            try : nv = VolumeViewer.volume.add_data_set ( ndata, None )
            except : nv = VolumeViewer.volume.volume_from_grid_data ( ndata )

            suff = "_MC_%d_%d_%d_%d.mrc" % (li, lj, lk, n)

            nv.name = os.path.splitext (dmap.name) [0] + suff
            nv.openState.xform = dmap.openState.xform

            path = os.path.splitext (dmap.data.path) [0] + suff
            nv.write_file ( path, "mrc" )
        



def show_extract_region_dialog ( closeOld = True ):

	from chimera import dialogs
	
	d = dialogs.find ( "extract region", create=False )
	if d :
		if closeOld :
			d.toplevel_widget.update_idletasks ()
			d.Close()
			d.toplevel_widget.update_idletasks ()
		else :
			# is there a way to bring it to front?
			return d
	
	dialogs.register (Extract_Region_Dialog.name, Extract_Region_Dialog, replace = True)
	
	d = dialogs.find ( "extract region", create=True )
	d.toplevel_widget.update_idletasks ()
	d.enter()
	
	return d


def dialog ():

	from chimera import dialogs
	return dialogs.find ( "extract region", create=False )
	
	

def MapSS ( dmap, d ) :

    ndata = dmap.data
    
    st1 = ndata.step[0] / d
    st2 = ndata.step[1] / d
    st3 = ndata.step[2] / d

    n1 = int ( numpy.ceil ( ndata.size[0] * ( ndata.step[0] / st1 ) ) )
    n2 = int ( numpy.ceil ( ndata.size[1] * ( ndata.step[1] / st2 ) ) )
    n3 = int ( numpy.ceil ( ndata.size[2] * ( ndata.step[2] / st3 ) ) )

    print " - new dimensions: %d %d %d" % (n1, n2, n3)

    # make a new matrix with the desired voxel size
    nmat = numpy.zeros ( (n3,n2,n1), numpy.float32 )
    ndata = VolumeData.Array_Grid_Data ( nmat, ndata.origin, (st1,st2,st3), ndata.cell_angles )

    npoints = VolumeData.grid_indices ( (n1,n2,n3), numpy.single)  # i,j,k indices
    _contour.affine_transform_vertices ( npoints, ndata.ijk_to_xyz_transform )

    dvals = dmap.interpolated_values ( npoints, None, method = 'linear' )  # dmap.openState.xform )
    #dvals = numpy.where ( dvals > threshold, dvals, numpy.zeros_like(dvals) )
    #nze = numpy.nonzero ( dvals )

    nmat = dvals.reshape( (n3,n2,n1) )
    #f_mat = fmap.data.full_matrix()
    #f_mask = numpy.where ( f_mat > fmap.surface_levels[0], numpy.ones_like(f_mat), numpy.zeros_like(f_mat) )
    #df_mat = df_mat * f_mask

    ndata = VolumeData.Array_Grid_Data ( nmat, ndata.origin, (st1,st2,st3), ndata.cell_angles )

    return ndata



def place_map_resample_R ( fmap, dmap, useThr = True ) :

    # get bounds of points above threshold
    fpoints = VolumeData.grid_indices (fmap.data.size, numpy.single)  # i,j,k indices
    _contour.affine_transform_vertices ( fpoints, fmap.data.ijk_to_xyz_transform )
    
    bound = 0

    if useThr :
        mat = fmap.data.full_matrix ()
        fpoint_weights = numpy.ravel(mat).astype(numpy.single)
        threshold = 0.01 # fmap.surface_levels[0]
        ge = numpy.greater_equal(fpoint_weights, threshold)
        fpoints = numpy.compress(ge, fpoints, 0)
        fpoint_weights = numpy.compress(ge, fpoint_weights)
        nz = numpy.nonzero( fpoint_weights )[0]
        print " - %d above %f in %s" % (len(nz), threshold, fmap.name)
        #print "points: ", fpoints
        #print "weights: ", fpoint_weights
        bound = 2

    _contour.affine_transform_vertices ( fpoints, Matrix.xform_matrix( fmap.openState.xform ) )
    _contour.affine_transform_vertices ( fpoints, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
    _contour.affine_transform_vertices ( fpoints, dmap.data.xyz_to_ijk_transform )
    #print "points in %s ref:" % dmap.name, fpoints

    li,lj,lk = numpy.min ( fpoints, axis=0 ) - (bound, bound, bound)
    hi,hj,hk = numpy.max ( fpoints, axis=0 ) + (bound, bound, bound)

    n1 = hi - li + 1
    n2 = hj - lj + 1
    n3 = hk - lk + 1

    print " - bounds - %d %d %d --> %d %d %d --> %d %d %d" % ( li, lj, lk, hi, hj, hk, n1,n2,n3 )

    #nmat = numpy.zeros ( (n1,n2,n3), numpy.float32 )
    #dmat = dmap.full_matrix()

    nn1 = int ( round (dmap.data.step[0] * float(n1) / fmap.data.step[0]) )
    nn2 = int ( round (dmap.data.step[1] * float(n2) / fmap.data.step[1]) )
    nn3 = int ( round (dmap.data.step[2] * float(n3) / fmap.data.step[2]) )

    O = dmap.data.origin
    print " - %s origin:" % dmap.name, O
    nO = ( O[0] + float(li) * dmap.data.step[0],
           O[1] + float(lj) * dmap.data.step[1],
           O[2] + float(lk) * dmap.data.step[2] )
    
    print " - new map origin:", nO

    nmat = numpy.zeros ( (nn1,nn2,nn3), numpy.float32 )
    ndata = VolumeData.Array_Grid_Data ( nmat, nO, fmap.data.step, dmap.data.cell_angles )

    print " - fmap grid dim: ", numpy.shape ( fmap.full_matrix() )
    print " - new map grid dim: ", numpy.shape ( nmat )

    npoints = VolumeData.grid_indices ( (nn1, nn2, nn3), numpy.single)  # i,j,k indices
    _contour.affine_transform_vertices ( npoints, ndata.ijk_to_xyz_transform )

    dvals = fmap.interpolated_values ( npoints, dmap.openState.xform )
    #dvals = numpy.where ( dvals > threshold, dvals, numpy.zeros_like(dvals) )
    #nze = numpy.nonzero ( dvals )

    nmat = dvals.reshape( (nn3,nn2,nn1) )
    #f_mat = fmap.data.full_matrix()
    #f_mask = numpy.where ( f_mat > fmap.surface_levels[0], numpy.ones_like(f_mat), numpy.zeros_like(f_mat) )
    #df_mat = df_mat * f_mask

    ndata = VolumeData.Array_Grid_Data ( nmat, nO, fmap.data.step, dmap.data.cell_angles )
    try : nv = VolumeViewer.volume.add_data_set ( ndata, None )
    except : nv = VolumeViewer.volume.volume_from_grid_data ( ndata )

    return nv




def place_map_resample ( densitiesFromMap, toGridOfMap, mask = False ) :

    #print "Taking densities from %s with grid of %s" % ( densitiesFromMap.name, toGridOfMap.name )

    #mmc = fmap.writable_copy ( require_copy = True )
    #mmc.name = rname
    #print " - cloned", fmap.name

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
    
    
    ndata = VolumeData.Array_Grid_Data ( df_mat, fmap.data.origin, fmap.data.step, fmap.data.cell_angles )
    try : nv = VolumeViewer.volume.add_data_set ( ndata, None )
    except : nv = VolumeViewer.volume.volume_from_grid_data ( ndata )

    return nv



