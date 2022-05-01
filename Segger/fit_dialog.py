
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
import numpy

from VolumeData import grid_indices, zone_masked_grid_data, interpolate_volume_data
from _multiscale import get_atom_coordinates
from _contour import affine_transform_vertices as transform_vertices
from Matrix import xform_matrix, multiply_matrices, chimera_xform, identity_matrix, invert_matrix, shift_and_angle
from VolumeViewer import volume_from_grid_data
from VolumeViewer.volume import Rendering_Options
from time import clock
from random import random as rand
import FitMap
import VolumeViewer
import Segger.quaternion
import Matrix
import VolumeData

from axes import prAxes
from regions import mask_volume, regions_radius

from segment_dialog import current_segmentation, segmentation_map
from Segger import dev_menus, timing, seggerVersion

OML = chimera.openModels.list

SAF_DVOL = 0.75
SAF_DBRAD = 0.3
SAF_LS_DEPTH = 4
SAF_LS_NGROUPS = 1000

REG_OPACITY = 0.45
MAX_NUM_GROUPS = 1000



def umsg ( txt ) :
    print txt
    status ( txt )


def status ( txt ) :
    txt = txt.rstrip('\n')
    msg.configure(text = txt)
    msg.update_idletasks()


#from fit_devel import Fit_Devel

class Fit_Segments_Dialog ( chimera.baseDialog.ModelessDialog ):

    title = "SegFit"
    name = "fit segments"
    #buttons = ( 'SMS', 'Scores', 'Fit', 'Options', "Close")
    #buttons = ( 'Place', 'Fit', 'Options', "Close")
    buttons = ( 'Clear', 'Fit', 'Stop', 'Options', "Close")
    help = 'https://github.com/gregdp/segger'

    def fillInUI(self, parent):

        import Tkinter
        from CGLtk import Hybrid

        tw = parent.winfo_toplevel()
        self.toplevel_widget = tw
        tw.withdraw()

        parent.columnconfigure(0, weight = 1)

        row = 1

        menubar = Tkinter.Menu(parent, type = 'menubar', tearoff = False)
        tw.config(menu = menubar)

        self.UseAllMods = Tkinter.IntVar()
        self.UseAllMods.set ( 0 )

        fit_menu_entries = [
            ('Delete selected fits from list', self.delete_fit_cb),
            ('Delete ALL fits from list', self.delete_all_fit_cb),
            'separator',
            ('Place molecule copies', self.place_copies_cb),
            ('Place map copies', self.place_map_copies_cb),
            #('Cube map', self.extractCubeMap),
            ('Close placed copies', self.close_copies_cb),

            'separator',
            ("Save chosen fit molecules", self.SaveStrucFit),

            'separator',
            ('Place selected map relative to segmented map', self.save_map_resample),

            'separator']

        if dev_menus :
            fit_menu_entries = fit_menu_entries + [
                ('Group regions by SS in visible (Molecule) models', self.GroupRegionsBySS),
                ('Mask map with selection', self.MaskWithSel)
            ]

        fit_menu_entries = fit_menu_entries + [
            ('Group regions by visible (Molecule) models', self.GroupRegionsByMols),
            ('Group regions by chains in visible (Molecule) models', self.GroupRegionsByChains),

            'separator',
            ("Show molecule axes", self.StrucShowAxes),
            ("Hide molecule axes", self.StrucHideAxes),
            ("Show overlapping regions", self.ShowOverlappingRegions),

            'separator',
            ("Export fit scores", self.ExportFitScores),
            ("Plot fit scores", self.PlotFitScores),
            ( "Score Visible", self.VisiScores )
            ]

        if dev_menus :
            fit_menu_entries = fit_menu_entries + [
                'separator',
                ('Difference map', self.DifferenceMap),
                ('Intersection map', self.IntersectionMap),
                ('Shape match', self.ShapeMatch),

                'separator',
                ('Fit all visible maps to selected', self.FitAllVisMaps),
                ('Make average map of visible fitted maps', self.AvgFMaps2),
                ('Make difference map of visible fitted maps', self.DifFMaps2),
                ('Take fitted map densities into segmented map', self.TakeFMap_with_DMap0),
                ('Take fitted map densities into segmented map + noise', self.TakeFMap_with_DMapN),
                ('Take fitted map densities into segmented map (Shrink/grow)', self.TakeFMap_with_DMap),
                ('Save visible fitted maps in segmented map grid', self.TakeFMapsVis),
                ('Take segmented map densities into fitted map', self.TakeDMap_with_FMap),
                'separator',
                ('Group regions by chains in visible molecules', self.GroupRegionsByChains),
                ('Group regions by visible molecules', self.GroupRegionsByMols),
                ('Group regions by selected fitted molecules', self.GroupRegionsByFittedMols),
                ('Group regions by visible maps', self.GroupRegionsByVisiMaps),
                'separator',
                #('0 map with selection', self.ZeroMapBySel),
                #('0 map with visible molecules', self.ZeroMapByMols),
                ('0 map with selected fitted molecules', self.ZeroMapFittedMols),
                ('0 map with visible molecules', self.ZeroMapVisMols),
                ('Values in map', self.ValuesInMap),
                ('Mask with selected map/model', self.MaskMapWithSel)
                 ]


        fmenu = Hybrid.cascade_menu(menubar, 'Fit', fit_menu_entries)
        #self.add_devel_menus(fmenu)

        from chimera.tkgui import aquaMenuBar
        aquaMenuBar(menubar, parent, row = 0, columnspan=3)


        f = Tkinter.Frame(parent)
        f.grid(column=0, row=row, sticky='ew')
        row += 1

        l = Tkinter.Label(f, text='Structure or Map to fit')
        l.grid(column=0, row=0, sticky='w')

        self.struc = Tkinter.StringVar(parent)
        self.strucMB  = Tkinter.Menubutton ( f, textvariable=self.struc, relief=Tkinter.RAISED )
        self.strucMB.grid (column=1, row=0, sticky='we', padx=5)
        self.strucMB.menu  =  Tkinter.Menu ( self.strucMB, tearoff=0, postcommand=self.StrucMenu )
        self.strucMB["menu"]  =  self.strucMB.menu


        h = '%10s %10s %10s %10s %10s %10s %10s %10s' % ('Corr.', 'At. Incl.', 'BB Incl.', 'Clashes', 'Dens. Occ.', 'Molecule', 'Map', 'Region')
        fl = Hybrid.Scrollable_List(parent, h, 8, self.fit_selection_cb)
        self.fit_listbox = fl.listbox
        self.list_fits = []
        fl.frame.grid(row = row, column = 0, sticky = 'news')
        parent.rowconfigure(row, weight = 1)
        row += 1
        self.fit_listbox.bind('<KeyPress-Delete>', self.delete_fit_cb)


        op = Hybrid.Popup_Panel(parent)
        opf = op.frame
        opf.grid(row = row, column = 0, sticky = 'news')
        opf.grid_remove()
        opf.columnconfigure(0, weight=1)
        self.optionsPanel = op.panel_shown_variable
        row += 1
        orow = 0

        cb = op.make_close_button(opf)
        cb.grid(row = orow, column = 0, sticky = 'e')

        l = Tkinter.Label(opf, text='Fitting Options', font = 'TkCaptionFont')
        l.grid(column=0, row=orow, sticky='w', pady=5)
        orow += 1

        fopt = Tkinter.Frame(opf)
        fopt.grid(column=0, row=orow, sticky='ew', padx=10)
        orow += 1
        forow = 0


        if 0 :
            oft = Hybrid.Checkbutton(fopt, 'Treat all sub-models as one structure', False)
            oft.button.grid(row = forow, column = 0, sticky = 'w')
            self.lump_subids = oft.variable
            forow += 1

        f = Tkinter.Frame(fopt)
        f.grid(column=0, row=forow, sticky='w')

        l = Tkinter.Label(f, text='Density map resolution:')
        l.grid(column=0, row=0, sticky='w')

        self.simRes = Tkinter.StringVar(fopt)
        e = Tkinter.Entry(f, width=4, textvariable=self.simRes)
        e.grid(column=1, row=0, sticky='w', padx=5)

        l = Tkinter.Label(f, text='grid spacing:')
        l.grid(column=2, row=0, sticky='w')

        self.simGridSp = Tkinter.StringVar(fopt)
        e = Tkinter.Entry(f, width=4, textvariable=self.simGridSp)
        e.grid(column=3, row=0, sticky='w', padx=5)

        b = Tkinter.Button(f, text="Calculate Map", command=self.GenStrucMap)
        b.grid (column=4, row=0, sticky='w', padx=5)


        forow += 1

        #l = Tkinter.Label(fopt, text='Fit to:')
        #l.grid(column=0, row=forow, sticky='w')

        forow += 1

        if 1 :
            f = Tkinter.Frame(fopt)
            f.grid(column=0, row=forow, sticky='w')

            l = Tkinter.Label(f, text='Fit to:')
            l.grid(column=0, row=0, sticky='w')

            self.alignTo = Tkinter.StringVar()
            self.alignTo.set ( 'combined_selected_regions' )

            #l = Tkinter.Label(f, text=' ', width=5)
            #l.grid(column=0, row=0, sticky='w')

            c = Tkinter.Radiobutton(f, text="Combined selected regions", variable=self.alignTo, value = 'combined_selected_regions')
            c.grid (column=1, row=0, sticky='w')

            c = Tkinter.Radiobutton(f, text="Each selected region", variable=self.alignTo, value = 'each_selected_region')
            c.grid (column=2, row=0, sticky='w')

            #c = Tkinter.Radiobutton(f, text="Groups of regions including selected region(s)", variable=self.alignTo, value = 'around_selected')
            #c.grid (column=1, row=2, sticky='w')

            #c = Tkinter.Radiobutton(f, text="Groups of regions including all regions", variable=self.alignTo, value = 'all_groups')
            #c.grid (column=1, row=3, sticky='w')

        forow += 1

        if 1 :
            f = Tkinter.Frame(fopt)
            f.grid(column=0, row=forow, sticky='w')

            l = Tkinter.Label(f, text='Fit by:')
            l.grid(column=0, row=0, sticky='w')

            self.rotaSearch = Tkinter.IntVar()
            self.rotaSearch.set ( 0 )

            #l = Tkinter.Label(f, text=' ', width=5)
            #l.grid(column=0, row=0, sticky='w')

            c = Tkinter.Radiobutton(f, text="PCA (faster)", variable=self.rotaSearch, value = 0)
            c.grid (column=1, row = 0, sticky='w')

            #l = Tkinter.Label(f, text=' ', width=5)
            #l.grid(column=0, row=0, sticky='w')

            c = Tkinter.Radiobutton(f, text="Centers +", variable=self.rotaSearch, value = 1)
            c.grid (column=2, row = 0, sticky='w')

            self.rotaSearchNum = Tkinter.StringVar(f, "100")
            e = Tkinter.Entry(f, width=5, textvariable=self.rotaSearchNum)
            e.grid(column=3, row=0, sticky='w', padx=5)

            l = Tkinter.Label(f, text='rotations (more accurate)')
            l.grid(column=4, row=0, sticky='w')


        forow += 1

        oft = Hybrid.Checkbutton(fopt, 'Mask map with region(s) to prevent large drifts', False)
        oft.button.grid(row = forow, column = 0, sticky = 'w')
        self.mask_map_when_fitting = oft.variable


        forow += 1

        if 1 :
            oft = Hybrid.Checkbutton(fopt, 'Use Laplacian filter', False)
            #oft.button.grid(row = forow, column = 0, sticky = 'w')
            self.useLaplace = oft.variable


        forow += 1

        oft = Hybrid.Checkbutton(fopt, 'Optimize fits', True)
        oft.button.grid(row = forow, column = 0, sticky = 'w')
        self.optimize_fits = oft.variable


        forow += 1

        f = Tkinter.Frame(fopt)
        f.grid(column=0, row=forow, sticky='w')

        oft = Hybrid.Checkbutton(f, 'Cluster fits that are <', True)
        oft.button.grid(row = 0, column = 0, sticky = 'w')
        self.doClusterFits = oft.variable

        self.positionTolString = Tkinter.StringVar(f, "5.0")
        e = Tkinter.Entry(f, width=3, textvariable=self.positionTolString)
        e.grid(column=1, row=0, sticky='w', padx=5)

        l = Tkinter.Label(f, text='Angstroms and <')
        l.grid(column=2, row=0, sticky='w')

        self.angleTolString = Tkinter.StringVar(f, "3.0")
        e = Tkinter.Entry(f, width=3, textvariable=self.angleTolString)
        e.grid(column=3, row=0, sticky='w', padx=5)

        l = Tkinter.Label(f, text='degrees apart' )
        l.grid(column=4, row=0, sticky='w')


        forow += 1

        f = Tkinter.Frame(fopt)
        f.grid(column=0, row=forow, sticky='w')

        l = Tkinter.Label(f, text='Add top')
        l.grid(column=0, row=0, sticky='w')

        self.numFitsToAdd = Tkinter.StringVar(f, "1")
        e = Tkinter.Entry(f, width=5, textvariable=self.numFitsToAdd)
        e.grid(column=1, row=0, sticky='w', padx=5)

        l = Tkinter.Label(f, text='fit(s) to list (empty to add all fits to list)')
        l.grid(column=2, row=0, sticky='w')


        forow += 1

        f = Tkinter.Frame(fopt)
        f.grid(column=0, row=forow, sticky='w')

        oft = Hybrid.Checkbutton(f, 'Clashes using symmetry:', False)
        oft.button.grid(row = 0, column = 0, sticky = 'w')
        self.calcSymmetryClashes = oft.variable

        self.symmetryString = Tkinter.StringVar(f)
        e = Tkinter.Entry(f, width=10, textvariable=self.symmetryString)
        e.grid(column=1, row=0, sticky='w', padx=5)

        b = Tkinter.Button(f, text="Detect", command=self.DetectSym)
        b.grid (column=2, row=0, sticky='w', padx=5)

        b = Tkinter.Button(f, text="Show", command=self.PlaceSym)
        b.grid (column=3, row=0, sticky='w', padx=5)

        b = Tkinter.Button(f, text="Place", command=self.PlaceSym2)
        b.grid (column=4, row=0, sticky='w', padx=5)


        dummyFrame = Tkinter.Frame(parent, relief='groove', borderwidth=1)
        Tkinter.Frame(dummyFrame).pack()
        dummyFrame.grid(row=row,column=0,columnspan=7, pady=7, sticky='we')

        row = row + 1

        global msg
        msg = Tkinter.Label(parent, width = 60, anchor = 'w', justify = 'left', fg="red")
        msg.grid(column=0, row=row, sticky='ew')
        self.msg = msg
        row += 1

        umsg ( 'To cite Segger in your paper, please press the Help button for more information.')

        self.SetResolution()

        chimera.openModels.addRemoveHandler(self.ModelClosed, None)

        mlist = OML(modelTypes = [chimera.Molecule])
        if mlist:
            #self.struc.set(self.menu_name(mlist[0]))
            self.cur_mol = mlist[0]
            print " - set fit mol/map: %s" % self.cur_mol.name
            self.struc.set ( self.cur_mol.name + " (%d)" % self.cur_mol.id)

        if dev_menus :
            self.optionsPanel.set(True)


        self.saveFrames = False
        self.frameAt = 0




    def PlotFitScores ( self ) :


        # self.list_fits.append((fmap, dmap, fmap.M, corr, atomI, bbI, bbC, hdo, regions))

        print "Plotting %d fits:" % len ( self.list_fits )

        if len ( self.list_fits ) == 0 :
            umsg ( "No fits in list to plot" )
            return


        minCorr = 1e9
        maxCorr = 0
        for fmap, dmap, M, corr, atomI, bbI, bbC, hdo, regions in self.list_fits :
            if maxCorr < corr : maxCorr = corr
            if minCorr > corr : minCorr = corr

        print " - maxCorr: %.3f minCorr: %.3f" % (maxCorr, minCorr)

        minCorr, maxCorr = 0, 1

        w = 600
        h = 400
        import PIL
        from PIL import Image, ImageDraw

        im = PIL.Image.new ( 'RGB', (w,h), (255,255,255) )

        #im.putpixel ( (i,j), (fclr[0], fclr[1], fclr[2]) )

        chartW = w - 40
        chartH = h - 40
        chartX = 20
        chartY = 20

        draw = ImageDraw.Draw(im) # Create a draw object

        lineClr = (120,120,120)
        draw.rectangle((10, h-10, w-10, h-10), fill=lineClr, outline=lineClr)
        draw.rectangle((10, 10, 10, h-10), fill=lineClr, outline=lineClr)

        xAt = chartX
        for fmap, dmap, M, corr, atomI, bbI, bbC, hdo, regions in self.list_fits :

            barWidth = int ( float(chartW) / float(len(self.list_fits)) )
            xPos = int ( xAt + barWidth/2 )

            x1 = xAt
            x2 = xAt + barWidth
            if ( barWidth > 3 ) :
                x1 = x1 + 1
                x2 = x2 - 2

            yTop = int ( h - ( chartY + (corr - minCorr) * chartH / (maxCorr - minCorr) ) )
            yBot = int ( h - chartY )

            lineClr = ( int(rand()*255.0), int(rand()*255.0), int(rand()*255.0) )
            draw.rectangle((x1, yTop, x2, yBot), fill=lineClr, outline=lineClr)

            print "x:%d barW %.2f cor %.2f height %.3f yTop %d yBot %d" % ( xAt, barWidth, corr, chartH, yTop, yBot )

            xAt = xAt + barWidth


        #im.save ( "plot.png", 'PNG' )

        def save ( okay, dialog ):
            if okay:
                paths = dialog.getPaths ( )
                if paths:
                    path = paths[0]
                    umsg ( "Saved to " + path )
                    im.save ( path, 'PNG' )

        idir = None
        ifile = None

        from OpenSave import SaveModeless
        SaveModeless ( title = 'Save Plot',
                       filters = [('PNG', '*.png', '.png')],
                       initialdir = idir, initialfile = ifile, command = save )






    def PlotFitScores_Experimental ( self ) :

        print "Plotting fits:"

        N = int ( self.numFitsToAdd.get() )

        totAngles = 0
        minCorr = 1e9
        maxCorr = 0
        for corr, M, regions, stats in self.cfits [0 : N-1] :
            print " - #fits:%d maxAngle:%.1f maxShift:%.2f maxHeight:%.2f" % ( stats['numFits'], stats['maxAngle'], stats['maxShift'], stats['maxHeight'] )
            totAngles = totAngles + stats['maxAngle']
            if maxCorr < corr : maxCorr = corr
            if minCorr > corr-stats['maxHeight'] : minCorr = corr-stats['maxHeight']

        minCorr, maxCorr = 0, 1

        print " - totAngles: %.2f, maxCorr: %.3f minCorr: %.3f" % (totAngles, maxCorr, minCorr)

        w = 600
        h = 400
        import PIL
        from PIL import Image, ImageDraw

        im = PIL.Image.new ( 'RGB', (w,h), (255,255,255) )

        #im.putpixel ( (i,j), (fclr[0], fclr[1], fclr[2]) )

        chartW = w - 40
        chartH = h - 40
        chartX = 20
        chartY = 20

        draw = ImageDraw.Draw(im) # Create a draw object

        lineClr = (120,120,120)
        draw.rectangle((10, h-10, w-10, h-10), fill=lineClr, outline=lineClr)
        draw.rectangle((10, 10, 10, h-10), fill=lineClr, outline=lineClr)

        xAt = chartX
        for corr, M, regions, stats in self.cfits [0 : N-1] :

            barWidth = int ( max ( 1, numpy.floor ( stats['maxAngle']  * float(chartW) / totAngles ) ) )
            xPos = int ( xAt + barWidth/2 )

            x1 = xAt
            x2 = xAt + barWidth
            if ( barWidth > 3 ) :
                x1 = x1 + 1
                x2 = x2 - 2

            yTop = int ( h - ( chartY + (corr - minCorr) * chartH / (maxCorr - minCorr) ) )
            yBot = int ( h - ( chartY + (corr - stats['maxHeight'] - minCorr) * chartH / (maxCorr - minCorr) ) )

            lineClr = ( int(rand()*255.0), int(rand()*255.0), int(rand()*255.0) )
            draw.rectangle((x1, yTop, x2, yBot), fill=lineClr, outline=lineClr)

            print "maxA %.2f x:%d barW %.2f cor %.2f height %.3f yTop %d yBot %d" % ( stats['maxAngle'], xAt, barWidth, corr, stats['maxHeight'], yTop, yBot )

            xAt = xAt + barWidth


        #im.save ( "plot.png", 'PNG' )

        def save ( okay, dialog, lfits = lfits ):
            if okay:
                paths = dialog.getPaths ( )
                if paths:
                    path = paths[0]
                    umsg ( "Saved to " + path )
                    im.save ( path, 'PNG' )

        idir = None
        ifile = None

        from OpenSave import SaveModeless
        SaveModeless ( title = 'Save Plot',
                       filters = [('PNG', '*.png', '.png')],
                       initialdir = idir, initialfile = ifile, command = save )




    def ExportFitScores ( self ) :

        num = self.fit_listbox.size()
        if num == 0 :
            umsg ( "No fits to export" )
            return

        scores = []

        # (fmap, dmap, fmap.M, corr, atomI, bbI, regions)

        for lf in self.list_fits :
            scores.append ( [ lf[3], lf[4], lf[5], lf[6], lf[7] ] )

        def ZZ ( scs ) :
            if len(scs) < 3 :
                return 0.0

            best_score = scs[0]
            other_scores = scs[1:14]
            avg = numpy.average ( other_scores )
            stdev = numpy.std ( other_scores )
            return [(best_score - avg) / stdev, best_score, avg, stdev]

        def save ( okay, dialog, scores = scores ):
            if okay:
                paths = dialog.getPaths ( )
                if paths:
                    path = paths[0]
                    f = open ( path, "a" )

                    f.write ( "%s\t%s\t%s\t%s\t%s\n" % (
                        "Cross-correlation",
                        "Atom Inclusion",
                        "Backbone-Atom Inclusion",
                        "Clash score",
                        "Density occupancy" ) )

                    s1, s2, s3, s4, s5 = [], [], [], [], []
                    for s in scores :
                        f.write ( "%f\t%f\t%f\t%f\t%f\n" % (s[0], s[1], s[2], s[3], s[4]) )
                        s1.append ( s[0] )
                        s2.append ( s[1] )
                        s3.append ( s[2] )
                        s4.append ( s[3] )
                        s5.append ( s[4] )

                    Z1, Z2, Z3, Z4, Z5 = ZZ(s1), ZZ(s2), ZZ(s3), ZZ(s4), ZZ(s5)
                    #f.write ( "Zscores: %f\t%f\t%f\t%f\t%f\n" % (Z1[0], Z2[0], Z3[0], Z4[0], Z5[0]) )
                    f.write ( "Score\tZ-score\tTop score\tMean\tSTDev\n" )
                    f.write ( "Cross-correlation:\t%f\t%f\t%f\t%f\n" % (Z1[0], Z1[1], Z1[2], Z1[3]) )
                    f.write ( "Atom Inclusion:\t%f\t%f\t%f\t%f\n" % (Z2[0], Z2[1], Z2[2], Z2[3]) )
                    f.write ( "Density occupancy:\t%f\t%f\t%f\t%f\n" % (Z5[0], Z5[1], Z5[2], Z5[3]) )
                    f.write ( "Clash score:\t%f\t%f\t%f\t%f\n" % (Z4[0], Z4[1], Z4[2], Z4[3]) )


                    f.close ()
                    umsg ( "Wrote %d fits to %s" % ( len(scores), path ) )


        idir = None
        ifile = None

        first_fit = self.list_fits[0]
        fit_map = first_fit[0]
        ref_map = first_fit[1]
        regions = first_fit[8]

        fit_path = ""
        try : fit_path = fit_map.mols[0].openedAs[0]
        except : fit_path = fit_map.data.path

        import os.path
        idir, ifile = os.path.split(fit_path)
        base, suf = os.path.splitext(ifile)
        map_base, map_suf = os.path.splitext( ref_map.name )
        ifile = base + "_fits_in_%s_regs" % map_base

        for r in regions :
            ifile = ifile + ("_%d" % r.rid)

        from OpenSave import SaveModeless
        SaveModeless ( title = 'Save Fit Scores',
                       filters = [('TXT', '*.txt', '.txt')],
                       initialdir = idir, initialfile = ifile, command = save )


    def add_fit (self, fmap, dmap):

        corr = fmap.fit_score
        regions = fmap.fit_regions
        atomI = fmap.atomInclusion
        bbI = fmap.bbAtomInclusion
        bbC = fmap.bbClashes
        hdo = fmap.hdoScore

        ids = ','.join(['%d' % r.rid for r in regions])
        line = '%8.4f %8.4f %8.4f %8.4f %8.4f %10s %10s %10s' % (corr, atomI, bbI, bbC, hdo, fmap.struc_name, dmap.name, ids)
        self.list_fits.append((fmap, dmap, fmap.M, corr, atomI, bbI, bbC, hdo, regions))
        self.fit_listbox.insert('end', line)

    def fit_selection_cb (self, event) :

        lfits = self.selected_listbox_fits()
        if len(lfits) == 0:
            return

        fmap, dmap, mat, corr, aI, bI, bC, hdo, regions = lfits[0]
        for mol in fmap.mols :
            if mol.__destroyed__:
                umsg('Fit molecule was closed')
        else:
            self.place_molecule(fmap, mat, dmap)
        self.make_regions_transparent(regions)


    def place_molecule(self, fmap, mat, dmap):

        import numpy
        tf = numpy.array(mat[:3,:])
        try :
            xf = dmap.openState.xform
        except :
            print "Reference map no longer open"
            return

        com = chimera_xform(tf).getTranslation()
        q = Segger.quaternion.Quaternion()
        q.fromXform ( chimera_xform(tf) )

        print "COM: ", com, "Q: %.6f" % q.s, q.v

        xf.multiply(chimera_xform(tf))

        try :
            fmap.openState.xform = xf
        except :
            print "Fitted map no longer open"
            return

        fmap.M = mat

        for mol in fmap.mols :
            mol.openState.xform = xf

    def make_regions_transparent(self, regions):

        for r in regions:
            if r.has_surface():
                c = r.color
                r.surface().color = ( c[0], c[1], c[2], REG_OPACITY )

    def selected_listbox_fits(self):

        return [self.list_fits[int(i)] for i in self.fit_listbox.curselection()]


    def delete_all_fit_cb ( self ) :

        num = self.fit_listbox.size()
        indices = range ( num )
        indices.reverse()

        for i in indices:
            self.fit_listbox.delete(i)
            del self.list_fits[i]


    def delete_fit_cb(self):

        indices = [int(i) for i in self.fit_listbox.curselection()]
        if len(indices) == 0:
            status('No fits chosen from list')
            return
        indices.sort()
        indices.reverse()
        for i in indices:
            self.fit_listbox.delete(i)
            del self.list_fits[i]

        status('Deleted %d fits' % len(indices))


    def place_copies_cb(self):

        lfits = self.selected_listbox_fits()
        if len(lfits) == 0:
            status('No fits chosen from list')
            return


        fmap, dmap = lfits[0][0], lfits[0][1]
        dmapM = xf_2_MM ( dmap.openState.xform )

        nPlaced = 0
        for fmap, dmap, mat, corr, aI, bI, bC, bO, regions in lfits:
            if len ( fmap.mols ) > 0 :
                self.PlaceCopy(fmap.mols, mat, dmap, (rand(),rand(),rand(),1) )
                nPlaced += 1

        status('Placed %d models' % nPlaced)




    def place_map_copies_cb ( self ) :

        lfits = self.selected_listbox_fits()
        if len(lfits) == 0 :
            status('No fits chosen from list')
            return

        fmap, dmap = lfits[0][0], lfits[0][1]
        dmapM = xf_2_MM ( dmap.openState.xform )

        for fmap, dmap, mat, corr, aI, bI, bC, bO, regions in lfits:
            self.place_molecule(fmap, mat, dmap)
            sf = "_F2Rid%d.mrc" % regions[0].rid
            pmap = place_map_resample ( fmap, dmap, sf )

            try :
                self.fitted_mols = self.fitted_mols + [pmap]
            except :
                self.fitted_mols = [pmap]


        status('Placed %d models' % len(lfits))


    def save_map_resample ( self ) :

        dmap = segmentation_map()
        if dmap == None : umsg ( "%s is not open" % self.dmap.get() ); return

        fmap = self.cur_mol
        if fmap == None :
            print "Select the map to save"
            return

        if type(fmap) != VolumeViewer.volume.Volume :
            umsg ( "Please select a map to save (molecule is selected)" )

        else :
            place_map_resample ( fmap, dmap, "_F2Rid.mrc" )


    def extractCubeMap ( self ) :

        dmap = segmentation_map()
        if dmap == None : umsg ( "%s is not open" % self.dmap.get() ); return

        fmap = self.cur_mol
        if fmap == None :
            print "Select the map to save"
            return

        print "Saving ", fmap.name


        npoints = grid_indices ( dmap.data.size, numpy.single)  # i,j,k indices
        transform_vertices ( npoints, dmap.data.ijk_to_xyz_transform )

        dvals = fmap.interpolated_values ( npoints, dmap.openState.xform )
        #dvals = numpy.where ( dvals > threshold, dvals, numpy.zeros_like(dvals) )
        #nze = numpy.nonzero ( dvals )

        nmat = dvals.reshape( dmap.data.size )
        #f_mat = fmap.data.full_matrix()
        #f_mask = numpy.where ( f_mat > fmap.surface_levels[0], numpy.ones_like(f_mat), numpy.zeros_like(f_mat) )
        #df_mat = df_mat * f_mask

        ndata = VolumeData.Array_Grid_Data ( nmat, dmap.data.origin, dmap.data.step, dmap.data.cell_angles )
        try : nv = VolumeViewer.volume.add_data_set ( ndata, None )
        except : nv = VolumeViewer.volume.volume_from_grid_data ( ndata )

        fmap_base = os.path.splitext (fmap.name)[0]
        dmap_base = os.path.splitext (dmap.name)[0]
        fmap_path = os.path.splitext (fmap.data.path)[0]
        dmap_path = os.path.splitext (dmap.data.path)[0]

        nv.name = "emd_1093_62.mrc"
        nv.openState.xform = dmap.openState.xform




    def close_copies_cb ( self ) :

        try :
            len ( self.fitted_mols )
        except :
            umsg ( "No fitted molecules found" )
            return

        chimera.openModels.close ( self.fitted_mols )



    def Options(self):

        self.optionsPanel.set(not self.optionsPanel.get())


    def ModelClosed(self, trigger, n, mlist):

        # Clear molecule menu if selected molecule is closed.

        mvar = self.struc
        #if len( mvar.get() ) > 0 and len( self.StructuresToFit() ) == 0:
        #    mvar.set('')

        found = False
        for m in chimera.openModels.list() :
            if m.name + " (%d)" % m.id == mvar.get() :
                found = True

        if not found :
            #print " - closed model - selected model not found", mvar.get()
            mvar.set('')
            self.cur_mol = None


    def SetResolution ( self ):

        dmap = segmentation_map()
        if dmap == None : return

        if len ( self.simRes.get() ) == 0:
            res = min(dmap.data.step) * 3
            self.simRes.set ( '%.3g' % res )
            self.simGridSp.set ( '%.3g' % (res/3.0) )


    def CurrentSegmentation ( self, warn = True ):

        return current_segmentation(warn)


    def StrucMenu0 ( self ) :

        self.strucMB.menu.delete ( 0, 'end' )   # Clear menu

        self.strucMB.menu.add_radiobutton ( label="Open structures:" )
        self.strucMB.menu.add_separator()

        id_struc = {}
        open_mols = {}
        for m in OML() :
            if type(m) == chimera.Molecule or type(m) == VolumeViewer.volume.Volume:
                try : id_struc [ m.id ].append ( m )
                except : id_struc [ m.id ] = [m]
                open_mols[m.name] = True
            else :
                #print type(m)
                pass

        if 1 :
            cur_sel_found = False
            for mid, mols in id_struc.iteritems() :

                if len(mols) == 1 or self.lump_subids.get () :
                    mol = mols[0]
                    label = self.menu_name(mol)
                    self.strucMB.menu.add_radiobutton ( label= label,
                                                        variable=self.struc,
                                                        command = lambda mol=mol: self.StrucSelected(mol) )
                    if label == self.struc.get() :
                        cur_sel_found = True

                else :
                    for mol in mols :
                        label = self.menu_name(mol)
                        self.strucMB.menu.add_radiobutton ( label=label,
                                                            variable=self.struc,
                                                            command = lambda mol=mol: self.StrucSelected(mol) )
                        if label == self.struc.get() :
                            cur_sel_found = True


        if not cur_sel_found :
            self.struc.set ( "" )
            self.cur_mol = None
            print " - set fit mol/map: ?"

        if 0 :
            dmap = segmentation_map()
            if dmap == None : return

            path = os.path.dirname ( dmap.data.path ) + os.path.sep
            files = os.listdir ( path );
            mols_in_path = []
            for f in files :
                if f.find ( ".pdb" ) >= 0 and open_mols.has_key(f) == False :
                    mols_in_path.append ( f )

            if len ( mols_in_path ) == 0 : return

            self.strucMB.menu.add_separator()
            self.strucMB.menu.add_radiobutton ( label="In %s:" % path )
            self.strucMB.menu.add_separator()

            for fm in mols_in_path :
                self.strucMB.menu.add_radiobutton ( label=fm, variable=self.struc,
                                                    command = self.StrucSelected )



    def StrucMenu ( self ) :

        self.strucMB.menu.delete ( 0, 'end' )   # Clear menu

        self.strucMB.menu.add_radiobutton ( label="Open structures:" )
        self.strucMB.menu.add_separator()

        id_struc = {}
        open_mols = {}
        for m in OML() :
            if type(m) == chimera.Molecule or type(m) == VolumeViewer.volume.Volume:

                self.strucMB.menu.add_radiobutton ( label= m.name + " (%d)" % m.id,
                                                    variable = self.struc,
                                                    command = lambda m=m: self.StrucSelected(m, None) )

                open_mols[m.name] = True

        if 0 :
            dmap = segmentation_map()
            if dmap == None : return

            path = os.path.dirname ( dmap.data.path ) + os.path.sep
            files = os.listdir ( path );
            mols_in_path = []
            for f in files :
                if f.find ( ".pdb" ) >= 0 and open_mols.has_key(f) == False :
                    mols_in_path.append ( f )

            if len ( mols_in_path ) == 0 : return

            self.strucMB.menu.add_separator()
            self.strucMB.menu.add_radiobutton ( label="In %s:" % path )
            self.strucMB.menu.add_separator()

            for fm in mols_in_path :
                self.strucMB.menu.add_radiobutton ( label=fm, variable=self.struc,
                                                    command = lambda fm=fm: self.StrucSelected(None, fm) )



    def StrucSelected ( self, mol, molNameToLoad ) :

        # Check if selected entry is an existing open molecule.

        if mol != None :
            self.cur_mol = mol

        if molNameToLoad != None :

            dmap = segmentation_map()
            if dmap == None : print " - no map selected"; return

            path = os.path.dirname(dmap.data.path) + os.path.sep

            print " - opening: %s" % molNameToLoad
            fmol = chimera.openModels.open ( path + molNameToLoad )[0]

            self.struc.set(fmol.name + " (%d)" % fmol.id)
            self.cur_mol = fmol

        print " - set fit mol/map: %s" % self.cur_mol.name


    def StructuresToFit ( self ) :
        #t = self.struc.get()
        #return [m for m in OML() if self.menu_name(m) == t]
        return [self.cur_mol]


    #def menu_name(self, mol):
    #    show_subid = not self.lump_subids.get() and mol.subid != 0
    #    id =  '%d.%d' % (mol.id, mol.subid) if show_subid else '%d' % mol.id
    #    mname = "%s (%s)" % (mol.name, id)
    #    return mname


    def Place ( self ) :

	    self.save_map_resample ();

    def Clear ( self ) :
        self.delete_all_fit_cb ()


    def Fit ( self ) :


        dmap = segmentation_map()
        if dmap == None :
            umsg ( "Density map not found or not selected" )
            return

        if self.cur_mol == None :
            umsg ( "Please select a structure or map to fit" )
            return



        from chimera import tasks, CancelOperation
        task = tasks.Task('Fitting %s to %s' % (self.cur_mol.name, dmap.name), modal = True)
        print "Started task..."

        try:
            smod = self.Fit_(dmap, self.cur_mol, task)
        except CancelOperation:
            umsg('Cancelled fitting')
            return None
        finally:
            task.finished()



    def Fit_ ( self, dmap, fmol, task=None ) :

        self.doStop = False

        #if type(fmol) == VolumeViewer.volume.Volume :
        if type(fmol) == chimera.Molecule :

            # this property added to self will indicate to functions called
            # below whether we are fitting a map instead of a molecule
            self.map_to_fit = None

            if fmol.openState is dmap.openState :
                umsg('Molecule cannot be moved relative to map\nbecause they have the same model id number.')
                return

        else :
            # looks like we are fitting a map, not a molecule
            fmap = fmol

            self.map_to_fit = fmap
            fmap.mols = []
            fmap.struc_name = fmap.name

            points, weights = fit_points ( fmap )
            # print "Points : ", points

            fmap.COM, fmap.U, fmap.S, fmap.V = prAxes ( points )
            print "COM : ", fmap.COM
            print "U : ", fmap.U

            toCOM = numpy.matrix ( [
                [ 1, 0, 0, -fmap.COM[0] ],
                [ 0, 1, 0, -fmap.COM[1] ],
                [ 0, 0, 1, -fmap.COM[2] ],
                [ 0, 0, 0,      1 ]  ] )

            mR = numpy.matrix ( [
                [ fmap.V[0,0], fmap.V[0,1], fmap.V[0,2], 0 ],
                [ fmap.V[1,0], fmap.V[1,1], fmap.V[1,2], 0 ],
                [ fmap.V[2,0], fmap.V[2,1], fmap.V[2,2], 0 ],
                [            0,            0,            0, 1 ] ] )

            # this matrix centers the map and aligns its principal axes
            # to the x-y-z axes
            fmap.preM = mR * toCOM


        alignTo = self.alignTo.get()

        if alignTo == "combined_selected_regions" :
            self.FitMapToSelRGroup ( dmap, task )

        elif alignTo == "each_selected_region" :
            self.FitToEachRegion ( dmap, task )

        elif alignTo == "around_selected" :
            self.FitMapToRegionsAroundSel ( dmap, task )

        elif alignTo == "all_groups" :
            self.FitMapToRGroups ( dmap, task )

        dmap = segmentation_map()
        #if dmap:
        #    dmap.display = False


    def Stop (self) :
        print "will stop..."
        self.doStop = True



    def DetectSym ( self ) :

        dmap = segmentation_map()

        if dmap == None:
            umsg ( "Please select a map in the Segment Map dialog" )
            return []

        print "Symmetry for", dmap.name

        from Measure.symmetry import find_point_symmetry

        syms, msg = find_point_symmetry ( dmap, nMax=8 )

        if syms is None :
            umsg ( "No symmetry detected for %s" % dmap.name )
            self.symmetryString.set ( "No symmetry detected" )
            return []

        umsg ( msg )
        start = msg.find(': ')+2
        end = msg.find (', center')
        self.symmetryString.set ( msg [start : end] )

        for i, sym in enumerate ( syms ) :
            #print i, " -> ", sym
            pass

        return syms


    def PlaceSym ( self ) :

        fmap = self.MoleculeMap()
        dmap = segmentation_map()

        if fmap == None or dmap == None:
            umsg ( "Please select an open structure to fit" )
            return

        from Measure.symmetry import centers_and_points

        syms = []
        esym = self.symmetryString.get()
        if 1 or len (esym) == 0 :
            syms = self.DetectSym ()
        else :
            print "Custom sym:", esym
            if ( esym == "C3" ) :

                print " - dmap: ", dmap.name
                mpoints, mpoint_weights = fit_points(dmap)
                COM, U, S, V = prAxes ( mpoints )
                print "COM: ", COM
                print "U: ", U
                print "S: ", S

                ax = chimera.Vector ( 0, 1, 0 )
                #ax = dmap.openState.xform.inverse().apply ( ax )

                syms.append ( Matrix.identity_matrix () )
                rm1 = Matrix.rotation_transform ( (ax.x,ax.y,ax.z), 360.0/3.0, COM )
                print rm1
                syms.append ( rm1 )
                #syms.append ( Matrix.rotation_transform ( (1.0,0.0,0.0), 2.0*360.0/3.0 ) )

                #centers, xyz, w = centers_and_points(dmap)
                #print " - center:", centers
                #ctf = Matrix.translation_matrix([-x for x in COM[0]])
                #syms = Matrix.coordinate_transform_list(syms, ctf)


        smols = []

        for si, sym in enumerate ( syms [1 : ] ) :

            T = numpy.array ( sym )
            #print "\nSym %d\n" % si, T

            xf = chimera.Xform.xform ( T[0,0], T[0,1], T[0,2], T[0,3], T[1,0], T[1,1], T[1,2], T[1,3], T[2,0], T[2,1], T[2,2], T[2,3] )
            M = xf_2_MM ( xf )

            #mols = self.PlaceCopy (fmap.mols, M*fmap.M, dmap, (0,0,0,1) )
            mols = self.PlaceCopy (fmap.mols, M*fmap.M, dmap, (.4, .8, .4, 1) )
            for m in mols : m.openState.xform = dmap.openState.xform
            smols = smols + mols

        return smols



    def PlaceSym2 ( self ) :

        #fmap = self.MoleculeMap()
        dmap = segmentation_map()


        label = self.struc.get()
        sel_str = "#" + label [ label.rfind("(")+1 : label.rfind(")") ]
        fmol = None
        try :
            fmol = chimera.selection.OSLSelection(sel_str).molecules()[0]
        except :
            umsg ( "%s not open - " % self.struc.get() ); return


        if fmol == None or dmap == None:
            umsg ( "Please select an open structure and/or map" )
            return

        from Measure.symmetry import centers_and_points

        syms = []
        esym = self.symmetryString.get()
        if 1 or len (esym) == 0 :
            syms = self.DetectSym ()
        else :
            print "Custom sym:", esym
            if ( esym == "C3" ) :

                print " - dmap: ", dmap.name
                mpoints, mpoint_weights = fit_points(dmap)
                COM, U, S, V = prAxes ( mpoints )
                print "COM: ", COM
                print "U: ", U
                print "S: ", S

                ax = chimera.Vector ( 0, 1, 0 )
                #ax = dmap.openState.xform.inverse().apply ( ax )

                syms.append ( Matrix.identity_matrix () )
                rm1 = Matrix.rotation_transform ( (ax.x,ax.y,ax.z), 360.0/3.0, COM )
                print rm1
                syms.append ( rm1 )
                #syms.append ( Matrix.rotation_transform ( (1.0,0.0,0.0), 2.0*360.0/3.0 ) )

                #centers, xyz, w = centers_and_points(dmap)
                #print " - center:", centers
                #ctf = Matrix.translation_matrix([-x for x in COM[0]])
                #syms = Matrix.coordinate_transform_list(syms, ctf)


        from SWIM import SetBBAts
        SetBBAts ( fmol )

        if 0 :
            smols = []

            #mol = fmap.mols[0]
            cid = fmol.residues[0].id.chainId
            print "Symming %s, chain %s" % (fmol.name, cid)
            nmol = CopyChain ( fmol, None, cid, cid, dmap.openState.xform.inverse() )
            chimera.openModels.add ( [nmol] )

            chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZ01234567890abcdefghijklmnopqrstuvwxyz"
            atCi = 0

            for si, sym in enumerate ( syms [1 : ] ) :

                T = numpy.array ( sym )
                #print "\nSym %d\n" % si, T

                xf = chimera.Xform.xform ( T[0,0], T[0,1], T[0,2], T[0,3], T[1,0], T[1,1], T[1,2], T[1,3], T[2,0], T[2,1], T[2,2], T[2,3] )
                #M = xf_2_MM ( xf )

                if chains[atCi] == cid :
                    atCi += 1
                ncid = chains[atCi]

                xf1 = dmap.openState.xform.inverse()
                xf1.premultiply (xf)

                print " - %d - %s" % (atCi, ncid)

                CopyChain ( fmol, nmol, cid, ncid, xf1 )

                atCi += 1
                #mols = self.PlaceCopy (fmap.mols, M*fmap.M, dmap, (0,0,0,1) )
                #mols = self.PlaceCopy (fmap.mols, M*fmap.M, dmap, (.4, .8, .4, 1) )

                #for m in mols :
                #    m.openState.xform = dmap.openState.xform
                #smols = smols + mols

                #break

            return smols

        else :

            cmap = {}
            for r in fmol.residues :
                cmap[r.id.chainId] = 1

            chains = cmap.keys()

            for si, sym in enumerate ( syms [1 : ] ) :

                T = numpy.array ( sym )
                xf = chimera.Xform.xform ( T[0,0], T[0,1], T[0,2], T[0,3], T[1,0], T[1,1], T[1,2], T[1,3], T[2,0], T[2,1], T[2,2], T[2,3] )
                xf1 = dmap.openState.xform.inverse()
                xf1.premultiply (xf)

                nmol = CopyMolX ( fmol, xf1 )
                nmol.name = fmol.name + "__sym%d" % si
                chimera.openModels.add ( [nmol] )
                print ".",

            print ""




    def PlaceSymOld ( self ) :

        fmap = self.MoleculeMap()
        dmap = segmentation_map()

        if fmap == None or dmap == None:
            print "Fit or segmentation map not found"
            return

        fpoints = grid_indices(dmap.data.size, numpy.single)  # i,j,k indices
        transform_vertices( fpoints, dmap.data.ijk_to_xyz_transform )
        mat = dmap.data.full_matrix()
        fpoint_weights = numpy.ravel(mat).astype(numpy.single)

        threshold = dmap.surface_levels[0]
        ge = numpy.greater_equal(fpoint_weights, threshold)
        fpoints = numpy.compress(ge, fpoints, 0)
        fpoint_weights = numpy.compress(ge, fpoint_weights)
        nz = numpy.nonzero( fpoint_weights )[0]

        print "%d above %f in %s\n" % (len(nz), threshold, dmap.name)

        COM, U, S, V = prAxes ( fpoints )

        print "COM: ", COM
        print "U: ", U
        print "S: ", S

        T0 = numpy.matrix ( [
            [ 1, 0, 0, -COM[0] ],
            [ 0, 1, 0, -COM[1] ],
            [ 0, 0, 1, -COM[2] ],
            [ 0, 0, 0,      1 ]  ] )

        T = numpy.matrix ( [
            [ 1, 0, 0, COM[0] ],
            [ 0, 1, 0, COM[1] ],
            [ 0, 0, 1, COM[2] ],
            [ 0, 0, 0,      1 ]  ] )


        fmapM = xf_2_MM ( fmap.openState.xform )
        dmapM = xf_2_MM ( dmap.openState.xform )


        smols = []

        if 0 :
            M = xf_2_MM ( chimera.Xform.rotation( 0, 0, 1, 360.0/7.0 ) )
            mols = self.PlaceCopy (fmap.mols, M*fmap.M, dmap, (0,0,0,1) )
            for m in mols : m.openState.xform = dmap.openState.xform
            smols = smols + mols

            M = xf_2_MM ( chimera.Xform.rotation( 0, 0, 1, -360.0/7.0 ) )
            mols = self.PlaceCopy (fmap.mols, M*fmap.M, dmap, (0,0,0,1))
            for m in mols : m.openState.xform = dmap.openState.xform
            smols = smols + mols

            M1 = xf_2_MM ( chimera.Xform.rotation( 1, 0, 0, 180.0 ) )
            M = xf_2_MM ( chimera.Xform.rotation( 0, 0, 1, 2.0*360.0/7.0 ) )
            mols = self.PlaceCopy (fmap.mols, M*M1*fmap.M, dmap, (0,0,0,1))
            for m in mols : m.openState.xform = dmap.openState.xform
            smols = smols + mols

            M1 = xf_2_MM ( chimera.Xform.rotation( 0, 0, 1, 180.0 ) )
            M = xf_2_MM ( chimera.Xform.rotation( 0, 0, 1, 3.0*360.0/7.0 ) )
            mols = self.PlaceCopy (fmap.mols, M*M1*fmap.M, dmap, (0,0,0,1) )
            for m in mols : m.openState.xform = dmap.openState.xform
            smols = smols + mols


        M = xf_2_MM ( chimera.Xform.rotation( U[0,2], U[1,2], U[2,2], 360.0/7.0 ) )
        mols = self.PlaceCopy (fmap.mols, T*M*T0*fmap.M, dmap, (0,0,0,1) )
        for m in mols : m.openState.xform = dmap.openState.xform
        smols = smols + mols

        M = xf_2_MM ( chimera.Xform.rotation( U[0,2], U[1,2], U[2,2], -360.0/7.0 ) )
        mols = self.PlaceCopy (fmap.mols, T*M*T0*fmap.M, dmap, (0,0,0,1) )
        for m in mols : m.openState.xform = dmap.openState.xform
        smols = smols + mols

        return smols


    def StrucCenter ( self ) :

        label = self.struc.get()
        sel_str = "#" + label [ label.rfind("(")+1 : label.rfind(")") ]
        mols = centerMol ( sel_str )

        mm = segmentation_map()
        if mm :
            for mol in mols :
                mol.openState.xform = mm.openState.xform



    def StrucShowAxes ( self ) :

        if len ( self.struc.get() ) == 0 :
            print "Please select a structure"
            return

        label = self.struc.get()
        sel_str = "#" + label [ label.rfind("(")+1 : label.rfind(")") ]
        mols = centerMol ( sel_str )

        if len(mols) == 0:
            print self.struc.get(), "not open";
            return

        fmol = mols[0]
        print "Showing axes for", fmol.name
        print " - COM:", fmol.COM
        print " - extents:", fmol.Extents

        try :
            chimera.openModels.close ( fmol.axes )
            fmol.axes = None
        except :
            pass

        import axes
        fmol.axes =  axes.AxesMod ( Extents = fmol.Extents, rad = 1.0,
                                    alignTo = fmol )
        fmol.axes.name = os.path.splitext (fmol.name)[0] + "_axes"



    def StrucHideAxes ( self ) :

        label = self.struc.get()
        sel_str = "#" + label [ label.rfind("(")+1 : label.rfind(")") ]
        mols = centerMol ( sel_str )

        if len(mols) == 0 :
            print self.struc.get(), "not open"; return

        fmol = mols[0]

        try :
            chimera.openModels.close ( fmol.axes )
            fmol.axes = None
        except :
            pass


    def GenStrucMap ( self, show = True ) :

        self.SetResolution()

        try : res = float ( self.simRes.get() )
        except :
            umsg ( "Invalid resolution entered, please enter a number" )
            return

        try : grid = float ( self.simGridSp.get() )
        except :
            umsg ( "Invalid grid spacing entered, using resolution/3.0" )
            grid = res/3.0

        umsg ( "Simulating map res %.3f, grid %.3f" % (res, grid) )

        if len(self.struc.get()) == 0 :
            umsg ( "Please select a Structure to fit in the field above" )
            return

        label = self.struc.get()
        sel_str = label [ label.rfind("(")+1 : label.rfind(")") ]
        mols = centerMol ( "#" + sel_str )
        if len(mols) == 0 :
            umsg ( "%s not open" % self.struc.get() )
            return

        mol = mols[0]

        base = os.path.splitext(mol.name)[0]
        mname = base + "_" + sel_str + "_r%.1f_sp%.1f" % (res, grid)
        #if self.useLaplace.get() :
        #    mname = mname + "_L"
        mname = mname + ".mrc"

        mv = getMod ( mname )
        if mv != None :
            print "Found", mname
            return mv

        print "Generating", mname

        #cmd = "molmap #%s:.C-D@CA %f sigmaFactor 0.187 gridSpacing %f replace false" % ( sel_str, res, grid )
        cmd = "molmap #%s %f sigmaFactor 0.187 gridSpacing %f replace false" % ( sel_str, res, grid )
        print " -", cmd
        chimera.runCommand ( cmd )

        mv = None
        for mod in chimera.openModels.list() :
            ts = mod.name.split()
            if len(ts) > 1 and mod.name.find("map") >=0 and mod.name.find("res") >=0 :
                print " - found", mod.name
                mv = mod
                mv.name = mname
                break

        if mv == None :
            umsg ("Map not generated - molmap command did not produce expected result.")
            return

        if 0 or self.useLaplace.get() :
            umsg ("Generating Laplacian...")
            from VolumeFilter import laplacian
            mvl = laplacian ( mv )
            chimera.openModels.close ( [mv] )
            mv = mvl
            mv.name = mname

        mv.display = show
        clr = mol.color.rgba()
        mv.surfacePieces[0].color = ( clr[0], clr[1], clr[2], 1.0 )
        mv.mols = mols
        mv.struc_name = label

        # for consistency when fitting maps, which need this pre-transform
        # since they don't get transformed like the molecules
        mv.preM = numpy.matrix ( [
            [ 1, 0, 0, 0 ],
            [ 0, 1, 0, 0 ],
            [ 0, 0, 1, 0 ],
            [ 0, 0, 0, 1 ]  ] )

        return mv




    def GenChMaps ( self, m ) :

        try : res = float ( self.simRes.get() )
        except :
            umsg ( "Invalid resolution entered, please enter a number" )
            return

        try : grid = float ( self.simGridSp.get() )
        except :
            umsg ( "Invalid grid spacing entered, using resolution/3.0" )
            grid = res/3.0




        if len(self.struc.get()) == 0 :
            umsg ( "Please select a Structure to fit in the field above" )
            return

        label = self.struc.get()
        sel_str = label [ label.rfind("(")+1 : label.rfind(")") ]
        mols = centerMol ( "#" + sel_str )
        if len(mols) == 0 :
            umsg ( "%s not open" % self.struc.get() )
            return

        mol = mols[0]

        base = os.path.splitext(mol.name)[0]
        mname = base + "_" + sel_str + "_r%.1f_sp%.1f" % (res, grid)
        if 0 and self.useLaplace.get() :
            mname = mname + "_L"
        mname = mname + ".mrc"

        mv = getMod ( mname )
        if mv != None :
            print "Found", mname
            return mv

        print "Generating", mname

        #cmd = "molmap #%s:.C-D@CA %f sigmaFactor 0.187 gridSpacing %f replace false" % ( sel_str, res, grid )
        cmd = "molmap #%s %f sigmaFactor 0.187 gridSpacing %f replace false" % ( sel_str, res, grid )
        print " -", cmd
        chimera.runCommand ( cmd )

        mv = None
        for mod in chimera.openModels.list() :
            ts = mod.name.split()
            if len(ts) > 1 and mod.name.find("map") >=0 and mod.name.find("res") >=0 :
                print " - found", mod.name
                mv = mod
                mv.name = mname
                break

        if mv == None :
            umsg ("Map not generated - molmap command did not produce expected result.")
            return

        if 0 and self.useLaplace.get() :
            umsg ("Generating Laplacian...")
            from VolumeFilter import laplacian
            mvl = laplacian ( mv )
            chimera.openModels.close ( [mv] )
            mv = mvl
            mv.name = mname

        mv.display = show
        clr = mol.color.rgba()
        mv.surfacePieces[0].color = ( clr[0], clr[1], clr[2], 1.0 )
        mv.mols = mols
        mv.struc_name = label

        # for consistency when fitting maps, which need this pre-transform
        # since they don't get transformed like the molecules
        mv.preM = numpy.matrix ( [
            [ 1, 0, 0, 0 ],
            [ 0, 1, 0, 0 ],
            [ 0, 0, 1, 0 ],
            [ 0, 0, 0, 1 ]  ] )

        return mv





    def SaveStrucFit ( self ) :

        lfits = self.selected_listbox_fits()
        if len(lfits) == 0:
            status('No fits chosen from list')
            return

        def save ( okay, dialog, lfits = lfits ):
            if okay:
                paths = dialog.getPaths ( )
                if paths:
                    path = paths[0]
                    import Midas
                    if len(lfits) > 1 and path.find('%d') == -1:
                        base, suf = os.path.splitext(path)
                        path = base + '_fit%d' + suf
                    for i, (fmap, dmap, mat, corr, aI, bI, bC, hdo, regions) in enumerate(lfits):
                        p = path if len(lfits) == 1 else path % (i+1)
                        self.place_molecule(fmap, mat, dmap)
                        Midas.write(fmap.mols, relModel = dmap, filename = p)

        mol = lfits[0][0].mols[0]
        if hasattr(mol, 'openedAs'):
            import os.path
            idir, ifile = os.path.split(mol.openedAs[0])
            base, suf = os.path.splitext(ifile)
            if len(lfits) > 1:
                ifile = base + '_fit%d' + suf
            else:
                ifile = base + '_fit' + suf
        else:
            idir = None
            ifile = None

        from OpenSave import SaveModeless
        SaveModeless ( title = 'Save Fit Molecules',
                       filters = [('PDB', '*.pdb', '.pdb')],
                       initialdir = idir, initialfile = ifile, command = save )

    def SaveFit ( self, fmap, clr=None ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        mols = self.PlaceCopy(fmap.mols, fmap.M, dmap, clr)
        path = os.path.dirname ( dmap.data.path ) + os.path.sep
        print "Saving:"
        for mol in mols :
            print " - %s %d.%d" % ( path + mol.name, mol.id, mol.subid )
        chimera.PDBio().writePDBfile ( mols, path + mols[0].name )
        umsg ( "Saved fit (%d structures)" % len(mols) )

        return mols



    def PlaceCopy(self, molecules, mat, dmap, clr=None):

        try : fit_m_at = len ( self.fitted_mols ) + 1
        except : fit_m_at = 1; self.fitted_mols = []

        new_mols = []

        for molecule in molecules :

            mol = CopyMol ( molecule )
            mol.name = os.path.splitext ( mol.name )[0] + "_f%d.pdb" % fit_m_at

            if clr :
                r, g, b, a = clr
                mclr = chimera.MaterialColor ( r, g, b, a )
            else : mclr = molecule.residues[0].ribbonColor

            if mat != None : molApplyT ( mol, mat )
            new_mols.append ( mol )

            for res in mol.residues :
                res.ribbonDisplay = True
                res.ribbonDrawMode = 2
                res.ribbonColor = mclr
                for at in res.atoms : at.display = False

            mol.display = True

        if dmap :
            self.fitted_mols = self.fitted_mols + new_mols

        chimera.openModels.add ( new_mols, noprefs = True )

        return new_mols



    def FitToEachRegion ( self, dmap, task=None ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        smod = self.CurrentSegmentation()
        if smod is None : return

        fmap = None
        if self.map_to_fit :
            fmap = self.map_to_fit
        else :
            fmap = self.MoleculeMap()
            #if not hasattr(fmap.mols[0], 'centered'):
            self.StrucCenter()

        regs = smod.selected_regions()
        if len(regs) == 0 :
            umsg ( "Please select one or more regions to align the structure to" )
            return

        path = os.path.dirname ( dmap.data.path ) + os.path.sep

        for reg in regs :

            self.fits = []
            reportFitRegions(fmap.name, regs)
            scores = []
            corrs = []

            sp = reg.surface_piece
            sp.display = True
            clr = sp.region.color
            sp.color = ( clr[0], clr[1], clr[2], REG_OPACITY )

            #fmap.display = False
            for mol in fmap.mols : mol.display = True

            fmap.fit_regions = [reg]

            to_map = dmap
            reg_map = None
            if 1 :
                reg_map = mask_volume( [reg], dmap )

            bCloseMap = False

            if 0 and self.useLaplace.get () :
                umsg ("Generating Laplace version of " + dmap.name)
                from VolumeFilter import laplacian
                to_map = laplacian ( dmap )
                bCloseMap = True

            elif self.mask_map_when_fitting.get() :
                to_map = reg_map
                bCloseMap = True
                if to_map is None:
                    umsg ('Could not create masked map')
                    return

            tpoints = reg.map_points()
            if self.rotaSearch.get () :
                self.saFitMapToPoints_byRot ( fmap, tpoints, to_map )
            else :
                self.saFitMapToPoints ( fmap, tpoints, to_map )

            scores.append ( 0 )
            corrs.append ( fmap.fit_score )

            umsg ( "Cross-correlation of fit: %f" % fmap.fit_score )

            self.cfits = self.ClusterFits ( self.fits )
            self.cfits.sort ( reverse=True, key=lambda x: x[0] )
            #self.cfits.sort()
            #self.cfits.reverse()

            #frame_at = 0

            try : nToAdd = int ( self.numFitsToAdd.get() )
            except : nToAdd = len (self.cfits)
            for corr, M, regions, stats in self.cfits [ 0 : nToAdd ] :
                print " -- clustered Fit -- #fits:%d maxAngle:%.1f maxShift:%.2f maxHeight:%.2f" % ( stats['numFits'], stats['maxAngle'], stats['maxShift'], stats['maxHeight'] )
                fmap.fit_score, fmap.M, fmap.fit_regions = corr, M, regions
                fmap.atomInclusion, fmap.bbAtomInclusion, fmap.bbClashes, fmap.hdoScore = self.FitScores ( fmap, reg_map )
                #frame_at = frame_at + 1
                self.add_fit (fmap, dmap)

            # FitScores fn moves the model, so move it back to top fit
            move_fit_models(fmap, self.cfits[0][1], dmap.openState.xform)

            #umsg ( "Cross-correlation: %.4f\n" % (fmap.fit_score) )
            self.ReportZScore ( self.cfits )

            # close masked map if it was created
            #if bCloseMap : chimera.openModels.close ( [to_map] )
            if reg_map :
                chimera.openModels.close ( [reg_map] )


    def saFitMapToPoints ( self, fmap, points, dmap, task=None ) :

        print "fitting %s in map %s, to %d points" % (fmap.name, dmap.name, len(points))

        fpoints, fpoint_weights = fit_points(fmap)

        # the 4 alignments to try...
        flips = [ (1,1,1), (-1,-1,1), (-1,1,-1), (1,-1,-1) ]
        #flips = [ (1,-1,-1)  ]
        mlist = principle_axes_alignments ( points, flips, fmap.preM )

        best = (-2, None, None)
        names =  ['%.1f*X %.1f*Y %.1f*Z' % f for f in flips]

        optimize = self.optimize_fits.get()

        fits = optimize_fits(fpoints, fpoint_weights, mlist, dmap, names, None, optimize)
        corr, Mfit, i = self.make_best_fit(fits, fmap, dmap)

        f = flips[i]
        print " - best fit: %f for %.1f*X %.1f*Y %.1f*Z" % (
            corr, f[0], f[1], f[2] )


    def ReportZScore ( self, fits ) :

        fit_scores = [c for c,M,regs,stats in fits]
        if ( len(fit_scores) > 3 ) :
            fit_scores.sort ()
            fit_scores.reverse ()
            best_score = fit_scores[0]
            other_scores = fit_scores[1:14]
            print "Best score: ", best_score
            print "Next scores: ", other_scores
            avg = numpy.average ( other_scores )
            stdev = numpy.std ( other_scores )
            self.zscore = (best_score - avg) / stdev
            umsg ( "Top score: %.5f, z-score: %.5f (avg: %.4f, stdev: %.4f)" % (
                best_score, self.zscore, avg, stdev) )



    def ClusterFits ( self, fits ) :

        class ClusterEntry :
            def __init__ (self, T, corr, regs, stats) :
                self.M = T
                xf = chimera.Xform.xform ( T[0,0], T[0,1], T[0,2], T[0,3], T[1,0], T[1,1], T[1,2], T[1,3], T[2,0], T[2,1], T[2,2], T[2,3] )
                self.COM = xf.getTranslation()
                self.Q = Segger.quaternion.Quaternion()
                self.Q.fromXform ( xf )
                self.corr = corr
                self.regs = regs
                self.stats = stats

        class Cluster :
            def __init__ ( self, e ) :
                self.entries = [ e ]
                self.COM = e.COM
                self.Q = e.Q
                self.corr = e.corr
                self.M = e.M
                self.regs = e.regs
                self.stats = {}
                self.stats['maxAngle'] = self.sumAngles = e.stats['totAngle']
                self.stats['maxShift'] = self.sumShifts = e.stats['totShift']
                self.stats['maxHeight'] = self.sumHeights = e.stats['difCC']
                self.stats['numFits'] = 1

            def AddEntry ( self, new_e ) :
                self.entries.append  ( new_e )
                if new_e.corr > self.corr :
                    self.corr = new_e.corr
                    self.M = new_e.M
                    self.regs = new_e.regs

                totAngle, totShift, difCC = new_e.stats['totAngle'], new_e.stats['totShift'], new_e.stats['difCC']
                if totAngle > self.stats['maxAngle'] : self.stats['maxAngle'] = totAngle
                if totShift > self.stats['maxShift'] : self.stats['maxShift'] = totShift
                if difCC > self.stats['maxHeight'] : self.stats['maxHeight'] = difCC

                self.stats['numFits'] = self.stats['numFits'] + 1

                # compute the new averages
                self.COM = chimera.Vector (0,0,0)
                self.Q = Segger.quaternion.Quaternion(0, chimera.Vector(0,0,0))

                for e in self.entries :
                    self.COM = self.COM + e.COM
                    self.Q = self.Q + e.Q
                    self.sumAngles = self.sumAngles + totAngle
                    self.sumShifts = self.sumShifts + totShift
                    self.sumHeights = self.sumHeights + difCC

                self.COM = self.COM / float ( len(self.entries) )
                self.Q.normalize()

                self.avgAngle = self.sumAngles / float ( len(self.entries) )
                self.avgShift = self.sumShifts / float ( len(self.entries) )
                self.avgHeight = self.sumHeights / float ( len(self.entries) )


            def SimilarTo ( self, e, posTol=5.0, angleTol=5.0 ) :
                if ( (self.COM - e.COM).length < posTol and
                     self.Q.angleTo ( e.Q ) * 180.0 / numpy.pi < angleTol ) :
                    return True
                return False

        posTol = float ( self.positionTolString.get() )
        angleTol = float ( self.angleTolString.get() )
        if self.doClusterFits.get () :
            print "Clustering %d fits..." % len(fits)
            print " - distance < ", posTol
            print " - angle < ", angleTol

        clusters = []
        for corr, M, regs, stats in fits :

            e = ClusterEntry ( M, corr, regs, stats )

            bAdded = False

            if self.doClusterFits.get () :
                for c in clusters :
                    if c.SimilarTo ( e, posTol, angleTol ) :
                        c.AddEntry ( e )
                        bAdded = True
                        break

            if bAdded == False :
                clusters.append ( Cluster (e) )

        print "%d clusters" % len(clusters)

        cfits = []
        for c in clusters :
            cfits.append ( [c.corr, c.M, c.regs, c.stats] )

        return cfits




    def saFitMapToPoints_byRot ( self, fmap, points, dmap, task=None, N=10, M=10 ) :

        print "fitting %s in map %s, to %d points, by rotation" % (fmap.name, dmap.name, len(points))

        num = float ( self.rotaSearchNum.get() )
        N = int ( numpy.floor ( numpy.sqrt ( num ) ) )
        M = int ( numpy.floor ( num / N ) )

        fpoints, fpoint_weights = fit_points ( fmap, (self.useLaplace.get()==False) )

        print "%d fits - rotations %d axes, %d angles" % (num, N, M)
        alist = uniform_rotation_angles(N, M)

        COM, U, S, V = prAxes ( points )
        comT = numpy.matrix ( [
            [ 1, 0, 0, COM[0] ],
            [ 0, 1, 0, COM[1] ],
            [ 0, 0, 1, COM[2] ],
            [ 0, 0, 0,      1 ]  ] )

        mlist = [comT*rotation_from_angles(*angles)*fmap.preM for angles in alist]

        from math import pi
        names = ['theta %.0f, phi %.0f, rot %.0f'
                 % tuple([a*180/pi for a in a3]) for a3 in alist]
        status_text = 'Rotational fit'

        optimize = self.optimize_fits.get()

        fits = optimize_fits(fpoints, fpoint_weights, mlist, dmap, names, status_text, optimize, False, task)
        corr, Mfit, i = self.make_best_fit(fits, fmap, dmap)

        print " - best fit: %f\n" % ( corr, )


    def make_best_fit(self, fits, fmap, dmap):

        i = numpy.argmax([c for Mf,c,stats in fits])
        Mfit, corr, stats = fits[i]       # highest correlation fit
        fmap.fit_score = corr
        fmap.M = Mfit

        list_fits = [(c, Mf, fmap.fit_regions, stats) for Mf,c,stats in fits]
        self.fits.extend(list_fits)

        move_fit_models(fmap, Mfit, dmap.openState.xform)

        return corr, Mfit, i


    def MoleculeMap ( self, create = True, warn = True ) :

        label = self.struc.get()
        sel_str = "#" + label [ label.rfind("(")+1 : label.rfind(")") ]
        try :
            fmol = chimera.selection.OSLSelection(sel_str).molecules()[0]
        except :
            umsg ( "%s not open - " % self.struc.get() ); return

        try :
            fmol.fitting_map.name
            return fmol.fitting_map
        except : pass
        if create:
            fmol.fitting_map = self.GenStrucMap(show = False)
            return fmol.fitting_map
        return None




    def GroupAroundReg ( self, smod, regs, target_volume, bRad=-1.0 ) :

        dv_rgroups = []
        maxDepthReached = 0

        stack = [([reg], reg.enclosed_volume(), 0) for reg in regs]

        status ( "Making groups around %d regions" % len(regs) )

        while len(stack) > 0 :

            regs, vol_at, depth_at = stack.pop()

            if depth_at > maxDepthReached : maxDepthReached = depth_at

            if depth_at >= SAF_LS_DEPTH : continue
            if vol_at >= target_volume * (1.0 + SAF_DVOL)  : continue

            dv = abs ( vol_at - target_volume ) / target_volume
            dv_rgroups.append ( [dv, regs] )

            if len(dv_rgroups) % 100 == 0 :
                status ( "Making groups around %d - %d groups" % (reg.rid, len(dv_rgroups)) ),

            reg_at = regs[0]
            for cr in reg_at.contacting_regions():
                if regs.count ( cr ) != 0 : continue
                if cr.placed : continue
                vol = vol_at + cr.enclosed_volume()
                stack.insert ( 0, [ [cr]+regs, vol, depth_at+1 ] )

        dv_rgroups_f = self.FilterGroups ( dv_rgroups, bRad )

        print " - %d groups --> %d filtered groups" % (len(dv_rgroups), len(dv_rgroups_f))
        return [dv_rgroups_f, maxDepthReached]





    def FilterGroups ( self, dv_rgroups, bRad = -1.0 ) :

        dv_rgroups_f = []
        len_groups = {}
        inc_regs_map = {}

        gi, ngroups = 0, len(dv_rgroups)

        for dv, regs in dv_rgroups :

            gi = gi + 1
            if gi % 100 == 0 :
                status ( "Filtering group %d/%d" % (gi,ngroups) )

            if dv > SAF_DVOL : continue

            included = False

            regs.sort()
            inc_regs_at = inc_regs_map
            for reg in regs :
                try : included, inc_regs_at = inc_regs_at[reg]
                except : included = False; break

            if included : continue

            if bRad > 0.0 :
                regs_bRad = regions_radius(regs)

                brad_d = abs ( regs_bRad - bRad ) / bRad

                if brad_d > SAF_DBRAD : continue

            dv_rgroups_f.append ( [dv, regs] )

            inc_regs_at = inc_regs_map
            last_arr = None
            for reg in regs :
                try :
                    last_arr = inc_regs_at[reg]
                    inc_regs_at = last_arr[1]
                except :
                    last_arr = [ False, {} ]
                    inc_regs_at[reg] = last_arr
                    inc_regs_at = last_arr[1]

            last_arr[0] = True


        return dv_rgroups_f




    def GroupAllRegions ( self, smod, target_volume, bRad=-1.0) :

        dv_rgroups = []
        maxDepthReached = 0

        print "Grouping %d regions in %s, target volume %.2f, bounding radius %.2f" % (
            len(smod.surfacePieces), smod.name, target_volume, bRad )

        ri, nregs = 0, len(smod.regions)

        for reg in smod.regions :

            ri = ri + 1

            dv_rgroupsR, maxDepthReachedR = self.GroupAroundReg ( smod, [reg], target_volume, bRad )

            dv_rgroups = dv_rgroups + dv_rgroupsR

            if maxDepthReachedR > maxDepthReached : maxDepthReached = maxDepthReachedR


        print "\n - max depth reached: %d" % maxDepthReached

        print " - filtering %d groups..." % len( dv_rgroups )
        return self.FilterGroups ( dv_rgroups, bRad )



    def FitMapToSelRGroup ( self, dmap = None, task = None ) :

        if dmap is None:
            dmap = segmentation_map()
            if dmap == None : print "No segmentation map"; return

        fmap = None
        if self.map_to_fit : fmap = self.map_to_fit
        else : fmap = self.MoleculeMap()

        thr = fmap.surface_levels[0]
        mm = fmap.data.matrix()
        mmab = numpy.where ( mm > thr, numpy.ones_like(mm), numpy.zeros_like(mm) )
        nz = numpy.shape ( numpy.nonzero ( mmab ) )[1]
        tvol = float(nz) * fmap.data.step[0] * fmap.data.step[1] * fmap.data.step[2]
        print "%s - %d above %f, volume %.3f" % (fmap.name, nz, thr, tvol)

        smod = self.CurrentSegmentation()
        if smod == None : return

        regs = smod.selected_regions()
        if len(regs)==0 : umsg ( "Please select a region to fit to" ); return

        reportFitRegions(fmap.name, regs)

        for reg in regs:
            clr = reg.color
            reg.surface_piece.color = ( clr[0], clr[1], clr[2], REG_OPACITY )
        print ""

        points = numpy.concatenate ( [r.map_points()
                                      for r in regs], axis=0 )
        regions_vol = sum([r.enclosed_volume() for r in regs])

        dv = abs(regions_vol - tvol) / tvol
        print " - regions volume: %.2f - dv %.5f" % (regions_vol, dv)

        self.fits = []
        fmap.fit_regions = regs

        to_map = dmap
        reg_map = None
        if 1 :
            reg_map = mask_volume( regs, dmap )


        bCloseMap = False

        if self.useLaplace.get () :
            umsg ("Generating Laplace version of " + dmap.name)
            from VolumeFilter import laplacian
            to_map = laplacian ( dmap )
            bCloseMap = True

        elif self.mask_map_when_fitting.get() :
            to_map = reg_map
            bCloseMap = True
            if to_map is None:
                umsg ('Could not create masked map')
                return

        if self.rotaSearch.get () :
            self.saFitMapToPoints_byRot ( fmap, points, to_map, task )
        else :
            self.saFitMapToPoints ( fmap, points, to_map, task )

        self.cfits = self.ClusterFits ( self.fits )
        self.cfits.sort ( reverse=True, key=lambda x: x[0] )
        #cfits.sort()
        #cfits.reverse()

        try : nToAdd = int ( self.numFitsToAdd.get () )
        except : nToAdd = len (self.cfits)
        for corr, M, regions, stats in self.cfits [ 0 : nToAdd ] :
            print " -- clustered Fit -- #fits:%d maxAngle:%.1f maxShift:%.2f maxHeight:%.2f" % ( stats['numFits'], stats['maxAngle'], stats['maxShift'], stats['maxHeight'] )
            fmap.fit_score, fmap.M, fmap.fit_regions = corr, M, regions
            fmap.atomInclusion, fmap.bbAtomInclusion, fmap.bbClashes, fmap.hdoScore  = self.FitScores ( fmap, reg_map )
            self.add_fit (fmap, dmap)

        move_fit_models(fmap, self.cfits[0][1], dmap.openState.xform)

        #umsg ( "Cross-correlation: %.4f\n" % (fmap.fit_score) )
        self.ReportZScore ( self.cfits )

        # close masked map if it was created
        #if bCloseMap : chimera.openModels.close ( [to_map] )
        if reg_map :
            chimera.openModels.close ( [reg_map] )

        for m in fmap.mols :
            m.display = True





    def MapIndexesInMap ( self, ref_map, mask_map ) :

        thr = mask_map.surface_levels[0]
        mm = mask_map.data.matrix()
        mm = numpy.where ( mm > thr, mm, numpy.zeros_like(mm) )

        nze = numpy.nonzero ( mm )

        # copy is needed! transform_vertices requires contiguous array
        points =  numpy.empty ( (len(nze[0]), 3), numpy.float32)
        points[:,0] = nze[2]
        points[:,1] = nze[1]
        points[:,2] = nze[0]

        print "Making map indices for %s in %s" % ( mask_map.name, ref_map.name )
        print " - %d points above %.3f" % ( len(points), thr )

        # transform to index reference frame of ref_map
        f1 = mask_map.data.ijk_to_xyz_transform
        f2 = xform_matrix ( mask_map.openState.xform )
        f3 = xform_matrix ( ref_map.openState.xform.inverse() )
        f4 = ref_map.data.xyz_to_ijk_transform

        tf = multiply_matrices( f2, f1 )
        tf = multiply_matrices( f3, tf )
        tf = multiply_matrices( f4, tf )
        transform_vertices ( points, tf )

        imap = set()
        for fi, fj, fk in points :
            for i in [ int(numpy.floor(fi)), int(numpy.ceil(fi)) ] :
                for j in [ int(numpy.floor(fj)), int(numpy.ceil(fj)) ] :
                    for k in [ int(numpy.floor(fk)), int(numpy.ceil(fk)) ] :
                        imap.add((i,j,k))

        return imap


    def ZeroMatWitMap ( self, ref_mat, ref_map, mask_map ) :

        thr = mask_map.surface_levels[0]
        mm = mask_map.data.matrix()
        mm = numpy.where ( mm > thr, mm, numpy.zeros_like(mm) )

        nze = numpy.nonzero ( mm )

        # copy is needed! transform_vertices requires contiguous array
        points =  numpy.empty ( (len(nze[0]), 3), numpy.float32)
        points[:,0] = nze[2]
        points[:,1] = nze[1]
        points[:,2] = nze[0]

        print "Making map indices for %s in %s" % ( mask_map.name, ref_map.name )
        print " - %d points above %.3f" % ( len(points), thr )

        # transform to index reference frame of ref_map
        f1 = mask_map.data.ijk_to_xyz_transform
        f2 = xform_matrix ( mask_map.openState.xform )
        f3 = xform_matrix ( ref_map.openState.xform.inverse() )
        f4 = ref_map.data.xyz_to_ijk_transform

        tf = multiply_matrices( f2, f1 )
        tf = multiply_matrices( f3, tf )
        tf = multiply_matrices( f4, tf )
        transform_vertices ( points, tf )

        imap = set()
        for fi, fj, fk in points :
            for i in [ int(numpy.floor(fi)), int(numpy.ceil(fi)) ] :
                for j in [ int(numpy.floor(fj)), int(numpy.ceil(fj)) ] :
                    for k in [ int(numpy.floor(fk)), int(numpy.ceil(fk)) ] :
                        #imap.add((i,j,k))
                        try :
                            ref_mat[k,j,i] = 0
                        except :
                            pass



    def ZeroMatWitMol ( self, M, dmap, mol ) :

        import _multiscale
        points = _multiscale.get_atom_coordinates ( mol.atoms, transformed = False )

        # transform to index reference frame of ref_map
        f1 = xform_matrix ( mol.openState.xform )
        f2 = xform_matrix ( dmap.openState.xform.inverse() )
        f3 = dmap.data.xyz_to_ijk_transform

        tf = multiply_matrices( f2, f1 )
        tf = multiply_matrices( f3, tf )

        transform_vertices ( points, tf )

        print " - ", len(points), "points"

        l = 3.0 / dmap.data.step[0]
        l2 = l*l
        print l

        nz = 0

        if 0 :
            for fi, fj, fk in points :
                for i in [ int(numpy.floor(fi-l)), int(numpy.ceil(fi+l))+1 ] :
                    for j in [ int(numpy.floor(fj-l)), int(numpy.ceil(fj+l))+1 ] :
                        for k in [ int(numpy.floor(fk-l)), int(numpy.ceil(fk+l))+1 ] :
                            di, dj, dk = i-fi, j-fj, k-fk
                            d = di*di + dj*dj + dk*dk
                            if d < l2 :
                                try :
                                    print k, j, i
                                    M[k,j,i] = 0
                                    nz += 1
                                    if nz > 20 :
                                        return M
                                except :
                                    pass
        for fi, fj, fk in points :
            for i in range ( int(numpy.floor(fi-l)), int(numpy.ceil(fi+l))+1 ) :
                for j in range ( int(numpy.floor(fj-l)), int(numpy.ceil(fj+l))+1 ) :
                    for k in range ( int(numpy.floor(fk-l)), int(numpy.ceil(fk+l))+1 ) :
                        di, dj, dk = i-fi, j-fj, k-fk
                        d = di*di + dj*dj + dk*dk
                        if d < l2 :
                            try :
                                M[k,j,i] = 0
                                nz += 1
                            except :
                                pass

        print nz
        return M






    def OverlappingRegions ( self, dmap, fmap, smod, hide_others = True ) :

        imap = self.MapIndexesInMap ( dmap, fmap )
        print 'imap', len(imap)

        try : fmap.COM
        except : fmap.COM = fmap.mols[0].COM; fmap.bRad = fmap.mols[0].BoundRad

        p = numpy.array ( [ fmap.COM ], numpy.float32 )
        transform_vertices ( p, xform_matrix( fmap.openState.xform ) )
        transform_vertices ( p, xform_matrix( dmap.openState.xform.inverse() ) )
        f_COM = chimera.Vector ( *p[0] )
        f_bRad = fmap.bRad
        print " center", f_COM, "brad", f_bRad

        jregs = []

        for r in smod.regions :

            if r.placed:
                continue

            ipoints = r.points()
            noverlap = 0
            for i,j,k in ipoints :
                if (i,j,k) in imap:
                    noverlap += 1

            ov = float(noverlap) / len(ipoints)

            if ov > .8 :
                jregs.append ( r )

            try :
                if ov > r.max_ov :
                    r.max_ov = ov
                    r.max_ov_cid = fmap.chain_id
                    r.max_ov_bioM = smod.bio_mt_at
            except :
                pass


        oregs = jregs

        umsg ( "Found %d regions overlapping" % len(oregs) )
        sel_sps = []

        for sp in smod.surfacePieces :
            try : clr = sp.region.color
            except : continue
            if oregs.count ( sp.region ) > 0 :
                sp.color = ( clr[0], clr[1], clr[2], REG_OPACITY )
                sp.display = True
                sel_sps.append ( sp )
            else :
                sp.color = ( clr[0], clr[1], clr[2], 1.0 )
                sp.display = not hide_others

        return jregs





    def ShowOverlappingRegions ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        smod = self.CurrentSegmentation()
        if smod is None : return

        fmap = self.MoleculeMap()
        if fmap == None : return


        oregs = self.OverlappingRegions ( dmap, fmap, smod )

        return oregs




    def FitMapToGroupsAround ( self, fmap, smod, regs, dmap, bFirst = True ) :

        bestFitScore = -1e99
        bestFitM = None
        bestFitGroup = None
        bestFitRegions = None
        fmap.fit_score = None

        tvol = self.MapVolume ( fmap )
        bRad = -1.0 # fmap.mol.BoundRad / float(dmap.data.step[0]); # self.MapBoundingRad ( fmap )
        bRad = fmap.mols[0].BoundRad

        print "\nMaking groups around %d regions - target vol %.3f, b-Rad %.3f" % ( len(regs), tvol, bRad )
        smod.rgroups, maxDepthReached = self.GroupAroundReg ( smod, regs, tvol, bRad )
        smod.rgroups.sort()
        print " - depth reached: %d" % maxDepthReached

        if len(smod.rgroups) == 0 : umsg ( "No groups found" ); return -1e99

        nsearchgrps = min ( len(smod.rgroups), SAF_LS_NGROUPS )
        print "________________________________________________________________________"
        umsg ( "Fitting %s to %d/%d groups..." % ( fmap.name, nsearchgrps, len(smod.rgroups) ) )
        print "________________________________________________________________________"


        self.fits = []

        for i, dv_regs in enumerate ( smod.rgroups[0:nsearchgrps] ) :

            dv, regs = dv_regs

            umsg ( "Fitting to group %d/%d, dVolume %.4f, %d regions" % (i+1, nsearchgrps, dv, len(regs) ) )
            print " - regions:",
            for r in regs : print r.rid,
            print ""

            fmap.fit_regions = regs

            points = numpy.concatenate ( [r.map_points() for r in regs], axis=0 )

            if self.rotaSearch.get () :
                self.saFitMapToPoints_byRot ( fmap, points, dmap, task )
            else :
                self.saFitMapToPoints ( fmap, points, dmap, task )

            print ""

            if fmap.fit_score > bestFitScore :
                bestFitScore = fmap.fit_score
                bestFitM = fmap.M
                bestFitRegions = regs

        umsg ( "Best cross-correlation: %.4f\n\n" % ( bestFitScore ) )

        fmap.fit_score = bestFitScore
        fmap.M = bestFitM
        fmap.fit_regions = bestFitRegions

        tXO, tXR = xf_2_M ( dmap.openState.xform )
        T = tXO * tXR * fmap.M
        xfA = chimera.Xform.xform ( T[0,0], T[0,1], T[0,2], T[0,3], T[1,0], T[1,1], T[1,2], T[1,3], T[2,0], T[2,1], T[2,2], T[2,3] )
        fmap.openState.xform = xfA
        for mol in fmap.mols : mol.openState.xform = xfA

        self.cfits = self.ClusterFits ( self.fits )
        #cfits.sort ( reverse=True, key=lambda x: x[0] )
        self.cfits.sort()
        self.cfits.reverse()

        try : nToAdd = int ( self.numFitsToAdd.get () )
        except : nToAdd = len (self.cfits)

        for corr, M, regions, stats in self.cfits [ 0 : nToAdd ] :
            fmap.fit_score, fmap.M, fmap.fit_regions = corr, M, regions
            self.add_fit (fmap, dmap)
            # TODO - add atom inclusion comp

        #umsg ( "Cross-correlation: %.4f\n" % (fmap.fit_score) )
        self.ReportZScore ( self.cfits )


    def FitMapToRegionsAroundSel ( self, task=None ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        smod = self.CurrentSegmentation()
        if smod is None : return

        regs = smod.selected_regions()

        fmap = self.MoleculeMap()
        if fmap == None : return

        if timing: t0 = clock()
        self.FitMapToGroupsAround ( fmap, smod, regs, dmap )
        if timing:
            t1 = clock()
            print "Time: %.1f sec" % (t1-t0,)

        if fmap.fit_score is None:
            umsg('No groups of regions meet size requirement')
            return

        oregs = self.OverlappingRegions ( dmap, fmap, smod, hide_others = False )

        if len(oregs) == 1 :
            oregs[0].placed = True


    def FitOpenMapsToGroupsAround ( self, smod, reg, dmap, bFirst=True ) :

        bestFitScore = -1e99
        bestFitMap = None
        bestFitM = None

        fmaps = []

        # TODO: Don't use "centered" to decide what maps to fit.
        for fmap in OML() :
            try : fmap.mols[0].centered
            except : continue
            #fmap.display = False
            for mol in fmap.mols : mol.display = False
            fmaps.append ( fmap )


        for fmap in fmaps :

            print "\n****************************************************"
            print "Fitting: ", fmap.name
            print "****************************************************\n"

            for mol in fmap.mols : mol.display = True

            self.FitMapToGroupsAround ( fmap, smod, reg, dmap, bFirst )

            if fmap.fit_score == None :
                print " - no fits for map", fmap.name

            elif fmap.fit_score > bestFitScore :
                bestFitScore = fmap.fit_score
                bestFitMap = fmap
                bestFitM = fmap.M

            for mol in fmap.mols : mol.display = False


        if bestFitM == None :
            print "No best fit recorded, perhaps there were not groups to start with"
            return None

        fmap = bestFitMap
        for mol in fmap.mols : mol.display = True

        print "\n**********************************************************"
        print "Best fit score was %.4f for %s" % (bestFitScore, fmap.name)
        print "**********************************************************\n"
        fmap.M = bestFitM

        tXO, tXR = xf_2_M ( dmap.openState.xform )
        T = tXO * tXR * fmap.M
        xfA = chimera.Xform.xform ( T[0,0], T[0,1], T[0,2], T[0,3], T[1,0], T[1,1], T[1,2], T[1,3], T[2,0], T[2,1], T[2,2], T[2,3] )
        fmap.openState.xform = xfA
        for mol in fmap.mols : mol.openState.xform = xfA

        oregs = self.OverlappingRegions ( dmap, fmap, smod, hide_others = False  )

        if oregs.count ( reg ) == 0 :
            print "Overlapping regions not inclusive"
            return False


        if len(oregs) > 1 :
            jreg = smod.join_regions ( oregs )
            jreg.placed = True
            self.ReportRegionCount(smod)

            jsp = jreg.surface_piece
            clr = jreg.color
            jsp.color = ( clr[0], clr[1], clr[2], REG_OPACITY )
            jsp.display = False

        elif len(oregs) == 1 :
            oregs[0].placed = True
            oregs[0].surface_piece.display = False

        else :
            for mol in fmap.mols : mol.display = False
            return False

        for mol in fmap.mols : mol.display = False

        self.add_fit(fmap, dmap)

        return True



    def FitMapsToRGroups ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        smod = self.CurrentSegmentation()
        if smod == None : return

        for sp in smod.surfacePieces :
            clr = sp.region.color
            sp.color = ( clr[0], clr[1], clr[2], REG_OPACITY )
            #sp.display = False


        if timing: t0 = clock()

        while 1 :

            sp = None
            for spi in smod.surfacePieces :

                try :
                    spi.region.placed.display = False

                    clr = spi.region.color
                    spi.color = ( clr[0], clr[1], clr[2], REG_OPACITY )
                    spi.display = True
                    spi.region.placed.display = True
                    spi.display = False

                except :
                    pass

                if spi.region.placed == False :
                    sp = spi
                    break

            if sp == None :
                print "\nAll regions have maps placed"
                break

            clr = sp.region.color
            sp.color = ( clr[0], clr[1], clr[2], REG_OPACITY )
            sp.display = True

            if self.FitOpenMapsToGroupsAround ( smod, sp.region, dmap, False ) == False :
                print "___ No map fit for region! ___", sp.region.rid, sp.region.enclosed_volume()
                sp.failed_fit = True

            sp.region.placed = True


        if timing:
            t1 = clock()
            print "Time: %.1f sec" % (t1-t0)

        for sp in smod.surfacePieces :
            try : sp.failed_fit
            except : continue
            print "Region %d failed fit and still in model" % sp.region.rid



    def MapVolume ( self, fmap ) :

        thr = fmap.surface_levels[0]
        mm = fmap.data.matrix()
        mmab = numpy.where ( mm > thr, numpy.ones_like(mm), numpy.zeros_like(mm) )
        nz = numpy.shape ( numpy.nonzero ( mmab ) )[1]
        vvol = fmap.data.step[0] * fmap.data.step[1] * fmap.data.step[2]
        tvol = vvol * float(nz)
        print "%s - %d above %f, VOLUME %.3f" % (fmap.name, nz, thr, tvol)
        return tvol


    def Scores ( self ) :

        fmap = self.MoleculeMap()
        if fmap == None : return

        self.FitScores ( fmap )


    def SMS ( self ) :

        dmap = segmentation_map()
        if dmap == None :
            print "No segmentation map";
            return

        fmap = self.MoleculeMap()
        if fmap == None : return

        sel_str = "#%d@C,N,CA" % fmap.mols[0].id
        sel = chimera.selection.OSLSelection (sel_str)
        backbone_atoms = sel.atoms()

        sms = ShapeMatchScore ( backbone_atoms, dmap )



    def VisiScores ( self ) :

        print "Visi scores..."
        dmap = segmentation_map()
        if dmap == None :
            return

        print " - in map: " + dmap.name

        molmap = None
        mols = []

        for m in chimera.openModels.list() :
            if m.display == False :
                continue
            if type(m) == chimera.Molecule :
                print " - mol: ", m.name
                mols.append ( m )
            if type(m) == VolumeViewer.Volume :
                print " - map: ", m.name
                molmap = m

        if molmap != None :
            fpoints, fpoint_weights = fit_points(molmap, True)
            map_values = dmap.interpolated_values ( fpoints, molmap.openState.xform )
            olap, corr = overlap_and_correlation ( fpoint_weights, map_values )
            print " - Overlap: %f, Cross-correlation: %f" % (olap, corr)



        otherPoints = None
        otherAtoms = []


        for m in mols :
            otherAtoms = otherAtoms + m.atoms
            mpoints = get_atom_coordinates ( m.atoms, transformed = True )
            if otherPoints == None : otherPoints = mpoints
            else : otherPoints = numpy.concatenate ( [otherPoints, mpoints], axis=0 )

        # print "Doing tree with %d %d" % ( len(otherPoints), len(otherAtoms) )

        print " - making tree, %d atoms" % len(otherAtoms)

        from CGLutil.AdaptiveTree import AdaptiveTree
        searchTreeAll = AdaptiveTree (otherPoints.tolist(), otherAtoms, 4.0)

        print " - checking clashes, %d atoms" % len(otherAtoms)

        numClash = 0.0
        for at in otherAtoms :
            nearby = searchTreeAll.searchTree ( at.xformCoord().data(), 3.0 )
            for nb in nearby :
                if nb.molecule != at.molecule :
                    numClash = numClash + 1.0
                    break

        bbClashes = numClash / float ( len(otherAtoms) )

        print " - clashes: %.0f/%.0f = %0.3f clash-free" % (numClash, len(otherAtoms), 1.0-bbClashes);






    def FitScores ( self, fmap, regionMap = None ) :

        dmap = segmentation_map()
        if dmap == None :
            print "No segmentation map";
            return [0.0, 0.0, 0.0, 0.0]

        print "Fit scores for", fmap.name, "in", dmap.name

        import numpy
        import _contour

        # move fmap (and structures) to fmap.M if it's there
        # (it's not there if we do scores for a selected structure
        # that hasn't been fit yet)
        try : fmap.M; move = True
        except : move = False

        if move :
            tf = numpy.array(fmap.M)
            xf = dmap.openState.xform
            xf.multiply(chimera_xform(tf))
            fmap.openState.xform = xf
            for mol in fmap.mols :
                mol.openState.xform = xf


        # ---------------------------------------------------------------
        # Cross-correlation
        # ---------------------------------------------------------------
        fpoints, fpoint_weights = fit_points(fmap)
        map_values = dmap.interpolated_values ( fpoints, fmap.openState.xform )
        olap, corr = overlap_and_correlation ( fpoint_weights, map_values )
        print " - Overlap: %f, Cross-correlation: %f" % (olap, corr)


        # ---------------------------------------------------------------
        # By-residue cross-correlation
        # ---------------------------------------------------------------
        if 0 :
            cc_by_residue ( fmap, dmap, 16 )


        # ---------------------------------------------------------------
        # Atom inclusion -- all atoms
        # ---------------------------------------------------------------
        backbone_atoms = []
        allIncl = 0.0
        bbIncl = 0.0
        bbClashes = 0.0

        all_atoms = []
        for mol in fmap.mols : all_atoms = all_atoms + mol.atoms
        numAllAtoms = float ( len(all_atoms) )

        if len(all_atoms) == 0 :
            return [0.0, 0.0, 0.0, 0.0]

        #dmapXfInv = xform_matrix( dmap.openState.xform.inverse() )
        #transform_vertices( points, dmapXfInv )
        points = get_atom_coordinates ( all_atoms, transformed = True )
        dvals = dmap.interpolated_values ( points, chimera.Xform() )
        min_d = dmap.surface_levels[0]
        dvals = numpy.where ( dvals > min_d, dvals, numpy.zeros_like(dvals) )
        nze = numpy.nonzero ( dvals )
        allIn = float(len(nze[0]))
        allIncl = allIn / numAllAtoms

        print " - Atom inclusion: %.0f/%.0f = %.3f" % ( allIn, numAllAtoms, allIncl )

        # ---------------------------------------------------------------
        # Atom inclusion -- backbone atoms
        # ---------------------------------------------------------------
        bbIncl = 0.0

        sel_str = "#%d@C,N,CA" % fmap.mols[0].id
        sel = chimera.selection.OSLSelection (sel_str)
        backbone_atoms = sel.atoms()

        if len(backbone_atoms) > 0 :

            points = get_atom_coordinates ( backbone_atoms, transformed = True )
            #dmapXfInv = xform_matrix( dmap.openState.xform.inverse() )
            #transform_vertices( points, dmapXfInv )
            dvals = dmap.interpolated_values ( points, chimera.Xform() )
            min_d = dmap.surface_levels[0]
            dvals = numpy.where ( dvals > min_d, dvals, numpy.zeros_like(dvals) )
            nze = numpy.nonzero ( dvals )
            bbIn = float(len(nze[0]))
            numBBAtoms = float(len(backbone_atoms))
            bbIncl = bbIn / numBBAtoms
            print " - BB Atom inclusion: %.0f/%.0f = %.3f" % (bbIn, numBBAtoms, bbIncl );


        if 0 :
            sms = ShapeMatchScore ( backbone_atoms, dmap )


        # ---------------------------------------------------------------
        # Coverage of high-density areas - Density Occupancy
        # ---------------------------------------------------------------

        hdo = 0.0
        if regionMap :
            #fpoints, fpoint_weights = fit_points ( regionMap )
            #nz_fpoints = len(fpoints)

            nz_fpoints = len ( numpy.nonzero ( regionMap.data.full_matrix() )[0] )

            _contour.affine_transform_vertices ( points, Matrix.xform_matrix( regionMap.openState.xform.inverse() ) )
            s = regionMap.data.step[0]
            mdata = VolumeData.zone_masked_grid_data ( regionMap.data, points, numpy.sqrt(3*s*s) )
            #gv = VolumeViewer.volume.volume_from_grid_data ( mdata )
            #gv.openState.xform = dmap.openState.xform
            #gv.name = "Masked"

            mat = mdata.full_matrix()
            nz_mdata = len ( numpy.nonzero ( mat )[0] )
            if nz_fpoints > 0 :
                hdo = float (nz_mdata) / float(nz_fpoints)
                print " - Density Occupancy: %d / %d grid points above %.3f occupied (%.4f)" % (
                    nz_mdata, nz_fpoints, dmap.surface_levels[0], hdo )
        else :
            print " - not computing density occupancy"


        # ---------------------------------------------------------------
        # Clashes with symmetric copies
        # ---------------------------------------------------------------
        if self.calcSymmetryClashes.get() :

          symMols = self.PlaceSym ()
          if symMols:

            otherPoints = None
            otherAtoms = []

            for m in symMols :
                otherAtoms = otherAtoms + m.atoms
                mpoints = get_atom_coordinates ( m.atoms, transformed = True )
                if otherPoints == None : otherPoints = mpoints
                else : otherPoints = numpy.concatenate ( [otherPoints, mpoints], axis=0 )

            # print "Doing tree with %d %d" % ( len(otherPoints), len(otherAtoms) )

            from CGLutil.AdaptiveTree import AdaptiveTree
            searchTreeAll = AdaptiveTree (otherPoints.tolist(), otherAtoms, 4.0)

            #if len ( backbone_atoms ) == 0 :
            #    sel_str = "#%d@C,N,CA" % fmap.mols[0].id
            #    sel = chimera.selection.OSLSelection (sel_str)
            #    backbone_atoms = sel.atoms()

            numClash = 0.0
            for at in all_atoms :
                nearby = searchTreeAll.searchTree ( at.xformCoord().data(), 3.0 )
                if len(nearby) > 0 :
                    numClash = numClash + 1.0

            bbClashes = numClash / numAllAtoms

            chimera.openModels.close ( symMols )

            print " - Clashes with symmetric copies: %.0f/%.0f = %0.3f" % (numClash, numAllAtoms, bbClashes);

        print fmap.name, corr, allIncl, bbClashes, hdo

        #for i in range ( 1 ) :
        #    print (frame_at+i),
        #    chimera.printer.saveImage ( "./frames/%06d.png" % (frame_at + i) )
        #print ""

        return [allIncl, bbIncl, bbClashes, hdo]


    def FitMapToRGroups ( self, task=None ) :

        print "_______________________________________________________________"

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        fmap = self.MoleculeMap()
        if fmap == None : return
        tvol = self.MapVolume ( fmap )

        smod = self.CurrentSegmentation()
        if smod is None : return

        print "---"

        if timing: t0 = clock()

        bRad = fmap.mols[0].BoundRad

        smod.rgroups = self.GroupAllRegions ( smod, tvol, bRad )
        smod.rgroups.sort()

        bestFitScore = -1e99
        bestFitM = None

        print "Got %d groups..." % (len(smod.rgroups) )

        dmap_name = os.path.splitext ( dmap.name )[0]
        path = os.path.dirname ( dmap.data.path ) + os.path.sep

        nsearchgrps = min ( MAX_NUM_GROUPS, len(smod.rgroups) )
        if nsearchgrps == 0:
            umsg('No groups of regions meet size requirement')
            return

        self.fits = []

        for i, dv_regs in enumerate ( smod.rgroups [0:nsearchgrps] ) :

            dv, regs = dv_regs

            umsg ( "Fitting to group %d/%d, dVolume %.4f, %d regions" % (i+1, nsearchgrps, dv, len(regs) ) )
            print " - regions:",
            for r in regs : print r.rid,
            print ""

            for sp in smod.surfacePieces :
                if regs.count ( sp.region ) > 0 :
                    sp.display = True
                    clr = sp.region.color
                    sp.color = ( clr[0], clr[1], clr[2], REG_OPACITY )

                else : sp.display = False

            fmap.fit_regions = regs

            # TODO: points need to be in dmap coordinates.
            points = numpy.concatenate ( [r.map_points()
                                          for r in regs], axis=0 )
            if self.rotaSearch.get () :
                self.saFitMapToPoints_byRot ( fmap, points, dmap )
            else :
                self.saFitMapToPoints ( fmap, points, dmap )

            if fmap.fit_score > bestFitScore :
                bestFitScore = fmap.fit_score
                bestFitM = fmap.M
                bestFitRegs = regs

        umsg ( "Best cross-correlation: %.4f\n\n" % ( bestFitScore ) )

        fmap.fit_score = bestFitScore
        fmap.M = bestFitM

        tXO, tXR = xf_2_M ( dmap.openState.xform )
        T = tXO * tXR * fmap.M
        xfA = chimera.Xform.xform ( T[0,0], T[0,1], T[0,2], T[0,3], T[1,0], T[1,1], T[1,2], T[1,3], T[2,0], T[2,1], T[2,2], T[2,3] )
        fmap.openState.xform = xfA
        for mol in fmap.mols : mol.openState.xform = xfA

        if timing:
            t1 = clock()
            print "Time: %.1f sec" % (t1-t0)

        oregs = self.OverlappingRegions ( dmap, fmap, smod, hide_others = False  )

        self.cfits = self.ClusterFits ( self.fits )
        self.cfits.sort ( reverse=True, key=lambda x: x[0] )
        #cfits.sort()
        #cfits.reverse()

        try : nToAdd = int ( self.numFitsToAdd.get () )
        except : nToAdd = len (self.cfits)

        for corr, M, regions, stats in self.cfits [ 0 : nToAdd ] :
            fmap.fit_score, fmap.M, fmap.fit_regions = corr, M, regions
            self.add_fit (fmap, dmap)
            # TODO - add atom inclusion comp

        #umsg ( "Cross-correlation: %.4f\n" % (fmap.fit_score) )
        self.ReportZScore ( self.cfits )




    def GetMapFromMolRes ( self, mol, cid, rStart, rEnd ) :

        sel_str = "#%d:%d-%d.%s" % (mol.id, rStart, rEnd, cid)
        print "[%s]" % (sel_str),


        res = float ( self.simRes.get() )
        grid = float ( self.simGridSp.get() )

        cmd = "molmap %s %f sigmaFactor 0.187 gridSpacing %f replace false" % ( sel_str, res, grid )
        chimera.runCommand ( cmd )

        mv = None
        for mod in chimera.openModels.list() :
            ts = mod.name.split()
            if len(ts) > 1 and mod.name.find("map") >=0 and mod.name.find("res") >=0 :
                #print " - found", mod.name
                mv = mod
                mv.name = "_" + sel_str
                break

        if mv == None :
            umsg (" - error - could not find chain map")

        return mv



    def GetMapFromMolRanges ( self, mol, cid, ranges ) :

        sel_str = "#%d:%s" % (mol.id, ranges)
        print "[%s]" % (sel_str),


        res = float ( self.simRes.get() )
        grid = float ( self.simGridSp.get() )

        cmd = "molmap %s %f sigmaFactor 0.187 gridSpacing %f replace false" % ( sel_str, res, grid )
        chimera.runCommand ( cmd )

        mv = None
        for mod in chimera.openModels.list() :
            ts = mod.name.split()
            if len(ts) > 1 and mod.name.find("map") >=0 and mod.name.find("res") >=0 :
                #print " - found", mod.name
                mv = mod
                mv.name = "_" + sel_str
                break

        if mv == None :
            umsg (" - error - could not find chain map")

        return mv



    def GroupRegionsBySS ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        #fmap = self.MoleculeMap()
        #if fmap == None : return

        smod = self.CurrentSegmentation()
        if smod is None : return

        print "---"

        #mol = fmap.mols[0]
        #print mol.name

        #chain_colors = RandColorChains ( mol )

        chain_maps = []

        res = float ( self.simRes.get() )
        grid = float ( self.simGridSp.get() )

        print "_____________ res %2f _______ grid %.2f _________________________________" % (res, grid)


        for mol in chimera.openModels.list() :

            if type(mol) != chimera.Molecule or mol.display == False : continue

            basename = os.path.splitext ( mol.name )[0]
            #chain_colors = RandColorChains ( mol )


            chainsRes = {}
            for res in mol.residues :
                try :
                    chainsRes[res.id.chainId].append ( res )
                except :
                    chainsRes[res.id.chainId] = [res]

            chainsList = chainsRes.keys()
            chainsList.sort()


            for chainId in chainsList :

                residues = chainsRes[chainId]

                print " - chain " + chainId + ", %d " % len(residues) + " residues"

                ss, rStart = "", 0
                rI = 0

                oRanges = ""

                while 1 :
                    res = residues[rI]

                    if rStart == 0 :
                        print "  - at first res %d " % rI + ", pos: %d " % res.id.position,

                        rStart = res.id.position
                        if res.isHelix :
                            print " - H"
                            ss = "H"
                        else :
                            print ""
                            ss = ""

                    else :
                        #print "  - at res %d " % rI + ", pos: %d " % res.id.position,

                        if res.isHelix :
                            #print " - H "
                            if ss != "H" :
                                print "  - _->H - at res %d " % rI + ", pos: %d " % res.id.position
                                #mv = self.GetMapFromMolRes ( mol, chainId, rStart, res.id.position-1 )
                                #chain_maps.append ( [mv, self.MapIndexesInMap ( dmap, mv )] )
                                #mv.chain_id = basename + "_" + chainId + "_H%d" % rStart
                                if len(oRanges) > 0 : oRanges = oRanges + ","
                                oRanges = oRanges + "%d-%d.%s" % (rStart, res.id.position-1,chainId)
                                rStart = res.id.position
                                ss = "H"
                        else :
                            #print ""
                            if ss == "H" :
                                print "  - H->_ - at res %d " % rI + ", pos: %d " % res.id.position
                                mv = self.GetMapFromMolRes ( mol, chainId, rStart, res.id.position-1 )
                                chain_maps.append ( [mv, self.MapIndexesInMap ( dmap, mv )] )
                                mv.chain_id = basename + "_" + chainId + "_%d" % rStart
                                rStart = res.id.position
                                ss = ""

                    rI += 1
                    if rI >= len(residues) :
                        print "  - done chain " + chainId + " - at res %d " % rI + ", pos: %d " % res.id.position,

                        if res.isHelix :
                            print " - H "
                            mv = self.GetMapFromMolRes ( mol, chainId, rStart, res.id.position )
                            chain_maps.append ( [mv, self.MapIndexesInMap ( dmap, mv )] )
                            mv.chain_id = basename + "_" + chainId + "_" + ss + "%d" % rStart
                        else :
                            print ""
                            if len(oRanges) > 0 : oRanges = oRanges + ","
                            oRanges = oRanges + "%d-%d.%s" % (rStart, res.id.position, chainId)


                        mv = self.GetMapFromMolRanges ( mol, chainId, oRanges )
                        chain_maps.append ( [mv, self.MapIndexesInMap ( dmap, mv )] )
                        mv.chain_id = basename + "_" + chainId

                        break;


            #break

        # print chain_maps

        rgroups = {}

        print " - %d regions" % len(smod.regions)

        for ri, reg in enumerate ( smod.regions ) :

            if ri % 100 == 0 :
              print " %d/%d " % (ri+1, len(smod.regions) )

            max_ov = 0.0
            max_ov_chm = None
            for chmImap in chain_maps :
                chm, imap = chmImap
                ipoints = reg.points()
                noverlap = 0
                for i,j,k in ipoints :
                    if (i,j,k) in imap:
                        noverlap += 1

                #print " - ", chm.name, noverlap

                ov = float(noverlap) / reg.point_count()
                if ov > max_ov :
                    max_ov = ov
                    max_ov_chm = chm

            if max_ov_chm :
                try : rgroups[max_ov_chm.chain_id]
                except : rgroups[max_ov_chm.chain_id] = []
                rgroups[max_ov_chm.chain_id].append ( reg )


        import regions

        for chid, regs in rgroups.iteritems () :
            print "Chain %s - %d regions" % (chid, len(regs))

            jregs = regions.TopParentRegions(regs)
            jreg = smod.join_regions ( jregs )
            jreg.make_surface(None, None, smod.regions_scale)


        for chmImap in chain_maps :
          chimera.openModels.close ( chmImap )



    def GroupRegionsByChains ( self ) :

        dmap = segmentation_map()
        if dmap == None :
            umsg ( "Please choose map in Segment Dialog" )
            return

        #fmap = self.MoleculeMap()
        #if fmap == None : return

        smod = self.CurrentSegmentation()
        if smod is None :
            umsg ( "Please select a segmentation in Segment Dialog" )
            return

        print "---"

        #mol = fmap.mols[0]
        #print mol.name

        #chain_colors = RandColorChains ( mol )


        umsg ( "Grouping with chains... making chain maps..." )


        chain_maps = []

        res = float ( self.simRes.get() )
        grid = float ( self.simGridSp.get() )

        print "_____________ res %2f _______ grid %.2f _________________________________" % (res, grid)


        mols = []

        for mol in chimera.openModels.list() :

            if type(mol) != chimera.Molecule or mol.display == False :
                continue

            mols.append ( mol )

        nchains = 0
        from random import random as rand

        mol_ch_colors = {}

        for i, mol in enumerate (mols) :

            chain_colors = {} # RandColorChains ( mol )
            for r in mol.residues:
                if hasattr ( r, 'ribbonColor' ) and r.ribbonColor != None :
                    chain_colors[r.id.chainId] = r.ribbonColor.rgba()
                    mol_ch_colors[mol.name + "_" + r.id.chainId] = r.ribbonColor.rgba()
                else :
                    if not r.id.chainId in chain_colors :
                        clr = ( rand()*.7, rand()*.7, rand()*.7, 1.0 )
                        chain_colors[r.id.chainId] = clr
                        mol_ch_colors[mol.name + "_" + r.id.chainId] = clr


            ci = 1

            for cid, clr in chain_colors.iteritems() :

                umsg ( "Grouping with chains... making map for chain %d/%d of mol %d/%d" % (ci,len(chain_colors),i+1,len(mols)) )
                ci += 1
                nchains += 1

                basename = os.path.splitext ( mol.name )[0]
                #cname = basename + "_" + cid
                cname = basename + "_" + cid
                #cname = basename.split ("__")[-1]

                sel_str = "#%d:.%s" % (mol.id, cid)
                print "%s [%s]" % (cname, sel_str),

                cmd = "molmap %s %f sigmaFactor 0.187 gridSpacing %f replace false" % ( sel_str, res, grid )
                chimera.runCommand ( cmd )

                if cid.lower() == cid :
                    cid = "_" + cid

                cname = basename + "_" + cid

                mv = None
                for mod in chimera.openModels.list() :
                  ts = mod.name.split()
                  if len(ts) > 1 and mod.name.find("map") >=0 and mod.name.find("res") >=0 :
                      #print " - found", mod.name
                      mv = mod
                      mv.name = cname
                      break

                if mv == None :
                  umsg (" - error - could not find chain map")
                  return

                imap = self.MapIndexesInMap ( dmap, mv )

                chain_maps.append ( [mv, imap] )
                mv.mol_name = mol.name
                mv.chain_id = cid # cname

            #break

        # print chain_maps

        rgroups = {}

        print " - %d regions" % len(smod.regions)

        for ri, reg in enumerate ( smod.regions ) :

            if ri % 100 == 0 :
              umsg ( "Grouping regions... %d/%d " % (ri+1, len(smod.regions) ) )

            max_ov = 0.0
            max_ov_chm = None
            for chmImap in chain_maps :
                chm, imap = chmImap
                ipoints = reg.points()
                noverlap = 0
                for i,j,k in ipoints :
                    if (i,j,k) in imap:
                        noverlap += 1

                #print " - ", chm.name, noverlap

                ov = float(noverlap) / reg.point_count()
                if ov > max_ov :
                    max_ov = ov
                    max_ov_chm = chm

            if max_ov_chm :
                if not max_ov_chm in rgroups : # max_ov_chm.chain_id?
                    rgroups[max_ov_chm] = []
                rgroups[max_ov_chm].append ( reg )


        import regions
        from Segger.extract_region_dialog import dialog as exdialog

        base = ""
        if exdialog() != None :
            base = exdialog().saveMapsBaseName.get()

        for mv, regs in rgroups.iteritems () :

            #cid = chid.split("_")[-1]

            print "%s.%s - %d regions" % (mv.mol_name, mv.chain_id, len(regs))

            jregs = regions.TopParentRegions(regs)
            jreg = smod.join_regions ( jregs )
            jreg.color = (.7,.7,.7,1)
            mvId = mv.mol_name + "_" + mv.chain_id
            if mvId in mol_ch_colors :
                jreg.color = mol_ch_colors[mvId]
            jreg.make_surface(None, None, smod.regions_scale)
            jreg.chain_id = mv.chain_id

            #jreg.chain_id = chid

            if 0 and exdialog() != None :
                exdialog().saveMapsBaseName.set( base % cid )
                exdialog().Extract2 ( dmap, dmap, smod, [jreg] )


        for chmImap in chain_maps :
          chimera.openModels.close ( chmImap )


        umsg ( "Done - total %d chains in %d visible Molecules" % (nchains,len(mols)) )


    def MaskWithSel ( self ) :

        selats = chimera.selection.currentAtoms()
        print "%d selected atoms" % len(selats)

        dmap = None
        for m in chimera.openModels.list() :
            if m.display == True and type(m) == VolumeViewer.volume.Volume :
                dmap = m
                break

        if dmap == None :
            return

        print "map: %s" % dmap.name


        import _multiscale
        points = _multiscale.get_atom_coordinates ( selats, transformed = True )

        import _contour
        _contour.affine_transform_vertices ( points, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )

        s = dmap.data.step[0]
        s2 = numpy.sqrt ( s*s + s*s + s*s )
        mdata = VolumeData.zone_masked_grid_data ( dmap.data, points, numpy.sqrt(s2) )

        from VolumeFilter import gaussian
        gvm = gaussian.gaussian_convolution ( mdata.full_matrix(), (.1,.1,.1) )
        #gvm = gvol.full_matrix()

        gdata = VolumeData.Array_Grid_Data ( gvm, dmap.data.origin, dmap.data.step, dmap.data.cell_angles, name = dmap.name + "_m" )
        nvg = VolumeViewer.volume.volume_from_grid_data ( gdata )
        nvg.name = dmap.name + "___"





    def GroupRegionsByMols ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        #fmap = self.MoleculeMap()
        #if fmap == None : return

        smod = self.CurrentSegmentation()
        if smod is None : return

        print "---"

        #mol = fmap.mols[0]
        #print mol.name

        #chain_colors = RandColorChains ( mol )

        chain_maps = []

        res = float ( self.simRes.get() )
        grid = float ( self.simGridSp.get() )

        print "_____________ res %2f _______ grid %.2f ______________" % (res, grid)


        for mol in chimera.openModels.list() :

            if type(mol) != chimera.Molecule or mol.display == False : continue

            chain_colors = RandColorChains ( mol )

            basename = os.path.splitext ( mol.name )[0]
            cname = basename
            sel_str = "#%d" % (mol.id)
            print "%s [%s]" % (mol.name, sel_str),

            cmd = "molmap %s %f sigmaFactor 0.187 gridSpacing %f replace false" % ( sel_str, res, grid )
            chimera.runCommand ( cmd )

            mv = None
            for mod in chimera.openModels.list() :
              ts = mod.name.split()
              if len(ts) > 1 and mod.name.find("map") >=0 and mod.name.find("res") >=0 :
                  #print " - found", mod.name
                  mv = mod
                  mv.name = cname
                  break

            if mv == None :
              umsg (" - error - could not find chain map")
              return

            imap = self.MapIndexesInMap ( dmap, mv )

            chain_maps.append ( [mv, imap] )
            mv.chain_id = cname

            #break

        # print chain_maps

        rgroups = {}

        print " - %d regions" % len(smod.regions)

        for ri, reg in enumerate ( smod.regions ) :

            if ri % 100 == 0 :
              status ( " %d/%d " % (ri+1, len(smod.regions) ) )

            max_ov = 0.0
            max_ov_chm = None
            for chmImap in chain_maps :
                chm, imap = chmImap
                ipoints = reg.points()
                noverlap = 0
                for i,j,k in ipoints :
                    if (i,j,k) in imap:
                        noverlap += 1

                #print " - ", chm.name, noverlap

                ov = float(noverlap) / reg.point_count()
                if ov > 0.8 and ov > max_ov :
                    max_ov = ov
                    max_ov_chm = chm

            if max_ov_chm :
                try : rgroups[max_ov_chm.chain_id]
                except : rgroups[max_ov_chm.chain_id] = []
                rgroups[max_ov_chm.chain_id].append ( reg )


        import regions

        for chid, regs in rgroups.iteritems () :
            print "Chain %s - %d regions" % (chid, len(regs))

            jregs = regions.TopParentRegions(regs)
            jreg = smod.join_regions ( jregs )
            jreg.make_surface(None, None, smod.regions_scale)


        for chmImap in chain_maps :
          chimera.openModels.close ( chmImap )



    def GroupRegionsByFittedMols ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        #fmap = self.MoleculeMap()
        #if fmap == None : return

        smod = self.CurrentSegmentation()
        if smod is None : return

        print "---"

        #mol = fmap.mols[0]
        #print mol.name

        #chain_colors = RandColorChains ( mol )

        chain_maps = []

        res = float ( self.simRes.get() )
        grid = float ( self.simGridSp.get() )

        print "_____________ res %2f _______ grid %.2f _________________________________" % (res, grid)


        #lfits = self.selected_listbox_fits()

        lfits = self.list_fits

        if len(lfits) == 0:
            umsg ( "No selected fitted molecules" )
            return

        umsg('Looking at %d fitted molecules' % len(lfits))

        fit_i = 1
        for fmap, dmap, mat, corr, aI, bI, bC, bO, regions in lfits:

            for mol in fmap.mols :
                if mol.__destroyed__:
                    umsg('Fit molecule was closed - ')
                    return

            self.place_molecule(fmap, mat, dmap)

            mol = fmap.mols[0]

            basename = os.path.splitext ( mol.name )[0]
            cname = basename
            sel_str = "#%d" % (mol.id)
            print "%s [%s]" % (mol.name, sel_str),

            cmd = "molmap %s %f sigmaFactor 0.187 gridSpacing %f replace false" % ( sel_str, res, grid )
            chimera.runCommand ( cmd )

            mv = None
            for mod in chimera.openModels.list() :
              ts = mod.name.split()
              if len(ts) > 1 and mod.name.find("map") >=0 and mod.name.find("res") >=0 :
                  #print " - found", mod.name
                  mv = mod
                  mv.name = cname
                  break

            if mv == None :
              umsg (" - error - could not find chain map")
              return

            imap = self.MapIndexesInMap ( dmap, mv )

            chain_maps.append ( [fit_i, imap] )
            #mv.chain_id = cname
            fit_i += 1
            chimera.openModels.close ( mv )

            #break

        # print chain_maps

        rgroups = {}

        print " - %d regions" % len(smod.regions)

        for ri, reg in enumerate ( smod.regions ) :

            if ri % 1000 == 0 :
              #print " %d/%d " % (ri, len(smod.regions) )
              status ( " %d/%d " % (ri, len(smod.regions) ) )

            max_ov = 0.1
            max_ov_chm = 0
            for chmImap in chain_maps :
                fit_i, imap = chmImap
                ipoints = reg.points()
                noverlap = 0
                for i,j,k in ipoints :
                    if (i,j,k) in imap:
                        noverlap += 1


                ov = float(noverlap) / reg.point_count()

                #print " - fit %d to reg %d, num ov %d, ov %.2f " % (fit_i, ri, noverlap, ov)

                if ov > max_ov :
                    max_ov = ov
                    max_ov_chm = fit_i

            if max_ov_chm > 0 :
                try : rgroups[max_ov_chm]
                except : rgroups[max_ov_chm] = []
                rgroups[max_ov_chm].append ( reg )


        import regions

        print  len( rgroups.keys() ), "groups"

        sregs = []

        for chid, regs in rgroups.iteritems () :
            print "Fit %d - %d regions" % (chid, len(regs))

            jregs = regions.TopParentRegions(regs)
            jreg = smod.join_regions ( jregs )
            sregs.append ( jreg )
            jreg.make_surface(None, None, smod.regions_scale)


        return

        #sel_regs = set ( smod.selected_regions() )
        surfs = [r.surface_piece for r in sregs
                 if not r in sel_regs and r.surface_piece]

        chimera.selection.clearCurrent ()
        chimera.selection.addCurrent ( surfs )

        smod.remove_regions ( regs, remove_children = True )




    def GroupRegionsByVisiMaps ( self ) :

        umsg ( "Grouping by visible maps..." )

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        #fmap = self.MoleculeMap()
        #if fmap == None : return

        smod = self.CurrentSegmentation()
        if smod is None : return

        print "---"

        #mol = fmap.mols[0]
        #print mol.name

        #chain_colors = RandColorChains ( mol )

        chain_maps = []

        res = float ( self.simRes.get() )
        grid = float ( self.simGridSp.get() )

        print "_____________ res %2f _______ grid %.2f _________________________________" % (res, grid)


        for mmap in chimera.openModels.list() :

            if type(mmap) != VolumeViewer.volume.Volume or mmap.display == False : continue

            print " -- map: ", mmap.name

            imap = self.MapIndexesInMap ( dmap, mmap )
            chain_maps.append ( [mmap, imap] )
            mmap.chain_id = mmap.name

            #break

        # print chain_maps

        rgroups = {}

        print " - %d regions" % len(smod.regions)

        for ri, reg in enumerate ( smod.regions ) :

            if ri % 100 == 0 :
              status ( " %d/%d " % (ri+1, len(smod.regions) ) )
              print ",",

            max_ov = 0.0
            max_ov_chm = None
            for chmImap in chain_maps :
                chm, imap = chmImap
                ipoints = reg.points()
                noverlap = 0
                for i,j,k in ipoints :
                    if (i,j,k) in imap:
                        noverlap += 1

                #print " - ", chm.name, noverlap

                ov = float(noverlap) / reg.point_count()
                if ov > max_ov :
                    max_ov = ov
                    max_ov_chm = chm

            if max_ov_chm :
                try : rgroups[max_ov_chm.chain_id]
                except : rgroups[max_ov_chm.chain_id] = []
                rgroups[max_ov_chm.chain_id].append ( reg )


        import regions
        print "."

        for chid, regs in rgroups.iteritems () :
            print "Chain %s - %d regions" % (chid, len(regs))

            jregs = regions.TopParentRegions(regs)
            jreg = smod.join_regions ( jregs )
            jreg.make_surface(None, None, smod.regions_scale)

        umsg ( "Done grouping by visible maps" )


        #for chmImap in chain_maps :
        #  chimera.openModels.close ( chmImap )



	# -----------------------------------------------------------------------------------------------------


    def ZeroMapBySel ( self ) :

        print "0"


    def ZeroMapByMols ( self ) :

        print "0"



    def ZeroMapFittedMols ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        vmat = dmap.full_matrix().copy()

        print "---"

        res = float ( self.simRes.get() )
        grid = float ( self.simGridSp.get() )

        print "_____________ res %2f _______ grid %.2f _________________________________" % (res, grid)


        #lfits = self.selected_listbox_fits()

        lfits = self.list_fits

        if len(lfits) == 0:
            umsg ( "No selected fitted molecules" )
            return

        umsg('Looking at %d fitted molecules' % len(lfits))


        for fmap, dmap, mat, corr, aI, bI, bC, bO, regions in lfits:

            for mol in fmap.mols :
                if mol.__destroyed__:
                    umsg('Fit molecule was closed - ')
                    return

            self.place_molecule(fmap, mat, dmap)

            mol = fmap.mols[0]

            basename = os.path.splitext ( mol.name )[0]
            cname = basename
            sel_str = "#%d" % (mol.id)
            print "%s [%s]" % (mol.name, sel_str),

            cmd = "molmap %s %f sigmaFactor 0.187 gridSpacing %f replace false" % ( sel_str, res, grid )
            chimera.runCommand ( cmd )

            mv = None
            for mod in chimera.openModels.list() :
              ts = mod.name.split()
              if len(ts) > 1 and mod.name.find("map") >=0 and mod.name.find("res") >=0 :
                  #print " - found", mod.name
                  mv = mod
                  mv.name = cname
                  break

            if mv == None :
				umsg (" - error - could not find chain map")
				return

            self.ZeroMatWitMap ( vmat, dmap, mv )
            chimera.openModels.close ( mv )

            #break

        # print chain_maps

        nname = os.path.splitext(dmap.name)[0] + "_zeroed"

        from VolumeData import Array_Grid_Data
        mgrid = Array_Grid_Data ( vmat, dmap.data.origin, dmap.data.step, dmap.data.cell_angles, name=nname)
        import VolumeViewer
        nv = VolumeViewer.volume_from_grid_data ( mgrid, show_data = False, show_dialog = False )
        nv.name = nname
        #nv.copy_settings_from(volume)
        nv.show()



    def ZeroMapVisMols ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        vmat = dmap.full_matrix().copy()

        for m in chimera.openModels.list() :
            if type(m) == chimera.Molecule and m.display == True :
                print m.name
                vmat = self.ZeroMatWitMol ( vmat, dmap, m  )


        nname = os.path.splitext(dmap.name)[0] + "_zeroed"

        from VolumeData import Array_Grid_Data
        mgrid = Array_Grid_Data ( vmat, dmap.data.origin, dmap.data.step, dmap.data.cell_angles, name=nname)

        import VolumeViewer
        #nv = VolumeViewer.volume_from_grid_data ( mgrid, show_data = False, show_dialog = False )

        try : df_v = VolumeViewer.volume.add_data_set ( mgrid, None )
        except : df_v = VolumeViewer.volume.volume_from_grid_data ( mgrid )

        df_v.name = nname
        #nv.copy_settings_from(volume)
        df_v.show()




    def ValuesInMap ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        dmat = dmap.full_matrix().copy()

        chain_maps = []

        res = float ( self.simRes.get() )
        grid = float ( self.simGridSp.get() )

        print "_____________ res %2f _______ grid %.2f _________________________________" % (res, grid)


        for mmap in chimera.openModels.list() :

            if type(mmap) != VolumeViewer.volume.Volume or mmap.display == False : continue

            print " -- map: ", mmap.name

            imap = self.MapIndexesInMap ( dmap, mmap )
            chain_maps.append ( [mmap, imap] )
            mmap.chain_id = mmap.name

            #break

        sumd = 0
        n = 0
        imap = set()
        values = []

        for cm, points in chain_maps :
            print " --- at map --- : " + cm.name
            #n += len(points)
            for fi, fj, fk in points :
                val = dmat[fk,fj,fi]
                if val < 17 :
                    sumd += dmat[fk,fj,fi]
                    values.append ( dmat[fk,fj,fi] )
                    #print dmat[fk,fj,fi]
                    n += 1.0

        avg = sumd / len(points)

        print " Average value: %.3f"%avg + " at %d"%len(points) + " points"

        return;

        import numpy
        min = numpy.min(values)
        amax = numpy.max(values)
        max = numpy.min ( [amax, 16] )
        d = 1
        print "Min: %.2f, max: %.2f (%.2f), step %.2f, avg %.2f" % (min, amax, max, d, numpy.average(values))
        buckets = numpy.zeros ( (max - min) / d )
        for v in values :
            if v > max :
                continue
            r = (v - min) / (max - min)
            bi = int ( numpy.floor ( r * (len(buckets)-1) ) )
            buckets[bi] += 1

        for bi, bnum in enumerate ( buckets ) :
            bv = float (bi) / float(len(buckets)) * (max - min) + min
            print "%.1f,%d" % (bv,bnum)





    def MaskMapWithSel ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        fmap = None
        [fmol, fmap2] = self.MolOrMapSelected ();
        if fmol != None :
            print " - got molecule map"
            fmap = self.MoleculeMap()
        elif fmap2 != None :
            print " - got map map"
            fmap = fmap2
        else :
            umsg ('Please select an open molecule or map in the field above')
            return

        df_mat = self.Map2Map ( fmap, dmap )

        print " - using surf level %.5f for mask" % fmap.surface_levels[0]

        if 1 :

            try :
                res = float ( self.simRes.get() )
            except :
                umsg ( "Invalid resolution entered, please enter a number" )
                return


            s = dmap.data.step # A/pixel
            diag_l = numpy.sqrt ( s[0]*s[0] + s[1]*s[1] + s[2]*s[2] ) # A/pixel
            num_it = res / diag_l # how many iterations will reach the desired width
            numit = int ( numpy.ceil ( num_it ) )

            print " - using res %.4f for dropoff, diag is %.3f, #it: %d" % (res, diag_l, numit)

            in_mask = numpy.where ( df_mat > fmap.surface_levels[0], numpy.ones_like(df_mat), numpy.zeros_like(df_mat) )
            out_mask = numpy.where ( in_mask > 0, numpy.zeros_like(in_mask), numpy.ones_like(in_mask) )
            gvm = in_mask.copy();
            for i in range (numit) :
                nv_1 = numpy.roll(gvm, 1, axis=0)
                nv_2 = numpy.roll(gvm, -1, axis=0)
                nv_3 = numpy.roll(gvm, 1, axis=1)
                nv_4 = numpy.roll(gvm, -1, axis=1)
                nv_5 = numpy.roll(gvm, 1, axis=2)
                nv_6 = numpy.roll(gvm, -1, axis=2)
                gvm = 1.0/6.0 * ( nv_1 + nv_2 + nv_3 + nv_4 + nv_5 + nv_6 )
                gvm = out_mask * gvm + in_mask
                umsg ("Adding drop-off - iteration %d" % i)

            df_mat = gvm


        mmat = dmap.data.full_matrix() * df_mat

        df_data = VolumeData.Array_Grid_Data ( df_mat, dmap.data.origin, dmap.data.step, dmap.data.cell_angles, name=(dmap.name + "__MaskedWith__" + fmap.name) )

        try : df_v = VolumeViewer.volume.add_data_set ( df_data, None )
        except : df_v = VolumeViewer.volume.volume_from_grid_data ( df_data )


        if 0 :
            mapMean, mapStDev = MapStats ( df_v )
            df_v = AddNoiseToMap ( df_v, mapMean, mapStDev / 3.0 )


        df_v.name = dmap.name + "__MaskedWith__" + fmap.name
        df_v.openState.xform = dmap.openState.xform





	# -----------------------------------------------------------------------------------------------------




    def MolOrMapSelected ( self ) :

        label = self.struc.get()

        if len(label) == 0 :
            umsg ( "No structure selected" )
            return [None, None]

        mod_num = label [ label.rfind("(")+1 : label.rfind(")") ]

        if len(mod_num) == 0 :
            # this isn't possible given a (#) is added to each name...
            umsg ( "An internal error that shouldn't happen did." )
            return [None, None]

        sel_str = "#" + mod_num

        fmol = None
        try :
            fmol = chimera.selection.OSLSelection ( sel_str ).molecules()[0]
        except :
            print ( "Selected model is not a molecule..." )

        fmap = None
        try :
            fmap = chimera.selection.OSLSelection ( sel_str ).models()[0]
        except :
            print ( "Selected model is not a map..." )

        return [fmol, fmap]



    def Map2Map ( self, densitiesFromMap, toGridOfMap, mask = False ) :

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

        return df_mat




    def Map2MapResize (self, fmap, dmap) :

        import axes
        fpoints, weights = axes.map_points ( fmap )
        print "Fit map - got %d points in contour" % len (fpoints)

        from _contour import affine_transform_vertices as transform_vertices
        #print "Fit map - xf: ", fmap.openState.xform
        transform_vertices( fpoints,  Matrix.xform_matrix( fmap.openState.xform ) )
        #print "Seg map - xf: ", dmap.openState.xform
        transform_vertices( fpoints,  Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
        transform_vertices ( fpoints, dmap.data.xyz_to_ijk_transform )
        #print "points in %s ref:" % dmap.name, fpoints

        bound = 5
        li,lj,lk = numpy.min ( fpoints, axis=0 ) - (bound, bound, bound)
        hi,hj,hk = numpy.max ( fpoints, axis=0 ) + (bound, bound, bound)

        n1 = hi - li + 1
        n2 = hj - lj + 1
        n3 = hk - lk + 1

        #print " - bounds - %d %d %d --> %d %d %d --> %d %d %d" % ( li, lj, lk, hi, hj, hk, n1,n2,n3 )

        #nmat = numpy.zeros ( (n1,n2,n3), numpy.float32 )
        #dmat = dmap.full_matrix()

        nstep = (fmap.data.step[0], fmap.data.step[1], fmap.data.step[2] )
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

        ox = round ( nO[0]/dmap.data.step[0] ) * dmap.data.step[0]
        oy = round ( nO[1]/dmap.data.step[1] ) * dmap.data.step[1]
        oz = round ( nO[2]/dmap.data.step[2] ) * dmap.data.step[2]

        nO = ( ox, oy, oz )

        print " - new map origin:", nO


        nmat = numpy.zeros ( (nn1,nn2,nn3), numpy.float32 )
        ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles )

        #print " - fmap grid dim: ", numpy.shape ( fmap.full_matrix() )
        #print " - new map grid dim: ", numpy.shape ( nmat )

        npoints = grid_indices ( (nn1, nn2, nn3), numpy.single)  # i,j,k indices
        transform_vertices ( npoints, ndata.ijk_to_xyz_transform )

        dvals = fmap.interpolated_values ( npoints, dmap.openState.xform )
        #dvals = numpy.where ( dvals > threshold, dvals, numpy.zeros_like(dvals) )
        #nze = numpy.nonzero ( dvals )

        nmat = dvals.reshape( (nn3,nn2,nn1) )
        #f_mat = fmap.data.full_matrix()
        #f_mask = numpy.where ( f_mat > fmap.surface_levels[0], numpy.ones_like(f_mat), numpy.zeros_like(f_mat) )
        #df_mat = df_mat * f_mask

        ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles )
        try : nv = VolumeViewer.volume.add_data_set ( ndata, None )
        except : nv = VolumeViewer.volume.volume_from_grid_data ( ndata )


        fmap_base = os.path.splitext(fmap.name)[0]
        dmap_base = os.path.splitext(dmap.name)[0]
        fmap_path = os.path.splitext (fmap.data.path)[0]
        dmap_path = os.path.splitext (dmap.data.path)[0]

        nv.name = fmap_base + "__in__" + dmap_base
        nv.openState.xform = dmap.openState.xform

        #npath = dmap_path + fnamesuf
        #nv.write_file ( npath, "mrc" )
        #print "Wrote ", npath

        return nv



    def TakeDMap_with_FMap ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        fmap = None
        [fmol, fmap2] = self.MolOrMapSelected ();
        if fmol != None :
            print " - got molecule map"
            fmap = self.MoleculeMap()
        elif fmap2 != None :
            print " - got map map"
            fmap = fmap2
        else :
            umsg ('Please select an open molecule or map in the field above')
            return


        df_mat = self.Map2Map ( dmap, fmap )
        df_data = VolumeData.Array_Grid_Data ( df_mat, fmap.data.origin, fmap.data.step, fmap.data.cell_angles )

        try : df_v = VolumeViewer.volume.add_data_set ( df_data, None )
        except : df_v = VolumeViewer.volume.volume_from_grid_data ( df_data )
        df_v.name = dmap.name + "_in_" + fmap.name
        df_v.openState.xform = fmap.openState.xform



    def TakeFMap_with_DMap0 ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        fmap = None
        [fmol, fmap2] = self.MolOrMapSelected ();
        if fmol != None :
            print " - got molecule map"
            fmap = self.MoleculeMap()
        elif fmap2 != None :
            print " - got map map"
            fmap = fmap2
        else :
            umsg ('Please select an open molecule or map in the field above')
            return

        df_mat = self.Map2Map ( fmap, dmap )
        df_data = VolumeData.Array_Grid_Data ( df_mat, dmap.data.origin, dmap.data.step, dmap.data.cell_angles )

        try : df_v = VolumeViewer.volume.add_data_set ( df_data, None )
        except : df_v = VolumeViewer.volume.volume_from_grid_data ( df_data )


        if 0 :
            mapMean, mapStDev = MapStats ( df_v )
            df_v = AddNoiseToMap ( df_v, mapMean, mapStDev / 3.0 )


        df_v.name = fmap.name + "_in_" + dmap.name
        df_v.openState.xform = dmap.openState.xform



    def TakeFMap_with_DMapN ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        fmap = None
        [fmol, fmap2] = self.MolOrMapSelected ();
        if fmol != None :
            print " - got molecule map"
            fmap = self.MoleculeMap()
        elif fmap2 != None :
            print " - got map map"
            fmap = fmap2
        else :
            umsg ('Please select an open molecule or map in the field above')
            return

        df_mat = self.Map2Map ( fmap, dmap )
        df_data = VolumeData.Array_Grid_Data ( df_mat, dmap.data.origin, dmap.data.step, dmap.data.cell_angles )

        try : df_v = VolumeViewer.volume.add_data_set ( df_data, None )
        except : df_v = VolumeViewer.volume.volume_from_grid_data ( df_data )


        if 1 :
            mapMean, mapStDev = MapStats ( df_v )
            df_v = AddNoiseToMap ( df_v, mapMean, mapStDev / 3.0 )


        df_v.name = fmap.name + "_in_" + dmap.name
        df_v.openState.xform = dmap.openState.xform




    def TakeFMap_with_DMap ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        fmap = None
        [fmol, fmap2] = self.MolOrMapSelected ();
        if fmol != None :
            print " - got molecule map"
            fmap = self.MoleculeMap()
        elif fmap2 != None :
            print " - got map map"
            fmap = fmap2
        else :
            umsg ('Please select an open molecule or map in the field above')
            return

        nv = self.Map2MapResize ( fmap, dmap )




    def FitAllVisMaps ( self ) :

        from VolumeViewer import Volume
        from chimera import Molecule
        mlist = OML(modelTypes = [Volume,Molecule])
        for m in mlist :
            label = m.name + " (%d)" % m.id
            print "---------------------", label, "---------------------"

            if m.display == False :
                continue

            self.struc.set(label)
            self.cur_mol = m
            self.Fit()



    def AvgFMaps ( self ) :

        from VolumeViewer import Volume
        mlist = OML(modelTypes = [Volume])

        fmap = None
        avgMat = None
        N = 0.0

        for m in mlist :
            if m.display == True :
                print m.name

                if avgMat == None :
                    avgMat = m.data.full_matrix()
                    fmap = m
                    N = 1.0
                else :
                    avgMat = avgMat + m.data.full_matrix()
                    N = N + 1.0

        avgMat = avgMat / N

        df_data = VolumeData.Array_Grid_Data ( avgMat, fmap.data.origin, fmap.data.step, fmap.data.cell_angles )

        try : df_v = VolumeViewer.volume.add_data_set ( df_data, None )
        except : df_v = VolumeViewer.volume.volume_from_grid_data ( df_data )
        df_v.name = "Avg"
        df_v.openState.xform = fmap.openState.xform

        return


    def DifFMaps2 ( self ) :

        from VolumeViewer import Volume
        mlist = OML(modelTypes = [Volume])

        smaps = []
        for m in mlist :
            if m.display == True :
                print " - ", m.name
                smaps.append ( m )


        if len(smaps) != 2 :
            umsg ( "Need only 2 maps visible" )
            return


        m1 = smaps[0]
        m2 = smaps[1]


        m2_mat = self.Map2Map ( m2, m1 )

        difMat = NormalizeMat ( m1.data.full_matrix() ) - NormalizeMat ( m2_mat )
        df_data = VolumeData.Array_Grid_Data ( difMat, m1.data.origin, m1.data.step, m1.data.cell_angles )
        try : df_v = VolumeViewer.volume.add_data_set ( df_data, None )
        except : df_v = VolumeViewer.volume.volume_from_grid_data ( df_data )
        df_v.name = m1.name + "__-__" + m2.name
        df_v.openState.xform = m1.openState.xform

        difMat = NormalizeMat ( m2_mat ) - NormalizeMat ( m1.data.full_matrix() )
        df_data = VolumeData.Array_Grid_Data ( difMat, m1.data.origin, m1.data.step, m1.data.cell_angles )
        try : df_v = VolumeViewer.volume.add_data_set ( df_data, None )
        except : df_v = VolumeViewer.volume.volume_from_grid_data ( df_data )
        df_v.name = m2.name + "__-__" + m1.name
        df_v.openState.xform = m1.openState.xform


        return



    def AvgFMaps2 (self) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        import extract_region_dialog
        reload ( extract_region_dialog )

        superSampleBy = 1

        if superSampleBy > 1 :
            ndata = extract_region_dialog.MapSS ( dmap, superSampleBy )
            try : nmap = VolumeViewer.volume.add_data_set ( ndata, None )
            except : nmap = VolumeViewer.volume.volume_from_grid_data ( ndata )
            nmap.name = dmap.name + "_M%d" % n
            nmap.openState.xform = dmap.openState.xform
            dmap = nmap


        from VolumeViewer import Volume
        mlist = OML(modelTypes = [Volume])

        if 0 :
            print " -- making base map -- "

            bmap = None
            m0 = None
            for m in mlist :
                if m.display == True :
                    print m.name
                    if bmap == None :
                        bmap = m
                        m0 = m
                    else :
                        bmap0 = bmap
                        bmap = self.Map2MapResize (m, bmap)
                        if bmap0 != m0 :
                            chimera.openModels.close ( [bmap0] )


            bmap.name = "base"
            dmap = bmap

        if 0 :
            print " -- finding base map --- "
            largestMap = None
            maxD = 0
            for m in mlist :
                if m.display == True :
                    d = numpy.sum ( m.data.size )
                    if d > maxD :
                        maxD = d
                        largestMap = m

            print " - largest map: ", largestMap.name



        fmap = None
        [fmol, fmap2] = self.MolOrMapSelected ();
        if fmol != None :
            umsg ('Please select the base map in the field at the top - molecule found')
            return
        elif fmap2 != None :
            print " - got map map"
            fmap = fmap2
        else :
            umsg ('Please select the base map in the field at the top')
            return

        dmap = fmap2
        #dmap.display = False
        umsg ( "Using as base map: %s" % dmap.name )

        avgMat = dmap.full_matrix()
        fmap = dmap

        weights = avgMat.ravel()
        smin = numpy.min (weights)
        sdev = numpy.std (weights)
        savg = numpy.average(weights)
        smax = numpy.max (weights)

        print " - (%.4f,%.4f) |%.4f| +/- %.4f" % (smin, smax, savg, sdev)
        avgMat = avgMat * (1.0 / smax)

        N = 1.0


        #fmap = None
        #avgMat = None
        #N = 0.0

        print " ----------- Averaging... ---------------------"

        for m in mlist :
            if m.display == True and m != dmap :
            #if m.display == True :
                print m.name

                df_mat = self.Map2Map ( m, dmap )
                m.display = False

                weights = df_mat.ravel()
                smin = numpy.min (weights)
                sdev = numpy.std (weights)
                savg = numpy.average(weights)
                smax = numpy.max (weights)

                thr = 0 # m.surface_levels[0]

                print "%s - (%.4f,%.4f) |%.4f| +/- %.4f -- %.4f" % (m.name, smin, smax, savg, sdev, thr)

                N = N + 1.0
                #df_mat = df_mat - ( numpy.ones_like(df_mat)*thr )

                #df_mat = numpy.where ( df_mat > thr, df_mat, numpy.zeros_like(df_mat) )

                df_mat = df_mat * (1.0 / smax)
                #df_mat = df_mat + ( numpy.ones_like(df_mat) * 10.0 )

                if 0 :
                    imhist,bins = numpy.histogram ( df_mat.flatten(), 20, normed=True )
                    print " ------- Histogram:"
                    print imhist
                    print " ------- Bins:"
                    print bins

                    cdf = imhist.cumsum() #cumulative distribution function
                    cdf = 10.0 * cdf / cdf[-1] #normalize

                    print cdf

                    #use linear interpolation of cdf to find new pixel values
                    #df_mat = numpy.interp ( df_mat.flatten(), bins[:-1], cdf )
                    #df_mat = df_mat.reshape(dmap.data.full_matrix().shape)




                if avgMat == None :
                    avgMat = df_mat
                    fmap = m
                else :
                    avgMat = avgMat + df_mat

        print " ----------- n=%f ---------------------" % N

        avgMat = avgMat / N
        df_data = VolumeData.Array_Grid_Data ( avgMat, dmap.data.origin, dmap.data.step, dmap.data.cell_angles, name="avg" )
        try : df_v = VolumeViewer.volume.add_data_set ( df_data, None )
        except : df_v = VolumeViewer.volume.volume_from_grid_data ( df_data )
        df_v.name = "Avg"
        df_v.openState.xform = dmap.openState.xform

        nv = self.ShrinkMap ( df_v, 1e-3 )
        chimera.openModels.close ( [df_v] )

        if 0 :
            stdMat = None
            N = 0.0

            for m in mlist :
                if m.display == True :
                    print m.name

                    if m.name == "Avg" :
                        print "skipping avg vol"
                        continue

                    df_mat = self.Map2Map ( m, dmap )
                    N = N + 1.0

                    print " - sub from avg..."
                    d = numpy.power ( df_mat - avgMat, 2 )
                    if stdMat == None :
                        stdMat = d
                    else :
                        stdMat = stdMat + d

            stdMat = numpy.power ( stdMat / N, 0.5 )
            df_data = VolumeData.Array_Grid_Data ( stdMat, dmap.data.origin, dmap.data.step, dmap.data.cell_angles )
            try : df_v = VolumeViewer.volume.add_data_set ( df_data, None )
            except : df_v = VolumeViewer.volume.volume_from_grid_data ( df_data )
            df_v.name = "Stdev"
            df_v.openState.xform = dmap.openState.xform



        # chimera.openModels.close ( dmap )

    def ShrinkMap ( self, dmap, thr ) :

        import axes
        dmap.surface_levels[0] = thr
        fpoints, weights = axes.map_points ( dmap )
        #print "%s / %f - %d points in contour" % (dmap.name, thr, len (fpoints))

        from _contour import affine_transform_vertices as transform_vertices
        #print "Fit map - xf: ", fmap.openState.xform
        #transform_vertices( fpoints,  Matrix.xform_matrix( fmap.openState.xform ) )

        #print "Seg map - xf: ", dmap.openState.xform
        #transform_vertices( fpoints,  Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
        transform_vertices ( fpoints, dmap.data.xyz_to_ijk_transform )
        #print "points in %s ref:" % dmap.name, fpoints

        bound = 4
        li,lj,lk = numpy.min ( fpoints, axis=0 ) - (bound, bound, bound)
        hi,hj,hk = numpy.max ( fpoints, axis=0 ) + (bound, bound, bound)

        n1 = int(hi - li + 1)
        n2 = int(hj - lj + 1)
        n3 = int(hk - lk + 1)

        #print " - bounds - %d %d %d --> %d %d %d --> %d %d %d" % ( li, lj, lk, hi, hj, hk, n1,n2,n3 )

        #nmat = numpy.zeros ( (n1,n2,n3), numpy.float32 )
        #dmat = dmap.full_matrix()

        O = dmap.data.origin
        #print " - %s origin:" % dmap.name, O
        #print " - %s step:" % dmap.name, dmap.data.step
        nO = ( O[0] + float(li) * dmap.data.step[0],
               O[1] + float(lj) * dmap.data.step[1],
               O[2] + float(lk) * dmap.data.step[2] )

        #print " - new map origin:", nO

        nmat = numpy.zeros ( (n1,n2,n3), numpy.float32 )
        ndata = VolumeData.Array_Grid_Data ( nmat, nO, dmap.data.step, dmap.data.cell_angles )

        #print " - new map grid dim: ", numpy.shape ( nmat )

        npoints = grid_indices ( (n1, n2, n3), numpy.single)  # i,j,k indices
        transform_vertices ( npoints, ndata.ijk_to_xyz_transform )

        dvals = dmap.interpolated_values ( npoints, dmap.openState.xform )
        #dvals = numpy.where ( dvals > threshold, dvals, numpy.zeros_like(dvals) )
        #nze = numpy.nonzero ( dvals )

        nmat = dvals.reshape( (n3,n2,n1) )
        #f_mat = fmap.data.full_matrix()
        #f_mask = numpy.where ( f_mat > fmap.surface_levels[0], numpy.ones_like(f_mat), numpy.zeros_like(f_mat) )
        #df_mat = df_mat * f_mask

        ndata = VolumeData.Array_Grid_Data ( nmat, nO, dmap.data.step, dmap.data.cell_angles, name = dmap.name )
        try : nv = VolumeViewer.volume.add_data_set ( ndata, None )
        except : nv = VolumeViewer.volume.volume_from_grid_data ( ndata )

        return nv




    def TakeFMapsVis ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        from VolumeViewer import Volume
        mlist = OML(modelTypes = [Volume])
        for m in mlist :
            if m.display == True :

                df_mat = self.Map2Map ( m, dmap )
                df_data = VolumeData.Array_Grid_Data ( df_mat, dmap.data.origin, dmap.data.step, dmap.data.cell_angles )
                try : df_v = VolumeViewer.volume.add_data_set ( df_data, None )
                except : df_v = VolumeViewer.volume.volume_from_grid_data ( df_data )
                df_v.openState.xform = dmap.openState.xform

                mdir, mfile = os.path.split(m.data.path)
                df_v.name = "f_" + mfile

                print m.name, "->", df_v.name

                dpath = mdir + "/" + df_v.name
                df_v.write_file ( dpath, "mrc" )



    def DifferenceMap ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        fmap = None
        [fmol, fmap2] = self.MolOrMapSelected ();
        if fmol != None :
            print " - using molecule map"
            fmap = self.MoleculeMap()
        elif fmap2 != None :
            print " - using map map"
            fmap = fmap2
        else :
            umsg ('Please select an open molecule or map in the field above')
            return


        print "\n\nDiff map ", dmap.name, " <=> ", fmap.name

        closeDMap = False

        smod = self.CurrentSegmentation()
        regs = smod.selected_regions()
        if len(regs) > 0 :
            dmap = mask_volume( regs, dmap )
            closeDMap = True



        df_mat = self.Map2Map ( fmap, dmap )
        df_data = VolumeData.Array_Grid_Data ( df_mat, dmap.data.origin, dmap.data.step, dmap.data.cell_angles )

        #MapStats ( dmap )
        #MapDataStats ( df_data )

        print ""
        print "Normalizing", dmap.name
        dmap_data_n = NormalizeData ( dmap.data )
        #dmap_data_n = dmap.data

        if 0 :
            try : nv = VolumeViewer.volume.add_data_set ( dmap_data_n, None )
            except : nv = VolumeViewer.volume.volume_from_grid_data ( dmap_data_n )
            nv.name = os.path.splitext(dmap.name)[0] + "_norm.mrc"
            nv.openState.xform = dmap.openState.xform
        #fmapn = NormalizeMap ( fmap )

        print ""
        print "Normalizing transferred fit map"
        df_data_n = NormalizeData (  df_data )
        #df_data_n = df_data

        if 0 :
            try : nv = VolumeViewer.volume.add_data_set ( df_data_n, None )
            except : nv = VolumeViewer.volume.volume_from_grid_data ( df_data_n )
            nv.name = os.path.splitext(dmap.name)[0] + "_fmap_norm.mrc"
            nv.openState.xform = dmap.openState.xform


        diff_mat = numpy.fabs ( dmap_data_n.full_matrix () - df_data_n.full_matrix () )
        weights = diff_mat.ravel()
        smin = numpy.min (weights)
        sdev = numpy.std (weights)
        savg = numpy.average(weights)
        smax = numpy.max(weights)

        print ""
        print "Difference map:"
        #print " -", len(nz), " nonzero"
        print " - range: %.3f -> %.3f, avg=%.3f, sdev=%.3f" % (smin, smax, savg, sdev)

        diff_data = VolumeData.Array_Grid_Data ( diff_mat, df_data.origin, df_data.step, df_data.cell_angles )

        try : nv = VolumeViewer.volume.add_data_set ( diff_data, None )
        except : nv = VolumeViewer.volume.volume_from_grid_data ( diff_data )

        nv.name = os.path.splitext(dmap.name)[0] + "_--_" + fmap.name
        nv.openState.xform = dmap.openState.xform

        if closeDMap :
            chimera.openModels.close ( [closeDMap] )



    def IntersectionMap ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        fmap = None
        [fmol, fmap2] = self.MolOrMapSelected ();
        if fmol != None :
            print " - using molecule map"
            fmap = self.MoleculeMap()
        elif fmap2 != None :
            print " - using map map"
            fmap = fmap2
        else :
            umsg ('Please select an open molecule or map in the field above')
            return


        print "\n\nDiff map ", dmap.name, " <=> ", fmap.name

        closeDMap = False

        smod = self.CurrentSegmentation()
        regs = smod.selected_regions()
        if len(regs) > 0 :
            dmap = mask_volume( regs, dmap )
            closeDMap = True


        df_mat = self.Map2Map ( fmap, dmap )
        df_data = VolumeData.Array_Grid_Data ( df_mat, dmap.data.origin, dmap.data.step, dmap.data.cell_angles )

        #MapStats ( dmap )
        #MapDataStats ( df_data )

        print ""
        print "Normalizing", dmap.name
        dmap_data_n = NormalizeData ( dmap.data )

        if 0 :
            try : nv = VolumeViewer.volume.add_data_set ( dmap_data_n, None )
            except : nv = VolumeViewer.volume.volume_from_grid_data ( dmap_data_n )
            nv.name = os.path.splitext(dmap.name)[0] + "_norm.mrc"
            nv.openState.xform = dmap.openState.xform
        #fmapn = NormalizeMap ( fmap )

        print ""
        print "Normalizing transferred fit map"
        df_data_n = NormalizeData (  df_data )

        if 0 :
            try : nv = VolumeViewer.volume.add_data_set ( df_data_n, None )
            except : nv = VolumeViewer.volume.volume_from_grid_data ( df_data_n )
            nv.name = os.path.splitext(dmap.name)[0] + "_fmap_norm.mrc"
            nv.openState.xform = dmap.openState.xform


        diff_mat = numpy.fabs ( dmap_data_n.full_matrix () - df_data_n.full_matrix () )
        weights = diff_mat.ravel()
        smin = numpy.min (weights)
        sdev = numpy.std (weights)
        savg = numpy.average(weights)
        smax = numpy.max(weights)

        print ""
        print "Difference map:"
        #print " -", len(nz), " nonzero"
        print " - range: %.3f -> %.3f, avg=%.3f, sdev=%.3f" % (smin, smax, savg, sdev)

        diff_data = VolumeData.Array_Grid_Data ( diff_mat, df_data.origin, df_data.step, df_data.cell_angles )

        try : nv = VolumeViewer.volume.add_data_set ( diff_data, None )
        except : nv = VolumeViewer.volume.volume_from_grid_data ( diff_data )

        nv.name = os.path.splitext(dmap.name)[0] + "_--_" + fmap.name
        nv.openState.xform = dmap.openState.xform

        if closeDMap :
            chimera.openModels.close ( [closeDMap] )


    def ShapeMatch ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        fmap = None
        [fmol, fmap2] = self.MolOrMapSelected ();
        if fmol != None :
            print " - using molecule map"
            fmap = self.MoleculeMap()
        elif fmap2 != None :
            print " - using map map"
            fmap = fmap2
        else :
            umsg ('Please select an open molecule or map in the field above')
            return


        print "\n\nDiff map ", dmap.name, " <=> ", fmap.name

        closeDMap = False
        realDMap = dmap

        smod = self.CurrentSegmentation()
        if smod != None :
            regs = smod.selected_regions()
            if len(regs) > 0 :
                dmap = mask_volume( regs, dmap )
                closeDMap = True


        df_mat = self.Map2Map ( fmap, dmap )
        #df_data = VolumeData.Array_Grid_Data ( df_mat, dmap.data.origin, dmap.data.step, dmap.data.cell_angles )


        thr = dmap.surface_levels[0]

        umsg ("Generating 1/-1 map for " + dmap.name + " thr: %.3f" % thr)

        m0 = dmap.data.full_matrix()
        m1 = numpy.where ( m0 > realDMap.surface_levels[0], numpy.ones_like(m0)*1, numpy.zeros_like(m0) )

        m2 = numpy.where ( df_mat > fmap.surface_levels[0], numpy.ones_like(df_mat)*1, numpy.zeros_like(df_mat) )

        mi = m1 * m2
        mu = m1 + m2

        if 0 :
            mid = VolumeData.Array_Grid_Data ( mi, realDMap.data.origin, realDMap.data.step, realDMap.data.cell_angles, name="inter" )
            mud = VolumeData.Array_Grid_Data ( mu, realDMap.data.origin, realDMap.data.step, realDMap.data.cell_angles, name="union" )

            nv = VolumeViewer.volume_from_grid_data ( mid )
            nv = VolumeViewer.volume_from_grid_data ( mud )


        nz_int =  numpy.shape ( (mi).nonzero () )[1]
        nz_uni =  numpy.shape ( (mu).nonzero () )[1]

        sm_score = float(nz_int) / float (nz_uni)

        print " - intersection %d, union %d - sm: %.3f" % (nz_int, nz_uni, sm_score)


        ndata = VolumeData.Array_Grid_Data ( mi, dmap.data.origin, dmap.data.step, dmap.data.cell_angles )

        try : nv = VolumeViewer.volume.add_data_set ( ndata, None )
        except : nv = VolumeViewer.volume.volume_from_grid_data ( ndata )

        nv.name = os.path.splitext(dmap.name)[0] + "_--_" + fmap.name
        nv.openState.xform = dmap.openState.xform


        if closeDMap :
            chimera.openModels.close ( [dmap] )





def fit_segments_dialog ( create=False ) :

  from chimera import dialogs
  return dialogs.find ( Fit_Segments_Dialog.name, create=create )


def close_fit_segments_dialog ():

    from chimera import dialogs
    d = fit_segments_dialog ()
    if d :
        d.toplevel_widget.update_idletasks ()
        d.Close()
        d.toplevel_widget.update_idletasks ()

def show_fit_segments_dialog ():

    from chimera import dialogs
    d = fit_segments_dialog ( create = True )
    # Avoid transient dialog resizing when created and mapped for first time.
    d.toplevel_widget.update_idletasks ()
    d.enter()
    return d


def new_fit_segments_dialog ( closeExisting = True ):

    if closeExisting : close_fit_segments_dialog ()
    show_fit_segments_dialog ()


# -----------------------------------------------------------------------------
#
from chimera import dialogs
dialogs.register (Fit_Segments_Dialog.name, Fit_Segments_Dialog,
                  replace = True)



# -----------------------------------------------------------------------------
#
def optimize_fits(fpoints, fpoint_weights, mlist, dmap,
                  names = None, status_text = None,
                  optimize = True, use_threads = False, task=None):

    from time import time
    c0 = time()

    darray = dmap.data.matrix()
    xyz_to_ijk_tf = dmap.data.xyz_to_ijk_transform

    if 0 or use_threads:
        # TODO: report status messages.
        print " - in parallel!"
        fits = parallel_fitting(fpoints, fpoint_weights,
                                mlist, darray, xyz_to_ijk_tf, optimize)
    else:
        fits = []
        for i, Mi in enumerate(mlist):
            #if names:
            #    print "%d/%d : %s" % ( i+1, len(mlist), names[i] )
            if task :
                task.updateStatus ( "Fit %d/%d" % (i+1, len(mlist)) )
            if status_text:
                status ( "%s %d/%d" % (status_text, i+1, len(mlist)) )
            Mfit, corr, stats = FitMap_T(fpoints, fpoint_weights, Mi, darray, xyz_to_ijk_tf, optimize = optimize)
            #print "Fit ", i, ":", "Shift: ", stats['totShift'], "Angle:", stats['totAngle'], "height", stats['difCC'], "Final", corr
            fits.append((Mfit, corr, stats))

    c1 = time()
    print '%d fits took %.2f seconds' % (len(fits), c1-c0)

    return fits


# -----------------------------------------------------------------------------
#
def parallel_fitting(fpoints, fpoint_weights, mlist, darray, xyz_to_ijk_tf,
                     optimize = True):

    #
    # Choose number of threads to match number of cores.  Using more threads
    # creates large inefficiency (2x slower) in Python 2.7 due to context
    # switching overhead (Dave Beazley lecture).
    #
    # System usually reports twice actual number of cores due to hyperthreading.
    # Hyperthreading doesn't help if fitting tests so half that number.
    #
    import multiprocessing
    threads = multiprocessing.cpu_count()
    print 'parallel fitting using %d threads' % threads

    # Avoid periodic Python context switching.
    import sys
    original_check_interval = sys.getcheckinterval()
    sys.setcheckinterval(1000000000)

    # Define thread class for fitting.
    from threading import Thread
    class Fit_Thread(Thread):
        def __init__(self, mlist):
            Thread.__init__(self)
            self.mlist = mlist
        def run(self):
            self.fits = [FitMap_T(fpoints, fpoint_weights, m, darray,
                                  xyz_to_ijk_tf, optimize = optimize)
                         for m in self.mlist]

    # Starts threads with each calculating an equal number of fits.
    n  = len(mlist)
    g = [mlist[(n*c)/threads:(n*(c+1))/threads] for c in range(threads)]
    threads = [Fit_Thread(ml) for ml in g]
    for t in threads:
        t.start()

    # Wait for all threads to finish
    for t in threads:
        t.join()

    # Restore periodic context switching.
    sys.setcheckinterval(original_check_interval)

    # Collect fit results from all threads.
    fits = []
    for t in threads:
        for Mfit, corr, stats in t.fits:
            fits.append((Mfit, corr, stats))

    return fits


# -----------------------------------------------------------------------------
#
def FitMap_T ( fpoints, fpoint_weights, M, darray, xyz_to_ijk_transform,
               bTrans=True, bRot=True, optimize=True ) :

    xyz_to_ijk_tf = multiply_matrices(xyz_to_ijk_transform, M.tolist())

    if optimize:
        from FitMap import locate_maximum
        totShift = 0.0
        totAngle = 0.0
        map_values, outside = interpolate_volume_data(fpoints, xyz_to_ijk_tf, darray)
        initOlap, initCC = overlap_and_correlation ( fpoint_weights, map_values )

        for i in range (5) :
            move_tf, stats = locate_maximum(fpoints, fpoint_weights,
                                            darray, xyz_to_ijk_tf,
                                            max_steps = 1000,
                                            ijk_step_size_min = 0.01,
                                            ijk_step_size_max = 0.5,
                                            optimize_translation = bTrans,
                                            optimize_rotation = bRot,
                                            metric = 'sum product',
                                            request_stop_cb = None)

            xT, xR = xf_2_M ( chimera_xform ( move_tf ) )
            M = M * xT * xR
            corr = stats['correlation']

            #print ' \t%d steps: d %.3g, r %.3g, cor %f' % (stats['steps'], stats['shift'], stats['angle'], corr )

            totShift = totShift + stats['shift']
            totAngle = totAngle + stats['angle']

            if ( stats['shift'] < 0.1 and stats['angle'] < 0.1 ) :
                break

            xyz_to_ijk_tf = multiply_matrices(xyz_to_ijk_transform, M.tolist())

        stats['totAngle'] = totAngle
        stats['totShift'] = totShift
        stats['difCC'] = corr - initCC

    else:
        map_values, outside = interpolate_volume_data(fpoints, xyz_to_ijk_tf,
                                                      darray )
        olap, corr = overlap_and_correlation ( fpoint_weights, map_values )
        stats = {}

        stats['totAngle'] = 0.0
        stats['totShift'] = 0.0
        stats['difCC'] = 0.0

    return M, corr, stats


def molApplyT ( mol, T ) :

    xf = chimera.Xform.xform ( T[0,0], T[0,1], T[0,2], T[0,3], T[1,0], T[1,1], T[1,2], T[1,3], T[2,0], T[2,1], T[2,2], T[2,3], True )
    # print xf

    mol.COM = chimera.Vector (0,0,0)

    for at in mol.atoms :
        c = xf.apply ( at.coord() )
        at.setCoord ( c )
        mol.COM = mol.COM + c.toVector()

    mol.COM = mol.COM / float ( len(mol.atoms) )

#
# Change atom coordinates so the center of mass is at the origin and
# the principal axes are x, y and z.  Make a compensating transformation
# of the molecule coordinate system so the molecule does not move in the
# graphics window.
#
def centerMol ( sel_str ):

    sel = chimera.selection.OSLSelection (sel_str)
    mols = sel.molecules ()
    atoms = sel.atoms()

    if len(mols) == 0 :
        print "Failed to center molecule"
        return []

    if hasattr(mols[0], 'centered') :
        return mols

    umsg ( "Centering %d structures, %d atoms" % (len(mols), len(atoms) ) )

    points = get_atom_coordinates ( atoms, transformed = False )
    COM, U, S, V = prAxes ( points )

    # move COM to origin and align pr. axes with XYZ
    tAO = numpy.matrix ( [
        [ 1, 0, 0, -COM[0] ],
        [ 0, 1, 0, -COM[1] ],
        [ 0, 0, 1, -COM[2] ],
        [ 0, 0, 0,       1 ]  ] )

    tAR = numpy.matrix ( [
        [ V[0,0], V[0,1], V[0,2], 0 ],
        [ V[1,0], V[1,1], V[1,2], 0 ],
        [ V[2,0], V[2,1], V[2,2], 0 ],
        [      0,      0,      0, 1 ]  ] )

    # Adjust coordinate system so molecule does not appear to move.
    tf = invert_matrix((tAR*tAO).tolist()[:3])

    for fmol in mols :

        fmol.COM, fmol.U, fmol.S, fmol.V = COM, U, S, V

        print "Mol %s .%d" % (fmol.name, fmol.subid)
        molApplyT ( fmol, tAO )
        print " - COM after translation:", fmol.COM
        molApplyT ( fmol, tAR )
        print " - COM after rotation:", fmol.COM

        fmol.openState.localXform(chimera_xform(tf))

        fmol.mT = numpy.matrix ( [
            [ 1, 0, 0, 0 ],
            [ 0, 1, 0, 0 ],
            [ 0, 0, 1, 0 ],
            [ 0, 0, 0, 1 ]  ] )

        fmol.mR = numpy.matrix ( [
            [ 1, 0, 0, 0 ],
            [ 0, 1, 0, 0 ],
            [ 0, 0, 1, 0 ],
            [ 0, 0, 0, 1 ]  ] )

        fmol.M = fmol.mT * fmol.mR

    points = get_atom_coordinates ( atoms, transformed = False )

    for fmol in mols :
        fmol.COM, fmol.U, fmol.S, fmol.V = prAxes ( points )

        ppoints = points * fmol.U

        fmol.BoundRad = numpy.sqrt ( numpy.max ( numpy.sum ( numpy.square (ppoints), 1 ) ) )
        fmol.Extents = numpy.asarray ( numpy.max ( numpy.abs ( ppoints ), 0 ) )[0]

        fmol.Extents[0] = fmol.Extents[0] + 5.0
        fmol.Extents[1] = fmol.Extents[1] + 5.0
        fmol.Extents[2] = fmol.Extents[2] + 5.0

        fmol.centered = True

        umsg ( "Centered %s .%d (radius %.2fA, extents %.2fA %.2fA %.2fA)" % (
            fmol.name, fmol.subid, fmol.BoundRad, fmol.Extents[0], fmol.Extents[1], fmol.Extents[2] ) )

    return mols


def fit_points(fmap, useThreshold = True):

    mat = fmap.data.full_matrix()
    threshold = fmap.surface_levels[0]

    if useThreshold == False :
        threshold = -1e9
        print " - not using threshold"

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



def move_fit_models(fmap, M, dmap_xform):

        tXO, tXR = xf_2_M ( dmap_xform )
        T = tXO * tXR * M
        xfA = chimera.Xform.xform ( T[0,0], T[0,1], T[0,2], T[0,3], T[1,0], T[1,1], T[1,2], T[1,3], T[2,0], T[2,1], T[2,2], T[2,3] )
        fmap.openState.xform = xfA

        print "moving %d mols to fitted position" % len(fmap.mols)
        for mol in fmap.mols :
            mol.openState.xform = xfA


def principle_axes_alignments ( points, flips, preM ):

        COM, U, S, V = prAxes ( points )

        comT = numpy.matrix ( [
            [ 1, 0, 0, COM[0] ],
            [ 0, 1, 0, COM[1] ],
            [ 0, 0, 1, COM[2] ],
            [ 0, 0, 0,      1 ]  ] )

        mlist = []
        for j in range( len(flips) ) :

            af = flips[j]

            mR = numpy.matrix ( [
                [ af[0]*U[0,0], af[1]*U[0,1], af[2]*U[0,2], 0 ],
                [ af[0]*U[1,0], af[1]*U[1,1], af[2]*U[1,2], 0 ],
                [ af[0]*U[2,0], af[1]*U[2,1], af[2]*U[2,2], 0 ],
                [            0,            0,            0, 1 ] ] )

            M = comT * mR * preM
            mlist.append(M)

        return mlist

#
# Return list of rotation xforms uniformly distributed rotating about
# N axis vectors and M angles about each axis.
#
# http://www.math.niu.edu/~rusin/known-math/97/spherefaq
#
def uniform_rotation_angles(N, M) :

    thetas, phis = [], []
    from math import acos, sin, cos, sqrt, pi
    for k in range ( 1, N+1 ) :
        h = -1.0 + ( 2.0*float(k-1)/float(N-1) )
        phis.append ( acos(h) )
        thetas.append ( 0 if k == 1 or k == N else
                        (thetas[k-2] + 3.6/sqrt(N*(1.0-h**2.0))) % (2*pi) )

    ralist = []
    for theta, phi in zip(thetas, phis):
        for m in range ( M ) :
            rot = 2*pi*float(m)/float(M)
            ralist.append((theta,phi,rot))

    return ralist


def rotation_from_angles(theta, phi, rot) :

    from math import sin, cos, pi
    v = chimera.Vector (sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi))
    xfR = chimera.Xform.rotation ( v, rot*180/pi )
    Mt, Mr = xf_2_M ( xfR )
    return Mr


def xf_2_M (xf) :

    X = ( numpy.matrix (xf.getOpenGLMatrix()) ).reshape([4,4]).transpose()

    tXO = numpy.matrix ( [
        [ 1, 0, 0, X[0,3] ],
        [ 0, 1, 0, X[1,3] ],
        [ 0, 0, 1, X[2,3] ],
        [ 0, 0, 0,      1 ]  ] )

    tXR = numpy.matrix ( [
        [ X[0,0], X[0,1], X[0,2], 0 ],
        [ X[1,0], X[1,1], X[1,2], 0 ],
        [ X[2,0], X[2,1], X[2,2], 0 ],
        [      0,      0,      0, 1 ]  ] )

    return [tXO, tXR]


def xf_2_MM (xf) :

    X = ( numpy.matrix (xf.getOpenGLMatrix()) ).reshape([4,4]).transpose()

    tXO = numpy.matrix ( [
        [ 1, 0, 0, X[0,3] ],
        [ 0, 1, 0, X[1,3] ],
        [ 0, 0, 1, X[2,3] ],
        [ 0, 0, 0,      1 ]  ] )

    tXR = numpy.matrix ( [
        [ X[0,0], X[0,1], X[0,2], 0 ],
        [ X[1,0], X[1,1], X[1,2], 0 ],
        [ X[2,0], X[2,1], X[2,2], 0 ],
        [      0,      0,      0, 1 ]  ] )

    return tXO * tXR




def place_map_resample ( fmap, dmap, fnamesuf ) :

    # get bounds of points above threshold
    fpoints = grid_indices (fmap.data.size, numpy.single)  # i,j,k indices
    transform_vertices ( fpoints, fmap.data.ijk_to_xyz_transform )
    mat = fmap.data.full_matrix ()
    fpoint_weights = numpy.ravel(mat).astype(numpy.single)
    threshold = fmap.surface_levels[0]
    ge = numpy.greater_equal(fpoint_weights, threshold)
    fpoints = numpy.compress(ge, fpoints, 0)
    fpoint_weights = numpy.compress(ge, fpoint_weights)
    nz = numpy.nonzero( fpoint_weights )[0]
    print " - %d above %f in %s" % (len(nz), threshold, fmap.name)
    #print "points: ", fpoints
    #print "weights: ", fpoint_weights

    transform_vertices ( fpoints, Matrix.xform_matrix( fmap.openState.xform ) )
    transform_vertices ( fpoints, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
    transform_vertices ( fpoints, dmap.data.xyz_to_ijk_transform )
    #print "points in %s ref:" % dmap.name, fpoints

    bound = 2
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

    npoints = grid_indices ( (nn1, nn2, nn3), numpy.single)  # i,j,k indices
    transform_vertices ( npoints, ndata.ijk_to_xyz_transform )

    dvals = fmap.interpolated_values ( npoints, dmap.openState.xform )
    #dvals = numpy.where ( dvals > threshold, dvals, numpy.zeros_like(dvals) )
    #nze = numpy.nonzero ( dvals )

    nmat = dvals.reshape( (nn3,nn2,nn1) )
    #f_mat = fmap.data.full_matrix()
    #f_mask = numpy.where ( f_mat > fmap.surface_levels[0], numpy.ones_like(f_mat), numpy.zeros_like(f_mat) )
    #df_mat = df_mat * f_mask



    fmap_base = os.path.splitext(fmap.name)[0]
    dmap_base = os.path.splitext(dmap.name)[0]
    fmap_path = os.path.splitext (fmap.data.path)[0]
    dmap_path = os.path.splitext (dmap.data.path)[0]


    ndata = VolumeData.Array_Grid_Data ( nmat, nO, fmap.data.step, dmap.data.cell_angles, name=(dmap_base + fnamesuf) )
    try : nv = VolumeViewer.volume.add_data_set ( ndata, None )
    except : nv = VolumeViewer.volume.volume_from_grid_data ( ndata )


    nv.name = dmap_base + fnamesuf
    nv.openState.xform = dmap.openState.xform

    npath = dmap_path + fnamesuf
    nv.write_file ( npath, "mrc" )
    print "Wrote ", npath

    return nv




def CopyMol ( mol ) :

    nmol = chimera.Molecule()
    nmol.name = mol.name

    aMap = dict()
    clr = ( rand(), rand(), rand() )

    for res in mol.residues :
        nres = nmol.newResidue (res.type, chimera.MolResId(res.id.chainId, res.id.position))
        # print "New res: %s %d" % (nres.id.chainId, nres.id.position)
        for at in res.atoms :
            nat = nmol.newAtom (at.name, chimera.Element(at.element.number))
            # todo: handle alt
            aMap[at] = nat
            nres.addAtom( nat )
            nat.setCoord ( at.coord() )
            nat.drawMode = nat.Sphere
            nat.color = chimera.MaterialColor( clr[0], clr[1], clr[2], 1.0 )
            nat.display = True
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

    for bond in mol.bonds :
        nb = nmol.newBond ( aMap[bond.atoms[0]], aMap[bond.atoms[1]] )
        nb.display = nb.Smart

    return nmol


def CopyMolX ( mol, xf ) :

    nmol = chimera.Molecule()
    nmol.name = mol.name

    aMap = dict()
    from random import random as rand
    clr = ( rand(), rand(), rand() )

    for res in mol.residues :
        #nres = nmol.newResidue (res.type, chimera.MolResId(res.id.chainId, res.id.position))
        nres = nmol.newResidue (res.type, chimera.MolResId(res.id.chainId, res.id.position))
        # print "New res: %s %d" % (nres.id.chainId, nres.id.position)
        for at in res.atoms :
            nat = nmol.newAtom (at.name, chimera.Element(at.element.number))
            aMap[at] = nat
            nres.addAtom( nat )
            nat.setCoord ( xf.apply(at.xformCoord()) )
            nat.altLoc = at.altLoc
            nat.occupancy = at.occupancy
            nat.bfactor = at.bfactor
            if res.isProt or res.isNA :
                nat.display = False
            else :
                nat.display = True
                nat.radius=1.46
            nat.color = chimera.MaterialColor( clr[0], clr[1], clr[2], 1.0 )
            nat.drawMode = nat.EndCap

        nres.isHelix = res.isHelix
        nres.isHet = res.isHet
        nres.isSheet = res.isSheet
        nres.isStrand = res.isStrand
        nres.ribbonDisplay = True
        nres.ribbonDrawMode = 2
        nres.ribbonColor = chimera.MaterialColor( clr[0], clr[1], clr[2], 1.0 );

    for bond in mol.bonds :
        nb = nmol.newBond ( aMap[bond.atoms[0]], aMap[bond.atoms[1]] )
        nb.display = nb.Smart

    return nmol


def CopyChain ( mol, nmol, cid, ncid, xf  ) :

    if nmol == None :
        nmol = chimera.Molecule()
        nmol.name = mol.name + "_sym"

    aMap = dict()
    clr = ( rand(), rand(), rand() )

    from SWIM import SetBBAts
    SetBBAts ( nmol )

    ligandAtoms = []
    for res in nmol.residues :
        if not res.isProt and not res.isNA :
            for at in res.atoms :
                ligandAtoms.append ( at )

    print " %d ligats" % len(ligandAtoms)


    for res in mol.residues :

        if res.id.chainId == cid :

            isDuplicate = False
            if not res.isProt and not res.isNA :
                for at in res.atoms :
                    atP = xf.apply(at.xformCoord())
                    for ligAt in ligandAtoms :
                        v = ligAt.coord() - atP
                        if v.length < 0.2 :
                            isDuplicate = True
                            break
                    if isDuplicate :
                        break

            if isDuplicate :
                continue


            nres = nmol.newResidue (res.type, chimera.MolResId(ncid, res.id.position))
            # print "New res: %s %d" % (nres.id.chainId, nres.id.position)
            for at in res.atoms :
                nat = nmol.newAtom (at.name, chimera.Element(at.element.number))
                # todo: handle alt
                aMap[at] = nat
                nres.addAtom( nat )
                nat.setCoord ( xf.apply(at.xformCoord()) )
                nat.altLoc = at.altLoc
                nat.occupancy = at.occupancy
                nat.bfactor = at.bfactor
                if res.isProt or res.isNA :
                    nat.display = False
                else :
                    nat.display = True
                    nat.radius=1.46
                nat.color = chimera.MaterialColor( clr[0], clr[1], clr[2], 1.0 )
                nat.drawMode = nat.EndCap

            nres.isHelix = res.isHelix
            nres.isHet = res.isHet
            nres.isSheet = res.isSheet
            nres.isStrand = res.isStrand
            nres.ribbonDisplay = True
            nres.ribbonDrawMode = 2
            nres.ribbonColor = chimera.MaterialColor( clr[0], clr[1], clr[2], 1.0 );

    for bond in mol.bonds :
        at1, at2 = bond.atoms
        if at1 in aMap and at2 in aMap :
            nb = nmol.newBond ( aMap[at1], aMap[at2] )
            nb.display = nb.Smart

    return nmol





def map_overlap_and_correlation (map1, map2, above_threshold):

    import FitMap
    olap, cor = FitMap.map_overlap_and_correlation ( v1, v2, above_threshold )[:2]
    return olap, cor



def overlap_and_correlation ( v1, v2 ):

    import FitMap
    olap, cor = FitMap.overlap_and_correlation ( v1, v2 )[:2]
    return olap, cor



def reportFitRegions(map_name, regs):

    r = ', '.join(str(reg.rid) for reg in regs[:5])
    if len(regs) > 5:
        r += '...'
    umsg ( "Fitting %s to %d regions (%s)" % ( map_name, len(regs), r ) )


def getMod ( name ) :

    import chimera
    mlist = chimera.openModels.list ()
    for mol in mlist :
        if mol.name == name :
            return mol
    return None




def ShapeMatchScore ( atoms, dmap, bPrint=False ) :

    #fmol = fmap.mol
    #print "atoms from", fmol.name
    #points = get_atom_coordinates ( fmol.atoms, transformed = True )

    print "shape match of %d atoms with map %s" % (len(atoms), dmap.name)
    points = get_atom_coordinates ( atoms, transformed = True )
    transform_vertices ( points, xform_matrix ( dmap.openState.xform.inverse() ) )
    points0 = points.copy()
    transform_vertices ( points, dmap.data.xyz_to_ijk_transform )
    #print "points in %s ref:" % dmap.name, fpoints

    bound = int ( numpy.ceil (3.0 * max(dmap.data.step)) ) + 2
    print " - bound:", bound
    lo = numpy.floor ( numpy.min ( points, axis=0 ) ) - (bound, bound, bound)
    hi = numpy.ceil  ( numpy.max ( points, axis=0 ) ) + (bound, bound, bound)
    print " - min:", lo
    print " - max:", hi


    O = list ( dmap.data.origin )
    n = list ( dmap.data.size )
    s = dmap.data.step
    print " - dmap size:", n
    print " - dmap O:", O

    for i in (0,1,2) :
        if lo[i] < 0 :
            n[i] -= lo[i]
            O[i] += lo[i]*s[i]

    for i in (0,1,2) :
        if hi[i] > n[i] :
            n[i] = hi[i]

    print " - dmap size:", n
    print " - dmap O:", O

    nmat = numpy.ones ( (n[2], n[1], n[0]) )
    eps = 0.5 * numpy.sqrt ( (s[0] * s[0]) + (s[1] * s[1]) + (s[2] * s[2]) )
    ndata = VolumeData.Array_Grid_Data ( nmat, O, s, dmap.data.cell_angles )
    amap_data = VolumeData.zone_masked_grid_data ( ndata, points0, max(3.0, eps) )
    amat = amap_data.full_matrix()
    if 0 :
        amap = VolumeViewer.volume_from_grid_data ( amap_data )
        amap.name = dmap.name + "_()_" + atoms[0].molecule.name
        amap.openState.xform = dmap.openState.xform

    npoints = grid_indices ( (int(n[0]), int(n[1]), int(n[2])), numpy.single)  # i,j,k indices
    transform_vertices ( npoints, ndata.ijk_to_xyz_transform )
    dvals = dmap.interpolated_values ( npoints, dmap.openState.xform )
    #dvals = numpy.where ( dvals > threshold, dvals, numpy.zeros_like(dvals) )
    #nze = numpy.nonzero ( dvals )

    nmat = dvals.reshape( (n[2], n[1], n[0]) )
    nmatm = numpy.where ( nmat > dmap.surface_levels[0], numpy.ones_like(nmat), numpy.zeros_like(nmat) )
    #df_mat = df_mat * f_mask

    if 0 :
        ndata = VolumeData.Array_Grid_Data ( nmatm, O, s, dmap.data.cell_angles )
        nmap = VolumeViewer.volume_from_grid_data ( ndata )
        nmap.name = dmap.name + "_(2)"
        nmap.openState.xform = dmap.openState.xform

    nmatm = nmatm.astype ( numpy.int )
    amat = amat.astype ( numpy.int )
    imat = nmatm & amat
    umat = nmatm | amat

    if 0 :
        ndata = VolumeData.Array_Grid_Data ( umat, O, s, dmap.data.cell_angles )
        nmap = VolumeViewer.volume_from_grid_data ( ndata )
        nmap.name = dmap.name + "_(U)_" + atoms[0].molecule.name
        nmap.openState.xform = dmap.openState.xform

        ndata = VolumeData.Array_Grid_Data ( imat, O, s, dmap.data.cell_angles )
        nmap = VolumeViewer.volume_from_grid_data ( ndata )
        nmap.name = dmap.name + "_(I)_" + atoms[0].molecule.name
        nmap.openState.xform = dmap.openState.xform


    nz_int =  numpy.shape ( (imat).nonzero () )[1]
    nz_uni =  numpy.shape ( (umat).nonzero () )[1]

    sm_score = float(nz_int) / float (nz_uni)

    print " - intersection %d, union %d - sm: %.3f" % (nz_int, nz_uni, sm_score)

    return sm_score



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


def cc_by_residue ( fmap, dmap, w ) :

    rccs = []
    rmap = None
    rmap_pos = None
    rpoints, rpoint_weights = None, None
    if hasattr ( fmap, "mols" ) :
        for mol in fmap.mols :
            for ri, res in enumerate ( mol.residues ) :

                try :
                    cat = res.atomsMap["CA"][0]
                except :
                    continue

                xf = None
                if rmap == None :
                    rmap = makeMap ( "#%d:%d@CA" % (mol.id, res.id.position)
                                     , 16.0, 1.0, (.5, .5, .5, 1.0), "resmap" )
                    rmap_pos = cat.coord()
                    rpoints, rpoint_weights = fit_points(rmap)
                    xf = rmap.openState.xform

                else :
                    #new_rmap_pos = cat.coord()
                    d = cat.coord() - rmap_pos
                    xf = rmap.openState.xform
                    xf.multiply ( chimera.Xform.translation ( d ) )

                rmap_values = dmap.interpolated_values ( rpoints, xf )
                olap, corr = overlap_and_correlation ( rpoint_weights, rmap_values )
                #print " - overlap: %f, cross-correlation: %f" % (olap, corr)
                #chimera.openModels.close ( rmap )
                rccs.append ( corr )
                #print corr,

        fp = open ( "ff_prcc_w%d.txt" % w, "a" )
        fp.write ( "%s" % fmap.mols[0].name )
        for i, cc in enumerate ( rccs ) :
            if w == 1 :
                fp.write ( "\t%f" % cc )
            else :
                wscores = rccs [ max(i-w, 0) : min(i+w,len(rccs)) ]
                wscore = float ( sum ( wscores ) ) / float( len(wscores) )
                fp.write ( "\t%f" % wscore )
        fp.write ( "\n" )
        fp.close ()
        if rmap : chimera.openModels.close ( rmap )




def RandColorChains ( m ) :

    ct = {}
    for r in m.residues: ct[r.id.chainId] = 1
    clist = ct.keys()
    clist.sort()
    chains_clrs = {}
    cnames = ""

    for ci, cid in enumerate ( clist ) :
        clr = ( rand()*.7, rand()*.7, rand()*.7 )
        #print "- %s: clr(%.2f, %.2f, %.2f)" % (cid, clr[0], clr[1], clr[2])
        chains_clrs[cid] = chimera.MaterialColor ( clr[0], clr[1], clr[2], 1.0 )
        cnames = cnames + cid

    print "%s - color ribbon for %d chains -" % ( m.name, len(cnames) ), cnames

    # color atoms
    for r in m.residues :
        clr = chains_clrs[r.id.chainId]
        r.ribbonDrawMode = 2
        r.ribbonColor = clr
        r.ribbonDisplay = True
        for at in r.atoms :
            at.display = False
            at.color = clr

    return chains_clrs



def MapStats ( dmap, aboveZero = True ) :

    print "map: %s" % (dmap.name)

    MapDataStats ( dmap.data )


def MapDataStats ( data, aboveZero = True ) :

    mat = data.full_matrix ()

    if aboveZero :
        mat = numpy.where ( mat >= 0.0, mat, numpy.zeros_like(mat) )

    weights = mat.ravel()
    #ge = numpy.greater_equal(weights, 0.0)
    #weights = numpy.compress(ge, weights)
    #nz = numpy.nonzero( weights )[0]

    smin = numpy.min (weights)
    sdev = numpy.std (weights)
    savg = numpy.average(weights)
    smax = numpy.max(weights)

    #print " -", len(nz), " nonzero"
    print " - range: %.3f -> %.3f, avg=%.3f, sdev=%.3f" % (smin, smax, savg, sdev)



def NormalizeMap ( dmap ) :

    print "Normalizing map: %s" % (dmap.name)

    ndata = NormalizeData ( dmap.data )

    try : nv = VolumeViewer.volume.add_data_set ( ndata, None )
    except : nv = VolumeViewer.volume.volume_from_grid_data ( ndata )

    nv.name = os.path.splitext(dmap.name)[0] + "_norm.mrc"
    nv.openState.xform = dmap.openState.xform

    return nv



def NormalizeData ( data ) :

    O = data.origin
    mat = data.full_matrix ()

    mat = numpy.where ( mat >= 0.0, mat, numpy.zeros_like(mat) )

    weights = mat.ravel()
    #ge = numpy.greater_equal(weights, 0.0)
    #weights = numpy.compress(ge, weights)
    #nz = numpy.nonzero( weights )[0]

    smin = numpy.min (weights)
    sdev = numpy.std (weights)
    savg = numpy.average(weights)
    smax = numpy.max(weights)

    print " - initial - range: %.3f -> %.3f, avg=%.3f, sdev=%.3f" % (smin, smax, savg, sdev)

    #mat0 = mat0 - savg
    mat0 = mat / sdev
    #mat0 = mat / smax

    weights = mat0.ravel()
    #ge = numpy.greater_equal(weights, 0.0)
    #weights = numpy.compress(ge, weights)
    smin = numpy.min (weights)
    sdev = numpy.std (weights)
    savg = numpy.average(weights)
    smax = numpy.max(weights)

    print " - normalized - range: %.3f -> %.3f, avg=%.3f, sdev=%.3f" % (smin, smax, savg, sdev)

    return VolumeData.Array_Grid_Data ( mat0, O, data.step, data.cell_angles )


def NormalizeMat ( mat ) :

    #mat = numpy.where ( mat >= 0.0, mat, numpy.zeros_like(mat) )

    weights = mat.ravel()
    #ge = numpy.greater_equal(weights, 0.0)
    #weights = numpy.compress(ge, weights)
    #nz = numpy.nonzero( weights )[0]

    smin = numpy.min (weights)
    sdev = numpy.std (weights)
    savg = numpy.average(weights)
    smax = numpy.max(weights)

    print " - initial - range: %.3f -> %.3f, avg=%.3f, sdev=%.3f" % (smin, smax, savg, sdev)

    #mat0 = mat0 - savg
    mat0 = mat / sdev
    #mat0 = mat / smax

    weights = mat0.ravel()
    #ge = numpy.greater_equal(weights, 0.0)
    #weights = numpy.compress(ge, weights)
    smin = numpy.min (weights)
    sdev = numpy.std (weights)
    savg = numpy.average(weights)
    smax = numpy.max(weights)

    print " - normalized - range: %.3f -> %.3f, avg=%.3f, sdev=%.3f" % (smin, smax, savg, sdev)

    return mat0



def OneMinusOneMap ( dmap ) :

    thr = dmap.surface_levels[0]

    umsg ("Generating 1/-1 map for " + dmap.name + " thr: %.3f" % thr)

    m2 = None


    if 0 :
        m1 = dmap.data.full_matrix()
        m2 = numpy.where ( m1 > thr, numpy.ones_like(m1)*1, numpy.ones_like(m1)*-1.0 )
    else :

        m0 = dmap.data.full_matrix()
        inside_start = numpy.where ( m0 > thr, numpy.ones_like(m0)*1, numpy.zeros_like(m0) )
        outside_mask = numpy.where ( m0 < thr, numpy.ones_like(m0)*1, numpy.zeros_like(m0) )

        gvm = inside_start.copy();
        for i in range (numit) :
            nv_1 = numpy.roll(gvm, 1, axis=0)
            nv_2 = numpy.roll(gvm, -1, axis=0)
            nv_3 = numpy.roll(gvm, 1, axis=1)
            nv_4 = numpy.roll(gvm, -1, axis=1)
            nv_5 = numpy.roll(gvm, 1, axis=2)
            nv_6 = numpy.roll(gvm, -1, axis=2)
            gvm = 1.0/6.0 * ( nv_1 + nv_2 + nv_3 + nv_4 + nv_5 + nv_6 )
            gvm = outside_mask * gvm + inside_start


    from VolumeData import Array_Grid_Data
    mgrid = Array_Grid_Data ( m2, dmap.data.origin, dmap.data.step, dmap.data.cell_angles, name="map_one_minus_one")

    import VolumeViewer
    #nv = VolumeViewer.volume_from_grid_data ( mgrid, show_data = False, show_dialog = False )
    return VolumeViewer.volume_from_grid_data ( mgrid )



def MapStats ( dmap, aboveZero = False ) :

    print "Map Stats: %s" % (dmap.name)

    mat = dmap.data.full_matrix ()

    if aboveZero :
        mat = numpy.where ( mat >= 0.0, mat, numpy.zeros_like(mat) )

    weights = mat.ravel()
    #ge = numpy.greater_equal(weights, 0.0)
    #weights = numpy.compress(ge, weights)
    #nz = numpy.nonzero( weights )[0]

    smin = numpy.min (weights)
    sdev = numpy.std (weights)
    savg = numpy.average(weights)
    smax = numpy.max(weights)

    #print " -", len(nz), " nonzero"
    print " - range: %.3f -> %.3f, avg=%.3f, sdev=%.3f" % (smin, smax, savg, sdev)

    return savg, sdev



def AddNoiseToMap ( mv, mean, stdev ) :

    print "\n---adding noise mean:",mean, " stdev:", stdev, "---\n"

    nvm = mv.full_matrix()
    #f_mask = numpy.where ( nvm > 0, numpy.zeros_like(nvm), numpy.ones_like(nvm) )

    from numpy.random import standard_normal as srand
    s=mv.data.size

    noisem = srand ( (s[2],s[1],s[0]) ) * stdev - (numpy.ones_like(nvm) * mean)
    ngvm = noisem + nvm

    ndata = VolumeData.Array_Grid_Data ( ngvm, mv.data.origin, mv.data.step, mv.data.cell_angles )
    try : nvg = VolumeViewer.volume.add_data_set ( ndata, None )
    except : nvg = VolumeViewer.volume.volume_from_grid_data ( ndata )
    nvg.name = mv.name

    chimera.openModels.close ( [mv] )
    return nvg
