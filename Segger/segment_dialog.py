
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
import graph; reload(graph)
from Segger import dev_menus, timing, seggerVersion, showDevTools

#dev_menus = False
#showDevTools = True




OML = chimera.openModels.list

REG_OPACITY = 0.45


def umsg ( txt ) :
    print txt
    status ( txt )

def status ( txt ) :
    txt = txt.rstrip('\n')
    msg.configure(text = txt)
    msg.update_idletasks()


class Volume_Segmentation_Dialog ( chimera.baseDialog.ModelessDialog ):

    title = "Segger (v" + seggerVersion + ")"
    name = "segment map"
    #buttons = ('Segment', 'Group', 'Ungroup', 'Options', 'Shortcuts', "Tools", "Close")
    buttons = ('Options', 'Shortcuts', "Tools", "Log", "Close")
    help = 'https://github.com/gregdp/segger'

    def fillInUI(self, parent):

        self.group_mouse_mode = None

        tw = parent.winfo_toplevel()
        self.toplevel_widget = tw
        tw.withdraw()

        parent.columnconfigure(0, weight = 1)

        row = 1

        menubar = Tkinter.Menu(parent, type = 'menubar', tearoff = False)
        tw.config(menu = menubar)

        file_menu_entries = (
            ('Open segmentation...', self.OpenSegmentation),
            ('Save segmentation', self.SaveSegmentation),
            ('Save segmentation as...', self.SaveSegmentationAs),
            ("Save selected regions to .mrc file...", self.WriteSelRegionsMRCFile),
            ("Save all regions to .mrc file...", self.WriteAllRegionsMRCFile),
            ("Save each region to .mrc file...", self.WriteEachRegionMRCFile),
            ("Close segmentation", self.CloseSeg),
            ("Close all segmentations except displayed", self.CloseHiddenSeg),
            ("Close all segmentations", self.CloseAll),
            ("Associate Selected", self.Associate),
            )

        fmenu = Hybrid.cascade_menu(menubar, 'File', file_menu_entries)

        import attributes
        regions_menu_entries = (
            'separator',
            ("Show all", self.RegSurfsShowAll),
            ("Show only selected", self.RegSurfsShowOnlySelected),
            ("Show adjacent", self.RegSurfsShowAdjacent),
            ('Show grouping', self.ShowUngroupedSurfaces),
            ('Unshow grouping', self.ShowGroupSurfaces),
            ("Hide", self.RegSurfsHide),
            ("Make transparent", self.RegSurfsTransparent),
            ("Make opaque", self.RegSurfsOpaque),
            ('Color density map', self.ColorDensity),
            'separator',
            ('Select groups', self.SelectGroups),
            ('Select boundary regions', self.SelectBoundaryRegions),
            ("Invert selection", self.Invert),
            ("Regions overlapping current selection", self.Overlapping),
            'separator',
            #("Group selected", self.JoinSelRegs),
            #("Ungroup selected", self.UngroupSelRegs),
            #("Smooth and group", self.SmoothAndGroupOneStep),
            ("Delete selected regions", self.DelSelRegs),
            ("Delete all except selected", self.DelExcSelRegs),
            "separator",
            ("Enclosed volume", self.RegionsVolume),
            ("Mean and SD", self.RegionMeanAndSD),
            #("Mask map with selected", self.MaskMapWRegions),
            #("Mask another map with selected (shrink map)", self.MaskAnotherMapWRegionsShrink),
            #("Mask another map with selected (keep map dimensions)", self.MaskAnotherMapWRegions),
            ("Extract densities...", self.ExtractDensities),
            ("Subtract selected from map", self.SubtractRegionsFromMap),
            ("Show axes for selected", self.ShowRegionAxesSelected),
            ("Hide all axes", self.HideRegionAxes),
            "separator",
            ("Attributes table...", attributes.show_region_attributes_dialog),
            ("How many sub-regions", self.ShowNumSubRegs)
            )

        self.regsVisMode = Tkinter.StringVar()
        self.regsVisMode.set ( 'Voxel_Surfaces' )

        rmenu = Hybrid.cascade_menu(menubar, 'Regions', regions_menu_entries)

        if dev_menus:
            rmenu.add_separator()
            for lbl, var, val, cmd in (
                ("Surfaces around voxels", self.regsVisMode,
                 'Voxel_Surfaces', self.RegsDispUpdate),
                ("Map iso-surfaces", self.regsVisMode,
                 'Iso_Surfaces', self.RegsDispUpdate),
                #(" - delete files", self.deleteCloseFiles, 1, None),
                ):
                rmenu.add_radiobutton(label = lbl, variable = var, value = val,
                                      command = cmd)
            #rmenu.add_separator()
            #for lbl, cmd in (#("Adjacency graph",self.RegionsAdjGraph),
                #("Group connected", self.GroupConnectedRegs),
                #("Select non-placed", self.SelectNonPlacedRegions),
                #("Apply threshold", self.RegsDispThr),
                #("Group connected", self.GroupConnectedRegs),
                #("Group by contacts", self.GroupByContacts),
                #("Group using all fits", self.GroupUsingFits),
                #("Ungroup ALL", self.UngroupAllRegs),
                #("Reduce map", self.ReduceMap),
                #("Close", self.CloseRegions),
                #):
                #rmenu.add_command(label = lbl, command = cmd)
            #self.deleteCloseFiles = Tkinter.IntVar()
            #self.deleteCloseFiles.set ( 1 )

            #rmenu.add_separator()
            #for lbl, cmd in (
            #    ("Mask map with selected (Cube Result)", self.MaskMapWRegionsCube),
            #    ("Extract map cube around selected", self.ExtractMapWRegionsCube),
            #    ):
            #    rmenu.add_command(label = lbl, command = cmd)

        if dev_menus:
            graph_menu_entries = (
                #('Save graph', self.SaveGraph),
                #('Load graph', self.LoadGraph),
                #'separator',
                ('Create graph with uniform link radii', self.Graph),
                #(' - Use maximum density between regions', self.GraphMaxD),
                #(' - Use average density between regions', self.GraphAvgD),
                (' - Use number of contacting voxel pairs', self.GraphN),
                'separator',
                ('Close graph', self.CloseGraph),
                #('Break selected links', graph.break_selected_links),
                #('Link selected', graph.link_selected),
                #('Group regions using skeleton', self.GroupBySkeleton),
                #('Join selected', self.GroupBySkeleton)
                )
            smenu = Hybrid.cascade_menu(menubar, 'Graph',
                                        graph_menu_entries)

        self.UseAllMods = Tkinter.IntVar()
        self.UseAllMods.set ( 0 )

        from chimera.tkgui import aquaMenuBar
        aquaMenuBar(menubar, parent, row = 0, columnspan=3)

        #umsg ( '')
        f = Tkinter.Frame(parent)
        f.grid(column=0, row=row, sticky='ew')
        row += 1

        l = Tkinter.Label(f, text=' Map:')
        l.grid(column=0, row=0, sticky='w')

        self.cur_dmap = None
        self.dmap = Tkinter.StringVar(parent)

        self.mb  = Tkinter.Menubutton ( f, textvariable=self.dmap, relief=Tkinter.RAISED )
        self.mb.grid (column=1, row=0, sticky='we', padx=5)
        self.mb.menu  =  Tkinter.Menu ( self.mb, tearoff=0, postcommand=self.MapMenu )
        self.mb["menu"]  =  self.mb.menu

        if 1:
            b = Tkinter.Button(f, text="Center", command=self.MapCOM)
            b.grid (column=6, row=0, sticky='w', padx=5)

        b = Tkinter.Button(f, text="Segment", command=self.Segment)
        b.grid (column=7, row=0, sticky='w', padx=0)


        f = Tkinter.Frame(parent)
        f.grid(column=0, row=row, sticky='ew')
        row += 1

        l = Tkinter.Label(f, text=' Segmentation:')
        l.grid(column=0, row=0, sticky='w')

        self.cur_seg = None
        self.regions_file = Tkinter.StringVar(parent)

        rmb  = Tkinter.Menubutton ( f, textvariable=self.regions_file,
                                    relief=Tkinter.RAISED )
        rmb.grid (column=1, row=0, sticky='we', padx=5)
        rmb.menu  =  Tkinter.Menu ( rmb, tearoff=0,
                                    postcommand=self.FillSegmentationMenu )
        self.mbSegmentationMenu = rmb.menu
        rmb["menu"]  = rmb.menu

        rc = Tkinter.Label(f, text='')
        rc.grid(column=2, row=0, sticky='w')
        self.regionCount = rc

        b = Tkinter.Button(f, text="Group", command=self.Group)
        b.grid (column=7, row=0, sticky='w', padx=0)

        b = Tkinter.Button(f, text="Ungroup", command=self.Ungroup)
        b.grid (column=8, row=0, sticky='w', padx=0)




        if dev_menus:
            cp = Hybrid.Popup_Panel(parent)
            cpf = cp.frame
            cpf.grid(row = row, column = 0, sticky = 'news')
            cpf.grid_remove()
            cpf.columnconfigure(0, weight=1)
            self.contactsPanel = cp.panel_shown_variable
            row += 1
            orow = 0

            cb = cp.make_close_button(cpf)
            cb.grid(row = orow, column = 0, sticky = 'e')

            l = Tkinter.Label(cpf, text='Contact Grouping', font = 'TkCaptionFont')
            l.grid(column=0, row=orow, sticky='w', pady=5)
            orow += 1

            s = Hybrid.Scale(cpf, 'Coloring level ', 1, 100, 1, 100)
            s.frame.grid(row = orow, column = 0, sticky = 'ew')
            orow += 1
            self.colorLevel = s
            s.callback(self.SetColorLevel)

            b = Tkinter.Button(cpf, text = 'Set Grouping',
                               command = self.SetContactGrouping)
            b.grid(row = orow, column = 0, sticky = 'w')
            orow += 1



        # --- Options Frame ----------------------------------------------------------------------

        op = Hybrid.Popup_Panel(parent)
        opf = op.frame
        opf.grid(row = row, column = 0, sticky = 'news')
        opf.grid_remove()
        opf.columnconfigure(0, weight=1)
        self.optionsPanel = op.panel_shown_variable
        row += 1
        orow = 0

        dummyFrame = Tkinter.Frame(opf, relief='groove', borderwidth=1)
        Tkinter.Frame(dummyFrame).pack()
        dummyFrame.grid(row=orow,column=0,columnspan=7, pady=1, sticky='we')
        orow += 1

        cb = op.make_close_button(opf)
        cb.grid(row = orow, column = 1, sticky = 'e')

        l = Tkinter.Label(opf, text='Segmenting Options', font = 'TkCaptionFont')
        l.grid(column=0, row=orow, sticky='w', pady=5)
        orow += 1

        sopt = Tkinter.Frame(opf)
        sopt.grid(column=0, row=orow, sticky='ew', padx=10)
        orow += 1

        sorow = 0


        if 1 :
            f = Tkinter.Frame(sopt)
            f.grid(column=0, row=sorow, sticky='w')
            sorow += 1

            l = Tkinter.Label(f, text='Display at most')
            l.grid(column=0, row=0, sticky='w')

            self.maxNumRegions = Tkinter.StringVar(sopt)
            self.maxNumRegions.set ( '6000' )
            e = Tkinter.Entry(f, width=5, textvariable=self.maxNumRegions)
            e.grid(column=1, row=0, sticky='w', padx=5)
            e.bind('<KeyPress-Return>', self.NewMaxRegions)

            l = Tkinter.Label(f, text='regions, granularity')
            l.grid(column=2, row=0, sticky='w')


            self.surfaceGranularity = Tkinter.StringVar(sopt)
            self.surfaceGranularity.set ( '1' )
            e = Tkinter.Entry(f, width=5, textvariable=self.surfaceGranularity)
            e.bind('<KeyPress-Return>', self.NewSurfaceResolution)
            e.grid(column=3, row=0, sticky='w', padx=5)

            l = Tkinter.Label(f, text='voxels')
            l.grid(column=4, row=0, sticky='w')



        if 1 :

            f = Tkinter.Frame(sopt)
            f.grid(column=0, row=sorow, sticky='w')
            sorow += 1

            l = Tkinter.Label(f, text="Keep regions having at least")
            l.grid(column=0, row=0, sticky='w')

            self.minRegionSize = Tkinter.StringVar(parent)
            self.minRegionSize.set ( '1' )
            e = Tkinter.Entry(f, width=5, textvariable=self.minRegionSize)
            e.grid(column=1, row=0, sticky='w', padx=5)
            e.bind('<KeyPress-Return>', lambda e: self.RemoveSmallRegions())

            l = Tkinter.Label(f, text='voxels')
            l.grid(column=2, row=0, sticky='w')

            l = Tkinter.Label(f, text=', ')
            l.grid(column=4, row=0, sticky='w')

            self.minContactSize = Tkinter.StringVar(parent)
            self.minContactSize.set ( '0' )
            e = Tkinter.Entry(f, width=5, textvariable=self.minContactSize)
            e.grid(column=5, row=0, sticky='w', padx=5)
            e.bind('<KeyPress-Return>', lambda e: self.RemoveContactRegions())

            l = Tkinter.Label(f, text='contact voxels')
            l.grid(column=6, row=0, sticky='w')


        if 1 :

            f = Tkinter.Frame(sopt)
            f.grid(column=0, row=sorow, sticky='w')
            sorow += 1

            self.groupMode = Tkinter.StringVar()
            self.groupMode.set ( 'cons' )


            c = Tkinter.Radiobutton(f, text="Group by Smoothing", variable=self.groupMode, value = 'smooth')
            c.grid (column=0, row=0, sticky='w')

            self.numSteps = Tkinter.StringVar(sopt)
            self.numSteps.set ( '4' )
            e = Tkinter.Entry(f, width=5, textvariable=self.numSteps)
            e.grid(column=1, row=0, sticky='w', padx=5)

            l = Tkinter.Label(f, text='steps of size')
            l.grid(column=2, row=0, sticky='w')

            self.stepSize = ss = Tkinter.StringVar(sopt)
            ss.set('1')
            e = Tkinter.Entry(f, width=2, textvariable=self.stepSize)
            e.grid(column=3, row=0, sticky='w', padx=5)

            l = Tkinter.Label(f, text=', stop at')
            l.grid(column=4, row=0, sticky='w')

            self.targNRegions = Tkinter.StringVar(sopt)
            self.targNRegions.set ( '1' )
            e = Tkinter.Entry(f, width=2, textvariable=self.targNRegions)
            e.grid(column=5, row=0, sticky='w', padx=5)

            l = Tkinter.Label(f, text='regions')
            l.grid(column=6, row=0, sticky='w')


        if 1 :

            f = Tkinter.Frame(sopt)
            f.grid(column=0, row=sorow, sticky='w')
            sorow += 1

            c = Tkinter.Radiobutton(f, text="Group by Connectivity", variable=self.groupMode, value = 'cons')
            c.grid (column=0, row=0, sticky='w')

            self.numStepsCon = Tkinter.StringVar(sopt)
            self.numStepsCon.set ( '20' )
            e = Tkinter.Entry(f, width=5, textvariable=self.numStepsCon)
            e.grid(column=1, row=0, sticky='w', padx=5)

            l = Tkinter.Label(f, text='steps, stop at')
            l.grid(column=4, row=0, sticky='w')

            self.targNRegionsCon = Tkinter.StringVar(sopt)
            self.targNRegionsCon.set ( '1' )
            e = Tkinter.Entry(f, width=2, textvariable=self.targNRegionsCon)
            e.grid(column=5, row=0, sticky='w', padx=5)

            l = Tkinter.Label(f, text='regions')
            l.grid(column=6, row=0, sticky='w')

            oft = Hybrid.Checkbutton(f, 'only visible', False )
            oft.button.grid(column=7, row=0, sticky='w')
            self.groupByConsOnlyVis = oft.variable


        if 0:
            f = Tkinter.Frame(sopt)
            f.grid(column=0, row=sorow, sticky='w')
            sorow += 1

            l = Tkinter.Label(f, text='Minimum connecting voxels ')
            l.grid(column=0, row=0, sticky='w')

            self.minConnection = Tkinter.StringVar(sopt)
            self.minConnection.set ( '0' )
            e = Tkinter.Entry(f, width=5, textvariable=self.minConnection)
            e.grid(column=1, row=0, sticky='w', padx=5)


            mmf = Tkinter.Frame(sopt)
            mmf.grid(row = sorow, column = 0, sticky = 'ew')
            sorow += 1

            mg = Hybrid.Checkbutton(mmf, 'Group with mouse ', False)
            mg.button.grid(row = 0, column = 0, sticky = 'w')
            self.mouse_group = mg.variable
            mg.callback(self.mouse_group_cb)

            mgb = Hybrid.Option_Menu(mmf, '', 'button 1', 'button 2',
                                     'button 3', 'ctrl button 1',
                                     'ctrl button 2', 'ctrl button 3')
            mgb.variable.set('button 3')
            mgb.frame.grid(row = 0, column = 1, sticky = 'w')
            mgb.add_callback(self.mouse_group_button_cb)
            self.mouse_group_button = mgb



        if 1 :
            f = Tkinter.Frame(sopt)
            f.grid(column=0, row=sorow, sticky='w')

            oft = Hybrid.Checkbutton(f, '', False )
            #oft.button.grid(row = 0, column = 0, sticky = 'w')
            self.useSymmetry = oft.variable
            oft.button.grid(column=0, row=0, sticky='w')

            l = Tkinter.Label(f, text='Symmetry:')
            l.grid(column=1, row=0, sticky='w')

            self.symmetryString = Tkinter.StringVar(f)
            e = Tkinter.Entry(f, width=15, textvariable=self.symmetryString)
            e.grid(column=2, row=0, sticky='w', padx=0)

            b = Tkinter.Button(f, text="Detect", command=self.DetectSym)
            b.grid (column=3, row=0, sticky='w', padx=0)

            b = Tkinter.Button(f, text="Show Sel", command=self.ShowSelSymm)
            b.grid (column=4, row=0, sticky='w', padx=0)

            l = Tkinter.Label(f, text='Color:')
            l.grid(column=5, row=0, sticky='w')

            b = Tkinter.Button(f, text="Same", command=self.ColorSymmSame)
            b.grid (column=6, row=0, sticky='w', padx=0)

            b = Tkinter.Button(f, text="Diff", command=self.ColorSymmDiff)
            b.grid (column=7, row=0, sticky='w', padx=0)

            #sorow += 1



        # --- Shortcuts Frame ----------------------------------------------------------------------

        sc = Hybrid.Popup_Panel(parent)
        scf = sc.frame
        scf.grid(row = row, column = 0, sticky = 'news')
        scf.grid_remove()
        scf.columnconfigure(0, weight=1)
        self.shortcutsPanelShownVar = sc.panel_shown_variable
        row += 1
        orow = 0

        dummyFrame = Tkinter.Frame(scf, relief='groove', borderwidth=1)
        Tkinter.Frame(dummyFrame).pack()
        dummyFrame.grid(row=orow,column=0,columnspan=7, pady=1, sticky='we')
        orow += 1

        cb = sc.make_close_button(scf)
        cb.grid(row = orow, column = 1, sticky = 'e')

        l = Tkinter.Label(scf, text='Shortcuts', font = 'TkCaptionFont')
        l.grid(column=0, row=orow, sticky='w', pady=5)
        orow += 1

        sopt = Tkinter.Frame(scf)
        sopt.grid(column=0, row=orow, sticky='ew', padx=10)
        orow += 1
        sorow = 0


        f = Tkinter.Frame(sopt)
        f.grid(column=0, row=sorow, sticky='w')
        sorow += 1
        if 1 :
            b = Tkinter.Label(f, text="Select regions: ", width=15, anchor=Tkinter.E)
            b.grid (column=0, row=0, sticky='w', padx=0)

            self.overlappingPercentage = Tkinter.StringVar(f)
            self.overlappingPercentage.set ( "50" )
            #e = Tkinter.Entry(f, width=4, textvariable=self.overlappingPercentage)
            #e.grid(column=1, row=0, sticky='w', padx=5)
            #l = Tkinter.Label(f, text="%")
            #l.grid(column=2, row=0, sticky='w')


            b = Tkinter.Button(f, text="All", command=self.SelectAllRegions)
            b.grid (column=1, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="Over Sel. Atoms", command=self.Overlapping)
            b.grid (column=2, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="Invert", command=self.Invert)
            b.grid (column=3, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="Not-", command=self.SelectNotGrouped)
            b.grid (column=4, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="Grouped", command=self.SelectGrouped)
            b.grid (column=5, row=0, sticky='w', padx=5)


        f = Tkinter.Frame(sopt)
        f.grid(column=0, row=sorow, sticky='w')
        sorow += 1
        if 1 :
            l = Tkinter.Label(f, text='Selected regions: ', width=15, anchor=Tkinter.E)
            l.grid(column=0, row=0)

            #b = Tkinter.Button(f, text="Group", command=self.Group)
            #b.grid (column=1, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(f, text="Ungroup", command=self.Ungroup)
            #b.grid (column=2, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="Hide", command=self.RegSurfsHide)
            b.grid (column=1, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="Show", command=self.RegSurfsShow)
            b.grid (column=2, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="Delete", command=self.DelSelRegs)
            b.grid (column=3, row=0, sticky='w', padx=5)


        #f = Tkinter.Frame(sopt)
        #f.grid(column=0, row=sorow, sticky='w')
        #sorow += 1
        #if 1 :
            #l = Tkinter.Label(f, text=' ', width=15)
            #l.grid(column=0, row=0, sticky='w')

            b = Tkinter.Button(f, text="Tr.", command=self.RegSurfsTransparent)
            b.grid (column=4, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="Opaque", command=self.RegSurfsOpaque)
            b.grid (column=5, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="Mesh", command=self.RegSurfsMesh)
            b.grid (column=6, row=0, sticky='w', padx=5)


            #b = Tkinter.Button(f, text="Invert Selection", command=self.Invert)
            #b.grid (column=3, row=0, sticky='w', padx=5)


        f = Tkinter.Frame(sopt)
        f.grid(column=0, row=sorow, sticky='w')
        sorow += 1
        if 1 :
            l = Tkinter.Label(f, text='Show regions: ', width=15, anchor=Tkinter.E)
            l.grid(column=0, row=0, sticky='w')

            b = Tkinter.Button(f, text="None", command=self.RegSurfsShowNone)
            b.grid (column=1, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="All", command=self.RegSurfsShowAll)
            b.grid (column=2, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="Only Sel.", command=self.RegSurfsShowOnlySelected)
            b.grid (column=3, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="Adj.", command=self.RegSurfsShowAdjacent)
            b.grid (column=4, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="Not-", command=self.RegSurfsShowNotGrouped)
            b.grid (column=5, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="Grouped", command=self.RegSurfsShowGrouped)
            b.grid (column=6, row=0, sticky='w', padx=5)


            if 0 :
                b = Tkinter.Button(f, text="Axes", command=self.ShowRegionAxesSelected)
                b.grid (column=4, row=0, sticky='w', padx=5)

                self.axesFactor = Tkinter.StringVar(f)
                self.axesFactor.set ( "3" )
                e = Tkinter.Entry(f, width=4, textvariable=self.axesFactor)
                e.grid(column=5, row=0, sticky='w', padx=5)

        # ---------- end of shortcuts frame  ------------------------------------------------------



        # --- Tools Frame ----------------------------------------------------------------------

        sc = Hybrid.Popup_Panel(parent)
        scf = sc.frame
        scf.grid(row = row, column = 0, sticky = 'news')
        scf.grid_remove()
        scf.columnconfigure(0, weight=1)
        self.toolsPanelShownVar = sc.panel_shown_variable
        row += 1
        orow = 0

        dummyFrame = Tkinter.Frame(scf, relief='groove', borderwidth=1)
        Tkinter.Frame(dummyFrame).pack()
        dummyFrame.grid(row=orow,column=0,columnspan=7, pady=1, sticky='we')
        orow += 1

        cb = sc.make_close_button(scf)
        cb.grid(row = orow, column = 1, sticky = 'e')

        l = Tkinter.Label(scf, text=' Tools', font = 'TkCaptionFont')
        l.grid(column=0, row=orow, sticky='w', pady=1)
        orow += 1

        sopt = scf
        #sopt = Tkinter.Frame(scf)
        #sopt.grid(column=0, row=orow, sticky='ew', padx=10)
        #orow += 1
        sorow = orow


        if 1 :
            # flat, groove, raised, ridge, solid, or sunken
            dummyFrame = Tkinter.Frame(sopt, relief='flat', borderwidth=1)
            Tkinter.Frame(dummyFrame).pack()
            dummyFrame.grid(row=sorow,column=0,columnspan=7, pady=3, sticky='we')
            #sorow += 1

            f = Tkinter.Frame(sopt)
            f.grid(column=0, row=sorow, sticky='w')
            sorow += 1

            l = Tkinter.Label(f, text=' ', width=3, anchor=Tkinter.E)
            l.grid(column=0, row=0)

            b = Tkinter.Button(f, text="Extract", command=self.ExtractDensities)
            b.grid (column=1, row=0, sticky='w', padx=2)

            b = Tkinter.Button(f, text="SegFit", command=self.FitDialog)
            b.grid (column=2, row=0, sticky='w', padx=2)

            b = Tkinter.Button(f, text="rSeg", command=self.RSeg)
            b.grid (column=3, row=0, sticky='w', padx=2)

            b = Tkinter.Button(f, text="iSeg", command=self.ISeg)
            b.grid (column=4, row=0, sticky='w', padx=2)

            #b = Tkinter.Button(f, text="SegLoop", command=self.SegLoop)
            #b.grid (column=5, row=0, sticky='w', padx=2)

            b = Tkinter.Button(f, text="ProMod", command=self.ProMod)
            b.grid (column=5, row=0, sticky='w', padx=2)

            #b = Tkinter.Button(f, text="ModelZ", command=self.ModelZ)
            #b.grid (column=7, row=0, sticky='w', padx=2)

            b = Tkinter.Button(f, text="Movie", command=self.BioMovie)
            b.grid (column=6, row=0, sticky='w', padx=2)

            b = Tkinter.Button(f, text="MapQ", command=self.MapQ)
            b.grid (column=7, row=0, sticky='w', padx=2)

            b = Tkinter.Button(f, text="SWIM", command=self.SWIM)
            b.grid (column=8, row=0, sticky='w', padx=2)



            #b = Tkinter.Button(f, text="Frk", command=self.Frankensteinify)
            #b.grid (column=9, row=0, sticky='w', padx=2)



        if showDevTools or dev_menus :
            f = Tkinter.Frame(sopt)
            f.grid(column=0, row=sorow, sticky='w')
            sorow += 1

            l = Tkinter.Label(f, text=' ', width=15, anchor=Tkinter.E)
            l.grid(column=0, row=0)

            #b = Tkinter.Button(f, text="GeoSeg", command=self.GeoSegDialog)
            #b.grid (column=1, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(f, text="SSE", command=self.SSE)
            #b.grid (column=2, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="Animate", command=self.Animate)
            b.grid (column=3, row=0, sticky='w', padx=5)


            b = Tkinter.Button(f, text="FlexFit", command=self.FlexFit)
            b.grid (column=5, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="Geo", command=self.GeoSeg)
            b.grid (column=6, row=0, sticky='w', padx=5)


            b = Tkinter.Button(f, text="M", command=self.MDFF)
            b.grid (column=9, row=0, sticky='w', padx=2)


            f = Tkinter.Frame(sopt)
            f.grid(column=0, row=sorow, sticky='w')
            sorow += 1

            l = Tkinter.Label(f, text=' ', width=15, anchor=Tkinter.E)
            l.grid(column=0, row=0)

            b = Tkinter.Button(f, text="Mod", command=self.SegMod)
            b.grid (column=1, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="NA", command=self.SegNA)
            b.grid (column=2, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="Ar", command=self.Ar)
            b.grid (column=3, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="VR", command=self.Vr)
            b.grid (column=4, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="Mono", command=self.CamMono)
            b.grid (column=5, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="SBS", command=self.CamSBS)
            b.grid (column=6, row=0, sticky='w', padx=5)

            b = Tkinter.Button(f, text="Spr", command=self.Spr)
            b.grid (column=7, row=0, sticky='w', padx=5)


        # ---------- end of tools frame  ------------------------------------------------------


        dummyFrame = Tkinter.Frame(parent, relief='groove', borderwidth=1)
        Tkinter.Frame(dummyFrame).pack()
        dummyFrame.grid(row=row,column=0,columnspan=7, pady=7, sticky='we')


        row += 1

        f = Tkinter.Frame(parent)
        f.grid(column=0, row=row, sticky='ew')
        row += 1

        l = Tkinter.Label(f, text='To cite Segger or learn more about it press the Help button', fg="blue")
        l.grid(column=0, row=0, sticky='w')

        dummyFrame = Tkinter.Frame(parent, relief='groove', borderwidth=1)
        Tkinter.Frame(dummyFrame).pack()
        dummyFrame.grid(row=row,column=0,columnspan=7, pady=3, sticky='we')
        row += 1

        global msg
        msg = Tkinter.Label(parent, width = 60, anchor = 'w', justify = 'left', fg="red")
        msg.grid(column=0, row=row, sticky='ew')
        self.msg = msg

        umsg ( 'Select an open density map in the field above and press Segment!' )
        row += 1

        vlist = VolumeViewer.volume_list()
        if vlist:
            self.SetMapMenu(vlist[0])

        for m in regions.segmentations() :

            v = m.volume_data()
            if v and m.display :
                self.SetCurrentSegmentation ( m )
                try : self.SetMapMenu ( v )
                except : pass

        chimera.openModels.addRemoveHandler(self.ModelClosed, None)

        if dev_menus :
            self.optionsPanel.set(True)
            self.shortcutsPanelShownVar.set(True)
            self.toolsPanelShownVar.set(True)


    def SetColorLevel(self):

        smod = self.CurrentSegmentation()
        if smod is None:
            return

        s = self.colorLevel
        lev = int(s.value())

        regions = [r for r in smod.all_regions()
                   if hasattr(r, 'color_level') and r.color_level >= lev and
                   (r.preg is None or r.preg.color_level < lev)]
        smod.color_density(regions)
        return

        # TODO: Unused code adjusts region surface colors.

        if not hasattr(smod, 'contact_grouping'):
            cg = smod.contact_grouping = regions.contact_grouping(smod)
            smod.region_count = len(smod.childless_regions())
        cg = smod.contact_grouping
        range = (smod.region_count - len(cg), smod.region_count)
        if s.range() != range:
            s.set_range(range[0], range[1], step = 1)
        p = max(0, smod.region_count - int(s.value()))

        cs, pairs = regions.connected_subsets(cg[:p])

        # Reset colors
        for r in smod.regions:
            sp = r.surface_piece
            if sp:
                sp.color = r.color

        # Color groups
        for rlist in cs:
            r0 = rlist[0]
            sp0 = r0.surface_piece
            if sp0:
                c = sp0.color
                for r in rlist[1:]:
                    sp = r.surface_piece
                    if sp:
                        sp.color = c

    def SetColorLevelRange(self):

        smod = self.CurrentSegmentation()
        if smod is None:
            return

        clevels = [r.color_level for r in smod.all_regions()
                   if hasattr(r, 'color_level')]
        clmin = min(clevels)
        clmax = max(clevels)
        cl = self.colorLevel
        cl.set_range(clmin, clmax, step = 1)
        cl.set_value(clmin)

    def SetContactGrouping(self):

        smod = self.CurrentSegmentation()
        if smod is None:
            return

        s = self.colorLevel
        lev = int(s.value())

        regions = [r for r in smod.all_regions()
                   if hasattr(r, 'color_level') and r.color_level < lev]
        smod.remove_regions(regions, update_surfaces = True)
        self.RegsDispUpdate()


    def ColorDensity(self):

        smod = self.CurrentSegmentation()
        if smod is None:
            return

        smod.color_density()




    def Options(self) :
        self.optionsPanel.set (not self.optionsPanel.get())

    def Shortcuts (self) :
        print "shortcuts"
        self.shortcutsPanelShownVar.set ( not self.shortcutsPanelShownVar.get() )

    def Tools (self) :
        print "tools"
        self.toolsPanelShownVar.set ( not self.toolsPanelShownVar.get() )

    def Log ( self ) :
        import Idle
        Idle.start_shell()


    def RSeg ( self ) :
        import Segger.rseg_dialog
        reload ( Segger.rseg_dialog )
        Segger.rseg_dialog.show_dialog()

    def ISeg ( self ) :
        import Segger.iseg_dialog
        reload ( Segger.iseg_dialog )
        Segger.iseg_dialog.show_dialog()

    def SSE ( self ) :
        # self.ssePanelShownVar.set ( not self.ssePanelShownVar.get() )
        import Segger.sse_dialog
        reload ( Segger.sse_dialog )
        Segger.sse_dialog.show_sse_dialog()

    def SegLoop ( self ) :
        # self.ssePanelShownVar.set ( not self.ssePanelShownVar.get() )
        import Segger.segloop_dialog
        reload ( Segger.segloop_dialog )
        Segger.segloop_dialog.show_dialog()

    def SegMod ( self ) :
        # self.ssePanelShownVar.set ( not self.ssePanelShownVar.get() )
        import Segger.segmod_dialog
        reload ( Segger.segmod_dialog )
        Segger.segmod_dialog.show_dialog()

    def SegNA ( self ) :
        # self.ssePanelShownVar.set ( not self.ssePanelShownVar.get() )
        import Segger.segna_dialog
        reload ( Segger.segna_dialog )
        Segger.segna_dialog.show_dialog()


    def Ar ( self ) :
        # self.ssePanelShownVar.set ( not self.ssePanelShownVar.get() )
        import Segger.ar_dialog
        reload ( Segger.ar_dialog )
        Segger.ar_dialog.show_dialog()


    def Spr ( self ) :
        # self.ssePanelShownVar.set ( not self.ssePanelShownVar.get() )
        import Segger.spr_dialog
        reload ( Segger.spr_dialog )
        Segger.spr_dialog.show_dialog()


    def Vr ( self ) :
        # self.ssePanelShownVar.set ( not self.ssePanelShownVar.get() )
        import Segger.vr_dialog
        reload ( Segger.vr_dialog )
        Segger.vr_dialog.show_dialog()


    def CamMono ( self ) :
        chimera.viewer.camera.setMode ( "mono" )

    def CamSBS ( self ) :
        chimera.viewer.camera.setMode ( "DTI side-by-side stereo" )


    def ProMod ( self ) :
        # self.ssePanelShownVar.set ( not self.ssePanelShownVar.get() )
        import Segger.promod_dialog
        reload ( Segger.promod_dialog )
        Segger.promod_dialog.show_dialog()


    def ModelZ ( self ) :
        # self.ssePanelShownVar.set ( not self.ssePanelShownVar.get() )
        import Segger.modelz
        reload ( Segger.modelz )
        Segger.modelz.show_dialog()

    def MapQ ( self ) :
        # self.ssePanelShownVar.set ( not self.ssePanelShownVar.get() )
        import Segger.mapq
        reload ( Segger.mapq )
        Segger.mapq.show_dialog()



    def BioMovie ( self ) :
        # self.ssePanelShownVar.set ( not self.ssePanelShownVar.get() )
        import Segger.biomovie2
        reload ( Segger.biomovie2 )
        Segger.biomovie2.show_dialog()


    def SWIM ( self ) :
        # self.ssePanelShownVar.set ( not self.ssePanelShownVar.get() )
        import Segger.SWIM
        reload ( Segger.SWIM )
        Segger.SWIM.show_dialog()

    def MDFF ( self ) :
        # self.ssePanelShownVar.set ( not self.ssePanelShownVar.get() )
        import Segger.mdff_dialog
        reload ( Segger.mdff_dialog )
        Segger.mdff_dialog.show_dialog()

    def Animate ( self ) :
        # self.ssePanelShownVar.set ( not self.ssePanelShownVar.get() )
        import Segger.animate_dialog
        reload ( Segger.animate_dialog )
        Segger.animate_dialog.close_animate_dialog ()
        Segger.animate_dialog.show_dialog()

    def FlexFit ( self ) :
        # self.ssePanelShownVar.set ( not self.ssePanelShownVar.get() )
        import Segger.flexfit_dialog
        reload ( Segger.flexfit_dialog )
        Segger.flexfit_dialog.show_dialog()


    def Tomolog ( self ) :
        # self.ssePanelShownVar.set ( not self.ssePanelShownVar.get() )
        import Segger.tomolog_dialog
        reload ( Segger.tomolog_dialog )
        Segger.tomolog_dialog.show_dialog()

    def GeoSeg ( self ) :
        # self.ssePanelShownVar.set ( not self.ssePanelShownVar.get() )
        import Segger.geoseg_dialog
        reload ( Segger.geoseg_dialog )
        Segger.geoseg_dialog.show_dialog()


    def MapCOM ( self ) :

        dmap = self.SegmentationMap()

        import axes
        pts, weights = axes.map_points ( dmap )
        if len(pts) == 0 :
            print " - no pts at this threshold?"
            return

        COM, U, S, V = axes.prAxes ( pts )
        print "com:", COM

        #chimera.viewer.camera.center = chimera.Point ( COM[0], COM[1], COM[2] )
        #xf = chimera.Xform.translation ( chimera.Vector( -COM[0], -COM[1], -COM[2] ) )
        #dmap.openState.xform = xf

        p = chimera.Point ( COM[0], COM[1], COM[2] )
        chimera.openModels.cofr = dmap.openState.xform.apply ( p )

        moveCam = 1
        if moveCam :
            p0 = numpy.array ( chimera.viewer.camera.center )
            p1 = numpy.array ( chimera.openModels.cofr )
            for i in range (10) :
                f = float(i) / 9.0
                f1, f2 = 2.0*f*f*f-3.0*f*f+1.0, 3*f*f-2*f*f*f
                P = p0 * f1 + p1 * f2
                chimera.viewer.camera.center = (P[0],P[1],P[2])
                print ".",
            print ""




    def OpenSegmentation(self):

        dmap = self.SegmentationMap()
        dir = os.path.dirname(dmap.data.path) if dmap else None
        import segfile
        segfile.show_open_dialog(dir, self.OpenSegFiles)



    def OpenSegFiles(self, paths_and_types, open = True):

        from chimera import tasks, CancelOperation
        task = tasks.Task('Loading segmentation', modal = True)
        smods = []
        try:
            import segfile
            reload (segfile)
            for path, ftype in paths_and_types:
                if ftype == 'Segmentation':
                    try:
                        smod = segfile.read_segmentation(path, open, task)
                    except CancelOperation:
                        break
                elif ftype == 'Old regions file':
                    dmap = self.SegmentationMap()
                    if dmap is None:
                        from chimera.replyobj import error
                        from os.path import basename
                        error('Segmentation map must be open before opening old-style segmentation file\n\n\t%s\n\nbecause file does not contain grid size and spacing.' % basename(path))
                        return
                    import regionsfile
                    smod = regionsfile.ReadRegionsFile ( path, dmap )
                smods.append(smod)

            if len(smods) == 0:
                umsg ( "No segmentation was loaded." )
                return

            for smod in smods:
                smod.open_map()

            # TODO: Can't control whether marker model is opened.
            smod = smods[-1]
            self.SetCurrentSegmentation(smod)
            v = smod.volume_data()
            if v:
                self.SetMapMenu(v)
            else :
                umsg ( "Volume data not found" )
            try:
                self.RegsDispUpdate (task)
            except CancelOperation:
                pass

        finally:
            task.finished()

        for s in smods:
            mname = os.path.basename(getattr(s, 'map_path', 'unknown'))
            umsg('Opened segmentation %s of map %s, grid size (%d,%d,%d)'
                 % ((s.name, mname) + tuple(s.grid_size())))


        return smods


    def SaveSegmentation(self):

        smod = self.CurrentSegmentation()
        if smod:

            map = smod.volume_data()
            if map == None :
                umsg ( "Map not found - please associate a map first" )
                return

            if hasattr(smod, 'path') and smod.path:
                import segfile
                segfile.write_segmentation(smod, smod.path)
                umsg ( "Saved" )

            else:
                self.SaveSegmentationAs()
        else :
            umsg ( "No segmentation selected" )



    def SaveSegmentationAs(self):

        smod = self.CurrentSegmentation()
        if smod:
            import segfile
            segfile.show_save_dialog(smod, self.path_changed_cb)

    def path_changed_cb(self, seg):

        if seg is self.CurrentSegmentation():
            seg.name = os.path.basename(seg.path)
            self.regions_file.set(seg.name)

    def ModelClosed(self, trigger, n, mlist):

        # Clear menus that are showing closed models.
        if self.cur_dmap in mlist:
            self.SetMapMenu(None)
        if self.cur_seg in mlist:
            self.cur_seg = None
            self.regions_file.set('')

    def MapMenu ( self ) :

        self.mb.menu.delete ( 0, 'end' )        # Clear menu
        from VolumeViewer import Volume
        mlist = OML(modelTypes = [Volume])
        for m in mlist :
            self.mb.menu.add_radiobutton ( label=m.name, variable=self.dmap,
                                           command=lambda m=m: self.MapSelected(m) )

    def SetMapMenu (self, dmap):

        mname = dmap.name if dmap else ''
        self.dmap.set(mname)
        self.cur_dmap = dmap
        #print "Set map menu to ", dmap.name

    def MapSelected ( self, dmap ) :

        self.cur_dmap = dmap
        #if dmap:
        #    dmap.display = True

    def SegmentationMap(self):

        return self.cur_dmap

    def FillSegmentationMenu ( self ) :

        menu = self.mbSegmentationMenu
        menu.delete ( 0, 'end' )      # Clear menu

        open_names = [(m.name, m) for m in regions.segmentations()]
        if len(open_names ) > 0 :
            menu.add_radiobutton ( label="Open regions files:" )
            menu.add_separator()
            open_names.sort()
            for name, smod in open_names:
                menu.add_radiobutton (label=name, variable=self.regions_file,
                          command=lambda smod=smod: self.RFileSelected(smod) )

        smm = self.SegmentationMap()
        if smm == None :
            self.SetCurrentSegmentation(None)
            self.SetMapMenu(None)
            return

        path = os.path.dirname ( smm.data.path ) + os.path.sep
        bname = os.path.splitext ( smm.name ) [0]

        files = os.listdir ( path );
        names_in_path = []
        for f in files :
            if f.find ( bname ) < 0 or f.find('.seg') < 0 : continue
            if f.find ( ".txt" ) >= 0 : continue
            if f.find ( ".mrc" ) >= 0 : continue
            if not f in open_names :
                names_in_path.append ( f )

        if len ( names_in_path ) == 0 : return

        if len(open_names ) > 0 :
            menu.add_separator()

        menu.add_radiobutton ( label="In %s:" % path )
        menu.add_separator()

        for f in names_in_path :
            menu.add_radiobutton (
                label=f, variable=self.regions_file, command=self.RFileSelected )



    def RFileSelected ( self, rmod = None ) :

        if rmod is None:
            mm = self.SegmentationMap()
            if mm == None :
                print self.dmap.get(), "not open";
                return
            path = os.path.dirname(mm.data.path) + os.path.sep
            rfile = self.regions_file.get()
            print " - opening seg file " + rfile + " for map " + mm.name
            rmod = self.OpenSegFiles ( [(path + rfile, 'Segmentation')] )[-1]
        else:
            rmod.display = True
            if rmod.adj_graph : rmod.adj_graph.show_model(True)
            self.ReportRegionCount(rmod)
			# hiding all other segmentations can be annoying sometimes
            if 0 :
	            for m in regions.segmentations():
	                if m != rmod:
	                    m.display = False
	                    if m.adj_graph : m.adj_graph.show_model(False)
            rfile = rmod.name

        self.cur_seg = rmod
        self.SetMapMenu(rmod.volume_data())

        umsg ( "Showing %s - %d regions, %d surfaces" %
               (rfile, len(rmod.regions), len(rmod.surfacePieces)) )


    def CurrentSegmentation ( self, warn = True ):

        if warn and self.cur_seg is None:
            umsg ( "No segmentation chosen" )
        return self.cur_seg

    def SetCurrentSegmentation ( self, smod ):

        self.cur_seg = smod
        self.regions_file.set(smod.name if smod else '')
        if smod:
            self.SetMapMenu(smod.volume_data())

    def NewSurfaceResolution( self, event = None ):

        smod = self.CurrentSegmentation()
        if smod is None:
            return

        self.SetSurfaceGranularity(smod)

    def SetSurfaceGranularity ( self, smod ):

        g = self.surfaceGranularity.get()
        try:
            res = float(g)
        except:
            umsg ('Surface granularity "%s" is not a number' % g)
            return

        if res <= 0:
            return

        if smod.regions:
            from chimera import tasks, CancelOperation
            task = tasks.Task('Changing surface resolution', modal = True)
            try:
                smod.change_surface_resolution(res, task)
            except CancelOperation:
                pass
            finally:
                task.finished()
        else:
            smod.change_surface_resolution(res)


    def NewMaxRegions ( self, event = None ):

        from chimera import tasks, CancelOperation
        task = tasks.Task('Redisplaying regions', modal = True)
        try:
            self.RegsDispUpdate (task)
        except CancelOperation:
            pass
        finally:
            task.finished()


    def CloseHiddenSeg ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        for m in regions.segmentations () :
            if not m.display:
                umsg ( "Closed %s" % m.name )
                m.close()



    def SaveRegsToMRC ( self, regs, dmap, path = None ) :

        segs = set([r.segmentation for r in regs])
        for s in segs:
            if tuple(s.grid_size()) != tuple(dmap.data.size):
                from chimera import replyobj
                replyobj.error('Cannot mask map.\n\n'
                               'Map %s grid size (%d,%d,%d) does not match '
                               'segmentation %s grid size (%d,%d,%d).'
                               % ((dmap.name,) + tuple(dmap.data.size) +
                                  (s.name,) + tuple(s.grid_size())))
                return

        if path is None:
            # Show file chooser dialog.
            fprefix = os.path.splitext(dmap.name)[0]
            if len(regs) == 1 :
                fname = fprefix + "_region_%d.mrc" % regs[0].rid
            else :
                fname = fprefix + "_%d_regions.mrc" % len(regs)
            dir = os.path.dirname ( dmap.data.path )
            import OpenSave
            d = OpenSave.SaveModal ( title = "Save Masked Map",
                                     initialdir = dir, initialfile = fname,
                                     filters = [('MRC map', '*.mrc', '.mrc')] )
            paths_and_types = d.run ( self.toplevel_widget )
            if paths_and_types:
                path = paths_and_types[0][0]
            else:
                return

        (li,lj,lk), (hi,hj,hk) = regions.region_bounds(regs)

        bound = 2
        li = li - bound; lj = lj - bound; lk = lk - bound
        hi = hi + bound; hj = hj + bound; hk = hk + bound

        n1 = hi - li + 1
        n2 = hj - lj + 1
        n3 = hk - lk + 1

        print "Bounds - %d %d %d --> %d %d %d --> %d %d %d" % ( li, lj, lk, hi, hj, hk, n1,n2,n3 )

        umsg ( "Saving %d regions to mrc file..." % len(regs) )

        nmat = numpy.zeros ( (n3,n2,n1), numpy.float32 )
        dmat = dmap.full_matrix()

        #regs_name = ""
        for reg in regs :
            p = reg.points()
            i,j,k = p[:,0],p[:,1],p[:,2]
            nmat[k-lk,j-lj,i-li] = dmat[k,j,i]

        O = dmap.data.origin
        print "origin:", O
        nO = ( O[0] + float(li) * dmap.data.step[0],
               O[1] + float(lj) * dmap.data.step[1],
               O[2] + float(lk) * dmap.data.step[2] )

        print "new origin:", nO

        ndata = VolumeData.Array_Grid_Data ( nmat, nO, dmap.data.step, dmap.data.cell_angles )
        try : nv = VolumeViewer.volume.add_data_set ( ndata, None )
        except : nv = VolumeViewer.volume.volume_from_grid_data ( ndata )

        nv.name = os.path.basename ( path )

        nv.openState.xform = dmap.openState.xform

        nv.write_file ( path, "mrc" )

        if [s.volume_data() for s in segs] == [dmap]:
            umsg ( "Wrote %s" % ( nv.name, ) )
        else:
            umsg ( "Masked map %s, wrote %s" % ( dmap.name, nv.name ) )



    def WriteSelRegionsMRCFile ( self ) :

        if len(self.dmap.get()) == 0 : umsg ("Please select a map from which density will be taken" ); return
        dmap = self.SegmentationMap()
        if dmap == None : umsg ( "%s is not open" % self.dmap.get() ); return

        smod = self.CurrentSegmentation()
        if smod == None : return

        regs = smod.selected_regions()
        if len(regs)==0 :
            umsg ( "Please select one ore more regions to save to .mrc file" )
            return

        self.SaveRegsToMRC ( regs, dmap )


    def WriteAllRegionsMRCFile ( self ) :

        if len(self.dmap.get()) == 0 : umsg ("Please select a map from which density will be taken" ); return
        dmap = self.SegmentationMap()
        if dmap == None : umsg ( "%s is not open" % self.dmap.get() ); return

        smod = self.CurrentSegmentation()
        if smod == None : return

        regs = [ sp.region for sp in smod.surfacePieces
                 if hasattr(sp, 'region')]

        self.SaveRegsToMRC ( regs, dmap )




    def WriteEachRegionMRCFile ( self ) :

        if len(self.dmap.get()) == 0 : umsg ("Please select a map from which density will be taken" ); return
        dmap = self.SegmentationMap()
        if dmap == None : umsg ( "%s is not open" % self.dmap.get() ); return

        smod = self.CurrentSegmentation()
        if smod is None : return

        regs = smod.selected_regions()
        if len(regs)==0 :
            regs = smod.regions

        # Choose file path.
        dir = os.path.dirname ( dmap.data.path )
        fprefix = os.path.splitext ( dmap.name ) [0]
        fname = fprefix + "_region_%d.mrc"
        import OpenSave
        d = OpenSave.SaveModal ( title = "Save Masked Map",
                                 initialdir = dir, initialfile = fname,
                                 filters = [('MRC map', '*.mrc', '.mrc')] )
        paths_and_types = d.run ( self.toplevel_widget )
        if paths_and_types:
            path = paths_and_types[0][0]
        else:
            return
        if not '%d' in path:
            umsg ( "Must include '%d' in map file name for region number" )
            return

        print "Saving each of %d regions to .mrc files" % len(regs)

        for reg in regs :
            self.SaveRegsToMRC ( [reg], dmap, path % (reg.rid,) )



    def MaskMapWRegions ( self ) :

        if len(self.dmap.get()) == 0 : umsg ("Please select a map from which density will be taken" ); return
        dmap = self.SegmentationMap()
        if dmap == None : umsg ( "%s is not open" % self.dmap.get() ); return

        smod = self.CurrentSegmentation()
        if smod is None : return

        regs = smod.selected_regions()
        if len(regs)==0 :
            umsg ( "Please select one ore more regions" )
            return

        nv = regions.mask_volume(regs, dmap)
        if nv is None:
            umsg ('Map size %d,%d,%d is incompatible with mask size %d,%d,%d'
                  % (tuple(dmap.data.size) + tuple(smod.grid_size())))
            return

        return nv


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

        _contour.affine_transform_vertices ( points, smod.seg_map.data.ijk_to_xyz_transform )
        _contour.affine_transform_vertices ( points, Matrix.xform_matrix( smod.openState.xform ) )
        _contour.affine_transform_vertices ( points, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )

        sg = VolumeData.zone_masked_grid_data ( dmap.data, points, smod.seg_map.data.step[0] )

        try : gv = VolumeViewer.volume.add_data_set ( sg, None )
        except : gv = VolumeViewer.volume.volume_from_grid_data ( sg )
        gv.openState.xform = dmap.openState.xform
        #chimera.openModels.add ( [gv] )
        gv.name = "Masked"


    def ExtractDensities ( self ) :

        import Segger.extract_region_dialog
        reload ( Segger.extract_region_dialog )

        Segger.extract_region_dialog.show_extract_region_dialog()


    def Frankensteinify ( self ) :

        import Segger.frankensteinify
        reload ( Segger.frankensteinify )

        Segger.frankensteinify.show_dialog()


    def FitDialog ( self ) :

        import Segger.fit_dialog
        Segger.fit_dialog.close_fit_segments_dialog();
        reload(Segger.fit_dialog);
        Segger.fit_dialog.new_fit_segments_dialog()


    def GeoSegDialog ( self ) :

        import Segger.geoseg;
        reload(Segger.geoseg);
        Segger.geoseg.show_dialog();


    def MaskAnotherMapWRegionsShrink ( self ) :

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

        _contour.affine_transform_vertices ( points, smod.seg_map.data.ijk_to_xyz_transform )
        _contour.affine_transform_vertices ( points, Matrix.xform_matrix( smod.openState.xform ) )
        _contour.affine_transform_vertices ( points, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )

        sg = VolumeData.zone_masked_grid_data ( dmap.data, points, smod.seg_map.data.step[0] )
        regsm = sg.matrix()

        nze = numpy.nonzero ( regsm )

        print nze

        li = numpy.min ( nze[0] )
        lj = numpy.min ( nze[1] )
        lk = numpy.min ( nze[2] )

        hi = numpy.max ( nze[0] )
        hj = numpy.max ( nze[1] )
        hk = numpy.max ( nze[2] )

        bound = 2
        li = li - bound; lj = lj - bound; lk = lk - bound
        hi = hi + bound; hj = hj + bound; hk = hk + bound

        n1 = hi - li + 1
        n2 = hj - lj + 1
        n3 = hk - lk + 1

        print "Bounds - %d %d %d --> %d %d %d --> %d %d %d" % ( li, lj, lk, hi, hj, hk, n1,n2,n3 )

        umsg ( "Saving %d regions to mrc file..." % len(regs) )

        nmat = numpy.zeros ( (n1,n2,n3), numpy.float32 )
        #dmat = dmap.full_matrix()

        print "map grid dim: ", numpy.shape ( dmap.full_matrix() )
        print "masked grid dim: ", numpy.shape ( regsm )
        print "new map grid dim: ", numpy.shape ( nmat )


        #regs_name = ""
        for ii in range ( len(nze[0]) ) :
            i,j,k = nze[0][ii], nze[1][ii], nze[2][ii]
            #nmat[k-lk,j-lj,i-li] = regsm[k,j,i]
            nmat[i-li,j-lj,k-lk] = regsm[i,j,k]

        O = dmap.data.origin
        print "origin:", O
        nO = ( O[0] + float(lk) * dmap.data.step[0],
               O[1] + float(lj) * dmap.data.step[1],
               O[2] + float(li) * dmap.data.step[2] )

        print "new origin:", nO

        ndata = VolumeData.Array_Grid_Data ( nmat, nO, dmap.data.step, dmap.data.cell_angles )
        try : nv = VolumeViewer.volume.add_data_set ( ndata, None )
        except : nv = VolumeViewer.volume.volume_from_grid_data ( ndata )

        nv.name = "Masked"

        nv.openState.xform = dmap.openState.xform


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

            _contour.affine_transform_vertices ( points, smod.seg_map.data.ijk_to_xyz_transform )
            _contour.affine_transform_vertices ( points, Matrix.xform_matrix( smod.openState.xform ) )
            _contour.affine_transform_vertices ( points, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )

            sg = VolumeData.zone_masked_grid_data ( dmap.data, points, smod.seg_map.data.step[0] )
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

            import os.path
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

            import os.path
            nv.name = os.path.splitext (dmap.name) [0] + suff
            nv.openState.xform = dmap.openState.xform

            path = os.path.splitext (dmap.data.path) [0] + suff
            nv.write_file ( path, "mrc" )


    def SubtractRegionsFromMap ( self ) :

        if len(self.dmap.get()) == 0 : umsg ("Please select a map from which density will be taken" ); return
        dmap = self.SegmentationMap()
        if dmap == None : umsg ( "%s is not open" % self.dmap.get() ); return

        smod = self.CurrentSegmentation()
        if smod is None : return

        regs = smod.selected_regions()
        if len(regs)==0 :
            umsg ( "Please select one ore more regions" )
            return

        nv = regions.remove_mask_volume(regs, dmap)
        if nv is None:
            umsg ('Map size %d,%d,%d is incompatible with mask size %d,%d,%d'
                  % (tuple(dmap.data.size) + tuple(smod.grid_size())))
            return

        return nv



    def DetectSym ( self ) :

        self.syms = None

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

        from Measure.symmetry import centers_and_points

        centers, xyz, w = centers_and_points(dmap)
        print "Centers: ", centers
        tcenters = numpy.array(centers, numpy.float32)
        Matrix.transform_points ( tcenters, dmap.data.xyz_to_ijk_transform )


        print "TCenters: ", tcenters

        self.syms = syms
        self.scenters = tcenters
        return syms



    def SetSymColors ( self, regions ) :

        if not hasattr (self, 'syms') :
            return

        from random import random as rand

        # set same color for each symm matrix
        if not hasattr ( self, 'sym_colors' ) or self.sym_type != self.symmetryString.get() :
            print " - making symm colors..."
            self.sym_colors = {}
            self.sym_type = self.symmetryString.get()
            for si, smat in enumerate ( self.syms ) :
                self.sym_colors[si] = ( rand(), rand(), rand(), 1 )


        # in case region color were changed by user:
        if regions :
            clr = None
            for reg in regions :
                if reg.surface_piece.vertexColors is not None :
                    reg.surface_piece.vertexColors = None
                clr = reg.surface_piece.color
                reg.set_color ( clr )
            if clr :
                for reg in regions :
                    reg.set_color ( clr )
                self.sym_colors[0] = regions[0].color



    def ShowSelSymm ( self ) :

        dmap = segmentation_map()

        if dmap == None:
            umsg ( "Please select a map..." )
            return

        smod = self.CurrentSegmentation()
        if smod is None :
            umsg ( "Select a segmentation..." )
            return

        regions = smod.selected_regions()
        if len(regions)==0 :
            umsg ( "Select one or more regions..." )
            return

        if not hasattr(self,'syms') or self.syms == None :
            umsg ( "No symmetry? Press Detect first..." )
            return


        syms = self.syms
        centers = self.scenters

        print "Showing %d symmetric copies..." % len(syms)
        print "Centers:", centers

        com = centers[0]
        t_0_com = ( (1.0,0.0,0.0,-com[0]),
                    (0.0,1.0,0.0,-com[1]),
                    (0.0,0.0,1.0,-com[2]) )
        t_to_com = ( (1.0,0.0,0.0,com[0]),
                    (0.0,1.0,0.0,com[1]),
                    (0.0,0.0,1.0,com[2]) )

        ptf = smod.point_transform()
        #print ptf

        surf = _surface.SurfaceModel ()
        rname = ""

        self.SetSymColors (regions)

        #for si, smat in enumerate ( syms [1 : ] ) :
        for si, smat in enumerate ( syms ) :

            clr = self.sym_colors[si]

            rname = ""
            for i, reg in enumerate (regions) :

                tf = Matrix.multiply_matrices( t_to_com, smat, t_0_com )

                points = numpy.array ( reg.points(), numpy.float32 )
                _contour.affine_transform_vertices ( points, tf )
                #print points

                from MultiScale.surface import surface_points
                vertices, triangles, normals = \
                    surface_points ( points,
                                     resolution = smod.surface_resolution,
                                     density_threshold = 0.1,
                                     smoothing_factor = .25,
                                     smoothing_iterations = 5 )

                _contour.affine_transform_vertices ( vertices, ptf )

                nsp = surf.addPiece ( vertices, triangles, clr )
                nsp.oslName = "Reg_%d_sym_%d" % (reg.rid, si)

                rname = rname + "_%d" % reg.rid



        nn = os.path.splitext(dmap.name)[0]

        surf.name = nn + "_sym_%s" % self.symmetryString.get()  + rname
        chimera.openModels.add ( [surf] )


    def ColorSymmSame ( self ) :

        print "Color symm:"

        smod = self.CurrentSegmentation()
        if smod is None :
            umsg ( "Please segment first..." )
            return

        if not hasattr(self,'syms') or self.syms == None :
            umsg ( "No symmetry? Press Detect first..." )
            return

        regions = smod.selected_regions()
        if len(regions)==0 :
            print " - using all regions"
            regions = None
        else :
            print " - using %d selected regions" % len(regions)

        syms = self.syms
        centers = self.scenters

        if 1 :
            from chimera import tasks, CancelOperation
            task = tasks.Task('Coloring symmetries', modal = True)
            try:
                smod.find_sym_regions2 ( [centers, syms], symColors=None, regs=regions, task=None )
            except CancelOperation:
                umsg('Cancelled coloring symmetries')
            finally:
                task.finished()

        else :
            #smod.calculate_watershed_regions ( mm, thrD, csyms, task )
            smod.find_sym_regions2 ( [centers, syms], symColors=None, regs=regions, task=None )



    def ColorSymmDiff ( self ) :

        print "Color symm:"

        smod = self.CurrentSegmentation()
        if smod is None :
            umsg ( "Please segment first..." )
            return

        if not hasattr(self,'syms') or self.syms == None :
            umsg ( "No symmetry? Press Detect first..." )
            return

        regions = smod.selected_regions()
        if len(regions)==0 :
            print " - using all regions"
            regions = None
        else :
            print " - using %d selected regions" % len(regions)

        syms = self.syms
        centers = self.scenters

        self.SetSymColors (regions)


        if 1 :
            from chimera import tasks, CancelOperation
            task = tasks.Task('Coloring symmetries', modal = True)
            try:
                #smod = self.SegmentAndGroup(show, group, task)
                smod.find_sym_regions2 ( [centers, syms], symColors=self.sym_colors, regs=regions, task=task )
            except CancelOperation:
                umsg('Cancelled coloring symmetries')
            finally:
                task.finished()

        else :
            #smod.calculate_watershed_regions ( mm, thrD, csyms, task )
            smod.find_sym_regions2 ( [centers, syms], symColors=self.sym_colors, regs=regions, task=None )


    def RegsDispUpdate ( self, task = None ) :

        smod = self.CurrentSegmentation()
        if smod is None :
            print " - regs disp update - no smod"
            return

        if smod.volume_data() is None:
            print " - regs disp update - no smod"
            smod.set_volume_data(self.SegmentationMap())

        maxnr = self.MaximumRegionsToDisplay()

        if maxnr > 0 and len(smod.regions) >= maxnr  :
            umsg('Only showing %d of %d regions.' % (maxnr, len(smod.regions)))

        smod.display_regions(self.regsVisMode.get(), maxnr, task)

        if maxnr >= len(smod.regions):
            umsg ( "Showing %d region surfaces" % len(smod.regions) )
        else:
            umsg ( "Showing %d of %d region surfaces" %
                   (maxnr, len(smod.regions)) )

        self.ReportRegionCount(smod)

    def MaximumRegionsToDisplay ( self ) :

        try:
            maxnr = int(self.maxNumRegions.get())
        except:
            maxnr = 0
        return maxnr

    def ReportRegionCount ( self, smod ):

        if smod is None:
            s = ''
        else:
            s = "%s regions" % "{:,}".format( len(smod.regions) )
        self.regionCount["text"] = s


    def RegsDispThr ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        print "%s - thresholding %d regions" % ( smod.name, len(smod.regions) )

        dmap = self.SegmentationMap()
        if dmap == None : print "Map %s not open" % self.dmap.get(); return
        print " - using map:", dmap.name
        dthr = dmap.surface_levels[0]

        for r in smod.regions : r.remove_surface()

        maxnr = self.MaximumRegionsToDisplay()
        for reg in smod.regions :

            if maxnr > 0 and len(smod.surfacePieces) >= maxnr  :
                umsg('Only showing %d of %d regions.' %
                     (len(smod.surfacePieces), len(smod.regions)))
                break
            reg.make_surface()


    def RegsPrint ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        for reg in smod.regions :

            print "%d - %d %d %d" % ( reg.rid, reg.max_point[0], reg.max_point[1], reg.map_point[2] )


    def ReduceMap ( self ) :

        if len(self.dmap.get()) == 0 : umsg ("Please select a map first" ); return
        mm = self.SegmentationMap()
        if mm == None : umsg ( "%s is not open" % self.dmap.get() ); return

        path = os.path.dirname ( mm.data.path ) + os.path.sep
        mname = os.path.splitext ( mm.name )[0]

        d = mm.data
        m2 = d.matrix ( ijk_step=(2,2,2) )
        step2 = ( d.step[0]*2.0, d.step[1]*2.0, d.step[2]*2.0 )
        ld = VolumeData.Array_Grid_Data(m2, d.origin, step2, d.cell_angles, d.rotation,
                           name = mname + '_s2.mrc')

        gv = VolumeViewer.volume.add_data_set ( ld, None )
        gv.name = mname + '_s2.mrc'

        print "writing", path + gv.name
        mod.write_file ( path + gv.name, "mrc" )


    def GetUseSymmetry ( self ) :

        csyms = None
        err_msg = None

        if self.useSymmetry.get () :

            sstring = self.symmetryString.get ()
            if len ( sstring ) == 0 :
                umsg ("Detecting symmetry...")
                self.DetectSym ()
                sstring = self.symmetryString.get ()

            if len ( sstring ) == 0 :
                umsg ( "Enter a symmetry string, e.g. D8" )
                return [None, "No symmetry specified or found"]

            print "Using symmetry:", sstring

            from Measure.symmetry import centers_and_points

            dmap = segmentation_map()
            centers, xyz, w = centers_and_points(dmap)
            #print "Centers: ", centers
            tcenters = numpy.array(centers, numpy.float32)
            Matrix.transform_points ( tcenters, dmap.data.xyz_to_ijk_transform )
            #print "Centers in ijk coords: ", tcenters

            import Symmetry

            if sstring[0] == "D" :
                print "Dihedral syms"
                syms = Symmetry.dihedral_symmetry_matrices ( int(sstring[1]) )
                csyms = [tcenters, syms]

            elif sstring[0] == "C" :
                print "Cyclic syms"
                syms = Symmetry.cyclic_symmetry_matrices ( int(sstring[1]) )
                csyms = [tcenters, syms]

            else :
                err_msg = "Symmetry string not recognized"


        return [ csyms, err_msg ]



    def Segment ( self, show = True, group = True ) :

        smod = self.CurrentSegmentation(warn = False)
        if smod :
            chimera.openModels.close ( [smod] )

        if self.cur_dmap :
            mname, mext = os.path.splitext ( self.cur_dmap.name )
            print " - current map: %s" % self.cur_dmap.name
            remm = []
            for m in chimera.openModels.list() :
                if ".seg" in m.name and mname in m.name :
                    print " - closing %s" % m.name
                    remm.append ( m )
            if len(remm) > 0 :
                chimera.openModels.close ( remm )


        from chimera import tasks, CancelOperation
        task = tasks.Task('Segmenting %s' % self.dmap.get(), modal = True)

        try:
            smod = self.SegmentAndGroup(show, group, task)
        except CancelOperation:
            umsg('Cancelled segmentation')
            return None
        finally:
            task.finished()

        return smod



    def SegmentAndGroup ( self, show = True, group = True, task = None ) :

        if len(self.dmap.get()) == 0 :
            umsg ("Select a density map in the Segment map field" );
            return

        mm = self.SegmentationMap()
        if mm == None : umsg ( "%s is not open" % self.dmap.get() ); return

        thrD = mm.surface_levels[0]
        mm.segmentThreshold = thrD
        print "\n___________________________"
        umsg ( "Segmenting %s, density threshold %f" % (mm.name, thrD) )


        csyms = None

        smod = self.CurrentSegmentation(warn = False)
        if smod is None or smod.volume_data() != mm:
            mbase, msuf = os.path.splitext ( mm.name )
            msp = msuf.find(' ')
            mend = '' if msp == -1 else msuf[msp:]
            segname = mbase + mend + '.seg'
            smod = regions.Segmentation(segname, mm)
            self.SetSurfaceGranularity(smod)
            self.SetCurrentSegmentation(smod)

        if timing: t0 = clock()
        smod.calculate_watershed_regions ( mm, thrD, csyms, task )

        if timing: t1 = clock()
        self.RemoveSmallRegions(smod, task)
        self.RemoveContactRegions(smod, task)
        nwr = len(smod.regions)

        if timing: t2 = clock()

        if group:
            if self.groupMode.get() == 'smooth' :
                self.SmoothAndGroup ( smod, task )
            else :
                self.GroupByCons ( smod, task )


        # Undisplay other segmentations
        if timing: t3 = clock()
        for m in regions.segmentations() :
            if m != smod:
                m.display = False

        self.RegsDispUpdate ( task )	 # Display region surfaces
#        mm.display = False              # Undisplay map

        if timing :
            t4 = clock()
            print "Time %.2f sec: watershed %.2f sec, small %.2f, group %.2f sec, display %.2f sec" % (t4-t0, t1-t0, t2-t1, t3-t2, t4-t3)

        umsg ( '%d watershed regions, grouped to %d regions' % ( nwr, len(smod.regions)) )

        return smod



    def RemoveSmallRegions(self, smod = None, task = None):

        if smod is None:
            smod = self.CurrentSegmentation()
            if smod is None:
                return

        mrs = self.minRegionSize.get()
        try:
            minsize = int ( mrs )
        except:
            print 'Minimum region size "%s" is not an integer' % mrs
            minsize = 1

        if minsize <= 1:
            return

        if task is None:
            from chimera import tasks, CancelOperation
            task = tasks.Task('Removing small regions', modal = True)
            try:
                smod.remove_small_regions(minsize, task)
                self.RegsDispUpdate(task)
            except CancelOperation:
                umsg('Cancelled removing small regions')
                return
            finally:
                task.finished()
        else:
            smod.remove_small_regions(minsize, task)
            self.RegsDispUpdate(task)

        self.ReportRegionCount(smod)


    def RemoveContactRegions(self, smod = None, task = None):

        if smod is None:
            smod = self.CurrentSegmentation()
            if smod is None:
                return

        mrs = self.minContactSize.get()
        try:
            minsize = int ( mrs )
        except:
            print 'Minimum contact size "%s" is not an integer' % mrs
            minsize = 1

        if minsize <= 0:
            return

        if task is None:
            from chimera import tasks, CancelOperation
            task = tasks.Task('Removing contact regions', modal = True)
            try:
                smod.remove_contact_regions(minsize, task)
                self.RegsDispUpdate(task)
            except CancelOperation:
                umsg('Cancelled removing small regions')
                return
            finally:
                task.finished()
        else:
            smod.remove_contact_regions(minsize, task)
            self.RegsDispUpdate(task)

        self.ReportRegionCount(smod)



    def GroupConnectedRegs ( self ) :

        if 0 :
            min_contact = int(self.minConnection.get())

            smod = self.CurrentSegmentation()
            if smod == None : return

            regions = smod.selected_regions()
            if len(regions)==0 :
                regions = smod.regions

            smod.group_connected ( regions, min_contact)

            self.RegsDispUpdate()
            self.ReportRegionCount(smod)
            if regions:
                from regions import TopParentRegions
                nsurfs = [r.surface_piece for r in TopParentRegions(regions)]
                chimera.selection.setCurrent(nsurfs)

        else :

            print " - connection grouping step"


    def GroupByContacts ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        from chimera import tasks, CancelOperation
        task = tasks.Task('Contact grouping', modal = True)
        try:
            regions.group_by_contacts(smod, task)
            self.RegsDispUpdate(task)
        except CancelOperation:
            umsg('Cancelled contact grouping')
        finally:
            task.finished()

        self.SetColorLevelRange()
        self.ReportRegionCount(smod)

        self.contactsPanel.set(True)

    def CloseAll ( self ) :

        umsg ( "Closing all segmentations." )

        dmap = self.SegmentationMap()
        if dmap != None :
            dmap.display = True

        for m in regions.segmentations ():
            print 'Closed', m.name
            m.close()

        self.SetCurrentSegmentation ( None )
        self.ReportRegionCount(None)


    def Associate ( self ) :

        seg = self.CurrentSegmentation ()
        if seg :
            print " - seg: ", seg.name

            if self.cur_dmap :
                print " - map: " + self.cur_dmap.name
                seg.set_volume_data ( self.cur_dmap )
                umsg ( "Map %s is now associated with %s" % (self.cur_dmap.name, seg.name) )
            else :
                umsg ( "No map selected" )


    def SmoothAndGroup ( self, smod, task = None ) :

        try :
            numit = int ( self.numSteps.get() )
            sdev = float ( self.stepSize.get() )
        except :
            umsg ( "Enter an integer for # smoothing steps, float for step size" )
            return

        try :
            targNRegs = int(self.targNRegions.get())
        except :
            umsg ( "Enter an integer for target # of regions" );
            return


        csyms = None
        if self.useSymmetry.get() :
            print "Using symmetry..."
            self.DetectSym ()
            csyms = [self.scenters, self.syms]


        if targNRegs <= 0 :
            umsg ( "# of regions" )
            return

        smod.smooth_and_group(numit, sdev, targNRegs, csyms, task)

        self.ReportRegionCount(smod)




    def GroupByCons ( self, smod, task = None ) :

        try :
            numit = int ( self.numStepsCon.get() )
            #sdev = float ( self.stepSize.get() )
        except :
            umsg ( "Enter an integer for # steps" )
            return

        try :
            targNRegs = int(self.targNRegionsCon.get())
        except :
            umsg ( "Enter an integer for target # of regions" );
            return


        #csyms, sym_err = self.GetUseSymmetry ()
        #if sym_err :
        #    umsg ( sym_err )
        #    return

        csyms = None
        if self.useSymmetry.get() :
            print "Using symmetry..."
            self.DetectSym ()
            csyms = [self.scenters, self.syms]


        if targNRegs <= 0 :
            umsg ( "Enter an integer > 0 for target # of regions" )
            return

        print " - grouping %d steps, target %d" % (numit, targNRegs)

        #smod.smooth_and_group(numit, sdev, targNRegs, csyms, task)
        smod.group_connected_n ( numit, targNRegs, None, csyms, task )



        self.ReportRegionCount(smod)



    def GroupByConsOneStep ( self, task = None  ) :

        smod = self.CurrentSegmentation()
        if smod is None:
            return

        if smod.volume_data() is None:
            umsg ('Segmentation map not opened')
            return

        if len(smod.regions) <= 1:
            umsg ('%s has %d regions' % (smod.name, len(smod.regions)))
            return


        csyms = None
        if self.useSymmetry.get() :
            print "Using symmetry..."
            self.DetectSym ()
            csyms = [self.scenters, self.syms]


        regions = None
        if  self.groupByConsOnlyVis.get() :
            regions = smod.visible_regions()
            if len(regions) == 0 :
                umsg ("Grouping by connections: no visible regions found or they are from a different model" )
                return

            umsg ("Grouping by connections: applying only to %d regions visible" % len(regions) )


        if 1 :
            from chimera import tasks, CancelOperation
            task = tasks.Task('Group by connections', modal = True)
            try :
                newRegs, removedRegs = smod.group_connected_n ( 1, 1, regions, csyms, task )
                #self.RegsDispUpdate ( task )	 # Display region surfaces
            except CancelOperation :
                umsg('Cancelled group by connections')
                return
            finally:
                task.finished()
        else :
            newRegs, removedRegs = smod.group_connected_n ( 1, 1, regions, csyms, task )

        for r in newRegs : r.make_surface (None, None, smod.regions_scale)
        print " - removig %d surfs" % len(removedRegs)
        for r in removedRegs : r.remove_surface()

        self.ReportRegionCount(smod)

        if smod.adj_graph :
            graph.create_graph ( smod, smod.graph_links )

        umsg ( "Got %d regions after grouping by connections" % (len(smod.regions)) )


    def SmoothAndGroupOneStep ( self ) :

        smod = self.CurrentSegmentation()
        if smod is None:
            return

        if smod.volume_data() is None:
            umsg ('Segmentation map not opened')
            return

        if len(smod.regions) <= 1:
            umsg ('%s has %d regions' % (smod.name, len(smod.regions)))
            return

        try :
            step = float ( self.stepSize.get() )
        except :
            umsg ( "Enter <float> for step size" )
            return

        sdev = step + smod.smoothing_level

        csyms = None
        if self.useSymmetry.get() :
            print "Using symmetry..."
            self.DetectSym ()
            csyms = [self.scenters, self.syms]


        umsg ( "Smoothing and grouping, standard deviation %.3g voxels" % sdev)


        from chimera import tasks, CancelOperation
        task = tasks.Task('Smooth and group', modal = True)
        try:
            for i in range ( 10 ) :
                new_regs = len(smod.smooth_and_group(1, sdev, 1, csyms, task))

                # if symmetry is being used we should stop after one step
                # since symmetry can block regions from joining indefinitely
                if new_regs > 0 : break

                umsg ('No new groups smoothing %.3g voxels' % sdev)
                sdev += step
            self.RegsDispUpdate ( task )	 # Display region surfaces

        except CancelOperation:
            umsg('Cancelled smooth and group')
            return
        finally:
            task.finished()

        self.ReportRegionCount(smod)

        if smod.adj_graph :
            graph.create_graph ( smod, smod.graph_links )

        umsg ( "Got %d regions after smoothing %.3g voxels." %
               (len(smod.regions), sdev) )


    def Overlapping ( self ) :

        dmap = self.SegmentationMap()
        if dmap == None :
            umsg ( "No map selected" )
            return

        smod = self.CurrentSegmentation()
        if smod == None :
            umsg ( "No segmentation selected" )
            return

        if len(smod.regions) == 0 :
            umsg ( "No regions found in %s" % smod.name )
            return

        selatoms = chimera.selection.currentAtoms()
        spoints = None

        if len ( selatoms ) > 0 :
            spoints = _multiscale.get_atom_coordinates ( selatoms, transformed = True )

        else :
            mods = chimera.selection._currentSelection.models()
            if len(mods) == 1 :
                mod = mods[0]
                print "Using for selection:", mod.name

                import axes
                spoints, weights = axes.map_points ( mod, True )
                print " - map - got %d points in contour" % len (spoints)
                from _contour import affine_transform_vertices as transform_vertices
                transform_vertices( spoints,  Matrix.xform_matrix( mod.openState.xform ) )
            else :
                umsg ("0 or more than 1 model selected")
                return


        simap = self.PointIndexesInMap ( spoints, dmap )

        umsg ( "Overlapping %d atoms with %d regions" % (
            len(selatoms), len(smod.regions) ) )

        ovRatio = float ( self.overlappingPercentage.get() ) / 100.0
        print " - overlap ratio: %f" % ovRatio

        oregs = []
        for ri, r in enumerate ( smod.regions ) :
            ipoints = r.points()
            noverlap = 0
            for i,j,k in ipoints :
                try : simap[i][j][k]
                except: continue
                noverlap += 1
            ov = float ( noverlap ) / float ( len(ipoints) )
            if ov > ovRatio : oregs.append ( r )
            #if noverlap > 0 : oregs.append ( r )
        regions.select_regions ( oregs )

        umsg ( "Selected %d regions" % ( len(oregs) ) )




    def GroupUsingFits ( self ) :

        dmap = self.SegmentationMap()
        if dmap == None : print "Map %s not open" % self.dmap.get(); return

        smod = self.CurrentSegmentation()
        if smod == None : return

        if len(smod.regions) == 0 : print "No regions in", smod.name; return
        try : dmap.fitted_mols
        except : dmap.fitted_mols = []
        if len(dmap.fitted_mols) == 0 : print "No fits found for", dmap.name; return

        print "Grouping %d regions by overlap to %d fitted structures" % (
            len(smod.regions), len(dmap.fitted_mols) )

        dmap.chain_maps = []

        for mol in dmap.fitted_mols :
            try : mol.fmap.imap
            except : mol.fmap.imap = self.MapIndexesInMap ( dmap, mol.fmap )
            from random import random as rand
            mol.fmap.surf_color = ( rand(), rand(), rand(), 1 )
            dmap.chain_maps.append ( mol.fmap )

#        self.SegAccuracy ( "_fits_acc", True )



    def RegSurfsShowNone ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        for reg in smod.regions :
            if reg.surface_piece:
                reg.surface_piece.display = False


    def RegSurfsShowAll ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        from chimera import tasks, CancelOperation
        task = tasks.Task('Showing all regions', modal = True)
        try:
            self.RegsDispUpdate(task)
        except CancelOperation:
            pass
        finally:
            task.finished()


    def RegSurfsShowOnlySelected ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        regions.show_only_regions(smod.selected_regions())


    def RegSurfsHide ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        sregs = smod.selected_regions()
        #if len(sregs) == 0 : sregs = smod.all_regions()

        for r in sregs : r.hide_surface()


    def RegSurfsShow ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        sregs = smod.selected_regions()
        #if len(sregs) == 0 : sregs = smod.all_regions()

        for r in sregs : r.show_surface()



    def RegSurfsShowAdjacent ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        sregs = smod.selected_regions()
        if len(sregs) == 0 :
            return

        cr = set()
        for r in sregs :
            cr.update(r.contacting_regions())

        umsg ( "Region has %d adjacent regions" % len(cr) )
        for r in cr :
            r.show_surface()






    def RegSurfsShowNotGrouped ( self ) :

        print "Showing not-grouped regions..."

        smod = self.CurrentSegmentation()
        if smod == None : return

        for reg in smod.regions :
            if len(reg.cregs) == 0 :
                if reg.surface_piece:
                    reg.surface_piece.display = True
            else :
                if reg.surface_piece:
                    reg.surface_piece.display = False



    def SelectGrouped ( self ) :

        print "Selecting grouped regions..."

        smod = self.CurrentSegmentation()
        if smod == None : return

        surfs = []
        for reg in smod.regions :
            if len(reg.cregs) > 0 :
                if reg.surface_piece:
                    surfs.append ( reg.surface_piece )

        chimera.selection.clearCurrent ()
        chimera.selection.addCurrent ( surfs )




    def SelectNotGrouped ( self ) :

        print "Showing not-grouped regions..."

        smod = self.CurrentSegmentation()
        if smod == None : return

        surfs = []
        for reg in smod.regions :
            if len(reg.cregs) == 0 :
                if reg.surface_piece:
                    surfs.append ( reg.surface_piece )

        chimera.selection.clearCurrent ()
        chimera.selection.addCurrent ( surfs )



    def RegSurfsShowGrouped ( self ) :

        print "Showing grouped regions..."

        smod = self.CurrentSegmentation()
        if smod == None : return

        sregs = smod.grouped_regions()
        if len(sregs) == 0 :
            umsg ( "No grouped regions" )
            return

        umsg ( "Showing %d grouped regions" % len(sregs) )

        regions.show_only_regions(sregs)



    def RegSurfsTransparent ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        sregs = smod.selected_regions()
        if len(sregs) == 0 : sregs = smod.all_regions()

        for r in sregs :
            if r.has_surface():
                cr,cg,cb = r.surface_piece.color[:3] #r.color[:3]
                r.surface_piece.color = ( cr, cg, cb, REG_OPACITY )
                r.surface_piece.displayStyle = r.surface_piece.Solid


    def RegSurfsOpaque ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        sregs = smod.selected_regions()
        if len(sregs) == 0 : sregs = smod.all_regions()

        for r in sregs :
            if r.has_surface():
                cr,cg,cb = r.surface_piece.color[:3] #r.color[:3]
                r.surface_piece.color = ( cr, cg, cb, 1.0 )
                r.surface_piece.displayStyle = r.surface_piece.Solid


    def RegSurfsMesh ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        sregs = smod.selected_regions()
        if len(sregs) == 0 : sregs = smod.all_regions()

        for r in sregs :
            if r.has_surface():
                cr,cg,cb = r.surface_piece.color[:3] #r.color[:3]
                r.surface_piece.color = ( cr, cg, cb, 1.0 )
                r.surface_piece.displayStyle = r.surface_piece.Mesh
                r.surface_piece.lineThickness = 1.0



    def SelectAllRegions ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        sel_regs = set ( smod.selected_regions() )
        surfs = [r.surface_piece for r in smod.regions
                 if 1]

        chimera.selection.clearCurrent ()
        chimera.selection.addCurrent ( surfs )


    def Invert ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        sel_regs = set ( smod.selected_regions() )
        surfs = [r.surface_piece for r in smod.regions
                 if not r in sel_regs and r.surface_piece]

        chimera.selection.clearCurrent ()
        chimera.selection.addCurrent ( surfs )


    def Group ( self ):


        if NothingSelected():

            if self.groupMode.get() == 'smooth' :
                self.SmoothAndGroupOneStep()
            else :
                self.GroupByConsOneStep()

        else:
            self.JoinSelRegs()


    def JoinSelRegs ( self ) :

        smod = self.CurrentSegmentation()
        if smod is None : return

        regs = smod.selected_regions()
        if len(regs)==0 :
            umsg ( "No regions selected" )
            return
        regs = regions.TopParentRegions(regs)

        jreg = smod.join_regions ( regs )
        jreg.make_surface(None, None, smod.regions_scale)

        if smod.adj_graph :
            graph.create_graph ( smod, smod.graph_links )

        chimera.selection.setCurrent([jreg.surface_piece])

        self.ReportRegionCount(smod)
        umsg ( "Grouped %d regions" % len(regs) )




    def DelSelRegs ( self ) :

        smod = self.CurrentSegmentation()
        if smod is None :
            umsg ( "No segmentation selected..." )
            return

        regs = smod.selected_regions()
        if len(regs)==0 :
            umsg ( "Select one or more regions to delete" )
            return

        smod.remove_regions ( regs, update_surfaces = True, remove_children = True )

        self.ReportRegionCount(smod)
        umsg ( "Deleted %d regions" % len(regs) )


    def DelExcSelRegs ( self ) :

        smod = self.CurrentSegmentation()
        if smod is None :
            umsg ( "No segmentation selected..." )
            return

        sel_regs = smod.selected_regions()
        if len(sel_regs)==0 :
            umsg ( "No regions selected..." )
            return

        dregs = [r for r in smod.regions
                 if not r in sel_regs]

        smod.remove_regions ( dregs, update_surfaces = True, remove_children = True )

        self.ReportRegionCount(smod)
        umsg ( "Deleted %d regions" % len(dregs) )


    def Ungroup ( self ):

        if NothingSelected():
            self.UngroupLastSmoothing()
        else:
            self.UngroupSelRegs()




    def SafeCreateSurfsForRegs ( self, smod, rlist, rregs ) :

        maxnr = self.MaximumRegionsToDisplay()

        nsurfs = 0
        for r in smod.regions :
            if r.has_surface() :
                nsurfs += 1

        print " - %d surfs have pieces before" % nsurfs

        # surfs that will go away...
        for r in rregs :
            if r.has_surface() :
                nsurfs -= 1

        print " - %d surfs will have pieces after removing selected" % nsurfs


        if nsurfs >= maxnr :
            umsg('Ungrouped to %d regions, but did not show their surfaces, see Options' % len(rlist) )

        else :
            canshow = maxnr - nsurfs

            if canshow < len(rlist) :
                umsg('Ungrouped to %d regions, but did not show all surfaces, see Options' % len(rlist) )
            else :
                umsg('Ungrouped to %d regions' % len(rlist) )

            from chimera import tasks, CancelOperation
            task = tasks.Task('Adding surfaces', modal = True)
            try:
                for ri, reg in enumerate ( rlist ) :
                    if ri >= canshow :
                        break
                    reg.make_surface(None, None, smod.regions_scale)
            except CancelOperation:
                pass
            finally:
                task.finished()



    def ShowNumSubRegs ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        sregs = smod.selected_regions()
        if len(sregs) == 0 :
            umsg ( "No regions selected" )
            return
        sregs = regions.TopParentRegions(sregs)

        num = 0
        for r in sregs :
            if len(r.cregs) == 0 :
                pass
            else :
                num += len(r.cregs)

        umsg ( "selected regions have %d total sub regions" % num )



    def UngroupSelRegs ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        sregs = smod.selected_regions()
        if len(sregs) == 0 :
            umsg ( "No regions selected" )
            return
        sregs = regions.TopParentRegions(sregs)

        chimera.selection.clearCurrent ()

        [rlist, removedRegs] = smod.ungroup_regions ( sregs )
        self.SafeCreateSurfsForRegs ( smod, rlist, removedRegs )
        for r in removedRegs : r.remove_surface()

        print " - now %d regions" % len(smod.regions)

        if smod.adj_graph :
            graph.create_graph ( smod, smod.graph_links )

        chimera.selection.setCurrent ( [r.surface_piece for r in rlist if (hasattr(r,'surface_piece') and r.surface_piece != None)] )

        self.ReportRegionCount(smod)



    def UngroupAllRegs ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        rlist = list(smod.regions)
        [rlist2, removedRegs] = smod.ungroup_regions(rlist)

        self.SafeCreateSurfsForRegs ( smod, rlist2, removedRegs )
        for r in removedRegs : r.remove_surface()

        self.ReportRegionCount(smod)



    def UngroupLastSmoothing ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        levels = [r.smoothing_level for r in smod.regions]
        if len(levels) == 0:
            return
        slev = max(levels)

        rlev = [r for r in smod.regions if r.smoothing_level == slev]
        rlist2 = []
        removedRegs = []

        from chimera import tasks, CancelOperation
        task = tasks.Task('Ungrouping', modal = True)
        try:
            [rlist2, removedRegs] = smod.ungroup_regions(rlev, task)
        except CancelOperation:
            pass
        finally:
            task.finished()

        self.SafeCreateSurfsForRegs ( smod, rlist2, removedRegs )
        for r in removedRegs : r.remove_surface()

        levels = [r.smoothing_level for r in smod.regions]
        smod.smoothing_level = max(levels)

        if smod.adj_graph :
            graph.create_graph ( smod, smod.graph_links )

        #umsg ( "Ungrouped to %.3g voxel smoothing, %d regions" % (smod.smoothing_level, len(smod.regions)) )
        self.ReportRegionCount(smod)


    def CloseRegions ( self ) :

        smod = self.CurrentSegmentation()
        if smod is None : return

        chimera.openModels.remove ( smod )
        self.SetCurrentSegmentation(None)
        self.ReportRegionCount(None)

        if smod.adj_graph : smod.adj_graph.close()


    def CloseSeg ( self ) :

        smod = self.CurrentSegmentation()
        if smod is None : return

        smod.close()
        self.SetCurrentSegmentation(None)



    def RegionsVolume ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        sregs = smod.selected_regions()
        print "%d selected regions" % len(sregs)

        if len(sregs) == 0 :
            sregs = smod.regions

        if len(sregs) == 0 :
            umsg ( "No regions found in %s" % smod.name )
            return

        tvol = sum([reg.enclosed_volume() for reg in sregs])
        pcount = sum([reg.point_count() for reg in sregs])

        rw = "region"
        if len(sregs) > 1 : rw = "regions"
        umsg ( "Volume of %d %s: %.3g Angstroms^3, %d points" % ( len(sregs), rw, tvol, pcount ) )

    def RegionMeanAndSD ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        sregs = smod.selected_regions()
        if len(sregs) == 0 :
            umsg ( "No regions selected in %s" % smod.name )
            return

        v = self.SegmentationMap()
        if v is None:
            v = smod.volume_data()
            if v is None:
                umsg ( 'No map specified' )
                return

        means, sdevs = regions.mean_and_sd(sregs, v)
        for r, m, sd in zip(sregs, means, sdevs):
            umsg ( 'Region %d mean %.5g, SD %.5g' % (r.rid, m, sd) )


    def Graph ( self ) :

        smod = self.CurrentSegmentation()
        if smod:
            graph.create_graph(smod,"uniform")

    def GraphAvgD ( self ) :

        smod = self.CurrentSegmentation()
        if smod:
            graph.create_graph(smod,"avgd")

    def GraphMaxD ( self ) :

        smod = self.CurrentSegmentation()
        if smod:
            graph.create_graph(smod,"maxd")

    def GraphN ( self ) :

        smod = self.CurrentSegmentation()
        if smod:
            graph.create_graph(smod,"N")


    def LoadGraph ( self ) :

        smod = self.CurrentSegmentation()
        if smod is None:
            graph.open_skeleton(smod)

    def SaveGraph ( self ) :

        smod = self.CurrentSegmentation()
        if smod:
            graph.read_graph(smod)

    def CloseGraph ( self ) :

        smod = self.CurrentSegmentation()
        if smod:
            graph.close(smod)
            smod.display = True

    def GroupBySkeleton ( self ) :

        smod = self.CurrentSegmentation()
        if smod:
            skeleton.group_by_skeleton(smod)
            smod.display = True
            smod.display_regions()
            self.ReportRegionCount(smod)

    def RemoveGraphLinks ( self ) :

        graph.remove_graph_links()


    def ShowRegionsAxes ( self, regs ) :

        smod = self.CurrentSegmentation()
        if smod is None: return

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
            sp.Extents = numpy.asarray ( numpy.max ( numpy.abs ( ppoints ), 0 ) )[0]

            sp.Extents[0] += 5.0
            sp.Extents[1] += 5.0
            sp.Extents[2] += 5.0

            import axes
            reload (axes)

            if 0 :
            # for ribosome direction
                sp.Extents[1] = sp.Extents[1] * float(self.axesFactor.get())

                sp.axes = axes.AxesMod ( sp.COM, sp.U, sp.Extents, 6, 1.0, alignTo = sp.model )
            else :
                sp.axes = axes.AxesMod ( sp.COM, sp.U, sp.Extents, 1.0, 1.1, alignTo = sp.model )

            sp.axes.name = "region_%d_axes" % r.rid



    def ShowRegionAxesSelected ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        sregs = smod.selected_regions()
        if len(sregs)==0 : print "no selected regions found"; return

        self.ShowRegionsAxes ( sregs )


    def HideRegionAxes ( self ) :

        print "hiding axes"
        for m in OML() :
            t = m.name.split ("_")
            if t[0] == "region" and t[2] == "axes" :
                print "- removing", m.name
                chimera.openModels.close( m )




    def PointIndexesInMap ( self, points, ref_map ) :

        print "Making map indices for %d points in %s" % ( len(points), ref_map.name )

        _contour.affine_transform_vertices ( points, Matrix.xform_matrix ( ref_map.openState.xform.inverse() ) )
        _contour.affine_transform_vertices ( points, ref_map.data.xyz_to_ijk_transform )

        imap = {}
        for fi, fj, fk in points :
            if 0 :
                i, j, k = int(numpy.round(fi)), int(numpy.round(fj)), int(numpy.round(fk))
                try : mi = imap[i]
                except : mi = {}; imap[i] = mi
                try : mij = mi[j]
                except : mij = {}; mi[j] = mij
                mij[k] = 1
                continue
            for i in [ int(numpy.floor(fi)), int(numpy.ceil(fi)) ] :
                for j in [ int(numpy.floor(fj)), int(numpy.ceil(fj)) ] :
                    for k in [ int(numpy.floor(fk)), int(numpy.ceil(fk)) ] :
                        try : mi = imap[i]
                        except : mi = {}; imap[i] = mi
                        try : mij = mi[j]
                        except : mij = {}; mi[j] = mij
                        mij[k] = 1


        return imap #, C, bRad


    def MapIndexesInMap ( self, ref_map, mask_map ) :

        thr = mask_map.surface_levels[0]
        mm = mask_map.data.matrix()
        mm = numpy.where ( mm > thr, mm, numpy.zeros_like(mm) )

        nze = numpy.nonzero ( mm )
        nzs = numpy.array ( [nze[2], nze[1], nze[0]] )

        # the copy is needed! otherwise the _contour.afine_transform does not work for some reason
        points =  numpy.transpose ( nzs ).astype(numpy.float32)

        #points = numpy.zeros ( ( len(nze[0]), 3 ), numpy.float32 )
        #for ei, i in enumerate ( nze[0] ) :
        #    j = nze[1][ei]
        #    k = nze[2][ei]
        #    points[ei][0], points[ei][1], points[ei][2] = float (k), float(j), float(i)

        #print points[0]

        print "Making map indices for %s in %s" % ( mask_map.name, ref_map.name )
        print " - %d points above %.3f" % ( len(points), thr )

        # transform to index reference frame of ref_map
        f1 = mask_map.data.ijk_to_xyz_transform
        f2 = Matrix.xform_matrix ( mask_map.openState.xform )
        f3 = Matrix.xform_matrix ( ref_map.openState.xform.inverse() )
        f4 = ref_map.data.xyz_to_ijk_transform

        tf = Matrix.multiply_matrices( f2, f1 )
        tf = Matrix.multiply_matrices( f3, tf )
        tf = Matrix.multiply_matrices( f4, tf )
        _contour.affine_transform_vertices ( points, tf )

        #_contour.affine_transform_vertices ( points, f1 )
        #_contour.affine_transform_vertices ( points, f2 )
        #_contour.affine_transform_vertices ( points, f3 )

        #print points[0]

        #com = numpy.sum (points, axis=0) / len(points)
        #C = chimera.Vector ( com[0], com[1], com[2] )
        #comv = numpy.ones_like ( points ) * com
        #points_v = points - comv
        #bRad = numpy.sqrt ( numpy.max ( numpy.sum ( numpy.square (points_v), 1 ) ) )

        # transform points to indexes in reference map
        # _contour.affine_transform_vertices ( points, ref_map.data.xyz_to_ijk_transform )

        imap = {}
        for fk, fj, fi in points :

            for i in [ int(numpy.floor(fi)), int(numpy.ceil(fi)) ] :
                for j in [ int(numpy.floor(fj)), int(numpy.ceil(fj)) ] :
                    for k in [ int(numpy.floor(fk)), int(numpy.ceil(fk)) ] :

                        try : mi = imap[i]
                        except : mi = {}; imap[i] = mi
                        try : mij = mi[j]
                        except : mij = {}; mi[j] = mij
                        mij[k] = 1


        return imap #, C, bRad


    def ShowGroupSurfaces ( self ) :

        smod = self.CurrentSegmentation()
        if smod is None : return

        sregs = smod.selected_regions()
        regs = set([r.top_parent() for r in sregs])
        if len(regs)==0 :
            regs = smod.regions

        for r in regs:
            for c in r.all_children():
                c.remove_surface()
            r.make_surface()

        if sregs:
            surfs = [r.surface() for r in regs if r.has_surface()]
            from chimera import selection
            selection.setCurrent(surfs)



    def ShowUngroupedSurfaces ( self ) :

        smod = self.CurrentSegmentation()
        if smod is None : return

        sregs = smod.selected_regions()
        regs = set([r.top_parent() for r in sregs])
        if len(regs)==0 :
            regs = smod.regions

        from chimera import tasks, CancelOperation
        task = tasks.Task('Showing ungrouped surfaces', modal = True)
        try:
            for i, r in enumerate(regs):
                if task and i % 100 == 0:
                    task.updateStatus('region %d of %d' % (i, len(regs)))
                if r.has_children():
                    r.remove_surface()
                    for c in r.childless_regions():
                        c.make_surface()
                else:
                    r.make_surface()
        except CancelOperation:
            pass
        finally:
            task.finished()

        if sregs:
            from regions import all_regions
            surfs = [r.surface() for r in all_regions(regs) if r.has_surface()]
            from chimera import selection
            selection.setCurrent(surfs)


    def SelectGroups ( self ) :

        smod = self.CurrentSegmentation()
        if smod is None : return

        rlist = [r for r in smod.regions if r.has_children()]
        from regions import all_regions
        surfs = [r.surface() for r in all_regions(rlist) if r.has_surface()]
        from chimera import selection
        selection.setCurrent(surfs)



    # Regions that have voxels on the mask boundary.
    def SelectBoundaryRegions ( self, pad = 3 ) :

        smod = self.CurrentSegmentation()
        if smod is None : return

        m = smod.mask
        if m is None:
            return

        import _segment
        b = _segment.region_bounds(m)

        rset = set()
        kmax, jmax, imax = [(s-1)-pad for s in m.shape]
        for r in smod.childless_regions():
            i = r.rid
            if (i < len(b) and b[i,6] > 0 and
                (b[i,0] <= pad or b[i,1] <= pad or b[i,2] <= pad or
                 b[i,3] >= kmax or b[i,4] >= jmax or b[i,5] >= imax)):
                rset.add(r.top_parent())
        from regions import all_regions, select_regions
        select_regions(all_regions(rset))


    def SelectNonPlacedRegions ( self ) :

        if len(self.dmap.get()) == 0 : umsg ("Please select a density map"); return
        dmap = self.SegmentationMap()
        if dmap == None : umsg ( "%s is not open" % self.dmap.get() ); return

        smod = self.CurrentSegmentation()
        if smod is None : return

        nsel = 0
        chimera.selection.clearCurrent ()
        for sp in smod.surfacePieces :
            try :
                sp.region.placed
            except :
                chimera.selection.addCurrent ( sp )
                nsel = nsel + 1

        print "%d non-placed regions" % nsel


    def mouse_group_cb(self):

        gmm = self.group_mouse_mode
        if self.mouse_group.get():
            if gmm is None:
                import mousemode
                gmm = mousemode.Group_Connected_Mouse_Mode()
                self.group_mouse_mode = gmm
            button, modifiers = self.mouse_button_spec()
            gmm.bind_mouse_button(button, modifiers)
        elif gmm:
            gmm.unbind_mouse_button()

    def mouse_group_button_cb(self):

        if self.mouse_group.get() and self.group_mouse_mode:
            button, modifiers = self.mouse_button_spec()
            self.group_mouse_mode.bind_mouse_button(button, modifiers)

    def mouse_button_spec(self):

        name = self.mouse_group_button.variable.get()
        name_to_bspec = {'button 1':('1', []), 'ctrl button 1':('1', ['Ctrl']),
                         'button 2':('2', []), 'ctrl button 2':('2', ['Ctrl']),
                         'button 3':('3', []), 'ctrl button 3':('3', ['Ctrl'])}
        bspec = name_to_bspec[name]
        return bspec


def NothingSelected():

    from chimera import selection
    return selection.currentEmpty()


def volume_segmentation_dialog ( create=False ) :

  from chimera import dialogs
  return dialogs.find ( Volume_Segmentation_Dialog.name, create=create )


def show_volume_segmentation_dialog ():

    print "hi"

    from chimera import dialogs
    d = volume_segmentation_dialog ( create = True )
    # Avoid transient dialog resizing when created and mapped for first time.
    #d.toplevel_widget.update_idletasks ()
    #d.enter()

    d = dialogs.display(Volume_Segmentation_Dialog.name)

    #from Accelerators import add_accelerator
    #add_accelerator('gg', 'Group regions', d.JoinSelRegs )
    #add_accelerator('uu', 'Ungroup regions', d.UngroupSelRegs )
    #add_accelerator('rh', 'Hide regions', d.RegSurfsHide )
    #add_accelerator('rs', 'Show regions', d.RegSurfsShowOnlySelected )
    #add_accelerator('dd', 'Delete regions', d.DelSelRegs )
    #print "gg - groups regions"
    #print "uu - ungroup regions"

    return d



def current_segmentation(warn = True):

    d = volume_segmentation_dialog()
    if d:
        return d.CurrentSegmentation(warn)
    elif warn:
        umsg ( "No segmentation opened" )
    return None

def segmentation_map():

    d = volume_segmentation_dialog()
    if d:
        return d.SegmentationMap()
    return None




# -----------------------------------------------------------------------------
#
from chimera import dialogs
dialogs.register (Volume_Segmentation_Dialog.name, Volume_Segmentation_Dialog, replace = True)
