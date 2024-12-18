
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
from Segger import showDevTools, timing, seggerVersion
from CGLutil.AdaptiveTree import AdaptiveTree

from chimera.resCode import nucleic3to1
from chimera.resCode import protein3to1

devMenus = True

import qscores
reload (qscores)

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

from chimera.resCode import nucleic3to1
from chimera.resCode import protein3to1, protein1to3
protein3to1['HSD'] = protein3to1['HIS']
protein3to1['HSE'] = protein3to1['HIS']


# protein-dna
xSeq = "MTAGTETDTQPAQLCAADSHDMIRVHGARENNLKNVQVEIPKRRLTVFTGVSGSGKSSLVFDTIAAESQRLINETYSAFI"\
        "QGFMPTLARPEVDVLDGLTTAILVDQQPMGTSLRSTVGTATDAGTLLRILFSRLAKPYIGTQKAFAFNVASADASGVLVV"\
        "NGKKIEKGFSVVGGMCLACEGIGSVSDIDPAQLFDASKSLADGAITVPGWKPDGWVVQSFTESGFFDPHKAIRDYTEQER"\
        "HGFLHGDPVKVKVKGVNTTYEGLLARVRKSFLSKDKETLQPHIRAFVDRAVTFSACSECHGTRLSETARSAKIDGLSIAD"\
        "ASAMQISDLAAWIRGLTDPSVTTLLTVLGQTLESFVQIGLGYLSLDRSSSTLSGGEAQRVKMVRHLGSALTDVTYVFDEP"\
        "TVGLHPHDIQRMNELLLRLRDKGNTVLVVEHKPETIVIADHVVDLGPLAGTKGGEVVFEGTVEGLRASGTVTGRHLDDRA"\
        "SLKPSVRQRTGVVEVRGADAHNLRDVDVDIPLGVLTVVTGVAGSGKSSLIHGSVAGRDGVVTVDQSPIKGSRRSNPATYT"\
        "GMLEPIRKTFAKANGVKPALFSPNSEGACPTCKGAGVIYTDLAIMAGVATTCEDCGGKRFQPSVLQYRVGGRDISEVFAM"\
        "PVAEAAEFFRTGEARTPAACTVLDRLAEVGLGYLSLGQPLTTLSGGERQRLKLAGHMGGAGSVYILDEPTSGLHLADVEQ"\
        "LLRLLDRLVDSGKTVIVVEHHQAVMAHADWIIDLGPGAGHDGGRVVFEGTPADLVAARSTLTGEHLAQYVGA"

# clc-2
xSeq = "MAAAAAEEGMEPRALQYEQTLMYGRYTQDLGAFAKEEAARIRLGGPEPWKGPPSSRAAPELLEYGRSRCARCRVCSVRCH"\
        "KFLVSRVGEDWIFLVLLGLLMALVSWVMDYAIAACLQAQQWMSRGLNTSILLQYLAWVTYPVVLITFSAGFTQILAPQAV"\
        "GSGIPEMKTILRGVVLKEYLTLKTFIAKVIGLTCALGSGMPLGKEGPFVHIASMCAALLSKFLSLFGGIYENESRNTEML"\
        "AAACAVGVGCCFAAPIGGVLFSIEVTSTFFAVRNYWRGFFAATFSAFIFRVLAVWNRDEETITALFKTRFRLDFPFDLQE"\
        "LPAFAVIGIASGFGGALFVYLNRKIVQVMRKQKTINRFLMRKRLLFPALVTLLISTLTFPPGFGQFMAGQLSQKETLVTL"\
        "FDNRTWVRQGLVEELEPPSTSQAWNPPRANVFLTLVIFILMKFWMSALATTIPVPCGAFMPVFVIGAAFGRLVGESMAAW"\
        "FPDGIHTDSSTYRIVPGGYAVVGAAALAGAVTHTVSTAVIVFELTGQIAHILPVMIAVILANAVAQSLQPSLYDSIIRIK"\
        "KLPYLPELGWGRHQQYRVRVEDIMVRDVPHVALSCTFRDLRLALHRTKGRMLALVESPESMILLGSIERSQVVALLGAQL"\
        "SPARRRQHMQERRATQTSPLSDQEGPPSPEASVCFQVNTEDSAFPAARGETHKPLKPALKRGPSVTRNLGESPTGSAESA"\
        "GIALRSLFCGSPPPEAASEKLESCEKRKLKRVRISLASDADLEGEMSPEEILEWEEQQLDEPVNFSDCKIDPAPFQLVER"\
        "TSLHKLRKAIEGSVTAQGVKVRPPLASFRDSATSSSDTETTEVHALWGPHSRHGLPREGSPSDSDDKCQ"

# clc-0
xSeq = "MSHEKNEASGYPEAQSWKSQEAMLGARTEVSRWRAVKNCLYRHLVKVLGEDWIFLLLLGALMALVSWAMDFIGSRGLRFYKY"\
        "LFALVEGNIGLQYLVWVCYPLALILFSSLFCQIVSPQAVGSGIPELKTIIRGAVLHEYLTLRTFVAKTVGLTVALSAGFP"\
        "LGKEGPFVHIASICATLLNQLLCFISGRREEPYYLRADILTVGCALGISCCFGTPLAGVLFSIEVTCSHFGVRSYWRGFL"\
        "GGAFSAFIFRVLSVWVKDTVTLTALFKTNFRGDIPFDLQEMPAFAIIGIASGFFGALFVYLNRQIIVFMRKKNFVTKIL"\
        "KKQRLIYPAVVTFVLATLRFPPGVGQFFGAGLMPRETINSLFDNYTWTKTIDPRGLGNSAQWFIPHLNIFIVMALYFVM"\
        "HFWMAALAVTMPVPCGAFVPVFNLGAVLGRFVGELMALLFPDGLVSNGNLYHILPGEYAVIGAAAMTGAVTHAVSTAVIC"\
        "FELTGQISHVLPMMVAVILANMVAQGLQPSLYDSIIQIKKLPYLPELSWSSANKYNIQVGDIMVRDVTSIASTSTYGDLL"\
        "HVLRQTKLKFFPFVDTPDTNTLLGSIDRTEVEGLLQRRISAYRRQPAAAAEADEEGRNGETGASFTGEAESSFAYIDQED"\
        "AEGQQREGLEAVKVQTEDPRPPSPVPAEEPTQTSGIYQKKQKGTGQVASRFEEMLTLEEIYRWEQREKNVVVNFETCRI"\
        "DQSPFQLVEGTSLQKTHTLFSLLGLDRAYVTSMGKLVGVVALAEIQAAIEGSYQKGFRLPPPLASFRDVKHARNSGRTA"\
        "TSNSSGK"

xSeq = "KPATLRCVPKIRPTTHPTNASHQLTILPTNSSFTSFLTISPKMREKRNLLETRLNVSDTVTLPTAPNMNSEPTLQPQTGEITNRMMDLTLNSSTATPVSPGSVTKETTTVIVTTTKSLPSDQVMLVYDQQEVE"

xSeq = "KPATLRCVPKIRPTTHPTNASHQLTILPTNSSFTSFLTISPKMREKRNLLETRLNVSDTVTLPTAPNMNSEPTLQPQTGEITNRMMDLTLNSSTATPVSPGSVTKETTTVIVTTTKSLPSDQVMLVYDQQEVERSSKPTCPPPELLGGPSVFIFPPKPKDTLMISRTPEVTCVVVDVSQDDPEVQFTWYINNEQVRTARPPLREQQFNSTIRVVSTLPIAHQDWLRGKEFKCKVHNKALPAPIEKTISKARGQPLEPKVYTMGPPREELSSRSVSLTCMINGFYPSDISVEWEKNGKAEDNYKTTPAVLDSDGSYFLYSKLSVPTSEWQRGDVFTCSVMHEALHNHYTQKSISRSPGK"

#xSeq = "MADDKVAILTDDEEEQKRKYVLADPFNGISREPEPPSNETPSSTETSAIPEEEIDWIEKHCVKINNDLLISKVFYFFFYSAYGSLYPLLPVY"\
#        "YKQLGMSPSQSGLLVGIRYFIEFCSAPFWGVVADRFKKGKIVLLFSLLCWVLFNLGIGFVKPATLRCVPKIRPTTHPTNASHQLTILPTNSS"\
#        "FTSFLTISPKMREKRNLLETRLNVSDTVTLPTAPNMNSEPTLQPQTGEITNRMMDLTLNSSTATPVSPGSVTKETTTVIVTTTKSLPSDQVML"\
#        "VYDQQEVEAIFLVILVVVIIGEFFSASSVTIVDTVTLQYLGKHRDRYGLQRMWGSLGWGLAMLSVGIGIDYTHIEVLIDGKGCKPPEYRNYQI"\
#        "VFIVFGVLMTMALIVATQFRFRYNHFKNDDSKGKEVEIPQVERNNSTESSEETPTTTSHSQAFNFWDLIKLLCSVQYGSVLFVAWFMGFGYGFV"\
#        "FTFLYWHLEDLNGTTTLFGVCSVLSHVSELTAYFFSHKLIELIGHIRVLYIGLACNTARYIYISYLENAWTVLPMEVLQGVTHAAIWAACISYLS"\
#        "AAVPPELRTSAQGILQGLHLGLGRGCGAMIGGVLVNYFGAAATFRGIGMACLVILLLFALIQWLAVPDEEEDKTMLAERIPVPSSPVPIATIDLV"\
#        "QQQTEDVMPRIEPRLPPKKTKHQEEQEDVNKPAWGVSSSPWVTFVYALYQIKEMMQLTRDNRASEIQPLQGTNENRENSPAGRAQPVPCETHSDPS"\
#        "RNQPSPDAAASQTQTSPAHPSVDPCTEESEEQQAQLAAGGH"


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

            b = Tkinter.Button(ff, text="_", command=self.Model)
            b.grid (column=5, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="-b", command=self.FindBonds)
            b.grid (column=6, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="-g", command=self.FindBonds2)
            b.grid (column=7, row=0, sticky='w', padx=0, pady=1)

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

            b = Tkinter.Button(ff, text="Select BB", command=self.SelectSelBB)
            b.grid (column=2, row=1, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="All", command=self.SelectAll)
            b.grid (column=3, row=1, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="Show", command=self.ShowSel)
            b.grid (column=4, row=1, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="Hide", command=self.HideSel)
            b.grid (column=5, row=1, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="Only", command=self.ShowSelOnly)
            b.grid (column=6, row=1, sticky='w', padx=0, pady=1)

            b = Tkinter.Button(ff, text="All", command=self.ShowAll)
            b.grid (column=7, row=1, sticky='w', padx=0, pady=1)


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

            b = Tkinter.Button(ff, text=" ", command=self.ZoneClear)
            b.grid (column=10, row=0, sticky='w', padx=1, pady=1)


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

            b = Tkinter.Button(ff, text="Add Sel", command=self.AddSelRes)
            b.grid (column=0, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="C", command=self.AddSelComp)
            b.grid (column=1, row=0, sticky='w', padx=1)

            b = Tkinter.Label(ff, text=" To Chain")
            b.grid (column=2, row=0, sticky='w', padx=0, pady=1)

            self.addToChain = Tkinter.StringVar(ff)
            self.addToChain.set ( "" )
            e = Tkinter.Entry(ff, width=5, textvariable=self.addToChain)
            e.grid(column=3, row=0, sticky='w', padx=1, pady=1)

            um = Hybrid.Checkbutton(ff, 'at end', False)
            um.button.grid(column = 10, row=0, sticky = 'w', padx=1)
            self.addAtEnd = um.variable

            #l = Tkinter.Label(ff, text=' At' )
            #l.grid(column=10, row=0, sticky='w')

            if 0 :
                self.addAtVar = Tkinter.StringVar(ff)
                self.addAtVar.set ( 'end' )

                c = Tkinter.Radiobutton(ff, text="At End", variable=self.addAtVar, value = 'end')
                c.grid (column=11, row=0, sticky='w')

                c = Tkinter.Radiobutton(ff, text="Pos:", variable=self.addAtVar, value = 'pos')
                c.grid (column=12, row=0, sticky='w')

            um = Hybrid.Checkbutton(ff, 'start at', False)
            um.button.grid(column = 11, row=0, sticky = 'w', padx=1)
            self.addAt = um.variable

            self.addAtPos = Tkinter.StringVar(ff)
            self.addAtPos.set ( "" )
            e = Tkinter.Entry(ff, width=4, textvariable=self.addAtPos)
            e.grid(column=12, row=0, sticky='w', padx=1, pady=1)

            if devMenus :
                um = Hybrid.Checkbutton(ff, 'rename', False)
                um.button.grid(column = 20, row=0, sticky = 'w', padx=1)
                self.renameAdd = um.variable


            #b = Tkinter.Button(ff, text="dif", command=self.DiffSelRes)
            #b.grid (column=11, row=0, sticky='w', padx=1)

            if 0 and devMenus :
                b = Tkinter.Button(ff, text="RC", command=self.RenameChains)
                b.grid (column=30, row=0, sticky='w', padx=1)

            if 0 :
                b = Tkinter.Button(ff, text="M", command=self.ModSel)
                b.grid (column=40, row=0, sticky='w', padx=1)

                b = Tkinter.Button(ff, text="S", command=self.SegToggle)
                b.grid (column=41, row=0, sticky='w', padx=1)


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

            #b = Tkinter.Button(ff, text="S", command=self.AddSheet)
            #b.grid (column=5, row=0, sticky='w', padx=1)

            #b = Tkinter.Button(ff, text="Ca", command=self.CaBlam)
            #b.grid (column=10, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="Seq", command=self.ResSeq)
            b.grid (column=15, row=0, sticky='w', padx=1)

            self.seqPos = Tkinter.StringVar(ff)
            self.seqPos.set ( "" )
            e = Tkinter.Entry(ff, width=4, textvariable=self.seqPos)
            e.grid(column=16, row=0, sticky='w', padx=1, pady=1)


        if devMenus :

            orow += 1
            ff = Tkinter.Frame(cpf)
            ff.grid(column=0, row=orow, sticky='w')

            l = Tkinter.Label(ff, text=' - rotamers:' )
            l.grid(column=0, row=0, sticky='w')

            b = Tkinter.Button(ff, text="R_", command=self.ResRota_)
            b.grid (column=11, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="R0", command=self.ResRota0)
            b.grid (column=12, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="R", command=self.ResRota)
            b.grid (column=13, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="Rc", command=self.ResRota2)
            b.grid (column=14, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="Rx_", command=self.ResRotaX_)
            b.grid (column=15, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="Rx", command=self.ResRotaX)
            b.grid (column=16, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="Ra", command=self.ResRotaAla)
            b.grid (column=17, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="Rq", command=self.ResQ)
            b.grid (column=18, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="Sw", command=self.SwapSelRes)
            b.grid (column=19, row=0, sticky='w', padx=1)

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
            self.addMolName.set ( "NAG" )
            e = Tkinter.Entry(ff, width=5, textvariable=self.addMolName)
            e.grid(column=1, row=0, sticky='w', padx=1, pady=1)

            b = Tkinter.Button(ff, text="+", command=self.AddLigand)
            b.grid (column=2, row=0, sticky='w', padx=1)

            l = Tkinter.Label(ff, text=' Sel:' )
            l.grid(column=5, row=0, sticky='w')

            b = Tkinter.Button(ff, text="Con", command=self.SelLigand)
            b.grid (column=6, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="All", command=self.SelLigands)
            b.grid (column=7, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="Col", command=self.ColorLigands)
            b.grid (column=8, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="QGly", command=self.QGlycans)
            b.grid (column=9, row=0, sticky='w', padx=1)


        if 1 :
            orow += 1
            ff = Tkinter.Frame(cpf)
            ff.grid(column=0, row=orow, sticky='w')

            l = Tkinter.Label(ff, text=' Glycans:' )
            l.grid(column=0, row=0, sticky='w')

            b = Tkinter.Button(ff, text="High Mannose", command=self.AddGlyMan)
            b.grid (column=10, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="Hybrid", command=self.AddGlyHybrid)
            b.grid (column=11, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="Complex", command=self.AddGlyComplex)
            b.grid (column=12, row=0, sticky='w', padx=1)


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

            self.randSearchSize = Tkinter.StringVar(ff, "18")
            e = Tkinter.Entry(ff, width=5, textvariable=self.randSearchSize)
            e.grid(column=2, row=0, sticky='w', padx=1)

            #l = Tkinter.Label(ff, text='step')
            #l.grid(column=3, row=0, sticky='w')

            b = Tkinter.Button(ff, text="Bi", command=self.TorFitBi)
            b.grid (column=8, row=0, sticky='w', padx=1)

            #b = Tkinter.Button(ff, text="<", command=self.TorFitBack)
            #b.grid (column=9, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="Ex", command=self.TorFitEx)
            b.grid (column=11, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="E", command=self.TorFitEnergy)
            b.grid (column=10, row=0, sticky='w', padx=1)

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


        if 1 :

            orow += 1
            ff = Tkinter.Frame(cpf)
            ff.grid(column=0, row=orow, sticky='w')

            l = Tkinter.Label(ff, text=' Rot:')
            l.grid(column=1, row=0, sticky='w')

            self.rotValue = Tkinter.DoubleVar()
            s = Tkinter.Scale (ff, from_=-360, to=360, length=300, showvalue=0, width=15, orient='horizontal', sliderlength = 30, command=self.RotBond, variable=self.rotValue)
            s.grid (column=2, row=0, sticky='ew', padx=1 )
            #self.rotValue = s.get()

            #ff.columnconfigure(2, weight=1)

            b = Tkinter.Button(ff, text="Start", command=self.StartRot)
            b.grid (column=3, row=0, sticky='w', padx=1)


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


            b = Tkinter.Button(ff, text="CpMod", command=self.CpMod)
            b.grid (column=27, row=0, sticky='w', padx=2)


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


    def GetSelAtoms ( self, bbOnly=False ) :

        atoms = []
        for to in self.tree.selection () :

            if to in self.toChain :
                #print " -- Chain:", self.toChain[to]
                for res in self.cur_mol.residues :
                    if res.id.chainId == self.toChain[to] :
                        if bbOnly :
                            atoms.extend ( res.bbAtoms )
                        else :
                            atoms.extend ( res.atoms )
            elif to in self.toRes :
                res = self.toRes[to]
                #print " -- Res: %d.%s.%s" % (res.id.position, res.type, res.id.chainId)
                if bbOnly :
                    atoms.extend ( res.bbAtoms )
                else :
                    atoms.extend ( res.atoms )
            elif to in self.toRess :
                ress = self.toRess[to]
                #print " -- %d res" % len(ress)
                for res in ress :
                    try :
                        if bbOnly :
                            atoms.extend ( res.bbAtoms )
                        else :
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


    def SelectSelBB ( self ) :
        print " - selecting..."

        if self.cur_mol == None :
            umsg ( "Select a molecule first" )
            return

        ats = self.GetSelAtoms ( bbOnly=True )
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
        for c in chimera.selection.currentContents()[0] + chimera.selection.currentContents()[1] :
            if type(c) == _surface.SurfacePiece :
                #print " - sp",
                selSp = c
                if hasattr ( selSp, 'region' ) :
                    selRegs.append ( selSp.region )
                    #selReg = selSp.region
                    #print " - reg: %d" % selSp.region.rid
                #else :
                #    print "?"
            elif type(c) == _molecule.Atom :
                selAt = c
                #print " - atom: %s" % selAt.name
            elif type(c) == _molecule.Bond :
                #print "."
                selAt = c.atoms[0]

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

        if self.cur_mol == None :
            umsg ( "Select a mol..." )
            return

        regs, selAt = self.GetSelRegsAt()

        if len(regs) == 0 :
            umsg ( "No regions selected?" )
            return

        dmap = regs[0].segmentation.seg_map
        rdata = molbuild.RegsData  ( regs )
        rmat = rdata.matrix()
        #nv = VolumeViewer.volume.volume_from_grid_data ( rdata )
        #nv.name = "regs"
        apoints = numpy.zeros ( [1, 3], numpy.float32 )

        mols = []
        #print "atoms near in:"
        #for m in chimera.openModels.list() :
        #    if type(m) == chimera.Molecule and m.display == True :

        for m in [self.cur_mol] :
                print " - in:", m.name
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


    def FindBonds ( self ) :

        print self.cur_mol.name
        mol = self.cur_mol
        allAts = [at for at in mol.atoms if not at.element.name == "H"]
        print "%d atoms" % len(allAts)
        from gridm import Grid
        atGrid = Grid ()
        atGrid.FromAtomsLocal ( allAts, 1.6 )

        for at in mol.atoms :
            if at.element.name == "H" :
                continue
            nearAts = atGrid.AtsNearPtLocal ( at.coord() )
            for nat, v in nearAts :
                if nat == at :
                    continue
                if at in nat.bondsMap :
                    continue

                if at.residue.type in protein3to1 and nat.residue.type in protein3to1 :
                    continue

                nb = mol.newBond ( at, nat )
                print "%s - %s" % ( self.Atom(at), self.Atom(nat) )
                nb.display = nb.Smart
                nb.drawMode = nb.Stick



    def FindBonds2 ( self ) :

        selStr = ""
        for r in chimera.selection.currentResidues() :
            if r.type == "ASN" :
                print "%d\t%s" % (r.id.position, r.id.chainId)
                selStr += "%d.%s," % (r.id.position, r.id.chainId)

        print selStr


    def Atom ( self, at ) :
        return "%s (%s.%d.%s)" % (at.name, at.residue.type, at.residue.id.position, at.residue.id.chainId)



    def Model_0 ( self ) :

        print self.cur_mol.name
        mol = self.cur_mol
        allAts = [at for at in mol.atoms if not at.element.name == "H"]
        print "%d atoms" % len(allAts)
        from gridm import Grid
        atGrid = Grid ()
        atGrid.FromAtomsLocal ( allAts, 4.5 )

        sasRess = []

        for r in chimera.selection.currentResidues () :
            sasRess.append ( r )

        if len(sasRess) == 0 :
            numGlyRes, numProtRes = 0, 0
            for ri, res in enumerate ( mol.residues ) :
                # is a surface residue if at least one atom is a surface atom

                if not res.type in protein3to1 :
                    continue

                if res.id.chainId == "A" :
                    isSurfaceRes = False
                    numProtRes += 1
                    isShielded = False
                    for at in res.atoms :
                        if at.areaSAS > 0.1 :
                            isSurfaceRes = True
                    if isSurfaceRes :
                        sasRess.append ( res )

                if res.type == "ASN" :
                    ndAt = res.atomsMap["ND2"][0]
                    hasGlycan = False
                    for bond in ndAt.bonds :
                        if bond.otherAtom ( ndAt ).name == "C1" :
                            hasGlycan = True
                    if hasGlycan == True :
                        numGlyRes += 1

            print " - %d/%d gly res" % (numGlyRes, numProtRes)


        numSasRess = len(sasRess)
        #print "%d sas res" % numSasRess

        numShielded = 0
        shieldRes = {}

        for res in sasRess :

            isShielded = False
            for at in res.atoms :
                nearAts = atGrid.AtsNearPtLocal ( at.coord() )
                for nat, v in nearAts :
                    if nat.residue.type in [ "NAG", "MAN", "BMA", "FUC", "GAL" ] :
                        shieldRes[nat.residue] = 1
                        isShielded = True
                        break
                    if nat.residue.type in [ "AFUC", "AGAN", "AMAN", "BGAL", "BGLN", "BMAN", "ANE5" ] :
                        shieldRes[nat.residue] = 1
                        isShielded = True
                        break
                #if isShielded :
                #    break

            if isShielded :
                numShielded += 1

        print " - %d sas, %d shielded, %.2f %%" % ( numSasRess, numShielded, 100.0*float(numShielded)/float(numSasRess))

        print " - %d gly res in shield" % len(shieldRes.keys())

        #chimera.selection.clearCurrent ()
        #chimera.selection.addCurrent ( shieldRes.values() )

        for r in self.cur_mol.residues :
            if r.type in [ "NAG", "MAN", "BMA", "FUC", "GAL" ] :
                if not r in shieldRes :
                    for at in r.atoms :
                        at.display = False
            if r.type in [ "AFUC", "AGAN", "AMAN", "BGAL", "BGLN", "BMAN" ] :
                if not r in shieldRes :
                    for at in r.atoms :
                        at.display = False



    #def isSurfPts ( at, atGrid, fabMol=None ) :

    def Model ( self ) :

        print self.cur_mol.name
        mol = self.cur_mol

        SetBBAts ( self.cur_mol )

        #protAts = [at for at in mol.atoms if (not at.element.name == "H" and at.residue.isProt and at.name=="CA")]
        protAts = []
        for res in mol.residues :
            if "CA" in res.atomsMap :
                protAts.append ( res.atomsMap["CA"][0] )
                protAts.append ( res.atomsMap["C"][0] )
                protAts.append ( res.atomsMap["N"][0] )

        glyAts = [at for at in mol.atoms if (not at.element.name == "H" and not at.residue.isProt)]

        print "%d protein (Ca) atoms, %d glycan atoms" % ( len(protAts), len(glyAts) )
        #from gridm import Grid
        import gridm
        reload(gridm)
        atGrid = gridm.Grid ()
        glyGrid = gridm.Grid ()

        import qscores
        vwRad = { "O" : 1.52, "N" : 1.55, "C" : 1.7, "S" : 1.8, "Cl" : 1.75, "H" : 1.2, "P" : 1.8, "Fl" : 1.47 }

        probeRad = 15.0
        #numSphPts = 50000
        numSphPts = 10000
        #atGrid.FromAtomsLocal ( allAts, 5.0 )
        atGrid.FromAtomsLocal ( protAts, probeRad )
        glyGrid.FromAtomsLocal ( glyAts, probeRad )

        from time import time

        print " -- probe rad %.3f, num sphere points %d" % ( probeRad, numSphPts )

        if 0 :
            fabMol = None
            fabCom = None
            for m in chimera.openModels.list() :
                if "FAB" in m.name :
                    fabMol = m
                    from _multiscale import get_atom_coordinates
                    points = get_atom_coordinates ( fabMol.atoms, transformed = True )
                    com = numpy.sum(points, axis=0) / len(points)
                    fabCom = chimera.Point ( com[0], com[1], com[2] )
                    print " - found fab mol: %s, com " % fabMol.name, com

            at = chimera.selection.currentAtoms()[0]
            probePts = qscores.SpherePts ( at.coord(), vwRad[at.element.name] + probeRad, numSphPts )
            print " %s - %.3f + probe rad %.3f" % (at.name, vwRad[at.element.name], probeRad)

            qscores.AddSpherePts ( probePts, (1,0,0,1), 0.2, mname = "Probe points" )

            isSurfacePt = False
            numSph = 0
            sphPts = []
            t0 = time ()
            for pt in probePts :
                isSurfacePt = True

                if 1 :
                    nearAts = atGrid.AtsNearPtLocal ( pt )
                    #print pt, len(nearAts)
                    for nat, v in nearAts :
                        if nat == at :
                            print " - found actual at"
                            continue
                        #print "  %s:%.3f / %.3f" % (nat.name, v.length, vwRad[nat.element.name] + 1.4)
                        #if v.length < vwRad[nat.element.name] + probeRad :
                        isSurfacePt = False
                        #break

                #if isSurfacePt == True :
                #isSurfacePt = not atGrid.AtsInSphere ( pt, probeRad )
                if isSurfacePt :
                    if fabCom :
                        cpt = chimera.Point ( pt[0], pt[1], pt[2] )
                        cpt = mol.openState.xform.apply ( cpt )
                        d = (cpt - fabCom).length
                        sphPts.append ( [d, pt] )
                    numSph += 1
                    #print "-"
                    if not fabCom and numSph == 1 :
                        qscores.AddSpherePts ( [pt], (0,1,0,.3), probeRad, mname = "Surf pt" )

            print ( " - %d allowed spheres, %d, %.2f" % (numSph, len(sphPts), time()-t0) )
            sphPts.sort () # = sorted(sphPts, key=lambda r: r[0], reverse=False)
            d, pt = sphPts[0]
            qscores.AddSpherePts ( [pt], (0,1,0,.3), probeRad, mname = "Surf pt" )
            print " - center %f,%f,%f radius %f" % (pt[0], pt[1], pt[2], probeRad)

            return

        numGlycans = 0
        glyRes = {}
        for res in mol.residues :
            if res.type == "ASN" :
                ndAt = res.atomsMap["ND2"][0]
                hasGlycan = False
                for bond in ndAt.bonds :
                    if bond.otherAtom ( ndAt ).name == "C1" :
                        hasGlycan = True
                if hasGlycan == True :
                    numGlycans += 1
                    glyRes["%s_%d"%(res.id.chainId, res.id.position)] = 1

        if 0 :
            for m in chimera.openModels.list() :
                if m != mol :
                    print " -- %s" % m.name
                    for res in m.residues :
                        if res.type == "ASN" :
                            ndAt = res.atomsMap["ND2"][0]
                            hasGlycan = False
                            for bond in ndAt.bonds :
                                if bond.otherAtom ( ndAt ).name == "C1" :
                                    hasGlycan = True
                            if hasGlycan == True :
                                if not "%s_%d"%(res.id.chainId, res.id.position) in glyRes :
                                    print " - not found in %s - %s_%d" % (mol.name, res.id.chainId, res.id.position)

        protRes = [res for res in mol.residues if res.isProt]
        glyRes = [res for res in mol.residues if not res.isProt]
        print " - %d prot, %d gly res in %d glycans" % (len(protRes), len(glyRes), numGlycans)

        #return

        totAcRes, totGacRes = 0, 0
        t0 = time()
        #for ai, at in enumerate ( mol.atoms ) :
        for ri, res in enumerate ( protRes [0:] ) :
            #if not "CA" in res.atomsMap : continue
            at = res.atomsMap["CA"][0]
            probePts = qscores.SpherePts ( at.coord(), vwRad[at.element.name] + probeRad, numSphPts )
            isAc, isGac = False, False
            for pt in probePts :
                #nearAts = atGrid.AtsNearPtLocal ( pt )
                if 1 :
                    nearAts = atGrid.AtsNearPtLocal ( pt )
                    if len(nearAts) == 0 :
                        isAc = True
                        nearGlyAts = glyGrid.AtsNearPtLocal ( pt )
                        if len(nearGlyAts) == 0 :
                            isGac = True
                            break
                else :
                    if not atGrid.HasNearPtLocal ( pt ) :
                        isAc = True
                        if not glyGrid.HasNearPtLocal ( pt ) :
                            isGac = True
                            break

            if isGac :
                at.residue.ribbonColor = chimera.MaterialColor ( .4, .9, .4, 1.0 )
                totGacRes += 1
                totAcRes += 1
                for at in res.atoms :
                    at.occupancy = 1.0
                    at.bfactor = 1.0
            elif isAc :
                at.residue.ribbonColor = chimera.MaterialColor ( .9, .4, .4, 1.0 )
                totAcRes += 1
                for at in res.atoms :
                    at.occupancy = 0.0
                    at.bfactor = 1.0
            else :
                at.residue.ribbonColor = chimera.MaterialColor ( .7, .7, .7, 1.0 )
                for at in res.atoms :
                    at.occupancy = 0.0
                    at.bfactor = 0.0

            if ri % 10 == 0 :
                t1 = time()
                if ri > 0 :
                    timePerAt = (t1-t0) / float(ri)
                    timeLeft = (len(protRes) - ri) * timePerAt
                    print "%d - %.2f, left %.1f min - %d ac, %d gac, %d prot res" % (ri, (t1-t0)/60.0, timeLeft/60.0, totAcRes, totGacRes, len(protRes))


        print ""
        print "%d protein residues" % len(protRes)
        print "%d saccharides in %d glycans" % ( len(glyRes), numGlycans )

        aden = totAcRes / float( len(protRes) )
        print "%d residues accessible without glycans: %.0f%%" % (totAcRes, aden*100.0)
        gden = totGacRes / float( len(protRes) )
        print "%d residues accessible with glycans: %.0f%%" % (totGacRes, gden*100.0)
        den = totGacRes / float(totAcRes)
        #print " %d acc, %d gacc -> %.2f%%" % (totAcRes, totGacRes, den*100.0)
        #print "%.1f %%" % ( numGlyRes * 100.0/float(numSurfRes) )
        den = (totAcRes-totGacRes)/float(totAcRes)
        print "%d of %d accessible residues shielded: %.0f%%" % (totAcRes-totGacRes, totAcRes, den*100.0  )
        print ""




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
        selResSorted = sorted(selRes, key=lambda r: r.id.position, reverse=False)

        print " - making rmap..."
        rmap = {}
        maxRi = 0
        for r in self.cur_mol.residues :
            if r.id.chainId == toChain :
                maxRi = max ( maxRi, r.id.position )
                rmap[r.id.position] = r

        print " - adding ress..."
        atPosI = None
        keepGaps = False
        #if self.addAtVar.get() == "end" :
        if self.addAt.get() :
            addAtStr = self.addAtPos.get()
            if len(addAtStr) > 0 :
                try :
                    atPosI = int ( addAtStr )
                    print " - starting at pos %d" % atPosI
                except :
                    atPosI = None
                    #print " - keeping same res #"
                    umsg ( "enter a number for start pos" )
                    return
            else :
                umsg ( "enter a number for start pos" )
                return
            if not self.addAtEnd.get() :
                print " - add at %d, keeping gaps" % atPosI
                keepGaps = True
        elif self.addAtEnd.get() :
            atPosI = maxRi + 1
            print " - add at end, pos %d" % atPosI
        else :
            print " - add with same res #"

        firstResI = selResSorted[0].id.position
        print " - first res pos is %d" % firstResI
        for res in selResSorted :

            rid = None
            if keepGaps :
                rid = res.id.position - firstResI + atPosI
            elif atPosI == None :
                # keep residue number, unles already exists...
                rid = res.id.position
            else :
                rid = atPosI
                atPosI += 1

            if rid in rmap :
                print " - trying to add res %d, but alraedy exists in chain, skipping" % rid
                rres = rmap[rid]
                if rres.type != res.type :
                    print " - sel res %s,%d != %s,%d" % (res.type, res.id.position, rres.type, rres.id.position)
                continue


            xf = res.molecule.openState.xform
            xf.premultiply ( toMol.openState.xform.inverse() )
            nres = molref.AddResToMol ( res, toMol, toChain, xf, withoutAtoms=[], rid=rid, asType=asType )

            if nres.type in protein3to1 or nres.type in nucleic3to1 :
                nres.ribbonDisplay = True
                for at in nres.atoms :
                    at.display = False

            status ( "Added res %s.%d.%s to %s chain %s position %d" % (res.type, res.id.position, res.id.chainId, toMol.name, toChain, nres.id.position) )

            rmap[nres.id.position] = nres
            if nres.type in protein3to1 :
                nres.ribbonColor = res.ribbonColor
                if nres.id.position-1 in rmap :
                    pres = rmap[nres.id.position-1]
                    if pres.type in protein3to1 :
                        nb = self.cur_mol.newBond ( pres.atomsMap['C'][0], nres.atomsMap['N'][0] )
                if nres.id.position+1 in rmap :
                    fres = rmap[nres.id.position+1]
                    if fres.type in protein3to1 :
                        nb = self.cur_mol.newBond ( nres.atomsMap['C'][0], fres.atomsMap['N'][0] )
            elif nres.type in nucleic3to1 :
                nres.ribbonColor = res.ribbonColor
                if nres.id.position-1 in rmap :
                    pres = rmap[nres.id.position-1]
                    if pres.type in nucleic3to1 :
                        nb = self.cur_mol.newBond ( pres.atomsMap["O3'"][0], nres.atomsMap['P'][0] )
                if nres.id.position+1 in rmap :
                    fres = rmap[nres.id.position+1]
                    if fres.type in nucleic3to1 :
                        nb = self.cur_mol.newBond ( nres.atomsMap['P'][0], fres.atomsMap["O3'"][0] )


        self.RefreshTree ()





    def AddSelComp ( self ) :

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
        selResSorted = sorted(selRes, key=lambda r: r.id.position, reverse=False)


        print " - finding connections"
        rCons = {}
        def AddCon (r1, r2) :
            if not r1 in rCons : rCons[r1] = {}
            rCons[r1][r2] = 1
            if not r2 in rCons : rCons[r2] = {}
            rCons[r2][r1] = 1

        cBonds = {}
        def AddBond (a1, a2) :
            if not a1 in cBonds : cBonds[a1] = {}
            cBonds[a1][a2] = 1
            if not a2 in cBonds : cBonds[a2] = {}
            cBonds[a2][a1] = 1

        leftRes = {}
        for r in selRes :
            leftRes[r] = 1
            for at in r.atoms :
                for nat in at.neighbors :
                    if nat.residue != r :
                        AddCon ( nat.residue, r )
                        AddBond ( at, nat )


        print " - finding connected components"
        compRes = []
        while len(leftRes) > 0 :
            r1 = leftRes.keys()[0]
            del leftRes[r1]
            Q = rCons[r1].keys()
            ress = [ r1 ]
            while len(Q) > 0 :
                r2 = Q.pop()
                if r2 in leftRes :
                    ress.append ( r2 )
                    del leftRes[r2]
                    Q = Q + rCons[r2].keys()
            compRes.append ( ress )

        chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        def NextCh ( chAt ) :
            i1, i2 = chars.index( chAt[0] ), chars.index( chAt[1] )
            i2 = i2 + 1
            if i2 >= len(chars) :
                i2 = 0
                i1 = i1 + 1
            return chars[i1] + chars[i2]

        atMap = {}
        chAt = "GA"
        print " - %d connected components" % len (compRes)
        for cli, cress in enumerate ( compRes ) :

            print "%d -- %s:" % (cli+1, chAt),
            for ri, res in enumerate ( cress ) :
                if res.type == "ASN" : continue
                print "%d.%s|%s" % (res.id.position, res.id.chainId, res.type),
                xf = res.molecule.openState.xform
                xf.premultiply ( toMol.openState.xform.inverse() )
                nres = molref.AddResToMol ( res, toMol, chAt, xf, withoutAtoms=[], rid=(ri+1) )
                for at in res.atoms : atMap[at] = nres.atomsMap[at.name][0]

            print ""
            chAt = NextCh ( chAt )

        molAts = {}
        for at in toMol.atoms :
            aid = "%s_%d_%s" % (at.residue.id.chainId, at.residue.id.position, at.name)
            molAts [ aid ] = at

        for a1, m in cBonds.iteritems() :
            for a2 in m.keys () :
                if cBonds[a1][a2] == 1 :
                    if a1 in atMap and a2 in atMap :
                        nb = toMol.newBond ( atMap[a1], atMap[a2] )
                        nb.display = nb.Smart
                        nb.drawMode = nb.Stick
                        cBonds[a2][a1] = 0
                    else :
                        if a1 in atMap :
                            aid = "%s_%d_%s" % (a2.residue.id.chainId, a2.residue.id.position, a2.name)
                            if aid in molAts :
                                nb = toMol.newBond ( atMap[a1], molAts[aid] )
                                nb.display = nb.Smart
                                nb.drawMode = nb.Stick
                                cBonds[a2][a1] = 0
                        if a2 in atMap :
                            aid = "%s_%d_%s" % (a1.residue.id.chainId, a1.residue.id.position, a1.name)
                            if aid in molAts :
                                nb = toMol.newBond ( atMap[a2], molAts[aid] )
                                nb.display = nb.Smart
                                nb.drawMode = nb.Stick
                                cBonds[a2][a1] = 0



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



    def BondGo ( self, at0 ) :

        from chimera.resCode import protein3to1, protein1to3

        visAts = { at0:1 }
        Q = [at0]

        while len(Q) > 0 :

            at = Q.pop(0)
            #print "%s " % at.name
            visAts[at] = 1

            if 0 :
                for nat, v in g1.AtsNearPt ( at.xformCoord() ) :
                    if nat.residue.type in protein3to1 :
                        if nat.residue in selResM :
                            return True

            for at2 in at.neighbors :
                #print " -> %s " % (at2.name),
                if not at2 in visAts :
                    if at2.residue.type in protein3to1 :
                        visAts[at2] = 1
                        continue
                    Q.append ( at2 )
                    #print " > "

        return visAts.keys()


    def SelLigand ( self ) :

        selAts = chimera.selection.currentAtoms()

        newSelAts = {}
        bonds = {}
        for at in selAts :
            ats = self.BondGo ( at )
            for at in ats :
                newSelAts[at] = 1
            for at0 in ats :
                for bond in at0.bonds :
                    a1, a2 = bond.atoms
                    if a1.residue.type in protein3to1 :
                        continue
                    if a2.residue.type in protein3to1 :
                        continue
                    bonds[bond] = 1

        chimera.selection.clearCurrent ()
        chimera.selection.addCurrent ( newSelAts.keys() + bonds.keys() )


    def SelLigands ( self ) :


        selAts = chimera.selection.currentAtoms()

        toSel = []
        bonds = {}
        for mol in chimera.openModels.list(modelTypes = [chimera.Molecule]) :
            for res in mol.residues :
                if res.type in protein3to1 or res.type in nucleic3to1 :
                    pass
                else :
                    toSel.extend ( res.atoms )
                    for at in res.atoms :
                        for bond in at.bonds :
                            a1, a2 = bond.atoms
                            if a1.residue.type in protein3to1 or a1.residue.type in nucleic3to1 :
                                continue
                            if a2.residue.type in protein3to1 or a2.residue.type in nucleic3to1 :
                                continue
                            bonds[bond] = 1

        chimera.selection.clearCurrent ()
        chimera.selection.addCurrent ( toSel + bonds.keys() )


    def ColorLigands ( self ) :

        for mol in chimera.openModels.list(modelTypes = [chimera.Molecule]) :
            for res in mol.residues :
                if res.type == "NAG" or res.type == "BGLN" or res.type == "AGLN" :
                    for at in res.atoms :
                        at.color = chimera.MaterialColor (0.0,0.27,0.98)
                elif res.type == "BMA" or res.type == "MAN" :
                    for at in res.atoms :
                        at.color = chimera.MaterialColor (0.20,0.59,0.20)
                elif res.type == "GAL" :
                    for at in res.atoms :
                        at.color = chimera.MaterialColor (1.00,0.97,0.21)
                elif res.type == "FUC" :
                    for at in res.atoms :
                        at.color = chimera.MaterialColor (0.99,0.18,0.10)


    def ConGlys ( self, G, conGly ) :

        visGly = {}
        glys = []

        Q = [G]
        while len(Q) > 0 :
            g = Q.pop()
            visGly[g] = 1
            glys.append ( g )
            for cg in conGly[g].keys() :
                if not cg in visGly and not cg.type == "ASN" :
                    Q.append ( cg )

        return glys


    def QGlycans ( self ) :

        print self.cur_mol.name

        conGly = {}
        def AddCon (r1, r2) :
            if not r1 in conGly :
                conGly[r1] = {}
            conGly[r1][r2] = 1

        for r in self.cur_mol.residues :
            if r.type == "NAG" or r.type == "BMA" or r.type == "GAL" or r.type == "FUC" or r.type == "MAN" :
                for at in r.atoms :
                    for nat in at.neighbors :
                        if nat.residue != r :
                            AddCon ( r, nat.residue )
                            AddCon ( nat.residue, r )

        print "%d gly & asn" % len(conGly.keys())


        import mapq
        mapqDlg = mapq.getdialog()
        sigma = float ( mapqDlg.sigma.get() )
        res = float ( mapqDlg.mapRes.get() )
        expQ, eqn = qscores.ExpectedQScore ( res, sigma )

        print " - map: %s" % self.cur_dmap.name
        print " - mol: %s" % self.cur_mol.name
        print " - sigma: %.2f" % sigma
        print " - res: %.2f, epected Q: %.2f" % (res, expQ)

        minD, maxD = qscores.MinMaxD ( self.cur_dmap )
        print " - mind: %.3f, maxd: %.3f" % (minD, maxD)

        ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]

        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        import gridm
        reload(gridm)
        ptGrid = gridm.Grid()
        ptGrid.FromPoints ( points, 3.0 )
        print " - %d pts grid" % len(points)


        def QForAtoms ( ats, forceCalc = True ) :
            sumQ, sumN = 0.0, 0.0
            for at in ats :
                if at.element.name != "H" :
                    if forceCalc or not hasattr ( at, 'Q' ) :
                        at.Q = qscores.QscorePt3 ( at.coord(), xfI, self.cur_dmap, sigma, ptGrid=ptGrid, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
                    sumQ += at.Q
                    sumN += 1.0
            return sumQ / sumN


        SetBBAts ( self.cur_mol )

        xfI = self.cur_dmap.openState.xform


        fname = self.cur_mol.openedAs[0]
        fname = os.path.splitext ( fname )[0] + "_glycan_Q-scores_.txt"
        print " -> %s" % fname
        fp = open ( fname, "w" )

        for r in self.cur_mol.residues :

            if not r.id.chainId == "A" :
                continue

            if r.type == "ASN" :

                if not r in conGly :
                    # does not have connected glycan
                    continue

                resQ = QForAtoms ( r.atoms )
                #print "%s\t%d\t%.3f\t%.3f" % (r.id.chainId, r.id.position, sumQ/sumN, expQ)
                fp.write ( "%s\t%s\t%d\t%.3f\t%.3f" % (r.type, r.id.chainId, r.id.position, expQ, resQ)  )

                g1 = conGly[r].keys()[0]
                #Q1 = QForAtoms ( g1.atoms )
                #fp.write ( "%s\t%s\t%d\t%.3f\t%.3f\n" % (g1.type, g1.id.chainId, g1.id.position, Q1, expQ)  )

                #if r.id.position % 20 == 0 :
                print "%d.%s" % (r.id.position, r.id.chainId),

                glys = self.ConGlys ( g1, conGly )
                allAts = []
                for g in glys :
                    #print " - %s %d" % (g.type, g.id.position)
                    allAts.extend ( g.atoms )

                Qall = QForAtoms ( allAts )
                fp.write ( "\t%.3f" % Qall  )

                for g in glys :
                    Qg = QForAtoms ( g.atoms )
                    fp.write ( "\t%.3f" % Qg  )

                fp.write ( "\n" )
                #break

        fp.close()
        umsg ( "done -> %s" % fname )





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

        elif 0 :
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



    def AddHighMannose ( self, res ) :

        atN = None
        if res.type == "ASN" : atN = res.atomsMap["ND2"][0]
        elif res.type == "SER" : atN = res.atomsMap["OG"][0]
        else : return

        res = molref.AddGly ( "NAG", atN, None, [] )

        atN = res.atomsMap["O4"][0]
        res = molref.AddGly ( "NAG", atN, None, [] )

        atN = res.atomsMap["O4"][0]
        bma = molref.AddGly ( "BMA", atN, None, [] )

        atN = bma.atomsMap["O3"][0]
        res = molref.AddGly ( "MAN", atN, None, [] )

        atN = bma.atomsMap["O6"][0]
        man = molref.AddGly ( "MAN", atN, None, [] )

        atN = man.atomsMap["O6"][0]
        res = molref.AddGly ( "MAN", atN, None, [] )

        atN = man.atomsMap["O3"][0]
        res = molref.AddGly ( "MAN", atN, None, [] )


    def AddHybrid ( self, res ) :

        atN = None
        if res.type == "ASN" : atN = res.atomsMap["ND2"][0]
        elif res.type == "SER" : atN = res.atomsMap["OG"][0]
        else : return

        res = molref.AddGly ( "NAG", atN, None, [] )

        atN = res.atomsMap["O4"][0]
        res = molref.AddGly ( "NAG", atN, None, [] )

        atN = res.atomsMap["O4"][0]
        bma = molref.AddGly ( "BMA", atN, None, [] )

        atN = bma.atomsMap["O6"][0]
        man = molref.AddGly ( "MAN", atN, None, [] )

        atN = man.atomsMap["O6"][0]
        res = molref.AddGly ( "MAN", atN, None, [] )

        atN = man.atomsMap["O3"][0]
        res = molref.AddGly ( "MAN", atN, None, [] )

        atN = bma.atomsMap["O3"][0]
        man = molref.AddGly ( "MAN", atN, None, [] )

        atN = man.atomsMap["O4"][0]
        nag = molref.AddGly ( "NAG", atN, None, [] )

        atN = nag.atomsMap["O4"][0]
        gal = molref.AddGly ( "GAL", atN, None, [] )



    def AddComplex ( self, res ) :

        atN = None
        if res.type == "ASN" : atN = res.atomsMap["ND2"][0]
        elif res.type == "SER" : atN = res.atomsMap["OG"][0]
        else : return

        res = molref.AddGly ( "NAG", atN, None, [] )

        atN = res.atomsMap["O4"][0]
        res = molref.AddGly ( "NAG", atN, None, [] )

        atN = res.atomsMap["O4"][0]
        bma = molref.AddGly ( "BMA", atN, None, [] )


        if 1 :
            atN = bma.atomsMap["O6"][0]
            man = molref.AddGly ( "MAN", atN, None, [] )

            atN = man.atomsMap["O4"][0]
            nag = molref.AddGly ( "NAG", atN, None, [] )

            atN = nag.atomsMap["O4"][0]
            gal = molref.AddGly ( "GAL", atN, None, [] )

            atN = gal.atomsMap["O4"][0]
            nag = molref.AddGly ( "NAG", atN, None, [] )

            atN = nag.atomsMap["O4"][0]
            gal = molref.AddGly ( "GAL", atN, None, [] )

        if 1 :
            atN = bma.atomsMap["O3"][0]
            man = molref.AddGly ( "MAN", atN, None, [] )

            atN = man.atomsMap["O4"][0]
            nag = molref.AddGly ( "NAG", atN, None, [] )

            atN = nag.atomsMap["O4"][0]
            gal = molref.AddGly ( "GAL", atN, None, [] )

            atN = gal.atomsMap["O4"][0]
            fuc = molref.AddGly ( "FUC", atN, None, [] )

    def AddGlyMan ( self ) :

        ndone = 0
        for r in chimera.selection.currentResidues () :
            if r.type == "ASN" or r.type == "SER" :
                print ""
                print "----------------- %s - %d.%s --------------------" % (r.type, r.id.position, r.id.chainId)
                self.AddHighMannose ( r )
            ndone += 1
            umsg ( "done %d/%d res" % (ndone, len(chimera.selection.currentResidues ())) )


    def AddGlyHybrid ( self ) :

        ndone = 0
        for r in chimera.selection.currentResidues () :
            if r.type == "ASN" or r.type == "SER" :
                print ""
                print "----------------- %s - %d.%s --------------------" % (r.type, r.id.position, r.id.chainId)
                self.AddHybrid ( r )
            ndone += 1
            umsg ( "done %d/%d res" % (ndone, len(chimera.selection.currentResidues ())) )

    def AddGlyComplex ( self ) :

        ndone = 0
        for r in chimera.selection.currentResidues () :
            if r.type == "ASN" or r.type == "SER" :
                print ""
                print "----------------- %s - %d.%s --------------------" % (r.type, r.id.position, r.id.chainId)
                self.AddComplex ( r )
            ndone += 1
            umsg ( "done %d/%d res" % (ndone, len(chimera.selection.currentResidues ())) )



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
        from BuildStructure import changeResidueType

        for i, r in enumerate (ress) :
            if seq[i].upper() in protein1to3 :

                swap(r, protein1to3[seq[i]], preserve=False, bfactor=False)

                for at in r.atoms :
                    at.drawMode = at.EndCap
                    at.display = True # not showRibbon
                    at.color = atomColors[at.element.name if at.element.name in atomColors else " "]

                if r.type == "HIP" :
                	changeResidueType ( r, "HIS" )



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


    def SwapSelRes ( self ) :

        print "swapping _ "

        dmap = self.cur_dmap
        if dmap == None :
            umsg ( "Select a map?" )
            return

        selRes = chimera.selection.currentResidues()
        if len(selRes) == 0 :
            umsg ( "Select residue(s)" )
            return


        mol = selRes[0].molecule
        SetBBAts ( mol )
        ats = []
        for r in mol.residues :
            #ats.extend ( r.bbAtoms )
            ats.extend ( r.atoms )

        import gridm; reload(gridm)
        atGrid = gridm.Grid()
        atGrid.FromAtomsLocal ( ats, 3.0 )

        from molref import ConAtsAtDepth
        conAtsMap = {}
        for r in mol.residues :
            for at in r.atoms :
                atId = "%s_%d_%s" % (r.id.chainId, r.id.position, at.name)
                conAtsMap[atId] = ConAtsAtDepth ( at, 3 )


        selResSorted = sorted(selRes, key=lambda r: r.id.position, reverse=False)

        from molref import ResRotaD
        from mmcif import ColorRes


        from SwapRes import swap, SwapResError
        from mmcif import ColorRes
        from BuildStructure import changeResidueType

        atSeqI = 0
        for res in selResSorted :
            print "%s - %d.%s" % (res.type, res.id.position, res.id.chainId)
            swap ( res, res.type, preserve=False, bfactor=False )
            ColorRes ( res )
            if res.type == "HIP" : changeResidueType ( res, "HIS" )
            atSeqI += 1


    def ResRota_ ( self ) :

        print "rota _ "

        dmap = self.cur_dmap
        if dmap == None :
            umsg ( "Select a map?" )
            return

        selRes = chimera.selection.currentResidues()
        if len(selRes) == 0 :
            umsg ( "Select residue(s)" )
            return


        mol = selRes[0].molecule
        SetBBAts ( mol )
        ats = []
        for r in mol.residues :
            #ats.extend ( r.bbAtoms )
            ats.extend ( r.atoms )

        import gridm; reload(gridm)
        atGrid = gridm.Grid()
        atGrid.FromAtomsLocal ( ats, 3.0 )

        from molref import ConAtsAtDepth
        conAtsMap = {}
        for r in mol.residues :
            for at in r.atoms :
                atId = "%s_%d_%s" % (r.id.chainId, r.id.position, at.name)
                conAtsMap[atId] = ConAtsAtDepth ( at, 3 )


        selResSorted = sorted(selRes, key=lambda r: r.id.position, reverse=False)

        from molref import ResRotaD
        from mmcif import ColorRes

        atSeqI = 0
        for res in selResSorted :
            print "%s - %d.%s" % (res.type, res.id.position, res.id.chainId)
            #ResRotaD ( res, dmap, atGrid, conAtsMap )
            dscore, apos = ResRotaD ( res, dmap, atGrid )
            print " - score: %.4f " % ( dscore  )
            atSeqI += 1


    def ResRota0 ( self ) :

        print "rota..."

        dmap = self.cur_dmap
        if dmap == None :
            umsg ( "Select a map?" )
            return

        selRes = chimera.selection.currentResidues()
        if len(selRes) == 0 :
            umsg ( "Select residue(s)" )
            return


        mol = selRes[0].molecule
        SetBBAts ( mol )
        ats = []
        for r in mol.residues :
            #ats.extend ( r.bbAtoms )
            ats.extend ( r.atoms )

        import gridm; reload(gridm)
        atGrid = gridm.Grid()
        atGrid.FromAtomsLocal ( ats, 3.0 )

        from molref import ConAtsAtDepth
        conAtsMap = {}
        for r in mol.residues :
            for at in r.atoms :
                atId = "%s_%d_%s" % (r.id.chainId, r.id.position, at.name)
                conAtsMap[atId] = ConAtsAtDepth ( at, 3 )


        selResSorted = sorted(selRes, key=lambda r: r.id.position, reverse=False)

        from molref import ResRotaD2
        from mmcif import ColorRes

        atSeqI = 0
        for res in selResSorted :
            print "%s - %d.%s" % (res.type, res.id.position, res.id.chainId)
            #ResRotaD ( res, dmap, atGrid, conAtsMap )
            dscore, clashRess = ResRotaD2 ( res, dmap, atGrid, conAtsMap, checkSideChainClashes=False )
            print " - score: %.4f, %d clash res" % ( dscore, len(clashRess) )
            atSeqI += 1



    def ResRota ( self ) :

        print "rota..."

        dmap = self.cur_dmap
        if dmap == None :
            umsg ( "Select a map?" )
            return

        selRes = chimera.selection.currentResidues()
        if len(selRes) == 0 :
            umsg ( "Select residue(s)" )
            return

        mol = selRes[0].molecule
        SetBBAts ( mol )
        ats = []
        if 0 :
            for r in mol.residues :
                ats.extend ( r.bbAtoms )
                if hasattr ( r, 'optimizedRotamer') :
                    #print "/%d" % r.id.position
                    ats.extend ( r.scAtoms )
        else :
            ats = [at for at in mol.atoms if not at.element.name == "H"]

        from molref import ConAtsAtDepth
        conAtsMap = {}
        for r in mol.residues :
            for at in r.atoms :
                atId = "%s_%d_%s" % (r.id.chainId, r.id.position, at.name)
                conAtsMap[atId] = ConAtsAtDepth ( at, 3 )

        import gridm; reload(gridm)
        atGrid = gridm.Grid()
        atGrid.FromAtomsLocal ( ats, 3.0 )

        selResSorted = sorted(selRes, key=lambda r: r.id.position, reverse=False)

        from molref import ResRotaD2
        from mmcif import ColorRes

        atSeqI = 0
        for res in selResSorted :
            print "%s - %d.%s" % (res.type, res.id.position, res.id.chainId)
            #ResRotaD ( res, dmap, atGrid, conAtsMap )
            dscore, clashRess = ResRotaD2 ( res, dmap, atGrid, conAtsMap, checkSideChainClashes=True )
            print " - score: %.4f, %d clash res" % ( dscore, len(clashRess) )
            atSeqI += 1



    def ResRota2 ( self ) :

        print "rota..."

        dmap = self.cur_dmap
        if dmap == None :
            umsg ( "Select a map?" )
            return

        selRes = chimera.selection.currentResidues()
        if len(selRes) == 0 :
            umsg ( "Select residue(s)" )
            return

        mol = selRes[0].molecule
        SetBBAts ( mol )
        ats = [at for at in mol.atoms if not at.element.name == "H"]

        from molref import ConAtsAtDepth
        conAtsMap = {}
        for r in mol.residues :
            for at in r.atoms :
                atId = "%s_%d_%s" % (r.id.chainId, r.id.position, at.name)
                conAtsMap[atId] = ConAtsAtDepth ( at, 3 )

        import gridm; reload(gridm)
        atGrid = gridm.Grid()
        atGrid.FromAtomsLocal ( ats, 3.0 )

        selResSorted = sorted(selRes, key=lambda r: r.id.position, reverse=False)

        from molref import ResRotaD3
        from mmcif import ColorRes

        atSeqI = 0
        for res in selResSorted :
            print "%s - %d.%s" % (res.type, res.id.position, res.id.chainId)
            dscore = ResRotaD3 ( res, dmap, atGrid, conAtsMap )
            atSeqI += 1
            print " ___  d-score: %.4f (%.4f) ___ " % (dscore, res.dscore)




    def ResRotaAla ( self ) :

        print "rota...ala"

        dmap = self.cur_dmap
        if dmap == None :
            umsg ( "Select a map?" )
            return

        selRes = chimera.selection.currentResidues()
        if len(selRes) == 0 :
            umsg ( "Select residue(s)" )
            return

        mol = selRes[0].molecule
        SetBBAts ( mol )
        ats = [at for at in mol.atoms if not at.element.name == "H"]

        from molref import ConAtsAtDepth
        conAtsMap = {}
        for r in mol.residues :
            for at in r.atoms :
                atId = "%s_%d_%s" % (r.id.chainId, r.id.position, at.name)
                conAtsMap[atId] = ConAtsAtDepth ( at, 3 )

        import gridm; reload(gridm)
        atGrid = gridm.Grid()
        atGrid.FromAtomsLocal ( ats, 3.0 )

        selResSorted = sorted(selRes, key=lambda r: r.id.position, reverse=False)

        from molref import ResRotaD3
        from mmcif import ColorRes

        from SwapRes import swap, SwapResError
        from BuildStructure import changeResidueType

        atSeqI = 0
        for res in selResSorted :
            print "%s - %d.%s" % (res.type, res.id.position, res.id.chainId),
            if res.type != "ALA" :
                print " -> ALA"
                swap ( res, "ALA", preserve=False, bfactor=False )
                ColorRes ( res )
                if res.type == "HIP" : changeResidueType ( res, "HIS" )
            else :
                print ""


    def ResQ ( self ) :

        r = chimera.selection.currentResidues()[0]

        print "-"
        print r.type, r.id.position, r.id.chainId
        print "-"

        rlist = []
        for res in r.molecule.residues :
            if res.type == r.type and res.id.chainId == r.id.chainId :
                #print " - ", res.type, res.id.position, res.id.chainId, res.qResidue
                rlist.append ( [res.qResidue, res] )

        rlist = sorted(rlist, key=lambda x: x[0], reverse=False)
        for qr, res in rlist :
            print " - ", "%.3f"%res.qResidue, res.type, res.id.position, res.id.chainId



    def ResRotaX_ ( self ) :

        print "rota X_"

        dmap = self.cur_dmap
        if dmap == None :
            umsg ( "Select a map?" )
            return

        selRes = chimera.selection.currentResidues()
        if len(selRes) == 0 :
            umsg ( "Select residue(s)" )
            return

        mol = selRes[0].molecule
        SetBBAts ( mol )
        ats = [] # [at for at in mol.atoms if not at.element.name == "H"]

        for r in mol.residues :
            ats.extend ( r.bbAtoms )

        import gridm; reload(gridm)
        atGrid = gridm.Grid()
        atGrid.FromAtomsLocal ( ats, 3.0 )

        selResSorted = sorted(selRes, key=lambda r: r.id.position, reverse=False)

        from molref import ResRotaX
        from molref import ResRotaD
        from mmcif import ColorRes

        log = True
        if len ( selResSorted ) > 1 :
            log = False

        atSeqI = 0
        for res in selResSorted :
            if 0 or not hasattr ( res, 'SCzs' ) :
                print "%s - %d.%s" % (res.type, res.id.position, res.id.chainId)
                #dscore = ResRotaD ( res, dmap, atGrid, conAtsMap=None )
                #try :
                res.SCzs = ResRotaX ( res, dmap, atGrid, setBest=False, log=log )
                #except :
                #    print "?"

                if log :
                    chimera.selection.addCurrent ( res )

            atSeqI += 1
            #print " ___  d-score: %.4f (%.4f) ___ " % (dscore, res.dscore)


        seq = xSeq

        scores, scores2 = [], []
        for i in range ( len(seq) ) :
            score, ii = 0.0, i
            for res in selResSorted :
                if ii < len(seq) :
                    seqResType = protein1to3 [ seq[ii] ]
                    score += res.SCzs [ seqResType ]
                else :
                    break
                ii += 1

            #print "%d - %.3f" % (i+1, score)
            scores.append ( score )
            scores2.append ( [score, i+1] )

        scores2.sort ( reverse=True, key=lambda x: x[0] )
        minSc, maxSc, avgSc, stdSc = min(scores), max(scores), numpy.mean(scores), numpy.std(scores)

        for sc, seqI in scores2[0:20] :
            z = (sc - avgSc) / stdSc
            nn = int ( 30.0 * (sc - minSc) / (maxSc - minSc) ) + 1
            print "%d\t%.3f\t%.3f\t%s" % (seqI, sc, z, "*"*nn)# toSeq[atSeqI],



    def ResRotaX ( self ) :

        print "rota X"

        dmap = self.cur_dmap
        if dmap == None :
            umsg ( "Select a map?" )
            return

        selRes = chimera.selection.currentResidues()
        if len(selRes) == 0 :
            umsg ( "Select residue(s)" )
            return

        mol = selRes[0].molecule
        SetBBAts ( mol )
        ats = [] # [at for at in mol.atoms if not at.element.name == "H"]

        for r in mol.residues :
            ats.extend ( r.bbAtoms )

        import gridm; reload(gridm)
        atGrid = gridm.Grid()
        atGrid.FromAtomsLocal ( ats, 3.0 )

        selResSorted = sorted(selRes, key=lambda r: r.id.position, reverse=False)

        from molref import ResRotaX
        from molref import ResRotaD
        from mmcif import ColorRes

        log = True
        if len ( selResSorted ) > 1 :
            log = False

        atSeqI = 0
        for res in selResSorted :
            if 0 or not hasattr ( res, 'SCzs' ) :
                print "%s - %d.%s" % (res.type, res.id.position, res.id.chainId)
                #dscore = ResRotaD ( res, dmap, atGrid, conAtsMap=None )
                #try :
                res.SCzs = ResRotaX ( res, dmap, atGrid, setBest=False, log=log )
                #except :
                #    print "?"

                if log :
                    chimera.selection.addCurrent ( res )

            atSeqI += 1
            #print " ___  d-score: %.4f (%.4f) ___ " % (dscore, res.dscore)


        seq = xSeq

        scores, scores2 = [], []
        for i in range ( len(seq) ) :

            score, ii = 0.0, i
            for res in selResSorted :
                if ii < len(seq) :
                    seqResType = protein1to3 [ seq[ii] ]
                    score += res.SCzs [ seqResType ]
                else :
                    break
                ii += 1

            #print "%d - %.3f" % (i+1, score)
            scores.append ( score )
            scores2.append ( [score, i+1] )

        scores2.sort ( reverse=True, key=lambda x: x[0] )
        minSc, maxSc, avgSc, stdSc = min(scores), max(scores), numpy.mean(scores), numpy.std(scores)

        for sc, seqI in scores2[0:20] :
            z = (sc - avgSc) / stdSc
            nn = int ( 30.0 * (sc - minSc) / (maxSc - minSc) ) + 1
            print "%d\t%.3f\t%.3f\t%s" % (seqI, sc, z, "*"*nn)# toSeq[atSeqI],




    def ResSeq ( self ) :

        print "..."

        dmap = self.cur_dmap
        if dmap == None :
            umsg ( "Select a map?" )
            return

        selRes = chimera.selection.currentResidues()
        if len(selRes) == 0 :
            #umsg ( "Select residue(s)" )
            #return
            pass

        seq = xSeq

        if 0 :
            rtypes = {}
            for i in range ( len(seq) ) :
                rtypes [ protein1to3 [ seq[i] ] ] = 1

            for rtype in rtypes.keys() :
                print '"%s", ' % ( rtype  ),

            print ""

        startI = 1
        print seq

        if 0 :
            for i in range ( len(seq) ) :
                if i + 3 < len(seq) :
                    if seq[i] == "C" and seq[i+3] == "C" :
                        print i + startI, seq[i:i+4]

        if 0 :
            sstr = "TTCED"
            for i in range ( len(seq) ) :
                if seq[i:i+len(sstr)] == sstr :
                    print " - %s at %d" % (sstr, i+1)



    	from SwapRes import swap, SwapResError
        from mmcif import ColorRes
        from BuildStructure import changeResidueType

        if len ( self.seqPos.get() ) == 0 :
            return

        try :
            startSeqI = int ( self.seqPos.get() )
        except :
            umsg ( "enter an intger for (Seq)uence position" )
            return


        print " - starting at %d" % startSeqI
        #toSeq = "CLAC"
        toSeq = seq[startSeqI-1:]
        print toSeq[0:50]
        #if 1 : return

        if len(selRes) == 0 :
            return

        selRes = sorted(selRes, key=lambda r: r.id.position, reverse=False)

        from molref import ResRota

        atSeqI = 0
        for res in selRes :
            print "%d: %s %d:%d.%s " % (atSeqI, res.type, res.molecule.id, res.id.position, res.id.chainId),

            toType = protein1to3[ toSeq[atSeqI] ]
            #print " -> [%s] %s " % ( toSeq[atSeqI], toType ),

            if toType != res.type :
                #print " - swapping..."
                print " -> %s" % toType # toSeq[atSeqI],
                swap ( res, toType, preserve=False, bfactor=False )
                ColorRes ( res )
                if res.type == "HIP" : changeResidueType ( res, "HIS" )

            else :
                print ""

            ResRota ( res, dmap ) # quicker; remove for more thorough rota search
            atSeqI += 1

        return # quicker; remove for more thorough rota search

        mol = selRes[0].molecule
        SetBBAts ( mol )
        ats = [at for at in mol.atoms if not at.element.name == "H"]
        import gridm; reload(gridm)
        atGrid = gridm.Grid()
        atGrid.FromAtomsLocal ( ats, 3.0 )

        from molref import ConAtsAtDepth
        conAtsMap = {}
        for r in mol.residues :
            for at in r.atoms :
                atId = "%s_%d_%s" % (r.id.chainId, r.id.position, at.name)
                conAtsMap[atId] = ConAtsAtDepth ( at, 3 )

        from molref import ResRotaD3
        for ati, res in enumerate ( selRes ) :
            print "%d/%d: %s %d.%s" % (ati+1, len(selRes), res.type, res.id.position, res.id.chainId)
            #print " - checking rotamers - in %s" % dmap.name
            #ResRotaD3 ( res, dmap )
            dscore = ResRotaD3 ( res, dmap, atGrid, conAtsMap )
            #ColorRes ( res )

        totDScore = 0.0
        for res in selRes :
            totDScore += res.dscore

        print ""
        print "DScore for %d res: %s" % (len(selRes), totDScore)
        print ""


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
            dmap = molbuild.RegsToMap ( selRegs )
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


    def StartRot ( self ) :

        b = chimera.selection.currentBonds()
        if len(b) != 1 :
            umsg ( "select one bond" )
            return

        self.rotBond = b[0]
        self.lastVal = 0.0
        self.rotValue.set(0)

        atsF, atsB = molref.BondAts ( self.rotBond.atoms[0], self.rotBond.atoms[1] )
        if len(atsF) < len(atsB) :
            self.rotAtoms = atsF
        else :
            self.rotAtoms = atsB



    def RotBond ( self, val ) :

        if not hasattr ( self, 'rotBond' ) :
            return

        #print float(val)
        #print "%.1f" % val
        atVal = self.rotValue.get()
        toRot = atVal - self.lastVal
        self.lastVal = atVal
        molref.RotBond ( self.rotBond.atoms[0], self.rotBond.atoms[1], self.rotAtoms, toRot )



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


    def TorFitBack ( self ) :

        mol = chimera.selection.currentAtoms()[0].molecule

        for at in mol.atoms :
            at.setCoord ( at.coord_ )



    def TorFitBi ( self ) :

        print "TorFit - Bi"

        selBonds = chimera.selection.currentBonds()
        print " - %d selected bonds" % len(selBonds)

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
        print " - step size: %.2f ----" % stepSize

        #res = selAt.residue
        #print " - residue %s %d.%s" % (res.type, res.id.position, res.id.chainId)
        ress = chimera.selection.currentResidues()
        print " - %d selected residues" % len(ress)

        #for r in ress :
        #    if r.type in protein3to1 or r.type in nucleic3to1 :
        #        umsg ( "Please select only ligands (no protein or nucleic)" )
        #        #return


        #if len(selBonds) == 0 :
        #    molref.TorFitRand ( ress, dmap, stepSize, parent=self.parent )
        #else :

        molref.TorFitBiSel ( ress, selBonds, dmap, stepSize )


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



    def TorFitEx ( self ) :

        print "TorFit - Exhaustive"

        selRegs, selAt = self.GetSelRegsAt()

        #if selAt == None :
        #    umsg ( "Select an atom" )
        #    return

        dmap = self.cur_dmap
        if len(selRegs) > 0 :
            dmap = molbuild.RegsToMap ( selRegs )
            dmap.delAfter = True
            print " - in regs map: %s" % dmap.name

        stepSize = float(self.randSearchSize.get())
        print " - step size: %.2f ----" % stepSize

        #res = selAt.residue
        #print " - residue %s %d.%s" % (res.type, res.id.position, res.id.chainId)
        ress = chimera.selection.currentResidues()

        print " - %d selected residues" % len(ress)

        #for r in ress :
        #    if r.type in protein3to1 or r.type in nucleic3to1 :
        #        umsg ( "Please select only ligands (no protein or nucleic)" )
        #        #return

        selBonds = chimera.selection.currentBonds()
        print " - %d selected bonds" % len(selBonds)

        if len(selBonds) == 0 :
            molref.TorFitEx ( ress, dmap, stepSize, parent=self.parent )
        else :
            molref.TorFitExSel1 ( ress, selBonds, dmap, stepSize )


        if hasattr ( dmap, 'delAfter' ) :
            chimera.openModels.close ( dmap )



    def TorFitEnergy ( self ) :

        print "TorFit - Energy"

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
        print " - step size: %.2f ----" % stepSize

        #res = selAt.residue
        #print " - residue %s %d.%s" % (res.type, res.id.position, res.id.chainId)
        ress = chimera.selection.currentResidues()

        #for r in ress :
        #    if r.type in protein3to1 or r.type in nucleic3to1 :
        #        umsg ( "Please select only ligands (no protein or nucleic)" )
        #        #return

        molref.TorFitEnergy ( ress, dmap, stepSize, parent=self.parent )

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
            dmap = molbuild.RegsToMap ( selRegs )
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
                #return

        molref.TorFitGrads ( ress, dmap )

        if delAfter :
            chimera.openModels.close ( dmap )



    def TorFitGradsSel ( self ) :

        selRegs, selAts = self.GetSelRegsAts()

        dmap = self.cur_dmap
        if len(selRegs) > 0 :
            dmap = molbuild.RegsToMap ( selRegs )
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
            dmap = molbuild.RegsToMap ( selRegs )
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
            dmap = molbuild.RegsToMap ( selRegs )
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

        if saveFile == True :
            mdir, mfile = os.path.split(self.cur_dmap.data.path)
            dpath = os.path.join ( mdir, nname )
            print " -> ", dpath
            cmap.write_file ( dpath, "mrc" )


    def ZoneClear ( self ) :
        self.zoneMapName.set( "" )

    def CpMod ( self ) :

        print self.cur_mol.name

        chRess = {}
        for r in self.cur_mol.residues :
            if not r.id.chainId in chRess :
                chRess[r.id.chainId] = []
            chRess[r.id.chainId].append ( [r.id.position, r] )


        nmol = chimera.Molecule()
        nmol.name = self.cur_mol.name + " - copied"

        aMap = {}

        for cid, ress in chRess.iteritems() :

            ress.sort()

            for ri, res in ress :
                nres = nmol.newResidue(res.type, chimera.MolResId(cid, res.id.position))
                # print "New res: %s %d" % (nres.id.chainId, nres.id.position)
                for at in res.atoms :
                    nat = nmol.newAtom ( at.name, chimera.Element(at.element.number) )
                    aMap[at] = nat
                    nres.addAtom( nat )
                    #if xf : nat.setCoord ( xf.apply( at.coord() )  )
                    #else :
                    nat.setCoord ( at.coord() )
                    #nat.drawMode = nat.Sphere
                    # todo: handle alt
                    #nat.color = chimera.MaterialColor( clr[0], clr[1], clr[2], 1.0 )
                    #nat.display = False
                    nat.altLoc = at.altLoc
                    nat.occupancy = at.occupancy

                nres.isHelix = res.isHelix
                nres.isHet = res.isHet
                nres.isSheet = res.isSheet
                nres.isStrand = res.isStrand

        for bond in self.cur_mol.bonds :
            a1 = aMap[bond.atoms[0]]
            a2 = aMap[bond.atoms[1]]
            nb = nmol.newBond ( a1, a2 )
            #if a1.display == True and a2.display == True :
            #    nb.display = True
            #    nb.drawMode = nb.Stick
            #nb.display = True

        chimera.openModels.add ( [nmol] )

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
