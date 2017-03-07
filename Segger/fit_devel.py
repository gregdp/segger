from Segger import dev_menus, timing

class Fit_Devel:

    def add_devel_menus(self, fmenu):
    
        if dev_menus:
            fmenu.add_separator()
            for lbl, var, val in (
                #("By principal axes", self.rotaSearch, 0),
                #("By rotation", self.rotaSearch, 1),
                ):
                fmenu.add_radiobutton(label = lbl, variable = var, value = val)
            for lbl, var in (("Chains = models", self.UseAllMods),
                             ):
                fmenu.add_checkbutton(label = lbl, variable = var)
            for lbl, cmd in (
                ("Export fit scores", self.ExportFitScores),
                #("Replicate with Bio-Matrices", self.StrucBioMT),
                ("Sim chain maps", self.SimChainMaps),
                ("Make chain maps", self.StrucChainMaps),
                (" - show chain maps", self.StrucShowChainMaps),
                (" - hide chain maps", self.StrucHideChainMaps),
                (" - close chain maps", self.StrucCloseChainMaps),
                (" - overlap regions", self.StrucChMapsOvRegs),
                (" - close 1 chain map", self.StrucCloseChainMap),
                #(" - delete chain maps", self.StrucDelChainMaps),
                #(" - show chain/regions", self.ShowChRegs),
                ("Segmentation accuracy", self.SegAccuracy),
                ("Fit RMSD", self.FitRMSD),
                #("Sel Chains", self.GetSelectedChains),
                ("Align to selected chain", self.AlignToSel),
                ("Extract proteins", self.ExtractProteins),
                #("Shape match with selected regions", self.SelRegsShapeScore),
                #("Adjust threshold for best match", self.SelRegsOptimizeShapeScore),
                #("Best shape score", self.StrucBestShapeScore),
                ("Group regions", self.StrucGroupRegions),
                ("Next group", self.NextRGroup),
                #("Find group w/selected", self.FindGroupFromSelRegs),
                #("Get fits", self.GetFits),
                #("Next fit", self.NextFit),
                #("Best fit around selected", self.FitMapToRSelGroups),
                #("Fit ALL around selected", self.FitMapsToRegionsAroundSel),
                #("Fit ALL structures", self.FitMapsToRGroups),
                #("Displayed volume", self.StrucMapVolume),
                ("Fit to map (local)", self.FitSMapToDMap),
                ("Zero density map", self.ZeroDMap_with_FMap),
                ("Extract density", self.TakeDMap_with_FMap),
                ("Masked map", self.MaskedMap),
                ("Next bio Matrix", self.NextBioMt),
                ("All bio matrices", self.GoBioMt),
                             ):
                fmenu.add_command(label = lbl, command = cmd)

        

    def ExportFitScores ( self ) :

        num = self.fit_listbox.size()
        if num == 0 :
            umsg ( "No fits to export" )
            return
        
        ccs = []
        for i in range ( num ) :
            le = self.fit_listbox.get ( i )
            toks = le.split ( " " )
            cc = None
            for t in toks :
                try :
                    cc = float ( t )
                    #print cc,
                    break
                except :
                    #print "[" + t + "]",
                    pass
            #print ""
            ccs.append ( cc )

        def save ( okay, dialog, ccs = ccs ):
            if okay:
                paths = dialog.getPaths ( )
                if paths:
                    path = paths[0]
                    f = open ( path, "a" )
                    for cc in ccs :
                        f.write ( "%f\t" % cc )
                    f.write ( "\n" )
                    f.close ()
                    umsg ( "Wrote %d fits to %s" % ( len(ccs), path ) )


        idir = None
        ifile = None

        mol = self.list_fits[0][0].mols[0]
        mmap = self.list_fits[0][1]
        
        if hasattr(mol, 'openedAs'):
            import os.path
            idir, ifile = os.path.split(mol.openedAs[0])
            base, suf = os.path.splitext(ifile)
            map_base, map_suf = os.path.splitext( mmap.name )
            ifile = base + "_fits_in_%s" % map_base
                       
        from OpenSave import SaveModeless
        SaveModeless ( title = 'Save Fit Scores',
                       filters = [('TXT', '*.txt', '.txt')],
                       initialdir = idir, initialfile = ifile, command = save )            

    def SimChainMaps ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        map_name = os.path.splitext ( dmap.name ) [0]
        path = os.path.dirname ( dmap.data.path ) + os.path.sep

        res, grid = None, None

        self.SetResolution()

        try : res = float ( self.simRes.get() )
        except : print "Invalid number entered for resolution:", self.simRes.get(); return

        try : grid = float ( self.simGridSp.get() )
        except : grid = res / 3.0

        mol = getMod ( self.struc.get() )
        if mol == None : print "Structure", self.struc.get(), "not found"; return

        try : mol.chain_colors
        except : mol.chain_colors = RandColorChains ( mol )
            
        print "Simulating %d chain maps for %s, res %.3f, grid %.3f" % (
            len(mol.chain_colors.keys()), mol.name, res, grid)


        dmap.fitted_mols = []

        for cid, clr in mol.chain_colors.iteritems() :

            cname = map_name + "_" + cid + ".mrc"

            mv = getMod ( cname )

            if mv == None :

                sel_str = "#%d:.%s" % (mol.id, cid)
                print "%s [%s]" % (cname, sel_str),

                cmd = "molmap %s %f sigmaFactor 0.187 gridSpacing %f replace false" % ( sel_str, res, grid )
                print " -", cmd
                chimera.runCommand ( cmd )

                for mod in chimera.openModels.list() :
                    ts = mod.name.split()
                    if len(ts) > 1 and mod.name.find("map") >=0 and mod.name.find("res") >=0 :
                        print " - saving to:", path + cname
                        mod.write_file ( path + cname, "mrc" )
                        chimera.openModels.close ( mod )
                        break

                mv = VolumeViewer.open_volume_file ( path + cname )[0]
                print " - loaded:", mv.name

            class FakeMolecule:
                def __init__(self, fmap):
                    self.fmap = fmap

            fmol = FakeMolecule ( mv )
            
            dmap.fitted_mols.append ( fmol )

            #gv.imap = imap
            #gv.name = cname
            #gv.chain_id = cid
            #dmap.chain_maps.append ( gv )

        


    def StrucChainMaps ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        mols = []
        if self.UseAllMods.get() :
            for mod in chimera.openModels.list() :
                if type(mod) == chimera.Molecule :
                    mols.append ( mod )

        else :
            for m in self.StructuresToFit():
                mols.append ( m )
                m.chain_colors = RandColorChains ( m )

        if len(mols) == 0 :
            print "No structures"; return

        self.MakeChainMaps ( mols, dmap )
        print "- %d chain or unit maps" % len ( dmap.chain_maps )



    def StrucShowChainMaps ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        print "%s - showing %d chain maps" % (dmap.name, len(dmap.chain_maps))
        for chm in dmap.chain_maps :
            chm.display = True


    def StrucHideChainMaps ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        print "%s - hiding %d chain maps" % (dmap.name, len(dmap.chain_maps))
        for chm in dmap.chain_maps :
            chm.display = False


    def StrucCloseChainMaps ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return
        
        print "%s - closing %d chain maps" % (dmap.name, len(dmap.chain_maps))

        for chm in dmap.chain_maps :
            chm.close()

        dmap.chain_maps = []



    def StrucDelChainMaps ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        path = os.path.dirname ( dmap.data.path ) + os.path.sep

        print "%s - deleting %d chain maps" % (dmap.name, len(dmap.chain_maps))
        for chm in dmap.chain_maps :
            try : os.remove ( path + chm.name )
            except : print " - could not delete", chm.name
            chm.close()            

        dmap.chain_maps = []


    def StrucChMapsOvRegs ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return
        path = os.path.dirname ( dmap.data.path ) + os.path.sep

        print "%s - %d chain maps" % (dmap.name, len(dmap.chain_maps))

        smod = self.CurrentSegmentation()
        if smod == None : return
        if len(smod.regions) == 0 : print " - no regions!"; return


        for mn in ["2AW4.pdb", "2AVY.pdb"] :
            m = getMod (mn)
            if m :
                for res in m.residues :
                    res.ribbonDisplay = False


        for chm in dmap.chain_maps :
            try : chm.display = False
            except : continue


        for chm in dmap.chain_maps :

            try : chm.display = False
            except : continue

            if chm.chain_id == "A" : continue
            if chm.chain_id == "B" and chm.mol.name == "2AW4.pdb" : continue

            chm.display = True
            print "%s - %s" % (chm.name, chm.chain_id)

            oregs = self.ShowOverlappingRegions ()
            chimera.viewer.viewAll ()


            chm.display = False

            if len(oregs) == 0 :
                print " - NO overlapping regions"
                continue


            ress = chimera.selection.OSLSelection( "#%d:.%s" % (chm.mol.id, chm.chain_id) ).residues()
            for res in ress : res.ribbonDisplay = True

            mname = os.path.splitext ( chm.mol.name )[0]
            fname = mname + "_%s.pdb" % chm.chain_id
            fmol = getMod ( fname )
            if fmol == None :
                fmol = chimera.openModels.open ( path + fname )[0]
                fmol.ch_colors = RandColorChains ( fmol )

            self.struc.set ( fname )
            self.StrucCenter ()
            self.simRes.set ( "14" )
            self.simGridSp.set ( "2" )
            self.GenStrucMap ()
            sim_dmap = self.MoleculeMap()
            sim_dmap.chmap = chm

            self.rotaSearch.set ( 0 )
            self.FitMapToSelRGroup ()
            rmsd = RMSD ( ress, fmol )
            print "\nRMSD: %f " % rmsd,
            tag = chm.chain_id

            msk_dmap = None
            
            if rmsd > 15.0 :
                msk_dmap = self.MaskMapWRegions ()
                self.FitMapToSelRGroup ()
                rmsd = RMSD ( ress, fmol )
                print "\nRMSD: %f " % rmsd,
                tag = chm.chain_id + "m"

            if rmsd > 15.0 :
                self.rotaSearch.set ( 1 )
                self.FitMapToSelRGroup ()
                rmsd = RMSD ( ress, fmol )
                print "\nRMSD: %f " % rmsd,
                tag = chm.chain_id + "r"

            if rmsd > 15.0 :
                self.rotaSearch.set ( 1 )
                self.FitMapToSelRGroup ( msk_dmap )
                rmsd = RMSD ( ress, fmol )
                print "\nRMSD: %f " % rmsd,
                tag = chm.chain_id + "rm"


            if rmsd < 15.0 :

                print "_________ YES _______\n"

                self.SelRegsOptimizeShapeScore ()
                oregs = self.ShowOverlappingRegions ()
                self.SelRegsOptimizeShapeScore ()
                oregs = self.ShowOverlappingRegions ()

                path = os.path.dirname ( dmap.data.path ) + os.path.sep
                log_file = path + "ribo_fits_sms_corr.txt"
                print "SMS log to", log_file
                fpl = open ( log_file, 'a' )
                fpl.write ( "%s %f " % (tag, sim_dmap.sms) )
                fpl.close ()

                if len(oregs) == 1 :
                    oregs[0].placed = True
                if len(oregs) >= 2 :
                    nreg = smod.join_regions ( oregs )

            else :

                print "_________ NO _______\n"


            #sim_dmap.close()
            sim_dmap.display = False

            if msk_dmap : msk_dmap.close()

            self.StrucCloseChainMap()

            #return
            

    def StrucCloseChainMap ( self ) :

        chm = self.MoleculeMap(create = False)
        if chm :
            print " --- closing ", chm.mol.name
            chimera.openModels.close ( [chm.mol] )
            print " --- closing ", chm.chmap.name
            chm.chmap.close()
            print " --- closing ", chm.name
            chm.close()



    def GoBioMt ( self ) :

        fmap = self.MoleculeMap()
        if fmap == None : return False
        try : mol = fmap.mol
        except : "no mol"; return False

        smod = self.CurrentSegmentation()
        if smod == None : return

        smod.bioM_groups = {}

        mol.bio_mt_at = 0

        for r in smod.surfacePieces :
            r.max_ov = 0.0
            r.max_ov_cid = None
            r.max_ov_bioM = None

        smod.p_ctrs = []            

        while self.NextBioMt () :

            self.FitSMapToDMap ()

            smod.bio_mt_at = mol.bio_mt_at
            smod.bioM_groups[mol.bio_mt_at] = self.ShowChRegs ()

            #self.ZeroDMap_with_FMap ()


            if 1 :
                atoms = chimera.selection.OSLSelection( "#%d:.%s" % (fmap.mol.id, 'P') ).atoms()
                points = _multiscale.get_atom_coordinates ( atoms, transformed = False )
                print " - %d points" % len(points)
                com = numpy.sum(points, axis=0) / len(points)
                comv = numpy.array ( [ [ com[0], com[1], com[2] ] ], numpy.float32 )
                _contour.affine_transform_vertices( comv, Matrix.xform_matrix(fmap.openState.xform) )
                _contour.affine_transform_vertices( comv, Matrix.xform_matrix(smod.openState.xform.inverse()) )

                print com[0], com[1], com[2]
                print comv[0][0], comv[0][1], comv[0][2]

                C = chimera.Vector ( comv[0][0], comv[0][1], comv[0][2] )
                smod.p_ctrs.append ( C )

            #return


    def StrucBioMT ( self ) :

        if len(self.struc.get()) == 0 :
            umsg ( "Please select a structure first" )
            return
        mol = getMod ( self.struc.get() )
        if mol == None :
            umsg ( "%s not open" % self.struc.get() )
            return


        mats = BioMatrices ( mol )

        umsg ( "%s - %d matrices" % ( mol.name, len(mats) ) )




    def NextBioMt ( self ) :

        fmap = self.MoleculeMap()
        if fmap == None : return False

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        try : mol = fmap.mol
        except : "no mol"; return False

        try : mol.bio_mts
        except : mol.bio_mts = BioMatrices ( mol )

        print " - %d bio matrices" % len(mol.bio_mts)

        try : mol.bio_mt_at = mol.bio_mt_at + 1
        except : mol.bio_mt_at = 1

        print "_______________ at matrix %d _______________" % mol.bio_mt_at

        try : amat = mol.bio_mts[ mol.bio_mt_at ]
        except : print " - no such matrix"; return False

        fmap.M = fmap.M0 * am_2_M(amat)

        # orthogonalize M
        T = fmap.M; xf = chimera.Xform.xform ( T[0,0], T[0,1], T[0,2], T[0,3], T[1,0], T[1,1], T[1,2], T[1,3], T[2,0], T[2,1], T[2,2], T[2,3], True )
        mt, mr = xf_2_M ( xf )
        fmap.M = mt * mr

        tXO, tXR = xf_2_M ( dmap.openState.xform )
        T = tXO * tXR * fmap.M
        xfA = chimera.Xform.xform ( T[0,0], T[0,1], T[0,2], T[0,3], T[1,0], T[1,1], T[1,2], T[1,3], T[2,0], T[2,1], T[2,2], T[2,3], True )
        fmap.openState.xform = xfA
        fmap.mol.openState.xform = xfA

        try : mol.chain_maps
        except : mol.chain_maps = []
        for chm in mol.chain_maps : chm.openState.xform = xfA

        return True



    def ShowChRegs ( self ) :

        mol = getMod ( self.struc.get() )
        if mol == None : print self.struc.get(), "not open"; return
        print "Structure:", mol.name

        try : chain_colors = mol.chain_colors
        except :  chain_colors = RandColorChains ( mol )
        mol.chain_colors = chain_colors

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        try : mol.chain_maps
        except : self.MakeChainMaps ( mol, dmap )
        print "- %d chain maps" % len ( mol.chain_maps )

        smod = self.CurrentSegmentation()
        if smod == None : return
        print "- %d regions" % len ( smod.surfacePieces )


        # for sp in smod.surfacePieces : sp.display = False

        rgroups = {}
        rgroups['all'] = []
        ov_regs = []

        for chmap in mol.chain_maps :

            # if chmap.chain_id == 'A' or chmap.chain_id == 'B' : continue
            # chmap.display = True

            regs = self.OverlappingRegions ( dmap, chmap, smod )

            ov_regs = ov_regs + regs

            if len(regs) == 0 :
                print " - no regions found";
                return

            # rgroups[chmap.chain_id] = regs

            #sms = self.ShapeMatchScore ( chmap, dmap, regs )
            #print " - chain %s, map overlaps %d regions - %.4f" % ( chmap.chain_id, len(regs), sms )

            for r in regs :
                #r.display = True
                #r.color = chmap.surf_color
                rgroups['all'].append ( r )

            # break

        if 0 :
            # for overlapping regs, which chain do they overlap the most
            for r in ov_regs :
                max_ov = 0.0
                max_ov_chm = None
                for chm in mol.chain_maps :
                    imap = self.MapIndexesInMap ( dmap, chm ) 
                    ipoints = r.points()
                    noverlap = 0
                    for i,j,k in ipoints :
                        if (i,j,k) in imap:
                            noverlap += 1

                    ov = float(noverlap) / r.point_count()
                    if ov > max_ov :
                        max_ov = ov
                        max_ov_chm = chm

                try : rgroups[max_ov_chm.chain_id]
                except : rgroups[max_ov_chm.chain_id] = []
                rgroups[max_ov_chm.chain_id].append ( r )

        return rgroups


    def ExtractProteins ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        path = os.path.dirname ( dmap.data.path ) + os.path.sep

        if len(self.struc.get()) == 0 : umsg ("Please select a structure first"); return
        fmol = getMod ( self.struc.get() )
        if fmol == None : umsg ( "%s not open" % self.struc.get() ); return

        try : fmol.ch_colors
        except : fmol.ch_colors = RandColorChains ( fmol )

        mname = os.path.splitext ( fmol.name )[0]

        print "%s - %d chains" % (fmol.name, len(fmol.ch_colors.keys()))

        for cid, clr in fmol.ch_colors.iteritems () :

            print cid,
            cmol = copyMolChain ( None, fmol, cid, cid, None, clr.rgba() )
            cmol.name = "%s_%s.pdb" % (mname, cid)

            print " - %d atoms" % len(cmol.atoms)

            if len(cmol.atoms) > 0 :
                print " - writing", path + cmol.name
                chimera.PDBio().writePDBfile ( [cmol], path + cmol.name )


        print ""

    def ZeroDMap_with_FMap ( self ) :
        
        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        fmap = self.MoleculeMap()
        if fmap == None : print 'Choose a molecule'; return

        print "Taking %s away from %s" % ( fmap.name, dmap.name )
        rname = dmap.name + '_-_' + fmap.name

        mmc = getMod ( rname )

        if mmc == None :
            mmc = dmap.writable_copy ()
            mmc.name = rname
            print " - cloned", dmap.name
        else :
            print " - found", mmc.name

        self.ZeroOverlappingRegion ( mmc, fmap )





    def ZeroOverlappingRegion ( self, ref_map, mask_map ) :

        thr = mask_map.surface_levels[0]
        mm = mask_map.data.matrix()
        mm = numpy.where ( mm > thr, mm, numpy.zeros_like(mm) )

        nze = numpy.nonzero ( mm )
        nzs = numpy.array ( [nze[2], nze[1], nze[0]] )

        # the copy is needed! otherwise the _contour.afine_transform does not work for some reason
        points = numpy.transpose ( nzs ).astype(numpy.float32)

        print " - %s - %d points above %.3f" % ( mask_map.name, len(points), thr )

        # transform to index reference frame of ref_map
        f1 = mask_map.data.ijk_to_xyz_transform
        f2 = Matrix.xform_matrix ( mask_map.openState.xform )
        f3 = Matrix.xform_matrix ( ref_map.openState.xform.inverse() )
        f4 = ref_map.data.xyz_to_ijk_transform

        tf = Matrix.multiply_matrices( f2, f1 )
        tf = Matrix.multiply_matrices( f3, tf )
        tf = Matrix.multiply_matrices( f4, tf )

        _contour.affine_transform_vertices ( points, tf )


        mm = ref_map.data.matrix()
        n1, n2, n3 = ref_map.data.size[0], ref_map.data.size[1], ref_map.data.size[2]

        for fi, fj, fk in points :

            for i in [ int(numpy.floor(fi)), int(numpy.ceil(fi)) ] :
                for j in [ int(numpy.floor(fj)), int(numpy.ceil(fj)) ] :
                    for k in [ int(numpy.floor(fk)), int(numpy.ceil(fk)) ] :
                        if i < 0 or i >= n3 : continue
                        if j < 0 or j >= n2 : continue
                        if k < 0 or k >= n1 : continue
                        mm[k,j,i] = 0.0



    def TakeDMap_with_FMap ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        fmap = self.MoleculeMap()
        if fmap == None : print 'Choose a molecule'; return

        print "Taking densities from %s with %s" % ( dmap.name, fmap.name )
        rname = fmap.name + '_-_' + dmap.name

        #mmc = fmap.writable_copy ( require_copy = True )
        #mmc.name = rname
        #print " - cloned", fmap.name


        n1, n2, n3 = fmap.data.size[0], fmap.data.size[1], fmap.data.size[2]
        f_points = VolumeData.grid_indices( (n1,n2,n3), numpy.single )  # i,j,k indices
        _contour.affine_transform_vertices( f_points, fmap.data.ijk_to_xyz_transform )

        d_vals = dmap.interpolated_values ( f_points, fmap.openState.xform )
        df_mat = d_vals.reshape( (n3,n2,n1) )

        f_mat = fmap.data.full_matrix()
        f_mask = numpy.where ( f_mat > fmap.surface_levels[0], numpy.ones_like(f_mat), numpy.zeros_like(f_mat) )

        df_mat = df_mat * f_mask

        df_data = VolumeData.Array_Grid_Data ( df_mat, fmap.data.origin, fmap.data.step, fmap.data.cell_angles )
        
        try : df_v = VolumeViewer.volume.add_data_set ( df_data, None )
        except : df_v = VolumeViewer.volume.volume_from_grid_data ( df_data )

        df_v.name = rname
        df_v.openState.xform = fmap.openState.xform



    def NextFit ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        self.fit_at = self.fit_at + 1
        sms, o = self.fits[self.fit_at]

        print "Fit %d, %s, corr %.3f, sms %.3f, regions" % (self.fit_at, o.mname, o.cor, o.sms),

        smod = self.CurrentSegmentation()
        if smod == None : return

        for r in smod.surfacePieces : r.display = False
        for r in o.regs : r.display = True; print r.rid,

        print ""

        fmol = getMod ( o.mname )
        if fmol == None : print "Fitted mol %s not open" % o.mname; return

        fmap = None
        for m in chimera.openModels.list() :
            # TODO: Don't use "centered" to decide what maps to fit.
            try : m.mols[0].centered
            except : continue
            if m.mols[0].name == o.mname :
                fmap = m
                m.display, m.mols[0].display = True, True
            else :
                m.display, m.mols[0].display = False, False

        self.struc.set ( fmap.mols[0].name )

        if fmap == None :
            print "Fitted map corresponding to %s not open" % o.mname; return

        print "Found fitted map", fmap.name                

        fmap.M = o.M
        tXO, tXR = xf_2_M ( dmap.openState.xform )
        T = tXO * tXR * o.M
        xfA = chimera.Xform.xform ( T[0,0], T[0,1], T[0,2], T[0,3], T[1,0], T[1,1], T[1,2], T[1,3], T[2,0], T[2,1], T[2,2], T[2,3] )
        fmap.openState.xform = xfA
        for mol in fmap.mols : mol.openState.xform = xfA


    def BestFits0 ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return
        path = os.path.dirname ( dmap.data.path ) + os.path.sep
        print " - path:", path

        smod = self.CurrentSegmentation()
        if smod == None : return

        for r in smod.surfacePieces : r.display = True

        self.GetFits ()
        regs_claimed = {}
        dmap.fitted_mols = []

        for i, sms_o in enumerate ( self.fits ) :

            sms, o = sms_o

            print "%d - %s - corr %.3f, sm_corr %.3f, regions, dv %.3f" % (
                i+1, o.mname, o.cor, o.sms, o.dv),

            bRegionsClaimed = False
            for r in o.regs :
                print r.rid,
                if regs_claimed.has_key ( r.rid ) :
                    bRegionsClaimed = True
                    print "*",
                else :
                    r.display = False
                    regs_claimed[r.rid] = True

            print ""

            if bRegionsClaimed :
                print " - one or more regions claimed"
                break

            nsp = smod.join_regions ( o.regs )
            self.nRegions.set ( len(smod.surfacePieces) )

            fmol = getMod ( o.mname )
            if fmol == None : print " - fitted mol %s not open" % o.mname; return

            fmap = fitMap ( o.mname )
            if fmap == None : print " - fitted map not open"

            fmap.M = o.M
            tXO, tXR = xf_2_M ( dmap.openState.xform )
            T = tXO * tXR * o.M
            xfA = chimera.Xform.xform ( T[0,0], T[0,1], T[0,2], T[0,3], T[1,0], T[1,1], T[1,2], T[1,3], T[2,0], T[2,1], T[2,2], T[2,3] )
            fmap.openState.xform = xfA
            for mol in fmap.mols : mol.openState.xform = xfA

            if 0 :
                self.saFitMapToPoints ( fmap, nsp.tpoints, dmap )
                fmap.sms = self.ShapeMatchScore ( fmap.fmol.atoms, dmap, [nsp] )

            self.add_fit ( fmap, dmap )




    def AlignToSel ( self ) :

        if len(self.struc.get()) == 0 : print "Please select a structure first"; return
        m = getMod ( self.struc.get() )
        if m == None : print self.struc.get(), "not open"; return

        m_cid = None
        to_cid = None
        m_to = None

        sela = chimera.selection.currentContents()[0]
        for sel_at in sela :
            try : cid = sel_at.residue.id.chainId
            except : pass

            if m_cid == None and sel_at.molecule == m :
                m_cid = cid
                print "Found selected atom in chain %s in %s" % (cid, m.name)
            elif m_to == None and sel_at.molecule != m :
                m_to = sel_at.molecule
                to_cid = cid
                print "Found selected atom in chain %s in %s" % (cid, m_to.name)

        if m_cid == None or m_to == None or to_cid == None :
            print "Please select at least two atoms in the corresponding structures to align"
            return
        
        print "Aligning %s chain %s to %s chain %s" % (m.name, m_cid, m_to.name, to_cid)

        m_seq = m.sequences ( asDict=True ) [m_cid]
        to_seq = m_to.sequences ( asDict=True ) [to_cid]

        print " -- %d residues -- %d residues" % ( len(m_seq.residues), len(to_seq.residues) )

        m_ats, to_ats = AlignChains ( m_seq, to_seq )
        xf, rmsd = chimera.match.matchAtoms ( to_ats, m_ats )

        m_points = _multiscale.get_atom_coordinates ( m_ats, transformed = True )
        f_points = _multiscale.get_atom_coordinates ( to_ats, transformed = True )
        vss = numpy.square ( numpy.subtract ( m_points, f_points ) )
        sums = numpy.sum ( numpy.sum ( vss, axis=1 ) )
        armsd = numpy.sqrt ( sums / float ( len(m_points) ) )

        print " - %d aligned - RMSD: %.4f, RMSD as placed: %.4f" % ( len(m_ats), rmsd, armsd )

        mxf = m_to.openState.xform
        mxf.multiply ( xf )
        m.openState.xform = mxf


        fmap = None
        for om in chimera.openModels.list() :
            try : om_mol = om.mol
            except : continue
            if om_mol == m : fmap = om; break

        if fmap == None :
            print " - no map found for struc"
            return

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        print ""
        umsg ( "Fitting %s and %s into %s" % ( fmap.name, m.name, dmap.name ) )

        fxf = m.openState.xform
        dxf = dmap.openState.xform

        mxf = dxf.inverse()
        mxf.multiply ( fxf ) # fxf = dxf * mxf

        tXO, tXR = xf_2_M ( mxf )
        fmap.M = tXO * tXR

        tXO, tXR = xf_2_M ( dmap.openState.xform )
        T = tXO * tXR * fmap.M
        xfA = chimera.Xform.xform ( T[0,0], T[0,1], T[0,2], T[0,3], T[1,0], T[1,1], T[1,2], T[1,3], T[2,0], T[2,1], T[2,2], T[2,3] )
        fmap.openState.xform = xfA
        for mol in fmap.mols : mol.openState.xform = xfA

        fmap.M0 = fmap.M



    def FitRMSD ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        print "Map: %s" % dmap.name
        try : print " - %d fitted molecules" % len(dmap.fitted_mols)
        except : print " - no fitted molecules found"; return


        if len(self.struc.get()) == 0 : print "Please select a structure first"; return
        m = getMod ( self.struc.get() )
        if m == None : print self.struc.get(), "not open"; return

        # try : m.ch_colors
        # except : m.ch_colors = RandColorChains ( m )
        # print "%s - %d chains -" % ( m.name, len(m.ch_colors) ), m.ch_colors.keys()

        print " - aligning %s (%d chains) to %d fit mols" % (
            m.name, len(m.sequences()), len(dmap.fitted_mols) )


        mseqs = m.sequences()
        fmol_msi = {}
        m_atoms = []
        f_atoms = []

        for msi, m_seq in enumerate ( mseqs ) :

            min_fmi = None
            min_rmsd = 1e99
            min_m_atoms = None
            min_f_atoms = None

            print "Chain %s - %d residues" % ( m_seq.chain, len(m_seq.residues) )

            for fmi, fm in enumerate ( dmap.fitted_mols ) :

                try : fmol_msi[fmi]; continue
                except : pass

                f_seq = fm.sequences()[0]
                print " - %d/%d %s, %d residues" % (
                    fmi+1, len(dmap.fitted_mols), dmap.fitted_mols[fmi].name, len(f_seq.residues) ),

                m_ats, f_ats = AlignChains ( m_seq, f_seq )
                xf, rmsd = chimera.match.matchAtoms ( m_ats, f_ats )
                m_points = _multiscale.get_atom_coordinates ( m_ats, transformed = True )
                f_points = _multiscale.get_atom_coordinates ( f_ats, transformed = True )
                vss = numpy.square ( numpy.subtract ( m_points, f_points ) )
                sums = numpy.sum ( numpy.sum ( vss, axis=1 ) )
                armsd = numpy.sqrt ( sums / float ( len(m_points) ) )

                print " %d aligned - RMSD: %.4f, RMSD as placed: %.4f" % (
                    len(m_ats), rmsd, armsd )

                if armsd < min_rmsd :
                    min_rmsd = armsd
                    min_fmi = fmi
                    min_m_atoms = m_ats
                    min_f_atoms = f_ats

            if min_fmi != None :
                print " - min is %d - %s" % (min_fmi, dmap.fitted_mols[min_fmi].name)
                fmol_msi[min_fmi] = msi
                m_atoms = m_atoms + min_m_atoms
                f_atoms = f_atoms + min_f_atoms
            else :
                print " - no maps left to align!"

        print ""
        
        # m.chain_fitmol = {}
        for fmi, msi in fmol_msi.iteritems () :
            print "Chain %s -- fit mol %s" % ( mseqs[msi].chain, dmap.fitted_mols[fmi].name )
            # m.chain_fitmol [ mseqs[msi].chain ] = dmap.fitted_mols[fmi]


        print ""

        xf, rmsd = chimera.match.matchAtoms ( m_atoms, f_atoms )
        m_points = _multiscale.get_atom_coordinates ( m_atoms, transformed = True )
        f_points = _multiscale.get_atom_coordinates ( f_atoms, transformed = True )
        vss = numpy.square ( numpy.subtract ( m_points, f_points ) )
        sums = numpy.sum ( numpy.sum ( vss, axis=1 ) )
        armsd = numpy.sqrt ( sums / float ( len(m_points) ) )

        print "Total %d residues aligned - RMSD: %.4f, RMSD as placed: %.4f" % (
            len(m_atoms), rmsd, armsd )




    def FitRMSDNoRef ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        print "Map: %s" % dmap.name
        try : print " - %d fitted molecules" % len(dmap.fitted_mols)
        except : print " - no fitted molecules found"; return


        if len(self.struc.get()) == 0 : print "Please select a structure first"; return
        m = getMod ( self.struc.get() )
        if m == None : print self.struc.get(), "not open"; return

        # try : m.ch_colors
        # except : m.ch_colors = RandColorChains ( m )
        # print "%s - %d chains -" % ( m.name, len(m.ch_colors) ), m.ch_colors.keys()

        print " - aligning %s (%d chains) to %d fit mols" % (
            m.name, len(m.sequences()), len(dmap.fitted_mols) )


        mseqs = m.sequences()
        cmap = []
        for msi, m_seq in enumerate ( mseqs ) :

            min_fid = None
            min_rmsd = 1e99
            min_txf = None
            cmap.append ( '?' )

            for fmi, fm in enumerate ( dmap.fitted_mols ) :

                bIncluded = False
                for cfmi in cmap :
                    if cfmi == fmi : bIncluded = True; break
                if bIncluded : continue

                cmap[msi] = fmi
                # print cmap,

                mats,fats = [],[]
                for i in range ( len(cmap) ) :
                    cid = cmap[i]
                    # print cid, fmols[i].name,
                    for ca_posi, cat in ch_cas[cid].iteritems() :
                        mats.append ( cat.coord() )
                        fats.append ( f_ch_cas[i][ca_posi].coord() )

                # print "matching %d to %d ca atoms" % ( len(fats), len(mats) )
                # txf, rmsd = chimera.match.matchAtoms(fats, mats)
                txf, rmsd = chimera.match.matchPositions(fats, mats)
                # print "- rmsd:", rmsd
                if rmsd < min_rmsd :
                    min_rmsd = rmsd
                    min_cid = chains[i2]
                    min_txf = txf

            print "- pos %d is %s, rmsd %f" % (i1, min_cid, min_rmsd),
            cmap[i1] = min_cid
            print cmap
            # return

            for at in m.atoms :
                c = min_txf.apply ( at.coord() )
                at.setCoord ( c )

        print cmap




    def PlaceBestFits ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return
        path = os.path.dirname ( dmap.data.path ) + os.path.sep
        print " - path:", path

        map_name = os.path.splitext ( dmap.name )[0]

        smod = self.CurrentSegmentation()
        if smod == None : return

        self.GetFits ()

        dmap.fitted_mols = []

        scores = []

        for i, sms_o in enumerate ( self.fits ) :

            if i >= int ( self.numFits.get() ) : break

            sms, o = sms_o

            print "%d/%d - %s - corr %.3f, sm_corr %.3f, dv %.3f" % (
                i+1, len(self.fits), o.mname, o.cor, o.sms, o.dv)

            fmol = getMod ( o.mname )
            if fmol == None : print " - fitted mol %s not found" % o.mname; return

            fmap = fitMap ( o.mname )
            if fmap == None : print " - fitted map for %s not found" % fmol.name; return

            fmap.M = o.M
            tXO, tXR = xf_2_M ( dmap.openState.xform )
            T = tXO * tXR * o.M
            xfA = chimera.Xform.xform ( T[0,0], T[0,1], T[0,2], T[0,3], T[1,0], T[1,1], T[1,2], T[1,3], T[2,0], T[2,1], T[2,2], T[2,3] )
            fmap.openState.xform = xfA
            for mol in fmap.mols : mol.openState.xform = xfA

            self.SaveFit ( fmap )
            

            if 0 :
                oregs = self.OverlappingRegions ( dmap, fmap, smod )
                sms = self.ShapeMatchScore ( fmap.mol.atoms, dmap, oregs )
                print " ---- sms %.4f --- " % (sms)
                scores.append ( sms )

            # if joinRegs : smod.join_regions ( regs )


        if 0 :
            log_f = path + map_name + "_fits_acc.txt"
            print "\nAccuracies log file:", log_f

            try : fp = open ( log_f, "a" )
            except : print " - could not open log file"; return
            fp.write ( "%f %f %f " % ( min(scores), max(scores), sum(scores)/float(len(scores))) )
            for sm_score in scores : fp.write ( "%f " % sm_score )
            fp.write ( "\n" )
            fp.close ()


    def GetFits ( self ) :

        smod = self.CurrentSegmentation()

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        dmap_name = os.path.splitext ( dmap.name )[0]
        path = os.path.dirname ( dmap.data.path ) + os.path.sep

        log_f = path + "%s_fits.txt" % (dmap_name)
        print "Log file: %s" % log_f


        class ClusterEntry :
            def __init__ (self, COM, o) :
                self.COM = chimera.Vector ( COM[0], COM[1], COM[2] )
                self.o = o

        class Cluster :
            def __init__ (self) :
                self.entries = []
                self.COM = chimera.Vector (0,0,0)

            def AddEntry ( self, new_e ) :
                self.entries.append  ( [new_e.o.sms,  new_e] )

                # compute the new COM
                self.COM = chimera.Vector (0,0,0)
                for o, e in self.entries :
                    self.COM = self.COM + e.COM
                self.COM = self.COM / float ( len(self.entries) )

        class Fit:
            def __init__ ( self, mname, cor, sms, dv ):
                self.mname = mname
                self.cor = cor
                self.sms = sms
                self.dv = dv

        try : fp = open ( log_f, 'r' )
        except :
            print "Could not open log file:", log_f
            return

        self.clusters = []

        li = 0
        while 1 :
            li = li + 1
            line = fp.readline()
            if len(line) == 0 : break

            n = line.split()

            o = Fit ( n[0], float(n[1]), float(n[2]), float(n[3]) )

            # get M
            mi = 4
            M = numpy.matrix ( [ [ float(n[mi+ 0]), float(n[mi+ 1]), float(n[mi+ 2]), float(n[mi+ 3]) ],
                                 [ float(n[mi+ 4]), float(n[mi+ 5]), float(n[mi+ 6]), float(n[mi+ 7]) ],
                                 [ float(n[mi+ 8]), float(n[mi+ 9]), float(n[mi+10]), float(n[mi+11]) ],
                                 [ 0.0, 0.0, 0.0, 1.0 ]  ] )
            # orthogonalize M
            T = M ; xf = chimera.Xform.xform ( T[0,0], T[0,1], T[0,2], T[0,3], T[1,0], T[1,1], T[1,2], T[1,3], T[2,0], T[2,1], T[2,2], T[2,3], True )
            mt, mr = xf_2_M ( xf )
            o.M = mt * mr

            COM = chimera.Vector ( M[0,3], M[1,3], M[2,3] )
            #print COM

            # get MI
            mi = mi+12
            MI = numpy.matrix ( [ [ float(n[mi+ 0]), float(n[mi+ 1]), float(n[mi+ 2]), float(n[mi+ 3]) ],
                                  [ float(n[mi+ 4]), float(n[mi+ 5]), float(n[mi+ 6]), float(n[mi+ 7]) ],
                                  [ float(n[mi+ 8]), float(n[mi+ 9]), float(n[mi+10]), float(n[mi+11]) ],
                                  [ 0.0, 0.0, 0.0, 1.0 ]  ] )
            # orthogonalize MI
            T = MI ; xf = chimera.Xform.xform ( T[0,0], T[0,1], T[0,2], T[0,3], T[1,0], T[1,1], T[1,2], T[1,3], T[2,0], T[2,1], T[2,2], T[2,3], True )
            mt, mr = xf_2_M ( xf )
            MI = mt * mr

            mi = mi+12

            o.regs = []
            for rmi in range ( mi, len(n) ) :
                ri = int ( n[rmi] )
                o.regs.append(smod.id_to_region[ri])
            # print "%s - corr %.3f, shape score %.3f, regions:" % (o.mname, o.cor, o.sms), o.regs

            cluster = None
            for cl in self.clusters :
                v = cl.COM - COM
                if v.length < 6.0 :
                    cluster = cl

            if cluster == None :
                #print "+++ Creating new cluster: ",
                cluster = Cluster()
                self.clusters.append ( cluster )

            cluster.AddEntry ( ClusterEntry (COM, o) )
            #print "- cluster COM (%f, %f, %f)" % (cluster.COM[0], cluster.COM[1], cluster.COM[2])

            # self.fits.append ( [o.sms, o] )

        fp.close()

        print "____________ %d alignments, %d clusters ____________" % ( li, len(self.clusters) )

        self.fits = []

        for cl in self.clusters :
            cl.entries.sort()
            cl.entries.reverse()
            sms, cle = cl.entries[0]
            self.fits.append ( [cle.o.cor, cle.o] )

        self.fits.sort()
        self.fits.reverse()

        log_f = path + "%s_fits_sorted.txt" % (dmap_name)
        print "Writing fits to: %s" % log_f
        lfp = open ( log_f, "w" )

        for i, sms_o in enumerate ( self.fits ) :
            sms, o = sms_o

            if i < 20 :
                print "%d - %s - correlation %.3f, dv %.3f" % (
                    i+1, o.mname, o.cor, o.dv),
                #for r in o.regs : print r.rid,
                print ""

            lfp.write ( "%d - structure: %s, cross-correlation: %.3f, dVolume: %.3f\n" % (i+1, o.mname, o.cor, o.dv) )

        print ""

        lfp.close ()

        self.fit_at = -1


    def StrucMapVolume ( self ) :

        fmap = self.MoleculeMap()
        if fmap == None : return

        vol = _surface.enclosed_volume ( *(fmap.surfacePieces[0].geometry) )[0]
        print fmap.name + " volume: %f" % vol

        print " - surface levels",fmap.surface_levels
        thr = fmap.surface_levels[0]

        mm = fmap.data.matrix()
        mmab = numpy.where ( mm > thr, numpy.ones_like(mm), numpy.zeros_like(mm) )
        nz = numpy.shape ( numpy.nonzero ( mmab ) )[1]
        vvol = fmap.data.step[0] * fmap.data.step[1] * fmap.data.step[2]
        tvol = vvol * float(nz)
        print "%s - %d above %f, volume %.3f" % (fmap.name, nz, thr, tvol)

        vvol = fmap.data.step[0] * fmap.data.step[1] * fmap.data.step[2]

        tvol = vvol * float(nz)
        print " - %d above %f, volume %.3f" % (nz, thr, tvol)



    def MaskedMap ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "Please select a density map"; return

        fmap = self.MoleculeMap()
        if fmap == None : print 'Choose a molecule'; return


        imap, f_COM, f_bRad = self.MapMaskedMapIndexes ( dmap, fmap, True )


        

    def MapMaskedMapIndexes ( self, ref_map, mask_map, bAddMap=False ) :

        f1 = ref_map.data.ijk_to_xyz_transform
        f2 = xform_matrix( ref_map.openState.xform )
        tf = multiply_matrices( f2, f1 )

        n1, n2, n3 = ref_map.data.size[0], ref_map.data.size[1], ref_map.data.size[2]
        m_points = VolumeData.grid_indices( (n1,n2,n3), numpy.single )  # i,j,k indices
        transform_vertices( m_points, tf )

        a_vals = mask_map.interpolated_values ( m_points, chimera.Xform() )
        mm = a_vals.reshape( (n3,n2,n1) )

        thr = mask_map.surface_levels[0]
        mm = numpy.where ( mm > thr, mm, numpy.zeros_like(mm) )
        print "- in masked map, %d above %.3f" % ( numpy.shape(mm.nonzero())[1], thr )

        if bAddMap :
            mmd = VolumeData.Array_Grid_Data (mm, ref_map.data.origin, ref_map.data.step, ref_map.data.cell_angles, ref_map.data.rotation)
            gv = VolumeViewer.volume_from_grid_data ( mmd )
            gv.name = ref_map.name + "---" + mask_map.name
            #gv.openState.xform = dmap.openState.xform

        nze = numpy.nonzero ( mm )
        print " - index values of %d non-zero points" % len ( nze[0] )

        imap = {}
        for ei, i in enumerate ( nze[0] ) :
            j = nze[1][ei]
            k = nze[2][ei]

            try : mi = imap[i]
            except : mi = {}; imap[i] = mi

            try : mij = mi[j]
            except : mij = {}; mi[j] = mij

            mij[k] = 1


        nzs = numpy.array ( [nze[2], nze[1], nze[0]] )
        points = numpy.cast[numpy.float32] ( numpy.transpose (nzs) )
        transform_vertices( points, f1 )

        com = numpy.sum(points, axis=0) / len(points)
        C = chimera.Vector ( com[0], com[1], com[2] )
        comv = numpy.ones_like ( points ) * com
        points = points - comv
        bRad = numpy.sqrt ( numpy.max ( numpy.sum ( numpy.square (points), 1 ) ) )

        return imap, C, bRad



    def ShapeMaskedCorr ( self, fmap, dmap, regs ) :

        points = numpy.concatenate ( [r.map_points()
                                      for r in regs], axis=0 )

        #tf = xform_matrix ( dmap.openState.xform )
        #transform_vertices(points, tf)

        sg = VolumeData.zone_masked_grid_data ( dmap.data, points, 0.5 )

        #gv = volume_from_grid_data ( sg, None )
        #gv.name = 'regions_mask'
        #gv.openState.xform = dmap.openState.xform

        nz = numpy.shape ( numpy.nonzero ( sg.matrix() ) )[1]
        if len(points) != nz :
            print "mask failed [%d points, %d masked]" % ( len(points), nz )
            return 0.0

        r_mask_matrix = sg.matrix()

        fmapm = fmap.data.matrix()
        f_weights = numpy.ravel(fmapm).astype(numpy.single)

        size = list(fmapm.shape)
        size.reverse()
        f_points = VolumeData.grid_indices(size, numpy.single)        # i,j,k indices

        thr = fmap.surface_levels[0]
        ge = numpy.greater_equal(f_weights, thr)
        f_points = numpy.compress(ge, f_points, 0)
        f_weights = numpy.compress(ge, f_weights)

        from Matrix import multiply_matrices as MM
        from Matrix import xform_matrix as XFM

        p2m_transform = MM ( dmap.data.xyz_to_ijk_transform,
                        MM ( XFM ( dmap.openState.xform.inverse() ),
                        MM ( XFM ( fmap.openState.xform ), fmap.data.ijk_to_xyz_transform ) ) )

        #transform_vertices(f_points, tf)

        m_weights, outside = VolumeData.interpolate_volume_data (
            f_points, p2m_transform, r_mask_matrix, 'linear' )

        o, shape_cor = overlap_and_correlation ( f_weights, m_weights )

        print " * shape map-masked correlation %.3f" % shape_cor

        return shape_cor



    def OptShapeScore ( self, fmap, dmap, regs ) :

        thr_at = fmap.surface_levels[0]
        sms = self.ShapeMatchScore ( fmap, dmap, regs )

        umsg ( "Optimizing threshold for %s - %.2f/%.4f" % (fmap.name, thr_at, sms) )

        while 1 :
            fmap.surface_levels[0] = thr_at + .01
            new_sms = self.ShapeMatchScore ( fmap, dmap, regs )

            if new_sms > sms :
                thr_at = fmap.surface_levels[0]
                sms = new_sms
                #print "- %.2f/%.4f" % (thr_at, sms),
                
                ro = VolumeViewer.volume.Rendering_Options()
                fmap.update_surface ( False, ro )
                for sp in fmap.surfacePieces :
                    v, t = sp.geometry
                    if len(v) == 8 and len(t) == 12 : sp.display = False

                umsg ( "Optimizing threshold for %s - %.2f/%.4f" % (fmap.name, thr_at, sms) )

            else :
                fmap.surface_levels[0] = thr_at
                break

        while 1 :
            fmap.surface_levels[0] = thr_at - .01
            new_sms = self.ShapeMatchScore ( fmap, dmap, regs )
            if new_sms > sms :
                thr_at = fmap.surface_levels[0]
                sms = new_sms
                print "- %.2f/%.4f" % (thr_at, sms),

                ro = VolumeViewer.volume.Rendering_Options()
                fmap.update_surface ( False, ro )
                for sp in fmap.surfacePieces :
                    v, t = sp.geometry
                    if len(v) == 8 and len(t) == 12 : sp.display = False

                umsg ( "Optimizing threshold for %s - %.2f/%.4f" % (fmap.name, thr_at, sms) )

            else :
                fmap.surface_levels[0] = thr_at
                break

        #print "| %.2f/%.4f" % (thr_at, sms)
        return sms



    def OptShapeScoreRegs ( self, fmap, dmap, regs ) :

        thr_at = dmap.surface_levels[0]
        sms = self.ShapeMatchScore ( fmap, dmap, regs )

        print ", optz regs - %.2f/%.4f" % (thr_at, sms),

        while 1 :
            dmap.surface_levels[0] = thr_at + .01
            new_sms = self.ShapeMatchScore ( fmap, dmap, regs )
            if new_sms > sms :
                thr_at = dmap.surface_levels[0]
                sms = new_sms
                print "- %.2f/%.4f" % (thr_at, sms),
            else :
                break

            ro = VolumeViewer.volume.Rendering_Options()
            dmap.update_surface ( False, ro )
            for sp in dmap.surfacePieces :
                v, t = sp.geometry
                if len(v) == 8 and len(t) == 12 : sp.display = False

        while 1 :
            dmap.surface_levels[0] = thr_at - .01
            new_sms = self.ShapeMatchScore ( fmap, dmap, regs )
            if new_sms > sms :
                thr_at = dmap.surface_levels[0]
                sms = new_sms
                print "- %.2f/%.4f" % (thr_at, sms),
            else :
                break

            ro = VolumeViewer.volume.Rendering_Options()
            dmap.update_surface ( False, ro )
            for sp in dmap.surfacePieces :
                v, t = sp.geometry
                if len(v) == 8 and len(t) == 12 : sp.display = False

        dmap.surface_levels[0] = thr_at
        print "| %.2f/%.4f" % (thr_at, sms)
        return sms


    def StrucBestShapeScore ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        fmap = self.MoleculeMap()
        if fmap == None : return

        # first find out which region overlaps the most
        
        smod = self.CurrentSegmentation()
        if smod == None : return
        if len(smod.regions) == 0 : print " - no regions!"; return

        max_overlap_num = 0
        max_overlap_reg = None

        print "Region with highest overlap: ",

        for region in smod.regions :

            noverlap = 0
            for i,j,k in region.points() :
                try : noverlap += imap[k][j][i]
                except : continue

            if noverlap > max_overlap_num :
                max_overlap_num = noverlap
                max_overlap_reg = region

        if max_overlap_num == 0 :
            print "no region overlaps"
            return None

        else :
            print "%d (%d)" % (max_overlap_reg.rid, max_overlap_num)    
            sms = self.ShapeMatchScore ( fmap.fmol.atoms, dmap, [max_overlap_reg] )
            print " - shape match : %f" % sms

        return sms




    def SelRegsShapeScore ( self ) :


        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        smod = self.CurrentSegmentation()
        if smod is None : return

        fmap = self.MoleculeMap()
        if fmap == None:
            return


        regs = smod.selected_regions()
        if len(regs)==0 :
            umsg ( "Please select a region" );
            return

        print "Shape score of %s to %d regions:" % (fmap.name, len(regs)),
        for r in regs :
            print r.rid,
        print ""

        fmap.sms = self.ShapeMatchScore ( fmap.fmol.atoms, dmap, regs )
        ov, fmap.fit_score = map_overlap_and_correlation ( fmap, dmap, True )

        #self.ShapeMaskedCorr ( fmap, dmap, regs )
        #self.OptShapeScore ( fmap, dmap, regs )

        umsg ( "Shape match score: %.3f, cross-correlation: %.3f" % (fmap.sms, fmap.fit_score) )

        path = os.path.dirname ( dmap.data.path ) + os.path.sep
        log_file = path + smod.name + "_fits_sms.txt"
        print "SMS to", log_file

        fpl = open ( log_file, 'a' )
        fpl.write ( "%s %f %f\n" % (fmap.name, fmap.sms, fmap.fit_score) )
        fpl.close ()




    def SelRegsOptimizeShapeScore ( self ) :


        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        smod = self.CurrentSegmentation()
        if smod is None : return

        fmap = self.MoleculeMap()
        if fmap == None:
            return


        regs = smod.selected_regions()
        if len(regs)==0 :
            umsg ( "Please select a region" );
            return

        print "Shape score of %s to %d regions:" % (fmap.name, len(regs)),
        for r in regs :
            print r.rid,
        print ""

        # fmap.sms = self.ShapeMatchScore ( fmap, dmap, regs )
        # self.ShapeMaskedCorr ( fmap, dmap, regs )
        # fmap.sms = self.OptShapeScore ( fmap, dmap, regs )

        umsg ( "Optimized shape match score: %.3f" % fmap.sms )



    def SegAccuracy ( self, tag = "_acc", joinRegs=True ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return
        map_name = os.path.splitext ( dmap.name )[0]

        smod = self.CurrentSegmentation()
        if smod == None : return

        cmaps = dmap.chain_maps
        print "\nMapping %d regions -> %d chain maps..." % (
            len(smod.regions), len(cmaps) )

        if len(smod.regions) == 0 : print " - no regions!"; return
        if len(cmaps) == 0 : print " - no chain maps!"; return

        ch_regs = {}

        for reg_i, region in enumerate ( smod.regions ) :

            max_overlap_num = 0
            max_overlap_cmap = None

            ipoints = region.points()

            for cm in cmaps :

                noverlap = 0
                for i,j,k in ipoints :
                    try : noverlap += cm.imap[k][j][i]
                    except : continue

                if noverlap > max_overlap_num :
                    max_overlap_num = noverlap
                    max_overlap_cmap = cm

            if max_overlap_num == 0 :
                print "%d of %d - %d %d points - " % (reg_i+1, len(smod.regions), region.rid, len(ipoints) )

            else :
                # print "%d/%d - %d %d points - max overlap %d with %s" % (reg_i+1, len(smod.regions), region.rid, len(ipoints), max_overlap_num, max_overlap_cmap.name )
                c = max_overlap_cmap.surf_color
                region.color = c

                try : ch_regs[max_overlap_cmap].append ( region )
                except : ch_regs[max_overlap_cmap] = [ region ]

        thr = dmap.surface_levels[0]
        scores = []

        for cm in cmaps :
            print cm.name,
            if ch_regs.has_key ( cm ) == False :
                print " - no regions!"
                continue

            regs = ch_regs[cm]
            print " - %d regions" % len(regs),

            sms = self.ShapeMatchScore ( cm.atoms, dmap, regs )
            print " - sms %.4f" % (sms)

            if sms > 0.3 and len(regs) < 10 :
                scores.append ( sms )

            else :
                print "sms low: %.3f, nregs: %d" % (sms, len(regs))

            if joinRegs : smod.join_regions ( regs )


        map_name = os.path.splitext ( dmap.name ) [0]
        path = os.path.dirname ( dmap.data.path ) + os.path.sep

        log_f = path + map_name + tag + ".txt"
        print "\nAccuracies log file:", log_f

        res = NumberFromName ( dmap.name, 'r' );
        if res == None : res = 0.0
        print " - resolution: %.0f, %d scores" % ( res, len(scores) )

        try : fp = open ( log_f, "a" )
        except : print " - could not open log file"; return
        fp.write ( "%f %f %f %f " % ( res, min(scores), max(scores), sum(scores)/float(len(scores))) )
        for sm_score in scores : fp.write ( "%f " % sm_score )
        fp.write ( "\n" )
        fp.close ()



    def ShapeMatchScore ( self, atoms, dmap, regs, bPrint=False ) :

        #fmol = fmap.mol
        #print "atoms from", fmol.name
        #points = get_atom_coordinates ( fmol.atoms, transformed = True )

        print "shape match of %d atoms" % len(atoms)
        points = get_atom_coordinates ( atoms, transformed = True )

        tfd = xform_matrix( dmap.openState.xform.inverse() )
        #tff = xform_matrix( fmap.openState.xform.inverse() )
        #tf = multiply_matrices( tff, tfd )
        transform_vertices( points, tfd )

        sg_fmap = VolumeData.zone_masked_grid_data ( dmap.data, points, max(3.0,dmap.data.step[0]) )

        #gv = VolumeViewer.volume_from_grid_data ( sg_fmap )
        #gv.name = 'registered_structure_region'
        #gv.openState.xform = dmap.openState.xform


        points = numpy.concatenate ( [r.map_points() for r in regs], axis=0 )

        sg_dmap = VolumeData.zone_masked_grid_data ( dmap.data, points, dmap.data.step[0] / 2.0 )
        #gv = VolumeViewer.volume.volume_from_grid_data ( sg_dmap )
        #gv.name = 'regions_mask'

        regs_m = sg_dmap.matrix()
        regs_f = sg_fmap.matrix()
        regs_mz = numpy.where ( regs_m > dmap.surface_levels[0], numpy.ones_like(regs_m), numpy.zeros_like(regs_m) )
        regs_fz = numpy.where ( regs_f > dmap.surface_levels[0], numpy.ones_like(regs_f), numpy.zeros_like(regs_f) )
        
        nz = numpy.shape ( numpy.nonzero ( regs_mz ) )[1]
        print "regions %d nonzero, %d points" % (nz, len(points))

        if nz != len(points) :
            print "mask failed - %d of %d pts nonzero" % (nz, len(points))

        nz = numpy.shape ( numpy.nonzero ( regs_fz ) )[1]
        print "struct. %d nonzero, %d atoms" % (nz, len(atoms))

        nz_int =  numpy.shape ( (regs_mz * regs_fz).nonzero () )[1]
        nz_uni =  numpy.shape ( (regs_mz + regs_fz).nonzero () )[1]

        sm_score = float(nz_int) / float (nz_uni)

        if 1 or bPrint : print " - intersection %d, union %d" % (nz_int, nz_uni)

        print " * shape match score %.3f" % sm_score

        return sm_score




    def MakeChainMaps ( self, mols, dmap ) :

        try : print len(dmap.chain_maps), "chain maps so far"
        except : dmap.chain_maps = []


        if len(mols) == 1 :

            mol = mols[0]
            print "Making chain maps for %s in %s:" % (mol.name, dmap.name)

            for cid, clr in mol.chain_colors.iteritems() :

                map_name = os.path.splitext ( dmap.name )[0]
                cname = map_name + "_" + cid
                sel_str = "#%d:.%s" % (mol.id, cid)
                print "%s [%s]" % (cname, sel_str),

                gv = self.SelStrucMap ( sel_str, mol, dmap, clr.rgba() )
                if gv :
                    gv.name = cname
                    gv.chain_id = cid
                    dmap.chain_maps.append ( gv )
                    gv.mol = mol

                #return


        else :

            print "Making struc maps for %d strucs" % len(mols)

            for mol in mols :

                map_name = os.path.splitext ( dmap.name )[0]
                mol_name = os.path.splitext ( mol.name )[0]
                cname = map_name + "_" + mol_name
                sel_str = "#%d" % (mol.id)
                print "%s [%s]" % (cname, sel_str),
                
                clr = (rand(), rand(), rand(), 1.0)
                gv = self.SelStrucMap ( sel_str, mol, dmap, clr )
                gv.name = cname
                gv.chain_id = mol_name
                dmap.chain_maps.append ( gv )


    def SelStrucMap ( self, sel_str, mol, dmap, color ) :

        atoms = chimera.selection.OSLSelection( sel_str ).atoms()

        if len(atoms) == 0 :
            print "- empty"
            return None
        
        points = get_atom_coordinates ( atoms, transformed = True )
        transform_vertices ( points, xform_matrix(dmap.openState.xform.inverse()) )

        step = chimera.Vector ( dmap.data.step[0], dmap.data.step[1], dmap.data.step[2] )
        sg = zone_masked_grid_data ( dmap.data, points, 2.0 )

        sgm = sg.matrix()
        #sgmt = numpy.where ( sgm > 0, sgm, numpy.ones_like(sgm)*1e99 )
        #min_d = numpy.min ( sgmt )

        #sgmt = numpy.where ( sgm > 0, numpy.ones_like(sgm), numpy.zeros_like(sgm) )
        #sgd = VolumeData.Array_Grid_Data (sgmt, dmap.data.origin, dmap.data.step, dmap.data.cell_angles, dmap.data.rotation)

        gv = volume_from_grid_data ( sg )
        gv.openState.xform = dmap.openState.xform

        dvals = dmap.interpolated_values ( points, mol.openState.xform )
        min_d = numpy.min ( dvals ) / 2.0
        min_d = dmap.surface_levels[0]

        print " - min d: %.4f," % min_d,

        gv.region = ( gv.region[0], gv.region[1], [1,1,1] )

        gv.display = True
        gv.surface_levels = [ min_d ]
        #gv.surface_levels = [20.0]
        gv.surface_colors = [ color ]
        gv.surf_color = color
        ro = Rendering_Options()
        ro.surface_smoothing = True
        #ro.smoothing_factor = .25
        #ro.smoothing_iterations = 2
        gv.update_surface ( False, ro )
        gv.atoms = atoms

        for sp in gv.surfacePieces :
            v, t = sp.geometry
            if len(v) == 8 and len(t) == 12 : sp.display = False

        cm = gv.data.matrix()

        nze = numpy.nonzero ( cm )
        print "%d NZ \\" % len ( nze[0] ),

        cmt = numpy.where ( cm > min_d, cm, numpy.zeros_like(cm) )
        nze = numpy.nonzero ( cmt )
        print " %d NZ" % len ( nze[0] )
        gv.d_thr = min_d

        gv.imap = {}
        for ei, i in enumerate ( nze[0] ) :
            j = nze[1][ei]
            k = nze[2][ei]

            try : mi = gv.imap[i]
            except : mi = {}; gv.imap[i] = mi

            try : mij = mi[j]
            except : mij = {}; mi[j] = mij

            mij[k] = 1


        com = numpy.sum(points, axis=0) / len(points)
        gv.COM = chimera.Vector ( com[0], com[1], com[2] )
        comv = numpy.ones_like ( points ) * com
        points_v = points - comv
        gv.bRad = numpy.sqrt ( numpy.max ( numpy.sum ( numpy.square (points_v), 1 ) ) )

        return gv


    def FitSMapToDMap ( self ) :

        fmap = self.MoleculeMap()
        if fmap == None : return

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        print "Fitting %s to %s" % ( fmap.name, dmap.name )

        ov, corr = self.FitMapLocal ( fmap, dmap )

        tXO, tXR = xf_2_M ( dmap.openState.xform )
        T = tXO * tXR * fmap.M
        xfA = chimera.Xform.xform ( T[0,0], T[0,1], T[0,2], T[0,3], T[1,0], T[1,1], T[1,2], T[1,3], T[2,0], T[2,1], T[2,2], T[2,3] )
        fmap.openState.xform = xfA
        for mol in fmap.mols : mol.openState.xform = xfA

        try : fmap.mols.chain_maps
        except : fmap.mol.chain_maps = []
        for chm in fmap.mol.chain_maps : chm.openState.xform = xfA


    def FitMapLocal ( self, fmap, dmap ) :

        f_m = fmap.data.full_matrix(); size = list(f_m.shape); size.reverse()
        fmap.fpoints = grid_indices(size, numpy.single)        # i,j,k indices
        transform_vertices( fmap.fpoints, fmap.data.ijk_to_xyz_transform )
        fmap.fpoint_weights = numpy.ravel(f_m).astype(numpy.single)

        threshold = fmap.surface_levels[0]
        #threshold = .3 * max ( numpy.ravel(f_m).astype(numpy.single) )

        ge = numpy.greater_equal(fmap.fpoint_weights, threshold)
        fmap.fpoints = numpy.compress(ge, fmap.fpoints, 0)
        fmap.fpoint_weights = numpy.compress(ge, fmap.fpoint_weights)
        nz = numpy.nonzero( fmap.fpoint_weights )[0]

        print " - threshold %.4f, %d nonzero" % ( threshold, len(nz) )

        if len(nz) < len (fmap.fpoint_weights) :
            fmap.fpoints = numpy.take( fmap.fpoints, nz, axis=0 )
            fmap.fpoint_weights = numpy.take(fmap.fpoint_weights, nz, axis=0)

        mass = numpy.sum(fmap.fpoint_weights, dtype=numpy.single)
        fmap.rotation_center = numpy.dot(fmap.fpoint_weights,fmap.fpoints) / mass

        olap, cor = self.FitMap_T ( fmap, dmap, num_tries=1 )

        fmap.fpoints = None
        fmap.fpoint_weights = None
        fmap.rotation_center = None

        return olap, cor

    def StrucGroupRegions ( self ) :

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        fmap = self.MoleculeMap()
        if fmap == None : print 'Choose a molecule'; return

        smod = self.CurrentSegmentation()
        if smod == None : return


        regs = smod.selected_regions()

        if len(regs)==0 :
            print "\nGrouping all regions"
            tvol = self.MapVolume ( fmap )
            smod.rgroups = self.GroupAllRegions ( smod, tvol )

        elif len(regs) == 1 :
            tvol = self.MapVolume ( fmap )
            bRad = -1.0 # self.MapBoundingRad ( fmap )

            print "\nMaking groups around region %d - target vol %.3f, b-Rad %.3f" % (regs[0].rid, tvol, bRad)
            smod.rgroups, maxDepthReached = self.GroupAroundReg ( smod, regs[0], tvol, bRad )
            print " - depth reached: %d" % maxDepthReached


        else :
            print "Please select no regions (for global groups) or one region (for local groups)"
            return

        if len ( smod.rgroups ) == 0 :
            print "No groups result!"
            return

        smod.rgroups.sort()
        smod.rgroup_at = 0
        dv, regs = smod.rgroups[smod.rgroup_at]
        print "Group %d of %d - dv %f, regions:" % (smod.rgroup_at+1, len(smod.rgroups), dv),
        for r in regs : print r.rid,
        print ""

        for sp in smod.surfacePieces :
            if regs.count ( sp.region ) > 0 : sp.display = True
            else : sp.display = False



    def NextRGroup ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        smod.rgroup_at = smod.rgroup_at + 1

        if smod.rgroup_at >= len(smod.rgroups) : smod.rgroup_at = 0

        dv, regs = smod.rgroups[smod.rgroup_at]

        print "Group %d/%d - dv %f, regions:" % (smod.rgroup_at+1, len(smod.rgroups), dv),
        for r in regs : print r.rid,
        print ""

        for sp in smod.surfacePieces :
            if regs.count ( sp.region ) > 0 : sp.display = True
            else : sp.display = False



    def FindGroupFromSelRegs ( self ) :

        smod = self.CurrentSegmentation()
        if smod == None : return

        print "\nStructure:",
        fmap = self.MoleculeMap()
        if fmap == None : print 'Choose a molecule'; return
        tvol = self.MapVolume ( fmap )


        regs = smod.selected_regions()        
        if len(regs)==0 : print "no selected regions found"; return

        print "Finding group for %d selected regions:" % ( len(regs) ),
        regs_vol = 0.0
        for r in regs :
            print r.rid,
            regs_vol = regs_vol + r.enclosed_volume()
        print ""

        dv = abs ( tvol - regs_vol ) / tvol;
        print " - total regions volume: %.3f - DV %.5f" % (regs_vol, dv)

        points = numpy.concatenate ( [r.map_points()
                                      for r in regs], axis=0 )

        com = numpy.sum(points, axis=0) / len(points)
        C = chimera.Vector ( com[0], com[1], com[2] )
        comv = numpy.ones_like ( points ) * com
        points = points - comv
        regs_bRad = numpy.sqrt ( numpy.max ( numpy.sum ( numpy.square (points), 1 ) ) )

        print " - COM: %.3f %.3f %.3f, bounding rad %.3f" % ( C[0], C[1], C[2], regs_bRad )


        print "Searching %d groups" % len(smod.rgroups)
        gi = 0
        for dv, gregs in smod.rgroups :
            gi = gi + 1
            if len(regs) != len (gregs) : continue
            allIn = True
            for r in regs :
                if gregs.count ( r ) == 0 : allIn = False; break
            if allIn :
                print " - FOUND! %d - dv %f - regs" % (gi, dv),
                for r in gregs : print r.rid,
                print ""

        print " - done searching"



    def FitMapsToRegionsAroundSel ( self ) :

        print "_______________________________________________________________"

        dmap = segmentation_map()
        if dmap == None : print "No segmentation map"; return

        smod = self.CurrentSegmentation()
        if smod is None : return


        sregs = smod.selected_regions()        
        if len(sregs) != 1 : print "please selected 1 region"; return

        if timing: t0 = clock()

        for sp in smod.surfacePieces :
            if sp.region == sregs[0].region :
                sp.display = True
                clr = sp.region.color
                sp.color = ( clr[0], clr[1], clr[2], REG_OPACITY )
            else :
                sp.display = False


        self.FitOpenMapsToGroupsAround ( smod, sregs[0].region, dmap )

        if timing:
            t1 = clock()
            print "Time: %.1f sec" % (t1 - t0)


    def MapBoundingRad ( self, fmap ) :
        
        fmol = fmap.mol
        points = get_atom_coordinates ( fmol.atoms, transformed = False )

        print "%s (%s)\n - COM: %.3f %.3f %.3f, bounding rad %.3f" % ( fmap.name, fmol.name, C[0], C[1], C[2], bRad )

        return bRad



def fitMap(name):
    # GP - I'm not sure what this function is for - the centerMol function has changed
    # and the code below won't work
    for m in chimera.openModels.list() :
        if hasattr(m, 'mol') and m.mols[0].name == name :
            if not hasattr(m.mols[0], 'centered'):
                centerMol ( m.mols[0] )
            return m
    return None



def NumberFromName ( name, tag ) :

    #name = os.path.splitext ( name ) [0]
    if name.rfind('.mrc') >= 0 : name = name [ : name.rfind('.mrc') ]
    ts = name.split ("_")
    num = None

    for t in ts :
        if t[0:len(tag)] == tag :
            try :
                num = float ( t[len(tag):] )
                break
            except :
                continue

    return num



def BioMatrices ( m ) :

    matrices = {}

    for rm in m.pdbHeaders['REMARK'] :
        s = rm.split()
        if len(s) == 8 and s[2].find('BIOMT') == 0 :
            # print s
            mi = int ( s[3] )
            ri = int ( s[2][5] ) - 1
            try : matrices[mi]
            except : matrices[mi] = numpy.zeros ( (3,4), numpy.float32 )

            matrices[mi][ri][0] = float ( s[4] )
            matrices[mi][ri][1] = float ( s[5] )
            matrices[mi][ri][2] = float ( s[6] )
            matrices[mi][ri][3] = float ( s[7] )

    return matrices                                       



def RMSD ( self, ress, fmol ) :

    N, D, posr = 0.0, 0.0, {}
    for r in ress : posr[r.id.position] = r
    for r in fmol.residues :
        v = r.atomsMap["CA"][0].xformCoord() - posr[r.id.position].atomsMap["CA"][0].xformCoord()
        N, D = N + 1, D + (v.length * v.length)

    return numpy.sqrt ( D / N )


def AlignChains ( ref, match ) :

    import MatchMaker
    import MatchMaker.prefs
    import chimera.misc

    gapOpen = MatchMaker.prefs.defaults[MatchMaker.prefs.GAP_OPEN]
    gapExtend = MatchMaker.prefs.defaults[MatchMaker.prefs.GAP_EXTEND]
    ksdsspCache = set()
    ssFraction = MatchMaker.prefs.defaults[MatchMaker.prefs.SS_MIXTURE]
    ssMatrix = MatchMaker.prefs.defaults[MatchMaker.prefs.SS_SCORES]
    alignKw = { 'ssFraction': ssFraction, 'ssMatrix': ssMatrix, 'computeSS': False }

    # print "Aligning %d res to %d res" % ( len(ref), len(match) )

    score, s1, s2 = MatchMaker.align( ref, match, 'BLOSUM-62', "nw", gapOpen, gapExtend, ksdsspCache, **alignKw)

    refAtoms, matchAtoms = [], []

    for i in range( len(s1) ) :
        if s1[i] == "." or s2[i] == ".": continue
        refRes = s1.residues[s1.gapped2ungapped(i)]
        matchRes = s2.residues[s2.gapped2ungapped(i)]
        refAtom = chimera.misc.principalAtom(refRes)
        if not refAtom: continue
        matchAtom = chimera.misc.principalAtom(matchRes)
        if not matchAtom: continue
        if refAtom.name != matchAtom.name:
            # nucleic P-only trace vs. full nucleic
            if refAtom.name != "P":
                try : refAtom = refAtom.residue.atomsMap["P"][0]
                except KeyError: continue
            else:
                try : matchAtom = matchAtom.residue.atomsMap["P"][0]
                except KeyError: continue
        refAtoms.append(refAtom)
        matchAtoms.append(matchAtom)

    # print " - aligned %d atoms" % len (refAtoms)

    if len ( refAtoms ) < 3 :
        print "too few atoms to perform 3D alignment"
        return None

    return refAtoms, matchAtoms



def copyMolChain ( nmol, mol, cid, new_cid, xf, clr ) :

    if nmol == None :
        nmol = chimera.Molecule()
        nmol.name = "complex"

    print "Copying chain %s (%s) to %s" % (cid, mol.name, new_cid)
    aMap = dict()
    for res in mol.residues :
        if res.id.chainId == cid :
            nres = nmol.newResidue(res.type,
                    chimera.MolResId(new_cid, res.id.position))
            # print "New res: %s %d" % (nres.id.chainId, nres.id.position)
            for at in res.atoms :
                nat = nmol.newAtom (at.name, chimera.Element(at.element.number))
                aMap[at] = nat
                nres.addAtom( nat )
                if xf : nat.setCoord ( xf.apply( at.coord() )  )
                else : nat.setCoord ( at.coord() )
                nat.drawMode = nat.Sphere
                nat.color = chimera.MaterialColor( clr[0], clr[1], clr[2], 1.0 )
                nat.display = True

            nres.isHelix = res.isHelix
            nres.isHet = res.isHet
            nres.isSheet = res.isSheet
            nres.isStrand = res.isStrand
            nres.ribbonDisplay = True
            nres.ribbonDrawMode = 2
            nres.ribbonColor = chimera.MaterialColor( clr[0], clr[1], clr[2], 1.0 );

    for bond in mol.bonds :
        if (bond.atoms[0].residue.id.chainId == cid and
            bond.atoms[1].residue.id.chainId == cid) :
            nb = nmol.newBond ( aMap[bond.atoms[0]], aMap[bond.atoms[1]] )
            nb.display = True

    return nmol



def TransferChainMap ( a, m ) :

    print " - tr %s" % ( a.name ),

    new_v = m.writable_copy()
    new_v.name = a.name + "_tr"
    mm = new_v.data.full_matrix()

    n1, n2, n3 = m.data.size[0], m.data.size[1], m.data.size[2]
    an1, an2, an3 = a.data.size[0], a.data.size[1], a.data.size[2]

    f_i2s = m.data.ijk_to_xyz_transform
    f_t = xform_matrix( m.openState.xform )
    tf = multiply_matrices( f_t, f_i2s )

    m_points = grid_indices( (n1,n2,n3), numpy.single )  # i,j,k indices
    transform_vertices( m_points, tf )

    a_vals = a.interpolated_values ( m_points, chimera.Xform() )
    mm[:,:,:] = a_vals.reshape( mm.shape )

    gv = new_v
    gv.region = ( gv.region[0], gv.region[1], [1,1,1] )

    gv.display = True
    #gv.surface_levels = [ dmap.surface_levels[0] ]
    gv.surface_levels = [ a.d_thr ]
    gv.surf_color = ( rand(), rand(), rand(), 1.0 )
    gv.surface_colors = [ gv.surf_color ]
    ro = Rendering_Options()
    #ro.surface_smoothing = True
    #ro.smoothing_factor = .25
    #ro.smoothing_iterations = 10
    gv.update_surface ( False, ro )

    for sp in gv.surfacePieces :
        v, t = sp.geometry
        if len(v) == 8 and len(t) == 12 : sp.display = False


    mmN = numpy.where ( mm > a.d_thr, mm, numpy.zeros_like(mm) )
    nze = numpy.nonzero ( mmN )
    print "%d NZ" % len ( nze[0] )

    gv.imap = {}
    for ei, i in enumerate ( nze[0] ) :
        j = nze[1][ei]
        k = nze[2][ei]

        try : mi = gv.imap[i]
        except : mi = {}; gv.imap[i] = mi

        try : mij = mi[j]
        except : mij = {}; mi[j] = mij

        mij[k] = 1

    nzs = numpy.array ( [nze[2], nze[1], nze[0]] )
    points = numpy.transpose ( nzs ).astype(numpy.float32)
    transform_vertices ( points, m.data.ijk_to_xyz_transform )

    com = numpy.sum(points, axis=0) / len(points)
    gv.COM = chimera.Vector ( com[0], com[1], com[2] )
    comv = numpy.ones_like ( points ) * com
    points_v = points - comv
    gv.bRad = numpy.sqrt ( numpy.max ( numpy.sum ( numpy.square (points_v), 1 ) ) )

    return gv


def am_2_M ( X ) :

    M = numpy.matrix ( [
        [ X[0,0], X[0,1], X[0,2], X[0,3] ],
        [ X[1,0], X[1,1], X[1,2], X[1,3] ],
        [ X[2,0], X[2,1], X[2,2], X[2,3] ],
        [      0,      0,      0,      1 ]  ] )

    return M


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




    
