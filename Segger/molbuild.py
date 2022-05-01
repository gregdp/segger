
import chimera
import random
import numpy
import _surface
import _contour
import Matrix
import VolumeData
import VolumeViewer
from Matrix import multiply_matrices as MX
import FitMap

from chimera.resCode import nucleic3to1

from chimera.resCode import protein3to1


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



class Helix :

    def __init__ (self) :
        self.mod = None
        self.type = "Helix"
        self.id = random.random()
        self.ress = None


    def Make ( self, regs, mol, chainId, dmap, surfMod ) :

        self.FromRegs ( regs )
        self.MakeMod ( surfMod )
        self.BuildModel ( dmap, mol, chainId, regs )


    def FromRegs ( self, regs ) :

        print "Helix from %d regs" % len(regs)

        if hasattr ( regs[0], 'hx' ) :
            hx = regs[0].hx
            self.switch = hx.switch
            self.sps = hx.sps
            print " - found original helix, switch", self.switch
            try :
                DelRes ( hx.ress )
            except :
                pass


        regs[0].hx = self

        tpoints = numpy.concatenate ( [r.map_points() for r in regs], axis=0 )
        import axes; # reload (axes)
        self.COM, self.U, self.S, self.V = axes.prAxes ( tpoints )

        com = numpy.sum(tpoints, axis=0) / len(tpoints)
        comv = numpy.ones_like ( tpoints ) * com
        points = tpoints - comv
        ppoints = points * self.U
        ex_max = numpy.asarray( numpy.max ( ppoints, 0 ) )[0]
        ex_min = numpy.asarray( numpy.min ( ppoints, 0 ) )[0]

        self.dir = chimera.Vector ( self.U[0,2], self.U[1,2], self.U[2,2] )
        print " - dir:  %.2f %.2f %.2f" % (self.dir[0], self.dir[1], self.dir[2])
        zmax = ex_max[2]
        zmin = -ex_min[2]
        if zmax > zmin :
            move_d = (zmax+zmin)/2.0 - zmin
            #print " move + %.3f" % move_d
            self.COM = self.COM + self.dir * move_d
        elif zmax < zmin :
            move_d = (zmax+zmin)/2.0 - zmax
            #print " move - %.3f" % move_d
            self.COM = self.COM - self.dir * move_d

        #self.Extents = numpy.asarray ( numpy.max ( numpy.abs ( ppoints ), 0 ) )[0]
        self.length = (zmax+zmin)
        self.numRes = int ( numpy.round ( self.length / 1.5 ) )

        print " - length %.2f, ~%d res ~~" % (self.length, self.length/1.5)

        self.heightAdj = 0 # float ( self.helixLength.get() )
        self.width = 2.5 #float ( self.helixWidthF.get() )
        #self.regs = regs

        #if not hasattr ( self, 'switch' ) :
        #    self.switch = False




    def MakeMod ( self, surfMod ) :

        print " - helix make vmod"

        self.length = self.numRes * 1.5

        import axes; reload (axes)

        newMod = _surface.SurfaceModel()
        #newMod = axes.AddCylinderSolid ( chimera.Vector(0,0,0), chimera.Vector(0,0,1), self.length, (1,.2,.2,1), self.width, newMod )
        p = chimera.Vector(0,0,0)
        v = chimera.Vector(0,0,1)
        newMod = axes.AddArrow4 ( p, v, self.length, (1,.2,.2,1), self.width, newMod, 3.5, 4.0 )

        #AddCylinderSolid ( pos, v, d, clr=(0,1,1,1), rad=0.2, mol=None ) :
        #       AddArrow2 ( pos, v, d, clr=(0,1,1,1), rad=0.2, mol=None ) :


        import Matrix
        U = self.U
        COM = self.COM

        R = numpy.array([
                         [  U[0,0], U[0,1], U[0,2], 0.0    ],
                         [  U[1,0], U[1,1], U[1,2], 0.0    ],
                         [  U[2,0], U[2,1], U[2,2], 0.0    ]  ] )

        if hasattr (self, 'switch') and self.switch :
            R = numpy.array([
                             [  -U[0,0], U[0,1], -U[0,2], 0.0    ],
                             [  -U[1,0], U[1,1], -U[1,2], 0.0    ],
                             [  -U[2,0], U[2,1], -U[2,2], 0.0    ]  ] )

        T = Matrix.translation_matrix ( COM )
        Ti = Matrix.translation_matrix ( [0,0,-0.5*self.length] )

        M = MX ( T, MX(R,Ti) )

        #sp = newMod.surfacePieces[-1]
        #v, t = numpy.copy (sp.geometry[0]), numpy.copy(sp.geometry[1])
        #segMod.removePiece ( sp )

        if hasattr ( self, 'sps' ) :
            RemoveSurfPieces ( surfMod, self.sps )

        self.sps = []
        for sp in newMod.surfacePieces :
            v, t = numpy.copy (sp.geometry[0]), numpy.copy(sp.geometry[1])
            _contour.affine_transform_vertices( v, M )
            sp = surfMod.addPiece ( v, t, (1,.2,.2,1) )
            self.sps.append ( sp )
            sp.displayStyle = sp.Mesh
            sp.lineThickness = 2.0





    def BuildModel ( self, dmap, mol, chainId, regs ) :

        print " - build helix %d res, %d regs -- to chain %s" % (self.numRes, len(regs), chainId)
        seq = "C" * self.numRes

        phisPsis = []
        for i in range ( len(seq) ) :
            phisPsis.append ( [-57,-47] )

        startRi = 0
        for r in mol.residues :
            if r.id.position > startRi :
                startRi = r.id.position

        startRi += 5
        print "- mod helix - ri %d" % (startRi)

        self.ress = AddRess ( mol, chainId, seq, startRi, phisPsis, None, None, doRota=False )

        atoms = []
        #apoints = numpy.zeros ( [len(self.ress), 3], numpy.float32 )

        for ri, res in enumerate ( self.ress ) :
            res.isHelix, res.isStrand, res.isSheet, res.ribbonDisplay = True, False, False, True
            res.phi_at, res.phi_eq = -57, -57
            res.psi_at, res.psi_eq = -47, -47
            res.ribbonColor = chimera.MaterialColor ( 0.7, 0.3, 0.3, 1.0 )
            atoms.extend ( res.atoms )

        from _multiscale import get_atom_coordinates
        apoints = get_atom_coordinates ( atoms, transformed = False )

        import axes
        com, U, S, V = axes.prAxes ( apoints )
        #print "- mod helix COM:  %.2f %.2f %.2f" % (com[0], com[1], com[2])

        dir = chimera.Vector ( U[0,2], U[1,2], U[2,2] )

        # match dir to self.dir, com to self.COM
        Acom = Matrix.translation_matrix ( com*-1.0 )
        Hcom = Matrix.translation_matrix ( self.COM )

        rdata, rmat = None, None
        if 0 :
            print " - full reg map from:", dmap.name
            zoneR = 3.0 # numpy.min(dmap.data.step)
            rpoints = numpy.concatenate ( [r.map_points() for r in regs], axis=0 ).astype ( numpy.float32 )
            rdata = VolumeData.zone_masked_grid_data ( dmap.data, rpoints, zoneR )
            rmat = rdata.matrix()
        else :
            print " - part reg map from:", dmap.name
            rdata = RegsData  ( regs )
            rmat = rdata.matrix()


        #print rmat
        ##gdata = VolumeData.Array_Grid_Data ( ndata.full_matrix(), dmap.data.origin, dmap.data.step, dmap.data.cell_angles, name = "atom masked" )
        #nv = VolumeViewer.volume.volume_from_grid_data ( rdata )
        #nv.name = "helix mask small"

        #rmap_values = nv.interpolated_values ( apoints, mol.openState.xform )
        #olap, corr, other = overlap_and_correlation ( rpoint_weights, rmap_values )
        #avg0a = numpy.average ( rmap_values )

        maxXf = None
        maxAvg = -1e9
        useSwitch = None

        apoints0 = numpy.copy ( apoints )

        minD = numpy.min ( rmat [ numpy.where ( rmat > 0 ) ] )
        #minD = numpy.min ( rmat )
        print "min-d:", minD

        maxD = numpy.max ( rmat [ numpy.where ( rmat > 0 ) ] )
        print "max-d:", maxD

        switches = [False, True]
        if hasattr ( self, 'switch' ) :
            print " - old switch : ", self.switch
            self.switch = not self.switch
            print " - new switch : ", self.switch
            switches = [self.switch]

        self.dir = self.dir * -1.0 # to make ress increase in direction of arrow

        #N = 1
        #for angi in range (N) :
        for switch in switches :

            #R2 = Matrix.rotation_transform ( self.dir, angi * (360.0/N) )
            #M = MX(Hcom,MX(R2, MX(R,Acom)))

            if switch : dir = dir * -1.0
            rvec = chimera.cross ( dir, self.dir )
            ang = numpy.arccos ( dir*self.dir / (dir.length * self.dir.length) ) * 180.0 / numpy.pi
            #print " - rotation axis: %.2f, %.2f, %.2f -- %.2f" % (rvec[0], rvec[1], rvec[2], ang)
            R = Matrix.rotation_transform ( rvec, ang )
            M = MX(Hcom,MX(R,Acom))

            apoints = numpy.copy ( apoints0 )
            _contour.affine_transform_vertices ( apoints, M )

            values, outside = VolumeData.interpolate_volume_data ( apoints, rdata.xyz_to_ijk_transform, rmat )
            #avg0 = len ( numpy.where ( values >= minD )[0] ) / float( len(apoints) )
            avg0_ = numpy.average ( values )

            xf1 = FitPointsToData ( apoints, rdata )
            _contour.affine_transform_vertices ( apoints, Matrix.xform_matrix(xf1) )
            #print xf1

            values, outside = VolumeData.interpolate_volume_data ( apoints, rdata.xyz_to_ijk_transform, rdata.matrix() )
            avg1_ = numpy.average ( values )
            #print " -- %d -- in %.4f -> %.4f, avgd %.4f -> %.4f" % (angi, avg0, avg1, avg0_, avg1_)
            print " -- switch", switch, " %.3f -> %.3f" % (avg0_, avg1_)

            if maxXf == None or avg1_ > maxAvg :
                maxAvg = avg1_
                xf1.multiply ( Matrix.chimera_xform(M) )
                maxXf = xf1
                self.switch = switch

        print " -- max score: %.4f" % maxAvg, "switch", self.switch

        for res in self.ress :
            for at in res.atoms :
                at.setCoord ( maxXf.apply (at.coord()) )

        _contour.affine_transform_vertices ( apoints0, Matrix.xform_matrix(maxXf) )
        #values, outside = VolumeData.interpolate_volume_data ( apoints0, rdata.xyz_to_ijk_transform, rdata.matrix() )
        #avg1 = len ( numpy.where ( values >= minD )[0] ) / float( len(apoints0) )
        #print " -- max AI ?: %.3f" % avg1

        SetDrawMode ( self.ress )







class Sheet :

    def __init__ (self) :
        self.mod = None
        self.type = "Sheet"
        self.id = random.random()
        self.width = 1.0
        self.ress = None


    def FromRegs ( self, sregs, numStrands ) :

        self.numStrands = numStrands

        reg = sregs[0]
        #tpoints = reg.map_points()
        tpoints = numpy.concatenate ( [r.map_points() for r in sregs], axis=0 )

        import axes; # reload (axes)
        self.COM, self.U, self.S, self.V = axes.prAxes ( tpoints )

        print "Region %d - %f %f %f, %d points" % (reg.rid, self.S[0], self.S[1], self.S[2], len(tpoints) )

        com = numpy.sum(tpoints, axis=0) / len(tpoints)
        comv = numpy.ones_like ( tpoints ) * com
        points = tpoints - comv

        #print points

        ppoints = points * self.U

        ex_max = numpy.asarray( numpy.max ( ppoints, 0 ) )[0]
        ex_min = numpy.asarray( numpy.min ( ppoints, 0 ) )[0]

        print "Max:"
        print ex_max

        print "Min:"
        print ex_min

        self.dir = chimera.Vector ( self.U[0,2], self.U[1,2], self.U[2,2] )
        print "DIR:  %.2f %.2f %.2f" % (self.dir[0], self.dir[1], self.dir[2])

        zmax = ex_max[2]
        zmin = -ex_min[2]
        if zmax > zmin :
            move_d = (zmax+zmin)/2.0 - zmin
            print " move z + %.3f" % move_d
            self.COM = self.COM + self.dir * move_d
        elif zmax < zmin :
            move_d = (zmax+zmin)/2.0 - zmax
            print " move z - %.3f" % move_d
            self.COM = self.COM - self.dir * move_d


        ymax = ex_max[1]
        ymin = ex_min[1]
        self.height = (ymax - ymin) * 0.6

        xmax = ex_max[0]
        xmin = ex_min[0]


        #self.Extents = numpy.asarray ( numpy.max ( numpy.abs ( ppoints ), 0 ) )[0]
        self.length = (zmax+zmin)
        self.numRes = int ( numpy.round(self.length / 3.5) )

        print " - length %.2f, height %.2f ~%d res" % (self.length, self.height, self.numRes)

        self.rid = reg.rid
        self.regs = sregs

        return self.MakeMod ()



    def MakeMod ( self ) :

        self.length = self.numRes * 3.5

        segMod = GetSegMod ()

        import axes; reload (axes)

        print " - making sheet mod", self.id

        if hasattr ( self, 'sps' ) :
            RemoveSurfPieces ( surfMod, self.sps )

        if not hasattr ( self, "switches" ) :
            self.switches = []
            for i in range ( self.numStrands ) :
                self.switches.append ( i % 2 == 1 )



        self.sps = []
        gap_height = 2.0
        strand_height = 2.3
        self.height = self.numStrands * strand_height + (self.numStrands-1) * gap_height

        self.strands = []

        for strand_i in range ( self.numStrands ) :

            #str_h = ( self.height - (self.numStrands-1) * gap_h ) / self.numStrands
            str_y = -self.height/2.0 + float(strand_i) * (strand_height + gap_height)

            smod = Strand ()
            smod.switch = self.switches[strand_i] # strand_i % 2 == 1
            smod.FromSheet ( self.U, self.COM, str_y, strand_height, self )

            self.strands.append ( smod )

            print "%d - %.3f, %.2f per sheet" % (strand_i, str_y, strand_height)

            if 0 :
                import _surface
                newMod = _surface.SurfaceModel()
                clr = ( segmod_dialog().sheetBaseClr[0], segmod_dialog().sheetBaseClr[1], segmod_dialog().sheetBaseClr[2], 1 )
                axes.BoxMesh ( self.width, str_h, self.length, clr, None, newMod )

                U = self.U
                COM = self.COM
                #if U != None :

                R = numpy.array([
                                 [  U[0,0], U[0,1], U[0,2], 0.0    ],
                                 [  U[1,0], U[1,1], U[1,2], 0.0    ],
                                 [  U[2,0], U[2,1], U[2,2], 0.0    ]  ] )

                T = numpy.array([
                                 [  1.0, 0.0, 0.0, COM[0]   ],
                                 [  0.0, 1.0, 0.0, COM[1]   ],
                                 [  0.0, 0.0, 1.0, COM[2]   ]  ] )

                Ti = numpy.array([
                                  [  1.0, 0.0, 0.0, 0  ],
                                  [  0.0, 1.0, 0.0, str_p   ],
                                  [  0.0, 0.0, 1.0, 0   ]  ] )

                import Matrix
                M = Matrix.multiply_matrices ( R, Ti )
                M = Matrix.multiply_matrices ( T, M )

                sp = newMod.surfacePieces[-1]
                v, t = numpy.copy (sp.geometry[0]), numpy.copy(sp.geometry[1])
                #segMod.removePiece ( sp )

                _contour.affine_transform_vertices( v, M )
                sp = segMod.addPiece ( v, t, (1,.2,.2,1) )
                self.sps.append ( sp )
                sp.displayStyle = sp.Mesh
                sp.lineThickness = 2.0


        if not self in segMod.mods :
            segMod.mods.append ( self )
            umsg ( "Sheet, length %.2f, ~%.0f res, now %d parts" % ( self.length, self.numRes, len(segMod.mods) ) )
        else :
            umsg ( "Updated sheet --??-- length %.2f, ~%.0f res, %d parts" % ( self.length, self.numRes, len(segMod.mods) ) )


        return segMod







def FitPointsToData ( fpoints, data ) :

    #print " - fitting %s to %s" % (m2.name, dmap.name)

    #fpoints, fpoint_weights = fit_points ( m2, False )
    #_contour.affine_transform_vertices ( fpoints, Matrix.xform_matrix(m2.openState.xform) )

    #m1M = Matrix.multiply_matrices( data.xyz_to_ijk_transform, Matrix.xform_matrix(dmap.openState.xform.inverse()) )

    xf = chimera.Xform.identity()
    m1M = data.xyz_to_ijk_transform

    for i in range ( 1 ) :

        fpoint_weights = numpy.ones ( len(fpoints), numpy.float32 )
        move_tf, stats = FitMap.locate_maximum(fpoints, fpoint_weights,
                                        data.matrix(), m1M,
                                        max_steps = 1000,
                                        ijk_step_size_min = 0.01,
                                        ijk_step_size_max = 0.5,
                                        optimize_translation = True,
                                        optimize_rotation = True,
                                        metric = 'sum product',
                                        request_stop_cb = None)

        moveXf = Matrix.chimera_xform ( move_tf )
        xf.premultiply ( moveXf )

        m1M = Matrix.multiply_matrices( move_tf, m1M )

    return xf





def SetBBAts ( ress ) :
    for r in ress :
        try :
            r.CA = r.atomsMap["CA"][0]
            r.C = r.atomsMap["C"][0]
            r.N = r.atomsMap["N"][0]
            r.O = r.atomsMap["O"][0]
        except :
            #print " - res %d.%s has not O" % (r.id.position, r.type)
            pass
        try :
            r.CB = r.atomsMap["CB"][0]
        except :
            #print " - res %d.%s has not O" % (r.id.position, r.type)
            r.CB = None




def vnorm ( a, b, c ) :
    v1 = a - b
    v2 = c - b
    n = chimera.cross ( v1, v2 )
    n.normalize()
    return n


# using atoms
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

# using coords
def diha_ ( a1, a2, a3, a4 ) :
    n1 = vnorm ( a1, a2, a3 )
    n2 = vnorm ( a4, a3, a2 )
    return numpy.arccos ( n2 * n1 ) * 180.0 / numpy.pi - 180


def angle ( a1, a2, a3 ) :
    n1 = a1.coord() - a2.coord()
    n2 = a3.coord() - a2.coord()
    return numpy.arccos ( (n2/n1.length) * (n1/n2.length) )  * 180.0 / numpy.pi



def ptOnPlane ( O, N, p ) :

    d = (p - O) * N
    r = p - d * N
    #d2 = (r - O) * N
    #print d2
    return r





def BuildModLoop ( mol, startRi, endRi, seq, chainId ) :


    print " - build loop in %s, cid %s, start %d" % (mol.name, chainId, startRi)

    rids = {}
    for r in mol.residues :
        if r.id.chainId == chainId :
            rids[r.id.position] = r


    seq0 = seq
    startRes, endRes = None, None
    phisPsis = []
    startPhi, startPsi = random.random()*360.0, random.random()*360.0
    endPhi, endPsi = random.random()*360.0, random.random()*360.0

    if startRi-1 in rids :
        startRes = rids[ startRi-1 ]
        #print " - start res: %d - %s" % ( startRes.id.position, startRes.type )
        SetBBAts ( [startRes] )

        startPsi = diha ( startRes.N, startRes.CA, startRes.C, startRes.O ) + 180.0 # adding 180 because of measurement to O atom
        print " - start res: %d - %s - psi: %.3f" % ( startRes.id.position, startRes.type, startPsi )

        seq = protein3to1[startRes.type] + seq
        phisPsis.append ( [startPhi, startPsi] )
        startRi -= 1

    for i in range ( len(seq0) ) :
        phisPsis.append ( [random.random()*360.0,random.random()*360.0] )

    if endRi+1 in rids :
        endRes = rids[ endRi+1 ]
        #print "   - end res: %d - %s" % ( endRes.id.position, endRes.type )
        endPsi = diha ( endRes.N, endRes.CA, endRes.C, endRes.O ) + 180.0 # adding 180 because of measurement to O atom
        print " - end res: %d - %s - psi: %.3f" % ( endRes.id.position, endRes.type, endPsi )

        seq = seq + protein3to1[endRes.type]
        phisPsis.append ( [endPhi, endPsi] )
        endRi += 1

    print seq
    print phisPsis



    for ri in range ( startRi, endRi+1 ) :
        if ri in rids :
            res = rids[ri]
            print " - deleting residue %d - %s" % (ri, res.type)
            mol.deleteResidue ( res )


    nmol = chimera.Molecule()
    nmol.name = "%d-%s-%d" % (startRi, seq, endRi)
    chimera.openModels.add ( [nmol] )

    ress = AddRess ( nmol, chainId, seq, startRi, phisPsis, startRes, endRes, False )

    bMap = {}
    for at in nmol.atoms :
        at.display = True
        at.drawMode = at.EndCap
        #at.drawMode = at.Dot
    for b in nmol.bonds :
        b.drawMode = bd.Smart
        for at in b.atoms :
            if at in bMap : bMap[at].append ( b )
            else : bMap[at] = [b]

    for ri, res in enumerate ( ress ) :
        res.isHelix, res.isStrand, res.isSheet, res.ribbonDisplay = False, False, False, True


    #chimera.openModels.remove ( [tomol], destroying=True )
    #RefineLoop3 ( lmod, 1, True )




def AddRess_N ( seq, connectRes, helix=False ) :

    toMol, toChain = connectRes.molecule, connectRes.id.chainId
    rmap = ResMap ( toMol, toChain )
    fromRi = connectRes.id.position-len(seq)
    toRi = connectRes.id.position-1
    for ri in range ( fromRi, toRi+1 ) :
        if ri in rmap :
            umsg ( "Molecule already has residue at position %d" % ri  )
            return

    R = connectRes
    atN, atC, atCA, atO = R.atomsMap["N"][0], R.atomsMap["C"][0], R.atomsMap["CA"][0], R.atomsMap["O"][0]
    psi = 180.0 + diha ( atN, atCA, atC, atO )

    from random import random

    phisPsis = []
    for i in range ( len(seq) ) :
        if helix :
            phisPsis.append ( [-57,-47] )
        else :
            phisPsis.append ( [random()*360.0,random()*360.0] )

    phisPsis.append ( [0,psi] )
    seq = seq + "A"

    nmol = chimera.Molecule()
    nmol.name = "_add_%s" % seq
    ress = AddRess ( nmol, 'A', seq, 1, phisPsis, None, None, doRota=False )

    atoms = []
    #apoints = numpy.zeros ( [len(self.ress), 3], numpy.float32 )

    for ri, res in enumerate ( ress ) :
        res.isHelix, res.isStrand, res.isSheet, res.ribbonDisplay = helix, False, False, False
        res.ribbonColor = chimera.MaterialColor ( 0.7, 0.3, 0.3, 1.0 )
        for at in res.atoms :
            at.display = True
            at.drawMode = at.EndCap
            at.color = atomColors[at.element.name if at.element.name in atomColors else " "]

    for b in nmol.bonds :
        b.display = b.Smart
        b.drawMode = b.Stick

    print " - added %d residues, %d atoms, %d bonds" % (len(ress), len(nmol.atoms), len(nmol.bonds))

    res_ = ress[-1]
    atN_, atC_, atCA_, atO_ = res_.atomsMap["N"][0], res_.atomsMap["C"][0], res_.atomsMap["CA"][0], res_.atomsMap["O"][0]

    # match new res CA to apoints :)
    refAts = [atN, atC, atCA, atO]
    newAts = [atN_, atC_, atCA_, atO_]
    xf, min_rmsd = chimera.match.matchAtoms ( refAts, newAts )
    print " - align rmsd: ", min_rmsd
    for res in ress :
        for at in res.atoms :
            at.setCoord ( xf.apply (at.coord()) )

    #chimera.openModels.add ( [nmol] )


    #startI = res.id.position - len(seq)
    #print " - adding ress at %d" % startI

    lastRes = None
    if fromRi-1 in rmap :
        lastRes = rmap[fromRi-1]
        print " - connecting to previous res: %d" % (fromRi-1)

    atRi = fromRi
    for res in ress[0:-1] :
        nres = AddRes ( toMol, toChain, atRi, res )
        if lastRes != None :
            nb = toMol.newBond ( lastRes.atomsMap["C"][0], nres.atomsMap["N"][0] )
            nb.display = nb.Smart
            nb.drawMode = nb.Stick
        lastRes = nres
        atRi += 1

    nb = toMol.newBond ( lastRes.atomsMap["C"][0], connectRes.atomsMap["N"][0] )
    nb.display = nb.Smart
    nb.drawMode = nb.Stick




def AddRess_N0 ( seq, connectRes, helix=False ) :

    res = connectRes
    atN, atC, atCA, atO = res.atomsMap["N"][0], res.atomsMap["C"][0], res.atomsMap["CA"][0], res.atomsMap["O"][0]
    psi = 180.0 + diha ( atN, atCA, atC, atO )

    from random import random

    phisPsis = []
    for i in range ( len(seq) ) :
        if helix :
            phisPsis.append ( [-57,-47] )
        else :
            phisPsis.append ( [random()*360.0,random()*360.0] )

    phisPsis.append ( [0,psi] )
    seq = seq + "A"

    nmol = chimera.Molecule()
    nmol.name = "_add_%s" % seq
    ress = AddRess ( nmol, 'A', seq, 1, phisPsis, None, None, doRota=False )

    atoms = []
    #apoints = numpy.zeros ( [len(self.ress), 3], numpy.float32 )

    for ri, res in enumerate ( ress ) :
        res.isHelix, res.isStrand, res.isSheet, res.ribbonDisplay = helix, False, False, False
        res.ribbonColor = chimera.MaterialColor ( 0.7, 0.3, 0.3, 1.0 )
        for at in res.atoms :
            at.display = True
            at.drawMode = at.EndCap
            at.color = atomColors[at.element.name if at.element.name in atomColors else " "]

    for b in nmol.bonds :
        b.display = b.Smart
        b.drawMode = b.Stick

    print " - added %d residues, %d atoms, %d bonds" % (len(ress), len(nmol.atoms), len(nmol.bonds))

    res_ = ress[-1]
    atN_, atC_, atCA_, atO_ = res_.atomsMap["N"][0], res_.atomsMap["C"][0], res_.atomsMap["CA"][0], res_.atomsMap["O"][0]

    # match new res CA to apoints :)
    refAts = [atN, atC, atCA, atO]
    newAts = [atN_, atC_, atCA_, atO_]
    xf, min_rmsd = chimera.match.matchAtoms ( refAts, newAts )
    print " - align rmsd: ", min_rmsd
    for res in ress :
        for at in res.atoms :
            at.setCoord ( xf.apply (at.coord()) )

    #chimera.openModels.add ( [nmol] )


    #startI = res.id.position - len(seq)
    #print " - adding ress at %d" % startI

    rmap = {}
    maxRi = 0
    for r in connectRes.molecule.residues :
        if r.id.chainId == res.id.chainId :
            maxRi = max ( maxRi,  r.id.position )
            rmap[r.id.position] = r

    cress = []
    ri = connectRes.id.position
    while ri in rmap :
        cress.append ( rmap[ri] )
        ri += 1

    atRi = maxRi + 2
    print " - new start %d, + %d -> %d" % ( maxRi, cress[0].id.position, cress[-1].id.position )

    toMol = connectRes.molecule
    lastRes = None
    for res in ress[0:-1] :
        nres = AddRes ( toMol, connectRes.id.chainId, atRi, res )
        if lastRes != None :
            nb = toMol.newBond ( lastRes.atomsMap["C"][0], nres.atomsMap["N"][0] )
            nb.display = nb.Smart
            nb.drawMode = nb.Stick
        lastRes = nres
        atRi += 1

    delAts, delBonds = {}, {}
    for res in cress :
        nres = AddRes ( toMol, connectRes.id.chainId, atRi, res )
        for at in res.atoms :
            delAts[at] = 1
            for b in at.bondsMap.values() :
                delBonds[b] = 1
        if lastRes != None :
            nb = toMol.newBond ( lastRes.atomsMap["C"][0], nres.atomsMap["N"][0] )
            nb.display, nb.drawMode = nb.Smart, nb.Stick
        lastRes = nres
        atRi += 1

    if 1 :
        for b in delBonds.keys() :
            toMol.deleteBond ( b )
        for at in delAts.keys() :
            toMol.deleteAtom ( at )
        for res in cress :
            toMol.deleteResidue ( res )



def ConnectRess ( atN, atC ) :

    rmapN = {}
    for r in atN.molecule.residues :
        if r.id.chainId == atN.residue.id.chainId :
            rmapN[r.id.position] = 1

    rmapC = {}
    for r in atC.molecule.residues :
        if r.id.chainId == atC.residue.id.chainId :
            rmapC[r.id.position] = 1

    if atC.residue.id.position + 1 in rmapC :
        # N -> C
        if atN.residue.id.position + 1 in rmapN :
            print " - residue with N already has +1 connection"
            return

        print " N -> C"

    elif atN.residue.id.position + 1 in rmapN :
        # C -> N
        if atN.residue.id.position - 1 in rmapN :
            print " - residue with N already has -1 connection"
            return
        if atC.residue.id.position + 1 in rmapC :
            print " - residue with C already has +1 connection"
            return

        Connect_C_N ( atN, atC )


def ResMap ( mol, chainId ) :
    rmap = {}
    for res in mol.residues :
        if res.id.chainId == chainId :
            rmap[res.id.position] = res
    return rmap




def Connect_C_N ( atN, atC ) :

    print " C -> N "
    print "   %d -> %d " % (atC.residue.id.position, atN.residue.id.position)

    molC = atC.molecule
    chainC = atC.residue.id.chainId

    chainN = atN.residue.id.chainId
    molN = atN.molecule

    rmapN = ResMap ( molN, chainN )
    riN = atN.residue.id.position
    ressN = []
    rmapN_ = {}
    while riN in rmapN :
        res = rmapN[riN]
        ressN.append ( res )
        rmapN_[res] = 1
        riN += 1

    print " -N fragment has %d residues, chain %s" % (len(ressN), chainN)

    nmol = chimera.Molecule()
    nmol.name = molC.name + "_"

    # copy all other chains...

    # copy chainC, up to C atom residue
    rmapC = ResMap ( molC, chainC )
    riMin, riMax = min(rmapC.keys()), max(rmapC.keys())
    toRi = 5
    lastRes2, lastRes = None, None
    for ri in range (riMin, atC.residue.id.position+1) :
        if ri in rmapC :

            resToAdd = rmapC[ri]
            if resToAdd in rmapN_ :
                #print " -1- (skipping res %d in N-fragment...)" % (resToAdd.id.position)
                continue

            if lastRes2 != None and resToAdd.id.position > lastRes2.id.position + 1 :
                toRi += 7
            lastRes2 = resToAdd

            #print " -1- adding res %d.%s -> %d.%s" % (resToAdd.id.position, resToAdd.id.chainId, toRi, chainC)

            nres = AddRes ( nmol, chainC, toRi, resToAdd )
            if lastRes != None and nres.id.position == lastRes.id.position + 1 :
                nb = nmol.newBond ( lastRes.atomsMap["C"][0], nres.atomsMap["N"][0] )
                nb.display, nb.drawMode = nb.Smart, nb.Stick
            lastRes = nres
            toRi += 1


    # copy fragment starting at N
    for res in ressN :

        #print " -f- adding res %d.%s -> %d.%s" % (res.id.position, res.id.chainId, toRi, chainC)
        nres = AddRes ( nmol, chainC, toRi, res )
        if lastRes != None and nres.id.position == lastRes.id.position + 1 :
            nb = nmol.newBond ( lastRes.atomsMap["C"][0], nres.atomsMap["N"][0] )
            nb.display, nb.drawMode = nb.Smart, nb.Stick
        lastRes = nres
        toRi += 1

    # copy rest of chainC
    for ri in range ( atC.residue.id.position+1, riMax+1 ) :
        if ri in rmapC :

            resToAdd = rmapC[ri]
            if resToAdd in rmapN_ :
                #print " -2- (skipping res %d in N-fragment...)" % (resToAdd.id.position)
                continue

            if lastRes2 != None and resToAdd.id.position > lastRes2.id.position + 1 :
                toRi += 7
            lastRes2 = resToAdd

            #print " -2- adding res %d.%s -> %d.%s" % (resToAdd.id.position, resToAdd.id.chainId, toRi, chainC)

            nres = AddRes ( nmol, chainC, toRi, resToAdd )
            if lastRes != None and nres.id.position == lastRes.id.position + 1 :
                nb = nmol.newBond ( lastRes.atomsMap["C"][0], nres.atomsMap["N"][0] )
                nb.display, nb.drawMode = nb.Smart, nb.Stick
            lastRes = nres
            toRi += 1

    #chimera.openModels.close ( [molC] )
    chimera.openModels.add ( [nmol] )




def Connect_C_N_0_del ( atN, atC ) :

        print " C -> N "
        print "   %d -> %d " % (atC.residue.id.position, atN.residue.id.position)

        sameChainMol = True
        if 0 and atN.molecule != atC.molecule :
            print " C and N are in different molecules..."
            sameChainMol = False
            return

        if 0 and atN.residue.id.chainId != atC.residue.id.chainId :
            print " C and N are in different chains..."
            sameChainMol = False
            return


        chainN = atN.residue.id.chainId
        posN = atN.residue.id.position
        # copy fragment starting at N residue to new molecule
        Nmol = chimera.Molecule ()
        Nmol.name = " _ N - copy _ "
        rmap = ResMap ( atN.molecule, chainN )
        atRi = atN.residue.id.position
        while atRi in rmap :
            res = rmap[atRi]
            nres = AddRes ( Nmol, res.id.chainId, res.id.position, res )
            res.molecule.deleteResidue ( res )
            atRi += 1
        chimera.openModels.add ( [Nmol] )

        rmapN = ResMap ( Nmol, chainN )
        atN = rmapN[posN].atomsMap["N"][0]

        # copy (and delete) atC molecule - only residues after atC - to Cmol
        Cmol = chimera.Molecule()
        Cmol.name = " _ C - copy _ "
        rmap = {}
        for res in atC.molecule.residues [:] :
            if res.id.chainId == atC.residue.id.chainId :
                if res.id.position > atC.residue.id.position :
                    #if sameChainMol and res.id.position in rmapN :
                    #    continue
                    nres = AddRes ( Cmol, res.id.chainId, res.id.position, res )
                    rmap[nres.id.position] = nres
                    atC.molecule.deleteResidue (res)

        chimera.openModels.add ( [Cmol] )

        # add (and delete) fragment starting with atN
        atRi = atN.residue.id.position
        toRi = atC.residue.id.position+1
        lastRes = atC.residue
        while 1 :
            if atRi not in rmapN :
                print " . added residues up to %d -> %d" % (atRi, toRi)
                break

            resToAdd = rmapN[atRi]
            print " -1- adding res %d.%s -> %d.%s" % (resToAdd.id.position, resToAdd.id.chainId, toRi, atC.residue.id.chainId)

            nres = AddRes ( atC.molecule, atC.residue.id.chainId, toRi, resToAdd )
            nb = atC.molecule.newBond ( lastRes.atomsMap["C"][0], nres.atomsMap["N"][0] )
            nb.display, nb.drawMode = nb.Smart, nb.Stick
            lastRes = nres
            #resToAdd.molecule.deleteResidue (resToAdd)

            atRi += 1
            toRi += 1

        # add remaining residues/fragments in Cmol
        if len(Cmol.residues) > 0 :
            rmap = ResMap ( Cmol, atC.residue.id.chainId )
            firstRi = min ( [res.id.position for res in Cmol.residues] )
            lastRi = max ( [res.id.position for res in Cmol.residues] )

            toRi += 5
            print " . adding remaining %d residues (%d -> %d), starting at %d" % ( len(Cmol.residues), firstRi, lastRi, toRi )

            lastRes, lastResN = None, None
            ri = firstRi
            while ri <= lastRi :
                if ri in rmap :
                    resToAdd = rmap[ri]

                    if lastResN != None and resToAdd.id.position > lastResN.id.position+1 :
                        toRi += 5
                    lastResN = resToAdd

                    print " -2- adding res %d.%s -> %d.%s" % (resToAdd.id.position, resToAdd.id.chainId, toRi, atC.residue.id.chainId)

                    nres = AddRes ( atC.molecule, atC.residue.id.chainId, toRi, resToAdd )

                    if lastRes != None and nres.id.position == lastRes.id.position+1 :
                        nb = atC.molecule.newBond ( lastRes.atomsMap["C"][0], nres.atomsMap["N"][0] )
                        nb.display, nb.drawMode = nb.Smart, nb.Stick
                    lastRes = nres

                    toRi += 1

                ri += 1

        nmol = chimera.Molecule ()
        nmol.name = " _ copy _ "
        rmap = ResMap ( atC.molecule, atC.residue.id.chainId )
        firstRi = min ( rmap.keys() )
        lastRi = max ( rmap.keys() )
        lastRes = None
        for ri in range ( firstRi, lastRi+1 ) :
            if ri in rmap :
                res = rmap[ri]
                nres = AddRes ( nmol, res.id.chainId, res.id.position, res )
                if lastRes != None and nres.id.position == lastRes.id.position+1 :
                    nb = atC.molecule.newBond ( lastRes.atomsMap["C"][0], nres.atomsMap["N"][0] )
                    nb.display, nb.drawMode = nb.Smart, nb.Stick
                lastRes = nres
        chimera.openModels.add ( [nmol] )



def SetResI ( atRes, newResI ) :


    addN = newResI - atRes.id.position
    print " - setting res %d.%s to %d +%d" % (atRes.id.position, atRes.id.chainId, newResI, addN)

    chain = atRes.id.chainId
    mol = atRes.molecule

    nmol = chimera.Molecule()
    nmol.name = mol.name + "-"
    #chimera.openModels.close ( [molC] )
    chimera.openModels.add ( [nmol] )

    # add all other chains ?
    aMap = {}
    for res in mol.residues :
        if res.id.chainId != chain :
            nres = nmol.newResidue (res.type, chimera.MolResId(res.id.chainId, res.id.position))
            for at in res.atoms :
                nat = nmol.newAtom (at.name, chimera.Element(at.element.number))
                aMap[at] = nat
                nres.addAtom( nat )
                nat.setCoord ( at.coord() )
                nat.altLoc, nat.occupancy, nat.bfactor = at.altLoc, at.occupancy, at.bfactor
                nat.color, nat.display, nat.drawMode = at.color, at.display, nat.EndCap
            nres.isHelix, nres.isHet, nres.isSheet, nres.isStrand = res.isHelix, res.isHet, res.isSheet, res.isStrand
            nres.ribbonDisplay, nres.ribbonColor = res.ribbonDisplay, res.ribbonColor
            if nres.type in nucleic3to1 :
                nres.fillDisplay = False

    for bond in mol.bonds :
        a1, a2 = bond.atoms
        if a1 in aMap and a2 in aMap :
            nb = nmol.newBond ( aMap[bond.atoms[0]], aMap[bond.atoms[1]] )
            nb.display, nb.drawMode = nb.Smart, nb.Stick

    # copy ress...
    rmap = ResMap ( mol, chain )
    riMin, riMax = min(rmap.keys()), max(rmap.keys())
    print " - min %d max %d" % (riMin, riMax)
    toRi = riMin
    lastRes = None
    for ri in range (riMin, riMax+addN+1) :
        if not ri in rmap :
            continue

        toRi = ri
        if ri >= atRes.id.position :
            toRi = ri + addN

        resToAdd = rmap[ri]
        #print " %d -> %d" % (ri, toRi)

        nres = AddRes ( nmol, chain, toRi, resToAdd )
        isProt = nres.type in protein3to1
        if lastRes != None and not lastRes.type in protein3to1 :
            isProt = False
        if isProt and lastRes != None and nres.id.position == lastRes.id.position + 1 :
            nb = nmol.newBond ( lastRes.atomsMap["C"][0], nres.atomsMap["N"][0] )
            nb.display, nb.drawMode = nb.Smart, nb.Stick
        lastRes = nres


    from NucleicAcids.cmd import sidechain
    sidechain("atoms", sel="#%d" % nmol.id)
    for r in nmol.residues :
        if r.type in nucleic3to1 :
            r.fillDisplay = False
            for at in r.atoms :
                at.display = False
    mol.display = False


def AddRes ( nmol, chainId, resI, res ) :

    aMap = {}
    #nres = nmol.newResidue (res.type, chimera.MolResId(res.id.chainId, res.id.position))
    nres = nmol.newResidue (res.type, chimera.MolResId(chainId, resI))
    # print "New res: %s %d" % (nres.id.chainId, nres.id.position)
    for at in res.atoms :
        nat = nmol.newAtom (at.name, chimera.Element(at.element.number))
        aMap[at] = nat
        nres.addAtom( nat )
        nat.setCoord ( at.coord() )
        nat.altLoc = at.altLoc
        nat.occupancy = at.occupancy
        nat.bfactor = at.bfactor
        nat.color = at.color
        #nat.radius = at.radius
        nat.display = at.display
        nat.drawMode = nat.EndCap

    nres.isHelix, nres.isHet, nres.isSheet, nres.isStrand = res.isHelix, res.isHet, res.isSheet, res.isStrand
    nres.ribbonDisplay = res.ribbonDisplay
    #nres.ribbonDrawMode = 2
    nres.ribbonColor = res.ribbonColor

    for bond in res.molecule.bonds :
        a1, a2 = bond.atoms
        if a1 in aMap and a2 in aMap :
            nb = nmol.newBond ( aMap[bond.atoms[0]], aMap[bond.atoms[1]] )
            nb.display = nb.Smart
            nb.drawMode = nb.Stick

    return nres



def AddRess ( mol, chainID, sequence, startResI, phiPsis, startRes=None, endRes=None, doRota=True ) :

    from chimera.resCode import protein1to3
    from chimera import Point, Element

    rotlib = None
    prev = [None] * 3
    pos = startResI
    from Midas.addAA import DIST_N_C, DIST_CA_N, DIST_C_CA, DIST_C_O
    from chimera.molEdit import findPt, addAtom, addDihedralAtom
    serialNumber = None
    residues = []
    prevPsi = 0
    for c, phiPsi in zip(sequence, phiPsis):
        phi, psi = phiPsi

        #while model.findResidue(chimera.MolResId(chainID, pos)):
        #	pos += 1

        r = mol.newResidue(protein1to3[c], chainID, pos, ' ')
        residues.append(r)
        pos += 1
        #print " - added %s at %d" % (r.type, pos)

        for backbone, dist, angle, dihed in (
				('N', DIST_N_C, 116.6, prevPsi),
				('CA', DIST_CA_N, 121.9, 180.0),
				('C', DIST_C_CA, 110.1, phi)):
            if prev[0] == None : pt = Point(0.0, 0.0, 0.0)
            elif prev[1] == None : pt = Point(dist, 0.0, 0.0)
            elif prev[2] == None : pt = findPt(prev[0].coord(), prev[1].coord(), Point(0.0, 1.0, 0.0), dist, angle, 0.0)
            else :                 pt = findPt(prev[0].coord(), prev[1].coord(), prev[2].coord(), dist, angle, dihed )

            if pos-1 == startResI and startRes != None and backbone == "N" :
                bat = startRes.atomsMap["C"][0]
                a = addAtom(backbone, Element(backbone[0]), r, pt, serialNumber=serialNumber, bondedTo=bat)
                #print " - C-bond to res %d %s" % (startRes.id.position, startRes.type)
            elif pos == startResI + len(phiPsis) and endRes != None and backbone == "C" :
                bat = endRes.atomsMap["N"][0]
                a = addAtom(backbone, Element(backbone[0]), r, pt, serialNumber=serialNumber, bondedTo=bat)
                # for some reason bonding C -> next N with bondedTo param removes the C-CA bond - add it back :)
                nb = mol.newBond ( r.atomsMap["CA"][0], r.atomsMap["C"][0] )
                #print " - N-bond to res %d %s" % (endRes.id.position, endRes.type)
            else :
                a = addAtom(backbone, Element(backbone[0]), r, pt, serialNumber=serialNumber, bondedTo=prev[0])
            #serialNumber = a.serialNumber + 1
            prev = [a] + prev[:2]

        o = addDihedralAtom("O", Element("O"), prev[0], prev[1], prev[2], DIST_C_O, 120.4, 180.0 + psi, bonded=True)
        prevPsi = psi

	# C terminus O/OXT at different angle than mainchain O
	#model.deleteAtom(o)
	#addDihedralAtom("O", Element("O"), prev[0], prev[1],
	#		prev[2], DIST_C_O, 117.0, 180.0 + psi, bonded=True)
	#addDihedralAtom("OXT", Element("O"), prev[0], prev[1], prev[2],
	#				DIST_C_O, 117.0, psi, bonded=True)

    if doRota :
        from Rotamers import useBestRotamers
        # have to process one by one, otherwise side-chain clashes will occur
        for r in residues:
            if 1 and r.type == "PRO" :
                useBestRotamers("same", [r], criteria="cp", log=False, lib="Richardson.mode")
                try :
                    print "pro"
                except :
                    print " - failed to adjust rotamer for r %d %s" % (r.id.position, r.type)
            else :
                useBestRotamers("same", [r], criteria="cp", log=False, lib="Dunbrack")
                try :
                    print "",
                    #useBestRotamers("same", [r], criteria="cp", log=False, lib="Dunbrack")
                except :
                    print " - failed to adjust rotamer for r %d %s" % (r.id.position, r.type)

    else :
    	from SwapRes import swap, SwapResError
        from chimera.resCode import protein1to3
        for i, r in enumerate ( residues ) :
			swap(r, protein1to3[sequence[i]], preserve=False, bfactor=False)

    return residues





def Connect ( mol, molc, chainId ) :

    print ""

    rids = {}
    for r in mol.residues :
        if r.id.chainId == chainId :
            rids[r.id.position] = r

    for r in m2.residues :
        if not r.id.position in rids :
            print " - %d %s %s" % (r.id.position, r.type, r.id.chainId)
            chimera.selection.addCurrent ( r )




def ResRotaD ( r, asType, agrid, dmap ) :


    try :
        r.CB = r.atomsMap["CB"][0]
    except :
        r.CB = None

    if r.CB == None :
        return

    if r.type == "ALA" or r.type == "GLY" :
        return

    from Rotamers import getRotamers
    bbdep, rmols = getRotamers ( r, log=False )

    apos = _multiscale.get_atom_coordinates(r.atoms, transformed = False)
    #print apos

    min_ri, min_ri_d = -1, -1
    min_score, min_score_d = 1e99, 1e99

    for ri, rmol in enumerate ( rmols ) :

        rotres = rmol.residues[0]

        #print rotres.atomsMap

        #to_ats = [r.N, r.CA, r.CB]
        to_ats = [ r.atomsMap['N'][0],r.atomsMap['CA'][0],r.atomsMap['CB'][0] ]
        rot_ats = [ rotres.atomsMap['N'][0],rotres.atomsMap['CA'][0],rotres.atomsMap['CB'][0] ]

        xf, rmsd = chimera.match.matchAtoms ( to_ats, rot_ats )

        score = 0.0
        useRot = 1
        for ai, rat in enumerate ( rotres.atoms ) :

            if rat.name == "C" or rat.name == "N" or rat.name == "CA" or rat.name == "O" :
                continue

            apos[ai] = xf.apply (rat.coord()).data()

            #at = r.atomsMap[rat.name][0]
            #at.setCoord ( xf.apply (rat.coord()) )
            trPt = xf.apply (rat.coord())

            atsNear = atTree.searchTree ( trPt.data(), 5.0 )
            for at in atsNear :
                if at.residue.id.position == r.id.position :
                    continue

                v = trPt - at.coord()
                if v.length < 0.001 :
                    score += 1000.0
                    useRot = 0
                elif v.length < 1.25 :
                    score += 1.0 / v.length
                    useRot = 0
                    #print "%s.%d.%s -- %s.%d.%s@%s " % (r.type,r.id.position,r.id.chainId, at.residue.type,at.residue.id.position,at.residue.id.chainId,at.name)
                elif v.length < 5 :
                    score += 1.0 / v.length


        if len(rotres.atoms) != len(apos) :
            print " %d mol ats -- %d rot ats" % (len(apos), len(rotres.atoms)),

        dvals = dmap.interpolated_values ( apos, r.molecule.openState.xform )
        dscore = -numpy.average ( dvals )


        if useRot :
            print " -- %d -- %.4f %.5f" % (ri, score, dscore)
            if dscore < min_score_d :
                min_ri_d = ri
                min_score_d = dscore
        else :
            print " -- %d -- %.4f %.5f - X" % (ri, score, dscore)

        if score < min_score :
            min_ri = ri
            min_score = score


    use_ri = -1
    if min_ri_d >= 0 :
        print " - best rotamer %d, score density %.5f" % (min_ri_d, min_score_d)
        use_ri = min_ri_d
    else :
        print " - best rotamer %d, score clash %.5f" % (min_ri, min_score)
        use_ri = min_ri

    rotres = rmols[use_ri].residues[0]
    #to_ats = [r.N, r.CA, r.CB]
    to_ats = [ r.atomsMap['N'][0],r.atomsMap['CA'][0],r.atomsMap['CB'][0] ]
    rot_ats = [ rotres.atomsMap['N'][0],rotres.atomsMap['CA'][0],rotres.atomsMap['CB'][0] ]
    xf, rmsd = chimera.match.matchAtoms ( to_ats, rot_ats )

    for rat in rotres.atoms :
        at = r.atomsMap[rat.name][0]
        trP = xf.apply (rat.coord())
        at.setCoord ( trP )
        #allAtPos[at.allPtI] = trP.data()

        #lmod.allAtTree = AdaptiveTree ( lmod.allAtPos.tolist(), lmod.allAts, 4.0)





def ResRotaD ( r, atTree, dmap ) :


    try :
        r.CB = r.atomsMap["CB"][0]
    except :
        r.CB = None

    if r.CB == None :
        return

    if r.type == "ALA" or r.type == "GLY" :
        return

    from Rotamers import getRotamers
    bbdep, rmols = getRotamers ( r, log=False )

    apos = _multiscale.get_atom_coordinates(r.atoms, transformed = False)
    #print apos

    min_ri, min_ri_d = -1, -1
    min_score, min_score_d = 1e99, 1e99

    for ri, rmol in enumerate ( rmols ) :

        rotres = rmol.residues[0]

        #print rotres.atomsMap

        #to_ats = [r.N, r.CA, r.CB]
        to_ats = [ r.atomsMap['N'][0],r.atomsMap['CA'][0],r.atomsMap['CB'][0] ]
        rot_ats = [ rotres.atomsMap['N'][0],rotres.atomsMap['CA'][0],rotres.atomsMap['CB'][0] ]

        xf, rmsd = chimera.match.matchAtoms ( to_ats, rot_ats )

        score = 0.0
        useRot = 1
        for ai, rat in enumerate ( rotres.atoms ) :

            if rat.name == "C" or rat.name == "N" or rat.name == "CA" or rat.name == "O" :
                continue

            apos[ai] = xf.apply (rat.coord()).data()

            #at = r.atomsMap[rat.name][0]
            #at.setCoord ( xf.apply (rat.coord()) )
            trPt = xf.apply (rat.coord())

            atsNear = atTree.searchTree ( trPt.data(), 5.0 )
            for at in atsNear :
                if at.residue.id.position == r.id.position :
                    continue

                v = trPt - at.coord()
                if v.length < 0.001 :
                    score += 1000.0
                    useRot = 0
                elif v.length < 1.25 :
                    score += 1.0 / v.length
                    useRot = 0
                    #print "%s.%d.%s -- %s.%d.%s@%s " % (r.type,r.id.position,r.id.chainId, at.residue.type,at.residue.id.position,at.residue.id.chainId,at.name)
                elif v.length < 5 :
                    score += 1.0 / v.length


        if len(rotres.atoms) != len(apos) :
            print " %d mol ats -- %d rot ats" % (len(apos), len(rotres.atoms)),

        dvals = dmap.interpolated_values ( apos, r.molecule.openState.xform )
        dscore = -numpy.average ( dvals )


        if useRot :
            print " -- %d -- %.4f %.5f" % (ri, score, dscore)
            if dscore < min_score_d :
                min_ri_d = ri
                min_score_d = dscore
        else :
            print " -- %d -- %.4f %.5f - X" % (ri, score, dscore)

        if score < min_score :
            min_ri = ri
            min_score = score


    use_ri = -1
    if min_ri_d >= 0 :
        print " - best rotamer %d, score density %.5f" % (min_ri_d, min_score_d)
        use_ri = min_ri_d
    else :
        print " - best rotamer %d, score clash %.5f" % (min_ri, min_score)
        use_ri = min_ri

    rotres = rmols[use_ri].residues[0]
    #to_ats = [r.N, r.CA, r.CB]
    to_ats = [ r.atomsMap['N'][0],r.atomsMap['CA'][0],r.atomsMap['CB'][0] ]
    rot_ats = [ rotres.atomsMap['N'][0],rotres.atomsMap['CA'][0],rotres.atomsMap['CB'][0] ]
    xf, rmsd = chimera.match.matchAtoms ( to_ats, rot_ats )

    for rat in rotres.atoms :
        at = r.atomsMap[rat.name][0]
        trP = xf.apply (rat.coord())
        at.setCoord ( trP )
        #allAtPos[at.allPtI] = trP.data()

        #lmod.allAtTree = AdaptiveTree ( lmod.allAtPos.tolist(), lmod.allAts, 4.0)




def SetDrawMode ( ress, showRibbon = True ) :

    atMap = {}
    atI = 0
    #c1 = (1.0,0.0,0.0,1)
    #c1 = (1.0,0.0,0.0,1)
    for res in ress :
        for at in res.atoms :
            at.drawMode = at.EndCap
            at.display = True # not showRibbon
            at.color = atomColors[at.element.name if at.element.name in atomColors else " "]
            atMap[at] = 1

        res.ribbonDisplay, res.ribbonDrawMode = showRibbon, res.Ribbon_Round
        #f = float(atI) / float(len(ress)-1)
        #res.ribbonColor = chimera.MaterialColor( f*0.8+0.2, 0.02, (1-f)*0.8+0.2, 1.0 );
        #atI+=1

    for bond in ress[0].molecule.bonds :
        if bond.atoms[0] in atMap or bond.atoms[1] in atMap :
            bond.display = bond.Smart
            bond.drawMode = bond.Stick




def DelRes ( ress ) :

    mol = ress[0].molecule
    datoms, dbonds, dres = {}, {}, {}
    for r in ress :
        dres[r] = 1
        for at in r.atoms :
            datoms[at] = 1
    #for b in mol.bonds :
    #    if b.atoms[0] in datoms or b.atoms[1] in datoms :
    #        dbonds[b] = 1
    for at in datoms.keys() :
        mol.deleteAtom ( at )
    #for b in dbonds.keys() :
    #    mol.deleteBond ( b )
    for r in dres.keys() :
        mol.deleteResidue ( r )


def model_is_deleted(m):
    'Test if model has been deleted'
    try:
        m.display
    except:
        return True
    return False


def RemoveSurfPieces ( surfMod, sps ) :

    if model_is_deleted ( surfMod ) :
        return

    for sp in sps :
        try :
            surfMod.removePiece ( sp )
        except :
            pass



def RegsData ( regs ) :

    segMap = regs[0].segmentation.seg_map
    dmap = segMap

    points = regs[0].points().astype ( numpy.float32 )
    for r in regs[1:] :
        npoints = r.points().astype ( numpy.float32 )
        points = numpy.concatenate ( [points, npoints], axis=0 )

    import _contour
    _contour.affine_transform_vertices ( points, segMap.data.ijk_to_xyz_transform )
    _contour.affine_transform_vertices ( points, Matrix.xform_matrix( segMap.openState.xform ) )
    _contour.affine_transform_vertices ( points, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )

    #points1 = numpy.copy ( points )
    #_contour.affine_transform_vertices ( points1, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
    points0 = numpy.copy ( points )
    _contour.affine_transform_vertices ( points, dmap.data.xyz_to_ijk_transform )

    bound = 5
    li,lj,lk = numpy.min ( points, axis=0 ) - (bound, bound, bound)
    hi,hj,hk = numpy.max ( points, axis=0 ) + (bound, bound, bound)

    n1 = hi - li + 1
    n2 = hj - lj + 1
    n3 = hk - lk + 1

    nstep = (dmap.data.step[0], dmap.data.step[1], dmap.data.step[2] )

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

    npoints = VolumeData.grid_indices ( (nn1, nn2, nn3), numpy.single)  # i,j,k indices
    _contour.affine_transform_vertices ( npoints, ndata.ijk_to_xyz_transform )

    dvals = dmap.interpolated_values ( npoints, dmap.openState.xform )
    nmat = dvals.reshape( (nn3,nn2,nn1) )
    #f_mat = fmap.data.full_matrix()
    #f_mask = numpy.where ( f_mat > fmap.surface_levels[0], numpy.ones_like(f_mat), numpy.zeros_like(f_mat) )
    #df_mat = df_mat * f_mask

    ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles )
    mdata = VolumeData.zone_masked_grid_data ( ndata, points0, dmap.data.step[0] * 0.9 )
    return mdata




def RegsToMap ( regs ) :

    segMap = regs[0].segmentation.seg_map
    dmap = segMap

    points = regs[0].points().astype ( numpy.float32 )
    for r in regs[1:] :
        npoints = r.points().astype ( numpy.float32 )
        points = numpy.concatenate ( [points, npoints], axis=0 )

    import _contour
    _contour.affine_transform_vertices ( points, segMap.data.ijk_to_xyz_transform )
    _contour.affine_transform_vertices ( points, Matrix.xform_matrix( segMap.openState.xform ) )
    _contour.affine_transform_vertices ( points, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )

    #points1 = numpy.copy ( points )
    #_contour.affine_transform_vertices ( points1, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
    points0 = numpy.copy ( points )
    _contour.affine_transform_vertices ( points, dmap.data.xyz_to_ijk_transform )

    bound = 5
    li,lj,lk = numpy.min ( points, axis=0 ) - (bound, bound, bound)
    hi,hj,hk = numpy.max ( points, axis=0 ) + (bound, bound, bound)

    n1 = hi - li + 1
    n2 = hj - lj + 1
    n3 = hk - lk + 1

    nstep = (dmap.data.step[0], dmap.data.step[1], dmap.data.step[2] )

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

    npoints = VolumeData.grid_indices ( (nn1, nn2, nn3), numpy.single)  # i,j,k indices
    _contour.affine_transform_vertices ( npoints, ndata.ijk_to_xyz_transform )

    dvals = dmap.interpolated_values ( npoints, dmap.openState.xform )
    nmat = dvals.reshape( (nn3,nn2,nn1) )
    #f_mat = fmap.data.full_matrix()
    #f_mask = numpy.where ( f_mat > fmap.surface_levels[0], numpy.ones_like(f_mat), numpy.zeros_like(f_mat) )
    #df_mat = df_mat * f_mask

    ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles )

    mdata = VolumeData.zone_masked_grid_data ( ndata, points0, dmap.data.step[0] * 0.9 )
    gdata = VolumeData.Array_Grid_Data ( mdata.full_matrix(), nO, nstep, dmap.data.cell_angles, name = "region masked" )
    nv = VolumeViewer.volume.volume_from_grid_data ( gdata )
    nv.openState.xform = dmap.openState.xform

    nv.name = "region masked"
    #dmap.display = False
    nv.region = ( nv.region[0], nv.region[1], [1,1,1] )
    nv.surface_levels[0] = dmap.surface_levels[0]

    if 0 :
        M = dmap.data.full_matrix()
        sdev, avg, thr = numpy.std(M), numpy.average(M), nv.surface_levels[0]

        M = dmap.data.full_matrix()
        lsdev, lavg = numpy.std(nmat), numpy.average(nmat)

        #print "Avg: %.3f, sdev: %.3f, thr: %.4f [%.4f sdev above mean]" % (avg, sdev, thr, (thr-avg)/sdev)
        sigmaGlobal = (thr-avg)/sdev
        sigmaLocal = (thr-lavg)/lsdev
        #umsg ( "Contour level: %.4f, %.2f/%.2f sigma above average global/local" % (thr, sigmaGlobal, sigmaLocal)  )
        #print sigmaGlobal, sigmaLocal


    ro = VolumeViewer.volume.Rendering_Options()
    ro.smoothing_factor = .3
    ro.smoothing_iterations = 2
    ro.surface_smoothing = False
    ro.square_mesh = True
    ro.line_thickness = 2
    nv.update_surface ( False, ro )
    #setro (ro)
    for sp in nv.surfacePieces :
        v, t = sp.geometry
        if len(v) == 8 and len(t) == 12 :
            sp.display = False
        else :
            if 0 :
                sp.color = (.5, .5, .5, 1.0)
                sp.displayStyle = sp.Mesh
            else :
                sp.color = (0.7, 0.7, 0.7, 0.4)


    return nv




# end
