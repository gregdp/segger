

import chimera
import random
import numpy

mrAtoms = None
mrBonds = None
mrLinkBonds = None
mrAngles = None
mrLinkAngles = None
mrTorsions = None

M_PI = numpy.pi
TWOPI = 2.0 * numpy.pi

def GetAtom ( l ) :
    global mrAtoms
    ts = l.split()
    if len (ts) == 17 :
        #rtype, a1, a2, a3, angle, esd, exception, descr = ts

        rtype, atomId, el = ts[0], ts[1], ts[3],
        try :
            charge = float(ts[10])
        except :
            charge = 0

        try :
            x, y, z = float(ts[11]), float(ts[12]), float(ts[13])
        except :
            x, y, z = 0.0, 0.0, 0.0

        if not rtype in mrAtoms :
            mrAtoms[rtype] = []
        mrAtoms[rtype].append ( (rtype, atomId, el, charge, x, y, z) )


def GetBond ( l ) :
    global mrBonds
    ts = l.split()
    if len (ts) == 11 :
        rtype, a1, a2, order, aromatic, stereo, ordinal, dist, esd, exception, descr = ts
        dist = float ( dist )
        esd = float ( esd )

        #if rtype == "ALA" :
        #    print rtype, a1, a2, dist, esd

        descr = descr.replace ( '"', "" )
        if not rtype in mrBonds :
            mrBonds[rtype] = []
        mrBonds[rtype].append ( [a1, a2, dist, esd, descr] )


def GetLinkBond ( l ) :
    global mrLinkBonds
    ts = l.split()

    #print l,

    if len (ts) > 4 :
        a1, a2, dist, esd = ts[0], ts[1], ts[2], ts[3]
        dist = float ( dist )
        esd = float ( esd )

        descr = ""
        for i in range ( 4, len(ts) ) :
            descr = descr + ts[i] + " "
        descr = descr[0:-1].replace( '"', '' )

        #print descr

        mrLinkBonds[descr] = [a1, a2, dist, esd]

    else :
        #print "?"
        pass



def GetAngle ( l ) :
    global mrAngles
    ts = l.split()
    if len (ts) == 8 :
        rtype, a1, a2, a3, angle, esd, exception, descr = ts
        angle = float ( angle ) * numpy.pi / 180.0
        esd = float ( esd ) * numpy.pi / 180.0

        #if rtype == "ALA" :
        #    print rtype, a1, a2, a3, angle, esd

        descr = descr.replace ( '"', "" )
        if not rtype in mrAngles :
            mrAngles[rtype] = []
        mrAngles[rtype].append ( [a1, a2, a3, angle, esd, descr] )


def GetLinkAngle ( l ) :
    global mrLinkAngles
    ts = l.split()
    if len (ts) > 5 :
        a1, a2, a3, angle, esd, descr = ts[0], ts[1], ts[2], ts[3], ts[4], ""
        angle = float ( angle ) * numpy.pi / 180.0
        esd = float ( esd ) * numpy.pi / 180.0

        for i in range ( 5, len(ts) ) :
            descr = descr + ts[i] + " "

        descr = descr[0:-1].replace( '"', '' )

        if not descr in mrLinkAngles :
            mrLinkAngles[descr] = {}

        mrLinkAngles[descr]["%s_%s_%s"%(a1,a2,a3)] = [angle, esd]
        #print " - angle", "%s_%s_%s : %.3f, %.3f"%(a1,a2,a3,angle,esd)
        #mrLinkAngles[descr].append ( [a1, a2, a3, angle, esd] )


def GetTorsion ( l ) :
    global mrTorsions
    ts = l.split()
    if len (ts) == 9 :
        # LYS      chi1     N      CA     CB     CG       180.0     15.0     3
        rtype, torId, a1, a2, a3, a4, angle, esd, period = ts
        angle = float ( angle ) * numpy.pi / 180.0
        esd = float ( esd ) * numpy.pi / 180.0
        period = int ( period )

        #if rtype == "ARG" or rtype == "ALA" :
        #    print rtype, torId, a1, a2, a3, a4, angle, esd, period

        if not rtype in mrTorsions :
            mrTorsions[rtype] = []
        mrTorsions[rtype].append ( [a1, a2, a3, a4, angle, esd, period] )



def ReadGeo ( LOG = 0 ) :

    global mrAngles, mrBonds, mrAtoms, mrLinkBonds, mrLinkAngles, mrTorsions

    import os
    dir_path = os.path.dirname(os.path.realpath(__file__))
    if LOG : print dir_path

    li = 0
    currentHeader, currentType = [], None
    mrAngles, numAngles = {}, 0
    mrBonds, numBonds = {}, 0
    mrAtoms, numAtoms = {}, 0
    mrLinkBonds, numLinkBonds = {}, 0
    mrLinkAngles, numLinkAngles = {}, 0
    mrTorsions, numTorsions = {}, 0

    fp = open ( dir_path + "/standard_geometry.cif.txt")
    while 1 :

        l = fp.readline(); li+=1;
        if not l : break

        if l[0:len("loop_")] == "loop_" :

            if LOG : print "\n%d:" % li,
            currentHeader, currentType = [], None

            while 1 :
                l = fp.readline(); li+=1;
                if not l : break

                if l[0:1] == "_" :
                    currentType = l.split(".")[0].strip()
                    label = l.split(".")[1].strip()
                    label = label.replace ( "number", "#" )
                    currentHeader.append ( label )
                else :
                    if LOG :
                        print "[%s]" % currentType, "|".join(currentHeader)
                        print ""
                    break

        if not l : break

        if currentType == "_chem_comp_atom" :
            GetAtom ( l ); numAtoms += 1

        elif currentType == "_chem_comp_angle" :
            GetAngle ( l ); numAngles += 1

        elif currentType == "_chem_link_angle" :
            GetLinkAngle ( l ); numLinkAngles += 1

        elif currentType == "_chem_comp_bond" :
            GetBond ( l ); numBonds += 1

        elif currentType == "_chem_link_bond" :
            GetLinkBond ( l ); numLinkBonds += 1

        elif currentType == "_chem_comp_tor" :
            GetTorsion ( l ); numTorsions += 1

    fp.close()

    if LOG :
        print li, "lines"
        print numAtoms, "atoms ->", len(mrAtoms)
        print numBonds, "bonds ->", len(mrBonds)
        print numLinkBonds, "link bonds ->", len(mrLinkBonds)
        print numAngles, "angles ->", len(mrAngles)
        print numLinkAngles, "link angles ->", len(mrLinkAngles)
        print numTorsions, "torsions ->", len(mrTorsions)




def ResParams ( res, log=False ) :

    global mrBonds, mrAngles, mrTorsions
    if mrBonds == None :
        #print " - getting geo..."
        ReadGeo ()

    rBonds, rAngles, rTorsions = [], [], []

    if res.type in mrBonds :
        for b in mrBonds[res.type] :

            a1, a2, dist, esd, descr = b

            if res.type == "CYS" :
                if descr == "."  : pass
                elif descr == "SH" : pass
                else : continue

            if res.type == "HIS" :
                if descr == "." : pass
                elif descr == "HISE" : pass
                else : continue

            if res.type == "PRO" :
                if descr == "." : pass
                elif descr == "trans" : pass
                else : continue

            at1, at2 = None, None
            if a1 in res.atomsMap :
                at1 = res.atomsMap[a1][0]
            if a2 in res.atomsMap :
                at2 = res.atomsMap[a2][0]

            if at1 and at2 :
                rBonds.append ( [at1, at2, dist, esd] )
                if log : print " -b: ", at1.name, at2.name, dist, esd, descr

    else :
        print " - bonds for res %d.%s:%s not found" % (res.id.position, res.id.chainId, res.type)

    if res.type in mrAngles :
        for a in mrAngles[res.type] :
            a1, a2, a3, angle, esd, descr = a

            if res.type == "CYS" :
                if descr == "."  : pass
                elif descr == "SH" : pass
                else : continue

            if res.type == "HIS" :
                if descr == "." : pass
                elif descr == "HISE" : pass
                else : continue

            if res.type == "PRO" :
                if descr == "." : pass
                elif descr == "trans" : pass
                else : continue

            at1, at2, at3 = None, None, None
            if a1 in res.atomsMap :
                at1 = res.atomsMap[a1][0]
            if a2 in res.atomsMap :
                at2 = res.atomsMap[a2][0]
            if a3 in res.atomsMap :
                at3 = res.atomsMap[a3][0]
            if at1 and at2 and at3 :
                rAngles.append ( [at1, at2, at3, angle, esd] )
                if log : print " -a: ", at1.name, at2.name, at3.name, angle, esd, descr
    else :
        print " - angles for res %d.%s:%s not found" % (res.id.position, res.id.chainId, res.type)


    if res.type in mrTorsions :
        for a in mrTorsions[res.type] :
            a1, a2, a3, a4, angle, esd, period = a

            at1, at2, at3, at4 = None, None, None, None
            if a1 in res.atomsMap :
                at1 = res.atomsMap[a1][0]
            if a2 in res.atomsMap :
                at2 = res.atomsMap[a2][0]
            if a3 in res.atomsMap :
                at3 = res.atomsMap[a3][0]
            if a4 in res.atomsMap :
                at4 = res.atomsMap[a4][0]

            if at1 and at2 and at3 and at4 :
                rTorsions.append ( [at1, at2, at3, at4, angle, esd, period] )
                if log : print " -t: ", At(at1), at2.name, at3.name, at4.name, angle, esd, period
    else :
        print " - torsions for res %d.%s:%s not found" % (res.id.position, res.id.chainId, res.type)


    return rBonds, rAngles, rTorsions




def ConnectRes (r1, r2, allBonds, allAngles) :

    t1 = 'NOT GLY PRO'
    if r1.type == "PRO" or r1.type == "GLY" :
        t1 = r1.type
    a1, a2, dist, esd = mrLinkBonds[t1]
    at1 = r1.atomsMap['C'][0]
    at2 = r2.atomsMap['N'][0]
    allBonds.append ( [at1, at2, dist, esd] )
    #print " -bl: ", at1.name, at2.name, dist, esd

    t1, t2 = "NOT GLY PRO", "NOT PRO GLY"
    #if r1.type == "PRO" or pres.type == "PRO" or res.type == "GLY" or res.type == "GLY" :
    #    t1, t2 = "NOT GLY PRO", "NOT PRO GLY"

    angle, esd = mrLinkAngles[t1]["CA_C_N"]
    at1, at2, at3 = r1.atomsMap['CA'][0], r1.atomsMap['C'][0], r2.atomsMap['N'][0]
    allAngles.append ( [at1, at2, at3, angle, esd]  )

    angle, esd = mrLinkAngles[t2]["O_C_N"]
    at1, at2, at3 = r1.atomsMap['O'][0], r1.atomsMap['C'][0], r2.atomsMap['N'][0]
    allAngles.append ( [at1, at2, at3, angle, esd]  )

    angle, esd = mrLinkAngles[t1]["C_N_CA"]
    at1, at2, at3 = r1.atomsMap['C'][0], r2.atomsMap['N'][0], r2.atomsMap['CA'][0]
    allAngles.append ( [at1, at2, at3, angle, esd]  )



def Refine ( ress, dmap ) :

    global mrBonds, mrAngles

    if mrBonds == None :
        print " - getting geo..."
        ReadGeo ()

    #print "Refining %d residues" % len(ress)

    atoms = []
    atomsMap = {}
    molMap = {}
    resMap = {}
    refAtoms = {}
    for r in ress :
        molMap[r.molecule] = 1
        for at in r.atoms :
            atoms.append ( at )
            atomsMap[at] = 1
            refAtoms[at] = 1


    bonds = []
    mi_ci_ri = {}
    for mol in molMap.keys() :

        for bond in mol.bonds :
            if bond.atoms[0] in atomsMap and bond.atoms[1] in atomsMap :
                bonds.append ( bond )

        if not mol.id in mi_ci_ri :
            mi_ci_ri[mol.id] = {}

        for res in mol.residues :
            if not res.id.chainId in mi_ci_ri[mol.id] :
                mi_ci_ri[mol.id][res.id.chainId] = {}
            if res.id.position in mi_ci_ri[mol.id][res.id.chainId] :
                print " - multiple residues with position %d in %d.%s" % (res.id.position, mol.id, res.id.chainId)
            else :
                mi_ci_ri[mol.id][res.id.chainId][res.id.position] = res


    nextRes = {}
    refBonds, refAngles, refTorsions = [], [], []
    for res in ress :

        for at in r.atoms :
            at.M = 1.0

        rBonds, rAngles, rTorsions = ResParams ( res )
        refBonds.extend ( rBonds )
        refAngles.extend ( rAngles )
        refTorsions.extend ( rTorsions )
        #print " - %d.%s %s - %d bonds, %d angles" % (r.id.position, r.id.chainId, r.type, len(rBonds), len(rAngles))

        if res.id.position-1 in mi_ci_ri[mol.id][res.id.chainId] :
            pres = mi_ci_ri[mol.id][res.id.chainId][res.id.position-1]

            # add connection angles only once...
            if pres in nextRes : continue
            nextRes[pres] = res

            #print " - con <- %d %s" % (res.id.position-1, pres.type)
            ConnectRes ( pres, res, refBonds, refAngles )

            for at in pres.atoms :
                refAtoms[at] = 1
                at.M = 0.0


        if res.id.position+1 in mi_ci_ri[mol.id][res.id.chainId] :
            nres = mi_ci_ri[mol.id][res.id.chainId][res.id.position+1]

            # add connection angles only once...
            if res in nextRes : continue
            nextRes[res] = nres

            #print " - con -> %d %s" % (res.id.position+1, nres.type)
            ConnectRes ( res, nres, refBonds, refAngles )

            for at in nres.atoms :
                refAtoms[at] = 1
                at.M = 0.0


    #print " - %d bonds, %d angles" % (len(rBonds), len(rAngles))
    params = { 'doBonds' : 1,
                'doAngles' : 1,
                'doTorsions' : 1,
                'doMap' : 1 }


    if not hasattr ( dmap, 'maxd' ) :
        M = dmap.data.full_matrix()
        avg, std = numpy.average(M), numpy.std(M)
        #maxM = numpy.max(M)
        #minM = numpy.min(M)
        dmap.maxd = min ( avg+std*10.0, numpy.max(M) )
        dmap.mind = max ( avg-std*1.0, numpy.min(M) )


    #DoRef ( refAtoms.keys(), refBonds, refAngles, params, N=10, log=0 )
    DoRef ( dmap, refAtoms.keys(), refBonds, refAngles, refTorsions, params )





def DoRef ( dmap, atoms, allBonds, allAngles, allTorsions, params = {}, N=100, log=1 ) :

    doBonds = 'doBonds' in params and params['doBonds'] == 1
    doAngles = 'doAngles' in params and params['doAngles'] == 1
    doTorsions = 'doTorsions' in params and params['doTorsions'] == 1
    doMap = 'doMap' in params and params['doMap'] == 1

    print " - in map %s" % dmap.name


    for i in range ( N ) :

        for at in atoms :
            #at.G = chimera.Vector(0,0,0)
            at.P = at.coord()
            #at.P = chimera.Point(0,0,0)
            at.G = chimera.Vector(0,0,0)

        #if log : print ""

        if doBonds :
            for b in allBonds :
                at1, at2, D, esd = b
                BondG ( at1, at2, D, esd, F=0.1 )

        #if log : print ""

        if doAngles :
            for a in allAngles :
                at1, at2, at3, A, esd = a
                AngleG ( at1, at2, at3, A, esd, F=0.01 )
                #break

        if doTorsions :
            for t in allTorsions :
                at1, at2, at3, at4, angle, esd, period = t

                #if at1.name == "N" :
                    #print "tor: %s, %s, %s, %s - angle %.3f, per %d" % ( At(at1), at2.name, at3.name, at4.name, angle*180.0/numpy.pi, period )
                    #RotDih ( at1, at2, at3, at4, angle, esd, period, dmap )

                TorsionG ( at1, at2, at3, at4, angle, esd, period, F=0.01 )

        if doMap :
            for at in atoms :
                MapG ( at, dmap, F=0.01 )
                #break


        for at in atoms :
            #at.setCoord ( at.coord() + at.G )
            at.P = at.P + at.G * at.M
            at.setCoord ( at.P )
            #pass

        if i % 10 == 0 :
            print ".",
    print ""



    #if log : print ""

    dB = 0.0
    if doBonds :
        for b in allBonds :
            at1, at2, D, esd = b
            dB += BondG ( at1, at2, D, esd, F=0 )

            #if log : print " - b %s _ %s %.3f (%.3f)" % (At(at1), At(at2), D, d)

    #if log : print ""

    dA = 0.0
    if doAngles :
        for a in allAngles :

            at1, at2, at3, A, esd = a
            dA += AngleG ( at1, at2, at3, A, esd, F=0 )

            #if log : print " - a %s _ %s _ %s %.3f (%.3f)" % (At(at1), At(at2), At(at3), A, d)

    dTor = 0.0
    if doTorsions :
        for t in allTorsions :
            at1, at2, at3, at4, angle, esd, period = t

            #if at1.name == "N" :
            #    print "tor: %s, %s, %s, %s - angle %.3f, per %d" % ( At(at1), at2.name, at3.name, at4.name, angle*180.0/numpy.pi, period )
            dTor += TorsionG ( at1, at2, at3, at4, angle, esd, period, F=0 )

    dMap = 0.0
    if doMap :
        for at in atoms :
            dMap += MapG ( at, dmap, F=0 )


    RefDisp ( atoms, dmap )

    #if log : print ""

    print " - ref %d - bonds: %.3f (%d), angles: %.3f (%d), torsions: %.3f (%d), map %.4f" % (N, dB, len(allBonds), dA, len(allAngles), dTor, len(allTorsions), dMap )




def MapG ( at, dmap, F=0.1 ) :

    #print at.P, type(at.P)
    tp = at.P
    #tp = at.molecule.openState.xform.apply ( at.P )
    #tp = dmap.openState.xform.inverse().apply ( tp )
    P = tp.data()

    if F <= 0.0 :
        return dmap.interpolated_values ( [P], at.molecule.openState.xform ) [0]

    dx, dy, dz = dmap.data.step
    pts = numpy.array ( [P] * 6 )

    #print pts

    #dx *= 0.1

    pts[0][0] -= dx; pts[1][0] += dx
    pts[2][1] -= dx; pts[3][1] += dx
    pts[4][2] -= dx; pts[5][2] += dx

    #print pts, type(pts)

    dmap.ctr = numpy.array(dmap.data.origin) + numpy.array(dmap.data.size)/2.0*numpy.array(dmap.data.step)


    vs = dmap.interpolated_values ( pts, at.molecule.openState.xform )

    #print mapvs, type(mapvs)

    vs = (vs - dmap.mind)/(dmap.maxd-dmap.mind)

    #print mapvs, type(mapvs)

    dx2 = 1.0/(dx*2.0)

    G = chimera.Vector ( (vs[1] - vs[0])*dx2, (vs[3] - vs[2])*dx2, (vs[5] - vs[4])*dx2  )
    #print G

    #G = dmap.ctr - numpy.array( tp )
    #G = chimera.Vector ( *G )

    #G = dmap.openState.xform.apply ( G )
    #G = at.molecule.openState.xform.inverse().apply ( G )
    #print G

    at.dg = G
    at.G += G * F




def RefDisp ( ats, dmap, mname = "RefDisp" ) :

    import axes; reload ( axes )

    smod = None
    for m in chimera.openModels.list() :
        if m.name == mname :
            smod = m

    if smod :
        for sp in smod.surfacePieces :
            smod.removePiece ( sp )

    else :
        import _surface
        smod = _surface.SurfaceModel()
        smod.name = mname
        chimera.openModels.add ( [smod] )


    #axes.SphereMesh (r=0.3, div=10, color=(.9,.5,.5,1.0), pos=dmap.ctr, mol=smod)


    for at in ats :

        if hasattr ( at, 'dg' ) :

            I = smod.openState.xform.inverse()
            P = I.apply(at.xformCoord())
            G = I.apply(at.molecule.openState.xform.apply(at.dg))

            #axes.AddArrow4 ( P, G, G.length, clr=(.5,.8,.5,1), rad=0.1, mol=smod, hrad=0.15, hlen=0.05 )
            if at.dg.length > 1e-3 :
                axes.AddArrow4 ( at.P, at.dg, at.dg.length, clr=(.5,.8,.5,1), rad=0.1, mol=smod, hrad=0.15, hlen=0.05 )




def RotDih ( at1, at2, at3, at4, angle, esd, period, dmap ) :

    A = diha ( at1, at2, at3, at4 )
    print " - start: %.3f" % A

    for a in range ( 0, 360, 10 ) :

        print "%.1f\t" % a,

        TorsionG ( at1, at2, at3, at4, angle, esd, period, F=1.0 )

        bx = at3.P - at2.P; bx.normalize()
        xf = chimera.Xform.translation (at3.P.toVector()*1.0)
        xf.multiply ( chimera.Xform.rotation ( bx, 10.0 ) )
        xf.multiply ( chimera.Xform.translation (at3.P.toVector()*-1.0) )

        trp = xf.apply ( at4.P )
        #at4.setCoord ( trp )
        AddSpherePts ( [ trp.data() ], (0,1,0,1), .1, mname = "RAD points" )

        trp = at3.molecule.openState.xform.apply ( trp )

        mapv = dmap.interpolated_values ( [trp.data()], dmap.openState.xform )
        print "\t%.3f" % ( mapv[0] )



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


def At ( at ) :
    return "%d.%s(%s)_%s" % (at.residue.id.position, at.residue.id.chainId, at.residue.type, at.name)



def AddSpherePts ( pts, clr, rad, mname = "RAD points" ) :

    from chimera import elements, Coord, Atom, MolResId

    ptsMol = None
    for m in chimera.openModels.list() :
        if m.name == mname :
            ptsMol = m

    res = None
    if ptsMol == None:
        from chimera import Molecule, openModels
        ptsMol = Molecule()
        ptsMol.name = mname
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





def BondG ( at1, at2, D, esd, F=0.1 ) :

    v = at1.P - at2.P
    d = D - v.length

    #if log : print " - b0 %s-%s %.3f/%.3f" % (at1.name, at2.name, D, d)

    if F <= 0.0 :
        return d * d


    if v.length < 0.01 :
        v = chimera.Vector ( random.random(), random.random(),random.random() )

    v.normalize()
    at1.G += v * d * F
    at2.G -= v * d * F
    #at1.P = at1.P + v * d * 0.1 * at1.M
    #at2.P = at2.P - v * d * 0.1 * at2.M



def AngleG_ (at1, at2, at3, theta0, esd, F=0.001) :

    r12 = at1.P - at2.P; d12 = r12.length
    r32 = at3.P - at2.P; d32 = r32.length

    cos_theta = (r12 * r32)/(d12*d32);
    #theta0 = a->dEqAngle;
    theta = numpy.arccos(cos_theta);
    diff = theta - theta0


    #print " - a %s _ %s _ %s %.3f/%.3f (%.3f)" % (At(at1), At(at2), At(at3), theta0, theta, diff)

    #energy = k *diff*diff;
    #dEAngles += energy;

    if F <= 0.0 :
        return diff*diff

    d12inv = 1.0 / d12;
    d32inv = 1.0 / d32;

    #  Calculate constant factor 2k(theta-theta0)/sin(theta)
    sin_theta = numpy.sqrt(1.0 - cos_theta*cos_theta)
    #diff *= (-2.0 * k) / sin_theta;
    diff = diff * (-2.0) / sin_theta;
    c1 = diff * d12inv
    c2 = diff * d32inv

    #  Calculate the actual forces
    force1 = (r12*(d12inv*cos_theta) - r32*d32inv)*c1
    force2 = force1
    force3 = (r32*(d32inv*cos_theta) - r12*d12inv)*c2
    force2 += force3
    force2 *= -1;

    at1.G += force1 * F
    at2.G += force2 * F
    at3.G += force3 * F




def AngleG (at1, at2, at3, theta0, esd, F=0.001) :

    rij = at1.P - at2.P
    rik = at3.P - at2.P

    lij, lik, = rij.length, rik.length
    if lij < 1e-3 or lik < 1e-3 :
        return 0.0

    rij.normalize()
    rik.normalize()

    cost = rij * rik
    t = numpy.arccos(cost);
    d = t - theta0
    ge = 2.0 * d / numpy.sin(t);

    if F <= 0.0 :
        return d * d


    i0 = -rik[0]/lij - rij[0]/lik - cost * ( -rij[0]/lij - rik[0]/lik )
    i1 = -rik[1]/lij - rij[1]/lik - cost * ( -rij[1]/lij - rik[1]/lik )
    i2 = -rik[2]/lij - rij[2]/lik - cost * ( -rij[2]/lij - rik[2]/lik )

    j0 = rik[0]/lij - cost * ( -rij[0]/lij )
    j1 = rik[1]/lij - cost * (  rij[1]/lij )
    j2 = rik[2]/lij - cost * ( -rij[2]/lij )

    k0 = rij[0]/lik - cost * ( rik[0]/lik )
    k1 = rij[1]/lik - cost * ( rik[1]/lik )
    k2 = rij[2]/lik - cost * ( rik[2]/lik )

    fi, fj, fk = chimera.Vector(i0,i1,i2), chimera.Vector(j0,j1,j2), chimera.Vector(k0,k1,k2)

    at1.G += fj * (ge * F)
    at2.G += fi * (ge * F)
    at3.G += fk * (ge * F)



def TorsionG ( at1, at2, at3, at4, angle, esd, period, F=0.01 ) :

	# Calculate the vectors between atoms
	#Vector3d pos0 = d->a1->pos, pos1 = d->a2->pos, pos2 = d->a3->pos, pos3 = d->a4->pos;
	#Vector3d r12 = pos0-pos1, r23 = pos1-pos2, r34 = pos2-pos3;
    r12 = at1.P - at2.P; r23 = at2.P - at3.P; r34 = at3.P - at4.P

	# Calculate the cross products and distances
	#Vector3d A = r12.Cross(r23); double rA = A.Length();
	#Vector3d B = r23.Cross(r34); double rB = B.Length();
	#Vector3d C = r23.Cross(A); double rC = C.Length();
    A = chimera.cross ( r12, r23 ); rA = A.length
    B = chimera.cross ( r23, r34 ); rB = B.length
    C = chimera.cross ( r23, A ); rC = C.length

    if rA < 1e-3 or rB < 1e-3 or rC < 1e-3 :
        return 0.0

	# Calculate the sin and cos
	#double cos_phi = (A.Dot(B))/(rA*rB);
	#double sin_phi = (C.Dot(B))/(rC*rB);
    cos_phi = A * B / (rA * rB)
    sin_phi = C * B / (rC * rB)

	#double phi= -atan2(sin_phi,cos_phi);
    phi = - numpy.arctan2 ( sin_phi, cos_phi )


    K=0;		# energy
    K1=0;		# force

	# get the dihedral information
	#int multiplicity = (int) d->values.size();
    #multiplicity = 0

	#  Loop through the multiple parameter sets for this
	#  bond.  We will only loop more than once if this
	#  has multiple parameter sets from Charmm22
	#for (int mult_num=0; mult_num<multiplicity; mult_num++)
	#{
		# get angle information
		#double k = d->values[mult_num].dFC;
		#double delta = d->values[mult_num].dPhase;
		#int n = d->values[mult_num].dPeriodicity;

    delta = angle
    n = period
    k = 1.0

	# Calculate the energy
    if n > 0 :
    	# Periodicity is greater than 0, so use cos form
    	K += k*(1+numpy.cos(n*phi - delta));
    	K1 += -n*k*numpy.sin(n*phi - delta);

    else :
		# Periodicity is 0, so just use the harmonic form
        diff = phi - delta
        if diff < -M_PI :
            diff += TWOPI
        elif diff > M_PI :
            diff -= TWOPI;

        K += k*diff*diff;
        K1 += 2.0*k*diff;

    #print "tor: %s, %s, %s, %s - angle %.3f / %.3f, per %d" % ( At(at1), at2.name, at3.name, at4.name, angle, phi*180/numpy.pi, period ),
    #print "%.3f\t%.3f" % (phi*180/numpy.pi, (-K+1.0)),


    #dEDihedrals += K;
    #}

    if F <= 0.0 :
        return K

    # using F instead of k here
    K1 *= F

    #Vector3d f1,f2,f3;

    # Normalize B
    rB = 1.0/rB;
    B *= rB;

	# Next, we want to calculate the forces.  In order
	# to do that, we first need to figure out whether the
	# sin or cos form will be more stable.  For this,
	# just look at the value of phi
    #print "sin_phi - %.3f - %.3f" % (sin_phi, abs(sin_phi))
    if abs(sin_phi) > 0.1 :
		# use the sin version to avoid 1/cos terms
		# Normalize A
        rA = 1.0/rA;
        A *= rA;
        dcosdA = (cos_phi*A-B)*rA;
        dcosdB = (cos_phi*B-A)*rB;

        #print "sin"

        K1 = K1/sin_phi;

        f1, f2, f3 = chimera.Vector(), chimera.Vector(), chimera.Vector()

        f1.x = K1*(r23.y*dcosdA.z - r23.z*dcosdA.y)
        f1.y = K1*(r23.z*dcosdA.x - r23.x*dcosdA.z)
        f1.z = K1*(r23.x*dcosdA.y - r23.y*dcosdA.x)

        f3.x = K1*(r23.z*dcosdB.y - r23.y*dcosdB.z)
        f3.y = K1*(r23.x*dcosdB.z - r23.z*dcosdB.x)
        f3.z = K1*(r23.y*dcosdB.x - r23.x*dcosdB.y)

        f2.x = K1*(r12.z*dcosdA.y - r12.y*dcosdA.z + r34.y*dcosdB.z - r34.z*dcosdB.y)
        f2.y = K1*(r12.x*dcosdA.z - r12.z*dcosdA.x + r34.z*dcosdB.x - r34.x*dcosdB.z)
        f2.z = K1*(r12.y*dcosdA.x - r12.x*dcosdA.y + r34.x*dcosdB.y - r34.y*dcosdB.x)

    else :

        # This angle is closer to 0 or 180 than it is to
        # 90, so use the cos version to avoid 1/sin terms

        #print "cos"

        # Normalize C
        rC = 1.0/rC
        C *= rC
        dsindC = (sin_phi*C-B)*rC
        dsindB = (sin_phi*B-C)*rB

        K1 = -K1/cos_phi

        f1, f2, f3 = chimera.Vector(), chimera.Vector(), chimera.Vector()

        f1.x = K1*((r23.y*r23.y + r23.z*r23.z)*dsindC.x - r23.x*r23.y*dsindC.y - r23.x*r23.z*dsindC.z)
        f1.y = K1*((r23.z*r23.z + r23.x*r23.x)*dsindC.y - r23.y*r23.z*dsindC.z - r23.y*r23.x*dsindC.x)
        f1.z = K1*((r23.x*r23.x + r23.y*r23.y)*dsindC.z - r23.z*r23.x*dsindC.x - r23.z*r23.y*dsindC.y)

        #f3 = K1*dsindB.Cross(r23)
        f3 = chimera.cross ( dsindB, r23 ) * K1

        f2.x = K1*(-(r23.y*r12.y + r23.z*r12.z)*dsindC.x +(2.0*r23.x*r12.y - r12.x*r23.y)*dsindC.y
        	+(2.0*r23.x*r12.z - r12.x*r23.z)*dsindC.z +dsindB.z*r34.y - dsindB.y*r34.z)
        f2.y = K1*(-(r23.z*r12.z + r23.x*r12.x)*dsindC.y+(2.0*r23.y*r12.z - r12.y*r23.z)*dsindC.z
        	+(2.0*r23.y*r12.x - r12.y*r23.x)*dsindC.x+dsindB.x*r34.z - dsindB.z*r34.x)
        f2.z = K1*(-(r23.x*r12.x + r23.y*r12.y)*dsindC.z +(2.0*r23.z*r12.x - r12.z*r23.x)*dsindC.x
        	+(2.0*r23.z*r12.y - r12.z*r23.y)*dsindC.y+dsindB.y*r34.x - dsindB.x*r34.y)


    # store the forces
    at1.G += f1
    at2.G += f2 - f1
    at3.G += f3 - f2
    at4.G += -f3








def AddResN ( rtype, atRes ) :


    global mrBonds, mrAngles, mrAtoms

    if mrBonds == None :
        print " - getting geo..."
        ReadGeo ()


    if len(rtype) == 1 :
        from chimera.resCode import protein1to3
        try :
            rtype = protein1to3[rtype]
        except :
            print "Adding - %s - unknown"

    ats = mrAtoms[rtype]

    addAts = []
    for at in ats :
        rtype, atomId, el, charge, x, y, z = at
        if el == "H" :
            continue
        addAts.append (at)

    mol = atRes.molecule
    chainId = atRes.id.chainId
    atPos = atRes.id.position

    print "Adding - %s - %d atoms, to %s, chaind %s, ri %d" % (rtype, len(addAts), mol.name, chainId, atPos)

    ress = []
    for r in mol.residues :
        if r.id.chainId == chainId :
            ress.append ( [r.id.position, r] )

    ress.sort()
    ress.reverse()

    print " - %d res in chain, first %d, last %d" % (len(ress), ress[0][1].id.position, ress[-1][1].id.position)


    aMap = {}

    bonds = []
    for bond in mol.bonds :
        bonds.append ( [bond.atoms[0], bond.atoms[1]] )

    for i, riRes in enumerate(ress) :

        ri, res = riRes

        if res.id.position >= atPos :
            ri = ri + 1
        else :
            break

        print " - at %d -> %d" % (res.id.position, ri)

        nres = mol.newResidue (res.type, chimera.MolResId(chainId, ri))
        for at in res.atoms :
            nat = mol.newAtom (at.name, chimera.Element(at.element.number))
            aMap[at] = nat
            nres.addAtom( nat )
            nat.setCoord ( at.coord() )
            nat.drawMode = nat.EndCap
            nat.display = True

        nres.isHelix = res.isHelix
        nres.isHet = res.isHet
        nres.isSheet = res.isSheet
        nres.isStrand = res.isStrand
        nres.ribbonDisplay = True
        nres.ribbonDrawMode = 2

        mol.deleteResidue ( res )


    nbondsAdded = 0
    for a1, a2 in bonds :
        if not a1 in aMap :
            continue
        if not a2 in aMap :
            continue
        nat1, nat2 = aMap[a1], aMap[a2]
        if abs(nat1.residue.id.position - nat2.residue.id.position) <= 1 :
            nb = mol.newBond ( aMap[a1], aMap[a2] )
            nb.display = nb.Smart
            nb.drawMode = nb.Stick
            nbondsAdded += 1

    print " - added %d bonds" % nbondsAdded

























# end
