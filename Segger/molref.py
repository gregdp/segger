

import chimera
import random
import numpy

mrAtoms = None
mrBonds = None
mrLinkBonds = None
mrAngles = None
mrLinkAngles = None


def ReadGeo ( LOG = 0 ) :

    global mrAngles, mrBonds, mrAtoms, mrLinkBonds,  mrLinkAngles

    import os
    dir_path = os.path.dirname(os.path.realpath(__file__))
    if LOG : print dir_path

    f = open ( dir_path + "/standard_geometry.cif.txt")
    li = 0
    getHead, head, gtype = 0, [], None
    mrAngles, numAngles = {}, 0
    mrBonds, numBonds = {}, 0
    mrAtoms, numAtoms = {}, 0
    mrLinkBonds, numLinkBonds = {}, 0
    mrLinkAngles, numLinkAngles = {}, 0
    angMap = {}
    for l in f :

        if getHead :
            if l[0:1] == "_" :
                gtype = l.split(".")[0].strip()
                label = l.split(".")[-1].strip()
                label = label.replace ( "number", "#" )
                head.append ( label )
            else :
                if LOG : print "head-", gtype
                #for h in head :
                #    print h,
                #print ""
                #print ""
                getHead = 0

        if not getHead:

            if gtype == "_chem_comp_atom" :

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
                    numAtoms += 1


            if gtype == "_chem_comp_angle" :

                ts = l.split()
                if len (ts) == 8 :
                    rtype, a1, a2, a3, angle, esd, exception, descr = ts
                    angle = float ( angle )
                    esd = float ( esd )

                    #if rtype == "ALA" :
                    #    print rtype, a1, a2, a3, angle, esd

                    descr = descr.replace ( '"', "" )

                    if not rtype in mrAngles :
                        mrAngles[rtype] = []
                    mrAngles[rtype].append ( [a1, a2, a3, angle, esd, descr] )
                    numAngles += 1

            if gtype == "_chem_link_angle" :

                ts = l.split()
                if len (ts) > 5 :
                    a1, a2, a3, angle, esd, descr = ts[0], ts[1], ts[2], ts[3], ts[4], ""
                    angle = float ( angle )
                    esd = float ( esd )

                    for i in range ( 5, len(ts) ) :
                        descr = descr + ts[i] + " "

                    descr = descr[0:-1].replace( '"', '' )

                    if not descr in mrLinkAngles :
                        mrLinkAngles[descr] = {}

                    mrLinkAngles[descr]["%s_%s_%s"%(a1,a2,a3)] = [angle, esd]

                    #print " - angle", "%s_%s_%s : %.3f, %.3f"%(a1,a2,a3,angle,esd)

                    #mrLinkAngles[descr].append ( [a1, a2, a3, angle, esd] )
                    numLinkAngles += 1


            if gtype == "_chem_comp_bond" :

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
                    numBonds += 1

            if gtype == "_chem_link_bond" :

                ts = l.split()
                if len (ts) > 4 :
                    a1, a2, dist, esd = ts[0], ts[1], ts[2], ts[3]

                    dist = float ( dist )
                    esd = float ( esd )

                    descr = ""
                    for i in range ( 4, len(ts) ) :
                        descr = descr + ts[i] + " "
                    descr = descr[0:-1].replace( '"', '' )

                    mrLinkBonds[descr] = [a1, a2, dist, esd]
                    numLinkBonds += 1


        if l[0:len("loop_")] == "loop_" :
            if LOG : print li, # "loop",
            getHead, head, gtype = 1, [], None

        li += 1

    if LOG :
        print li, "lines"
        print numAngles, "angles, ", len(mrAngles), "res"
        print numLinkAngles, "link angles, ", len(mrLinkAngles), "res"
        print numBonds, "bonds, ", len(mrBonds), "res"
        print numLinkBonds, "link bonds, ", len(mrLinkBonds), "res"
        print numAtoms, "atoms, ", len(mrAtoms), "res"




def ResParams ( res, log=False ) :

    global mrBonds, mrAngles
    if mrBonds == None :
        print " - getting geo..."
        ReadGeo ( LOG = 0 )

    rBonds, rAngles = [], []

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
                if log :
                    print " -b: ", at1.name, at2.name, dist, esd, descr

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
                if log :
                    print " -a: ", at1.name, at2.name, at3.name, angle, esd, descr
    else :
        print " - angles for res %d.%s:%s not found" % (res.id.position, res.id.chainId, res.type)

    return rBonds, rAngles




def Refine ( ress ) :

    global mrBonds, mrAngles

    if mrBonds == None :
        print " - getting geo..."
        ReadGeo ( LOG = 0 )

    #print "Refining %d residues" % len(ress)

    atoms = []
    atomsMap = {}
    molMap = {}
    resMap = {}
    for r in ress :
        molMap[r.molecule] = 1
        for at in r.atoms :
            atoms.append ( at )
            atomsMap[at] = 1


    for mol in molMap.keys() :
        for at in mol.atoms :
            at.M = 0.0
            at.P = at.coord()

    for r in ress :
        for at in r.atoms :
            at.M = 1.0


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



    allBonds, allAngles = [], []
    for res in ress :

        rBonds, rAngles = ResParams ( res, log=False )
        allBonds.extend ( rBonds )
        allAngles.extend ( rAngles )
        #print " - %d.%s %s - %d bonds, %d angles" % (r.id.position, r.id.chainId, r.type, len(rBonds), len(rAngles))

        if res.id.position-1 in mi_ci_ri[mol.id][res.id.chainId] :
            pres = mi_ci_ri[mol.id][res.id.chainId][res.id.position-1]
            print " - con res %d %s" % (res.id.position-1, pres.type)

            a1, a2, dist, esd = mrLinkBonds['NOT GLY PRO']
            at1 = cres.atomsMap['C'][0]
            at2 = res.atomsMap['N'][0]
            allBonds.append ( [at1, at2, dist, esd] )
            #print " -bl: ", at1.name, at2.name, dist, esd

            t1, t2 = "NOT GLY PRO", "NOT PRO GLY"
            if res.type == "PRO" or pres.type == "PRO" or res.type == "GLY" or res.type == "GLY" :
                t1, t2 = "NOT GLY PRO", "NOT PRO GLY"

            angle, esd = mrLinkAngles[t1]["CA_C_N"]
            at1, at2, at3 = pres.atomsMap['CA'][0], pres.atomsMap['C'][0], res.atomsMap['N'][0]
            allAngles.append ( [at1, at2, at3, angle, esd]  )

            angle, esd = mrLinkAngles[t2]["O_C_N"]
            at1, at2, at3 = pres.atomsMap['O'][0], pres.atomsMap['C'][0], res.atomsMap['N'][0]
            allAngles.append ( [at1, at2, at3, angle, esd]  )

            angle, esd = mrLinkAngles[t1]["C_N_CA"]
            at1, at2, at3 = pres.atomsMap['C'][0], res.atomsMap['N'][0], res.atomsMap['CA'][0]
            allAngles.append ( [at1, at2, at3, angle, esd]  )


        if res.id.position+1 in mi_ci_ri[mol.id][res.id.chainId] :
            cres = mi_ci_ri[mol.id][res.id.chainId][res.id.position+1]
            print " - con res %d %s" % (res.id.position+1, cres.type)


    #print " - %d bonds, %d angles" % (len(rBonds), len(rAngles))
    DoRef ( atoms, allBonds, allAngles, doBonds=1, doAngles=1, N=100, log=1 )





def DoRef ( atoms, allBonds, allAngles, doBonds=1, doAngles=1, N=10, log=1 ) :


    for i in range ( N ) :

        for at in atoms :
            #at.G = chimera.Vector(0,0,0)
            at.P = at.coord()

        #if log : print ""

        if doBonds :
            for b in allBonds :
                at1, at2, D, esd = b
                v = at1.P - at2.P
                d = D - v.length

                #if log : print " - b0 %s-%s %.3f/%.3f" % (at1.name, at2.name, D, d)

                if v.length < 0.01 :
                    v = chimera.Vector ( random.random(), random.random(),random.random() )

                v.normalize()
                #at1.G += v * d * 0.1
                #at2.G -= v * d * 0.1
                at1.P = at1.P + v * d * 0.1 * at1.M
                at2.P = at2.P - v * d * 0.1 * at2.M

        #if log : print ""

        if doAngles :
            for a in allAngles :
                at1, at2, at3, A, esd = a
                #at1Pt, at2Pt, at3Pt = at1.coord(), at2.coord(), at3.coord()
                v1 = at1.P - at2.P; v1L = v1.length
                v3 = at3.P - at2.P; v3L = v3.length
                if v1L < 1e-5 : continue
                if v3L < 1e-5 : continue
                #v1.normalize()
                #v3.normalize()
                a = numpy.arccos ( v1 * v3 / (v1L * v3L) ) * 180.0 / numpy.pi
                d = a-A

                # if log : print " - a %s-%s-%s %.3f/%.3f" % (at1.name, at2.name, at3.name, A, d)

                X = chimera.cross ( v1, v3 )
                if X.length < 1e-5 :
                    X = chimera.Vector ( random.random(), random.random(),random.random() )
                X.normalize()

                #a1v = chimera.cross ( x, v1 )
                #a3v = chimera.cross ( v2, x )
                ad = abs(d)
                if at1.M > 0.0 :
                    at1.P = at2.P + chimera.Xform.rotation ( X, max(min(d * 0.1,5.0),-5.0) ).apply (v1)
                if at3.M > 0.0 :
                    at3.P = at2.P + chimera.Xform.rotation ( X, -max(min(d * 0.1,5.0),-5.0) ).apply (v3)

                #at1.G += a1v * d * 0.01
                #at3.G += a3v * d * 0.01

        for at in atoms :
            #at.setCoord ( at.coord() + at.G )
            at.setCoord ( at.P )

        if i % 10 == 0 :
            print ".",
    print ""



    if log : print ""

    dB = 0.0
    if doBonds :
        for b in allBonds :
            at1, at2, D, esd = b
            v = at1.coord() - at2.coord()
            d = D - v.length
            dB += abs(d)
            if log : print " - b %d.%s:%s-%d.%s:%s %.3f (%.3f)" % (at1.residue.id.position, at1.residue.id.chainId, at1.name, at2.residue.id.position, at2.residue.id.chainId, at2.name, D, d)

    if log : print ""

    dA = 0.0
    if doAngles :
        for a in allAngles :
            at1, at2, at3, A, esd = a
            v1 = at1.coord() - at2.coord();
            v3 = at3.coord() - at2.coord();
            if v1.length < 1e-5 : continue
            if v3.length < 1e-5 : continue
            v1.normalize()
            v3.normalize()
            a = numpy.arccos ( v1 * v3 ) * 180.0 / numpy.pi
            d = A - a
            dA += abs(d)
            if log : print " - a %s-%s-%s %.3f (%.3f)" % (at1.name, at2.name, at3.name, A, d)

    if log : print ""

    print " - ref %d - bonds: %.3f (%d), angles: %.3f (%d)" % (N, dB, len(allBonds), dA, len(allAngles) )





def AddResN ( rtype, atRes ) :


    global mrBonds, mrAngles, mrAtoms

    if mrBonds == None :
        print " - getting geo..."
        ReadGeo ( LOG = 0 )


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
