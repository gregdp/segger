

import chimera
import random
import numpy
import _multiscale
import VolumeData
import FitMap
import _contour
import Matrix
import _distances

from chimera.resCode import nucleic3to1
from chimera.resCode import protein3to1, protein1to3
protein3to1['HSD'] = protein3to1['HIS']
protein3to1['HSE'] = protein3to1['HIS']

nucleic1to3 = { 'G':'GUA', 'A':'ADE', 'U':'URA', 'C':'CYT' }


atomColors = {'C' : chimera.MaterialColor (0.565,0.565,0.565),
            'Cbb' : chimera.MaterialColor (0.2,0.6,0.2),
            'S' : chimera.MaterialColor (1.000,1.000,0.188),
            'O' : chimera.MaterialColor (1.000,0.051,0.051),
            'N' : chimera.MaterialColor (0.188,0.314,0.973),
            'P' : chimera.MaterialColor (1.0, 0.502, 0.0),
            'H' : chimera.MaterialColor (0.9,.9,.9),
            ' ' : chimera.MaterialColor (0.2,1,.2),
            "MG" : chimera.MaterialColor (0,1,0),
            "NA" : chimera.MaterialColor (.6,.3,.6),
            "CL" : chimera.MaterialColor (.2,.6,.2),
            "CA" : chimera.MaterialColor (.4,.4,.6),
            "ZN" : chimera.MaterialColor (.2,.8,.2),
            "MN" : chimera.MaterialColor (.4,.4,.6),
            "FE" : chimera.MaterialColor (.4,.4,.6),
            "CO" : chimera.MaterialColor (.4,.4,.6),
            "NI" : chimera.MaterialColor (.4,.4,.6)
}



M_PI = numpy.pi
TWOPI = 2.0 * numpy.pi

class Params :

    def __init__ ( self ) :

        self.atoms = {}
        self.bonds = {}
        self.linkBonds = {}
        self.angles = {}
        self.linkAngles = {}
        self.torsions = {}
        self.planes = {}

        print "init params"


class RefineParams :

    def __init__ ( self ) :

        self.bonds = []
        self.angles = []
        self.torsions = []
        self.planes = []


mrPar = None


def GetProps ( l, label, props, toMap ) :
    global mrPar

    ts = l.split()

    if 1 or len (ts) == len(props) :
        #rtype, a1, a2, a3, angle, esd, exception, descr = ts

        par = {}
        for i in range ( len(props) ) :

            prop = props[i]
            if i >= len(ts) :
                print " -x- %s : length mismatch" % label
                print l,
                print props
                break

            val = ""
            if ts[i][0] == '"' :
                while i < len(ts) :
                    if len(val) == 0 :
                        val = ts[i]
                    else :
                        val += " " + ts[i]
                    if val[-1] == '"' :
                        break
                    i += 1
                val = val.replace ( '"', '' )
                #print l,
                #print " - val:", val
            else :
                val = ts[i]

            par[prop] = val

            try :
                numI = int ( val )
                par[prop] = numI
            except :
                pass

            try :
                numF = float ( val )
                par[prop] = numF
            except :
                pass

        if "comp_id" in par :
            rtype = par["comp_id"]
            if not rtype in toMap :
                toMap[rtype] = []
            toMap[rtype].append ( par )
        elif "link_id" in par :
            #rtype = par["link_id"]
            linkI = "%d" % len(toMap)
            toMap[linkI] = par
        else :
            print " -x- %s : no comp_id" % label
            print l,
            print props

    else :
        print " -x- %s : length mismatch" % label
        print l,
        print props


def GetBond ( l ) :
    global mrPar
    ts = l.split()
    if len (ts) == 11 :
        rtype, a1, a2, order, aromatic, stereo, ordinal, dist, esd, exception, descr = ts
        dist = float ( dist )
        esd = float ( esd )

        #if rtype == "ALA" :
        #    print rtype, a1, a2, dist, esd

        descr = descr.replace ( '"', "" )
        if not rtype in mrPar.bonds :
            mrPar.bonds[rtype] = []
        mrPar.bonds[rtype].append ( [a1, a2, dist, esd, descr] )


def GetLinkBond ( l ) :
    global mrPar
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

        mrPar.linkBonds[descr] = [a1, a2, dist, esd]

    else :
        #print "?"
        pass



def GetAngle ( l ) :
    global mrPar
    ts = l.split()

    rtype, a1, a2, a3, angle, esd, exception, descr = [None] * 8

    if len (ts) == 6 :
        rtype, a1, a2, a3, angle, esd = ts
    elif len (ts) == 8 :
        rtype, a1, a2, a3, angle, esd, exception, descr = ts
        descr = descr.replace ( '"', "" )
    else :
        print " -x- angle:", l,
        return

    angle = float ( angle ) * numpy.pi / 180.0
    esd = float ( esd ) * numpy.pi / 180.0

    #if rtype == "ALA" :
    #    print rtype, a1, a2, a3, angle, esd

    if not rtype in mrPar.angles :
        mrPar.angles[rtype] = []
    mrPar.angles[rtype].append ( [a1, a2, a3, angle, esd, descr] )


def GetLinkAngle ( l ) :
    global mrPar
    ts = l.split()
    if len (ts) > 5 :
        a1, a2, a3, angle, esd, descr = ts[0], ts[1], ts[2], ts[3], ts[4], ""
        angle = float ( angle ) * numpy.pi / 180.0
        esd = float ( esd ) * numpy.pi / 180.0

        for i in range ( 5, len(ts) ) :
            descr = descr + ts[i] + " "

        descr = descr[0:-1].replace( '"', '' )

        if not descr in mrPar.linkAngles :
            mrPar.linkAngles[descr] = {}

        mrPar.linkAngles[descr]["%s_%s_%s"%(a1,a2,a3)] = [angle, esd]
        #print " - angle", "%s_%s_%s : %.3f, %.3f"%(a1,a2,a3,angle,esd)
        #mrPar.linkAngles[descr].append ( [a1, a2, a3, angle, esd] )


def GetTorsion ( l ) :
    global mrPar
    ts = l.split()
    if len (ts) == 9 :
        # LYS      chi1     N      CA     CB     CG       180.0     15.0     3
        rtype, torId, a1, a2, a3, a4, angle, esd, period = ts
        angle = float ( angle ) * numpy.pi / 180.0
        esd = float ( esd ) * numpy.pi / 180.0
        period = int ( period )

        #if rtype == "ARG" or rtype == "ALA" :
        #    print rtype, torId, a1, a2, a3, a4, angle, esd, period

        if not rtype in mrPar.torsions :
            mrPar.torsions[rtype] = []
        mrPar.torsions[rtype].append ( [a1, a2, a3, a4, angle, esd, period] )



def GetPlane ( l ) :

    global mrPar
    ts = l.split()

    if len (ts) == 4 :
        # ARG      1      CD        0.02
        rtype, planeId, a1, esd = ts
        esd = float ( esd )

        #if rtype == "ASP" or rtype == "ASP" :
        #    print rtype, planeId, a1, esd

        if not rtype in mrPar.planes :
            mrPar.planes[rtype] = []
        mrPar.planes[rtype].append ( [rtype, planeId, a1, esd] )




def ReadParams ( resType ) :

    global mrPar

    if mrPar == None :
        mrPar = Params ()

    import os
    dir_path = os.path.dirname(os.path.realpath(__file__))
    #print dir_path

    li = 0

    paramFileName = ""
    if resType in protein3to1 or resType in nucleic3to1 :
        paramFileName = "standard_geometry.cif"
    else :
        paramFileName = resType.lower() + ".cif"


    fpath = os.path.join ( dir_path, "_param" )
    fpath = os.path.join ( fpath, paramFileName )

    if not os.path.isfile ( fpath ) :
        print " - could not find param file: %s" % fpath
        return

    print " - params from: %s" % fpath

    ctype = None

    fp = open ( fpath )
    while 1 :

        l = fp.readline(); li+=1;
        if not l : break

        if l[0:len("loop_")] == "loop_" :

            #print "\n%d:" % li,
            props, ctype = [], None

            while 1 :
                l = fp.readline(); li+=1;
                if not l : break

                if l[0:1] == "_" :
                    ctype, prop = l.split(".")
                    props.append ( prop.strip() )
                else :
                    #print "%s:%s" % (ctype, ",".join(props))
                    #print ""
                    break

        if l[0] == "#" :
            continue

        if not l : break

        if ctype == "_chem_comp_atom" :
            GetProps ( l, ctype, props, mrPar.atoms )

        elif ctype == "_chem_comp_bond" :
            #GetBond ( l, props )
            GetProps ( l, ctype, props, mrPar.bonds )

        elif ctype == "_chem_link_bond" :
            #GetLinkBond ( l, props )
            GetProps ( l, ctype, props, mrPar.linkBonds )

        elif ctype == "_chem_comp_angle" :
            #GetAngle ( l, props )
            GetProps ( l, ctype, props, mrPar.angles )

        elif ctype == "_chem_link_angle" :
            #GetLinkAngle ( l, props )
            GetProps ( l, ctype, props, mrPar.linkAngles )


        elif ctype == "_chem_comp_tor" :
            #GetTorsion ( l, props )
            GetProps ( l, ctype, props, mrPar.torsions )

        elif ctype == "_chem_comp_plane_atom" :
            #GetPlane ( l, props )
            GetProps ( l, ctype, props, mrPar.planes )


    fp.close()

    if 0 :
        print li, "lines"
        print "atoms ->", len(mrPar.atoms)
        print "bonds ->", len(mrPar.bonds)
        print "link bonds ->", len(mrPar.linkBonds)
        print "angles ->", len(mrPar.angles)
        print "link angles ->", len(mrPar.linkAngles)
        print "torsions ->", len(mrPar.torsions)
        print "planes ->", len(mrPar.planes)




def ResParams ( res, refPar ) :

    global mrPar

    if mrPar == None :
        mrPar = Params ()

    if res.type not in mrPar.atoms :
        ReadParams ( res.type )


    if res.type in mrPar.bonds :
        for prop in mrPar.bonds[res.type] :

            a1, a2 = prop["atom_id_1"], prop["atom_id_2"]
            dist, esd = prop["value_dist"], prop["value_dist_esd"]
            #a1, a2, dist, esd, descr = b

            useIt = False
            descr = ""

            if not "description" in prop :
                #print " - descr not found in bond"
                useIt = True

            else :
                descr = prop["description"]

                if res.type == "CYS" :
                    if descr == "."  :
                        useIt = True
                    elif descr == 'SH' :
                        useIt = True

                elif res.type == "HIS" :
                    if descr == "." :
                        useIt = True
                    elif descr == 'HISE' :
                        useIt = True

                elif res.type == "PRO" :
                    if descr == "." :
                        useIt = True
                    elif descr == 'trans' :
                        useIt = True

                elif res.type == "U" or res.type == "G" or res.type == "A" or res.type == "C" :
                    if descr == "." or descr == "s" :
                        useIt = True
                    elif descr == "C3'-endo" : #RNA
                        useIt = True
                    elif descr == "(H)C3'-end" : # ?
                        useIt = True

                else :
                    useIt = True

            if useIt :

                if not a1 in res.atomsMap :
                    #print " - %s.%d.%s - bond: %s-%s - %s ! a1 not found" % (res.type, res.id.position, res.id.chainId, a1, a2, descr)
                    continue

                if not a2 in res.atomsMap :
                    #print " - %s.%d.%s - bond: %s-%s - %s ! a2 not found" % (res.type, res.id.position, res.id.chainId, a1, a2, descr)
                    continue

                at1 = res.atomsMap[a1][0]
                at2 = res.atomsMap[a2][0]

                #print " - %s.%d.%s - bond: %s-%s - %s - %.3f %.3f" % (res.type, res.id.position, res.id.chainId, a1, a2, descr, dist, esd)
                refPar.bonds.append ( [at1, at2, dist, esd] )

    else :
        print " - no bonds for res %d.%s:%s" % (res.id.position, res.id.chainId, res.type)


    if res.type in mrPar.angles :
        for prop in mrPar.angles[res.type] :

            a1, a2, a3 = prop["atom_id_1"], prop["atom_id_2"], prop["atom_id_3"]
            angle, esd = prop["value_angle"], prop["value_angle_esd"]

            angle_rad = angle * numpy.pi / 180.0
            esd_rad = esd * numpy.pi / 180.0

            useIt = False
            descr = ""

            if not "description" in prop :
                #print " - descr not found in angle prop"
                useIt = True

            else :
                descr = prop["description"]

                if res.type == "CYS" :
                    if descr == "."  :
                        useIt = True
                    elif descr == "SH" :
                        useIt = True

                elif res.type == "HIS" :
                    if descr == "." :
                        useIt = True
                    elif descr == "HISE" :
                        useIt = True

                elif res.type == "PRO" :
                    if descr == "." :
                        useIt = True
                    elif descr == "trans" :
                        useIt = True

                elif res.type == "U" or res.type == "G" or res.type == "A" or res.type == "C" :
                    if descr == "." or descr == "s" :
                        useIt = True
                    elif descr == "C3'-endo" : #RNA
                        useIt = True
                    #elif descr == "(P)O5'-C5'" : # ?
                    #    useIt = True
                    elif descr == "small" : # ?
                        useIt = True


                else :
                    useIt = True

            if useIt :

                if not a1 in res.atomsMap :
                    #print " - %s.%d.%s - angle: %s-%s-%s - %s ! a1 not found" % (res.type, res.id.position, res.id.chainId, a1, a2, a3, descr)
                    continue

                if not a2 in res.atomsMap :
                    #print " - %s.%d.%s - angle: %s-%s-%s - %s ! a2 not found" % (res.type, res.id.position, res.id.chainId, a1, a2, a3, descr)
                    continue

                if not a3 in res.atomsMap :
                    #print " - %s.%d.%s - angle: %s-%s-%s - %s ! a3 not found" % (res.type, res.id.position, res.id.chainId, a1, a2, a3, descr)
                    continue

                at1 = res.atomsMap[a1][0]
                at2 = res.atomsMap[a2][0]
                at3 = res.atomsMap[a3][0]
                refPar.angles.append ( [at1, at2, at3, angle_rad, esd_rad] )

                #print " - %s.%d.%s - angle: %s-%s-%s - %s - %.3f %.3f" % (res.type, res.id.position, res.id.chainId, a1, a2, a3, descr, angle, esd)

    else :
        print " - angles for res %d.%s:%s not found" % (res.id.position, res.id.chainId, res.type)


    if res.type in mrPar.torsions :
        for P in mrPar.torsions[res.type] :

            a1, a2, a3, a4 = P["atom_id_1"], P["atom_id_2"], P["atom_id_3"], P["atom_id_4"]
            angle, esd, period = P["value_angle"], P["value_angle_esd"], P["period"]

            #a1, a2, a3, a4, angle, esd, period = a
            angle_rad = angle * numpy.pi / 180.0
            esd_rad = esd * numpy.pi / 180.0

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
                refPar.torsions.append ( [at1, at2, at3, at4, angle_rad, esd_rad, period] )
                #print " -t: ", at1.name, at2.name, at3.name, angle_rad, esd_rad, period
    else :
        #print " - torsions for res %d.%s:%s not found" % (res.id.position, res.id.chainId, res.type)
        pass


    if res.type in mrPar.planes :

        planeAtoms = []

        for P in mrPar.planes[res.type] :
            #print a
            planeId, a1, esd = P["plane_id"], P["atom_id"], P["dist_esd"]

            at1 = None
            if a1 in res.atomsMap :
                at1 = res.atomsMap[a1][0]

            if at1 :
                planeAtoms.append ( [at1, esd] )
                #if log : print " -p: ", At(at1), planeId, esd
            else :
                #print " - atom %s not found for plane, in res %d %s chain %s" % (a1, res.id.position, res.type, res.id.chainId)
                pass

        if len (planeAtoms) >= 3 :
            refPar.planes.append ( planeAtoms )
        #print " res %d.%s:%s - %d atoms in plane " % (res.id.position, res.id.chainId, res.type, len(planeAtoms) )

    #else :
    #    print " - res %d.%s:%s - no planes" % (res.id.position, res.id.chainId, res.type)



def CheckConnectDiS ( at1, refPar ) :

    P = None
    for linkI, lpar in mrPar.linkBonds.iteritems() :
        if lpar['link_id'] == "SG" and lpar["atom_id_1"] == "SG" and lpar["link_id"] == "CYS" :
            P = lpar

    if P :
        a1, a2, dist, esd = P['atom_id_1'], P['atom_id_2'], P['value_dist'], P['value_dist_esd']

        if not a1 in r1.atomsMap :
            print " - %s-%s [%s] link at %s-%d-%s - %s at not found" % (a1, a2, lpar['link_id'], r1.type, r1.id.position, r1.id.chainId, a1)
            return

        if not a2 in r2.atomsMap :
            print " - %s-%s [%s] link at %s-%d-%s - %s at not found" % (a1, a2, lpar['link_id'], r2.type, r2.id.position, r2.id.chainId, a2)
            return

        print " - diS link bond %s-%s [%s] %d.%s - %d.%s" % (a1, a2, lpar['link_id'], r1.id.position, r1.id.chainId, r2.id.position, r2.id.chainId)

        at1 = r1.atomsMap[a1][0]
        at2 = r2.atomsMap[a2][0]

        refPar.bonds.append ( [at1, at2, dist, esd] )


def ConnectRes (r1, r2, refPar) :

    for linkI, lpar in mrPar.linkBonds.iteritems() :

        P = None
        if r2.type in protein3to1 and r1.type in protein3to1 :
            if r2.type == "GLY" :
                if lpar['link_id'] == "GLY" :
                    P = lpar
            elif r2.type == "PRO" :
                if lpar['link_id'] == "PRO" :
                    P = lpar
            else :
                if lpar['link_id'] == "NOT GLY PRO" :
                    P = lpar

        elif r2.type in nucleic3to1 and r1.type in nucleic3to1 :
            if lpar['link_id'] == "NUC-ACID-ALL" :
                P = lpar



        if P :
            a1, a2, dist, esd = P['atom_id_1'], P['atom_id_2'], P['value_dist'], P['value_dist_esd']

            if not a1 in r1.atomsMap :
                print " - %s-%s [%s] link at %s-%d-%s - %s at not found" % (a1, a2, lpar['link_id'], r1.type, r1.id.position, r1.id.chainId, a1)
                continue

            if not a2 in r2.atomsMap :
                print " - %s-%s [%s] link at %s-%d-%s - %s at not found" % (a1, a2, lpar['link_id'], r2.type, r2.id.position, r2.id.chainId, a2)
                continue

            #print " - link bond %s-%s [%s] link at %s-%d-%s" % (a1, a2, lpar['link_id'], r1.type, r1.id.position, r1.id.chainId)

            at1 = r1.atomsMap[a1][0]
            at2 = r2.atomsMap[a2][0]

            refPar.bonds.append ( [at1, at2, dist, esd] )


    for linkI, lpar in mrPar.linkAngles.iteritems() :

        a1, a2, a3 = lpar['atom_id_1'], lpar['atom_id_2'], lpar['atom_id_3']

        P = None

        if r2.type in protein3to1 and r1.type in protein3to1 :
            if r2.type == "GLY" :
                if lpar['link_id'] == "GLY" :
                    P = lpar
            elif r2.type == "PRO" :
                isTrans = True
                if lpar['link_id'] == "PRO" :
                    P = lpar
                else :
                    if isTrans :
                        if lpar['link_id'] == "PRO-TRANS" :
                            P = lpar
                    else :
                        if lpar['link_id'] == "PRO-CIS" :
                            #P = lpar
                            #print " - adding PRO C-N link at %s-%d-%s" % (r1.type, r1.id.position, r1.id.chainId)
                            pass
            else :
                if lpar['link_id'] == "NOT GLY PRO" :
                    P = lpar
                if lpar['link_id'] == "NOT PRO GLY" :
                    P = lpar

        elif r2.type in nucleic1to3 and r1.type in nucleic1to3 :
            if lpar['link_id'] == "NA s" :
                P = lpar
            if lpar['link_id'] == "NA small" :
                P = lpar

        if P :

            angle, esd = P['value_angle'], P['value_angle_esd']

            angle_rad = angle * numpy.pi / 180.0
            esd_rad = esd * numpy.pi / 180.0

            rr1, rr2, rr3 = r1, None, r2
            if a1 == "CA" or a1 == "O" :
                rr2 = r1
            else :
                rr2 = r2

            if not a1 in rr1.atomsMap :
                print " - %s-%s-%s [%s] at %s-%d-%s ! a1 not found" % (a1, a2, a3, lpar['link_id'], rr1.type, rr1.id.position, rr1.id.chainId)
                continue

            if not a2 in rr2.atomsMap :
                print " - %s-%s-%s [%s] at %s-%d-%s ! a2 not found" % (a1, a2, a3, lpar['link_id'], rr2.type, rr2.id.position, rr2.id.chainId)
                continue

            if not a3 in rr3.atomsMap :
                print " - %s-%s-%s [%s] at %s-%d-%s ! a3 not found" % (a1, a2, a3, lpar['link_id'], rr3.type, rr3.id.position, rr3.id.chainId)
                continue

            #print " - link angle %s-%s-%s [%s] at %s-%d-%s - %.3f" % (a1, a2, a3, lpar['link_id'], r1.type, r1.id.position, r1.id.chainId, angle_rad*180.0/numpy.pi)

            at1 = rr1.atomsMap[a1][0]
            at2 = rr2.atomsMap[a2][0]
            at3 = rr3.atomsMap[a3][0]

            refPar.angles.append ( [at1, at2, at3, angle_rad, esd_rad]  )



refPar = None
refAtoms = None
refPos = None
refGrad = None
refAtM = None
gMaxDs = [0.0,0.0,0.0]
gMaxG = None


def RefStart ( ress, dmap ) :

    global refPar
    global refAtoms, refPos, refGrad, refAtM

    #print "Refining %d residues" % len(ress)

    atomsMap = {}
    molMap = {}
    resMap = {}

    for r in ress :
        molMap[r.molecule] = 1
        resMap[r] = 0
        for at in r.atoms :
            atomsMap[at] = 1


    bonds = []
    # resisdues mapped by mol_id, chain_id, and residue id
    mi_ci_ri = {}
    for mol in molMap.keys() :

        #for bond in mol.bonds :
        #    if bond.atoms[0] in atomsMap and bond.atoms[1] in atomsMap :
        #        bonds.append ( bond )

        if not mol.id in mi_ci_ri :
            mi_ci_ri[mol.id] = {}

        for res in mol.residues :
            if not res.id.chainId in mi_ci_ri[mol.id] :
                mi_ci_ri[mol.id][res.id.chainId] = {}
            if res.id.position in mi_ci_ri[mol.id][res.id.chainId] :
                print "multiple residues with position %d in %d.%s" % (res.id.position, mol.id, res.id.chainId)
                print " %s - %s" % (mi_ci_ri[mol.id][res.id.chainId][res.id.position].type, res.type)
            else :
                mi_ci_ri[mol.id][res.id.chainId][res.id.position] = res


    # stores if the next res from any residue has had parameters added yet
    #addedParForNextRes = {}

    # all parameters
    refPar = RefineParams ()

    for res in ress :

        mol = res.molecule

        for at in res.atoms :
            atomsMap[at] = 1
            at.M = 1.0

        ResParams ( res, refPar )

        #print " - %d.%s %s - %d bonds, %d angles" % (r.id.position, r.id.chainId, r.type, len(rBonds), len(rAngles))

        if res.type in protein3to1 or res.type in nucleic1to3 :

            if res.id.position-1 in mi_ci_ri[mol.id][res.id.chainId] :
                prevRes = mi_ci_ri[mol.id][res.id.chainId][res.id.position-1]

                if prevRes in resMap :
                    # connects from prevRes to this one
                    pass
                else :
                    ConnectRes ( prevRes, res, refPar )

                    for at in prevRes.atoms :
                        atomsMap[at] = 1
                        at.M = 0.0


            if res.id.position+1 in mi_ci_ri[mol.id][res.id.chainId] :
                nextRes = mi_ci_ri[mol.id][res.id.chainId][res.id.position+1]

                ConnectRes ( res, nextRes, refPar )

                if not nextRes in resMap :
                    for at in nextRes.atoms :
                        atomsMap[at] = 1
                        at.M = 0.0

        if 0 and res.type == "CYS" :
            CheckConnectDiS ()
            atO1 = res.atomsMap["OG"][0]
            if len(atO1.bonds) == 2 :
                atO2 = None
                for b in atO1.bonds :
                    if b.otherAtom(atO1).name == "OG" :
                        atO2 = b.otherAtom(atO1)
                if atO2 :
                    add = False
                    if atO2.residue in resMap :
                        # add from lowest to highest to avoid adding twice
                        if atO2.residue.id.position < atO1.residue.id.position :
                            add = True
                    else :
                        add = True
                    if add :
                        pass


        if res.type == "NAG" :
            C1 = res.atomsMap["C1"][0]
            for b in C1.bonds :

                bondedAt = b.otherAtom(C1)
                bondedRes = bondedAt.residue

                if bondedRes.type == "ASN" and bondedAt.name == "ND2" :
                    refPar.bonds.append ( [C1, bondedAt, 1.45, 0.02] )
                    CG = bondedRes.atomsMap["CG"][0]
                    refPar.angles.append ( [CG, bondedAt, C1, 124.7*numpy.pi/180.0, 5.0*numpy.pi/180.0] )
                    print " - NAG-ASN(%s) bond" % bondedAt.name
                    if not bondedRes in resMap :
                        for at in bondedRes.atoms :
                            atomsMap[at] = 1
                            at.M = 0.0

                elif bondedRes.type == "NAG" and bondedAt.name == "O4" :
                    refPar.bonds.append ( [C1, bondedAt, 1.433, 0.02] )
                    C4 = bondedRes.atomsMap["C4"][0]
                    refPar.angles.append ( [C4, bondedAt, C1, 118.567*numpy.pi/180.0, 5.0*numpy.pi/180.0] )
                    print " - NAG-NAG (%s) bond" % bondedAt.name
                    if not bondedRes in resMap :
                        for at in bondedRes.atoms :
                            atomsMap[at] = 1
                            at.M = 0.0

        if res.type == "BMA" :
            C1 = res.atomsMap["C1"][0]
            for b in C1.bonds :

                bondedAt = b.otherAtom(C1)
                bondedRes = bondedAt.residue

                if bondedRes.type == "NAG" and bondedAt.name == "O4" :
                    refPar.bonds.append ( [C1, bondedAt, 1.433, 0.02] )
                    C4 = bondedRes.atomsMap["C4"][0]
                    refPar.angles.append ( [C4, bondedAt, C1, 109.147*numpy.pi/180.0, 5.0*numpy.pi/180.0] )
                    print " - BMA-NAG (%s) bond" % bondedAt.name
                    if not bondedRes in resMap :
                        for at in bondedRes.atoms :
                            atomsMap[at] = 1
                            at.M = 0.0


        if res.type == "MAN" :
            C1 = res.atomsMap["C1"][0]
            for b in C1.bonds :

                bondedAt = b.otherAtom(C1)
                bondedRes = bondedAt.residue

                if bondedRes.type == "BMA" and bondedAt.name == "O6" :
                    refPar.bonds.append ( [C1, bondedAt, 1.425, 0.02] )
                    C6 = bondedRes.atomsMap["C6"][0]
                    refPar.angles.append ( [C6, bondedAt, C1, 115.695*numpy.pi/180.0, 5.0*numpy.pi/180.0] )
                    print " - MAN-BMA (%s) bond" % bondedAt.name
                    if not bondedRes in resMap :
                        for at in bondedRes.atoms :
                            atomsMap[at] = 1
                            at.M = 0.0

                elif bondedRes.type == "BMA" and bondedAt.name == "O3" :
                    refPar.bonds.append ( [C1, bondedAt, 1.475, 0.02] )
                    C3 = bondedRes.atomsMap["C3"][0]
                    refPar.angles.append ( [C3, bondedAt, C1, 110.731*numpy.pi/180.0, 5.0*numpy.pi/180.0] )
                    print " - MAN-BMA (%s) bond" % bondedAt.name
                    if not bondedRes in resMap :
                        for at in bondedRes.atoms :
                            atomsMap[at] = 1
                            at.M = 0.0


    if dmap and not hasattr ( dmap, 'maxd' ) :
        M = dmap.data.full_matrix()
        avg, std = numpy.average(M), numpy.std(M)
        #maxM = numpy.max(M)
        #minM = numpy.min(M)
        dmap.maxd = min ( avg+std*10.0, numpy.max(M) )
        dmap.mind = max ( avg-std*1.0, numpy.min(M) )

    refAtoms = atomsMap.keys()
    refPos = numpy.zeros ( [len(refAtoms),3] )
    refGrad = numpy.zeros ( [len(refAtoms),3] )
    refAtM = numpy.zeros ( [len(refAtoms),3] )

    #global gAtId
    atId = GetLastId () + 1
    for i, at in enumerate(refAtoms) :

        at.i = i
        #refPos[i] = at.coord()
        refAtM[i] = [at.M, at.M, at.M]

        if not hasattr ( at, 'coords' ) :
            at.coords = {}
        #at.coordAt = len( at.coords )
        at.coords[atId] = at.coord()

    print "Coord id is %d" % atId
    global gAtId
    gAtId = atId

    #global refAtomsPos0
    #refAtomsPos0 = _multiscale.get_atom_coordinates ( refAtoms, transformed = False )


def GetLastId () :
    lastId = -1
    for m in chimera.openModels.list(modelTypes = [chimera.Molecule]) :
        for at in m.atoms :
            if hasattr ( at, 'coords' ) :
                lastId = max ( lastId, max ( at.coords.keys() ) )
    return lastId

gAtId = None


def RefBack () :

    global gAtId
    if gAtId == None :
        gAtId = GetLastId ()

    if gAtId >= 0 :
        for m in chimera.openModels.list(modelTypes = [chimera.Molecule]) :
            for at in m.atoms :
                if hasattr ( at, 'coords' ) :
                    if gAtId in at.coords :
                        at.setCoord ( at.coords[gAtId] )

        print "Set coord id %d" % gAtId
        gAtId = gAtId - 1




import threading
import time
import Queue

class RefThread ( threading.Thread ):
    def __init__(self, outQueue, inQueue, dmap, mapF):
        threading.Thread.__init__(self)
        self.dmap = dmap
        self.mapF = mapF
        self.inQueue = inQueue
        self.outQueue = outQueue

    def run(self):
        #print "Thread start"

        while 1 :

            try:
                msg = self.inQueue.get(0)

            except Queue.Empty:
                time.sleep(0.1)
                #print "Thread done"
                continue

            if msg == "stop" :
                print " - thread msg:", msg
                break
            elif msg == "go" :
                startt = time.time()
                RefStep ( self.dmap, self.mapF  )
                dur = time.time() - startt

                e = RefE ( self.dmap )
                msg2 = "%.3fs / " % dur + e
                self.outQueue.put( msg2 )
            else :
                print " - thread msg: ? ", msg
                break



def RefStep ( dmap, mapF  ) :

    global refPar
    global refAtoms, refPos, refGrad, refAtM

    #mapF = float ( mapFStr.get() )
    #print mapF
    #if dmap :
    #    print " - in map %s" % dmap.name
    #mapF = 0.1

    for at in refAtoms :
        #at.G = chimera.Vector(0,0,0)
        #-at.P = at.coord()
        #at.P = chimera.Point(0,0,0)
        #-at.G = chimera.Vector(0,0,0)
        refPos[at.i] = at.coord()
        refGrad[at.i] = [0,0,0]


    #if log : print ""

    if 1 :
        for b in refPar.bonds :
            at1, at2, D, esd = b
            BondG ( at1, at2, D, esd, F=0.1 )

    #if log : print ""

    if 1 :
        for a in refPar.angles :
            at1, at2, at3, A, esd = a
            AngleG ( at1, at2, at3, A, esd, F=0.1 )
            #break

    if 1 :
        for t in refPar.torsions :
            at1, at2, at3, at4, angle, esd, period = t

            #if at1.name == "N" :
                #print "tor: %s, %s, %s, %s - angle %.3f, per %d" % ( At(at1), at2.name, at3.name, at4.name, angle*180.0/numpy.pi, period )
                #RotDih ( at1, at2, at3, at4, angle, esd, period, dmap )

            TorsionG ( at1, at2, at3, at4, angle, esd, period, F=0.001 )

    if 1 :
        for planeAtoms in refPar.planes :
            PlaneG ( planeAtoms, F=0.1 )
            pass

    if dmap :
        for at in refAtoms :
            MapG ( at, dmap, F=mapF )
            #break
            pass


    global gMaxDs, gMaxG

    maxD = 0.0

    #for at in refAtoms :
    #    #at.setCoord ( at.coord() + at.G )
    #    if not hasattr ( at, 'M' ) :
    #        print " - xMx - at %s" % At(at)
    #    else :
    #        #-gl = at.G.length
    #        #-if gl > gMaxG :
    #        #-    gMaxG = gl
    #        #-at.P = at.P + at.G * at.M
    #        #at.setCoord ( at.P )
    #        pass
    #    #pass

    gMaxDs[0] = gMaxDs[1]
    gMaxDs[1] = gMaxDs[2]
    gMaxDs[2] = gMaxG

    gMaxG = numpy.sqrt ( max ( numpy.sum ( refGrad*refGrad, axis=1 ) ) )

    scaleF = 1.0
    if gMaxG > 10.0 :
        scaleF = (10.0 / gMaxG)
        refGrad = refGrad * scaleF
        #gMaxG = numpy.sqrt ( max ( numpy.sum ( refGrad*refGrad, axis=1 ) ) )

    refPos = refPos + (refGrad * refAtM)



def RefPut () :
    global refAtoms, refPos
    for at in refAtoms :
        p = refPos[at.i]
        at.setCoord ( chimera.Point(*p) )



def RefE ( dmap ) :

    global refPar
    global refAtoms

    dB = 0.0
    if 1 :
        for b in refPar.bonds :
            at1, at2, D, esd = b
            dB += BondG ( at1, at2, D, esd, F=0 )

            #if log : print " - b %s _ %s %.3f (%.3f)" % (At(at1), At(at2), D, d)

    #if log : print ""

    dA = 0.0
    if 1 :
        for a in refPar.angles :

            at1, at2, at3, A, esd = a
            dA += AngleG ( at1, at2, at3, A, esd, F=0 )

            #if log : print " - a %s _ %s _ %s %.3f (%.3f)" % (At(at1), At(at2), At(at3), A, d)

    dTor = 0.0
    if 1 :
        for t in refPar.torsions :
            at1, at2, at3, at4, angle, esd, period = t

            #if at1.name == "N" :
            #    print "tor: %s, %s, %s, %s - angle %.3f, per %d" % ( At(at1), at2.name, at3.name, at4.name, angle*180.0/numpy.pi, period )
            dTor += TorsionG ( at1, at2, at3, at4, angle, esd, period, F=0 )

    dPlanes = 0.0
    if 1 :
        for planeAtoms in refPar.planes :
            dPlanes += PlaneG ( planeAtoms, F=0 )


    dMap = 0.0
    if dmap  :
        for at in refAtoms :
            dMap += MapG ( at, dmap, F=0 )

    #if dmap :
    #    RefDisp ( refAtoms, dmap )

    #if log : print ""
    #maxG = max ( numpy.sum ( refGrad*refGrad, axis=1 ) )
    global gMaxG

    estr =  "[%d] atoms, [%d] bonds:%.3f, [%d] angles:%.3f, [%d] torsions:%.3f, [%d] planes:%.3f, map:%.4f maxG:%.5f" % ( len(refAtoms),
                len(refPar.bonds), dB, len(refPar.angles), dA, len(refPar.torsions), dTor, len(refPar.planes), dPlanes, dMap, gMaxG )

    return estr




def MapG ( at, dmap, F=0.1 ) :

    global refAtoms, refPos, refGrad

    #print at.P, type(at.P)
    #tp = at.P
    ##tp = at.molecule.openState.xform.apply ( at.P )
    ##tp = dmap.openState.xform.inverse().apply ( tp )
    #P = tp.data()

    P = refPos[at.i]

    if F <= 0.0 :
        return dmap.interpolated_values ( [P], at.molecule.openState.xform ) [0]

    dx, dy, dz = dmap.data.step
    pts = numpy.array ( [P] * 6 )

    #print pts

    #dx *= 0.1

    pts[0][0] -= dx; pts[1][0] += dx
    pts[2][1] -= dy; pts[3][1] += dy
    pts[4][2] -= dz; pts[5][2] += dz

    #print pts, type(pts)

    dmap.ctr = numpy.array(dmap.data.origin) + numpy.array(dmap.data.size)/2.0*numpy.array(dmap.data.step)


    vs = dmap.interpolated_values ( pts, at.molecule.openState.xform )

    #print mapvs, type(mapvs)

    #vs = (vs - dmap.mind)/(dmap.maxd-dmap.mind)

    #print mapvs, type(mapvs)

    dx2 = 1.0/(dx*2.0)

    G = chimera.Vector ( (vs[1] - vs[0])*dx2, (vs[3] - vs[2])*dx2, (vs[5] - vs[4])*dx2  )
    #print G

    #G = dmap.ctr - numpy.array( tp )
    #G = chimera.Vector ( *G )

    #G = dmap.openState.xform.apply ( G )
    #G = at.molecule.openState.xform.inverse().apply ( G )
    #print G

    #at.dg = G
    #at.G += G * F
    refGrad[at.i] += G * F




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


        if hasattr ( at, 'toPlaneV' ) :

            #axes.AddArrow4 ( P, G, G.length, clr=(.5,.8,.5,1), rad=0.1, mol=smod, hrad=0.15, hlen=0.05 )
            if at.toPlaneV.length > 1e-3 :
                axes.AddArrow4 ( at.P, at.toPlaneV, at.toPlaneV.length, clr=(.5,.8,.5,1), rad=0.1, mol=smod, hrad=0.15, hlen=0.05 )
                #axes.AddArrow4 ( at.com, at.toPlaneV, at.toPlaneV.length, clr=(.7,.2,.5,1), rad=0.1, mol=smod, hrad=0.15, hlen=0.05 )



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


def dihap ( a1, a2, a3, a4 ) :
    #n1 = vnorm ( a1.coord(), a2.coord(), a3.coord() )
    #n2 = vnorm ( a2.coord(), a3.coord(), a4.coord() )
    #return numpy.arccos ( n2 * n1 * -1.0 ) * 180.0 / numpy.pi

    # http://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
    b1 = a2 - a1
    b2 = a3 - a2
    b3 = a4 - a3

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

    #v = at1.P - at2.P

    global refAtoms, refPos, refGrad
    pv = refPos[at1.i] - refPos[at2.i]
    v = chimera.Vector( *pv )

    d = D - v.length

    #if log : print " - b0 %s-%s %.3f/%.3f" % (at1.name, at2.name, D, d)

    if F <= 0.0 :
        return d * d


    if v.length < 0.01 :
        v = chimera.Vector ( random.random(), random.random(),random.random() )

    v.normalize()
    #at1.G += v * d * F
    #at2.G -= v * d * F
    #at1.P = at1.P + v * d * 0.1 * at1.M
    #at2.P = at2.P - v * d * 0.1 * at2.M
    refGrad[at1.i] += v * d * F
    refGrad[at2.i] -= v * d * F



def AngleG_ (at1, at2, at3, theta0, esd, F=0.001) :

    global refAtoms, refPos, refGrad
    pv = refPos[at1.i] - refPos[at2.i]

    r12 = refPos[at1.i] - refPos[at2.i]
    r32 = refPos[at3.i] - refPos[at2.i]

    r12 = chimera.Vector(*r12); d12 = r12.length
    r32 = chimera.Vector(*r32); d32 = r32.length

    #r12 = at1.P - at2.P; d12 = r12.length
    #r32 = at3.P - at2.P; d32 = r32.length

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

    #rij = at1.P - at2.P
    #rik = at3.P - at2.P

    global refAtoms, refPos, refGrad

    rij = refPos[at1.i] - refPos[at2.i]; rij = chimera.Vector (*rij)
    rik = refPos[at3.i] - refPos[at2.i]; rik = chimera.Vector (*rik)

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

    #at1.G += fj * (ge * F)
    #at2.G += fi * (ge * F)
    #at3.G += fk * (ge * F)

    refGrad[at1.i] += fj * (ge * F)
    refGrad[at2.i] += fi * (ge * F)
    refGrad[at3.i] += fk * (ge * F)



def TorsionG ( at1, at2, at3, at4, angle, esd, period, F=0.01 ) :

	# Calculate the vectors between atoms
	#Vector3d pos0 = d->a1->pos, pos1 = d->a2->pos, pos2 = d->a3->pos, pos3 = d->a4->pos;
	#Vector3d r12 = pos0-pos1, r23 = pos1-pos2, r34 = pos2-pos3;

    #r12 = at1.P - at2.P;
    #r23 = at2.P - at3.P;
    #r34 = at3.P - at4.P

    global refAtoms, refPos, refGrad
    r12 = refPos[at1.i] - refPos[at2.i]; r12 = chimera.Vector (*r12)
    r23 = refPos[at2.i] - refPos[at3.i]; r23 = chimera.Vector (*r23)
    r34 = refPos[at3.i] - refPos[at4.i]; r34 = chimera.Vector (*r34)


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
    diff = 0

	# Calculate the energy
    if n > 0 :
    	# Periodicity is greater than 0, so use cos form
    	#K += k*(1+numpy.cos(n*phi - delta));
    	#K1 += -n*k*numpy.sin(n*phi - delta);

        d = n*phi - delta
    	K += k*(1+numpy.cos(d))
    	K1 += -n*k*numpy.sin(d)

        #minDiff = 1e9
        #for i in range (-(n+1), n+1) :
        #    if i != 0 :
        #        d = i*phi - delta
        #        ad = abs(d)
        #        if ad < minDiff :
        #            minDiff = ad
        #            diff = d

    else :
		# Periodicity is 0, so just use the harmonic form
        diff = phi - delta
        if diff < -M_PI :
            diff += TWOPI
        elif diff > M_PI :
            diff -= TWOPI;

        K += k*diff*diff;
        K1 += 2.0*k*diff;

    #print "tor: %s, %s, %s, %s - angle %.3f / %.3f, per %d" % ( At(at1), at2.name, at3.name, at4.name, angle, phi*180/numpy.pi, period )
    #print "%.3f\t%.3f" % (phi*180/numpy.pi, (-K+1.0)),
    #print " - tor %s %s %s %s - per %d, phi %.2f (%.2f), diff %.2f -- %.3f" % (At(at1), at2.name, at3.name, at4.name, period, phi*180.0/numpy.pi, angle*180.0/numpy.pi, diff*180.0/numpy.pi, K)
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
    #at1.G += f1
    #at2.G += f2 - f1
    #at3.G += f3 - f2
    #at4.G += -f3

    refGrad[at1.i] += f1
    refGrad[at2.i] += f2 - f1
    refGrad[at3.i] += f3 - f2
    refGrad[at4.i] += -f3


def PlaneG ( planeAtoms, F=0.1 ) :

    global refAtoms, refPos, refGrad

    points = numpy.zeros ( (len(planeAtoms), 3) )
    i = 0
    for at, esd in planeAtoms :
        #points[i] = at.P
        points[i] = refPos[at.i]
        i += 1

    if 0 :
        com = numpy.sum(points, axis=0) / len(points)
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

        #if F > 0 :
            #print "ctr:", C
        #    print "U:", U
        #    print "S:", S
            #print "V:", V

        ni = numpy.argmin ( numpy.abs(S) )
        #N = numpy.array ( [U[0,ni], U[1,ni], U[2,ni]] )
        #N = numpy.array ( [V[ni,0], V[ni,1], V[ni,2]] )
        N = numpy.array ( [U[ni,0], U[ni,1], U[ni,2]] )

    # barycenter of the points
    # compute centered coordinates
    com = points.sum(axis=0) / points.shape[0]

    # run SVD
    u, s, vh = numpy.linalg.svd(points - com)

    #print com
    #print u
    #print s
    #print vh

    # unitary normal vector
    N = vh[2, :]

    #if F > 0 : print "N:", N

    i = 0
    sumL = 0.0
    for at, esd in planeAtoms :
        v = com - points[i]
        dot = numpy.dot ( v, N )
        #if F > 0 : print "%s -- %d : %.3f" % (at.name, i, dot)
        sumL += abs ( dot )
        toPlane = N * (dot * F)
        #at.G += chimera.Vector ( *toPlane )
        refGrad[at.i] += toPlane
        #at.toPlaneV = chimera.Vector ( *(N*dot) ) * 10.0
        #at.com = chimera.Point ( *com )
        i += 1

        #if hasattr ( at, 'toPlaneV' ) : del at.toPlaneV
        #if hasattr ( at, 'com' ) : del at.com

    #if F > 0 : print " - sumL: %.3f" % sumL

    if F <= 0 :
        return sumL






def AtomsNearRegs ( regs, maxD=4.0 ) :

    mols = []
    #print "atoms near in:"
    for m in chimera.openModels.list() :
        if type(m) == chimera.Molecule and m.display == True :
            mols.append ( m )

    import grid
    reload(grid)

    agrid = grid.Grid ()
    agrid.FromMols ( mols, maxD )

    n_regs = []

    points = regs[0].points().astype ( numpy.float32 )
    for r in regs[1:] :
        npoints = r.points().astype ( numpy.float32 )
        points = numpy.concatenate ( [points, npoints], axis=0 )

    segMap = regs[0].segmentation.seg_map
    _contour.affine_transform_vertices ( points, segMap.data.ijk_to_xyz_transform )
    _contour.affine_transform_vertices ( points, Matrix.xform_matrix( segMap.openState.xform ) )

    atMap = {}
    #resMap = {}
    for p in points :
        cp = chimera.Point(p[0],p[1],p[2])
        nearAts = agrid.AtsNearPt ( cp )
        for at in nearAts :
            atMap[at] = 1
            #resMap[at.residue] = 1

    return atMap.keys()




def AtomsNearAtoms ( atoms, maxD=4.0 ) :

    mols = {}
    #print "atoms near in:"
    #for m in chimera.openModels.list() :
    #    if type(m) == chimera.Molecule and m.display == True :
    #        mols.append ( m )

    for at in atoms :
        mols[at.molecule] = 1

    mols = mols.keys()

    import grid
    reload(grid)

    agrid = grid.Grid ()
    agrid.FromMols ( mols, maxD )

    atMap = {}
    #resMap = {}
    for at in atoms :
        nearAts = agrid.AtsNearPt ( at.xformCoord() )
        for at in nearAts :
            atMap[at] = 1
            #resMap[at.residue] = 1

    return atMap.keys()



def AddMol ( molName, selAt, inMap, regs, toMol=None, toChainId=None ) :

    if molName.lower() == "nag" :
        print " - adding nag"
        AddNAG ( selAt, inMap, regs )

    elif molName.lower() == "bma" :
        print " - adding bma"
        AddBMA ( selAt, inMap, regs )

    elif molName.lower() == "man" :
        print " - adding man"
        AddMAN ( selAt, inMap, regs )


    else :
        print " - adding %s" % molName.lower()
        AddMol2 ( molName.lower(), inMap, regs, toMol, toChainId )




def AddMol2 ( molName, inMap, regs, toMol, toChainId ) :

    fname = "/Users/greg/Dropbox/_mol/Segger/_param/%s.pdb" % molName

    from os import path
    if not path.isfile(fname) :
        print " - did not find %s" % fname
        return

    nmol = chimera.PDBio().readPDBfile ( fname )[0]
    print " - read %s - %d atoms - %d res" % ( nmol.name, len(nmol.atoms), len(nmol.residues) )
    addRes = nmol.residues[0]


    from axes import prAxes
    regsPoints = RegsPtsInMol ( regs, toMol )
    regsC, regsU, regsS, regsV = prAxes ( regsPoints )
    #print regsC, regsU

    molPoints = _multiscale.get_atom_coordinates ( addRes.atoms, transformed = False )
    molC, molU, molS, molV = prAxes ( molPoints )
    #print molC, molU

    import qscores
    reload ( qscores )

    xfs = uniform_rota_xfs ( 64 )
    score_xfs = []
    for rxf in xfs :

        xf = chimera.Xform.translation ( molC * -1.0 )
        xf.premultiply ( rxf )
        xf.premultiply ( chimera.Xform.translation ( regsC  ) )

        #nres = AddResToMol ( addRes, toMol, toChainId, xf, withoutAtoms=[] )

        cc0, cc1, xfm = FitAtomsXf ( addRes.atoms, xf, inMap )
        xf.premultiply ( xfm )

        takeIt = True
        rv, ang = xf.getRotation(); tr = xf.getTranslation()
        for cc1_, cc0_, xf0 in score_xfs :
            # todo: use rmsd instead!
            rv0, ang0 = xf0.getRotation(); tr0 = xf0.getTranslation()
            rd, rang, td = numpy.arccos(rv*rv0)*180.0/numpy.pi, abs(ang0-ang), (tr0-tr).length
            if rd < 5.0 and rang < 5.0 and td < 3.0 :
                takeIt = False

        if takeIt :
            molg = qscores.MyMolMapX2 ( addRes.atoms, 2.0, inMap.data.step[0], xf )
            fpoints, fpoint_weights = qscores.fit_points_g ( molg, 1e-2 )
            map_values = inMap.interpolated_values ( fpoints, toMol.openState.xform )
            #print map_values
            olap, CC, CCm = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
            print " - taking %.3f -> %.3f" % (cc1, CC)

            score_xfs.append ( [CC, cc0, xf] )


    score_xfs.sort ( reverse=True, key=lambda x: x[0] )
    print "%d unique" % len(score_xfs)

    #topCC, topCC0, xf = score_xfs[0]

    if 1 :
        for cc1, cc0, xf in score_xfs :
            nres = AddResToMol ( addRes, toMol, toChainId, xf, withoutAtoms=[] )
            nres.scoreCC = cc1
            print " - added %s - %d.%s -- %.4f" % (nres.type, nres.id.position, nres.id.chainId, cc1)

    else :
        q_xfs = []
        import qscores

        ats = [at for at in toMol.atoms if not at.element.name == "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(toMol.atoms) )
        pts = points.tolist()
        from CGLutil.AdaptiveTree import AdaptiveTree
        allPtTree = AdaptiveTree ( pts, pts, 1.0)
        #allAtTree = None
        minD, maxD = qscores.MinMaxD ( inMap )
        print " - minD %.3f, maxD %.3f" % (minD, maxD)


        for cc1, cc0, xf in score_xfs :
            #nres = AddResToMol ( addRes, toMol, toChainId, xf, withoutAtoms=[] )
            #nres.scoreCC = cc1
            #print " - added %s - %d.%s -- %.4f" % (nres.type, nres.id.position, nres.id.chainId, cc1)
            qavg, N = 0.0, 0.0
            for at in addRes.atoms :
                xfp = xf.apply ( at.coord() )
                #xfp = toMol.openState.xform.apply ( xfp )
                Q = qscores.QscorePt2 ( xfp, toMol.openState.xform, inMap, 0.6, allPtTree=allPtTree, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
                qavg += Q; N += 1.0
            Q = qavg / N

            q_xfs.append ( [Q, cc1, xf] )
            print " - cc %.3f -> Q %.3f" % (cc1, Q)

        q_xfs.sort ( reverse=True, key=lambda x: x[0] )

        for Q, cc1, xf in q_xfs :
            nres = AddResToMol ( addRes, toMol, toChainId, xf, withoutAtoms=[] )
            nres.scoreCC = cc1
            nres.scoreQ = Q
            print " - added %s - %d.%s -- cc %.3f, Q %.3f" % (nres.type, nres.id.position, nres.id.chainId, cc1, Q)




def SegFitRes ( res, inMap, regs, useAts=None ) :

    from axes import prAxes
    regsPoints = RegsPtsInMol ( regs, res.molecule )
    regsC, regsU, regsS, regsV = prAxes ( regsPoints )
    #print regsC, regsU

    ats = res.atoms if useAts == None else useAts
    points = _multiscale.get_atom_coordinates ( ats, transformed = False )
    molC, molU, molS, molV = prAxes ( points )
    #print molC, molU

    xfs = uniform_rota_xfs ( 64 )
    score_xfs = []
    for rxf in xfs :

        xf = chimera.Xform.translation ( molC * -1.0 )
        xf.premultiply ( rxf )
        xf.premultiply ( chimera.Xform.translation ( regsC  ) )

        #nres = AddResToMol ( addRes, toMol, toChainId, xf, withoutAtoms=[] )

        cc0, cc1, xfm = FitAtomsXf ( ats, xf, inMap )
        xf.premultiply ( xfm )

        rv, ang = xf.getRotation(); tr = xf.getTranslation()

        takeIt = True

        for cc1, cc0, xf0 in score_xfs :
            rv0, ang0 = xf0.getRotation(); tr0 = xf0.getTranslation()
            rd, rang, td = numpy.arccos(rv*rv0)*180.0/numpy.pi, abs(ang0-ang), (tr0-tr).length
            if rd < 5.0 and rang < 5.0 and td < 1.0 :
                takeIt = False

        if takeIt :
            score_xfs.append ( [cc1, cc0, xf] )


    score_xfs.sort ( reverse=True, key=lambda x: x[0] )
    print "%d unique" % len(score_xfs)

    #topCC, topCC0, xf = score_xfs[0]

    for cc1, cc0, xf in score_xfs :
        nres = AddResToMol ( res, res.molecule, res.id.chainId, xf, withoutAtoms=[] )
        nres.scoreCC = cc1
        print " - added %s - %d.%s -- %.4f" % (res.type, nres.id.position, nres.id.chainId, cc1)




def TorFit0 ( res, inMap ) :

    tors = FindTors ( [res] )
    print "%d tors" % len(tors)

    for tor in tors[1:2] :

        bond, ats1, ats2 = tor
        p2, p1 = bond.atoms[1].coord(), bond.atoms[0].coord()
        v = p2 - p1


        for i in range (180) :
            v.normalize()

            xf1 = chimera.Xform.translation ( p1.toVector() * -1 )
            xf1.premultiply ( chimera.Xform.rotation ( v, 1.0 ) )
            xf1.premultiply ( chimera.Xform.translation ( p1.toVector() ) )

            for at in ats1 :
                at.setCoord ( xf1.apply ( at.coord() ) )

            xf2 = chimera.Xform.translation ( p1.toVector() * -1 )
            xf2.premultiply ( chimera.Xform.rotation ( v, -1.0 ) )
            xf2.premultiply ( chimera.Xform.translation ( p1.toVector() ) )

            for at in ats2 :
                at.setCoord ( xf2.apply ( at.coord() ) )

            print ".",


        break




def TorFitGrads ( ress, inMap, useAtoms=None ) :

    tors = FindTors ( ress )
    print "%d tors" % len(tors)

    ressAtoms = []
    for r in ress :
        ressAtoms = ressAtoms + r.atoms

    conAts = ConAts ( ress )

    import grid
    reload(grid)

    if 0 :
        ats = res.atoms
        #ats = [at for at in atoms if not at.element.name == "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        #_contour.affine_transform_vertices ( fpoints,  Matrix.xform_matrix(xf0) )
        #_contour.affine_transform_vertices ( fpoints, inMap.openState.xform  )
        xyz_to_ijk_tf = inMap.data.xyz_to_ijk_transform
        matrix = inMap.data.full_matrix()
        from VolumeData import interpolate_volume_gradient
        gradients, outside = interpolate_volume_gradient(points, xyz_to_ijk_tf, matrix, 'linear')
        #print gradients
        for i, at in enumerate(ats) :
            at.i = i
            #print at.name, gradients[i]


    if useAtoms == None :
        useAtoms = ressAtoms
        print "using all %d atoms" % len(useAtoms)
    else :
        print "using %d atoms" % len(useAtoms)


    tf, rmat = inMap.data.xyz_to_ijk_transform, inMap.data.full_matrix()
    points = _multiscale.get_atom_coordinates ( useAtoms, transformed = False )
    last_avg = numpy.average ( VolumeData.interpolate_volume_data ( points, tf, rmat )[0] )
    print "%.5f -> " % last_avg,

    ijk_step_size_max = 0.5
    ijk_step_size_min = 0.01
    ijk_step_size = ijk_step_size_max

    agrid = None
    #agrid = grid.Grid ()
    #agrid.FromAtoms ( ressAtoms, 3.0 )

    for i in range ( 100 ) :
        for tor in tors :
            #bond, ats1, ats2 = tor
            bond, ats1, ats2 = tor
            p2, p1 = bond.atoms[1].coord(), bond.atoms[0].coord()
            v = p2 - p1
            AtTorques2 (bond, ats1, agrid, conAts, inMap.data, ijk_step_size)

            if 1 :
                AtTorques2 (bond, ats2, agrid, conAts, inMap.data, ijk_step_size)

            #break
        #break

        if 1 :
            cc, xf = FitAtoms ( useAtoms, inMap )
            for res in ress :
                for at in res.atoms :
                    at.setCoord ( xf.apply ( at.coord() ) )

        points = _multiscale.get_atom_coordinates ( useAtoms, transformed = False )
        avg1 = numpy.average ( VolumeData.interpolate_volume_data ( points, tf, rmat )[0] )
        print "%.3f" % avg1,

        if avg1 < last_avg :
            ijk_step_size = ijk_step_size / 2.0
            if ijk_step_size < ijk_step_size_min :
                print " ->| %.5f" % avg1
                #print " - reached min step size, stopping"
                break

        last_avg = avg1




def TorFitGrads_2sided ( ress, inMap, useAtoms=None ) :

    tors = FindTors ( ress )
    print "%d tors" % len(tors)
    ressAtoms = []
    for r in ress :
        ressAtoms = ressAtoms + r.atoms

    conAts = ConAts (res)

    import grid
    reload(grid)

    if 0 :
        ats = res.atoms
        #ats = [at for at in atoms if not at.element.name == "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        #_contour.affine_transform_vertices ( fpoints,  Matrix.xform_matrix(xf0) )
        #_contour.affine_transform_vertices ( fpoints, inMap.openState.xform  )
        xyz_to_ijk_tf = inMap.data.xyz_to_ijk_transform
        matrix = inMap.data.full_matrix()
        from VolumeData import interpolate_volume_gradient
        gradients, outside = interpolate_volume_gradient(points, xyz_to_ijk_tf, matrix, 'linear')
        #print gradients
        for i, at in enumerate(ats) :
            at.i = i
            #print at.name, gradients[i]


    if useAtoms == None :
        useAtoms = res.atoms
        print "using all %d atoms" % len(useAtoms)
    else :
        print "using %d atoms" % len(useAtoms)


    tf, rmat = inMap.data.xyz_to_ijk_transform, inMap.data.full_matrix()
    points = _multiscale.get_atom_coordinates ( useAtoms, transformed = False )
    last_avg = numpy.average ( VolumeData.interpolate_volume_data ( points, tf, rmat )[0] )
    print "%.5f -> " % last_avg,

    ijk_step_size_max = 0.5
    ijk_step_size_min = 0.01
    ijk_step_size = ijk_step_size_max

    for i in range ( 100 ) :

        agrid = grid.Grid ()
        agrid.FromAtoms ( res.atoms, 3.0 )

        for tor in tors :

            bond, ats1, ats2 = tor
            p2, p1 = bond.atoms[1].coord(), bond.atoms[0].coord()
            v = p2 - p1

            AtTorques2 (bond, ats1, agrid, conAts, inMap.data, ijk_step_size)
            AtTorques2 (bond, ats2, agrid, conAts, inMap.data, ijk_step_size)

            #break

        #break

        cc, xf = FitAtoms ( useAtoms, inMap )

        for at in res.atoms :
            at.setCoord ( xf.apply ( at.coord() ) )

        points = _multiscale.get_atom_coordinates ( useAtoms, transformed = False )
        avg1 = numpy.average ( VolumeData.interpolate_volume_data ( points, tf, rmat )[0] )
        print "%.3f" % avg1,

        if avg1 < last_avg :
            ijk_step_size = ijk_step_size / 2.0
            if ijk_step_size < ijk_step_size_min :
                print " -> %.5f" % avg1
                #print " - reached min step size, stopping"
                break

        last_avg = avg1




def TorFitGradsBonds ( selBonds, inMap, useAtoms=None ) :

    #tors = FindTors ( res )
    #print "%d tors" % len(tors)

    ress = []

    res = selBonds[0].atoms[0].residue
    mol = res.molecule
    SetBBAts ( mol )

    rmap = {}
    #for r in mol.residues :
    #    rmap[r.id.chainId + "%d"%r.id.position] = r
    for b in selBonds :
        rmap[b.atoms[0].residue] = 1
        rmap[b.atoms[1].residue] = 1

    if res.type in protein3to1 or res.type in nucleic1to3 :
        for r in res.molecule.residues :
            if r.id.chainId == res.id.chainId :
                ress.append ( r )
    else :
        ress = [res]

    tors = FindTors ( ress, selBonds )

    selTors = []
    for tor in tors :
        bond, ats1, ats2 = tor
        if bond in selBonds :
            selTors.append ( tor )

    print " - %d total tors, using %d" % ( len(tors), len(selTors) )

    if len(selTors) == 0 :
        return None


    scoreAtoms = []
    for r in rmap.keys() :
        #fitAtoms.extend ( r.atoms )
        scoreAtoms.extend ( r.bbAtoms )

    allAtoms = []
    for r in rmap.keys() :
        allAtoms.extend ( r.atoms )


    print " - %d atoms in %d res, score %d atoms - in %s" % (len(allAtoms), len(rmap.keys()), len(scoreAtoms), inMap.name)



    import grid
    reload(grid)


    if 0 :
        ats = res.atoms
        #ats = [at for at in atoms if not at.element.name == "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        #_contour.affine_transform_vertices ( fpoints,  Matrix.xform_matrix(xf0) )
        #_contour.affine_transform_vertices ( fpoints, inMap.openState.xform  )
        xyz_to_ijk_tf = inMap.data.xyz_to_ijk_transform
        matrix = inMap.data.full_matrix()
        from VolumeData import interpolate_volume_gradient
        gradients, outside = interpolate_volume_gradient(points, xyz_to_ijk_tf, matrix, 'linear')
        #print gradients
        for i, at in enumerate(ats) :
            at.i = i
            #print at.name, gradients[i]


    if useAtoms == None :
        useAtoms = res.atoms
        print "using all %d atoms" % len(useAtoms)
    else :
        print "using %d atoms" % len(useAtoms)

    useAtoms = scoreAtoms

    conAts = ConAts2(useAtoms)

    tf, rmat = inMap.data.xyz_to_ijk_transform, inMap.data.full_matrix()
    points = _multiscale.get_atom_coordinates ( useAtoms, transformed = False )
    last_avg = numpy.average ( VolumeData.interpolate_volume_data ( points, tf, rmat )[0] )
    print "%.5f -> " % last_avg,

    ijk_step_size_max = 0.5
    ijk_step_size_min = 0.01
    ijk_step_size = ijk_step_size_max

    for i in range ( 100 ) :

        agrid = grid.Grid ()
        agrid.FromAtoms ( useAtoms, 3.0 )

        for tor in selTors :

            bond, ats1, ats2 = tor
            p2, p1 = bond.atoms[1].coord(), bond.atoms[0].coord()
            v = p2 - p1

            AtTorques2 (bond, ats1, agrid, conAts, inMap.data, ijk_step_size)
            #AtTorques2 (bond, ats2, agrid, conAts, inMap.data, ijk_step_size)

            #break

        #break

        if 0 :
            cc, xf = FitAtoms ( useAtoms, inMap )
            for at in useAtoms :
                at.setCoord ( xf.apply ( at.coord() ) )

        #cc = FitScore ( useAtoms, inMap )

        points = _multiscale.get_atom_coordinates ( useAtoms, transformed = False )
        avg1 = numpy.average ( VolumeData.interpolate_volume_data ( points, tf, rmat )[0] )
        print "%.3f" % avg1,

        if avg1 < last_avg :
            ijk_step_size = ijk_step_size / 2.0
            if ijk_step_size < ijk_step_size_min :
                print " -> %.5f" % avg1
                #print " - reached min step size, stopping"
                break

        last_avg = avg1





def AtTorques2 (bond, atoms, agrid, conAts, mdata, step_size=0.5 ) :

    at1, at2 = bond.atoms
    center = at1.coord().data()
    rotVec = at2.coord() - at1.coord()
    rotVec.normalize()
    rotVecAr = numpy.array ( rotVec.data() )

    m1M = mdata.xyz_to_ijk_transform
    rmat = mdata.full_matrix()

    apoints = _multiscale.get_atom_coordinates ( atoms, transformed = False )
    point_weights = None
    #point_weights = numpy.ones ( len(apoints), numpy.float32 )

    gradients = VolumeData.interpolate_volume_gradient(apoints, m1M, rmat)[0]
    #maxG = numpy.max ( gradients, axis=0 )
    maxG = numpy.sqrt ( max ( numpy.sum ( gradients*gradients, axis=1 ) ) )
    #print "mg_%.3f" % ( maxG )

    #print gradients

    if agrid != None :
        grads2 = numpy.zeros ( [len(atoms), 3] )
        for i, at in enumerate ( atoms ) :
            if not at in conAts :
                print " at %s" % at.name
                continue
            nearAts = agrid.AtsNearPtLocal ( at.coord() )
            for atn, v in nearAts :
                if not atn in conAts[at] :
                    print "%s -- %s - %.3f" % (at.name, atn.name, v.length)
                    #D = v.length
                    v.normalize()
                    v = v * maxG * 3.0
                    grads2[i] -= v.data()

        #print grads2

        #maxG2 = numpy.sqrt ( max ( numpy.sum ( grads2*grads2, axis=1 ) ) )
        #print "max grad2:", maxG2

        gradients += grads2

    torque_axis = _distances.torque(apoints, point_weights, gradients, center )
    #print "torque axis 0:", torque_axis

    dt = numpy.dot ( torque_axis, rotVec.data() )
    torque_axis = rotVecAr if dt > 0 else rotVecAr * -1.0

    na = Matrix.norm(torque_axis)
    if na == 0 :
        #torque_axis = (0,0,1)
        pass
    else :
        torque_axis /= na
        angle = angle_step(torque_axis, apoints, center, m1M, step_size)
        #print "torque axis:", torque_axis, "angle:", angle
        move_tf = Matrix.rotation_transform ( torque_axis, angle, center )
        xf =  Matrix.chimera_xform ( move_tf )

        for at in atoms :
            at.setCoord ( xf.apply (at.coord()) )




def AtTorques1 (bond, atoms, mdata, step_size=0.5 ) :


    at1, at2 = bond.atoms
    center = at1.coord().data()
    rotVec = at2.coord() - at1.coord()
    rotVec.normalize()
    rotVecAr = numpy.array ( rotVec.data() )


    m1M = mdata.xyz_to_ijk_transform
    rmat = mdata.full_matrix()

    apoints = _multiscale.get_atom_coordinates ( atoms, transformed = False )
    point_weights = None
    #point_weights = numpy.ones ( len(apoints), numpy.float32 )

    gradients = VolumeData.interpolate_volume_gradient(apoints, m1M, rmat)[0]
    torque_axis = _distances.torque(apoints, point_weights, gradients, center )
    #print "torque axis 0:", torque_axis

    dt = numpy.dot ( torque_axis, rotVec.data() )
    torque_axis = rotVecAr if dt > 0 else rotVecAr * -1.0

    na = Matrix.norm(torque_axis)
    if na == 0 :
        #torque_axis = (0,0,1)
        pass
    else :
        torque_axis /= na
        angle = angle_step(torque_axis, apoints, center, m1M, step_size)
        #print "torque axis:", torque_axis, "angle:", angle
        move_tf = Matrix.rotation_transform ( torque_axis, angle, center )
        xf =  Matrix.chimera_xform ( move_tf )

        for at in atoms :
            at.setCoord ( xf.apply (at.coord()) )





def AtTorques (bond, atoms, mdata ) :

    apoints = _multiscale.get_atom_coordinates ( atoms, transformed = False )

    at1, at2 = bond.atoms
    center = at1.coord().data()
    rotVec = at2.coord() - at1.coord();
    rotVec.normalize()

    point_weights = None
    #point_weights = numpy.ones ( len(apoints), numpy.float32 )
    ijk_step_size_max = 0.5
    ijk_step_size_min = 0.01
    ijk_step_size = ijk_step_size_max

    #xyz_to_ijk_tf = inMap.data.xyz_to_ijk_transform
    #darray = inMap.data.full_matrix()
    #map_values, outside = VolumeData.interpolate_volume_data(fpoints, xyz_to_ijk_tf, darray)
    #olap0, cc0, other = FitMap.overlap_and_correlation ( fpoint_weights, map_values )

    m1M = mdata.xyz_to_ijk_transform
    rmat = mdata.full_matrix()

    values, outside = VolumeData.interpolate_volume_data ( apoints, m1M, rmat )
    avg0 = numpy.average ( values )
    last_avg = avg0
    print " - 0: %.3f" % last_avg

    xf = chimera.Xform.identity()

    rotVecAr = numpy.array ( rotVec.data() )
    totAngle = 0

    for i in range ( 100 ) :

        gradients = VolumeData.interpolate_volume_gradient(apoints, m1M, rmat)[0]
        torque_axis = _distances.torque(apoints, point_weights, gradients, center )
        #print "torque axis 0:", torque_axis

        dt = numpy.dot ( torque_axis, rotVec.data() )
        torque_axis = rotVecAr if dt > 0 else rotVecAr * -1.0

        na = Matrix.norm(torque_axis)
        if na == 0 :
            torque_axis = (0,0,1)
            angle = 0
            #print " - torque axis 0, stopping"
            break
        else :
            torque_axis /= na
            angle = angle_step(torque_axis, apoints, center, m1M, ijk_step_size)
            #print "torque axis:", torque_axis, "angle:", angle

        move_tf = Matrix.rotation_transform ( torque_axis, angle, center )
        _contour.affine_transform_vertices ( apoints, move_tf )

        values, outside = VolumeData.interpolate_volume_data ( apoints, m1M, rmat )
        avg1 = numpy.average ( values )

        print " %d - score: %.3f - step: %.3f, angle %.1f" % (i+1, avg1, ijk_step_size, angle)
        totAngle += angle

        xf.premultiply ( Matrix.chimera_xform ( move_tf ) )
        #m1M = Matrix.multiply_matrices( move_tf, m1M )

        if avg1 < last_avg :
            ijk_step_size = ijk_step_size / 2.0
            if ijk_step_size < ijk_step_size_min :
                print " - reached min step size, stopping"
                break

        last_avg = avg1

    #print " - dih %s(%d) - %s(%d), %d its, score: %.3f -> %.3f, step: %.3f, angle %.1f" % (at1.name, at1.residue.id.position, at2.name, at2.residue.id.position, i+1, avg0, last_avg, ijk_step_size, totAngle)

    #print " - moving atoms in res %d" % atoms[0].residue.id.position
    if 1 :
        for at in atoms :
            at.setCoord ( xf.apply (at.coord()) )

    return last_avg


# -----------------------------------------------------------------------------
# Return angle such that rotating point about given axis and center causes the
# largest motion in ijk space to equal ijk_step_size.
#
def angle_step(axis, points, center, xyz_to_ijk_transform, ijk_step_size):

    import Matrix as m
    tf = m.multiply_matrices(m.zero_translation(xyz_to_ijk_transform),
                             m.cross_product_transform(axis),
                             m.translation_matrix([-x for x in center]))

    import _distances as dist
    av = dist.maximum_norm(points, tf)

    if av > 0:
        from math import pi
        angle = (ijk_step_size / av) * 180.0/pi
    else:
        angle = 0
    return angle


def TorFitRandStep ( ress, inMap, stepSize, parent, data ) :

    tors, ressAtoms, cc0, atStep = data

    from random import random

    for tor in tors :

        #bond, ats1, ats2 = tor
        bond, ats1, ats2 = tor
        p2, p1 = bond.atoms[1].coord(), bond.atoms[0].coord()
        v = p2 - p1
        ang = (random()-0.5)*stepSize
        xf1 = chimera.Xform.translation ( p1.toVector() * -1 )
        xf1.premultiply ( chimera.Xform.rotation ( v, ang ) )
        xf1.premultiply ( chimera.Xform.translation ( p1.toVector() ) )
        for at in ats1 :
            at.setCoord ( xf1.apply ( at.coord() ) )

    cc, xf = FitAtoms ( ressAtoms, inMap )
    if cc < cc0 :
        #print "x"
        for at in ressAtoms :
            at.setCoord( at.coord0 )
    else :
        cc0 = cc
        print "%d_%.4f" % (atStep,cc),
        for at in ressAtoms :
            at.coord0 = xf.apply ( at.coord() )
            at.setCoord ( at.coord0 )

    if atStep >= 100 :
        print "done"
    else :
        print ".",
        atStep += 1
        data = [tors, ressAtoms, cc0, atStep]

        #parent.toplevel_widget.update_idletasks ()
        #from chimera import dialogs
        #dlg = dialogs.find ( "segment map", create=False )
        #dlg.toplevel_widget.update_idletasks ()
        #sz = chimera.viewer.windowSize
        #chimera.runCommand ( "windowsize %d %d" % (sz[0], sz[1]) )
        #dlg = dialogs.find ( "View Editor", create=False )
        #dlg.update()

        from chimera.tkgui import app
        # Make sure new dialogs are shown before returning focus to main window.
        app.update_idletasks()
        app.graphics.focus()

        parent.after(10, TorFitRandStep(ress, inMap, stepSize, parent, data))



def TorFitRand ( ress, inMap, stepSize, parent=None, task=None ) :

    tors = FindTors ( ress )
    print "%d tors" % len(tors)
    ressAtoms = []
    for r in ress :
        ressAtoms = ressAtoms + r.atoms

    from random import random

    cc0, xf = FitAtoms ( ressAtoms, inMap )
    for at in ressAtoms :
        at.setCoord ( xf.apply ( at.coord() ) )
        at.coord0 = at.coord()
    print "%.4f" % cc0,

    if 0 :
        atStep = 1
        data = [tors, ressAtoms, cc0, atStep]
        TorFitRandStep(ress, inMap, stepSize, parent, data)
        return


    for i in range ( 100 ) :

        for tor in tors :

            #bond, ats1, ats2 = tor
            bond, ats1, ats2 = tor
            p2, p1 = bond.atoms[1].coord(), bond.atoms[0].coord()
            v = p2 - p1
            ang = (random()-0.5)*stepSize
            xf1 = chimera.Xform.translation ( p1.toVector() * -1 )
            xf1.premultiply ( chimera.Xform.rotation ( v, ang ) )
            xf1.premultiply ( chimera.Xform.translation ( p1.toVector() ) )
            for at in ats1 :
                at.setCoord ( xf1.apply ( at.coord() ) )

        cc, xf = FitAtoms ( ressAtoms, inMap )
        if cc < cc0 :
            #print "x"
            for at in ressAtoms :
                at.setCoord( at.coord0 )
        else :
            cc0 = cc
            print "%d_%.4f" % (i,cc),
            for at in ressAtoms :
                at.coord0 = xf.apply ( at.coord() )
                at.setCoord ( at.coord0 )

        if task != None :
            task.updateStatus ( "Torsion fit random %d" % i )
            print "."

        #from chimera.tkgui import app
        #app.update_idletasks()
        #app.graphics.focus()






def TorFitRSel ( selBonds, inMap, stepSize, doRigidFit=False ) :

    ress = []

    res = selBonds[0].atoms[0].residue
    mol = res.molecule
    SetBBAts ( mol )

    rmap = {}
    #for r in mol.residues :
    #    rmap[r.id.chainId + "%d"%r.id.position] = r
    for b in selBonds :
        rmap[b.atoms[0].residue] = 1
        rmap[b.atoms[1].residue] = 1

    if res.type in protein3to1 or res.type in nucleic1to3 :
        for r in res.molecule.residues :
            if r.id.chainId == res.id.chainId :
                ress.append ( r )
    else :
        ress = [res]

    tors = FindTors ( ress, selBonds )
    print " - %d total tors" % len(tors)

    selTors = []
    for tor in tors :
        bond, ats1, ats2 = tor
        if bond in selBonds :
            selTors.append ( tor )

    if len(selTors) == 0 :
        print "sel tors not found"
        return None


    scoreAtoms = []
    for r in rmap.keys() :
        #fitAtoms.extend ( r.atoms )
        scoreAtoms.extend ( r.bbAtoms )

    allAtoms = []
    for r in rmap.keys() :
        allAtoms.extend ( r.atoms )

    print " - %d atoms in %d res, score %d atoms" % (len(allAtoms), len(rmap.keys()), len(scoreAtoms))

    from random import random

    cc0 = None
    if doRigidFit:
        cc0, xf = FitAtoms ( allAtoms, inMap )
        for at in allAtoms :
            at.setCoord ( xf.apply ( at.coord() ) )
            at.coord0 = at.coord()
    else :
        cc0 = FitScore ( scoreAtoms, inMap )
        for at in allAtoms :
            at.coord0 = at.coord()


    print "Fitting %d tors - in map: %s, cc: %.4f" % (len(selTors), inMap.name, cc0)

    for i in range ( 100 ) :

        for tor in selTors :

            bond, ats1, ats2 = tor
            p2, p1 = bond.atoms[1].coord(), bond.atoms[0].coord()
            v = p2 - p1

            ang = (random()-0.5)*stepSize

            xf1 = chimera.Xform.translation ( p1.toVector() * -1 )
            xf1.premultiply ( chimera.Xform.rotation ( v, ang ) )
            xf1.premultiply ( chimera.Xform.translation ( p1.toVector() ) )

            for at in ats1 :
                at.setCoord ( xf1.apply ( at.coord() ) )


        if doRigidFit :
            cc, xf = FitAtoms ( allAtoms, inMap )
            if cc < cc0 :
                print ".",
                for at in allAtoms :
                    at.setCoord( at.coord0 )
            else :
                #print "%d|%.4f" % (i,cc),
                print "%.4f" % (cc),
                cc0 = cc
                for at in allAtoms :
                    at.coord0 = xf.apply ( at.coord() )
                    at.setCoord ( at.coord0 )


        else :
            cc = FitScore ( scoreAtoms, inMap )
            if cc < cc0 :
                print ".",
                for at in allAtoms :
                    at.setCoord( at.coord0 )
            else :
                #print "%d|%.4f" % (i,cc),
                print "%.4f" % (cc),
                cc0 = cc
                for at in allAtoms :
                    at.coord0 = at.coord()



    print ""




def FindTors ( ress, selBonds=None ) :

    amap = {}
    for res in ress :
        for at in res.atoms :
            amap[at] = 1

    mol = ress[0].molecule
    doBonds = []
    for b in mol.bonds :
        at1, at2 = b.atoms
        if at1 in amap or at2 in amap :
            doBonds.append ( b )

    bonds = selBonds
    bonds = doBonds
    if bonds == None :
        # use all bonds - can be slow for large proteins/rna
        print " - using all bonds in res"
        bonds = []
        for b in mol.bonds :
            at1, at2 = b.atoms
            if at1 in amap or at2 in amap :
                bonds.append (b)

    print " - %d sel ress, %d/%d atoms, %d/%d bonds" % ( len(ress), len(amap), len(mol.atoms), len(bonds), len(mol.bonds) )

    tors = []
    for b in bonds :

        at1, at2 = b.atoms

        cycle, ats1 = BondGo ( at1, at2 )
        if cycle :
            #print "cycle"
            continue

        ats2 = {}
        cycle, ats2 = BondGo ( at2, at1 )
        if cycle :
            #print "cycle"
            continue

        if len(ats1) == 0 or len(ats2) == 0 :
            continue

        a1 = ats1 if len(ats1) < len(ats2) else ats2
        a2 = ats1 if a1 == ats2 else ats2
        tors.append ( [b, a1, a2] )
        #tors.append ( [b, a1] )

        if 0:
            print "bond %s-%s " % (at1.name, at2.name)
            print " -1 ",
            for at in a1 :
                print "%s.%d" % (at.name, at.residue.id.position),
            print ""
            if 0 :
                print " -2 ",
                for at in a2 :
                    print "%s.%d" % (at.name, at.residue.id.position),
                print ""

        #break

    return tors


def BondAts ( at1, at2 ) :

    cycle, ats1 = BondGo ( at1, at2 )
    if cycle :
        #print "cycle"
        return None, None

    ats2 = {}
    cycle, ats2 = BondGo ( at2, at1 )
    if cycle :
        #print "cycle"
        return None, None

    return ats1, ats2


def RotBond ( at1, at2, ats, ang ) :

    p2, p1 = at2.coord(), at1.coord()
    v = p2 - p1
    #ang = (random()-0.5)*stepSize
    xf1 = chimera.Xform.translation ( p1.toVector() * -1 )
    xf1.premultiply ( chimera.Xform.rotation ( v, ang ) )
    xf1.premultiply ( chimera.Xform.translation ( p1.toVector() ) )
    for at in ats :
        at.setCoord ( xf1.apply ( at.coord() ) )



def FindTorsDir ( ress, selBonds=None ) :

    # same as above but preserves direction

    bondedAtoms = {}
    amap = {}
    for res in ress :
        for at in res.atoms :
            amap[at] = 1
            bondedAtoms[at] = []

    for b in ress[0].molecule.bonds :
        at1, at2 = b.atoms
        if at1 in amap or at2 in amap :
            bondedAtoms[at1].append ( at2 )
            bondedAtoms[at2].append ( at1 )

    bonds = selBonds
    if selBonds == None :
        # use all bonds - can be slow for large proteins/rna
        bonds = []
        for b in ress[0].molecule.bonds :
            if at1 in amap or at2 in amap :
                bonds.append (b)

    print "%d atoms, %d bonds" % ( len(res.atoms), len(bonds) )

    tors = []
    for b in bonds :

        at1, at2 = b.atoms

        cycle, ats1 = BondGo ( at1, at2, bondedAtoms )
        if cycle :
            #print "cycle"
            continue

        ats2 = {}
        cycle, ats2 = BondGo ( at2, at1, bondedAtoms )
        if cycle :
            #print "cycle"
            continue

        if len(ats1) == 0 or len(ats2) == 0 :
            continue

        minDir = 1.0 if len(ats1) < len(ats2) else -1.0
        tors.append ( [b, ats1, ats2, minDir] )

        if 0:
            print "bond %s-%s " % (at1.name, at2.name),
            for at in ats :
                print at.name,
            print ""

        #break

    return tors



def BondGo ( at0, at1 ) :

    visAts = { at0:1, at1:1 }

    first = True
    cycle = False
    Q = [at1]

    while len(Q) > 0 :

        at = Q.pop(0)
        #print "%s " % at.name
        visAts[at] = 1

        for at2 in at.neighbors :
            #print " -> %s " % (at2.name),

            if not first and at2 == at0 :
                #print "cycle"
                cycle = True
                break

            if not at2 in visAts :
                Q.append ( at2 )
                #print " > "

        first = False
        if cycle :
            break

    del visAts[at0]
    del visAts[at1]
    return cycle, visAts.keys()


# uses atoms in residue only
def ConAts ( ress, maxDepth=3 ) :

    #bondedAtoms = {}
    amap = {}
    for res in ress :
        for at in res.atoms :
            amap[at] = 1
            #bondedAtoms[at] = []

    bonds = []
    for b in ress[0].molecule.bonds :
        at1, at2 = b.atoms
        if at1 in amap or at2 in amap :
            bonds.append (b)
            #bondedAtoms[at1].append ( at2 )
            #bondedAtoms[at2].append ( at1 )

    print " - conAts - %d atoms, %d bonds" % ( len(amap), len(bonds) )

    conAts = {}
    for res in ress :
        for at in res.atoms :
            conAts[at] = ConAtsGo ( at, maxDepth=3 )
            #print at.name, " : ",
            #for atc in conAts[at].keys() :
            #    print atc.name,
            #print ""


    return conAts


# uses list of atoms
def ConAts2 ( atoms ) :

    bondedAtoms = {}
    amap = {}
    for at in atoms :
        amap[at] = 1
        bondedAtoms[at] = []

    bonds = []
    for b in atoms[0].molecule.bonds :
        at1, at2 = b.atoms
        if at1 in amap or at2 in amap :
            bonds.append (b)
            if at1 in bondedAtoms :
                bondedAtoms[at1].append ( at2 )
            if at2 in bondedAtoms :
                bondedAtoms[at2].append ( at1 )

    print "ConAts - %d atoms, %d bonds" % ( len(atoms), len(bonds) )

    conAts = {}
    for at in atoms :

        conAts[at] = ConAtsGo ( at, bondedAtoms )

        #print at.name, " : ",
        #for atc in conAts[at].keys() :
        #    print atc.name,
        #print ""


    return conAts


def ConAtsGo ( startAt, maxDepth=3 ) :

    visAts = { startAt:1 }
    depthAt = { startAt:1 }
    Q = [startAt]

    while len(Q) > 0 :
        at = Q.pop(0)
        #print "%s " % at.name
        visAts[at] = 1
        if depthAt[at] > maxDepth :
            continue

        for at2 in at.neighbors :
            #print " -> %s " % (at2.name),
            if not at2 in visAts :
                Q.append ( at2 )
                depthAt[at2] = depthAt[at]+1
                #print " > "

    return visAts




def FitAtomsXf ( atoms, xf0, inMap, doTranslate = True, doRotate = True ) :

    #fpoints, fpoint_weights, darray = RessPtsInMap (ress, regs, segMap)

    #print inMap.name

    ats = [at for at in atoms if not at.element.name == "H"]
    fpoints = _multiscale.get_atom_coordinates ( ats, transformed = False )

    _contour.affine_transform_vertices ( fpoints,  Matrix.xform_matrix(xf0) )
    #_contour.affine_transform_vertices ( fpoints, inMap.openState.xform  )
    fpoint_weights = numpy.ones ( len(fpoints), numpy.float32 )

    #fpoints = numpy.array ( fpoints, dtype=numpy.float32 )

    xyz_to_ijk_tf = inMap.data.xyz_to_ijk_transform
    darray = inMap.data.full_matrix()
    map_values, outside = VolumeData.interpolate_volume_data(fpoints, xyz_to_ijk_tf, darray)
    olap0, cc0, other = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    #print cc0,

    move_tf, stats = FitMap.locate_maximum(fpoints, fpoint_weights,
                                    darray, xyz_to_ijk_tf,
                                    max_steps = 1000,
                                    ijk_step_size_min = 0.01,
                                    ijk_step_size_max = 0.5,
                                    optimize_translation = doTranslate,
                                    optimize_rotation = doRotate,
                                    metric = 'sum product', # 'correlation' or 'correlation about mean'
                                    request_stop_cb = None)

    xf = Matrix.chimera_xform ( move_tf )
    cc1 = float ( stats['correlation'] )
    #ApplyXf ( ress, xf )
    #print " -> ", cc1

    return cc0, cc1, xf




def FitAtoms ( atoms, inMap, doTranslate = True, doRotate = True ) :

    #fpoints, fpoint_weights, darray = RessPtsInMap (ress, regs, segMap)

    #ats = [at for at in atoms if not at.element.name == "H"]
    fpoints = _multiscale.get_atom_coordinates ( atoms, transformed = False )
    #_contour.affine_transform_vertices ( fpoints,  Matrix.xform_matrix(xf0) )
    fpoint_weights = numpy.ones ( len(fpoints), numpy.float32 )

    #fpoints = numpy.array ( fpoints, dtype=numpy.float32 )

    xyz_to_ijk_tf = inMap.data.xyz_to_ijk_transform
    darray = inMap.data.full_matrix()

    #map_values, outside = VolumeData.interpolate_volume_data(fpoints, xyz_to_ijk_tf, darray)
    #olap0, cc0, other = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    #print cc0,

    move_tf, stats = FitMap.locate_maximum(fpoints, fpoint_weights,
                                    darray, xyz_to_ijk_tf,
                                    max_steps = 1000,
                                    ijk_step_size_min = 0.01,
                                    ijk_step_size_max = 0.5,
                                    optimize_translation = doTranslate,
                                    optimize_rotation = doRotate,
                                    metric = 'sum product',
                                    request_stop_cb = None)

    xf = Matrix.chimera_xform ( move_tf )
    cc1 = float ( stats['correlation'] )
    #ApplyXf ( ress, xf )
    #print " -> ", cc1

    return cc1, xf



def FitScore ( atoms, inMap ) :

    fpoints = _multiscale.get_atom_coordinates ( atoms, transformed = False )
    #_contour.affine_transform_vertices ( fpoints,  Matrix.xform_matrix(xf0) )
    fpoint_weights = numpy.ones ( len(fpoints), numpy.float32 )

    xyz_to_ijk_tf = inMap.data.xyz_to_ijk_transform
    darray = inMap.data.full_matrix()

    map_values, outside = VolumeData.interpolate_volume_data(fpoints, xyz_to_ijk_tf, darray)
    olap0, cc0, other = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    return cc0




def RegsPtsInMol ( regs, toMol ) :

    regsPoints = regs[0].points().astype ( numpy.float32 )
    for r in regs[1:] :
        npoints = r.points().astype ( numpy.float32 )
        regsPoints = numpy.concatenate ( [regsPoints, npoints], axis=0 )

    segMap = regs[0].segmentation.seg_map
    _contour.affine_transform_vertices ( regsPoints, segMap.data.ijk_to_xyz_transform )
    _contour.affine_transform_vertices ( regsPoints, Matrix.xform_matrix( segMap.openState.xform ) )
    _contour.affine_transform_vertices ( regsPoints, Matrix.xform_matrix( toMol.openState.xform.inverse() ) )

    return regsPoints



def AlignResToRegs ( res, regs ) :

    regsPoints = RegsPtsInMol ( regs, res.molecule )

    from axes import prAxes

    regsC, regsU, regsS, regsV = prAxes ( regsPoints )
    #print regsC, regsU

    ats = [at for at in res.atoms if not at.element.name == "H"]
    molPoints = _multiscale.get_atom_coordinates ( ats, transformed = False )
    molC, molU, molS, molV = prAxes ( molPoints )
    #print molC, molU

    xf = chimera.Xform.translation ( molC * -1.0 )
    #xf.premultiply ( chimera.Xform.rotation(rax, ang) )
    xf.premultiply ( chimera.Xform.translation ( regsC  ) )

    for at in res.atoms :
        at.setCoord ( xf.apply ( at.coord() ) )




def uniform_rota_xfs ( num ) :

    N = int ( numpy.floor ( numpy.sqrt ( num ) ) )
    M = int ( numpy.floor ( num / N ) )

    thetas, phis = [], []
    from math import acos, sin, cos, sqrt, pi
    for k in range ( 1, N+1 ) :
        h = -1.0 + ( 2.0*float(k-1)/float(N-1) )
        phis.append ( acos(h) )
        thetas.append ( 0 if k == 1 or k == N else
                        (thetas[k-2] + 3.6/sqrt(N*(1.0-h**2.0))) % (2*pi) )

    xfs = []
    for theta, phi in zip(thetas, phis):
        for m in range ( M ) :
            rot = 2*pi*float(m)/float(M)
            #ralist.append((theta,phi,rot))
            v = chimera.Vector (sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi))
            xfR = chimera.Xform.rotation ( v, rot*180/pi )
            xfs.append ( xfR )

    return xfs






def AddNAG ( selAt, inMap, selReg ) :

    nmol = chimera.PDBio().readPDBfile ( "/Users/greg/Dropbox/_mol/Segger/_param/nag.pdb" )[0]
    print " - read %s - %d atoms" % ( nmol.name, len(nmol.atoms) )


    if selAt.residue.type == "ASN" :

        pCG = selAt.residue.atomsMap["CG"][0].coord()
        pN = selAt.residue.atomsMap["ND2"][0].coord()
        pO = selAt.residue.atomsMap["OD1"][0].coord()
        pO1_ = nmol.residues[0].atomsMap["O1"][0].coord()
        pC1_ = nmol.residues[0].atomsMap["C1"][0].coord()
        xf = ConnectXf ( pCG, pN, pO, 1.450, 124.669, pO1_, pC1_ )

        addRes = nmol.residues[0]
        toMol = selAt.molecule
        toChain = selAt.residue.id.chainId
        nres = AddResToMol ( addRes, toMol, toChain, xf, withoutAtoms=["O1"] )
        atN = selAt.residue.atomsMap["ND2"][0]
        atC1 = nres.atomsMap["C1"][0]
        nb = selAt.molecule.newBond ( atN, atC1 )
        nb.display = nb.Smart
        nb.drawMode = nb.Stick

        OptDihedral ( pN, atC1.coord(), nres.atoms, inMap, selReg  )

    elif selAt.residue.type == "NAG" :

        pC4 = selAt.residue.atomsMap["C4"][0].coord()
        pO4 = selAt.residue.atomsMap["O4"][0].coord()
        pC3 = selAt.residue.atomsMap["C3"][0].coord()
        pO1_ = nmol.residues[0].atomsMap["O1"][0].coord()
        pC1_ = nmol.residues[0].atomsMap["C1"][0].coord()
        xf = ConnectXf ( pC4, pO4, pC3, 1.433, 118.567, pO1_, pC1_ )

        addRes = nmol.residues[0]
        toMol = selAt.molecule
        toChain = selAt.residue.id.chainId
        nres = AddResToMol ( addRes, toMol, toChain, xf, withoutAtoms=["O1"] )
        atO4 = selAt.residue.atomsMap["O4"][0]
        atC1 = nres.atomsMap["C1"][0]
        nb = selAt.molecule.newBond ( atO4, atC1 )
        nb.display = nb.Smart
        nb.drawMode = nb.Stick

        SetDihedral ( pC3, pC4, pO4, atC1.coord(), 64.190, nres.atoms  )

        OptDihedral ( pO4, atC1.coord(), nres.atoms, inMap, selReg  )



def AddBMA ( selAt, inMap, selReg ) :

    nmol = chimera.PDBio().readPDBfile ( "/Users/greg/Dropbox/_mol/Segger/_param/bma.pdb" )[0]
    print " - read %s - %d atoms" % ( nmol.name, len(nmol.atoms) )


    if selAt.residue.type == "NAG" :

        pC4 = selAt.residue.atomsMap["C4"][0].coord()
        pO4 = selAt.residue.atomsMap["O4"][0].coord()
        pC3 = selAt.residue.atomsMap["C3"][0].coord()
        pO1_ = nmol.residues[0].atomsMap["O1"][0].coord()
        pC1_ = nmol.residues[0].atomsMap["C1"][0].coord()
        xf = ConnectXf ( pC4, pO4, pC3, 1.433, 109.147, pO1_, pC1_ )

        addRes = nmol.residues[0]
        toMol = selAt.molecule
        toChain = selAt.residue.id.chainId
        nres = AddResToMol ( addRes, toMol, toChain, xf, withoutAtoms=["O1"] )
        atO4 = selAt.residue.atomsMap["O4"][0]
        atC1 = nres.atomsMap["C1"][0]
        nb = selAt.molecule.newBond ( atO4, atC1 )
        nb.display = nb.Smart
        nb.drawMode = nb.Stick

        SetDihedral ( pC3, pC4, pO4, atC1.coord(), 137.239, nres.atoms  )

        OptDihedral ( pO4, atC1.coord(), nres.atoms, inMap, selReg  )



def AddMAN ( selAt, inMap, selReg ) :

    nmol = chimera.PDBio().readPDBfile ( "/Users/greg/Dropbox/_mol/Segger/_param/man.pdb" )[0]
    print " - read %s - %d atoms" % ( nmol.name, len(nmol.atoms) )


    if selAt.residue.type == "BMA" :

        if selAt.name == "O6" :

            pC6 = selAt.residue.atomsMap["C6"][0].coord()
            pO6 = selAt.residue.atomsMap["O6"][0].coord()
            pC5 = selAt.residue.atomsMap["C5"][0].coord()
            pO1_ = nmol.residues[0].atomsMap["O1"][0].coord()
            pC1_ = nmol.residues[0].atomsMap["C1"][0].coord()
            xf = ConnectXf ( pC6, pO6, pC5, 1.425, 115.695, pO1_, pC1_ )

            addRes = nmol.residues[0]
            toMol = selAt.molecule
            toChain = selAt.residue.id.chainId
            nres = AddResToMol ( addRes, toMol, toChain, xf, withoutAtoms=["O1"] )
            atO6 = selAt.residue.atomsMap["O6"][0]
            atC1 = nres.atomsMap["C1"][0]
            nb = selAt.molecule.newBond ( atO6, atC1 )
            nb.display = nb.Smart
            nb.drawMode = nb.Stick

            SetDihedral ( pC5, pC6, pO6, atC1.coord(), 177.537, nres.atoms  )

            OptDihedral ( pO6, atC1.coord(), nres.atoms, inMap, selReg  )


        if selAt.name == "O3" :

            pC3 = selAt.residue.atomsMap["C3"][0].coord()
            pO3 = selAt.residue.atomsMap["O3"][0].coord()
            pC4 = selAt.residue.atomsMap["C4"][0].coord()
            pO1_ = nmol.residues[0].atomsMap["O1"][0].coord()
            pC1_ = nmol.residues[0].atomsMap["C1"][0].coord()
            xf = ConnectXf ( pC3, pO3, pC4, 1.475, 110.731, pO1_, pC1_ )

            addRes = nmol.residues[0]
            toMol = selAt.molecule
            toChain = selAt.residue.id.chainId
            nres = AddResToMol ( addRes, toMol, toChain, xf, withoutAtoms=["O1"] )
            atO3 = selAt.residue.atomsMap["O3"][0]
            atC1 = nres.atomsMap["C1"][0]
            nb = selAt.molecule.newBond ( atO3, atC1 )
            nb.display = nb.Smart
            nb.drawMode = nb.Stick

            SetDihedral ( pC4, pC3, pO3, atC1.coord(), 144.288, nres.atoms  )

            OptDihedral ( pO3, atC1.coord(), nres.atoms, inMap, selReg  )






def AddRes ( resType, selAt, inMap, selReg ) :

    print "Adding:", resType
    for i in range ( len(resType) ) :

        RT = resType[i].upper()
        if RT in protein1to3 :
            rtype = protein1to3[RT]
            print " %s -prot> %s" % (RT, rtype)
            AddProtRes ( rtype, selAt, inMap, selReg )



def AddNuc ( resType, selAt, inMap, selReg ) :

    print "Adding:", resType
    for i in range ( len(resType) ) :

        RT = resType[i].upper()

        if RT in nucleic1to3 :
            rtype = nucleic1to3[RT]
            print " %s -nucleic> %s" % (RT, rtype)
            AddNucRes ( rtype, selAt, inMap, selReg )


def AddNucRes ( rtype, selAt, inMap, selReg ) :

    from mmcif import ParamPathPdb
    rpath = ParamPathPdb ( rtype )
    print " - res from %s" % rpath

    #nmol = chimera.PDBio().readPDBfile ( "/Users/greg/Dropbox/_mol/Segger/_param/%s.pdb" % rtype )[0]
    nmol = chimera.PDBio().readPDBfile ( rpath )[0]
    print " - read %s - %d atoms" % ( nmol.name, len(nmol.atoms) )

    #chimera.openModels.add ( [nmol] )

    addRes = nmol.residues[0]
    toMol = selAt.molecule
    toChain = selAt.residue.id.chainId
    rmap = {}
    for r in toMol.residues :
        if r.id.chainId == toChain :
            rmap[r.id.position] = r

    if selAt.name == "P" :

        rid = selAt.residue.id.position - 1
        if rid in rmap :
            umsg ( "Residue at position %d already exists" % rid )
            return

        pP = selAt.coord() #.residue.atomsMap["C"][0].coord()
        pO5p = selAt.residue.atomsMap["O5'"][0].coord()
        pC5p = selAt.residue.atomsMap["C5'"][0].coord()
        pO3p_ = addRes.atomsMap["O3'"][0].coord()
        pC3p_ = addRes.atomsMap["C3'"][0].coord()

        xf = ConnectXfR ( pO5p, pP, pC5p, pO3p_, pC3p_ ) #  pCA, pC, pO, pN_, pCA_

        nres = AddResToMol ( addRes, toMol, toChain, xf, withoutAtoms=[], rid=rid )
        atP = selAt.residue.atomsMap["P"][0]
        atO3p = nres.atomsMap["O3'"][0]
        nb = selAt.molecule.newBond ( atP, atO3p )
        nb.display = nb.Smart
        nb.drawMode = nb.Stick

        #OptDihedral ( pN, atC1.coord(), nres.atoms, inMap, selReg  )

        if 0 and rid + 1 in rmap :
            toRes = rmap[rid + 1]
            print "- connect to %s.%d" % (toRes.type, toRes.id.position)

            atC = nres.atomsMap["C"][0]
            atN = toRes.atomsMap["N"][0]
            nb = selAt.molecule.newBond ( atC, atN )
            nb.display = nb.Smart
            nb.drawMode = nb.Stick

    elif selAt.name == "O3'" :
        print " - todo - add to O3'"

    else :
        print " - replace %s -> %s" % ( selAt.residue.type, rtype )




def AddProtRes ( rtype, selAt, inMap, selReg ) :

    nmol = chimera.PDBio().readPDBfile ( "/Users/greg/Dropbox/_mol/Segger/_param/%s.pdb" % rtype.lower() )[0]
    print " - read %s - %d atoms" % ( nmol.name, len(nmol.atoms) )

    addRes = nmol.residues[0]
    toMol = selAt.molecule
    toChain = selAt.residue.id.chainId
    rmap = {}
    for r in toMol.residues :
        if r.id.chainId == toChain :
            rmap[r.id.position] = r

    if selAt.name == "C" :

        rid = selAt.residue.id.position + 1
        if rid in rmap :
            umsg ( "Residue at position %d already exists" % rid )
            return

        pC = selAt.coord() #.residue.atomsMap["C"][0].coord()
        pCA = selAt.residue.atomsMap["CA"][0].coord()
        pO = selAt.residue.atomsMap["O"][0].coord()
        pN_ = nmol.residues[0].atomsMap["N"][0].coord()
        pCA_ = nmol.residues[0].atomsMap["CA"][0].coord()
        xf = ConnectXfP ( pCA, pC, pO, pN_, pCA_ ) #  pCA, pC, pO, pN_, pCA_

        nres = AddResToMol ( addRes, toMol, toChain, xf, withoutAtoms=['OXT'], rid=rid )
        atC = selAt.residue.atomsMap["C"][0]
        atN = nres.atomsMap["N"][0]
        nb = selAt.molecule.newBond ( atC, atN )
        nb.display = nb.Smart
        nb.drawMode = nb.Stick

        #OptDihedral ( pN, atC1.coord(), nres.atoms, inMap, selReg  )


        if 0 and rid + 1 in rmap :
            toRes = rmap[rid + 1]
            print "- connect to %s.%d" % (toRes.type, toRes.id.position)

            atC = nres.atomsMap["C"][0]
            atN = toRes.atomsMap["N"][0]
            nb = selAt.molecule.newBond ( atC, atN )
            nb.display = nb.Smart
            nb.drawMode = nb.Stick






def SetDihedral ( p1, p2, p3, p4, toAngDeg, atoms ) :

    b1, b2, b3 = p2 - p1, p3 - p2, p4 - p3

    n1 = chimera.cross ( b1, b2 ); n1.normalize()
    n2 = chimera.cross ( b2, b3 ); n2.normalize()
    m1 = chimera.cross ( n1, b2 ); m1.normalize()

    x, y = n1 * n2, m1 * n2

    A = -1.0 * numpy.arctan2 (y, x) * 180.0 / numpy.pi

    print " - dih: %.3f -> %.3f" % (A, toAngDeg)

    V = p3.toVector()
    b2.normalize()
    xf = chimera.Xform.translation ( V * -1.0 )
    xf.premultiply ( chimera.Xform.rotation(b2, toAngDeg-A) )
    xf.premultiply ( chimera.Xform.translation ( V ) )

    for at in atoms :
        at.setCoord ( xf.apply(at.coord()) )





def ConnectXf ( P0, P1, P2, bondLength, angleDeg, M0, M1 ) :

    v1 = P1 - P0; v1.normalize()
    v2 = P2 - P0; v2.normalize()
    vA = chimera.cross ( v2, v1 ); vA.normalize()

    vP = chimera.Xform.rotation (vA, angleDeg) .apply ( v1*-1.0 )
    vP.normalize()
    pP = P1 + vP * bondLength

    vM = M1 - M0; vM.normalize()

    # align vM to vP, put M1 to pP

    rax = chimera.cross ( vM, vP ); rax.normalize()
    ang = numpy.arccos ( vM * vP ) * 180.0 / numpy.pi

    xf = chimera.Xform.translation ( M1.toVector() * -1.0 )
    xf.premultiply ( chimera.Xform.rotation(rax, ang) )
    xf.premultiply ( chimera.Xform.translation ( pP.toVector() ) )

    return xf



def ConnectXfP ( pCA, pC, pO, pN_, pCA_ ) :

    v1 = pCA - pC; v1.normalize()
    v2 = pO - pC; v2.normalize()
    vA = chimera.cross ( v2, v1 ); vA.normalize()

    vN = chimera.Xform.rotation (vA, 114.017) .apply ( v1 )
    vN.normalize()
    pN = pC + vN * 1.302

    vCA = chimera.Xform.rotation (vA*-1.0, 116.766) .apply ( vN * -1.0 )
    vCA.normalize()
    #pCA = pN + vCA * 1.373

    vCA_ = pCA_ - pN_; vCA_.normalize()

    # align vCA_ to vCA, put pN_ to pN

    rax = chimera.cross ( vCA_, vCA );
    if rax.length < 1e-4 :
        xf = chimera.Xform.translation ( pN_.toVector() * -1.0 )
        xf.premultiply ( chimera.Xform.translation ( pN.toVector() ) )
        return xf

    else :
        rax.normalize()
        ang = numpy.arccos ( vCA_ * vCA ) * 180.0 / numpy.pi
        xf = chimera.Xform.translation ( pN_.toVector() * -1.0 )
        xf.premultiply ( chimera.Xform.rotation(rax, ang) )
        xf.premultiply ( chimera.Xform.translation ( pN.toVector() ) )
        return xf



def ConnectXfR ( pO5p, pP, pC5p, pO3p_, pC3p_ ) :

    v1 = pP - pO5p; v1.normalize()
    v2 = pC5p - pO5p; v2.normalize()
    vA = chimera.cross ( v2, v1 ); vA.normalize()

    vN = chimera.Xform.rotation (vA, 104.0) .apply ( v1*-1.0 )
    vN.normalize()
    pO3p = pP + vN * 1.608

    vA = chimera.cross ( vN, v1 ); vA.normalize()
    #print vA
    v = chimera.Xform.rotation (vA*-1.0, 119.7) .apply ( vN*-1.0 )
    v.normalize()
    #print v

    v_ = pC3p_ - pO3p_; v_.normalize()

    # align v_ to v, put pO3p_ to pO3p
    rax = chimera.cross ( v_, v );
    if rax.length < 1e-3 :
        xf = chimera.Xform.translation ( pO3p_.toVector() * -1.0 )
        xf.premultiply ( chimera.Xform.translation ( pO3p.toVector() ) )
        return xf

    else :
        ang = numpy.arccos ( v * v_ ) * 180.0 / numpy.pi
        #print "ang0: %.3f" % ang
        rax.normalize()
        print ang
        xf = chimera.Xform.translation ( pO3p_.toVector() * -1.0 )
        xf.premultiply ( chimera.Xform.rotation(rax, ang) )
        xf.premultiply ( chimera.Xform.translation ( pO3p.toVector() ) )

        pO3p_ = xf.apply ( pO3p_ ); print "O3p: ", pO3p_
        pC3p_ = xf.apply ( pC3p_ ); print "C3p: ", pC3p_

        v_ = pC3p_ - pO3p_; v_.normalize()
        v = pP - pO3p_; v.normalize()
        rax = chimera.cross ( v, v_ );
        ang = numpy.arccos ( v * v_ ) * 180.0 / numpy.pi
        print "- new ang:", ang

        return xf



def AddResToMol ( res, toMol, toChain, xf, withoutAtoms, rid=None, asType=None ) :

    if rid == None :
        rid = 0
        for r in toMol.residues :
            if r.id.chainId == toChain :
                if r.id.position > rid :
                    rid = r.id.position
        rid += 1

    if asType == None :
        asType = res.type

    aMap = {}
    nres = toMol.newResidue ( asType, chimera.MolResId(toChain, rid))

    if hasattr ( res, 'isHelix' ) :
        nres.isHelix = res.isHelix
    if hasattr ( res, 'isSheet' ) :
        nres.isSheet = res.isSheet

    for at in res.atoms :
        if at.element.name == "H" :
            continue
        elif 1 and at.name in withoutAtoms :
            continue

        nat = toMol.newAtom (at.name, chimera.Element(at.element.number))
        aMap[at] = nat
        nres.addAtom( nat )
        nat.drawMode = nat.EndCap
        nat.setCoord ( xf.apply ( at.coord()) )
        nat.display = True
        if nat.element.name.upper() in atomColors : nat.color = atomColors[nat.element.name.upper()]

    for bond in res.molecule.bonds :
        if bond.atoms[0] in aMap and bond.atoms[1] in aMap :
            nb = toMol.newBond ( aMap[bond.atoms[0]], aMap[bond.atoms[1]] )
            nb.display = nb.Smart
            nb.drawMode = nb.Stick

    return nres



def OptDihedral ( P, P2, forAtoms, inMap, selRegs ) :

    V = P2 - P

    dmap, rdata, rmat = None, None, None

    if selRegs != None :
        dmap = selRegs[0].segmentation.seg_map
        print " - seg map:", dmap.name
        zoneR = dmap.data.step[0]/2.0
        rpoints = numpy.concatenate ( [r.map_points() for r in selRegs], axis=0 ).astype ( numpy.float32 )
        rdata = VolumeData.zone_masked_grid_data ( dmap.data, rpoints, zoneR )
        rmat = rdata.matrix()

    elif inMap != None :
        dmap = inMap
        print " - in map:", dmap.name
        rdata = dmap.data
        rmat = dmap.full_matrix()

    ##gdata = VolumeData.Array_Grid_Data ( ndata.full_matrix(), segMap.data.origin, segMap.data.step, segMap.data.cell_angles, name = "atom masked" )
    #nv = VolumeViewer.volume.volume_from_grid_data ( rdata )
    #nv.name = "helix mask 2"

    maxAng, maxD, angD = 0, -1e9, 1.0
    for ang in range ( 0, int(round(360.0/angD)), 1 ) :

        xf = chimera.Xform.translation ( P.toVector() * -1.0 )
        xf.premultiply ( chimera.Xform.rotation(V, angD) )
        xf.premultiply ( chimera.Xform.translation ( P.toVector() ) )

        for at in forAtoms :
            at.setCoord ( xf.apply(at.coord()) )

        if dmap != None :
            points = _multiscale.get_atom_coordinates ( forAtoms, transformed = True )
            _contour.affine_transform_vertices ( points, Matrix.xform_matrix(dmap.openState.xform.inverse()) )
            values, outside = VolumeData.interpolate_volume_data ( points, rdata.xyz_to_ijk_transform, rmat )
            #values = nv.interpolated_values ( points, selAt.molecule.openState.xform )
            #olap, corr, other = overlap_and_correlation ( rpoint_weights, rmap_values )
            avgD = numpy.average ( values )
            #print "%.1f\t%.4f" % (ang, avgD)

            if avgD > maxD :
                maxD = avgD
                maxAng = round(float(ang)/angD)


    print "Max ang: %.3f" % maxAng
    xf = chimera.Xform.translation ( P.toVector() * -1.0 )
    xf.premultiply ( chimera.Xform.rotation(V, maxAng) )
    xf.premultiply ( chimera.Xform.translation ( P.toVector() ) )

    for at in forAtoms :
        at.setCoord ( xf.apply(at.coord()) )




def AddResN ( rtype, atRes ) :


    global mrPar


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









def SetBBAts ( mol ) :

    #print " - setting bbAts in %s" % mol.name
    for r in mol.residues :

        r.isProt = r.type in protein3to1
        r.isNA = r.type in nucleic3to1

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

                a.isBB = n=="P" or n=="O1P" or n=="O2P" or n=="OP1" or n=="OP2" or n=="O5'" or n=="C5'" or n=="O3'" or n=="C3'" or n=="C4'"
                a.isSugar = n=="C1'" or n=="C2'" or n=="O4'" or n=="O2'" or n=="C3'" or n=="C4'"

                #a.isBB = a.isBB or a.isSugar
                a.isBase = not a.isBB
                a.isSC = a.isBase

                if a.isBB or a.isSugar :
                    r.bbAtoms.append ( a )
                else :
                    r.scAtoms.append ( a )


        else :
            for a in r.atoms :
                a.isBB, a.isSC, a.isSugar, a.isBase = False, False, False, False

















# end
