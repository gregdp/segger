
import sys, os
import numpy


chiPath = "/Users/greg/Desktop/Chimera.app/"
chiBin = "/Users/greg/Desktop/Chimera.app/Contents/MacOS/chimera"
inPath = "/Users/greg/_data/toxop_/SPM2"
mapPath = "/Users/greg/_data/toxop_/SPM2/_microtubule_crop.mrc"
fitsPath = "/Users/greg/_data/toxop_/SPM2/by_ccm"
numProc = 6

if 1 :
    chiPath = "/storage/UCSF-Chimera64-1.15rc"
    chiBin = "/storage/UCSF-Chimera64-1.15rc/bin/chimera"
    inPath = "/storage/gregp/toxop"
    mapPath = "/storage/gregp/toxop/maps/microtubule_crop.mrc"
    fitsPath = "/storage/gregp/toxop/maps/microtubule_crop__ez_fits_res700_step10_8rot/by_ccm"
    numProc = 64


inPath = "/Volumes/S1/toxop_/"


gSigma = 2.0

def CalcQ ( fitf ) :

    print ("")
    print " - calcq for [%s]" % fitf

    #scriptPath = os.path.dirname(os.path.realpath(__file__))
    scriptPath = inPath + "/_mapqScript.py"
    print (" - creating script -> %s" % scriptPath)

    molPath = fitsPath + "/" + fitf

    fp = open ( scriptPath, "w" )
    fp.write ( "import chimera\n" )
    fp.write ( "import VolumeViewer\n" )
    fp.write ( "from Segger import qscores\n" )
    fp.write ( "from Segger import mmcif\n" )
    #fp.write ( "dmap = VolumeViewer.open_volume_file ( '%s', 'ccp4')[0]\n" % mapPath )
    #fp.write ( "mol = chimera.openModels.open ( '%s', type='PDB' )[0]\n" % molPath )
    fp.write ( "mol = chimera.openModels.list(modelTypes = [chimera.Molecule])[0]\n" )
    fp.write ( "dmap = chimera.openModels.list(modelTypes = [VolumeViewer.volume.Volume])[0]\n" )
    fp.write ( "q = qscores.CalcQ2pn ( mol, 'All', dmap, %f, useOld=False, numProc=%d, chimeraPath='%s', closeMap=True )\n" % ( gSigma, numProc, chiBin) )
    fp.write ( 'print "Q--==>%.4f" % q\n' )
    fp.close ()

    args = [chiBin, "--nogui", "--silent", "--nostatus", mapPath, molPath, scriptPath]
    print args

    fout = open ( inPath + "/_mapqScript_out.txt", "w" )
    ferr = open ( inPath + "/_mapqScript_err.txt", "w" )
    from subprocess import Popen
    p = Popen(args, stdout=fout, stderr=ferr, cwd=inPath)
    p.wait()
    fout.close()
    ferr.close()

    #print "output:"
    Q = -1.0
    fp = open ( inPath + "/_mapqScript_out.txt" )
    for l in fp :
        #print " :-: ", l,
        qs = l.split("--==>")
        if len(qs) == 2 :
            Q = float ( qs[1] )
            print " - found q: %.4f" % Q

    return Q



def CalcSseQ ( fitf ) :

    print ("")
    print " - calc sseQ for [%s]" % fitf

    #scriptPath = os.path.dirname(os.path.realpath(__file__))
    scriptPath = inPath + "/_mapqScript.py"
    print (" - creating script -> %s" % scriptPath)

    molPath = fitsPath + "/" + fitf

    fp = open ( scriptPath, "w" )
    fp.write ( "import chimera\n" )
    fp.write ( "import VolumeViewer\n" )
    fp.write ( "from Segger import qscores\n" )
    fp.write ( "from Segger import mmcif\n" )
    #fp.write ( "dmap = VolumeViewer.open_volume_file ( '%s', 'ccp4')[0]\n" % mapPath )
    #fp.write ( "mol = chimera.openModels.open ( '%s', type='PDB' )[0]\n" % molPath )
    fp.write ( "mol = chimera.openModels.list(modelTypes = [chimera.Molecule])[0]\n" )
    fp.write ( "dmap = chimera.openModels.list(modelTypes = [VolumeViewer.volume.Volume])[0]\n" )
    fp.write ( "sseq, sseqmax = qscores.sseQscores ( mol, dmap, sigma=3.0 )\n"  )
    fp.write ( 'print "Q--==>%.4f" % sseq\n' )
    fp.close ()

    args = [chiBin, "--nogui", "--silent", "--nostatus", mapPath, molPath, scriptPath]
    print args

    fout = open ( inPath + "/_mapqScript_out.txt", "w" )
    ferr = open ( inPath + "/_mapqScript_err.txt", "w" )
    from subprocess import Popen
    p = Popen(args, stdout=fout, stderr=ferr, cwd=inPath)
    p.wait()
    fout.close()
    ferr.close()

    log = False
    if log : print "output:"
    Q = -1.0
    fp = open ( inPath + "/_mapqScript_out.txt" )
    for l in fp :
        if log : print " :-: ", l,
        qs = l.split("--==>")
        if len(qs) == 2 :
            Q = float ( qs[1] )
            print " - found q: %.4f" % Q

    return Q


def fitFilesQ ( fpath ) :


    outf = "ez_out_Q2.txt"
    #fout = open ( "/Users/greg/Desktop/ez_out.txt", "w" )
    fout = open ( outf, "w" )
    fout.write ( "kda all\tkda af55\tZ-score\tCCM\tQ2\tfit#\tseqLen0\tseqLen1\tname 1\tname 2\n" )
    fout.close()

    ati = 0
    for fname in os.listdir(fpath) :
        ati += 1

        print fname
        ts = fname.split ("_")

        kdaRange, Z, ccm, fitNum = ts[0].split("-"), float(ts[4])/100.0, float(ts[7]), int(ts[9])

        fname_, at = "", 11
        seqLen0, k0 = None, None
        seqLen1, k1 = None, None
        while 1 :
            l = ts[at]
            if l [len(l)-3:len(l)] == "kDa" :
                seqLen0,k0 = l[0:-3].split("-")
            if l [len(l)-4:len(l)] == ".pdb" :
                l = l[0:-4]
                seqLen1,k1 = l[0:-3].split("-")
                break
            elif l == "af55" :
                pass
            else :
                fname_ += "_" + l

            at += 1
            if at >= len(ts) :
                break

        fname_ = fname_.strip("_")
        print " -----%d----- [%s]" % (ati, fname_)
        print " - %s:%s, Z %.2f, ccm %.3f, fit %d, seql %s->%s, kda %s->%s" % (kdaRange[0], kdaRange[1], Z, ccm, fitNum, seqLen0, seqLen1, k0, k1)


        if 1 and fitNum == 1 :

            #if fname != "00020-00040_kDa__z_00337__ccm_0.266__01__B9Q7B8_TOXGV__emp24gp25lp24-family-protein__259-30kDa__af55__191-23kDa.pdb" :
            #    continue

            Q = -1.0
            try :
                Q = CalcQ ( fname )
                #Q = CalcSseQ ( fname )
            except :
                print " - could not get Q..."

            fout = open ( outf, "a" )
            #fout.write ( "%s\t%s\t%.2f\t%.3f\t%d\t%s\t%s\t%s\t%s\n" % (kdaRange[0], kdaRange[1], Z, ccm, fitNum, seqLen0, seqLen1, k0, k1) )
            fout.write ( "%s\t%s\t%.2f\t%.3f\t%.4f\t%d\t%s\t%s\t%s\t%s\n" % (k0, k1, Z, ccm, Q, fitNum, seqLen0, seqLen1, fname_, fname) )
            fout.close()

            #break



def fitFiles ( fpath ) :


    #fout = open ( "/Users/greg/Desktop/ez_out.txt", "w" )
    #toF = os.path.split(fpath)[0] + "/ez_out_3.txt"

    fpath0, fpath1 = os.path.split(fpath)
    toF = fpath0 + "/" + fpath1 + ".txt"
    print ( " -> %s" % toF )


    print ( " -> %s" % toF )
    fout = open ( toF, "w" )
    fout.write ( "mod_seqLen\tmod_kDa\tseqLen\tkDa\tZ-score\tCCM\tfit#\tname 1\tname 2\n" )

    percs, seqLens, modLens, seqWts, modWts = [], [], [], [], []

    for fname in os.listdir(fpath) :

        #print fname
        ts = fname.split ("_")

        # 00020-00040_kDa__z_00269__ccm_0.202__04__B9Q7C5_TOXGV__sphingolipid-4-desaturase__424-40kDa__af55__319-39kDa
        # 00000-00020_kDa__z_00165__ccm_0.811__04__A0A125YQH5_TOXGV__putative-transmembrane-protein__260-31kDa__af55__10-1kDa

        kdaRange, Z, ccm, fitNum = ts[0].split("-"), float(ts[4])/100.0, float(ts[7]), int(ts[9])

        fname_, at = "", 11
        seqLen0, k0 = None, None
        seqLen1, k1 = None, None
        while 1 :
            l = ts[at]
            if l [len(l)-3:len(l)] == "kDa" :
                seqLen0,k0 = l[0:-3].split("-")
            if l [len(l)-4:len(l)] == ".pdb" :
                l = l[0:-4]
                seqLen1,k1 = l[0:-3].split("-")
                break
            elif l == "af55" :
                pass
            else :
                fname_ += "_" + l

            at += 1
            if at >= len(ts) :
                break

        fname_ = fname_.strip("_")
        #print " - [%s]" % fname_
        #print " - %s:%s, Z %.2f, ccm %.3f, fit %d, seql %s->%s, kda %s->%s" % (kdaRange[0], kdaRange[1], Z, ccm, fitNum, seqLen0, seqLen1, k0, k1)

        #if fitNum == 1 and float(k0) <= 161 :
        if fitNum == 1 :
            perc = float (seqLen1) / float ( seqLen0 )
            percs.append ( perc )
            seqLens.append ( int(seqLen0) )
            modLens.append ( int(seqLen1) )
            seqWts.append ( float(k0) )
            modWts.append ( float(k1) )
            #fout.write ( "%s\t%s\t%.2f\t%.3f\t%d\t%s\t%s\t%s\t%s\n" % (k0, k1, Z, ccm, fitNum, seqLen0, seqLen1, fname_, fname) )
            fout.write ( "%s\t%s\t%s\t%s\t%.3f\t%.3f\t%d\t%s\t%s\n" % (seqLen1, k1, seqLen0, seqLen1, Z, ccm, fitNum, fname_, fname) )

    fout.close()

    inPath = os.path.split(fpath)[0]
    h = numpy.histogram ( percs, bins=10, range=None )
    fo = open ( fpath0 + "/ez_out_2_percs_10_1217.txt", "w" )
    for hi in range(len(h[0])) :
        fo.write ( "%.0f\t%d\n" % (h[1][hi], h[0][hi]) )
    fo.close()

    h = numpy.histogram ( seqWts, bins=20, range=None )
    fo = open ( fpath0 + "/ez_out_2_seqWt_20_1217.txt", "w" )
    for hi in range(len(h[0])) :
        fo.write ( "%.0f\t%d\n" % (h[1][hi], h[0][hi]) )
    fo.close()

    h = numpy.histogram ( modWts, bins=20, range=(min(seqWts), max(seqWts)) )
    fo = open ( fpath0 + "/ez_out_2_modWt_20_1217.txt", "w" )
    for hi in range(len(h[0])) :
        fo.write ( "%.0f\t%d\n" % (h[1][hi], h[0][hi]) )
    fo.close()



def fitFiles2 ( fpath ) :


    #fout = open ( "/Users/greg/Desktop/ez_out.txt", "w" )
    fpath0, fpath1 = os.path.split(fpath)
    toF = fpath0 + "/" + fpath1 + ".txt"
    print ( " -> %s" % toF )

    fout = open ( toF, "w" )
    fout.write ( "mod_seqLen\tmod_kDa\tseqLen\tkDa\tZ-score\tCCM\tfit#\tname 1\tname 2\n" )

    percs, seqLens, modLens = [], [], []

    for fname in os.listdir(fpath) :

        #print fname
        ts = fname.split ("__")

        # ccm_0.166813__ccm_0.166813__cc_0.193314__30__B9PJD4_TOXGV__tubulin-alpha-chain__453__52.99kDa__af55__432__51.09kDa__.pdb
        # ['ccm_0.166813', 'ccm_0.166813', 'cc_0.193314', '30', 'B9PJD4_TOXGV', 'tubulin-alpha-chain', '453', '52.99kDa', 'af55', '432', '51.09kDa', '.pdb']

        s1, score1 = ts[0].split("_")
        sCCm, CCm = ts[1].split("_")
        sCC, CC = ts[2].split("_")
        fitNum = int ( ts[3] )
        idStr = ts[4]
        molName = ts[5]
        seqLen = int(ts[6])
        seqWt = float(ts[7][0:-3])
        afThr = ts[8]
        afLen = int(ts[9])
        afWt = float(ts[10][0:-3])

        if 1 or fitNum == 1 :
            perc = afWt / seqWt
            percs.append ( perc )
            seqLens.append ( seqLen )
            modLens.append ( afLen )
            #fout.write ( "%s\t%s\t%.2f\t%.3f\t%d\t%s\t%s\t%s\t%s\n" % (k0, k1, Z, ccm, fitNum, seqLen0, seqLen1, fname_, fname) )
            fout.write ( "%d\t%.3f\t%d\t%.3f\t%s\t%s\t%s\t%d\t%s\t%s\n" % (afLen, afWt, seqLen, seqWt, score1, CCm, CC, fitNum, idStr, fname) )

    fout.close()



def fitFiles3 ( fpath=None, W=200, minZ=3.0 ) :

    from shutil import copyfile

    if fpath == None :
        if len ( sys.argv ) >= 2 :
            fpath = sys.argv[1]
            print sys.argv[1]
            if os.path.isdir ( fpath ) :
                print " - first param is a path"
            else :
                print " - first param should be a path to the fit files"
                return
        else :
            print " - no path specified"
            return

        print " - path specified: %s" % fpath
        fpath = os.path.abspath ( fpath )
        print " - absolute: %s" % fpath
        #return

        if len ( sys.argv ) >= 3 :
            try :
                W = int ( sys.argv[2] )
            except :
                print " - second param should be window size, e.g. 200"
                return
            print " - window size specified as %d" % W
        else :
            W = 200


    #fout = open ( "/Users/greg/Desktop/ez_out.txt", "w" )
    fpath0, fpath1 = os.path.split(fpath)
    fpath0, fpath1 = os.path.split(fpath0)

    toF = os.path.join ( fpath0, fpath1 + ".txt" )
    toF2 = os.path.join ( fpath0, fpath1 + "_eZ_%d.txt" % W )
    toF3 = os.path.join ( fpath0, fpath1 + "_eZ_%d_above_%.1f" % (W,minZ) )
    toF4 = os.path.join ( fpath0, fpath1 + "_eZ_%d_above_%.1f.txt" % (W,minZ) )

    print ( " all fits -> %s" % toF )
    print ( " all fits with eZ -> %s" % toF2 )
    print " all fits by eZ -> %s " % (toF3)
    print " all fits with eZ > %.1f -> %s " % (minZ, toF4)

    if os.path.isdir ( toF3 ) :
        print " - found %s" % toF3
    else :
        print " - making %s" % toF3
        os.mkdir(toF3)

    fout = open ( toF, "w" )
    fout.write ( "mod_seqLen\tmod_kDa\tseqLen\tkDa\tZ-score\tCC-mean\tCC\tfit#\tid1\tnfile\n" )

    percs, seqLens, modLens = [], [], []
    alls = []

    for fname in os.listdir(fpath) :

        #print fname
        ts = fname.split ("__")

        # lz_017.806__z_003.352__ccm_0.260757__01__A0A0F7V5W4_TOXGV__palmitoyltransferase__356__33.71kDa__af55__245__30.09kDa__.pdb
        # ['lz_017.806', 'z_003.352', 'ccm_0.260757', '01', 'A0A0F7V5W4_TOXGV', 'palmitoyltransferase', '356', '33.71kDa', 'af55', '245', '30.09kDa', '.pdb']

        #z_004.202__ccm_0.287509__cc_0.631451__01__B6K916_TOXGV__MORN_repeat_containing_protein__429__49.00kDa__af55__258__30.74kDa__
        # ['z_004.202', 'ccm_0.287509', 'cc_0.631451', '01', 'B6K916_TOXGV', 'MORN_repeat_containing_protein', '429', '49.00kDa', 'af55', '258', '30.74kDa', '']

        #slZ, lZ = ts[0].split("_")
        sZ, Z = ts[0].split("_")
        sCCm, CCm = ts[1].split("_")
        sCC, CC = ts[2].split("_")
        fitNum = int ( ts[3] )
        idStr = ts[4]
        molName = ts[5]
        seqLen = int(ts[6])
        seqWt = float(ts[7][0:-3])
        afThr = ts[8]
        afLen = int(ts[9])
        afWt = float(ts[10][0:-3])

        if fitNum == 1 :
            perc = afWt / seqWt
            percs.append ( perc )
            seqLens.append ( seqLen )
            modLens.append ( afLen )
            #fout.write ( "%s\t%s\t%.2f\t%.3f\t%d\t%s\t%s\t%s\t%s\n" % (k0, k1, Z, ccm, fitNum, seqLen0, seqLen1, fname_, fname) )
            fout.write ( "%d\t%.3f\t%d\t%.3f\t%s\t%s\t%s\t%d\t%s\t%s\n" % (afLen, afWt, seqLen, seqWt, Z, CCm, CC, fitNum, idStr, fname) )

            alls.append ( [afLen, afWt, afThr, seqLen, seqWt, float(CCm), fitNum, molName, idStr, fname] )

    fout.close()

    if 1 :
        print " - sorting %d entries" % len(alls)
        alls.sort ( reverse=False, key=lambda x: x[0] )

        allc = []
        for afLen, afWt, afThr, seqLen, seqWt, CCm, fitNum, molName, idStr, fname in alls :
            allc.append ( CCm )

        fout = open ( toF2, "w" )
        fout.write ( "mod_seqLen\tmod_kDa\tseqLen\tkDa\tCC-mean\teZ\tid1\tfile\n" )

        fout4 = open ( toF4, "w" )
        fout4.write ( "eZ\tCC-mean\tmod_seqLen\tmod_kDa\tseqLen\tkDa\tid1\tfile\n" )

        print " - writing and calculating local Z"
        for i, ts in enumerate ( alls ) :
            start_i = max ( 0, i-W )
            end_i = min ( len(alls), i+W )
            ccs = allc[start_i:end_i]
            mean, stdev = numpy.mean(ccs), numpy.std(ccs)

            afLen, afWt, afThr, seqLen, seqWt, CCm, fitNum, molName, idStr, fname = ts
            eZ = (float(CCm) - mean) / stdev
            fout.write ( "%d\t%.3f\t%d\t%.3f\t%.3f\t%.03f\t%s\t%s\n" % (afLen, afWt, seqLen, seqWt, CCm, eZ, idStr, fname) )

            if eZ >= minZ :
                fromf = os.path.join ( fpath, fname )
                tof = os.path.join ( toF3, "eZ_%07d__ccm__%.6f__%d__%s__%s__%d__%.3fkDa__%s__%d__%.3fkDa.pdb" % (int(eZ*1000), CCm, fitNum, idStr, molName, seqLen, seqWt, afThr, afLen, afWt) )
                if 0 :
                    print " --< %s" % fromf
                    print " --> %s" % tof
                    return
                copyfile ( fromf, tof )
                fout4.write ( "%.3f\t%.6f\t%d\t%.3f\t%d\t%.3f\t%s\t%s\n" % (eZ, CCm, afLen, afWt, seqLen, seqWt, idStr, fname) )

        fout.close()
        fout4.close()




def getAfMods () :

    gotIds = {}
    notIds = {}
    if 1 :
        if os.path.isfile ( inPath + "prots_af_stats/prots_af.txt" ) :
            fp1 = open ( inPath + "prots_af_stats/prots_af.txt" )
            for l in fp1 :
                ts = l.strip().split()
                if len(ts) == 5 :
                    id1, seqLen, modLen, modfLen, perc = ts
                    gotIds[id1] = ts
                elif len(ts) == 1 :
                    id1 = ts[0]
                    notIds[id1] = None
                else :
                    print " - blank line?"

            fp1.close()

    print ( " - got %d, no %d" % ( len(gotIds), len(notIds) ) )

    from Segger.fit import CopyAfMol

    #fp1 = open ( inPath + "prots_af.txt", "w" )
    #fp1.close()

    prots = getProts ()
    for N, P in enumerate ( prots ) :

        path = inPath + "prots_af/" + "%s.pdb" % P.id1
        #print "%d / %d -- %s" % (N+1, len(prots), path)

        if P.id1 in gotIds :
            #print ( " -- have it " )
            continue

        if P.id1 in notIds :
            #print ( " -- no pred " )
            continue

        if not os.path.isfile ( path ) :
            print " -- no file for %s" % P.id1
            continue

        try :
            mol = chimera.openModels.open ( path )[0]
        except :
            fp1 = open ( inPath + "prots_af.txt", "a" )
            fp1.write ( "%s\n" % (P.id1) )
            fp1.close()
            continue


        maf = CopyAfMol ( mol, chimera.Xform.identity(), bCutoff=55.0 )
        fp1 = open ( inPath + "prots_af.txt", "a" )
        fp1.write ( "%s\t%d\t%d\t%d\t%.2f\n" % ( P.id1, len(P.seq), len(mol.residues), len(maf.residues), perc ) )
        fp1.close()

        chimera.openModels.close ( [mol] )


    fp1.close()



def getAfMods2 () :

    percs, seqLens, modLens = [], [], []
    seqWeights, modWeights = [], []

    genSeq, totSeq, totModSeq = 0, 0, 0

    from shutil import copy2

    prots = getProts ()
    idProt = {}
    for N, P in enumerate ( prots ) :
        idProt[P.id1] = P


    fp1 = open ( inPath + "prots_af_stats/prots_af.txt" )
    gotIds = {}
    notIds = {}
    allIds = {}
    ati = 1
    for l in fp1 :
        ts = l.strip().split()
        if len(ts) == 5 :
            id1, seqLen, modLen, modfLen, perc = ts
            gotIds[id1] = ts
            allIds[id1] = None

            perc = float (modfLen) / float (seqLen) * 100.0
            percs.append ( perc )
            seqLens.append ( float(seqLen) )
            modLens.append ( float(modfLen) )
            seqWeights.append ( float(seqLen) * 0.1246 + 0.8389 )
            modWeights.append ( float(modfLen) * 0.1246 + 0.8389 )

            totSeq += float(seqLen)
            totModSeq += float(modfLen)

            if int(seqLen) > len(P.seq) :
                print " -> ? %d - %s - seq len %d/%d, model %d/%d" % (ati, P.id1, len(P.seq), int(seqLen), int(modLen), int(modfLen))
            #if P.seqLen < int(modfLen) :
            #    print " -> ?? %s - seq len %d, model %d" % (P.id1, P.seqLen, int(modfLen))
            ati += 1

            P = idProt[id1]
            if 0 :
                fromp = inPath + "prots_af/" + "%s.pdb" % id1
                top = inPath + "prots_af_/%s.pdb" % P.fname
                print " -> %s" % top
                copy2 ( fromp, top )

        elif len(ts) == 1 :
            id1 = ts[0]
            notIds[id1] = None
            allIds[id1] = None

        else :
            print (" - blank")

    numNF = 0
    for N, P in enumerate ( prots ) :
        genSeq += P.seqLen
        if P.id1 in gotIds or P.id1 in notIds :
            continue
        else :
            print " -? %s" % P.id1
            numNF += 1


    fp1.close()
    print ( " - got %d, not %d / %d | genome %d" % ( len(gotIds), len(notIds), len(allIds), len(prots) ) )

    print " - %d of genome %d -- %f%%" % ( totModSeq, genSeq, float(totModSeq)/float(genSeq)%100.0 )
    print " - %d in models %d -- %f%%" % ( totModSeq, totSeq, float(totModSeq)/float(totSeq)%100.0 )

    h = numpy.histogram ( percs, bins=10, range=None )
    fo = open ( inPath + "prots_af_stats/prots_af_percs_10.txt", "w" )
    for hi in range(len(h[0])) :
        fo.write ( "%.0f\t%d\n" % (h[1][hi], h[0][hi]) )
    fo.close()

    h = numpy.histogram ( seqLens, bins=20, range=None )
    fo = open ( inPath + "prots_af_stats/prots_af_seqWt_20.txt", "w" )
    for hi in range(len(h[0])) :
        fo.write ( "%.0f\t%d\n" % (h[1][hi], h[0][hi]) )
    fo.close()

    h = numpy.histogram ( modLens, bins=20, range=(min(seqLens), max(seqLens)) )
    fo = open ( inPath + "prots_af_stats/prots_af_modWt_20.txt", "w" )
    for hi in range(len(h[0])) :
        fo.write ( "%.0f\t%d\n" % (h[1][hi], h[0][hi]) )
    fo.close()



class Prot () :
    def __init__ (self) :
        descr = ""


def slugify ( str ) :
    from re import sub
    return sub('[^0-9a-zA-Z]+', '_', str)


def getProts () :

    fp = open ( "/Volumes/S1/toxop_/uniprot-proteome_UP000002226.fasta" )
    li = 0

    P = None
    types = {}
    proteins = []

    for l in fp :
        li += 1
        # if li > 100 : break
        if l[0] == ">" :
            if P != None :
                P.seqLen = len(P.seq)
                P.fname = "%s__%s__%d" % (P.id2, slugify(P.descr), P.seqLen)
            P = Prot()
            proteins.append ( P )
            ts = l.split()
            P.id = ts[0][1:]
            P.type, P.id1, P.id2 = P.id.split("|")
            P.header = l
            if P.type in types :
                types[P.type] += 1
            else :
                types[P.type] = 1
            desc = []
            ti = 1
            while 1 :
                if ts[ti][0:3] == "OS=" :
                    break
                desc.append ( ts[ti] )
                ti += 1
                if ti >= len(ts) :
                    break

            os = ts[ti]
            P.descr = " ".join(desc)
            P.seq = ""

        elif len(l.strip()) == 0 :
            print ( " -- %d" % li )
            if P != None :
                P.seqLen = len(P.seq)
                descr = P.descr.replace ( " ", "_" )
                P.fname = "%s__%s__%d" % (P.id2, descr, P.seqLen)
                #print ( P.seqLen, P.seq )
                #print ( "%4d %s %s |%s|" % (P.seqLen, P.type, P.id2, P.descr) )
        else :
            P.seq += l.strip()

    fp.close()

    if P != None :
        P.seqLen = len(P.seq)
        descr = P.descr.replace ( " ", "_" )
        P.fname = "%s__%s__%d" % (P.id2, descr, P.seqLen)

    print ( "\n" )
    print ( " - %d proteins" % len(proteins) )
    print ( "\n" )

    return proteins


# microtubule_crop__ez_fits_res700_step10_8rot
# lz_001.186__z_001.308__ccm_0.767795__01__B6KGT4_TOXGV__uncharacterized-protein__189__23.24kDa__af55__3__0.26kDa__

# conoid_inside_crop__ez_fits_res700_step15_6rot
# lz_000.589__z_001.121__ccm_0.695035__03__A0A125YV58_TOXGV__uncharacterized-protein__81__9.59kDa__af55__16__1.94kDa__

#fitFilesQ ( fitsPath )
#fitFiles ( "/Users/greg/_data/toxop_/CF2/by_ccm" )
#fitFiles ( "/Volumes/S1/toxop_/SPM2/by_ccm" )

#fitFiles2 ( "/Users/greg/Box Sync/_data/toxop_/SPM3/microtubule_crop__ez_fits_res700_step10_8rot__B9PJD4__tubulin/by_z-score" )
#fitFiles2 ( "/Users/greg/Box Sync/_data/toxop_/SPM3/microtubule_crop__ez_fits_res700_step10_8rot__A0A0F7UXJ4__nucleoredoxin/by_ccm" )
#fitFiles3 ( "/Volumes/S1/toxop_/SPM3/microtubule_crop__ez_fits_res700_step15_6rot/by_ccm" )
#fitFiles3 ( "/Volumes/S1/toxop_/CF3/conoid_inside_crop__ez_fits_res700_step15_6rot/by_ccm" )

#fitFiles3 ( "/Volumes/S1/toxop_/SPM3/microtubule_crop__ez_fits_res700_step15_6rot/by_ccm" )

#fitFiles3 ( "/Users/greg/_data/toxop__/SPM3/microtubule_crop__ez_fits_res700_step10_10rot/by_ccm" )
#fitFiles3 ( "/sdf/home/g/gregp/g24/toxop/maps2/microtubule_crop__ez_fits_res700_step10_10rot/by_ccm", W=400 )

fitFiles3 ()

#getAfMods ( )
#getAfMods2 ( )

#
