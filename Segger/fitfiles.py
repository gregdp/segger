
import sys, os



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
    fout = open ( "ez_out.txt", "w" )
    fout.write ( "kda all\tkda af55\tZ-score\tCCM\tfit#\tseqLen0\tseqLen1\tname 1\tname 2\n" )

    for fname in os.listdir(fpath) :

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
        print " - [%s]" % fname_
        print " - %s:%s, Z %.2f, ccm %.3f, fit %d, seql %s->%s, kda %s->%s" % (kdaRange[0], kdaRange[1], Z, ccm, fitNum, seqLen0, seqLen1, k0, k1)

        fout.write ( "%s\t%s\t%.2f\t%.3f\t%d\t%s\t%s\t%s\t%s\n" % (k0, k1, Z, ccm, fitNum, seqLen0, seqLen1, fname_, fname) )

    fout.close()

fitFilesQ ( fitsPath )
#fitFiles ( "/Users/greg/_data/toxop_/CF2/by_ccm" )


#
