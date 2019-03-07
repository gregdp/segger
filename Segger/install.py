import sys, os, shutil

#li = sys.argv.index ( "install.py" )
#print li


if len(sys.argv) != 2 :
    print ""
    print "Please add the path where Chimera is installed, e.g.:"
    print "   python install.py /home/greg/applications/Chimera"
    print ""

    exit()
    

print ""

opath1 = os.path.join ( sys.argv[1], "Contents" )
opath1 = os.path.join ( opath1, "Resources" )
opath1 = os.path.join ( opath1, "share" )

opath2 = os.path.join ( sys.argv[1], "share" )

didInstall = False

for opath in [opath1, opath2] :

    if os.path.isdir( opath ) :
        opath = os.path.join ( opath, "Segger" )
    
        if os.path.isdir( opath ) :
            print " - removing previous Segger:", opath
            try :
                shutil.rmtree(opath)
            except :
                pass
        
        #print " - copying from:", os.getcwd()
        print " - copying . ->", opath
    
        try :
            shutil.copytree ( os.getcwd(), opath )  
            didInstall = True
        except :
            print "Could not copy to ", opath
            print " 1. please check if you have write access"
            print " 2. try with sudo python install.py <path to Chimera>"
            print ""
            break

        didInstall = True

if didInstall :

    print ""
    print "Installation complete."
    print ""
    print "To use:"
    print " 1. Please restart Chimera."
    print " 2. Select Tools -> Volume Data -> Segment Map"
    print ' 3. Please note that on Mac OS, you may see the message "Chimera is damaged and cannot be opened." Please see the following link for the solution: https://www.santoshsrinivas.com/disable-gatekeeper-in-macos-sierra/'
    print ' 4. More info: https://cryoem.slac.stanford.edu/ncmi/resources/software/Segger'
    print ' 5. Tutorials: https://github.com/gregdp/segger/tree/master/tutorials'
    print ' 6. For other questions/comments/suggestions, please contact gregp@slac.stanford.edu'
    print ""
    
    #wh = os.path.join ( os.getcwd(), "install.html" )
    #import webbrowser
    #webbrowser.open( 'file://' + wh, new=2)


else :
    print ""
    print 'Chimera not found in "' + sys.argv[1] + '"'
    print " 1. please check the path"
    print " 2. remember you can auto-complete while typing the path with the tab key"
    print " 3. if issue persists, please report to gregp@slac.stanford.edu"
    print ""


