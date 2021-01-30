# Segger
Plugin to <a href="https://www.cgl.ucsf.edu/chimera/">UCSF Chimera</a> for segmenting 3D cryo-EM density maps, fitting atomic models to maps, and various other map-model related tools.

For details and tutorials, please see <a href="https://github.com/gregdp/segger/tree/master/tutorials">tutorials</a> folder.




To <strong>install</strong> the latest version:

1. First, <a href="https://www.cgl.ucsf.edu/chimera/download.html">download</a> and install Chimera. 
* Run it once before installing the plugin; on some platforms, e.g. MacOS, you may see a warning message on first run which you have to accept. This may prevent further issues after adding the plugin.
* On Windows, install to your home folder rather than to "Program Files". In the latter, the OS may not allow further modifications at a user level, i.e. adding this plugin.
2. <a href="https://github.com/gregdp/segger/tree/master/download">Download</a> latest version of Segger.
3. In a terminal, navigate to where the file was downloaded, then run the install script, e.g.:
* cd ~/Downloads
* unzip Segger_2_1.zip
* cd Segger_2_1
* python install.py ~/Desktop/Chimera.app
* (in the above, replace ~/Desktop/Chimera.app with the folder where Chimera was installed)
4. Note that on Windows, python may not be already installed; you may however use the python bundled with Chimera itself, e.g.:
* cd Downloads\Segger_2_4
* "C:\Users\greg\Chimera 1.14\bin\python.exe" install.py "C:\Users\greg\Chimera 1.14"
* (in the above, replace "C:\Users\greg\Chimera 1.14" with the directory where Chimera was installed)

To <strong>run</strong>:
1. (Re)start Chimera*
2. Start Segger: Tools -> Volume Data -> Segment Map
3. See [more details and tutorials](https://cryoem.slac.stanford.edu/ncmi/resources/software/segger)

\* On Mac OS, an error message may be shown on first run after installing, see [here](https://www.santoshsrinivas.com/disable-gatekeeper-in-macos-sierra/) for solution.

