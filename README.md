# Segger
Plugin to UCSF Chimera for segmenting Cryo-EM density maps, fitting atomic models to maps, and various other map-model related utilities.

To Install:

1. Please note that <a href="https://www.cgl.ucsf.edu/chimera/">UCSF Chimera</a> is required to use ModelZ
2. Download latest version of Segger from the <a href="https://github.com/gregdp/segger/tree/master/download">download</a> folder
3. In a terminal, navigate to where the file was downloaded, then run the following commands:
* unzip Segger_2_0.zip
* cd Segger_2_0
* python install.py [path to Chimera]

Note that on Windows, you may use the python bundled with Chimera itself, so the third command would be
* [path to Chimera]/bin/python install.py [path to Chimera]

To Run:
1. (Re)start Chimera*
2. Start Segger: Tools -> Volume Data -> Segment Map
3. See [more details and tutorials](https://cryoem.slac.stanford.edu/ncmi/resources/software/segger)

\* On Mac OS, an error message may be shown on first run after installing, see [here](https://www.santoshsrinivas.com/disable-gatekeeper-in-macos-sierra/) for solution.

