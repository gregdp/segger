# segger
Plugin to UCSF Chimera for segmenting Cryo-EM density maps, fitting atomic models to maps, and various other map-model related utilities.

To Install:

1. Please note that <a href="https://www.cgl.ucsf.edu/chimera/">UCSF Chimera</a> is required to use ModelZ
2. Download latest version of Segger file from the <a href="https://github.com/gregdp/segger/tree/master/download">download</a> folder
3. In a terminal, navigate to where the file was downloaded, then run the following commands:
* unzip Segger_2_0.zip
* cd Segger_2_0
* python install.py [path to Chimera]

  On Windows, you may use the Chimera python, e.g. [path to Chimera]/bin/python install.py [path to Chimera]

To Run:
1. (Re)start Chimera*
2. Start ModelZ: Tools -> Volume Data -> ModelZ
3. See [Tutorial](https://github.com/gregdp/modelz/blob/master/tutorials/Tutorial-ModelZ.pdf)

\* On Mac OS, an error message may be shown at first, see [here](https://www.santoshsrinivas.com/disable-gatekeeper-in-macos-sierra/) for solution.

