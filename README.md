# Segger
Plugin to <a href="https://www.cgl.ucsf.edu/chimera/">UCSF Chimera</a> for segmenting Cryo-EM density maps, fitting atomic models to maps, and various other map-model related utilities.

For technical details (and please cite :)
* (2010) Quantitative analysis of cryo-EM density map segmentation by watershed...<a href="https://pubmed.ncbi.nlm.nih.gov/20338243/" target="_blank">PubMed</a>
* (2012) Comparison of Segger and...<a href="https://pubmed.ncbi.nlm.nih.gov/22696409/" target="_blank">Pubmed</a>
* (2016) Resolution and Probabilistic Models of...<a href="https://pubmed.ncbi.nlm.nih.gov/26743049/" target="_blank">Pubmed</a>
* (2019) Segmentation and Comparative Modeling in...<a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6853598/" target="_blank">Pubmed</a>
* (2020) Measurement of atom resolvability in cryo-EM maps with Q-scores</a> <a href="https://www.nature.com/articles/s41592-020-0731-1" target="_blank">Nature Methods<a href="https://www.biorxiv.org/content/10.1101/722991v1" target="_blank">BioRXiv</a>

To install the latest version:

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

To Run:
1. (Re)start Chimera*
2. Start Segger: Tools -> Volume Data -> Segment Map
3. See [more details and tutorials](https://cryoem.slac.stanford.edu/ncmi/resources/software/segger)

\* On Mac OS, an error message may be shown on first run after installing, see [here](https://www.santoshsrinivas.com/disable-gatekeeper-in-macos-sierra/) for solution.

