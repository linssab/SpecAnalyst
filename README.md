# SpecAnalyst
Python script for quick evaluation of X-Ray fluorescence spectra. The script reads and plot the spectrum. The spectrum can be cleaned (have its continuum calculated and stripped) and saved as a txt file. Net peak areas can be calculated by inputting the desired chemical element. Preference is given to KA macro until 35 KeV, when its changed to LA macro.<br>
This script was tested with AMPTEK(C) *.mca spectra and XRMC spectra. If the files you are trying to read are not supported or present issues, please open an issue.<br>
An XRMC testfile is available [here](../master/newtest_1.txt)

## Installation:

This script has external dependencies.<br>
Xraylib package is optional to run the script, but strongly recommended. You can find it available at [this link](https://github.com/tschoonj/xraylib/wiki)<br>
Numpy, matplotlib and scipy packages are necessary.<br>
If you lack any of the three packages mentioned, simply install them via pip through your terminal:<br>
<br>
`pip install numpy`<br>
`pip install matplotlib`<br>
`pip install scipy`<br>
<br>
## Execution:

### Plot spectrum:
In the terminal type:<br>
<br>
`python SpecAnalysis.py -s "filename"`<br>
<br>
The spectrum file must be in the same directory as the script.<br>
If you wish to plot the spectrum without the background approximation line, type:<br>
<br>
`python SpecAnalysis.py -s "filename" -nobg`<br>
<br>
By default, the energy axis is set according to a list of anchors already defined inside the script. Obviously this calibration parameter does not fit other spectrum. To read the internal calibration parameters (in the case of *.mca files) add the `-fs` command in the end.<br>
<br>
`python SpecAnalysis.py -s "filename" -fs`<br>
<br>
### Wipe background and save:
You can save the spectrum file minus the background approximation by typing:<br>
<br>
`python SpecAnalysis.py -s "filename" -w "output file name"`<br>
<br>
Note that if you add the `-nobg` command, the output file will contain the same values as the input, since the background will not be calculated.<br>

### Get a peak area:
The script can extract the net area of a given peak. The input must be the chemical element name, e.g. "Fe", "Cu", "Ti", "Pb" and so on.<br>
To plot the spectrum with the net area, type:<br>
<br>
`python SpecAnalysis.py -s "filename" -e "element"`<br>
<br>
To get rid of the fill area, add the `-nofill` command. This opstion can be used together with `-nobg`, to get the peak total area, instead of the net area.<br>
