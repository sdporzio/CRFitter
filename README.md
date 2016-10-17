# CRFitter

CRFitter is a fitter for cosmic rays which calculates a *tolerance* parameter to better determine uncertainties associated with the fit parameters. This *tolerance* method is inspired by the work of Stump *et al.*, 2001 and details can be found in their papers (https://arxiv.org/abs/hep-ph/0101051, https://arxiv.org/abs/hep-ph/0101032).

CRFitter depends on Python, ROOT (version 6.4.4 has been tested) and PyROOT and its operation is based on two scripts found in the *Scripts/* directory. `runToleranceAnalysis.py` analyzes the datasets contained in the *Data/* directory and determines a *tolerance* factor. `fitWithTolerance.py` uses this *tolerance* to determine the uncertainties on the parameters of the fit, which are then used to draw the uncertainty corridors. 'fitWithTolerance.py' must be obviously executed after `runToleranceAnalysis.py`.

Both scripts rely on JSON configuration files which are located in *Data/SettingFiles/*. The JSON file must contain the path to all the datasets to be included in the analysis, plus additional settings for creating plots.

Deviation plots with inflation factors for uncertainties in the calorimeter region, as shown in Barr *et al.*, 2005 (https://arxiv.org/abs/astro-ph/0611266), can be produced using the specific argument `--Barr`.
******************************************
Usage:
```
fitWithTolerance.py [-h] [--tolerance TOLERANCE] [--data DATA]
                       [--primary PRIMARY] [--fname FNAME] [--Barr]
                       [--emin EMIN] [--emax EMAX] [--smin SMIN]
                       [--smax SMAX] [--outdir OUTDIR]

optional arguments:
-h, --help            show this help message and exit

--tolerance TOLERANCE     DAT file containing delta chi values from which to calculate the tolerance. (default: None)

--data DATA       JSON file containing locations of all of the data you want to include in the fit. (default: None)

--primary PRIMARY     Name of primary you are performing the global fit on. (default: Proton

--fname FNAME     Name of fitting function you want to use for the global fits. The choices are: * GSHL * AMSGSHL (GSHL with AMS modification) * Simple (Power Law) * H3a (default: GSHL)

--Barr        Flag if performing paper a la Barr et al. (default: False)

--emin EMIN       Minimum energy value for fitter. (default: 0.5)

--emax EMAX       Maximum energy value for fitter. (default: 250000)

--smin SMIN       Minimum solar modulation value for which to accept all data. (default: 416)

--smax SMAX       Maximum solar modulation value for which to accept all data. (default: 555)

--outdir OUTDIR       Directory to save output. Script is designed to create subdirectories and delete older plots if it needs to, so set this carefully! If none provided, all plots will be saved in CWD. (default: None)
```
Example usage:
```
python Scripts/runToleranceAnalysis.py --data Data/SettingsFiles/ProtonEnergyPerNucleonData.json --primary Proton --smin 385 --smax 564 --fname GSHL --Barr --outdir Outdir

python Scripts/fitWithTolerance.py --tolerance Outdir/GSHL/Proton/ToleranceMethod2/ToleranceData/Method2DeltaChi.dat --data Data/SettingsFiles/ProtonEnergyPerNucleonData.json --primary Proton --smin 385 --smax 564 --fname GSHL --Barr --outdir Outdir
```
Since the script takes a large number of arguments, it is possible to use a bash parser (`runCRFitter.sh`) which uses a setting file called `config.sh`.
Arguments for the python script can be provided to the `config.sh` and CRFitter can ran with:
```
./runCRFitter.sh
```
******************************************
If you have any questions, queries or comments please contact the authors:

* salvatore.porzio@postgrad.manchester.ac.uk

* steven.wren@icecube.wisc.edu

The authors make no guarrentee of the behaviour, stability or bug-free-ness of this code. This code is provide "as is", use is at own risk.

Copyright (2016) Salvatore Davide Porzio, Steven Wren
******************************************
