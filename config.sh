# ACTIVE SCRIPTS
# Select which scripts to run.
# If CALCULATE_TOLERANCE is set to True, runToleranceAnalysis.py will be executed.
# If FIT_WITH_TOLERANCE is set to True, fitWithTolerance.py will be executed.
# If OPEN_PLOTS is set to True, the final plots will be opened (this may not work on all operative systems).
CALCULATE_TOLERANCE=False
FIT_WITH_TOLERANCE=False
OPEN_PLOTS=False

# ARGUMENTS
# The arguments are identical to the one described in the README.md file.
# Suggested proton MV range: 416 - 555
# Suggested helium MV range: 372 - 555
PRIMARY="Proton"
FNAME="AMSGSHL"
BARR=True
DATA="Data/SettingsFiles/${PRIMARY}EnergyPerNucleonData.json"
OUTDIR="Outdir"
if [ ${PRIMARY} = "Proton" ]; then
  SMIN="416"
elif [ ${PRIMARY} = "Helium" ]; then
  SMIN="372"
else
  SMIN="0"
fi
SMAX="555"
