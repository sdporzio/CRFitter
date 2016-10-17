# Case settings
export DATA="Data/SettingsFiles/HeliumTestEnergyPerNucleonData.json"
export PRIMARY="Helium"
export FNAME="AMSGSHL"
export IFBARR="--Barr" # Else leave empty
export OUTDIR="TestOutdir"
export TOLERANCE="${OUTDIR}/${FNAME}/${PRIMARY}/ToleranceMethod2/ToleranceData/Method2DeltaChi.dat"
# export SMIN="416"
# export SMAX="555"
export SMIN="385"
export SMAX="564"

# Global Settings
CALC_TOL=False
FIT_TOL=True
OPEN_PLOT=True

# Function
if ${CALC_TOL}; then
  export CMD="python Scripts/runToleranceAnalysis.py --data ${DATA} --primary ${PRIMARY} --smin ${SMIN} --smax ${SMAX} --fname ${FNAME} ${IFBARR} --outdir ${OUTDIR}";
  echo; echo "> ${CMD}"; echo;
  eval ${CMD};
fi

if ${FIT_TOL}; then
  export CMD="python Scripts/fitWithTolerance_Format.py --tolerance ${TOLERANCE} --data ${DATA} --primary ${PRIMARY} --smin ${SMIN} --smax ${SMAX} --fname ${FNAME} ${IFBARR} --outdir ${OUTDIR}";
  echo; echo "> ${CMD}"; echo;
  eval ${CMD};
fi

if ${OPEN_PLOT}; then
  export CMD="open ${OUTDIR}/${FNAME}/${PRIMARY}/FinalPlots/*.pdf";
  echo; echo "> ${CMD}"; echo;
  eval ${CMD};
fi
