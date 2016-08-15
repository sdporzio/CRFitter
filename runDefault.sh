# Settings
export DATA="Data/SettingsFiles/ProtonEnergyPerNucleonData.json"
export PRIMARY="Proton"
export FNAME="GSHL"
export IFBARR="--Barr" # Else leave empty
export OUTDIR="Outdir"
export TOLERANCE="${OUTDIR}/${FNAME}/${PRIMARY}/ToleranceMethod2/ToleranceData/Method2DeltaChi.dat"
# export SMIN="385"
# export SMAX="564"
export SMIN="416"
export SMAX="555"

CALC_TOL=True
FIT_TOL=True
OPEN_PLOT=True

# Execution
if ${CALC_TOL}; then
  export CMD="python Scripts/runToleranceAnalysis.py --data ${DATA} --primary ${PRIMARY} --smin ${SMIN} --smax ${SMAX} --fname ${FNAME} ${IFBARR} --outdir ${OUTDIR}"
  echo; echo "> ${CMD}"; echo
  eval ${CMD}
fi

if ${FIT_TOL}; then
  export CMD="python Scripts/fitWithTolerance.py --tolerance ${TOLERANCE} --data ${DATA} --primary ${PRIMARY} --smin ${SMIN} --smax ${SMAX} --fname ${FNAME} ${IFBARR} --outdir ${OUTDIR}"
  echo; echo "> ${CMD}"; echo
  eval ${CMD}
fi

if ${OPEN_PLOT}; then
  export CMD="open ${OUTDIR}/${FNAME}/${PRIMARY}/FinalPlots/*.pdf"
  echo; echo "> ${CMD}"; echo
  eval ${CMD}
fi
