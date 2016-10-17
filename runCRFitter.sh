displayAndRun(){
  COLOR='\033[1;33m'
  DEFAULT='\033[0m'
  echo -e "${COLOR}-> ${1}${DEFAULT}";
  eval ${1};
}

source config.sh

if ${BARR}; then
  IFBARR="--Barr"
else
  IFBARR=""
fi
TOLERANCE="${OUTDIR}/${FNAME}/${PRIMARY}/ToleranceMethod2/ToleranceData/Method2DeltaChi.dat"

if ${CALC_TOL}; then
  CMD="python Scripts/runToleranceAnalysis.py --data ${DATA} --primary ${PRIMARY} --smin ${SMIN} --smax ${SMAX} --fname ${FNAME} ${IFBARR} --outdir ${OUTDIR}";
  displayAndRun "${CMD}"
fi

if ${FIT_TOL}; then
  CMD="python Scripts/fitWithTolerance.py --tolerance ${TOLERANCE} --data ${DATA} --primary ${PRIMARY} --smin ${SMIN} --smax ${SMAX} --fname ${FNAME} ${IFBARR} --outdir ${OUTDIR}";
  displayAndRun "${CMD}"
fi

if ${OPEN_PLOT}; then
  CMD="open ${OUTDIR}/${FNAME}/${PRIMARY}/FinalPlots/*.pdf";
  displayAndRun "${CMD}"
fi
