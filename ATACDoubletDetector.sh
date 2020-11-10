#!/bin/bash

BAMBC="CB"
BCIDX=0
CELLIDX=8
ISCELLIDX=9

ARGUMENTS=()

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --bambc) BAMBC="$2"; shift 2 ;;
        --bcidx) BCIDX=$2; shift 2 ;;
        --cellidx) CELLIDX=$2; shift 2 ;;
        --iscellidx) ISCELLIDX=$2; shift 2 ;;
        --*) echo "Unknown argument: $1"; exit 1 ;;
        *)  ARGUMENTS+=($1); shift ;;
    esac
done

if [[ ${#ARGUMENTS[@]} -ne 6 ]]; then
	echo "Usage: ATACDoubletDetector.sh bamfile barcodemap chromosomelist repeatfilter outputdirectory scriptpath";
	echo "Options: --bambc     Bamfile attribute used for the barcode. (Default=\"CB\")";
	echo "         --bcidx     The column index of the CSV for barcode. (Default: 0)";
	echo "         --cellidx   The column index of the CSV for cellid. (Default: 8)";
	echo "         --iscellidx The index for determining cells (selecting values=1). (Default: 9)";
	exit 0;
fi

BAMFILE=${ARGUMENTS[0]}
BCMAP=${ARGUMENTS[1]}
CHRLIST=${ARGUMENTS[2]}
REPFILTER=${ARGUMENTS[3]}
OUTDIR=${ARGUMENTS[4]}
SCRIPTPATH=${ARGUMENTS[5]}

java -jar ${SCRIPTPATH}/snATACOverlapCounter.jar --bambc ${BAMBC} --bcidx ${BCIDX} --cellidx ${CELLIDX} --iscellidx ${ISCELLIDX} ${BAMFILE} ${BCMAP} ${CHRLIST} ${OUTDIR}

python3 ${SCRIPTPATH}/ATACDoubletDetector.py --rfilter ${REPFILTER} ${OUTDIR}/Overlaps.txt ${OUTDIR}/OverlapSummary.txt ${OUTDIR}

