#!/bin/bash

EXPECTEDOVERLAP=2
BAMBC="CB"
BCIDX=0
CELLIDX=0
ISCELLIDX=9
MAPQTHRESH=30
MAXINSERT=900
FORCESORTED=""
STARTBASES="DEFAULT"
ENDBASES="DEFAULT"

ARGUMENTS=()

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --expectedoverlap) EXPECTEDOVERLAP="$2"; shift 2 ;;
        --bambc) BAMBC="$2"; shift 2 ;;
        --forcesorted) FORCESORTED="--forcesorted "; shift ;;
        --bcidx) BCIDX=$2; shift 2 ;;
        --cellidx) CELLIDX=$2; shift 2 ;;
        --iscellidx) ISCELLIDX=$2; shift 2 ;;
        --maxinsertsize) MAXINSERT=$2; shift 2 ;;
        --mapqthresh) MAPQTHRESH=$2; shift 2 ;;
        --startbases) STARTBASES="$2"; shift 2 ;;
        --endbases) ENDBASES="$2"; shift 2 ;;
        --*) echo "Unknown argument: $1"; exit 1 ;;
        *)  ARGUMENTS+=("$1"); shift ;;
    esac
done

if [[ ${#ARGUMENTS[@]} -ne 6 ]]; then
	echo "Usage: AMULET.sh bamfile barcodemap chromosomelist repeatfilter outputdirectory scriptpath";
	echo "Options: --expectedoverlap     Expected number of reads overlapping. (Default: 2)";
	echo "         --maxinsertsize The maximum insert size (in bp) between read pairs. (Default: 900)";
    echo "         --startbases The amount of bases add to the start position. (BAM Default: 4, Frag Default: 0)";
    echo "         --endbases The amount of bases to add to the end position. (BAM Default: -5, Frag Default: 0)";
	echo "         BAM Input Only Parameters:";
	echo "         --bambc     Bamfile attribute used for the barcode. (Default: \"CB\")";
	echo "         --forcesorted Forces the input bam file to be treated as sorted.";
	echo "         --bcidx     The column index of the CSV for barcode. (Default: 0)";
	echo "         --cellidx   The column index of the CSV for cellid. (Default: 0)";
	echo "         --iscellidx The index for determining cells (selecting values=1). (Default: 9)";
	echo "         --mapqthresh Threshold for filtering low map quality reads (<= comparison). (Default: 30)";
	echo "         Fragment file singlecell.csv input requires columns containing 'barcode' and 'is__cell_barcode':";

	exit 0;
fi

BAMFILE=${ARGUMENTS[0]}
BCMAP=${ARGUMENTS[1]}
CHRLIST=${ARGUMENTS[2]}
REPFILTER=${ARGUMENTS[3]}
OUTDIR=${ARGUMENTS[4]}
SCRIPTPATH=${ARGUMENTS[5]}

#Branch based on file extension. If BAM, use java code, if txt.gz/tsv.gz/txt/tsv, use python script
if [ "${BAMFILE: -4}" == ".bam" ]; then
	#These are defaults when the bam file is not shifted to account for the 9bp duplication
	#If the BAM file accounts for these, set these parameters to 0
    if [ ${STARTBASES} == "DEFAULT" ]; then
        STARTBASES=4
    fi
    if [ ${ENDBASES} == "DEFAULT" ]; then
        ENDBASES=-5
    fi

	#Another issue is duplicate marking. Better duplicate marking needs to account for the cell barcodes
	#Otherwise more reads are marked as duplicates when they are actually unique to that cell
    java -jar ${SCRIPTPATH}/snATACOverlapCounter.jar ${FORCESORTED}--expectedoverlap ${EXPECTEDOVERLAP} --bambc ${BAMBC} --bcidx ${BCIDX} --cellidx ${CELLIDX} --iscellidx ${ISCELLIDX} --mapqthresh ${MAPQTHRESH} --maxinsertsize ${MAXINSERT} --startbases ${STARTBASES} --endbases ${ENDBASES} ${BAMFILE} ${BCMAP} ${CHRLIST} ${OUTDIR}
    python3 ${SCRIPTPATH}/AMULET.py --rfilter ${REPFILTER} ${OUTDIR}/Overlaps.txt ${OUTDIR}/OverlapSummary.txt ${OUTDIR}

elif [ "${BAMFILE: -4}" == ".tsv" ] || [ "${BAMFILE: -4}" == ".txt" ] ||[ "${BAMFILE: -7}" == ".tsv.gz" ] || [ "${BAMFILE: -7}" == ".txt.gz" ]; then
    if [ ${STARTBASES} == "DEFAULT" ]; then
        STARTBASES=0
    fi
    if [ ${ENDBASES} == "DEFAULT" ]; then
        ENDBASES=0
    fi

    python3 ${SCRIPTPATH}/FragmentFileOverlapCounter.py --maxinsertsize ${MAXINSERT} --expectedoverlap ${EXPECTEDOVERLAP} --startbases ${STARTBASES} --endbases ${ENDBASES} ${BAMFILE} ${BCMAP} ${CHRLIST} ${OUTDIR}
    python3 ${SCRIPTPATH}/AMULET.py --rfilter ${REPFILTER} ${OUTDIR}/Overlaps.txt ${OUTDIR}/OverlapSummary.txt ${OUTDIR}
else
	echo "Unsupported file formatted.";
fi


