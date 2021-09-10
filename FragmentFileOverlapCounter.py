import numpy as np
import pandas as pd
import gzip
import time
import argparse

parser = argparse.ArgumentParser(description='OverlapCounter')
parser.add_argument("fragmentfile", help="Tab delimited file of fragments. (*.txt, *.tsv, *.txt.gz, *.tsv.gz)")
parser.add_argument("singlecellfile", help="Csv File containing 'barcode' and 'is__cell_barcode'.")
parser.add_argument("chrfile", help="List of chromosomes to read fragments from.")
parser.add_argument("outdir", help="The output directory.")

parser.add_argument('--expectedoverlap', dest='expectedoverlap', type=int, default=2,
                    help='Expected number of reads overlapping. (Default: 2)')

parser.add_argument('--maxinsertsize', dest='maxinsertsize', type=int, default=900,
                    help='The maximum insert size (in bp) between read pairs. (Default: 900)')
                    
                    
parser.add_argument('--startbases', dest='startbases', type=int, default=0,
                    help='The amount of bases add to the start position. (Default: 0)')

parser.add_argument('--endbases', dest='endbases', type=int, default=0,
                    help='The amount of bases to add to the end position. (Default: 0)')
                    

#Input arguments
args = parser.parse_args()
fragmentfile = args.fragmentfile
singlecellfile = args.singlecellfile
chrfile = args.chrfile
outdir = args.outdir

#optional arguments
expectedoverlap = args.expectedoverlap
maxinsert = args.maxinsertsize
startbases = args.startbases
endbases = args.endbases


chromosomes = list(pd.read_csv(chrfile, header=None).values.flatten())
chrdict = dict()
for curchr in chromosomes:
    chrdict[curchr] = True

def getOverlaps(reads, expectedoverlap):
    #If there are no reads, we are done
    if len(reads) <= expectedoverlap:
        return []

    #Create Overlap Index (1 if starting, -1 if ending) O(n)
    overlapindex = []
    for curread in reads:
        overlapindex.append([curread[1], 1])
        overlapindex.append([curread[2], -1])
    overlapindex = np.array(overlapindex)
    
    #Sort the overlap index by position O(nlogn)
    overlapindex = overlapindex[np.argsort(overlapindex[:,0]),:] 
    
    indexsize = len(overlapindex)
    
    #Calculate running sum O(n)
    runningsum = [overlapindex[0,1]]
    runningsumpos = [overlapindex[0,0]]
    for i in range(1, indexsize):
        previ = i-1
        cursum = runningsum[-1]+overlapindex[i,1]
        if overlapindex[previ,0] == overlapindex[i,0]:
            #Sum same positions together
            runningsum[-1] = cursum
        else:
            #Start a new position
            runningsum.append(cursum)
            runningsumpos.append(overlapindex[i,0])
    
    #Detect overlaps > the expected and report regions using the running sum O(n)
    rv = []
    chromosome = reads[0][0]
    withinsegment = False
    segmentstart = -1
    minoverlap = -1
    maxoverlap = -1
    for i in range(0, len(runningsum)):
        if withinsegment:
            if runningsum[i] <= expectedoverlap:
                rv.append([chromosome, segmentstart, runningsumpos[i], minoverlap, maxoverlap])
                withinsegment = False
                segmentstart = -1
                minoverlap = -1
                maxoverlap = -1
            else:
                maxoverlap = max(maxoverlap, runningsum[i])
        else:
            if runningsum[i] > expectedoverlap:
                segmentstart = runningsumpos[i]
                minoverlap = runningsum[i]
                maxoverlap = runningsum[i]
                withinsegment = True
    
    if withinsegment:
        rv.append([chromosome, segmentstart, runningsumpos[-1], minoverlap, maxoverlap])
    
    assignReadsWithinOverlaps(rv, reads)
    
    return rv

def assignReadsWithinOverlaps(overlaps, reads):
    numreads = len(reads)
    ri = 0
    multioverlaps = []
    
    for curoverlap in overlaps:
        curolstart = curoverlap[1]
        curolend = curoverlap[2]
        curstarts = []
        curends = []

        nextmultioverlaps = []
        for curmultioverlap in multioverlaps:
            curreadstart = curmultioverlap[1]
            curreadend = curmultioverlap[2]
            
            if curreadend >= curolstart and curolend >= curreadstart:
                curstarts.append(curreadstart)
                curends.append(curreadend)

            if curreadend >= curolend:
                nextmultioverlaps.append(curmultioverlap)
        
        multioverlaps = nextmultioverlaps

        while ri < numreads:
            curread = reads[ri]
            curreadstart = curread[1]
            curreadend = curread[2]
            
            if curreadend >= curolstart and curolend >= curreadstart:
                curstarts.append(curreadstart)
                curends.append(curreadend)
            
            if curreadend >= curolend:
                multioverlaps.append(curread)
            
            ri += 1
                
        startstring = ""
        endstring = ""
        for i in range(0, len(curstarts)):
            startstring += str(curstarts[i])+","
            endstring += str(curends[i])+","
            
        curoverlap.append(startstring)
        curoverlap.append(endstring)

def writeOverlaps(writer, overlaps, barcode):
    overlapout = []
    for curoverlap in overlaps:
        overlapout.append(curoverlap[0])
        overlapout.append("\t")
        overlapout.append(str(curoverlap[1]))
        overlapout.append("\t")
        overlapout.append(str(curoverlap[2]))
        overlapout.append("\t")
        overlapout.append(barcode)
        overlapout.append("\t")
        overlapout.append(str(curoverlap[3]))
        overlapout.append("\t")
        overlapout.append(str(curoverlap[4]))
        overlapout.append("\t")
        overlapout.append(".")
        overlapout.append("\t")
        overlapout.append(".")
        overlapout.append("\t")
        overlapout.append(".")
        overlapout.append("\t")
        overlapout.append(curoverlap[5])
        overlapout.append("\t")
        overlapout.append(curoverlap[6])
        overlapout.append("\n")
    overlapout = ''.join(overlapout)
    writer.write(overlapout)

def writeOverlapSummary(filepath, overlapcounts, vreadspercell, readspercell):
    
    colnames = ["Cell Id",
               "Number of Valid Reads",
               "Number of Overlaps",
               "Barcode",
               "Total Number of Reads"]
    
    
    table = np.empty(shape=(len(overlapcounts),len(colnames)), dtype=np.object)

    rowindex = 0
    for curbarcode in overlapcounts.keys():
        table[rowindex, 0] = curbarcode
        table[rowindex, 1] = vreadspercell[curbarcode]
        table[rowindex, 2] = overlapcounts[curbarcode]
        table[rowindex, 3] = curbarcode
        table[rowindex, 4] = readspercell[curbarcode]
        rowindex += 1
        
    pd.DataFrame(table, columns=colnames).to_csv(filepath, index=None, sep="\t")
    
def stringDecoder(line):
    return line

def gzipDecoder(line):
    return line.decode("utf-8")


def findOverlaps(fragmentfile, singlecellfile, expectedoverlap, chrlist, path, maxinsertsize=900):
    starttime = time.process_time()
    
    fragmentreader = []
    linedecoder = None

    if fragmentfile.endswith(".tsv.gz") or fragmentfile.endswith(".txt.gz"):
        fragmentreader = gzip.open(fragmentfile)
        linedecoder = gzipDecoder
    elif fragmentfile.endswith(".tsv") or fragmentfile.endswith(".txt"):
        fragmentreader = open(fragmentfile, 'r')
        linedecoder = stringDecoder
    else:
        print("Fragment file must be *.txt, *.txt.gz, *.tsv, or *.tsv.gz")
        return
    
    #Set up barcode maps
    sc_data = pd.read_csv(singlecellfile)
    sc_data = sc_data[sc_data['is__cell_barcode'] == 1]
    bc_map = dict()
    previous_reads = dict()
    previous_ends = dict()
    overlapcounts = dict()
    vreadspercell = dict()
    readspercell = dict()

    for curbarcode in sc_data['barcode']:
        bc_map[curbarcode] = []
        previous_reads[curbarcode] = []
        previous_ends[curbarcode] = -1
        overlapcounts[curbarcode] = 0
        vreadspercell[curbarcode] = 0
        readspercell[curbarcode] = 0

        
    #set up the overlap writer
    overlapwriter = open(path+"/Overlaps.txt", 'w')
    overlapwriter.write("chr\tstart\tend\tcell id\tMin Overlap Count\tMax Overlap Count\tMean Mapping Quality\tMin Mapping Quality\tMax Mapping Quality\tStarts\tEnds\n")
    
    #Variable to specify the current chromosome being processed
    chromosomeset = ""
    
    stats_readcount = 0
    stats_cellcount = len(bc_map)
    stats_overlapsizes = []
    
    #Loop through all reads to detect overlaps 
    for curline in fragmentreader:
        stats_readcount += 1
        
        decodedline = linedecoder(curline)

        if decodedline.strip().startswith("#"):
            continue

        split = decodedline.split("\t")

        curchr = split[0]
        curstart = int(split[1])+startbases
        curend = int(split[2])+endbases
        curbarcode = split[3]

        curlocation = [curchr, curstart, curend]

        insertsize = curend-curstart
        
        if curbarcode not in bc_map:
            continue

        readspercell[curbarcode] = readspercell[curbarcode]+1 
        
        if insertsize > maxinsertsize:
            continue #skip when the insert size is greater than the limit

        if curchr not in chrlist:
            continue #skip reads from chromosomes not in the provided chromosome list

        vreadspercell[curbarcode] = vreadspercell[curbarcode]+1 

        if curchr != chromosomeset:
            #Finish checking remaining overlaps and reset the running lists
            #This change indicates that we are done with the chromosome set as 'chromosomeset'

            #loop through all barcodes
            for curprevbarcode in previous_reads.keys():
                #1 - find overlaps with remaining lists
                curoverlaps = []
                if(len(previous_reads[curprevbarcode]) > 0):
                    curoverlaps = getOverlaps(previous_reads[curprevbarcode], expectedoverlap)
                    stats_overlapsizes.append(len(previous_reads[curprevbarcode]))
                #2 - writer overlaps
                writeOverlaps(overlapwriter, curoverlaps, curprevbarcode)
                #3 - Update overlap count for each
                overlapcounts[curprevbarcode] = overlapcounts[curprevbarcode]+len(curoverlaps)

                #set empty list
                previous_reads[curprevbarcode] = []
                #set prev end to -1
                previous_ends[curprevbarcode] = -1


            chromosomeset = curchr

        prevend = previous_ends[curbarcode]

        if prevend < curstart:
            #There can be no more overlaps for this segment.
            #Therefore, check for overlaps meeting the criteria and start a new list.

            #1: Find the overlaps
            curoverlaps = getOverlaps(previous_reads[curbarcode], expectedoverlap)
            stats_overlapsizes.append(len(previous_reads[curbarcode]))

            #2: Write the overlaps
            writeOverlaps(overlapwriter, curoverlaps, curbarcode)
            #3: Update the overlap count for the current barcode
            overlapcounts[curbarcode] = overlapcounts[curbarcode]+len(curoverlaps)
            #Start a new running list of reads to detect overlaps for the current barcode
            previous_reads[curbarcode] = [curlocation]


        else:
            #Add this read to the running list of overlaps for this barcode
            previous_reads[curbarcode].append(curlocation)

        #Assign a new endpoint for overlaps
        newend = max(prevend, curend);
        previous_ends[curbarcode] = newend
    fragmentreader.close()
    overlapwriter.close()
    writeOverlapSummary(path+"/OverlapSummary.txt", overlapcounts, vreadspercell, readspercell)
        
    endtime = time.process_time()
    elapsed = (endtime-starttime)
    print("Completed in "+str(elapsed)+" seconds.")
    print("Total number of reads: "+str(stats_readcount))
    print("Total number of cells: "+str(stats_cellcount))
    print("Total number of overlapping segments: "+str(len(stats_overlapsizes)))
    print("Average overlaps per segment: "+str(np.mean(stats_overlapsizes)))
    print("Max overlaps in segment: "+str(np.max(stats_overlapsizes)))
    
findOverlaps(fragmentfile, singlecellfile, expectedoverlap, chromosomes, outdir, maxinsert)
