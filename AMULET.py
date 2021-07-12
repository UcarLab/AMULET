import numpy as np
import pandas as pd
import peakoverlap as po
import scipy.stats as stats
import statsmodels.api as sm
import argparse


parser = argparse.ArgumentParser(description='AMULET: ATAC-seq MULtiplet Estimation Tool.')
parser.add_argument("overlaps")
parser.add_argument("overlapsummary")
parser.add_argument("outputdirectory")

parser.add_argument('--rfilter', dest='rfilter', help='Filepath of regions (e.g., known repetitive regions) to exclude.')

parser.add_argument('--q', dest='qthreshold', type=float, default=0.01,
                    help='FDR corrected probability threshold. (Default: 0.01)')

parser.add_argument('--qrep', dest='qrepthreshold', type=float, default=0.01,
                    help='FDR corrected probability threshold for inferring repetitive regions. (Default: 0.01)')

parser.add_argument('--expectedoverlap', dest='expectedoverlap', type=int, default=2,
                    help='Expected number of reads overlapping. (Default: 2)')

parser.add_argument('--minoverlap', dest='minoverlap', type=int, default=1,
                    help='The minimum length (in bp) of overlap to keep. (Default: 1)')

args = parser.parse_args()
qvalthresh = args.qthreshold
qvalrepthresh = args.qrepthreshold
outdir = args.outputdirectory
expectedoverlap = args.expectedoverlap
minoverlap = args.minoverlap-1

def generateMatrix(data, cellids, unionoverlaps):
    #Map cellids to integers
    celliddict = dict()
    rcelliddict = dict()
    for i in range(len(cellids)):
        celliddict[cellids[i]] = i
        rcelliddict[i] = cellids[i]
        
    sortedunionoverlaps = po.getChrStartSorted(unionoverlaps)

    regioninfo = dict()
    
    matrix = np.zeros((len(unionoverlaps), len(cellids)))
    for i in range(len(data)):
        curchr = data[i,0]
        curstart = data[i,1]
        curend = data[i,2]
        cellid = data[i,3]
        overlap = po.getOverlappingRegions(curchr, curstart, curend+1, sortedunionoverlaps, unionoverlaps)
        
        if cellid in celliddict:
            for oi in overlap:
                matrix[oi, celliddict[cellid]] = 1

                if oi not in regioninfo:
                    regioninfo[oi] = []
                mergedlength = unionoverlaps[oi][2]-unionoverlaps[oi][1]+1
                length = curend-curstart+1
                regioninfo[oi].append([length, length/mergedlength])
    
    return matrix, rcelliddict, regioninfo

def inferRepeats(matrix, unionoverlaps, threshold):
    rowsum = np.sum(matrix,axis=1)
    rep_probabilities = []
    rep_mean = np.mean(rowsum[:])
    
    for curval in rowsum:
        rep_probabilities.append(stats.poisson.sf(curval,rep_mean))
    
    rep_probabilities = np.array(rep_probabilities)
    corrected_rep_probabilities = sm.stats.multipletests(rep_probabilities, method='fdr_bh')
    
    rep_regions = unionoverlaps[corrected_rep_probabilities[1] < threshold]
    non_repregions = unionoverlaps[corrected_rep_probabilities[1] >= threshold]

    rep_probas = np.concatenate((unionoverlaps, rowsum[:,np.newaxis], rep_probabilities[:,np.newaxis], corrected_rep_probabilities[1][:,np.newaxis]), axis=1)
    return rep_regions, non_repregions, rep_probas

def getDoublets(matrix, unionoverlaps, rcelliddict):
    colsum = np.sum(matrix,axis=0)
    doublet_probabilities = []
    doublet_mean = np.mean(colsum[:])
    
    for curval in colsum:
        doublet_probabilities.append(stats.poisson.sf(curval,doublet_mean))
    
    doublet_probabilities = np.array(doublet_probabilities)
    corrected_doublet_probabilities = sm.stats.multipletests(doublet_probabilities, method='fdr_bh')

    doublets_probas = []
    
    index = 0
    for cur in corrected_doublet_probabilities[1]:
        doublets_probas.append([rcelliddict[index], doublet_probabilities[index], cur])
        index += 1

    return np.array(doublets_probas, dtype=np.object)



def getFilteredOverlaps(data, simplerepeats, expectedoverlap):
    rv = []
    sorted_repeats = po.getChrStartSorted(simplerepeats)
    for curoverlap in data:
        curchr = curoverlap[0]
        starts = np.array(str(curoverlap[-2]).split(",")[:-1],dtype=int)
        ends = np.array(str(curoverlap[-1]).split(",")[:-1],dtype=int)

        observedloci = dict()

        newstarts = []
        newends = []
        for i in range(len(starts)):
            key = str(starts[i])+"-"+str(ends[i])
            if key not in observedloci:
                observedloci[key] = True
                
                overlap = po.getOverlappingRegions(curchr, starts[i], ends[i], sorted_repeats, simplerepeats)
                if len(overlap) == 0:
                    newstarts.append(starts[i])
                    newends.append(ends[i])
                
        if len(newstarts) < len(starts):
            if len(newstarts) > expectedoverlap:
                #Recalculate overlaps

                #starts increment by 1
                counts = np.ones((len(starts), 2), dtype=np.object)
                counts[:,0] = starts
                counts

                #ends decrement by 1
                counts2 = -1*np.ones((len(starts), 2), dtype=np.object)
                counts2[:,0] = ends
                counts2

                #Combine the counts and sort them
                combinedcounts = np.concatenate((counts,counts2))
                countorder = np.argsort(combinedcounts[:,0])
                combinedcounts[countorder]

                #Scan through and maintain a running sum
                #when the running sum is > 2, continue until <= 2 and report that overlap
                runningsum = 0
                i = 0
                startoverlap = False
                startoverlapposition = 0

                while i < len(combinedcounts):
                    runningsum += combinedcounts[i][1]
                    j = i+1
                    while j < len(combinedcounts):
                        if combinedcounts[i,0] == combinedcounts[j,0]:
                            runningsum += combinedcounts[j][1]
                            j += 1
                        else:
                            break

                    if not startoverlap and runningsum > expectedoverlap:
                        startoverlap = True
                        startoverlapposition = combinedcounts[i][0]
                    elif startoverlap and runningsum <= expectedoverlap:
                        #append overlap
                        rv.append([curchr, startoverlapposition, combinedcounts[i][0], curoverlap[3]])
                        startoverlap= False

                    i = j
                if startoverlap:
                    rv.append([curchr, startoverlapposition, combinedcounts[-1][0], curoverlap[3]])
             
        else:
            rv.append(curoverlap[:4])
    rv = np.array(rv)
    return rv

######################
#Load Data/Preprocess#
######################


#Step 1: Load data & repeat region filter
data = pd.read_csv(args.overlaps, sep="\t").values
summarydata = pd.read_csv(args.overlapsummary, sep="\t").values
cellids = summarydata[:,0]

simplerepeats = np.zeros((0,3))
if args.rfilter:
    print("Filtering regions.")
    simplerepeats = po.getUnionPeaks([pd.read_csv(args.rfilter, sep="\t", header=None).values[:,0:3]])


#Step 2: Filter repetitive elementss
filtereddata = getFilteredOverlaps(data, simplerepeats, expectedoverlap)

#filter overlaps that are < a specified bp
lengths = filtereddata[:,2]-filtereddata[:,1]+1
filtereddata = filtereddata[lengths > minoverlap,:]
numfiltered = len(data)-len(filtereddata)

print("Number of regions filtered: "+str(numfiltered)+" ("+str(100*numfiltered/len(data))+"%)")

###################
#Doublet Detection#
###################

print("Detecting multiplets.")

#Step 3: Generate Matrix
unionoverlaps = po.getUnionPeaks([filtereddata])

matrix,rcelliddict,_ = generateMatrix(filtereddata, cellids, unionoverlaps)

#Step 4: Infer repetitive regions
repetitive,nonrepetitive, _ = inferRepeats(matrix, unionoverlaps, qvalrepthresh)

#Step 5: filter repetitive regions & run doublet detection
repfilterindex, _ = po.getOverlapCount(filtereddata, tuple([repetitive]))
repfiltereddata = filtereddata[repfilterindex == 0,:]
repfiltered_unionoverlaps = po.getUnionPeaks([repfiltereddata])
repfiltered_matrix, rcelliddict2, _ = generateMatrix(repfiltereddata, cellids, repfiltered_unionoverlaps)
doublets_with_prob = getDoublets(repfiltered_matrix, repfiltered_unionoverlaps, rcelliddict2)

summarydata_dict = dict()
for i in range(len(summarydata)):
    summarydata_dict[summarydata[i,0]] = summarydata[i,:]
    

doublets_cellids = []
doublets_barcodes = []
doublets_with_prob_barcode = []

for i in range(len(doublets_with_prob)):
    summaryrow = summarydata_dict[doublets_with_prob[i,0]]
    
    doublets_with_prob_barcode.append([doublets_with_prob[i,0], summaryrow[3], doublets_with_prob[i,1], doublets_with_prob[i,2]])
    
    if doublets_with_prob[i,2] < qvalthresh:
        doublets_cellids.append(doublets_with_prob[i,0])
        doublets_barcodes.append(summaryrow[3])

#Output doublets
pd.DataFrame(doublets_cellids).to_csv(outdir+"/MultipletCellIds_"+str(qvalthresh).split(".")[1]+".txt", header=None, index=None, sep="\t")
pd.DataFrame(doublets_barcodes).to_csv(outdir+"/MultipletBarcodes_"+str(qvalthresh).split(".")[1]+".txt", header=None, index=None, sep="\t")

pd.DataFrame(doublets_with_prob_barcode, columns=["cell_id", "barcode", "p-value", "q-value"]).to_csv(outdir+"/MultipletProbabilities.txt", index=None, sep="\t")

#Output stats
stats_numbercells = len(cellids)
stats_numberunionregions = len(unionoverlaps)
stats_numberdoublets = len(doublets_cellids)
stats_percentdoublets = stats_numberdoublets*100/stats_numbercells
doublet_stats = np.array([["Number of Cells", stats_numbercells], ["Number of Merged Regions", stats_numberunionregions], ["Number of Multiplets", stats_numberdoublets], ["Multiplet Percent", stats_percentdoublets]])

pd.DataFrame(doublet_stats).to_csv(outdir+"/MultipletSummary.txt", index=None, header=None, sep="\t")


print("Done.")
