import numpy as np
import pandas as pd

#Author: Asa Thibodeau

def convertToPositionFormatFromBED(data, startidx=1):
    """Converts BED format positions to position format.
        
        Parameters
        ----------
        data : array-like, shape (n_peaks, n_features)
        A numpy array containing peak data in BED format.
        
        startidx : int
        The start index of the array. Default: 1
        
        Returns
        -------
        posdata : array-like, shape (n_peaks, n_features)
        Returns a copy of the array incrementing the start
        position by one.
        """
    posdata = np.copy(data)
    posdata[:,startidx] = posdata[:,startidx]+1
    return posdata

def convertToBEDFormatFromPosition(data, startidx=1):
    """Converts position format positions to BED format.
        
        Parameters
        ----------
        data : array-like, shape (n_peaks, n_features)
        A numpy array containing peak data in position format.
        
        startidx : int
        The start index of the array. Default: 1
        
        Returns
        -------
        beddata : array-like, shape (n_peaks, n_features)
        Returns a copy of the array decrementing the start
        position by one.
        """
    beddata = np.copy(data)
    beddata[:,startidx] = beddata[:,startidx]-1
    return beddata


###################################################################
#The functions below assume position format. (1-base fully-closed)#
#This is only an issue if peaks from the same dataset overlap.    #
#If they do not, then this can be ignored.                        #
###################################################################


def getChrStartSorted(data, chridx=0, startidx=1):
    """Returns a dictionary by chromosome where each element in
        the dictionary contains an array containing start positions
        sorted in ascending order and the position of the original
        data element containing the start position.
        
        Time complexity: O(nlogn), n = # of peaks
        
        Parameters
        ----------
        data : array-like, shape (n_peaks, n_features)
        A numpy array containing peak data in position format.
        
        chridx : int
        The chromsome index of the array. Default: 0
        
        startidx : int
        The start index of the array. Default: 1
        
        Returns
        -------
        rv : dict
        A dictionary mapping each chromosome to an array:
        [0] = start position [1] = index in original data.
        """
    allchr = np.unique(data[:,chridx])
    rv = dict()
    for i in range(0, len(allchr)):
        idx = np.where(data[:,chridx] == allchr[i])[0]
        chrdata = data[idx,:]
        sidx = np.argsort(chrdata[:,startidx],kind="mergesort")
        rv[allchr[i]] = np.concatenate((np.transpose(chrdata[sidx,startidx][np.newaxis]), np.transpose(idx[sidx][np.newaxis])), axis=1)
    return rv


def getOverlappingRegions(chrom, start, end, chrstartsorted, data, eidx=2):
    """Returns the index of all regions that overlap the given
        chromosome, start, and end position.
        
        Time Complexity: O(logn), n = # of elements in the region list
        
        Parameters
        ----------
        chrom : str
        Chromsome of the position.
        
        start : int
        The start position.
        
        end : int
        The end position
        
        chrstartsorted : dict
        The chromosome start sorted dictionary of the regions
        to identify whether the given coordinates overlap.
        
        data : array-like, shape (n_peaks, n_features)
        The data corresponding to the sorted dictionary in
        positon format.
        
        eidx : int
        The end position indec within the data parameter. Default: 2
        
        Returns
        -------
        rv : tuple
        A tuple containing the index positions of data that
        overlap the given position.
        """
    try:
        startsorted = chrstartsorted[chrom]
    except:
        startsorted = []
    s = 0
    e = len(startsorted)
    while (e-s) > 1:
        mi = int(s+((e-s)/2))
        mstart = startsorted[mi,0]
        if mstart < start:
            s = mi
        elif mstart > start:
            e = mi
        else:
            s = mi
            e = mi
    #scan until starts are greater than end
    rv = []
    idx = s
    while idx < len(startsorted) and end > startsorted[idx,0]:
        didx = startsorted[idx,1]
        cstart = startsorted[idx,0]
        cend = data[didx,eidx]
        if start <= cend and end >= cstart: #position format comparison
            rv.append(didx)
        idx = idx+1
    return tuple(rv)


def getOverlapIndex(data, peakset, chridx=0, startidx=1, endidx=2, setchridx=0, setstartidx=1, setendidx=2):
    """Returns a boolean vector indicating whether or not the peak
        in the list overlaps a set of peaks.
        
        Time Complexity: O(nlogn), n = # of elements in the region list
        
        Parameters
        ----------
        data : array-like, shape (n_peaks, n_features)
        A numpy array containing peak data in position format.
        
        
        peakset : array-like, shape (n_peaks, n_features)
        A numpy array containing peak data in position format.
        
        
        chridx : int
        The chromsome index of the data parameter. Default: 0
        
        startidx : int
        The start index of the data parameter. Default: 1
        
        endidx : int
        The end index of the data parameter. Default: 2
        
        setchridx : int
        The chromosome index of the peakset parameter. Default: 0
        
        setstartidx : int
        The start index of the peakset parameter. Default: 1
        
        setendidx : int
        The end index of the peakset parameter. Default: 2
        
        Returns
        -------
        rv : arraylike, shape (n_peaks,)
        A boolean vector indicating whether the peaks in data overlap
        with the peaks in the peakset.
        """
    sortedconsensus = getChrStartSorted(peakset, setchridx, setstartidx)
    rv = np.zeros(len(data), dtype=bool)
    for i in range(0, len(data)):
        curchr = data[i, chridx]
        curstart = data[i, startidx]
        curend = data[i, endidx]
        if(len(getOverlappingRegions(curchr, curstart, curend, sortedconsensus, peakset, setendidx)) > 0):
            rv[i] = True
    return rv

def getOverlapCount(countdataset, datasets, chridx=0, startidx=1, endidx=2):
    """Counts how many peaks & in which dataset the current peak
        list overlaps over a set of peak lists.
        
        Time Complexity: O(m*nlogn)
        n = # of elements in the region list
        m = # of datasets
        
        Parameters
        ----------
        countdataset : array-like, shape (n_peaks, n_features)
        A numpy array containing peak data in position format.
        
        datasets : tuple
        A tuple containing multiple peak datasets in position format.
        
        chridx : int
        The chromsome index of the data parameter. Default: 0
        
        startidx : int
        The start index of the data parameter. Default: 1
        
        endidx : int
        The end index of the data parameter. Default: 2
        
        Returns
        -------
        overlapvector : arraylike, shape (n_peaks,)
        A vector indicating the number of datasets in datasets
        overlapping with the corresponding peak in countdataset.
        
        overlapmatrix : arraylike, shape (n_peaks, n_datasets)
        A boolean matrix indicating whether the peak countdataset
        overlaps with a peak in the corresponding dataset. Columns
        are ordered respective of the ordering in the dataset tuple.
        """
    overlapvector = np.zeros(len(countdataset))
    overlapmatrix = np.zeros((len(countdataset), len(datasets)))
    for i in range(0, len(datasets)):
        curv = getOverlapIndex(countdataset, datasets[i], chridx, startidx, endidx, chridx, startidx).astype(int)
        overlapmatrix[:,i] = curv
        overlapvector = overlapvector+curv
    return overlapvector, overlapmatrix



def getStrictConsensusPeaks(data, chridx=0, startidx=1, endidx=2):
    """Returns a strict set of consensus peaks requiring that every
        sample overlaps every other sample in the consensus peak region.
        
        Time Complexity: O(m^2*n + m*nlogn) m=#of datasets, n=#peaks
        
        Parameters
        ----------
        data : tuple
        A tuple of numpy arrays containing peak data in position format.
        
        chridx : int
        The chromsome index of the data parameter. Default: 0
        
        startidx : int
        The start index of the data parameter. Default: 1
        
        endidx : int
        The end index of the data parameter. Default: 2
        
        Returns
        -------
        rv : arraylike, shape (n_peaks,)
        A numpy array of genomic positions where peaks across
        all datasets overlap.
        
        Notes
        -----
        Chromosome start and end positions must be the same over all peak
        datasets.
        """
    
    sorteddata = dict()
    for i in range(0, len(data)):
        sorteddata[i] = getChrStartSorted(data[i], chridx, startidx)
    
    chromosomes = list(sorteddata[0])
    for i in range(1, len(data)):
        chromosomes = np.union1d(chromosomes, list(sorteddata[i]))
    
    allchrstrictpeaks = []
    for curchr in chromosomes:
        counters = np.zeros(len(data), dtype=np.int32)
        strictpeaks = []
        
        while True:
            curlist = []
            for i in range(0, len(data)):
                curdata = data[i]
                curindex = sorteddata[i][curchr][counters[i],1]
                curlist.append(curdata[curindex,:])
            curlist = np.array(curlist)
            
            #TODO get min and max positions
            minendidx = 0
            minend = curlist[0, endidx]
            maxstartidx = 0
            maxstart = curlist[0, startidx]
            for i in range(1, len(curlist)):
                curstart = curlist[i, startidx]
                if curstart > maxstart:
                    maxstart = curstart
                    maxstartidx = i
                curend = curlist[i, endidx]
                if curend < minend:
                    minend = curend
                    minendidx = i
            
            if maxstart < minend:
                #Add strict consensus peak
                strictpeaks.append([curchr, maxstart, minend])
            
            #Remove the least end
            counters[minendidx] = counters[minendidx]+1
            
            
            maxreached = False
            for i in range(0, len(counters)):
                if(counters[i] >= len(sorteddata[i][curchr])):
                    maxreached = True
                    break;
            
            if maxreached:
                break
        
        for cp in strictpeaks:
            allchrstrictpeaks.append(cp)
    
    return np.array(allchrstrictpeaks, dtype=object)


def getUnionPeaks(datasets, chridx=0, startidx=1, endidx=2):
    combineddata = np.concatenate(datasets)
    sortedlocations = getChrStartSorted(combineddata)
    rv = []
    for curchr in sortedlocations:
        locations = sortedlocations[curchr]
        curloci = combineddata[locations[0,1],:3]

        for i in range(1,len(locations)):
            nextloci = combineddata[locations[i,1],:3]
            if nextloci[1] > curloci[2]:
                rv.append([curchr, curloci[1], curloci[2]])
                curloci = nextloci
            else:
                curloci[2] = max(curloci[2], nextloci[2])

        rv.append([curchr, curloci[1], curloci[2]])
    return np.array(rv, dtype=np.object)


