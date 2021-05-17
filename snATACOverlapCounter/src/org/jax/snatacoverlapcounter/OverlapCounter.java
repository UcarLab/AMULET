package org.jax.snatacoverlapcounter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.TreeSet;

import org.jax.snatacoverlapcounter.util.Util;

import htsjdk.samtools.SAMFileHeader;
//import htsjdk.samtools.SAMFileHeader;
//import htsjdk.samtools.SAMFileWriter;
//import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class OverlapCounter {
	
	private String _barcodeattribute;
	private int _barcodeidx, _cellididx, _iscellidx, _mapqthreshold, _maxinsertsize;
	private boolean _forcesorted;
	
	public static void main(String[] args) {
		
		String[] parsedargs = new String[5];
		String barcodeattribute = "CB";
		String barcodeidx = "0";
		String cellidx = "0";
		String iscellidx = "9";
		String maxinsertsize = "900";
		String mapqthreshold = "30";
		boolean forcesorted = false;
		
		int argidx = 0;
		for(int i = 0; i < args.length; i++) {
			
			if(i < args.length-1) {
				switch(args[i]) {
					case "--bambc":  
						barcodeattribute = args[i+1];
						i++;
						break;
					case "--cellidx":  
						cellidx = args[i+1];
						i++;
						break;
					case "--bcidx":  
						barcodeidx = args[i+1];
						i++;
						break;
					case "--iscellidx":  
						iscellidx = args[i+1];
						i++;
						break;
					case "--maxinsertsize":  
						maxinsertsize = args[i+1];
						i++;
						break;
					case "--mapqthresh":  
						mapqthreshold = args[i+1];
						i++;
						break;
					case "--forcesorted":  
						forcesorted = true;
						break;
					default:
						if(args[i].startsWith("-")) {
							argidx = 5;
						}
						else {
							parsedargs[argidx++] = args[i];
						}
						break;
				}
			}
			else {
				parsedargs[argidx++] = args[i];
			}
			
			if(argidx > 4) {
				break;
			}
		}
		
		if(argidx != 4) {
			System.out.println("Usage: bamfile cellidbarcodemap chromosomelist outputdirectory");
			System.out.println("Options: --bambc     Bamfile attribute used for the barcode. (Default=\"CB\")");
			System.out.println("         --forcesorted Forces the input bam file to be treated as sorted.");
			System.out.println("         --bcidx     The column index of the CSV for barcode. (Default: 0)");
			System.out.println("         --cellidx   The column index of the CSV for cellid. (Default: 0)");
			System.out.println("         --iscellidx The index for determining cells (selecting values=1). (Default: 9)");
			System.out.println("         --mapqthresh Threshold for filtering low map quality reads (<= comparison). (Default: 30)");
			System.out.println("         --maxinsertsize The maximum insert size (in bp) between read pairs. (Default: 900)");

			System.exit(0);
		}
		long timestart = System.currentTimeMillis();
		
		String bamfile = parsedargs[0];
		String cellbarcodes = parsedargs[1];
		String chromlist = parsedargs[2];
		String outdir = parsedargs[3];
		
		OverlapCounter pc = new OverlapCounter(barcodeattribute, barcodeidx, cellidx, iscellidx, forcesorted, mapqthreshold, maxinsertsize);
		try {
			 pc.findOverlaps(bamfile, cellbarcodes, chromlist, outdir);
		} catch (IOException e) {
			e.printStackTrace();
			try {
				Util u = new Util();
				u.writeErrorFile(e, outdir);
			} catch (IOException e1) {
				e1.printStackTrace();
			}
			e.printStackTrace();
		}
		long timeend = System.currentTimeMillis();
		try {
			pc.writeRunTime(timestart, timeend, outdir);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public OverlapCounter(String bca, String bcidx, String cellidx, String iscellidx, boolean forcesorted, String mapqthresh, String maxinsertsize) {
		_barcodeattribute = bca;
		_forcesorted = forcesorted;
		try {
			_barcodeidx = Integer.parseInt(bcidx);
			_cellididx = Integer.parseInt(cellidx);
			_iscellidx = Integer.parseInt(iscellidx);
			_mapqthreshold = Integer.parseInt(mapqthresh);
			_maxinsertsize = Integer.parseInt(maxinsertsize);
		}
		catch(NumberFormatException e) {
			System.out.println("Please use integer values for index locations.");
		}
	}
	
	public void writeRunTime(long start, long end, String outdir) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outdir+"/RunTime.txt"));
		double duration = (double)(end-start)/1000;
		
		bw.write("Overlap counter finished in: "+Double.toString(duration)+" seconds.");
		bw.flush();
		bw.close();
	}
	
	public void findOverlaps(String bamfile, String cellbarcodes, String chromsizes, String outdir) throws IOException{

		SamReaderFactory factory = SamReaderFactory.makeDefault()
	              .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
	              .validationStringency(ValidationStringency.SILENT);
		
		final SamReader reader = factory.open(new File(bamfile));
		
		if(_forcesorted) {
			reader.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.coordinate);
		}
		
		if(!(reader.getFileHeader().getSortOrder().getComparatorInstance() instanceof SAMRecordCoordinateComparator)) {
			System.out.println("The input BAM file must be coordinate sorted. If you believe the bam file is sorted, use the --forcesorted option.");
			reader.close();
			System.exit(0);
		}
		
		SAMRecordIterator it = reader.iterator();
		Util u = new Util();
		TreeSet<String> chromsizesmap = u.readChromSizes(chromsizes);
		TreeMap<String, String> cellbarcodemap = u.readCellBarcodes(cellbarcodes, _barcodeidx, _cellididx, _iscellidx);

		
		System.out.println("Reading BAM file.");
		int totalreads = 0;
		int dupreads = 0;
		int lowmapq = 0;
		String curchr = null;

		
		TreeMap<String, LinkedList<Location>> previousreads = new TreeMap<String, LinkedList<Location>>();
		TreeMap<String, Integer> previousend = new TreeMap<String, Integer>();

		TreeMap<String, Incrementor> treadcounts = new TreeMap<String, Incrementor>();
		TreeMap<String, Incrementor> vreadcounts = new TreeMap<String, Incrementor>();
		TreeMap<String, Incrementor> overlapcounts = new TreeMap<String, Incrementor>();
		for(Iterator<String> cellidit = cellbarcodemap.values().iterator(); cellidit.hasNext();) {
			String cellid = cellidit.next();
			vreadcounts.put(cellid, new Incrementor(0));
			treadcounts.put(cellid, new Incrementor(0));
			overlapcounts.put(cellid, new Incrementor(0));
		}

		BufferedWriter bw = new BufferedWriter(new FileWriter(outdir+"/Overlaps.txt"));
		bw.write("chr");
		bw.write("\t");
		bw.write("start");
		bw.write("\t");
		bw.write("end");
		bw.write("\t");
		bw.write("cell id");
		bw.write("\t");
		bw.write("Min Overlap Count");
		bw.write("\t");
		bw.write("Max Overlap Count");
		bw.write("\t");
		bw.write("Mean Mapping Quality");
		bw.write("\t");
		bw.write("Min Mapping Quality");
		bw.write("\t");
		bw.write("Max Mapping Quality");
		bw.write("\t");
		bw.write("Starts");
		bw.write("\t");
		bw.write("Ends");
		bw.write("\n");
		
		
		int validreads = 0;
		long readlengthsum = 0;
		long insertsizesum = 0;
		
		while(it.hasNext()){
			if(totalreads % 10000000 == 0) {
				System.out.println(Integer.toString(totalreads));
			}
			
			SAMRecord next = it.next();
			
			totalreads++;
			if(next.getDuplicateReadFlag()) {
				dupreads++;
			}
			if(next.getMappingQuality() <= _mapqthreshold) {
				lowmapq++;
			}
			
			boolean isnegative = next.getReadNegativeStrandFlag();
			int insertsize = next.getInferredInsertSize();
			
			

			
			if(!isnegative && u.isValidRead(next) && insertsize > 0 && next.getMappingQuality() > _mapqthreshold) {
				
				String chr = next.getReferenceName();
				if(chromsizesmap.contains(chr)) {
					
					if(!chr.equals(curchr)) {
						//find/write overlaps
						for(Iterator<String> cellidit = cellbarcodemap.values().iterator(); cellidit.hasNext();) {
							String curcellid = cellidit.next();
							LinkedList<OverlapLocation> curol = findOverlaps(previousreads.get(curcellid), 3, curchr);
							writeOverlaps(bw, curcellid, curol);
							overlapcounts.get(curcellid).increment(curol.size());
						}
						
						//initialize observed read list, these will determine where the overlaps are
						previousreads = new TreeMap<String, LinkedList<Location>>();
						previousend = new TreeMap<String, Integer>();
						for(Iterator<String> cellidit = cellbarcodemap.values().iterator(); cellidit.hasNext();) {
							String cellid = cellidit.next();
							previousreads.put(cellid, new LinkedList<Location>());
							previousend.put(cellid, Integer.MIN_VALUE);
						}
						curchr = chr;
					}
					
					String barcode = u.getBarcode(next, _barcodeattribute);

					if(barcode == null) {
						continue;
					}
					
					String cellid = cellbarcodemap.get(barcode);
					
					if(cellid == null) {
						continue;
					}
					
					int start = next.getAlignmentStart();
					int end = start+insertsize-1;
					int endtoendinsertsize = end-start;
					
					if(end-start > _maxinsertsize) {
						continue;
					}
					
					
					readlengthsum += next.getReadLength();
					insertsizesum += endtoendinsertsize;
					//w.addAlignment(next);
					validreads += 1;
					
					int mappingquality = next.getMappingQuality();
					
					vreadcounts.get(cellid).increment();
					treadcounts.get(cellid).increment();
					
					Location readlocation = new Location(curchr, start, end, mappingquality);
					LinkedList<Location> curlist = previousreads.get(cellid);
					int endposition = previousend.get(cellid);
					
					if(endposition < start) {
						LinkedList<OverlapLocation> ol = findOverlaps(curlist, 3, curchr);
						writeOverlaps(bw, cellid, ol);
						overlapcounts.get(cellid).increment(ol.size());

						//Start a new list
						LinkedList<Location> newlist = new LinkedList<Location>();
						previousreads.put(cellid, newlist);						
					}
					else {
						//Current read overlaps with the previous reads and is added to the list
						curlist.add(readlocation);
					}
					int newend = Math.max(endposition, end);
					previousend.put(cellid, newend);
					
					
				}
				
			}
			else {
				String barcode = u.getBarcode(next, _barcodeattribute);

				if(barcode == null) {
					continue;
				}
				
				String cellid = cellbarcodemap.get(barcode);
				
				if(cellid == null) {
					continue;
				}
				
				treadcounts.get(cellid).increment();

			}
			
		}

		//w.close();
		
		bw.flush();
		bw.close();
		
		writeOverlapCounts(vreadcounts, overlapcounts, cellbarcodemap, treadcounts, outdir+"/OverlapSummary.txt");
		writeReadInsertStatistics(totalreads, dupreads, lowmapq, validreads, readlengthsum, insertsizesum, outdir+"/StatSummary.txt");
		
		System.out.println("Complete.");
	}
	
	private void writeReadInsertStatistics(int totalreads, int duplicates, int lowmapq, int validreads, long readlengthsum, long insertsizesum, String outfile) throws IOException {
		
		double meanreadlength = (double)readlengthsum/(double)validreads;
		double meaninsertsize = (double)insertsizesum/(double)validreads;
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		bw.write("Total Reads:\t"+Integer.toString(totalreads)+"\n");
		bw.write("Duplicate Reads:\t"+Integer.toString(duplicates)+"\n");
		bw.write("Low Mapping Quality Reads (<="+Integer.toString(_mapqthreshold)+"):\t"+Integer.toString(lowmapq)+"\n");
		bw.write("Valid Reads:\t"+Integer.toString(validreads)+"\n");
		bw.write("Total Reads:\t"+Integer.toString(totalreads)+"\n");

		bw.write("Mean Read Length:\t"+Double.toString(meanreadlength)+"\n");
		bw.write("Mean Insert Size:\t"+Double.toString(meaninsertsize)+"\n");

		bw.flush();
		bw.close();
		
	}
	
	private void writeOverlapCounts(TreeMap<String, Incrementor> numreads, TreeMap<String, Incrementor> overlaps, TreeMap<String, String> barcode, TreeMap<String, Incrementor> tnumreads, String outfile) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		
		bw.write("Cell Id");
		bw.write("\t");
		bw.write("Number of Valid Read Pairs");
		bw.write("\t");
		bw.write("Number of Overlaps");
		bw.write("\t");
		bw.write("Barcode");
		bw.write("\t");
		bw.write("Total Number of Reads");
		bw.write("\n");
		
		TreeMap<String, String> revbarcode = new TreeMap<String, String>();
		for (Iterator<Entry<String, String>> it = barcode.entrySet().iterator(); it.hasNext();) {
			Entry<String, String> next = it.next();
			String cellid = next.getValue();
			String curbarcode = next.getKey();
			revbarcode.put(cellid, curbarcode);
		}
		
		
		for(Iterator<String> it = numreads.keySet().iterator(); it.hasNext();) {
			String cellid = it.next();
			bw.write(cellid);
			bw.write("\t");
			bw.write(Integer.toString(numreads.get(cellid).getValue()));
			bw.write("\t");
			bw.write(Integer.toString(overlaps.get(cellid).getValue()));
			bw.write("\t");
			bw.write(revbarcode.get(cellid));
			bw.write("\t");
			bw.write(Integer.toString(tnumreads.get(cellid).getValue()));
			bw.write("\n");
		}
		
		bw.flush();
		bw.close();
		
	}
	
	private void writeOverlaps(BufferedWriter bw, String cellid, LinkedList<OverlapLocation> ol) throws IOException {

		
		for(Iterator<OverlapLocation> olit = ol.iterator(); olit.hasNext();) {
			OverlapLocation curoverlap = olit.next();
			
			bw.write(curoverlap.getChr());
			bw.write("\t");
			bw.write(Integer.toString(curoverlap.getStart()));
			bw.write("\t");
			bw.write(Integer.toString(curoverlap.getEnd()));
			bw.write("\t");
			bw.write(cellid);
			bw.write("\t");
			bw.write(Integer.toString(curoverlap.getMinOverlap()));
			bw.write("\t");
			bw.write(Integer.toString(curoverlap.getMaxOverlap()));
			bw.write("\t");
			bw.write(Double.toString(curoverlap.getMeanMappingQuality()));
			bw.write("\t");
			bw.write(Integer.toString(curoverlap.getMinMappingQuality()));
			bw.write("\t");
			bw.write(Integer.toString(curoverlap.getMaxMappingQuality()));
			bw.write("\t");
			bw.write(curoverlap.getStarts());
			bw.write("\t");
			bw.write(curoverlap.getEnds());
			bw.write("\n");
		}
	}

	private LinkedList<OverlapLocation> findOverlaps(LinkedList<Location> l, int overlapthresh, String chr) {
		if(l == null) {
			return new LinkedList<OverlapLocation>();
		}
		LinkedList<OverlapIndex> ol = new LinkedList<OverlapIndex>();
		
		for(Iterator<Location> it = l.iterator(); it.hasNext();) {
			Location next = it.next();
			ol.add(new OverlapIndex(next.getStart(),1));
			ol.add(new OverlapIndex(next.getEnd()+1,-1));
		}
		
		ol.sort(new OverlapIndexComparator());
		
		int runningsum = 0;
		int prevposition = -1;
		LinkedList<OverlapIndex> runningsumlist = new LinkedList<OverlapIndex>();
		
		for(Iterator<OverlapIndex> it = ol.iterator(); it.hasNext();){
			OverlapIndex next = it.next();
			runningsum += next.getValue();
			int position = next.getPosition();
			if(position == prevposition) {
				runningsumlist.getLast().setValue(runningsum);
			}
			else {
				runningsumlist.add(new OverlapIndex(runningsum, runningsum));
			}
			prevposition = position;
		}
		
		
		//Scan the running sum list and find all subintervals where the overlap is >= the overlapthresh,
		//breaking when dropping less than overlapthresh
		LinkedList<OverlapLocation> intervals = new LinkedList<OverlapLocation>();

		boolean start = false;
		int startposition = 0;
		int minoverlap = 0;
		int maxoverlap = 0;
		
		Iterator<OverlapIndex> positioniteror = ol.iterator();
		
		for(Iterator<OverlapIndex> it = runningsumlist.iterator(); it.hasNext();) {
			OverlapIndex next = it.next();
			OverlapIndex nextpos = positioniteror.next();
			
			int curposition = nextpos.getPosition();
			int curvalue = next.getValue();
			if(curvalue < overlapthresh && start) {
				intervals.add(new OverlapLocation(chr, startposition, curposition, minoverlap, maxoverlap));
				start = false;
				startposition = 0;
				minoverlap = 0;
				maxoverlap = 0;
			}
			else if(curvalue >= overlapthresh) {
				if(!start) {
					start = true;
					startposition = curposition;
					minoverlap = curvalue;
				}
				minoverlap = Math.min(minoverlap, curvalue);
				maxoverlap = Math.max(curvalue, maxoverlap);
			}
			
		}
		
		setOverlapInfo(intervals, l);
		return intervals;
		
	}
	
	private void setOverlapInfo(LinkedList<OverlapLocation> ol, LinkedList<Location> l) {
		Iterator<OverlapLocation> oli = ol.iterator();
		Iterator<Location> li = l.iterator();
		
		LinkedList<Location> multioverlaps = new LinkedList<Location>();
		while(oli.hasNext()) {
			OverlapLocation curol = oli.next();
			int curolstart = curol.getStart();
			int curolend = curol.getEnd();
			LinkedList<Location> curoverlaps = new LinkedList<Location>();

			//Check all previously checked locations that could be overlapping the current iteration
			LinkedList<Location> nextmultioverlaps = new LinkedList<Location>();
			for(Iterator<Location> poit = multioverlaps.iterator(); poit.hasNext();) {
				Location nextpo = poit.next();
				
				if(nextpo.getEnd() >= curolstart && curolend >= nextpo.getStart()) {
					curoverlaps.add(nextpo);
				}
				
				if(nextpo.getEnd() >= curolend) {
					nextmultioverlaps.add(nextpo);
				}
				
			}
			multioverlaps = nextmultioverlaps;
			
			//Check new locations until the next start > end
			while(li.hasNext()) {
				Location nextloci = li.next();
				
				if(nextloci.getEnd() >= curolstart && curolend >= nextloci.getStart()) {
					curoverlaps.add(nextloci);
				}
				
				if(nextloci.getEnd() >= curolend) {
					multioverlaps.add(nextloci);
				}
			}
			
			curol.setOverlapInfo(curoverlaps);
		}
		
	}
	
	private class OverlapIndex {
		
		private int _position;
		private int _value;
		
		public OverlapIndex(int position, int value) {
			_position = position;
			_value = value;
		}
		
		public int getPosition() {
			return _position;
		}
		
		public int getValue() {
			return _value;
		}
		
		public void setValue(int v) {
			_value = v;
		}
		
	}
	
	private class OverlapIndexComparator implements Comparator<OverlapIndex> {

		@Override
		public int compare(OverlapIndex o1, OverlapIndex o2) {
			if(o1.getPosition() < o2.getPosition()) {
				return -1;
			}
			else if(o1.getPosition() > o2.getPosition()) {
				return 1;
			}
			return 0;
		}
		
	}
	

	
}
