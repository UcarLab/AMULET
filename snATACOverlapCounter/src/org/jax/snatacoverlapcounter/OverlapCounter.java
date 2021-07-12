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
import org.jax.snatacoverlapcounter.util.ValidReadSummary;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class OverlapCounter {
	
	private String _barcodeattribute;
	private int _barcodeidx, _cellididx, _iscellidx, _mapqthreshold, _maxinsertsize, _forwardcorrection, _reversecorrection, _expectedoverlap;
	private boolean _forcesorted;
	
	public static void main(String[] args) {
		
		String[] parsedargs = new String[5];
		String expectedoverlap = "2";
		String barcodeattribute = "CB";
		String barcodeidx = "0";
		String cellidx = "0";
		String iscellidx = "9";
		String maxinsertsize = "900";
		String mapqthreshold = "30";
		String forwardcorrection = "4";
		String reversecorrection = "-5";
		boolean forcesorted = false;
		
		int argidx = 0;
		for(int i = 0; i < args.length; i++) {
			
			if(i < args.length-1) {
				switch(args[i]) {
					case "--expectedoverlap":  
						expectedoverlap = args[i+1];
						i++;
						break;
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
					case "--startbases":  
						forwardcorrection = args[i+1];
						i++;
						break;
					case "--endbases":  
						reversecorrection = args[i+1];
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
			System.out.println("Options: --expectedoverlap    Expected number of reads overlapping. (Default: 2)");
			System.out.println("         --bambc     Bamfile attribute used for the barcode. (Default: \"CB\")");
			System.out.println("         --forcesorted Forces the input bam file to be treated as sorted.");
			System.out.println("         --bcidx     The column index of the CSV for barcode. (Default: 0)");
			System.out.println("         --cellidx   The column index of the CSV for cellid. (Default: 0)");
			System.out.println("         --iscellidx The index for determining cells (selecting values=1). (Default: 9)");
			System.out.println("         --mapqthresh Threshold for filtering low map quality reads (<= comparison). (Default: 30)");
			System.out.println("         --maxinsertsize The maximum insert size (in bp) between read pairs. (Default: 900)");
			System.out.println("         --startbases The amount of bases add to the start position. (Default: 4)");
			System.out.println("         --endbases The amount of bases to add to the end position (can be negative). (Default: -5)");

			System.exit(0);
		}
		long timestart = System.currentTimeMillis();
		
		String bamfile = parsedargs[0];
		String cellbarcodes = parsedargs[1];
		String chromlist = parsedargs[2];
		String outdir = parsedargs[3];
		
		OverlapCounter pc = new OverlapCounter(expectedoverlap, barcodeattribute, barcodeidx, cellidx, iscellidx, forcesorted, mapqthreshold, maxinsertsize, forwardcorrection, reversecorrection);
		
		try {
			pc.writeParameters(outdir+"/OverlapCounter-LastRunParameters.txt");
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
	
	public OverlapCounter(String expectedoverlap, String bca, String bcidx, String cellidx, String iscellidx, boolean forcesorted, String mapqthresh, String maxinsertsize, String forwardcorrection, String reversecorrection) {
		_barcodeattribute = bca;
		_forcesorted = forcesorted;
		try {
			_expectedoverlap = Integer.parseInt(expectedoverlap);
			_barcodeidx = Integer.parseInt(bcidx);
			_cellididx = Integer.parseInt(cellidx);
			_iscellidx = Integer.parseInt(iscellidx);
			_mapqthreshold = Integer.parseInt(mapqthresh);
			_maxinsertsize = Integer.parseInt(maxinsertsize);
			_forwardcorrection = Integer.parseInt(forwardcorrection);
			_reversecorrection = Integer.parseInt(reversecorrection);
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
		int overlapthresh = _expectedoverlap+1; //If the expected is 2 then we are looking for overlaps of 3 or more
		
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
		ValidReadSummary vrs = new ValidReadSummary();
		TreeSet<String> chromsizesmap = u.readChromSizes(chromsizes);
		TreeMap<String, String> cellbarcodemap = u.readCellBarcodes(cellbarcodes, _barcodeidx, _cellididx, _iscellidx);

		
		System.out.println("Reading BAM file.");
		int totalreads = 0;
		
		int lowmapq = 0;
		int insertsizeflag = 0;
		int positive = 0;
		int negative = 0;
		int notinchromosome = 0;
		int validreads = 0;
		long readlengthsum = 0;
		long insertsizesum = 0;
		
		
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
		
		TreeSet<String> omittedchromosomes = new TreeSet<String>();

		BufferedWriter bw = getOverlapWriter(outdir+"/Overlaps.txt");
		

		while(it.hasNext()){
			if(totalreads % 10000000 == 0) {
				System.out.println(Integer.toString(totalreads));
			}
			
			SAMRecord next = it.next();
			
			totalreads++;
			int curmapq = next.getMappingQuality();
			if(curmapq <= _mapqthreshold) {
				lowmapq++;
			}
			
			boolean isnegative = next.getReadNegativeStrandFlag();
			if(isnegative) {
				negative++;
			}
			else {
				positive++;
			}
			
			int insertsize = next.getInferredInsertSize();
			
			
			if(insertsize > 0 && curmapq > _mapqthreshold && vrs.isValidRead(next)) {
				
				String chr = next.getReferenceName();
				if(chromsizesmap.contains(chr)) {
					
					if(!chr.equals(curchr)) {
						//find/write overlaps
						for(Iterator<String> cellidit = cellbarcodemap.values().iterator(); cellidit.hasNext();) {
							String curcellid = cellidit.next();
							LinkedList<OverlapLocation> curol = findOverlaps(previousreads.get(curcellid), overlapthresh, curchr);
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

					int orig_start = next.getAlignmentStart()-1;
					int start = orig_start+_forwardcorrection;
					int end = orig_start+insertsize+_reversecorrection;
					int endtoendinsertsize = end-start;
					
					
					if(endtoendinsertsize > _maxinsertsize || endtoendinsertsize <= 0) {
						insertsizeflag++;
						continue;
					}
					
					
					readlengthsum += next.getReadLength();
					insertsizesum += endtoendinsertsize;
					validreads += 1;
					
					int mappingquality = next.getMappingQuality();
					
					vreadcounts.get(cellid).increment();
					treadcounts.get(cellid).increment();
					
					Location readlocation = new Location(curchr, start, end, mappingquality);
					LinkedList<Location> curlist = previousreads.get(cellid);
					int endposition = previousend.get(cellid);
					
					if(endposition < start) {
						LinkedList<OverlapLocation> ol = findOverlaps(curlist, overlapthresh, curchr);
						writeOverlaps(bw, cellid, ol);
						overlapcounts.get(cellid).increment(ol.size());

						//Start a new list
						LinkedList<Location> newlist = new LinkedList<Location>();
						newlist.add(readlocation);
						previousreads.put(cellid, newlist);						
					}
					else {
						//Current read overlaps with the previous reads and is added to the list
						curlist.add(readlocation);
					}
					int newend = Math.max(endposition, end);
					previousend.put(cellid, newend);
					
					
				}
				else {
					String barcode = u.getBarcode(next, _barcodeattribute);

					if(barcode == null) {
						continue;
					}
					String cellid = cellbarcodemap.get(barcode);
					
					if(cellid != null) {
						treadcounts.get(cellid).increment();
						omittedchromosomes.add(chr);
						notinchromosome++;
					}
					
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

		bw.flush();
		bw.close();
		
		writeOverlapCounts(vreadcounts, overlapcounts, cellbarcodemap, treadcounts, outdir+"/OverlapSummary.txt");
		writeReadInsertStatistics(totalreads, validreads, positive, negative, notinchromosome, vrs, insertsizeflag, lowmapq, readlengthsum, insertsizesum, outdir+"/StatSummary.txt");

		if(omittedchromosomes.size() > 0) {
			System.out.println("The following chromosomes were omitted:");
			for(Iterator<String> ocit = omittedchromosomes.iterator(); ocit.hasNext();) {
				System.out.println(ocit.next());
			}
		}
		
		System.out.println("Complete.");
	}
	
	
	private BufferedWriter getOverlapWriter(String filepath) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(filepath));
		
		StringBuilder sb = new StringBuilder();
		sb.append("chr");
		sb.append("\t");
		sb.append("start");
		sb.append("\t");
		sb.append("end");
		sb.append("\t");
		sb.append("cell id");
		sb.append("\t");
		sb.append("Min Overlap Count");
		sb.append("\t");
		sb.append("Max Overlap Count");
		sb.append("\t");
		sb.append("Mean Mapping Quality");
		sb.append("\t");
		sb.append("Min Mapping Quality");
		sb.append("\t");
		sb.append("Max Mapping Quality");
		sb.append("\t");
		sb.append("Starts");
		sb.append("\t");
		sb.append("Ends");
		sb.append("\n");
		
		bw.write(sb.toString());
		return bw;
	}
	

	private void writeParameters(String outfile) throws IOException {
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		StringBuilder sb = new StringBuilder();
		sb.append("Parameter:\tValue\n");
		
		sb.append("expectedoverlap:\t");
		sb.append(Integer.toString(_expectedoverlap));
		sb.append("\n");
		
		sb.append("bambc:\t");
		sb.append(_barcodeattribute);
		sb.append("\n");
		
		sb.append("forcesorted:\t");
		sb.append(Boolean.toString(_forcesorted));
		sb.append("\n");
		
		sb.append("bcidx:\t");
		sb.append(Integer.toString(_barcodeidx));
		sb.append("\n");
		
		sb.append("cellidx:\t");
		sb.append(Integer.toString(_cellididx));
		sb.append("\n");
		
		sb.append("iscellidx:\t");
		sb.append(Integer.toString(_iscellidx));
		sb.append("\n");
		
		sb.append("mapqthresh:\t");
		sb.append(Integer.toString(_mapqthreshold));
		sb.append("\n");
		
		sb.append("maxinsertsize:\t");
		sb.append(Integer.toString(_maxinsertsize));
		sb.append("\n");
		
		sb.append("startbases:\t");
		sb.append(Integer.toString(_forwardcorrection));
		sb.append("\n");
		
		sb.append("endbases:\t");
		sb.append(Integer.toString(_reversecorrection));
		sb.append("\n");
		
		bw.write(sb.toString());
		bw.flush();
		bw.close();
		
	}
	
	private void writeReadInsertStatistics(int totalreads, int validreads, int positive, int negative, int notinchromosome, ValidReadSummary vrs, int insert, int mq, long readlengthsum, long insertsizesum, String outfile) throws IOException {
		
		double meanreadlength = (double)readlengthsum/(double)validreads;
		double meaninsertsize = (double)insertsizesum/(double)validreads;
		
		StringBuilder sb = new StringBuilder();
		
		sb.append("Total Reads:\t");
		sb.append(Integer.toString(totalreads));
		
		sb.append("\nValid Read Pairs:\t");
		sb.append(Integer.toString(validreads));

		sb.append("\nPositive:\t");
		sb.append(Integer.toString(positive));
		
		sb.append("\nNegative:\t");
		sb.append(Integer.toString(negative));
		
		sb.append("\nNot in chromosome list:\t");
		sb.append(Integer.toString(notinchromosome));

		sb.append("\nInsert size (<=0 or >");
		sb.append(Integer.toString(_maxinsertsize));
		sb.append("):\t");
		sb.append(Integer.toString(insert));

		sb.append("\nLow Mapping Quality Reads (<=");
		sb.append(Integer.toString(_mapqthreshold));
		sb.append("):\t");
		sb.append(Integer.toString(mq));

		sb.append("\nMean Read Length:\t");
		sb.append(Double.toString(meanreadlength));
		
		sb.append("\nMean Insert Size:\t");
		sb.append(Double.toString(meaninsertsize));

		//Write SAMTools Flag Checks
		sb.append("\nSAMTools Flag Checked:\t");
		sb.append(Integer.toString(vrs.getTotalReads()));
		
		sb.append("\nSAMTools Flag Passed:\t");
		sb.append(Integer.toString(vrs.getValidReads()));
		
		sb.append("\nSAMTools Flag Duplicate:\t");
		sb.append(Integer.toString(vrs.getDuplicate()));
		
		sb.append("\nSAMTools Flag Unpaired:\t");
		sb.append(Integer.toString(vrs.getUnpaired()));
		
		sb.append("\nSAMTools Flag Unmapped:\t");
		sb.append(Integer.toString(vrs.getUnmapped()));
		
		sb.append("\nSAMTools Flag Mate Unmapped:\t");
		sb.append(Integer.toString(vrs.getMateUnmapped()));
		
		sb.append("\nSAMTools Flag Secondary:\t");
		sb.append(Integer.toString(vrs.getSecondary()));
		
		sb.append("\nSAMTools Flag Reference Mismatch:\t");
		sb.append(Integer.toString(vrs.getRefMismatch()));
		sb.append("\n");
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		bw.write(sb.toString());
		bw.flush();
		bw.close();
		
	}
	
	private void writeOverlapCounts(TreeMap<String, Incrementor> numreads, TreeMap<String, Incrementor> overlaps, TreeMap<String, String> barcode, TreeMap<String, Incrementor> tnumreads, String outfile) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		
		StringBuilder sb = new StringBuilder();
		sb.append("Cell Id");
		sb.append("\t");
		sb.append("Number of Valid Read Pairs");
		sb.append("\t");
		sb.append("Number of Overlaps");
		sb.append("\t");
		sb.append("Barcode");
		sb.append("\t");
		sb.append("Total Number of Reads");
		sb.append("\n");
		
		bw.write(sb.toString());
		
		TreeMap<String, String> revbarcode = new TreeMap<String, String>();
		for (Iterator<Entry<String, String>> it = barcode.entrySet().iterator(); it.hasNext();) {
			Entry<String, String> next = it.next();
			String cellid = next.getValue();
			String curbarcode = next.getKey();
			revbarcode.put(cellid, curbarcode);
		}
		
		
		for(Iterator<String> it = numreads.keySet().iterator(); it.hasNext();) {
			String cellid = it.next();
			
			sb = new StringBuilder();
			sb.append(cellid);
			sb.append("\t");
			sb.append(Integer.toString(numreads.get(cellid).getValue()));
			sb.append("\t");
			sb.append(Integer.toString(overlaps.get(cellid).getValue()));
			sb.append("\t");
			sb.append(revbarcode.get(cellid));
			sb.append("\t");
			sb.append(Integer.toString(tnumreads.get(cellid).getValue()));
			sb.append("\n");
			bw.write(sb.toString());
		}
		
		bw.flush();
		bw.close();
		
	}
	
	private void writeOverlaps(BufferedWriter bw, String cellid, LinkedList<OverlapLocation> ol) throws IOException {

		StringBuilder sb = new StringBuilder();
		for(Iterator<OverlapLocation> olit = ol.iterator(); olit.hasNext();) {
			OverlapLocation curoverlap = olit.next();
			
			sb.append(curoverlap.getChr());
			sb.append("\t");
			sb.append(Integer.toString(curoverlap.getStart()));
			sb.append("\t");
			sb.append(Integer.toString(curoverlap.getEnd()));
			sb.append("\t");
			sb.append(cellid);
			sb.append("\t");
			sb.append(Integer.toString(curoverlap.getMinOverlap()));
			sb.append("\t");
			sb.append(Integer.toString(curoverlap.getMaxOverlap()));
			sb.append("\t");
			sb.append(Double.toString(curoverlap.getMeanMappingQuality()));
			sb.append("\t");
			sb.append(Integer.toString(curoverlap.getMinMappingQuality()));
			sb.append("\t");
			sb.append(Integer.toString(curoverlap.getMaxMappingQuality()));
			sb.append("\t");
			sb.append(curoverlap.getStarts());
			sb.append("\t");
			sb.append(curoverlap.getEnds());
			sb.append("\n");
		}
		
		bw.write(sb.toString());
	}

	private LinkedList<OverlapLocation> findOverlaps(LinkedList<Location> l, int overlapthresh, String chr) {
		if(l == null) {
			return new LinkedList<OverlapLocation>();
		}
		LinkedList<OverlapIndex> ol = new LinkedList<OverlapIndex>();
		
		for(Iterator<Location> it = l.iterator(); it.hasNext();) {
			Location next = it.next();
			ol.add(new OverlapIndex(next.getStart(),1));
			ol.add(new OverlapIndex(next.getEnd(),-1));
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
				runningsumlist.add(new OverlapIndex(position, runningsum));
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
		
		for(Iterator<OverlapIndex> it = runningsumlist.iterator(); it.hasNext();) {
			OverlapIndex next = it.next();
			int curposition = next.getPosition();
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
		
		//Include overlaps that extend for the rest of the list via same end points
		if(start) {
			intervals.add(new OverlapLocation(chr, startposition, runningsumlist.getLast().getPosition(), minoverlap, maxoverlap));
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
