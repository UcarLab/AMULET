package org.jax.snatacoverlapcounter.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.TreeMap;
import java.util.TreeSet;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;

public class Util {

	
	public String getBarcode(SAMRecord read, String attributeid) {
		List<SAMTagAndValue> attributes = read.getAttributes();
		for(Iterator<SAMTagAndValue> it = attributes.iterator(); it.hasNext();) {
			SAMTagAndValue next = it.next();
			if(next.tag.equals(attributeid)){
				return (String) next.value;
			}
		}
		return null;
	}
	
	public TreeSet<String> readChromSizes(String chromsizefile) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(chromsizefile));
		TreeSet<String> rv = new TreeSet<String>();
		while(br.ready()) {
			String line = br.readLine();
			String[] split = line.split("\t|,");
			String chr = split[0];
			rv.add(chr);
		}
		br.close();
		return rv;
	}
	
	public TreeMap<String, String> readCellBarcodes(String chromsizefile, int barcodeidx, int cellididx, int iscellbarcodeidx) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(chromsizefile));
		TreeMap<String, String> rv = new TreeMap<String, String>();
		while(br.ready()) {
			String line = br.readLine();
			String[] split = line.split(",");
			boolean iscellbarcode = split[iscellbarcodeidx].equals("1");
			if(iscellbarcode) {
				String barcode = split[barcodeidx];
				String cellid = split[cellididx];
				rv.put(barcode, cellid);
			}
		}
		br.close();
		return rv;
	}

	public boolean isValidRead(SAMRecord record){
		if ((!record.getReadPairedFlag() ||
                record.getReadUnmappedFlag() ||
                record.getMateUnmappedFlag() ||
                record.isSecondaryOrSupplementary() ||
                record.getDuplicateReadFlag())) {
			return false;
		}
		
		if(record.getReferenceIndex() != record.getMateReferenceIndex()){
			return false;
		}
		return true;
	}
	
	public String[] getFiles(String file) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(file));
		LinkedList<String> rv = new LinkedList<String>();
		while(br.ready()) {
			rv.add(br.readLine());
		}
		br.close();
		return rv.toArray(new String[0]);
	}
	
	public void writeCompletionFile(String outdir, String prefix, String output) throws IOException {
		String date = new SimpleDateFormat("yyyy-MM-dd", Locale.getDefault()). format(new Date());
		String time = new SimpleDateFormat("HHmmss", Locale.getDefault()). format(new Date());

		BufferedWriter bw = new BufferedWriter(new FileWriter(outdir+"/"+prefix+"_"+date+"T"+time+".txt"));
		bw.write(output);
		bw.flush();
		bw.close();
	}
	
	public void writeErrorFile(Exception e, String outdir) throws IOException {
		StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw);
		e.printStackTrace(pw);
		String stacktrace = sw.toString();
		writeCompletionFile(outdir, "Error", stacktrace);
	}
	
	public String[] inferSampleNames(String[] files) {
		String[] rv = new String[files.length];
		
		for (int i = 0; i < files.length ; i++){
			String curfilepath = files[i];
			File curfile = new File(curfilepath);

			if(curfile.isDirectory()){
				int index = curfilepath.lastIndexOf("/", curfilepath.length()-1);
				if(index > -1) {
					rv[i] = curfilepath.substring(index+1);
				}
				else {
					rv[i] = "sample"+i;
				}
			}
			else {
				int index1 = curfilepath.lastIndexOf("/");
				int index2 = curfilepath.lastIndexOf("/", index1-1);
				
				if(index2 > -1) {
					rv[i] = curfilepath.substring(index2+1, index1);
				}
				else {
					rv[i] = "sample"+i;
				}
			}
			

		}
		
		return rv;
	}
	
}
