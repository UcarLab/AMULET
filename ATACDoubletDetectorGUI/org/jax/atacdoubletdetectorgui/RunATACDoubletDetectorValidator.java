package org.jax.atacdoubletdetectorgui;



import java.io.File;

public class RunATACDoubletDetectorValidator implements InputValidator{

	private String _path;
	private String _bamfile;
	private String _csvfile;
	private String _outdir;
	private String _chrfile;
	private String _repeatfile;

	private String _bambarcode;
	private int _barcodeidx;
	private int _cellidx;
	private int _iscellidx;
	
	private boolean _haserror;
	private String _errormessages;
	
	public RunATACDoubletDetectorValidator(String path, String bamfile, String csvfile, String outdir, String chrfile, String repeatfile, String bambarcode, String bcidx, String cellidx, String iscell){
		StringBuilder errormessages = new StringBuilder();
		_haserror = false;
		
		
		_path = path.trim();
		if(!(new File(_path)).exists()){
			_haserror = true;
			errormessages.append("Incorrect file path set for ATAC-DoubletDetector.\n");
		}


		
		_bamfile = bamfile.trim();
		if(_bamfile.equals("") || !(new File(_bamfile)).exists()){
			_haserror = true;
			errormessages.append("BAM file does not exist.\n");
		}
		
		_csvfile = csvfile.trim();
		if(_csvfile.equals("") || !(new File(_csvfile)).exists()){
			_haserror = true;
			errormessages.append("CSV file does not exist.\n");
		}
		
		_outdir = outdir.trim();
		if(!(new File(_outdir)).isDirectory()){
			_haserror = true;
			errormessages.append("Output directory does not exist.\n");
		}
		if(!_outdir.endsWith("/")){
			_outdir = _outdir+"/";
		}
		
		_chrfile = chrfile.trim();
		if(_chrfile.equals("") || !(new File(_chrfile)).exists()){
			_haserror = true;
			errormessages.append("Chromosome list lfile does not exist.\n");
		}
		
		_repeatfile = repeatfile.trim();
		if(_repeatfile.equals("") || !(new File(_repeatfile)).exists()){
			_haserror = true;
			errormessages.append("Repeat region file does not exist.\n");
		}
		
		
		_bambarcode = bambarcode.trim();
		if(_bambarcode.equals("")) {
			_haserror = true;
			errormessages.append("Please specify the attribute code for barcodes in the bam file (BAM Barcode Attribute).\n");
		}
		
		try {
			_barcodeidx = Integer.parseInt(bcidx.trim());
		} catch(NumberFormatException e) {
			_haserror = true;
			errormessages.append("Barcode index must be set as an integer value.\n");
		}
		
		
		try {
			_cellidx = Integer.parseInt(cellidx.trim());
		} catch(NumberFormatException e) {
			_haserror = true;
			errormessages.append("Cell ID index must be set as an integer value.\n");
		}
		
		try {
			_iscellidx = Integer.parseInt(iscell.trim());
		} catch(NumberFormatException e) {
			_haserror = true;
			errormessages.append("Is cell index must be set as an integer value.\n");
		}
		
		
		_errormessages = errormessages.toString();
	}
	
	@Override
	public boolean hasError() {
		return _haserror;
	}

	@Override
	public String getErrorMessage() {
		return _errormessages;
	}

	public String getCommand(){
		StringBuilder sb = new StringBuilder();
		
		sb.append("\"");
		sb.append(_path);
		sb.append("ATACDoubletDetector.sh\"");
		sb.append(" ");
		
		sb.append("--bambc ");
		sb.append(_bambarcode);
		sb.append(" ");
		
		sb.append("--bcidx ");
		sb.append(_barcodeidx);
		sb.append(" ");
		
		sb.append("--cellidx ");
		sb.append(_cellidx);
		sb.append(" ");
		
		sb.append("--iscellidx ");
		sb.append(_iscellidx);
		sb.append(" ");
		
		sb.append("\"");
		sb.append(_bamfile);
		sb.append("\"");
		sb.append(" ");
		
		sb.append("\"");
		sb.append(_csvfile);
		sb.append("\"");
		sb.append(" ");
		
		sb.append("\"");
		sb.append(_chrfile);
		sb.append("\"");
		sb.append(" ");
		
		sb.append("\"");
		sb.append(_repeatfile);
		sb.append("\"");
		sb.append(" ");

		sb.append("\"");
		sb.append(_outdir);
		sb.append("\"");
		sb.append(" ");
		
		sb.append("\"");
		sb.append(_path);
		sb.append("\"");
		
		return sb.toString();
	}
	
	@Override
	public String[] getCommandArray() {
		return new String[] {"/bin/bash", "--login", "-c", getCommand()};
	}


}