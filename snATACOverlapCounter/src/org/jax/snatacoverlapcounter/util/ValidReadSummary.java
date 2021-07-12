package org.jax.snatacoverlapcounter.util;

import htsjdk.samtools.SAMRecord;

public class ValidReadSummary {

	
	private int _totalreads, _unpaired, _unmapped, _mateunmapped, _secondary, _duplicate, _refmismatch, _validreads;
	
	public ValidReadSummary() {
		_totalreads = 0;
		_unpaired = 0;
		_unmapped = 0;
		_mateunmapped = 0;
		_secondary = 0;
		_duplicate = 0;
		_refmismatch = 0;
		_validreads = 0;
	}
	
	public int getTotalReads() {
		return _totalreads;
	}
	
	public int getUnpaired() {
		return _unpaired;
	}
	
	public int getUnmapped() {
		return _unmapped;
	}
	
	public int getMateUnmapped() {
		return _mateunmapped;
	}
	
	public int getSecondary() {
		return _secondary;
	}
	
	public int getDuplicate() {
		return _duplicate;
	}
	
	public int getRefMismatch() {
		return _refmismatch;
	}
	
	public int getValidReads() {
		return _validreads;
	}
	
	public boolean isValidRead(SAMRecord record){
		//Assume the read is valid until it isn't
		_totalreads++;
		boolean isvalid = true;
		
		if(!record.getReadPairedFlag()) {
			_unpaired++;
			isvalid=false;
		}
		
		if(record.getReadUnmappedFlag()) {
			_unmapped++;
			isvalid=false;
		}
		
		if(record.getMateUnmappedFlag()) {
			_mateunmapped++;
			isvalid=false;
		}
		
		if(record.isSecondaryOrSupplementary()) {
			_secondary++;
			isvalid=false;
		}
		
		if(record.getDuplicateReadFlag()) {
			_duplicate++;
			isvalid=false;
		}
		
		
		if(record.getReferenceIndex() != record.getMateReferenceIndex()){
			_refmismatch++;
			isvalid=false;
		}
		
		if(isvalid) {
			_validreads++;
		}
		
		return isvalid;
	}
	

	
}
