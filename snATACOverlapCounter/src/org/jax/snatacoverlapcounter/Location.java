package org.jax.snatacoverlapcounter;

public class Location {
	
	private String _chr;
	private int _start;
	private int _end;
	private int _mappingquality;

	public Location(String chr, int start, int end, int mq) {
		_chr = chr;
		_start = start;
		_end = end;
		_mappingquality = mq;
	}

	public String getChr() {
		return _chr;
	}
	
	public int getStart() {
		return _start;
	}
	
	public int getEnd() {
		return _end;
	}
	
	public int getMappingQuality() {
		return _mappingquality;
	}
	
}
