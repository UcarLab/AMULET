package org.jax.snatacoverlapcounter;

import java.util.Iterator;
import java.util.LinkedList;

public class OverlapLocation extends Location {

	private int _minoverlap;
	private int _maxoverlap;
	
	private int _minmap;
	private int _maxmap;
	private double _meanmap;
	private LinkedList<Location> _overlaplocations;
	
	public OverlapLocation(String chr, int start, int end, int minoverlap, int maxoverlap) {
		super(chr, start, end, 0);
		_minoverlap = minoverlap;
		_maxoverlap = maxoverlap;
	}
	
	public void setOverlapInfo(LinkedList<Location> l) {
		_overlaplocations = l;
		_minmap = Integer.MAX_VALUE;
		_maxmap = Integer.MIN_VALUE;
		double runningmapping = 0;
		for(Iterator<Location> it = l.iterator(); it.hasNext();) {
			Location next = it.next();
			int mq = next.getMappingQuality();
			_minmap = Math.min(_minmap, mq);
			_maxmap = Math.max(_maxmap, mq);
			runningmapping += mq;
		}
		_meanmap = runningmapping/l.size();
	}
	
	public int getMinOverlap() {
		return _minoverlap;
	}
	
	public int getMaxOverlap() {
		return _maxoverlap;
	}
	
	public int getMinMappingQuality() {
		return _minmap;
	}
	
	public int getMaxMappingQuality() {
		return _maxmap;
	}
	
	public double getMeanMappingQuality() {
		return _meanmap;
	}
	
	public LinkedList<Location> getOverlapLocation() {
		return _overlaplocations;
	}
	
	public String getStarts() {
		StringBuilder sb = new StringBuilder();
		for(Iterator<Location> it = _overlaplocations.iterator(); it.hasNext();) {
			Location next = it.next();
			sb.append(next.getStart()+",");
		}
		return sb.toString();
	}

	public String getEnds() {
		StringBuilder sb = new StringBuilder();
		for(Iterator<Location> it = _overlaplocations.iterator(); it.hasNext();) {
			Location next = it.next();
			sb.append(next.getEnd()+",");
		}
		return sb.toString();
	}
}
