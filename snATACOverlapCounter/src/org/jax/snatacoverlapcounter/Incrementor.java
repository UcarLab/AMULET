package org.jax.snatacoverlapcounter;

public class Incrementor {

	private int _value;
	
	public Incrementor(int init) {
		_value = init;
	}
	
	public void increment() {
		_value++;
	}
	
	public void increment(int amt) {
		_value += amt;
	}
	
	
	public int getValue() {
		return _value;
	}
}
