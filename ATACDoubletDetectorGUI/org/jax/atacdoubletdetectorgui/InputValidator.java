package org.jax.atacdoubletdetectorgui;

public interface InputValidator {

	public boolean hasError();
	
	public String getErrorMessage();
	
	public String getCommand();
	
	public String[] getCommandArray();
	
}