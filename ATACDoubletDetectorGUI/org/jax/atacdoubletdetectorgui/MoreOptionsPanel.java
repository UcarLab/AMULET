package org.jax.atacdoubletdetectorgui;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JPanel;

@SuppressWarnings("serial")
public class MoreOptionsPanel extends JPanel{

	private String MOREOPTIONS = "Show Advanced Options";
	private String LESSOPTIONS = "Hide Advanced Options";
	private Component _component;
	
	public MoreOptionsPanel(Component c){
		_component = c;
		this.setLayout(new BorderLayout());
		this.add(c, BorderLayout.CENTER);
		c.setVisible(false);
		this.revalidate();
		this.repaint();
		
		JPanel buttonpanel = new JPanel();
		buttonpanel.setLayout(new BorderLayout());
		final JButton morelessbutton = new JButton(MOREOPTIONS);
		
		morelessbutton.addActionListener(new ActionListener(){

			@Override
			public void actionPerformed(ActionEvent e) {
				if(_component.isVisible()){
					morelessbutton.setText(MOREOPTIONS);
					_component.setVisible(false);
				}
				else{
					morelessbutton.setText(LESSOPTIONS);
					_component.setVisible(true);
				}
			}
			
		});
		buttonpanel.add(morelessbutton, BorderLayout.WEST);
		this.add(buttonpanel, BorderLayout.NORTH);
	}
	
	public Component getComponent(){
		return _component;
	}
	
}