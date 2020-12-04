package org.jax.atacdoubletdetectorgui;

import java.awt.BorderLayout;
import java.awt.Font;

import javax.swing.BoxLayout;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

@SuppressWarnings("serial")
public class CommandOutputPanel extends JPanel {
	
	private JTextArea _ota;

	public CommandOutputPanel(Font font, String label){
		setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
		
		JLabel sol = new JLabel(label);
		sol.setHorizontalAlignment(JLabel.LEFT);
		_ota = new JTextArea(15, 35);
		_ota.setEditable(false);
		_ota.setFont(font);

		JScrollPane osp = new JScrollPane(_ota);
		osp.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		osp.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
		JPanel opanel = new JPanel(new BorderLayout());
		opanel.add(sol, BorderLayout.NORTH);
		opanel.add(osp, BorderLayout.CENTER);

		add(opanel);
	}
	
	public void addLine(String s){
		_ota.append(s+"\n");
		revalidate();
		repaint();
	}
	
	public void reset() {
		_ota.setText("");
		revalidate();
		repaint();
	}
	
}
