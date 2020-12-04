package org.jax.atacdoubletdetectorgui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.GroupLayout;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.GroupLayout.ParallelGroup;
import javax.swing.GroupLayout.SequentialGroup;
import javax.swing.filechooser.FileFilter;

public class UIUtil {
	public JPanel getHorizontalField(JComponent ...c ){
		JPanel rvpanel = new JPanel();
		GroupLayout rvlayout = new GroupLayout(rvpanel);
		rvpanel.setLayout(rvlayout);
		SequentialGroup seqgroup = rvlayout.createSequentialGroup();
		ParallelGroup pgroup = rvlayout.createParallelGroup(GroupLayout.Alignment.BASELINE);
		for(int i = 0; i < c.length; i++){
			seqgroup.addComponent(c[i]);
			pgroup.addComponent(c[i]);
		}
		rvlayout.setHorizontalGroup(seqgroup);
		rvlayout.setVerticalGroup(rvlayout.createSequentialGroup().addGroup(pgroup));
		return rvpanel;
	}
	
	public void setFileAction(final JFileChooser fc, final JPanel comp, final JTextField filetext, final JButton filebutton, final FileFilter filter, final boolean directory){
		final ActionListener bamaction = new ActionListener(){
			@Override
			public void actionPerformed(ActionEvent e) {
				fc.resetChoosableFileFilters();
				fc.setFileFilter(filter);
				fc.setMultiSelectionEnabled(false);
				if(directory){
					fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				}
				else{
					fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
				}
				int result = fc.showOpenDialog(comp);
				if(result == JFileChooser.APPROVE_OPTION){
					filetext.setText(fc.getSelectedFile().getAbsolutePath());
				}
			}
		};
		filebutton.addActionListener(bamaction);
	}
	
	public void setFieldPanel(JPanel panel, JFileChooser fc, JLabel label, JTextField textfield, JButton button, FileFilter filter, boolean directory){
		textfield.setColumns(12);
		setFileAction(fc, panel, textfield, button, filter, directory);
		panel.add(getHorizontalField(label, textfield, button));
	}
	
	
}
