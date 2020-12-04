package org.jax.atacdoubletdetectorgui;

import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.SwingWorker;

public class CommandRunner {

	private JFrame _frame;
	private Font _labelfont;
	private Process _p;
	private String _title;
	private Function _postfunction;
	
	public CommandRunner(JFrame frame, Font lf, String title){
		_frame = frame;
		_labelfont = lf;
		_title = title;
	}
	
	public void setPostCommandFunction(Function f){
		_postfunction = f;
	}
	
	
	public void run(String[] command) {
		final JOptionPane pane = new JOptionPane();
		
		final CommandOutputPanel op = new CommandOutputPanel(_labelfont, "Command Line Output");
		
		final JProgressBar jpb = new JProgressBar();
		jpb.setSize(500, 10);
		jpb.setIndeterminate(true);
		
		final JButton cancel = new JButton("Cancel");
		final JButton finish = new JButton("Finish");

		
		final JPanel cp = new JPanel(new BorderLayout());
		cp.add(jpb, BorderLayout.CENTER);
		cp.add(cancel, BorderLayout.EAST);
		
		final JPanel combinedpanel = new JPanel(new BorderLayout());
		combinedpanel.add(op, BorderLayout.NORTH);
		combinedpanel.add(cp, BorderLayout.SOUTH);

		
		pane.setMessage(combinedpanel);
		pane.setOptions(new Object[]{});
		
		final JDialog dialog = pane.createDialog(_frame, _title);
		
		finish.addActionListener(new ActionListener(){
			@Override
			public void actionPerformed(ActionEvent e) {
				dialog.dispose();
			}
			
		});
		
		cancel.addActionListener(new ActionListener(){
			@Override
			public void actionPerformed(ActionEvent e) {
				if(_p != null){
					_p.destroy();
				}
				dialog.dispose();
			}
		});
		
		final ProcessBuilder pb = new ProcessBuilder(command);
		pb.redirectErrorStream(true);
		_p = null;

		SwingWorker<Void, Void> sw = new SwingWorker<Void, Void>(){

			@Override
			protected Void doInBackground() throws Exception {

				try {
					_p = pb.start();
					BufferedReader stdin = new BufferedReader(new InputStreamReader(_p.getInputStream()));
					
					String s = null;
					while((s = stdin.readLine()) != null){
						op.addLine(s);
						op.revalidate();
						op.repaint();
					}
					
					if(_postfunction != null){
						try {
							op.addLine(_postfunction.run());
							op.revalidate();
							op.repaint();
						} catch (IOException e1) {
							e1.printStackTrace();
						}
					}

					
				} catch (IOException e1) {
					e1.printStackTrace();
				}
				return null;
			}
			
			@Override 
			protected void done(){
				cp.remove(jpb);
				cp.remove(cancel);
				cp.add(finish, BorderLayout.EAST);
				dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
				dialog.revalidate();
				dialog.repaint();
			}
		
		};
		
		sw.execute();

		
		dialog.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
	    dialog.setVisible(true);
	    
	}

}