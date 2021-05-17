package org.jax.atacdoubletdetectorgui;

import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;

public class ATACDoubletDetectorGUI {

	private String _path = "./";
	private final JFileChooser _fc = new JFileChooser();
	private final FileFilter _txtfilter = new FileNameExtensionFilter("Text File (.txt)", "txt");
	private final FileFilter _bamfilter = new FileNameExtensionFilter("BAM File (.bam)", "bam");
	private final FileFilter _csvfilter = new FileNameExtensionFilter("CSV File (.csv)", "csv");
	private final FileFilter _bedfilter = new FileNameExtensionFilter("Chromosome Position File (.bed, .txt)", "bed", "txt");
	private final Font _commandfont = new Font("Courier", Font.PLAIN, 12);
	private JFrame _frame;
	
	public static void main(String[] args) {
		try {
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		} catch (UnsupportedLookAndFeelException e) {
		} catch (ClassNotFoundException e) {
		} catch (InstantiationException e) {
		} catch (IllegalAccessException e) {
		}
		new ATACDoubletDetectorGUI();

	}

	public ATACDoubletDetectorGUI(){
		_path = System.getProperty("user.dir");
		if(!_path.endsWith("/")){
			_path += "/";
		}
		
		_fc.setCurrentDirectory(new File(_path));
		
		//Set up the Frame
		_frame = new JFrame();
		_frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		/*JPanel featureextpanel = getFeatureExtractionPanel();
		JPanel featureannpanel = getFeatureAnnotationPanel();
		JPanel promoterpanelt = getPromoterTSSPanel();
		JPanel promoterpanelm = getPromoterModelPanel();
		JPanel enhpanel = getEnhancerPredictionPanel();
		JPanel trainpanel = getTrainModelPanel();
		JPanel configpanel = getConfigPanel();

		JPanel mainpanel = getMainPanel(featureextpanel,featureannpanel,promoterpanelt, promoterpanelm, enhpanel, trainpanel, configpanel);
		*/
		_frame.add(getMainPanel());
		_frame.setSize(560,500);
		_frame.setResizable(false);
		_frame.setVisible(true);
		_frame.setTitle("ATACDoubletDetector");
		
	}
	
	private JPanel getMainPanel() {
		UIUtil util = new UIUtil();

		JPanel rvpanel = new JPanel();
		rvpanel.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
		rvpanel.setLayout(new BoxLayout(rvpanel, BoxLayout.Y_AXIS));

		
		JPanel basicpanel = new JPanel();
		basicpanel.setLayout(new BoxLayout(basicpanel, BoxLayout.Y_AXIS));
		
		
		final JTextField bamfield = new JTextField();
		util.setFieldPanel(basicpanel, _fc, new JLabel("Position Sorted bam:"), bamfield, new JButton("..."), _bamfilter, false);
		
		final JTextField csvfield = new JTextField();
		util.setFieldPanel(basicpanel, _fc, new JLabel("Single cell csv:"), csvfield, new JButton("..."), _csvfilter, false);
		
		
		final JTextField outdirfield = new JTextField();
		util.setFieldPanel(basicpanel, _fc, new JLabel("Output Directory:"), outdirfield, new JButton("..."), null, true);
		
		
		final JTextField chromfield = new JTextField(_path+"human_autosomes.txt");
		util.setFieldPanel(basicpanel, _fc, new JLabel("Chromosome List:"), chromfield, new JButton("..."), _txtfilter, false);
		
		
		final JTextField repeatfield = new JTextField(_path+"blacklist_repeats_segdups_rmsk_hg38.bed");
		util.setFieldPanel(basicpanel, _fc, new JLabel("Repeat Regions:"), repeatfield, new JButton("..."), _bedfilter, false);
		
		
		
		
		JPanel advancedpanel = new JPanel();
		advancedpanel.setLayout(new BoxLayout(advancedpanel, BoxLayout.Y_AXIS));
		
		
		final JTextField bcfield = new JTextField("CB");
		advancedpanel.add(util.getHorizontalField(new JLabel("BAM Barcode Attribute:"), bcfield));
		
		final JTextField bcidx = new JTextField("0");
		advancedpanel.add(util.getHorizontalField(new JLabel("Barcode index (CSV):"), bcidx));
		
		final JTextField cellidfield = new JTextField("0");
		advancedpanel.add(util.getHorizontalField(new JLabel("Cell id index (CSV):"), cellidfield));
		
		final JTextField iscellfield = new JTextField("9");
		advancedpanel.add(util.getHorizontalField(new JLabel("Is cell index (CSV):"), iscellfield));
		
		MoreOptionsPanel advancedbuttonpanel = new MoreOptionsPanel(advancedpanel);

		
		CommandOutputPanel commandoutput = new CommandOutputPanel(_commandfont, "Commands/Errors:");
		rvpanel.add(basicpanel);
		rvpanel.add(advancedbuttonpanel);
		rvpanel.add(commandoutput);
		
		JPanel buttonpanel = new JPanel();
		buttonpanel.setLayout(new BoxLayout(buttonpanel, BoxLayout.X_AXIS));
		
		JButton getcommandbutton = new JButton("Get Shell Command");
		getcommandbutton.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				RunATACDoubletDetectorValidator validator = new RunATACDoubletDetectorValidator(
						_path, 
						bamfield.getText(), 
						csvfield.getText(), 
						outdirfield.getText(), 
						chromfield.getText(), 
						repeatfield.getText(), 
						bcfield.getText(), 
						bcidx.getText(), 
						cellidfield.getText(), 
						iscellfield.getText());
				
				if(validator.hasError()) {
					commandoutput.reset();
					commandoutput.addLine(validator.getErrorMessage());
				}
				else {
					commandoutput.reset();
					commandoutput.addLine(validator.getCommand());
				}
				
				
			}
			
		});
		
		JButton runcommandbutton = new JButton("Run Shell Script");
		CommandRunner cr = new CommandRunner(_frame, _commandfont, "Run ATAC-DoubletDetector");

		runcommandbutton.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				RunATACDoubletDetectorValidator validator = new RunATACDoubletDetectorValidator(
						_path, 
						bamfield.getText(), 
						csvfield.getText(), 
						outdirfield.getText(), 
						chromfield.getText(), 
						repeatfield.getText(), 
						bcfield.getText(), 
						bcidx.getText(), 
						cellidfield.getText(), 
						iscellfield.getText());
				
				if(validator.hasError()) {
					commandoutput.reset();
					commandoutput.addLine(validator.getErrorMessage());
				}
				else {
					commandoutput.reset();
					commandoutput.addLine(validator.getCommand());
					cr.run(validator.getCommandArray());
				}
				
				
			}
			
		});
		
		buttonpanel.add(getcommandbutton);
		buttonpanel.add(runcommandbutton);

		
		rvpanel.add(buttonpanel);
		
		
		return rvpanel;
	}

}
