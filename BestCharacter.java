package bestCharacter;

import java.awt.*; 
import java.awt.event.*;

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.text.BadLocationException;

import java.util.*;

/*                               Author: Nadia Talent, 2010-2014 
*                 The author waives copyright and places this work in the Public Domain.
*                        http://creativecommons.org/publicdomain/zero/1.0/
*
*                            Bestchar source code, Java version, 2014
*
* A program for calculating some representative 'Best character' coefficients for
* categorical data (the same procedure would be used for numeric ranges
* after they are converted to categories). A single character is processed.
*
*------------------------------------------------------------------------------------------------------
* The input text is a single column representing one character. 
* The lines list character states for each taxon, with '/' to separate states.
* Spaces and tab characters are ignored.
*
* For example, a character called 'Petal color' with the states for four taxa could be coded:
* red
* white/red
* white
* white
*
* Output appears in the graphical user interface, and some intermediate values appear on the console.
*------------------------------------------------------------------------------------------------------
* Coefficients based on information theory
* 1. Information coefficient reverse engineered from Intkey, logarithms to base 2 are used.
*    using {} to represent subscripts
*    H=sigma(H{i}), H{i}=-(k{i}/t)*log2(k{i}/t)*SIG{i}/k{i}
*    where t=number of remaining taxa,
*    j varies from 1 to t,
*    s is the total number of states allowed for this character (including all taxa),
*    i varies from 1 to s,
*    k{i} is the number of taxa that have state i, 
*    SIG{i}=sigma(1/n{j})
* 2. The 'Normalized information coefficient' as described by 
*    Talent, Dickinson, and Dickinson (2014, in Biodiversity Informatics); 
*    logarithms to base t are used.
*
* Coefficients that use pairwise comparisons of taxa, i.e. distance measures:
* 1. The Separation Coefficient counts only taxa that are completely separable, so if any state
*    is possible for two taxa they are not separable
* 2. Pairwise average Jaccard Coefficient */

public class BestCharacter extends JApplet {
	JButton jbtCompute = new JButton("Click to compute");
	TwoPanels canvas = new TwoPanels(); // TwoPanels class is declared below
	// data entered, treated as a queue, i.e., one line per taxon holds an array of allowed states
	ArrayList<String[]> totalTaxa = new ArrayList<String[]>();

	ArrayList<String> stateList = new ArrayList<String>(); // all the states that appear in the input
	ArrayList<Integer> listOfStateCounts = new ArrayList<Integer>(); // These are called k{i} above, in H{i}=log2(k{i}/t)*SIG{i}/t).
	ArrayList<Double> listOfSigmas = new ArrayList<Double>(); // These are called SIG{i} above.

	public BestCharacter() { // constructor for the applet class
		add(canvas, BorderLayout.CENTER);
		jbtCompute.setForeground(Color.BLUE); // text color of the button
		Border paneEdge = BorderFactory.createRaisedBevelBorder();
		jbtCompute.setBorder(paneEdge);
		JPanel panel = new JPanel(new FlowLayout()); // A panel to stop the button taking up the entire width
		panel.add(jbtCompute);
		add(panel, BorderLayout.SOUTH);		
		jbtCompute.addActionListener(new ComputeListener()); // Register the listener for the Compute button.		
	}

	/******************************** main method ***********************************/
	// This main method only enables the applet to run as an application, if required.
	// Otherwise, it is not used.
	public static void main(String[] args) { 
		JFrame frame = new JFrame("BestCharacter Applet"); // frame with title

		// Add an applet instance to the frame.
		frame.add(new BestCharacter(), BorderLayout.CENTER);

		// Display the frame
		frame.pack(); // better than setting the size
		frame.setLocationRelativeTo(null); // Centre the frame.
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setVisible(true);
	}

	/************************ inner class to detect button click **************************/
	class ComputeListener implements ActionListener { // listener for the Compute button
		@Override
		public void actionPerformed(ActionEvent e) {
			canvas.compute(); // Take appropriate action
		}
	}

	/************ inner class, a panel that contains input and output panels *****************/
	class TwoPanels extends JPanel { 
		JTextField characterNameField  = new JTextField ("");
		JPanel p1 = new JPanel(new BorderLayout()); // left outer panel, for the input fields and instructions
		JPanel p2 = new JPanel(new BorderLayout()); // right panel to contain the output text area
		JPanel p3 = new JPanel(new BorderLayout()); // panel to group the text entry fields and the input instructions
		
		JTextArea instructionTextArea = new JTextArea("Input is a single column representing one character, " + 
			"with character states\nfor one taxon per line, '/' to separate alternative states. Spaces and " +
			"tab\ncharacters are ignored. e.g.: " +
			"petal color for four taxa could be coded:\nred\nwhite/red\nwhite\nwhite\n");

		JTextArea inputTextArea = new JTextArea(25, 35); // Declare the input text area.
		JScrollPane scrollArea = new JScrollPane(inputTextArea); 
		JTextArea outputTextArea = new JTextArea(35, 35); // Declare the output text area.
			
		public TwoPanels() { // constructor
			Border greenLine = BorderFactory.createLoweredSoftBevelBorder();

			characterNameField.setBorder(new TitledBorder(greenLine, "Character name (optional, e.g., Petal Color)"));
			scrollArea.setBorder(new TitledBorder(greenLine, "Input, one taxon per line"));
			p3.add(characterNameField, BorderLayout.NORTH);
			p3.add(instructionTextArea, BorderLayout.SOUTH);
			instructionTextArea.setEditable(false);
			
			p1.add(p3, BorderLayout.NORTH);
			p1.add(scrollArea, BorderLayout.SOUTH);
			
			outputTextArea.setEditable(false);
			outputTextArea.setBorder(new TitledBorder(greenLine, "Results"));
			System.out.println("Bestchar program, Java version, diagnostic output\n");

			p2.add(new JScrollPane(outputTextArea)); // Add the output textarea within a scrollpane.
			add(p1, BorderLayout.WEST); // Add the left panel.
			add(p2, BorderLayout.EAST); // Add the right panel.
		}
		
	/*************** method for when the Compute button is clicked **************/
		public void compute() { // This is the top-level controller for most of the work.
			totalTaxa.clear();
			stateList.clear();
			listOfStateCounts.clear();
			listOfSigmas.clear();

			String charName = characterNameField.getText(); // optional field
			if (! charName.equals("")) {
				System.out.println("Character=" + charName + "\n");
			}
			
			inputText();
			int numTaxa = totalTaxa.size(); // called t in the comments above
			totalTaxaToString(totalTaxa); // show the data structure in the console log
			outputTextArea.append("Number of taxa= " + numTaxa + "\n");
			System.out.println("length=" + numTaxa);
			
			calculateStates();
			int totalNumStates = stateList.size();
			System.out.println("totalNumStates=" + totalNumStates);
			outputTextArea.append("\nNumber of states= " + totalNumStates);
			
			computeInformation(numTaxa, totalNumStates);
			
			computeOtherCoefficients(numTaxa, totalNumStates);

		}	

	/************************** Subsidiary methods ******************************/
		public void inputText() {
			String sReturn = inputTextArea.getText(); // This limits the amount of input to 65535 bytes of UTF-8 data.
			int numRows = inputTextArea.getLineCount(); // number of lines of data entered

			outputTextArea.setText(""); // Initialize.

			// Process the lines, one line per taxon.
			int j = 0; 
			int lineOffset = 0;
			for (int i=1; i<numRows; i++) {
				try {
					lineOffset = inputTextArea.getLineStartOffset(i);
					parseALine(sReturn.substring(j,lineOffset));
					j = lineOffset;
				}
				catch (BadLocationException ex) {
					outputTextArea.append("bad location exception caught\n");
				}
			}
			parseALine(sReturn.substring(lineOffset, sReturn.length())); // the last line
			System.out.println(); // Blank line at end of echoed input.
   		}
		
		void parseALine(String line) { // Parse one line of the data entered (one taxon).
			String line2 = line.replaceAll("\t",""); // Remove tab characters.
			String line3 = line2.replaceAll(" ",""); // Remove spaces.
			String line4 = line3.replaceAll("\n",""); // Remove end-of-line character.
			if (line4.length() > 0) { // an empty string if the input line was ignored
				System.out.println("input line=" + line4); // Echo the input.
				String[] elements = line4.split("/"); // Split up the entries by the slashes, if any.
				totalTaxa.add(elements); // Add the line for this taxon to the ArrayList.
			}
		}
		
		void totalTaxaToString (ArrayList<String[]> totalTaxa) { // Show the contents of the data structure on the console
			System.out.print("totalTaxa=");
			for (int i=0; i < totalTaxa.size(); i++) {
				arrayToString(totalTaxa.get(i));
				System.out.print(", ");
			}
			System.out.println("\n");
		}

		void arrayToString (String[] array) { // Show the contents of the data structure on the console
			for (int j=0; j < array.length; j++) {
				System.out.print(" " + array[j]);
			}
		}

/*------------------------------------------------------------------------------------------------------
* For the Information Statistic, build three parallel lists:
* stateList: a list of the possible states i, is used only during the first pass through the data
* listOfStateCounts: a list of tallies of taxa that allow each state (=k{i})
* listOfSigmas: a list of the sums for each state i, of 1/n{j} where j is the number of states
*       allowed for each taxon that allows state i
*-----------------------------------------------------------------------------------------------------*/
		void calculateStates() { // Loop through the taxa, processing the states specified for each.
			for (String[] taxon2 : totalTaxa) {
				int numStates = taxon2.length; // number of states for one taxon
				System.out.println("numStates= "+ numStates);
				for (String l : taxon2) { // each state for the taxon
					int index = containsString(stateList, l);
					System.out.println("index1= "+ index + " (-1 if not found)");
					if (index < 0) { // Returns -1 if not found, or the index if found.
						System.out.println ("appending new state=" + l);
						stateList.add(l);
						listOfStateCounts.add(1); // Initialize the counter to 1 (# of taxa that have this state).
						double temp = 1.0/numStates;
						listOfSigmas.add((Double)(temp)); // Initialize the SIG{i}.
						System.out.println ("numStates=" + numStates + " set sigma to " + temp);						
					}
					else { // This character state was previously found.
						int temp = listOfStateCounts.get(index);
						System.out.println ("incrementing at index=" + index + " (state " + l + ") count was " + temp);
						listOfStateCounts.set(index, ++temp);
						Double temp2 = listOfSigmas.get(index) + 1.0/numStates;
						listOfSigmas.set(index, temp2);
						System.out.println ("numStates=" + numStates + " added " + 1.0/numStates + " to sigma, giving " + temp2);
					}
				}
			}
			System.out.println("\nstatelist= " + stateList);
			System.out.println("listOfStateCounts= " + listOfStateCounts);
			System.out.println("listOfSigmas= " + listOfSigmas + "\n");
		}
		
		// Look for a character state (string) in stateList. This is used by the calculateStates method.
		int containsString(ArrayList<String> A, String s) { 
			int index = -1; // Return -1 if the character state was not found.
			for (int i = 0; i < A.size(); i ++) {
				if (s.equals(A.get(i))) { // The character state was found.
					index = i;
					break;
				}
			}
			return index;
		}
		
		// Calculate the Information Statistic and variants. Sends results directly to the output area.
		void computeInformation(int numTaxa, int totalNumStates) {
			double IntKeyH = 0.0;
			double PankhurstH = 0.0;
			double H = 0.0;
			int mindex = 0; // index for both lists
			for (int m : listOfStateCounts) {
				System.out.println("m= " + m);
				double mtemp = 1.00 * m/numTaxa;
				double m1 = Math.log(mtemp)/Math.log(2); // log2(k{i}/t)
				double m2 = Math.log(mtemp)/Math.log(totalNumStates); // log base totalNumStates
				double m3 = Math.log(mtemp)/Math.log(numTaxa); // log base numTaxa
				System.out.println("m1=" + m1 + "\nm2=" + m2 + "\nm3=" + m3);
				double n = listOfSigmas.get(mindex); // the equivalent entry in listofsigmas
				mindex++;
				double n1 = 1.00 * n/numTaxa;
				System.out.println("n1=" + n + "/" + numTaxa + " = " + n1);
				IntKeyH += -n1 * m1; // base 2 logarithms, as in IntKey
				System.out.println("IntKeyH=" + IntKeyH);
				PankhurstH += -n1 * m2; // base totalNumStates, as recommended by Pankhurst
				System.out.println("PankhurstH=" + PankhurstH);
				H += -n1 * m3;// base numtaxa (to normalize the range as the number of taxa (t) decreases)
				System.out.println("H=" + H + "\n");
			}
			outputTextArea.append(String.format("\n\nResults to two decimal places:\n"
				+ "Intkey-style information coefficient=%.2f\n"
				+ "Pankhurst's information coefficient=%.2f\n" 
				+ "Normalized information coefficient=%.2f\n", IntKeyH, PankhurstH, H));
			System.out.println(String.format("\nResults to two decimal places:\n"
				+ "Intkey-style information coefficient=%.2f\n"
				+ "Pankhurst's information coefficient=%.2f\n" 
				+ "Normalized information coefficient=%.2f\n", IntKeyH, PankhurstH, H));
		}

		// Calculate the Jaccard and Separation Coefficients. Sends results directly to the output area.
		void computeOtherCoefficients(int numTaxa, int totalNumStates) {
			// Calculate the number of possible pairs of taxa n!/2*(n-2)!=n*(n-1)/2.
			int totalNumPairs = numTaxa * (numTaxa-1) / 2;
			System.out.println("numTaxa=" + numTaxa + "\ntotalNumPairs=" + totalNumPairs + "\n");
			
			// Destructively compare taxon character-states from the top row to each of the other rows (taxa).
			String[] taxon1 = totalTaxa.get(0);
			totalTaxa.remove(0); // equivalent to "pop" in other languages
			System.out.println("popped taxon=" + taxon1);
			arrayToString(taxon1); // Display the contents of the array on the console
			System.out.println();
			
			int numDiffs = 0;
			double JaccardSimilarity = 0.0;
			while (totalTaxa.size() > 0) {
				System.out.println("Remaining taxa at start of loop " + totalTaxa);
				for (String[] taxon : totalTaxa) {
					// The separation coefficient counts only taxa that are completely separable,
					//so if any element of taxon1's list is included in taxon's list, they don't
					//count as separable
					Boolean match = false;
					//Average pairwise Jaccard distance uses the count of elements (character states) 
					//in the intersection (between two taxa) divided by the count of elements in the union.
					ArrayList<String> unionList = new ArrayList<String>();
					ArrayList<String> intersectionList = new ArrayList<String>();
					for (String k : taxon1) {
						if (! arrayContains(unionList, k))
							unionList.add(k);
						for (String j : taxon) {
							if (j.equals(k)) {
								match = true;
								intersectionList.add(k);
							} 
							else if (! arrayContains(unionList, j)) unionList.add(j);
								
						}
					}
					System.out.println("unionList=" + unionList);
					System.out.println("intersectionList=" + intersectionList);
					double JaccardCoefficient = 1.0 * intersectionList.size() / unionList.size();
					System.out.println("JaccardCoefficient=" + JaccardCoefficient + "\n");
					JaccardSimilarity += JaccardCoefficient;
					if (! match) numDiffs++;
				}
				taxon1 = totalTaxa.get(0);
				totalTaxa.remove(0); // equivalent to "pop" in other languages
				System.out.println("popped taxon=" + taxon1);
				arrayToString(taxon1); // Display the contents of the array on the console
				System.out.println();
			}
			// Separation coefficient = total separable pairs/total possible pairs.
			double separationCoefficient = numDiffs / totalNumPairs;
			System.out.println(String.format("separationCoefficient=%.2f", separationCoefficient));
			outputTextArea.append(String.format("\nSeparation coefficient=%.2f", separationCoefficient));

			// Average pairwise Jaccard distance
			System.out.println(String.format("JaccardSimilarity=" + JaccardSimilarity));
			double JaccardDistance = 1 - (JaccardSimilarity / totalNumPairs);
			System.out.println(String.format("Average pairwise Jaccard distance=" + JaccardDistance));
			outputTextArea.append(String.format("\nAverage pairwise Jaccard distance=%.2f", JaccardDistance));
		}
		
		// Look for a character state (string) in an arrayList. This is used by the computeOtherCoefficients method.
		Boolean arrayContains (ArrayList<String> array, String s) {
			Boolean returnValue = false;
			for (String t : array) {
				if (s.equals(t)) {
					returnValue = true;
				}
			}
			return returnValue;
		}
	}
}
