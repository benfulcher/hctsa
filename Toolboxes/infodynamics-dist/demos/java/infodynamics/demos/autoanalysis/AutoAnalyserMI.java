/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2015, Joseph T. Lizier
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package infodynamics.demos.autoanalysis;

import infodynamics.measures.continuous.ChannelCalculatorCommon;
import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.gaussian.MutualInfoCalculatorMultiVariateGaussian;
import infodynamics.measures.continuous.kernel.MutualInfoCalculatorMultiVariateKernel;
import infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov;
import infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1;
import infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov2;
import infodynamics.measures.discrete.MutualInformationCalculatorDiscrete;

import javax.swing.JOptionPane;
import javax.swing.event.DocumentListener;

import java.awt.event.ActionListener;
import java.awt.event.MouseListener;

/**
 * This class provides a GUI to build a simple mutual information calculation,
 *  and supply the code to execute it.
 * 
 * 
 * @author Joseph Lizier
 *
 */
public class AutoAnalyserMI extends AutoAnalyser
	implements ActionListener, DocumentListener, MouseListener {

	/**
	 * Need serialVersionUID to be serializable
	 */
	private static final long serialVersionUID = 1L;
	
	protected static final String DISCRETE_PROPNAME_TIME_DIFF = "time difference";
	
	protected static final String CALC_TYPE_KRASKOV_ALG1 = CALC_TYPE_KRASKOV + " alg. 1";
	protected static final String CALC_TYPE_KRASKOV_ALG2 = CALC_TYPE_KRASKOV + " alg. 2";
	
	/**
	 * Constructor to initialise the GUI for MI
	 */
	protected void makeSpecificInitialisations() {
		
		// Set up the properties for MI:
		measureAcronym = "MI";
		appletTitle = "JIDT Mutual Information Auto-Analyser"; 
		
		calcTypes = new String[] {
				CALC_TYPE_DISCRETE, CALC_TYPE_GAUSSIAN,
				CALC_TYPE_KRASKOV_ALG1, CALC_TYPE_KRASKOV_ALG2,
				CALC_TYPE_KERNEL};
		unitsForEachCalc = new String[] {"bits", "nats", "nats", "nats", "bits"};
		
		// Discrete:
		discreteClass = MutualInformationCalculatorDiscrete.class;
		discreteProperties = new String[] {
				DISCRETE_PROPNAME_BASE,
				DISCRETE_PROPNAME_TIME_DIFF
		};
		discretePropertyDefaultValues = new String[] {
				"2",
				"0",
		};
		discretePropertyDescriptions = new String[] {
				"Number of discrete states available for each variable (i.e. 2 for binary)",
				"Time-lag from source to dest to consider MI across; must be >= 0 (0 for standard MI)",
		};
		
		// Continuous:
		abstractContinuousClass = MutualInfoCalculatorMultiVariate.class;
		// Common properties for all continuous calcs:
		commonContPropertyNames = new String[] {
				MutualInfoCalculatorMultiVariate.PROP_TIME_DIFF
		};
		commonContPropertiesFieldNames = new String[] {
				"PROP_TIME_DIFF"
		};
		commonContPropertyDescriptions = new String[] {
				"Time-lag from source to dest to consider MI across; must be >= 0 (0 for standard MI)"
		};
		// Gaussian properties:
		gaussianProperties = new String[] {
		};
		gaussianPropertiesFieldNames = new String[] {
		};
		gaussianPropertyDescriptions = new String[] {
		};
		// Kernel:
		kernelProperties = new String[] {
				MutualInfoCalculatorMultiVariateKernel.KERNEL_WIDTH_PROP_NAME,
				MutualInfoCalculatorMultiVariateKernel.DYN_CORR_EXCL_TIME_NAME,
				MutualInfoCalculatorMultiVariateKernel.NORMALISE_PROP_NAME,			
		};
		kernelPropertiesFieldNames = new String[] {
				"KERNEL_WIDTH_PROP_NAME",
				"DYN_CORR_EXCL_TIME_NAME",
				"NORMALISE_PROP_NAME"			
		};
		kernelPropertyDescriptions = new String[] {
				"Kernel width to be used in the calculation. <br/>If the property " +
						MutualInfoCalculatorMultiVariateKernel.NORMALISE_PROP_NAME +
						" is set, then this is a number of standard deviations; " +
						"otherwise it is an absolute value.",
				"Dynamic correlation exclusion time or <br/>Theiler window (see Kantz and Schreiber); " +
						"0 (default) means no dynamic exclusion window",
				"(boolean) whether to normalise <br/>each incoming time-series to mean 0, standard deviation 1, or not  (recommended)",
		};
		// KSG (Kraskov):
		kraskovProperties = new String[] {
				MutualInfoCalculatorMultiVariateKraskov.PROP_NORMALISE,
				MutualInfoCalculatorMultiVariateKraskov.PROP_K,
				MutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE,
				MutualInfoCalculatorMultiVariateKraskov.PROP_DYN_CORR_EXCL_TIME,
				MutualInfoCalculatorMultiVariateKraskov.PROP_NORM_TYPE,
				MutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS,
		};
		kraskovPropertiesFieldNames = new String[] {
				"MutualInfoCalculatorMultiVariateKraskov.PROP_NORMALISE",
				"MutualInfoCalculatorMultiVariateKraskov.PROP_K",
				"MutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE",
				"MutualInfoCalculatorMultiVariateKraskov.PROP_DYN_CORR_EXCL_TIME",
				"MutualInfoCalculatorMultiVariateKraskov.PROP_NORM_TYPE",
				"MutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS",
		};
		kraskovPropertyDescriptions = new String[] {
				"(boolean) whether to normalise <br/>each incoming time-series to mean 0, standard deviation 1, or not (recommended)",
				"Number of k nearest neighbours to use <br/>in the full joint kernel space in the KSG algorithm",
				"Standard deviation for an amount <br/>of random Gaussian noise to add to each variable, " +
						"to avoid having neighbourhoods with artificially large counts. <br/>" +
						"(\"false\" may be used to indicate \"0\".). The amount is added in after any normalisation.",
				"Dynamic correlation exclusion time or <br/>Theiler window (see Kantz and Schreiber); " +
						"0 (default) means no dynamic exclusion window",
				"<br/>Norm type to use in KSG algorithm between the points in each marginal space. <br/>Options are: " +
						"\"MAX_NORM\" (default), otherwise \"EUCLIDEAN\" or \"EUCLIDEAN_SQUARED\" (both equivalent here)",
				"Number of parallel threads to use <br/>in computation: an integer > 0 or \"USE_ALL\" " +
						"(default, to indicate to use all available processors)",
		};
		
	}

	/**
	 * Method to assign and initialise our continuous calculator class
	 */
	protected ChannelCalculatorCommon assignCalcObjectContinuous(String selectedCalcType) throws Exception {
		if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_GAUSSIAN)) {
			return new MutualInfoCalculatorMultiVariateGaussian();
		} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KRASKOV_ALG1)) {
			return new MutualInfoCalculatorMultiVariateKraskov1();
		} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KRASKOV_ALG2)) {
			return new MutualInfoCalculatorMultiVariateKraskov2();
		} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KERNEL)) {
			return new MutualInfoCalculatorMultiVariateKernel();
		} else {
			throw new Exception("No recognised continuous calculator selected: " +
					selectedCalcType);
		}

	}

	/**
	 * Method to assign and initialise our discrete calculator class
	 */
	protected DiscreteCalcAndArguments assignCalcObjectDiscrete() throws Exception {
		String timeDiffPropValueStr, basePropValueStr;
		try {
			timeDiffPropValueStr = propertyValues.get(DISCRETE_PROPNAME_TIME_DIFF);
		} catch (Exception ex) {
			JOptionPane.showMessageDialog(this,
					ex.getMessage());
			resultsLabel.setText("Cannot find a value for property " + DISCRETE_PROPNAME_TIME_DIFF);
			return null;
		}
		try {
			basePropValueStr = propertyValues.get(DISCRETE_PROPNAME_BASE);
		} catch (Exception ex) {
			JOptionPane.showMessageDialog(this,
					ex.getMessage());
			resultsLabel.setText("Cannot find a value for property " + DISCRETE_PROPNAME_BASE);
			return null;
		}
		int timeDiff = Integer.parseInt(timeDiffPropValueStr);
		int base = Integer.parseInt(basePropValueStr);
		
		return new DiscreteCalcAndArguments(
				new MutualInformationCalculatorDiscrete(base, timeDiff),
				base,
				base + ", " + timeDiff);
	}

	protected void setObservations(ChannelCalculatorCommon calc,
			double[] source, double[] dest) throws Exception {
		// We know this is a MutualInfoCalculatorMultiVariate
		MutualInfoCalculatorMultiVariate miCalc = (MutualInfoCalculatorMultiVariate) calc;
		miCalc.setObservations(source, dest);
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new AutoAnalyserMI();
	}
}
